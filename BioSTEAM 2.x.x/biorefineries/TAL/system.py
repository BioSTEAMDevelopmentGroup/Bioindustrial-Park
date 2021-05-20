#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

Naming conventions:
    D = Distillation column
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specificiation, not physical units, but for adjusting streams

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Facilities

@author: sarangbhagwat
"""


# %% Setup


import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from math import exp as math_exp
from biosteam import main_flowsheet as F
from copy import deepcopy
from biosteam import System
from thermosteam import Stream
from biorefineries.cornstover import CellulosicEthanolTEA
from biorefineries.TAL import units, facilities
from biorefineries.TAL._process_specification import ProcessSpecification
from biorefineries.TAL.process_settings import price, CFs
from biorefineries.TAL.utils import find_split, splits_df, baseline_feedflow
from biorefineries.TAL.chemicals_data import TAL_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.TAL.tea import TALTEA
# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2
bst.speed_up()
flowsheet = bst.Flowsheet('TAL')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(TAL_chemicals)

# %% Utils
R = 8.314
TAL_Hm = 30883.66976 # by Dannenfelser-Yalkowsky method
TAL_Tm = TAL_chemicals['TAL'].Tm # 458.15 K
TAL_c = 6056.69421768496 # fitted parameter
TAL_c_by_R = TAL_c/R
TAL_Hm_by_R = TAL_Hm/R

def get_TAL_solubility_in_water(T): # mol TAL : mol (TAL+water)
    return math_exp(-(TAL_Hm_by_R) * (1/T - 1/TAL_Tm))/math_exp(TAL_c_by_R/T) 

def get_mol_TAL_dissolved(T, mol_water):
    TAL_x = get_TAL_solubility_in_water(T)
    return mol_water*TAL_x/(1-TAL_x)
# %% 

# =============================================================================
# Feedstock
# =============================================================================

feedstock = Stream('feedstock',
                    baseline_feedflow.copy(),
                    units='kg/hr',
                    price=price['Feedstock'])

U101 = units.FeedstockPreprocessing('U101', ins=feedstock)

# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0


# %% 

# =============================================================================
# Pretreatment streams
# =============================================================================

# For pretreatment, 93% purity
sulfuric_acid_T201 = Stream('sulfuric_acid_T201', units='kg/hr')
# To be mixed with sulfuric acid, flow updated in SulfuricAcidMixer
water_M201 = Stream('water_M201', T=300, units='kg/hr')
    
# To be used for feedstock conditioning
water_M202 = Stream('water_M202', T=300, units='kg/hr')

# To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
water_M203 = Stream('water_M203', phase='l', T=300, P=13.*101325, units='kg/hr')

# For neutralization of pretreatment hydrolysate
ammonia_M205 = Stream('ammonia_M205', phase='l', units='kg/hr')
# To be used for ammonia addition, flow updated by AmmoniaMixer
water_M205 = Stream('water_M205', units='kg/hr')


# =============================================================================
# Pretreatment units
# =============================================================================
H_M201 = bst.units.HXutility('H_M201', ins=water_M201,
                                 outs='steam_M201',
                                 T=99.+273.15, rigorous=True)

H_M201.heat_utilities[0].heat_transfer_efficiency = 1.
def H_M201_specification():
    T201._run()
    acid_imass = T201.outs[0].imass['SulfuricAcid']
    H_M201.ins[0].imass['Water'] = acid_imass / 0.05
    # H_M201.ins[0].imass['H2SO4'] = H_M201.ins[0].imass['Water']/1000.
    H_M201._run()
H_M201.specification = H_M201_specification
# H_M201._cost = lambda: None
# H_M201._design = lambda: None
# H_M201.heat_utilities[0].heat_exchanger = None
H_M202 = bst.units.HXutility('H_M202', ins=water_M202,
                                 outs='hot_water_M202',
                                 T=99.+273.15, rigorous=True)
H_M202.heat_utilities[0].heat_transfer_efficiency = 1.
def H_M202_specification():
    U101._run()
    H_M201.run()
    M201._run()
    feedstock, acid = U101.outs[0], M201.outs[0]
    recycled_water = H201.outs[0]
    mixture_F_mass = feedstock.F_mass + acid.F_mass
    mixture_imass_water = feedstock.imass['Water'] + acid.imass['Water'] + \
        recycled_water.imass['Water']
    total_mass = (mixture_F_mass - mixture_imass_water)/M202.solid_loading
    H_M202.ins[0].imass['Water'] = total_mass - mixture_F_mass
    # H_M202.ins[0].imass['H2SO4'] = H_M202.ins[0].imass['Water']/1000.
    H_M202._run()
H_M202.specification = H_M202_specification


# Prepare sulfuric acid
get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                      feedstock_dry_mass=get_feedstock_dry_mass())

M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, H_M201-0))
    
# Mix sulfuric acid and feedstock, adjust water loading for pretreatment
M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0, H_M202-0, ''))

# Mix feedstock/sulfuric acid mixture and steam
# M203 = units.SteamMixer('M203', ins=(M202-0, water_M203), P=5.5*101325)
M203 = bst.units.SteamMixer('M203', ins=(M202-0, water_M203), P=5.5*101325)
M203.heat_utilities[0].heat_transfer_efficiency = 1.
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0, outs=('R201_g', 'R201_l'))

# Pump bottom of the pretreatment products to the oligomer conversion tank
T202 = units.BlowdownTank('T202', ins=R201-1)
T203 = units.OligomerConversionTank('T203', ins=T202-0)
F201 = units.PretreatmentFlash('F201', ins=T203-0,
                               outs=('F201_waste_vapor', 'F201_to_fermentation'),
                               P=101325, Q=0)

M204 = bst.units.Mixer('M204', ins=(R201-0, F201-0))
H201 = bst.units.HXutility('H201', ins=M204-0,
                                 outs='condensed_pretreatment_waste_vapor',
                                 V=0, rigorous=True)

# Neutralize pretreatment hydrolysate
M205 = units.AmmoniaMixer('M205', ins=(ammonia_M205, water_M205))
def update_ammonia_and_mix():
    hydrolysate = F201.outs[1]
    # Load 10% extra
    ammonia_M205.imol['NH4OH'] = (2*hydrolysate.imol['H2SO4']) * 1.1
    M205._run()
M205.specification = update_ammonia_and_mix

T204 = units.AmmoniaAdditionTank('T204', ins=(F201-1, M205-0))
P201 = units.HydrolysatePump('P201', ins=T204-0)


# %% 

# =============================================================================
# Conversion streams
# =============================================================================

# Flow and price will be updated in EnzymeHydrolysateMixer
enzyme = Stream('enzyme', units='kg/hr', price=price['Enzyme'])
# Used to adjust enzymatic hydrolysis solid loading, will be updated in EnzymeHydrolysateMixer
enzyme_water = Stream('enzyme_water', units='kg/hr')

# Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
CSL = Stream('CSL', units='kg/hr')
# Lime for neutralization of produced acid
# fermentation_lime = Stream('fermentation_lime', units='kg/hr')

# For diluting concentrated, inhibitor-reduced hydrolysate
dilution_water = Stream('dilution_water', units='kg/hr')


# =============================================================================
# Conversion units
# =============================================================================

# Cool hydrolysate down to fermentation temperature at 50°C
H301 = bst.units.HXutility('H301', ins=P201-0, T=50+273.15)

# Mix enzyme with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, enzyme, enzyme_water))

# Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
M302 = bst.units.Mixer('M302', ins=(M301-0, ''))


# Saccharification and Cofermentation
# R301 = units.SaccharificationAndCoFermentation('R301', 
#                                                ins=(M302-0, CSL),
#                                                outs=('fermentation_effluent', 
#                                                      'sidedraw'))

# Saccharification
R301 = units.Saccharification('R301', 
                                ins=M302-0,
                                outs=('saccharification_effluent', 
                                      'sidedraw'))

M303 = bst.units.Mixer('M303', ins=(R301-0, ''))
M303_P = units.TALPump('M303_P', ins=M303-0)
# Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
S301_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
S301_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
S301_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
S301 = units.CellMassFilter('S301', ins=M303-0, outs=('solids', ''),
                            moisture_content=0.35,
                            split=find_split(S301_index,
                                              S301_cell_mass_split,
                                              S301_filtrate_split,
                                              chemical_groups))

# S302 = bst.units.Splitter('S302', ins=S301-1, outs=('to_cofermentation', 
#                                                     'to_evaporator'),
#                           split=0.2)



# F301 = bst.units.MultiEffectEvaporator('F301', ins=M304-0, outs=('F301_l', 'F301_g'),
#                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.9695)
# F301_H = bst.units.HXutility('F301_H', ins=F301-0, T=30+273.15)
# F301_H_P = units.TALPump('F301_H_P', ins=F301_H-0)



F301 = bst.units.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.7)
# # F301.V = 0.797 for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)


F301_P = units.TALPump('F301_P', ins=F301-0)


def adjust_M304_water():
    M304.ins[1].imol['Water'] = (M304.water_multiplier - 1) * M304.ins[0].imol['Water']
    M304._run()
    
M304 = bst.units.Mixer('M304', ins=(F301_P-0, dilution_water, ''))
M304.water_multiplier = 4.
M304.specification = adjust_M304_water
M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15)
M304_H_P = units.TALPump('M304_H_P', ins=M304_H-0)



# Cofermentation
R302 = units.CoFermentation('R302', 
                                ins=('', M304_H_P-0, CSL),
                                outs=('fermentation_effluent', 'CO2'))


# ferm_ratio is the ratio of conversion relative to the fermenter
R303 = units.SeedTrain('R303', ins=R301-1, outs=('seed',), ferm_ratio=0.9)

T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-M302)


# %% 

# =============================================================================
# Separation streams
# =============================================================================

# This flow will be automatically updated in CellMassFilter
# separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')

# # # To be mixed with sulfuric acid, will be updated in SulfuricAdditionTank
# # separation_acid_water = Stream('separation_acid_water', units='kg/hr')

# separation_DPHP = Stream('DPHP', DPHP =feedstock_dry_mass*22.1/1000*0.93,
#                                     H2O=feedstock_dry_mass*22.1/1000*0.07, units='kg/hr')

# # Ethanol for esterification reaction, will be updated in the EsterificationReactor
# separation_ethanol = Stream('separation_ethanol', Ethanol=feedstock_dry_mass*22.1/1000*0.93,
#                                     H2O=feedstock_dry_mass*22.1/1000*0.07, units='kg/hr')

# For ester hydrolysis
# separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')

Hexanol = Stream('Hexanol', units = 'kg/hr')
Heptane = Stream('Heptane', units = 'kg/hr')
Toluene = Stream('Toluene', units = 'kg/hr')
Hydrogen = Stream('Hydrogen', units = 'kg/hr')
KOH = Stream('KOH', units = 'kg/hr')
HCl = Stream('HCl', units = 'kg/hr')

# =============================================================================
# Separation units
# =============================================================================

# # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
S401 = bst.units.SolidsCentrifuge('S401', ins=R302-0, outs=('cell_mass', ''),
                            # moisture_content=0.50,
                            split=find_split(S401_index,
                                              S401_cell_mass_split,
                                              S401_filtrate_split,
                                              chemical_groups), solids =\
                                ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                  'Ash', 'Arabinan', 'Galactan', 'Mannan'])

def S401_TAL_split_spec():
    S401._run()
    S401_ins_0 = S401.ins[0]
    TOT_TAL = S401_ins_0.imol['TAL']
    dissolved_TAL = get_mol_TAL_dissolved(S401_ins_0.T, S401_ins_0.imol['Water'])
    
    S401.outs[0].imol['TAL'] = TOT_TAL - dissolved_TAL # crystallized TAL
    S401.outs[1].imol['TAL'] = dissolved_TAL

S401.specification = S401_TAL_split_spec

# S402 = units.TAL_Separation('S402', ins = S401-1, outs = ('', '')

# def adjust_S402():
#     S402._run()
#     TAL_mass = deepcopy(S402.ins[0].imass['TAL'])
#     S402.outs[0].mol = np.zeros(len(TAL_chemicals))
#     S402.outs[0].imass['TAL'] = TAL_mass
#     S402.outs[1] = S402.ins[0].copy()
#     S402.outs[1].imass['TAL'] = 0
#     S402-1-0-M501
    

# S402.specification = adjust_S402

M402 = bst.units.Mixer('M402', ins = (Toluene, 'recycled_toluene'), outs = 'toluene_solvent')
def adjust_M402_solvent():
    mixer = M402
    solvent = 'Toluene'
    ref_unit = S401
    reqd_solvent_frac = 0.8
    reqd_solvent_mass = reqd_solvent_frac * ref_unit.outs[1].F_mass
    mixer.ins[0].imass[solvent] = reqd_solvent_mass - mixer.ins[1].imass[solvent]
    mixer._run()
M402.specification = adjust_M402_solvent

S402 = bst.units.MultiStageMixerSettlers('S402', ins = (S401-1, M402-0), 
                                         outs = ('raffinate', 'extract'),
        partition_data={
        'K': np.array([.0001, 10000, .1, .1, .1, .1, 10000, 10000]),
        'IDs': ('Toluene', 'Water', 'Cellobiose', 'SolubleLignin', 'Enzyme', 
                'VitaminA', 'Glucose', 'TAL'),
        'phi' : 0.5,
        },
        N_stages = 2)
def S402_spec():
    M402.run()
    S402._run()
S402.specification = S402_spec

F402 = bst.units.MultiEffectEvaporator('F402', ins=S402-1, outs=('F402_l', 'F402_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.95)
                                            # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)


r_S402 = bst.units.Splitter('r_S402', ins = F402-0, 
                            outs = ('recycled_toluene', 'isolated_nonpolar_solids'), split = 0.5)
r_S402.line = 'Drum dryer'

def adjust_r_S402_recovery():
    splitter = r_S402
    solvent = 'Toluene'
    # solventobj = Toluene
    recovery = 0.97
    
    instream = splitter.ins[0]
    out0 = splitter.outs[0]
    out1 = splitter.outs[1]
    out0.T = out1.T = instream.T
    out0.P = out1.P = instream.P
    totsolvent = instream.imol[solvent]
    out0.imol[solvent] = recovery * totsolvent
    out1.copy_like(instream)
    out1.imol[solvent] = (1 - recovery) * totsolvent
    hu = bst.HeatUtility()
    tempstream = Stream('tempstream', Toluene = recovery * totsolvent, units = 'kmol/hr')
    tempstream.vle(T = instream.T, P = instream.P)
    H1 = tempstream.H
    tempstream.vle(T = TAL_chemicals[solvent].Tb, P = out0.P)
    H2 = tempstream.H
    duty = H2 - H1
    hu(duty, instream.T, TAL_chemicals[solvent].Tb)
    splitter._N_heat_utilities=1
    splitter.heat_utilities = tuple([hu])
    
r_S402.specification = adjust_r_S402_recovery

pre_M402 = bst.units.Mixer('pre_M402', ins=(F402-1, r_S402-0), outs='mixed_recycled_toluene')
pre_M402-0-1-M402



M403 = bst.units.Mixer('M403', ins = (Heptane, 'recycled_heptane'), outs = 'heptane_solvent')
def adjust_M403_solvent():
    mixer = M403
    solvent = 'Heptane'
    ref_unit = S402
    reqd_solvent_frac = 0.8
    reqd_solvent_mass = reqd_solvent_frac * ref_unit.outs[1].F_mass
    mixer.ins[0].imass[solvent] = reqd_solvent_mass - mixer.ins[1].imass[solvent]
    mixer._run()
M403.specification = adjust_M403_solvent

S403 = bst.units.MultiStageMixerSettlers('S403', ins = (S402-0, M403-0),
                                         outs = ('raffinate', 'extract'),
        partition_data={
        'K': np.array([.0001, 10000, .1, .1, .1, .1, .1, 10000]),
        'IDs': ('Heptane', 'Water', 'Cellobiose', 'SolubleLignin', 'Enzyme', 'VitaminD2', 'Glucose', 'TAL'),
        'phi' : 0.5,
        },
        N_stages = 2)
def S403_spec():
    M403.run()
    S403._run()
S403.specification = S403_spec


F403 = bst.units.MultiEffectEvaporator('F403', ins=S403-1, outs=('F403_l', 'F403_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.90)
                                            # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)


r_S403 = bst.units.Splitter('r_S403', ins = F403-0,
                            outs = ('recycled_heptane', 'isolated_nonpolar_solids'), split = 0.5)
r_S403.line = 'Drum dryer'

def adjust_r_S403_recovery():
    splitter = r_S403
    solvent = 'Heptane'
    # solventobj = Heptane
    recovery = 0.97
    instream = splitter.ins[0]
    out0 = splitter.outs[0]
    out1 = splitter.outs[1]
    out0.T = out1.T = instream.T
    out0.P = out1.P = instream.P
    totsolvent = instream.imol[solvent]
    out0.imol[solvent] = recovery * totsolvent
    out1.copy_like(instream)
    out1.imol[solvent] = (1 - recovery) * totsolvent
    hu = bst.HeatUtility()
    tempstream = Stream('tempstream', Heptane = recovery * totsolvent, units = 'kmol/hr')
    tempstream.vle(T = instream.T, P = instream.P)
    H1 = tempstream.H
    tempstream.vle(T = TAL_chemicals[solvent].Tb, P = out0.P)
    H2 = tempstream.H
    duty = H2 - H1
    hu(duty, instream.T, TAL_chemicals[solvent].Tb)
    splitter._N_heat_utilities=1
    splitter.heat_utilities = tuple([hu])
    
r_S403.specification = adjust_r_S403_recovery

pre_M403 = bst.units.Mixer('pre_M403', ins=(F403-1, r_S403-0), outs='mixed_recycled_toluene')
pre_M403-0-1-M403



M404 = bst.units.Mixer('M404', ins = (Hexanol, 'recycled_hexanol'), outs = 'hexanol_solvent')
def adjust_M404_solvent():
    mixer = M404
    solvent = 'Hexanol'
    ref_unit = S403
    reqd_solvent_frac = 1.
    reqd_solvent_mass = reqd_solvent_frac * ref_unit.outs[1].F_mass
    mixer.ins[0].imass[solvent] = reqd_solvent_mass - mixer.ins[1].imass[solvent]
    mixer._run()
M404.specification = adjust_M404_solvent

S404 = bst.units.MultiStageMixerSettlers('S404', ins = (S403-0, M404-0),
                                         outs = ('raffinate', 'extract'),
        partition_data={
        'K': np.array([1/7, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 1/10,
                       1/10, 1/10, 1/10, 1/10, 1/10, 1/10, 10000, .0001, 10000]),
        'IDs': ('TAL', 'Glucose', 'GlucoseOligomer', 'Xylose', 'XyloseOligomer', 
                'Protein', 'HMF', 'Mannose', 'Galactose', 'GalactoseOligomer',
                'Arabinose', 'ArabinoseOligomer', 'Furfural', 'AceticAcid', 'FermMicrobe',
                'Cellobiose', 'Hexanol', 'Water'),
        'phi' : 0.5,
        },
        N_stages = 4)


S405 = units.Adsorption_and_Centrifugation('S405', ins = S404-1, outs = ('diluteTAL', 'polarcompounds'))

def S405_spec():
    S405._run()
    sep_factor = 10
    polars = ('Glucose', 'GlucoseOligomer', 'Xylose', 'XyloseOligomer', 
                'Protein', 'HMF', 'Mannose', 'Galactose', 'GalactoseOligomer',
                'Arabinose', 'ArabinoseOligomer', 'Furfural', 'AceticAcid', 'FermMicrobe',
                'Cellobiose', 'Water')
    instream = S405.ins[0]
    out0 = S405.outs[0]
    out1 = S405.outs[1]
    out0.copy_like(instream)
    for polar in polars:
        molpolar = instream.imol[polar]
        out0.imol[polar] = instream.imol[polar]/sep_factor
        out1.imol[polar] = molpolar - out0.imol[polar]
S405.specification = S405_spec    

R401 = units.HydrogenationReactor('R401', ins = (S405-0, '', Hydrogen), outs = 'HMTHP',
                                  vessel_material='Stainless steel 316',)


R402 = units.DehydrationRingOpeningReactor('R402', ins = (R401-0, ''), outs = 'SA',
                                           vessel_material='Stainless steel 316',
                                           tau = 12)

R403 = units.HydrolysisReactor('R403', ins = (R402-0, '', KOH), outs = 'KSA',
                                           vessel_material='Stainless steel 316')

splits_S406 = np.zeros(len(TAL_chemicals))

splits_S406[TAL_chemicals.index('KSA')] = 0.98
splits_S406[TAL_chemicals.index('Water')] = 0.2

S406 = bst.units.SolidsCentrifuge('S406', ins=R403-0, outs=('K_sorbate', ''),
                            # moisture_content=0.50,
                            split=splits_S406,
                            solids = ['KSA'])
# !!! splits for moisture content (negligible in feed), hexanol content 
def S406_spec():
    try:
        S406._run()
    except:
        moisture_content = S406.moisture_content
        # S406.ins[0].imol['Water'] = 0.
        S406.moisture_content/= 5.
        S406._run()
        S406.moisture_content = moisture_content

S406.specification = S406_spec

S406-1-1-M404

S407 = units.Crystallization_Decantation('S407', ins = (S406-0, HCl, ''), outs = ('wet_SorbicAcid_crystals', 'KCl'))

def S407_spec():
    S407._run()
    # S408.run()
S407.specification = S407_spec

R404 = units.HClKOHRecovery('R404', ins = (S407-1, 'water'),
                            outs = ('HCl_recycle', 'KOH_recycle'))


R404-0-2-S407
R404-1-1-R403

S408 = bst.units.Flash('S408', ins = S407-0, outs = ('water', 'SorbicAcid_crystals'),
                       V = 1., P = 101325)
# !!! In vacuum


# def S408_spec():
#     instream = S408.ins[0]
#     S408.V = 0.99*instream.imol['H2O']/instream.F_mol
#     S408._run()
# S408.specification = S408_spec
S408.line = 'Drying'

# S415 = units.TAL_Separation('S415', ins = S407-0, outs = ('', ''))
# def adjust_S415():
#     S415._run()
#     SA_mass = deepcopy(S415.ins[0].imass['SA'])
#     S415.outs[0].mol = np.zeros(len(TAL_chemicals))
#     S415.outs[0].imass['SA'] = SA_mass
#     S415.outs[1] = S415.ins[0].copy()
#     S415.outs[1].imass['SA'] = 0
#     # S415-1-1-M404
    

# S415.specification = adjust_S415

# %% 

# =============================================================================
# Wastewater treatment streams
# =============================================================================

# For aerobic digestion, flow will be updated in AerobicDigestion
air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')

# To neutralize nitric acid formed by nitrification in aerobic digestion
# flow will be updated in AerobicDigestion
# The active chemical is modeled as NaOH, but the price is cheaper than that of NaOH
aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20+273.15, P=2*101325,
                          price=price['Caustics'])

# =============================================================================
# Wastewater treatment units
# =============================================================================

# Mix waste liquids for treatment
M501 = bst.units.Mixer('M501', ins=(F301-1,r_S402-1, r_S403-1, S404-0))

# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)

R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                outs=('biogas', 'anaerobic_treated_water', 
                                      'anaerobic_sludge'),
                                reactants=soluble_organics,
                                split=find_split(splits_df.index,
                                                 splits_df['stream_611'],
                                                 splits_df['stream_612'],
                                                 chemical_groups),
                                T=35+273.15)

get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185

# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
                              outs=('aerobic_vent', 'aerobic_treated_water'),
                              reactants=soluble_organics,
                              ratio=get_flow_tpd()/2205)

# Membrane bioreactor to split treated wastewater from R502
S501 = bst.units.Splitter('S501', ins=R502-1, outs=('membrane_treated_water', 
                                                    'membrane_sludge'),
                          split=find_split(splits_df.index,
                                           splits_df['stream_624'],
                                           splits_df['stream_625'],
                                           chemical_groups))

S501.line = 'Membrane bioreactor'

# Recycled sludge stream of memberane bioreactor, the majority of it (96%)
# goes to aerobic digestion and the rest to sludge holding tank then to BT
S502 = bst.units.Splitter('S502', ins=S501-1, outs=('to_aerobic_digestion', 
                                                    'to_boiler_turbogenerator'),
                          split=0.96)

M503 = bst.units.Mixer('M503', ins=(S502-0, 'centrate'), outs=1-M502)

# Mix anaerobic and 4% of membrane bioreactor sludge
M504 = bst.units.Mixer('M504', ins=(R501-2, S502-1))

# Sludge centrifuge to separate water (centrate) from sludge
S503 = bst.units.Splitter('S503', ins=M504-0, outs=(1-M503, 'sludge'),
                          split=find_split(splits_df.index,
                                           splits_df['stream_616'],
                                           splits_df['stream_623'],
                                           chemical_groups))
S503.line = 'Sludge centrifuge'

# Reverse osmosis to treat membrane separated water
S504 = bst.units.Splitter('S504', ins=S501-0, outs=('discharged_water', 'waste_brine'),
                          split=find_split(splits_df.index,
                                           splits_df['stream_626'],
                                           splits_df['stream_627'],
                                           chemical_groups))
S504.line = 'Reverse osmosis'

# Mix solid wastes to boiler turbogeneration
M505 = bst.units.Mixer('M505', ins=(S503-1, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')


# %% 

# =============================================================================
# Facilities streams
# =============================================================================

sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
# TCP_fresh = Stream('TCP_fresh',  price=price['TCP'])

ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
# lime_fresh = Stream('lime_fresh', price=price['Lime'])

HCl_fresh = Stream('HCl_fresh', price=price['HCl'])
hexanol_fresh = Stream('hexanol_fresh', price=price['Hexanol'])
heptane_fresh = Stream('heptane_fresh', price=price['Heptane'])
toluene_fresh = Stream('toluene_fresh', price=price['Toluene'])
hydrogen_fresh = Stream('hydrogen_fresh', price=price['Hydrogen'])
KOH_fresh = Stream('KOH_fresh', price=price['KOH'])
# S401_out1_F_mass = S401.outs[1].F_mass

# if not (S401_out1_F_mass == 0):
#     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M401.ins[3].imass['Ethanol']
#     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M401.ins[3].imass['Dipotassium hydrogen phosphate']
    
# else:
# ethanol_fresh = Stream('ethanol_fresh', Ethanol = get_feedstock_dry_mass()*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
# DPHP_fresh = Stream('DPHP_fresh', DPHP = get_feedstock_dry_mass()*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

# TAL stream
TAL = Stream('TAL', units='kg/hr', price=price['TAL'])
# SA product
SA = Stream('SA', units='kg/hr', price=price['SA'])
# Acetoin product
# Acetoin = Stream('Acetoin', units='kg/hr', price=price['Acetoin'])
# # Isobutyraldehyde product
# IBA = Stream('IBA', units='kg/hr', price=price['IBA'])
# Chemicals used/generated in BT
FGD_lime = Stream('FGD_lime')
ash = Stream('ash', price=price['Ash disposal'])
boiler_chems = Stream('boiler_chems', price=price['Boiler chems'])
baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
# Supplementary natural gas for BT if produced steam not enough for regenerating
# all steam streams required by the system
natural_gas = Stream('natural_gas', price=price['Natural gas'])

# Cooling tower chemicals
cooling_tower_chems = Stream('cooling_tower_chems', price=price['Cooling tower chems'])

# 145 based on equipment M-910 (clean-in-place system) in Humbird et al.
CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_tpd()/2205, units='kg/hr')

# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including enzyme production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                      N2=0.79*1372608*get_flow_tpd()/2205,
                      O2=0.21*1372608*get_flow_tpd()/2205)

# 8021 based on stream 713 in Humbird et al.
fire_water_in = Stream('fire_water_in', 
                       Water=8021*get_flow_tpd()/2205, units='kg/hr')

# =============================================================================
# Facilities units
# =============================================================================

T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh,
                                     outs=sulfuric_acid_T201)
T601.line = 'Sulfuric acid storage tank'
# S601 = bst.units.ReversedSplitter('S601', ins=T601-0, 
#                                   outs=(pretreatment_sulfuric_acid, 
#                                         ''))
# T608 = units.TCPStorageTank('T608', ins=TCP_fresh,
#                                      outs='TCP_catalyst')
# T608-0-3-R401
# T608.line = 'Tricalcium diphosphate storage tank'
#
T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia_M205)
T602.line = 'Ammonia storage tank'

T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
T603.line = 'CSL storage tank'

# DPHP storage
#!!! Yalin suggests to use BioSTEAM's storage tank, and maybe we don't need the ConveryingBelt
# (Yalin removed that from lactic acid biorefinery)
T604 = units.DPHPStorageTank('T604', ins=hexanol_fresh)
T604.line = 'Hexanol storage tank'
T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0, outs = Hexanol)

# 7-day storage time, similar to ethanol's in Humbird et al.
T605 = units.DPHPStorageTank('T605', ins=heptane_fresh)
T605.line = 'Heptane storage tank'
T605_P = units.TALPump('T605_P', ins=T605-0, outs = Heptane)

T606 = units.DPHPStorageTank('T606', ins=toluene_fresh)
T606.line = 'Toluene storage tank'
T606_P = units.TALPump('T606_P', ins=T606-0, outs = Toluene)


T607 = units.DPHPStorageTank('T607', ins=hydrogen_fresh, outs = Hydrogen)
T607.line = 'Hydrogen storage tank'

T608 = units.DPHPStorageTank('T608', ins=HCl_fresh, outs = HCl,
                             vessel_material = 'Stainless steel')
T608.line = 'HCl storage tank'

T609 = units.DPHPStorageTank('T609', ins=KOH_fresh, outs = KOH,
                             vessel_material = 'Stainless steel')
T609.line = 'KOH storage tank'

# T607_P = units.TALPump('T607_P', ins=T607-0, outs = Hydrogen)

# Connections to ATPE Mixer
# T604_P-0-1-M401
# T605_P-0-2-M401

# 7-day storage time, similar to ethanol's in Humbird et al.
T620 = units.TALStorageTank('T620', ins=S408-1, tau=7*24, V_wf=0.9,
                                      vessel_type='Floating roof',
                                      vessel_material='Stainless steel')



T620.line = 'SAStorageTank'
T620_P = units.TALPump('T620_P', ins=T620-0, outs=SA)


# # 7-day storage time, similar to ethanol's in Humbird et al.
# T607 = units.TALStorageTank('T607', ins=D402_H-0, tau=7*24, V_wf=0.9,
#                                       vessel_type='Floating roof',
#                                       vessel_material='Stainless steel')



# T607.line = 'AcetoinStorageTank'
# T607_P = units.TALPump('T607_P', ins=T607-0, outs=Acetoin)

# # 7-day storage time, similar to ethanol's in Humbird et al.
# T608 = units.TALStorageTank('T608', ins=D403_H-0, tau=7*24, V_wf=0.9,
#                                       vessel_type='Floating roof',
#                                       vessel_material='Stainless steel')



# T608.line = 'IBAStorageTank'
# T608_P = units.TALPump('T608_P', ins=T608-0, outs=IBA)


CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                     ratio=get_flow_tpd()/2205)


FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')

#!!! M304_H uses chilled water, thus requiring CWP
CWP = facilities.CWP('CWP', ins='return_chilled_water',
                     outs='process_chilled_water')

# M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
# BT = facilities.BT('BT', ins=(M505-0, R501-0, 
#                                           FGD_lime, boiler_chems,
#                                           baghouse_bag, natural_gas,
#                                           'BT_makeup_water'),
#                                 B_eff=0.8, TG_eff=0.85,
#                                 combustibles=combustibles,
#                                 side_streams_to_heat=(water_M201, water_M202, steam_M203),
#                                 outs=('gas_emission', ash, 'boiler_blowdown_water'))

BT = bst.facilities.BoilerTurbogenerator('BT',
                                                  ins=(M505-0,
                                                      R501-0, 
                                                      'boiler_makeup_water',
                                                      'natural_gas',
                                                      'lime',
                                                      'boilerchems'), 
                                                  outs=('gas_emission', 'boiler_blowdown_water', ash,),
                                                  turbogenerator_efficiency=0.85)

# BT = bst.BDunits.BoilerTurbogenerator('BT',
#                                    ins=(M505-0, R501-0, 'boiler_makeup_water', 'natural_gas', FGD_lime, boiler_chems),
#                                    boiler_efficiency=0.80,
#                                    turbogenerator_efficiency=0.85)

# Blowdown is discharged
CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                              'CT_makeup_water'),
                   outs=('process_cooling_water', 'cooling_tower_blowdown'))

# All water used in the system, here only consider water usage,
# if heating needed, then heeating duty required is considered in BT
process_water_streams = (water_M201, water_M202, water_M203, water_M205, 
                         enzyme_water,
                         aerobic_caustic, 
                         CIP.ins[-1], BT.ins[-1], CT.ins[-1])

PWC = facilities.PWC('PWC', ins=(system_makeup_water, S504-0),
                     process_water_streams=process_water_streams,
                     recycled_blowdown_streams=None,
                     outs=('process_water', 'discharged_water'))

# Heat exchange network
HXN = bst.facilities.HeatExchangerNetwork('HXN')


# HXN = HX_Network('HXN')

# %% 

# =============================================================================
# Complete system
# =============================================================================

# TAL_sys = System('TAL_sys',
#     [
#    # Feedstock preprocessing
#       U101,
      
#    # Pretreatment
#       T201, M201, # sulfuric acid mixing and addition
#       M202, # feedstock mixing
#       M203, R201, R201_H, # pretreatment 
#       T202, T203,# blowdown and oligomer conversion
#       F201, F201_H, # pretreatment flash and waste vapor condensation
#       M204, T204, T204_P, # ammonia addition
#       M205, M205_P, # waste vapor mixing and pumping
      
#    # Conversion
#       H301, # hydrolysate cooler
#       M301, # enzyme addition
#       System('fermentation_recycle',
#         [M302, R301, # simultaneous saccharification and co-fermentation
#           R302, T301], # seed train and seed holding tank
#         recycle=T301-0), # recycle seed
      
#    # Separation
#       S401, # cell mass filter
#       R401, R401_P, # acidulation
#       T401, T401_P, # sulfuric acid addition     
#       S402, # gypsum filter
#       F401, F401_H, F401_P, # separate water
#       D401, D401_H, D401_P, # separate other volatiles
#       System('esterification_recycle',
#         [System('outer_loop_acid_and_ester_recycle',
#             [System('inner_loop_ethanol_cycle',
#                 [R402, R402_P, # esterification of lactic acid
#                   D401, D401_H, D401_P], # separate out ethanol
#                 recycle=D401_H-0), # recycle ethanol
#               D401, D401_H, D401_P, S403], # separate out acid and ester
#             recycle=S403-0), # recycle acid and ester
#           System('hydrolysis_recycle',
#                 [R403, R403_P, # hydrolysis of ester
#                   D403, D403_H, D403_P, # separate out ethanol for recylcing
#                   F402, F402_H, F402_P], # separate out volatiles
#                 recycle=F402_H-0), # recycle ester
#           ],
#           recycle=D403_H-0), # recycle ethanol
#       D405, D405_H1, D405_H2, D405_P, # final purification of the acid product
      
#    # Wastewater treatment
#       M501, # mix all wastewater streams
#       WWT_cost, # total cost of wastewater treatment process
#       R501, # anaerobic digestion
#       System('wastewater_treatment_recycle',
#         [M502, R502, # aerobic digestion
#           S501, # membrane bioreactor
#           S502, M503], # sludge centrifuge
#         recycle=M503-0), # recycle sludge
#       M504, S503, # sludge centrifuge
#       S504, # reverse osmosis
#       M505, # sludge mixer
      
#    # Facilities
#       S601, T601, # sulfuric acid storage
#       T602, # ammonia storage
#       T603, # CSL storage
#       T604, T604_P, # lime storage
#       T605, T605_P, # ethanol storage
#       T606, T606_P], # lactic acid product storage
#     # facilities=(BT, CT, PWC, CIP, ADP, FWT))
#     facilities=(HXN, BT, CT, PWC, CIP, ADP, FWT))


#!!! Yalin strongly recommends reviewing the system path or manually set up the system
# for lactic acid, the automatically created system has bugs
# TAL_sys = bst.main_flowsheet.create_system(
#     'TAL_sys', feeds=[i for i in bst.main_flowsheet.stream
#                             if i.sink and not i.source])

f = bst.main_flowsheet

# TAL_sys = System('TAL_sys', path =
#     [U101,
#       T601,
#       T201,
#       M201,
#       M202,
#       M203,
#       R201,
#       T202,
#       T203,
#       F201,
#       T602,
#       M205,
#       T204,
#       P201,
#       H301,
#       M301,
#       System('a', path =
#         [M302,
#           R301,
#           R303,
#           T301],
#         recycle=T301.outs[0]),
#       M303,
#       S301,
#       F301,
#       F301_P,
#       M304,
#       M304_H,
#       M304_H_P,
#       T603,
#       R302,
#       S401,
#       System('b', path =
#         [S402,
#           r_S402,
#           T606,
#           T606_P,
#           M402],
#         recycle=f('toluene_solvent')),
#       System('c', path = 
#         [S403,
#           r_S403,
#           T605,
#           T605_P,
#           M403],
#         recycle=f('heptane_solvent')),
#       M501,
#       WWT_cost,
#       R501,
#       System('d', path =
#         [M502,
#           R502,
#           S501,
#           S502,
#           M503,
#           M504,
#           S503],
#         recycle=M503.outs[0]),
#       M505,
#       T604,
#       T604_P,
#       M404,
#       S404,
#       S405,
#       T607,
#       R401,
#       R402,
#       R403,
#       S406,
#       S407,
#       S408,
#       T620,
#       System('e', path =
#         [T609,
#           R403,
#           S406,
#           System('f', path =
#             [T608,
#               S407,
#               R404],
#             recycle=f('HCl_recycle'))],
#         recycle=f('KOH_recycle')),
#       T620, T620_P,
#       S504,
#       M204,
#       H201,
#       System('g', path =
#         [T609,
#           R403,
#           S406,
#           System('h', path =
#             [T608,
#               S407,
#               R404],
#             recycle=f('HCl_recycle'))],
#         recycle=f('KOH_recycle')),
#       System('i', path =
#         [M502,
#           R502,
#           S501,
#           S502,
#           M503,
#           M504,
#           S503],
#         recycle=M503.outs[0]),
#       ],
#     facilities = (BT, CT, FWT, PWC, ADP, CIP))

TAL_sys = bst.main_flowsheet.create_system('TAL_sys',
                            feeds=[i for i in bst.main_flowsheet.stream
                            if i.sink and not i.source])

BT_sys = System('BT_sys', path=(BT,))


TEA_feeds = set([i for i in TAL_sys.feeds if i.price]+ \
    [i for i in TAL_sys.feeds if i.price])

TEA_products = set([i for i in TAL_sys.products if i.price]+ \
    [i for i in TAL_sys.products if i.price]+[SA])
# %%
# =============================================================================
# TEA
# =============================================================================

#!!! Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legistration)
TAL_no_BT_tea = TALTEA(
        system=TAL_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(U101, WWT_cost,
                    T601, T602, T603, T604, T604_P, 
                    T605, T605_P, T606, T606_P,
                    T607,
                    T620, T620_P,
                    CWP, CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

TAL_no_BT_tea.units.remove(BT)

# # Removed because there is not double counting anyways.
# # Removes feeds/products of BT_sys from TAL_sys to avoid double-counting
# for i in BT_sys.feeds:
#     TAL_sys.feeds.remove(i)
# for i in BT_sys.products:
#     TAL_sys.products.remove(i)

# Boiler turbogenerator potentially has different depreciation schedule
# BT_tea = bst.TEA.like(BT_sys, TAL_no_BT_tea)
# BT_tea.labor_cost = 0

# Changed to MACRS 20 to be consistent with Humbird
# BT_tea.depreciation = 'MACRS20'
# BT_tea.OSBL_units = (BT,)

TAL_tea = CellulosicEthanolTEA(system=TAL_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(U101, WWT_cost,
                    T601, T602, T603, T604, T604_P, 
                    T605, T605_P, T606, T606_P,
                    T607,
                    T620, T620_P,
                    CWP, CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=BT)

TAL_sys._TEA = TAL_tea

# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

# def get_TAL_MPSP():
#     TAL_sys.simulate()
    
#     for i in range(3):
#         TAL.price = TAL_tea.solve_price(TAL, TAL_no_BT_tea)
#     return TAL.price

def get_SA_MPSP():
    for i in range(3):
        try:
            TAL_sys.simulate()
        except:
            pass
    for i in range(3):
        SA.price = TAL_tea.solve_price(SA)
    return SA.price

def get_titer():
    return R302.outs[0].imass['TAL']/R302.outs[0].F_vol

def set_titer(titer):
    M304.water_multiplier *= get_titer()/titer
    get_SA_MPSP()
    return get_titer()

# get_SA_MPSP()

# R301 = F('R301') # Fermentor
# yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = F301,
    mixer = M304,
    reactor=R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('TAL',),
    spec_1=100,
    spec_2=0.909,
    spec_3=18.5,
    path = (M304_H, M304_H_P),
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = None,
    evaporator_pump = F301_P)

# path = (F301, R302)
# @np.vectorize
# def calculate_titer(V):
#     F301.V = V
#     for i in path: i._run()
#     return spec._calculate_titer()

# @np.vectorize   
# def calculate_MPSP(V):
#     F301.V = V
#     TAL_sys.simulate()
#     MPSP = SA.price = TAL_tea.solve_price(SA, TAL_no_BT_tea)
#     return MPSP

# vapor_fractions = np.linspace(0.20, 0.80)
# titers = calculate_titer(vapor_fractions)
# MPSPs = calculate_MPSP(vapor_fractions)
# import matplotlib.pyplot as plt
# plt.plot(vapor_fractions, titers)
# plt.show()

# plt.plot(titers, MPSPs)
# plt.show()   

# %%

# =============================================================================
# Life cycle analysis (LCA), waste disposal emission not included
# =============================================================================

# 100-year global warming potential (GWP) from material flows
LCA_streams = TEA_feeds.copy()
LCA_stream = Stream('LCA_stream', units='kg/hr')
    
def get_material_GWP():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
    # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
    return chemical_GWP.sum()/SA.F_mass

# GWP from combustion of non-biogenic carbons
get_non_bio_GWP = lambda: (natural_gas.get_atomic_flow('C'))* TAL_chemicals.CO2.MW / SA.F_mass
                            # +ethanol_fresh.get_atomic_flow('C')) \
    

# GWP from electricity
get_electricity_use = lambda: sum(i.power_utility.rate for i in TAL_sys.units)
get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
    / SA.F_mass

# CO2 fixed in lactic acid product
get_fixed_GWP = lambda: \
    SA.get_atomic_flow('C')*TAL_chemicals.CO2.MW/SA.F_mass

# carbon_content_of_feedstock = 0
get_GWP = lambda: get_material_GWP()+get_non_bio_GWP()+get_electricity_GWP() 

# Fossil energy consumption (FEC) from materials
def get_material_FEC():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
    # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
    return chemical_FEC.sum()/SA.F_mass

# FEC from electricity
get_electricity_FEC = lambda: \
    get_electricity_use()*CFs['FEC_CFs']['Electricity']/SA.F_mass

# Total FEC
get_FEC = lambda: get_material_FEC()+get_electricity_FEC()

# get_SPED = lambda: BT.system_heating_demand*0.001/SA.F_mass
SA_LHV = 31.45 # MJ/kg SA

# %% Full analysis
def simulate_and_print():
    get_SA_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_SA_MPSP():.3f}/kg')
    # print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg SA')
    # print(f'FEC is {get_FEC():.2f} MJ/kg SA or {get_FEC()/SA_LHV:.2f} MJ/MJ SA')
    # print(f'SPED is {get_SPED():.2f} MJ/kg SA or {get_SPED()/SA_LHV:.2f} MJ/MJ SA')
    # print('--------------------\n')

# simulate_and_print()
TAL_sys.simulate()

# %% 

# =============================================================================
# For Monte Carlo and analyses
# =============================================================================

TAL_sub_sys = {
#     'feedstock_sys': (U101,),
#     'pretreatment_sys': (T201, M201, M202, M203, 
#                          R201, R201_H, T202, T203,
#                          F201, F201_H,
#                          M204, T204, T204_P,
#                          M205, M205_P),
#     'conversion_sys': (H301, M301, M302, R301, R302, T301),
    # 'separation_sys': (S401, M401, M401_P,
    #                     S402, 
    #                     # F401, F401_H, F401_P,
    #                     D401, D401_H, D401_P, S403,
    #                     M402_P, S403,
    #                     D403, D403_H, D403_P,
    #                     M501,
    #                     T606, T606_P, T607, T607_P)
                        # F402, F402_H, F402_P,
                        # D405, D405_H1, D405_H2, D405_P,
                        # M401, M401_P)
#     'wastewater_sys': (M501, WWT_cost, R501,
#                        M502, R502, S501, S502, M503,
#                        M504, S503, S504, M505),
#     'HXN': (HXN,),
#     'BT': (BT,),
#     'CT': (CT,),
#     'other_facilities': (T601, S601,
#                          T602, T603,
#                          T604, T604_P,
#                          T605, T605_P,
#                          T606, T606_P,
#                          PWC, CIP, ADP, FWT)
    }

# for unit in sum(TAL_sub_sys.values(), ()):
#     if not unit in TAL_sys.units:
#         print(f'{unit.ID} not in TAL_sys.units')

# for unit in TAL_sys.units:
#     if not unit in sum(TAL_sub_sys.values(), ()):
#         print(f'{unit.ID} not in TAL_sub_sys')









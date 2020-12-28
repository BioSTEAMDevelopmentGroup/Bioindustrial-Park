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
from biosteam import main_flowsheet as F

from biosteam import System
from thermosteam import Stream
from BDO import units, facilities
from BDO._process_specification import ProcessSpecification
from BDO.process_settings import price, CFs
from BDO.utils import find_split, splits_df, baseline_feedflow
from BDO.chemicals_data import BDO_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from BDO.tea import BDOTEA

# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2
bst.speed_up()
flowsheet = bst.Flowsheet('BDO')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(BDO_chemicals)



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
water_M201 = Stream('water_M201', T=114+273.15, units='kg/hr')

# To be used for feedstock conditioning
water_M202 = Stream('water_M202', T=95+273.15, units='kg/hr')

# To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
steam_M203 = Stream('steam_M203', phase='g',T=268+273.15, P=13*101325, units='kg/hr')

# For neutralization of pretreatment hydrolysate
ammonia_M205 = Stream('ammonia_M205', phase='l', units='kg/hr')
# To be used for ammonia addition, flow updated by AmmoniaMixer
water_M205 = Stream('water_M205', units='kg/hr')


# =============================================================================
# Pretreatment units
# =============================================================================

# Prepare sulfuric acid
get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                      feedstock_dry_mass=get_feedstock_dry_mass())

M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, water_M201))

# Mix sulfuric acid and feedstock, adjust water loading for pretreatment
M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0, water_M202))

# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, steam_M203), P=5.5*101325)
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
M303_P = units.BDOPump('M303_P', ins=M303-0)
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
# F301_H_P = units.BDOPump('F301_H_P', ins=F301_H-0)

F301 = bst.units.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                       P = (101325, 73581, 50892, 32777, 20000), V = 0.793)
# F301.V = 0.797 for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)


F301_P = units.BDOPump('F301_P', ins=F301-0)


def adjust_M304_water():
    M304.ins[1].imol['Water'] = (M304.water_multiplier - 1) * M304.ins[0].imol['Water']
    M304._run()
    
M304 = bst.units.Mixer('M304', ins=(F301_P-0, dilution_water, ''))
M304.water_multiplier = 1.
M304.specification = adjust_M304_water
M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15)
M304_H_P = units.BDOPump('M304_H_P', ins=M304_H-0)



# Cofermentation
R302 = units.CoFermentation('R302', 
                                ins=('', M304_H_P-0, CSL),
                                outs=('fermentation_effluent', 'CO2'))

# R302 = units.CoFermentation_original('R302', 
#                                 ins=('', M304_H_P-0, CSL),
#                                 outs=('fermentation_effluent', 'CO2'))

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

# separation_H2 = Stream('separation_H2', price = price['H2'])

# =============================================================================
# Separation units
# =============================================================================

# Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
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



M401 = bst.units.Mixer('M401', ins=(S401-1, '', ''))

def adjust_M401_DPHP():
    feed, feed_DPHP, recycle_DPHP = all_ins = M401.ins[0:3]
    # print(S401.outs[0].F_mass)
    # print(recycle_ethanol.imass['Ethanol'])
    tot_mass = sum([i.F_mass - i.imass['DPHP'] - i.imass['Ethanol'] for i in all_ins])
    # M401.ins[2].imass['Ethanol'] = max(0, 0.24 * tot_mass - recycle_ethanol.imass['Ethanol'])
    M401.ins[1].imass['Dipotassium hydrogen phosphate'] = max(0, 0.24 * tot_mass - recycle_DPHP.imass['DPHP'] - feed.imass['DPHP'])
    
    # print(M401.ins[2].imass['Ethanol'])
    M401._run()


M401.specification = adjust_M401_DPHP
M401_P = units.BDOPump('M401_P', ins=M401-0, outs='mixed_stream')

M401_P_H = bst.HXutility('M401_P_H', ins = M401_P-0, V = 0, rigorous = True)


M401_b = bst.units.Mixer('M401_b', ins=('', ''))

def adjust_M401_b_ethanol():
    # import pdb
    # pdb.set_trace()
    M401.specification()
    feed_ethanol, recycle_ethanol  = M401_b.ins 
    feed, feed_DPHP, recycle_DPHP = all_ins = M401.ins[0:3] # yes, M401
    tot_mass = sum([i.F_mass - i.imass['DPHP'] - i.imass['Ethanol'] for i in all_ins])
    M401_b.ins[0].imass['Ethanol'] = max(0, 0.25 * tot_mass - recycle_ethanol.imass['Ethanol'])
    M401_b._run()
    
M401_b.specification = adjust_M401_b_ethanol

S402 = bst.units.MultiStageMixerSettlers('S402', ins = (M401_P_H-0, M401_b-0),
                                         # model = 'partition_coefficients',
                                         partition_data={
        'K': np.array([1/28.34, 1/0.046, 1/10000, 10000, 1/28.34, 1/0.046]),
        'IDs': ('2,3-Butanediol', 'Glucose', 'Ethanol', 'Water', 'Acetoin', 'Xylose'),
        'phi' : 0.5,
        },
        N_stages = 5)

def adjust_S402_split():
    M401_b.specification()
    S402._run()

    
S402.specification = adjust_S402_split                                    
    # mat_a = np.array([[L, (1-L)], [Ks[0], -1]])
    # mat_b = np.array([[zs[0]], [0]])
    # mat_x = np.linalg.solve(mat_a, mat_b)
    
    # S402.isplit['2,3-Butanediol'] = L[0]
    # S402.isplit['Glucose'] = L[1]
    # S402.isplit['Xylose'] =  L[1]
    # S402.isplit['Ethanol'] = L[2]
    # S402.isplit['Water'] = L[3]
    # m_etoh = S402.ins[0].imass['Ethanol']
    # m_h2o = S402.ins[0].imass['Water']
    # S402.isplit['2,3-Butanediol'] = (k_23bdo * m_etoh / m_h2o)/(1+k_23bdo * m_etoh / m_h2o)
    # S402.isplit['Glucose'] = (k_glucose * m_etoh / m_h2o)/(1+k_glucose * m_etoh / m_h2o)
    # S402.isplit['Xylose'] =  S402.isplit['Glucose']
    # S402.isplit['Water'] = 0.01
    # S402.isplit['Ethanol'] = 0.99
    
    
    

# split = [0.001 for i in range(len(BDO_chemicals))]
# split[BDO_chemicals.index('Dipotassium hydrogen phosphate')] = 0

# S402 = bst.units.LiquidsSplitSettler('S402', ins = M401_P_H-0, split=split)
# S402.specification = adjust_S402_split


# Separate out the majority of water,
# no need to include agitator thus using biosteam Flash
# F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'),
#                                     # LHK=('AceticAcid', '2,3-Butanediol'),
#                                     # is_divided=True,
#                                     # product_specification_format='Recovery',
#                                     # Lr=0.8, Hr=0.8, k=1.2,
#                                     T = 379, P = 101325,
#                                     vessel_material = 'Stainless steel 316')


F401 = bst.units.Flash('F401', ins=S402-0, outs=('F401_g', 'F401_l'),
                                    # LHK=('AceticAcid', '2,3-Butanediol'),
                                    # is_divided=True,
                                    # product_specification_format='Recovery',
                                    # Lr=0.8, Hr=0.8, k=1.2,
                                    T = 379, P = 101325,
                                    vessel_material = 'Stainless steel 316')


# def adjust_F401_V():
#     F401.V = F401.ins[0].imol['H2O']/F401.ins[0].F_mol
#     F401._run()
    
# F401.specification = adjust_F401_V

# # Condense waste vapor for recycling
F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)
F401_P = units.BDOPump('F401_P', ins=F401-1)



# Placeholder separation operation for DPHP recovery; consider gravimetric separation
# (rho_DPHP ~= 2.5 g/cm3, rho_Glucose ~= 1.5 g/cm3)

S403 = bst.units.Splitter('S403', ins = F401_P-0,
                          outs = ('DPHP_recycle', 'to_M501'),
                          split = 0.98)

def adjust_S403_split():
    S403._run()
    separables = ('Extract', 'Arabinose', 'Xylose', 'XyloseOligomer',
                  'Mannose', 'Glucose', 'GlucoseOligomer', 'Galactose') # DPHP is insoluble in ethanol
    extract_mol = S403.ins[0].imol[separables]
    S403.outs[0].imol[separables] = 0
    S403.outs[1].imol[separables] = extract_mol
    
S403.specification = adjust_S403_split


# DPHP + sugars recycle to fermentor
# S403-0-0-R302
# DPHP recycle to mixer-settler
S403-0-2-M401



F402 = bst.units.Flash('F402', ins=S402-1,
                                    outs=('F402_g', 'F402_l'),
                                    T = 368, P = 101325,
                                    vessel_material = 'Stainless steel 316')
# D401.process_specification = adjust_D401_V

F402_H = bst.units.HXutility('F402_H', ins=F402-0, V=0, rigorous=True)
F402_P = units.BDOPump('F402_P', ins=F402-1)

D401 = bst.units.ShortcutColumn('D401', ins=F402_P-0,
                                    outs=('D401_g', 'D401_l'),
                                    LHK=('Ethanol', 'BDO'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')
D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0, rigorous=True)
D401_P = units.BDOPump('D401_P', ins=D401-1)

# ethanol recycle
M402 = bst.units.Mixer('M402', ins = (F402_H-0, D401_H-0), outs = 'ethanol_mixed')
M402_P = units.BDOPump('M402_P', ins=M402-0, outs='ethanol_recycle')
M402_P-0-1-M401_b


D402 = bst.units.ShortcutColumn('D402', ins=D401_P-0,
                                    outs=('D402_g', 'D402_l'),
                                    LHK=('Acetoin', 'BDO'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')

D402_H = bst.units.HXutility('D402_H', ins=D402-0, V=0, rigorous=True)
D402_P = units.BDOPump('D402_P', ins=D402-1)
# M402 = bst.units.Mixer('M402', ins=(D401_H-0))

# acetoin recycle to fermentor
D402_H-0-0-R302



R401 = units.DehydrationReactor('R401', ins = (D402_P-0, ''),
                                tau = 5,
                                vessel_material='Stainless steel 316')

D403 = bst.units.ShortcutColumn('D403', ins=R401-0,
                                    outs=('IBA', 'D403_l'),
                                    LHK=('IBA', 'MEK'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')
D403_H = bst.units.HXutility('D403_H', ins=D403-0, V=0, rigorous=True)  
D403_P = units.BDOPump('D403_P', ins=D403-1)

# F402 = bst.units.MultiEffectEvaporator('F402', ins=D403_P-0,
#                                     outs=('F402_l', 'MEK'),
#                                     P = (101325, 73581, 50892, 32777, 20000), V = 0.324)
# def adjust_F402_V():
#     instream = F402.ins[0]
#     F402.V = instream.imol['MEK']/instream.F_mol
#     F402._run()
# F402.specification = adjust_F402_V
# # F402_H = bst.units.HXutility('F402_H', ins=F402-0, V=0, rigorous=True)
# F402_P = units.BDOPump('F402_P', ins=F402-1)

D404 = bst.units.ShortcutColumn('D404', ins=D403_P-0,
                                    outs=('D404_g', 'D404_l'),
                                    LHK=('MEK', 'H2O'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')
D404_H = bst.units.HXutility('D404_H', ins=D404-0, V=0, rigorous=True)
D404_P = units.BDOPump('D404_P', ins=D404-1)



D405 = bst.units.ShortcutColumn('D405', ins=D404_P-0,
                                    outs=('D404_g', 'D404_l'),
                                    LHK=('H2O', 'BDO'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')
D405_H = bst.units.HXutility('D405_H', ins=D405-0, V=0, rigorous=True)
D405_P = units.BDOPump('D405_P', ins=D405-1)



# def adjust_D405_streams():
#     D405._run()
#     D405.outs[1].imass['Xylose'] = 0
    
# D405.specification = adjust_D405_streams

# D403 = bst.units.ShortcutColumn('D403', ins=D404_P-0,
#                                     outs=('D403_g', 'D403_l'),
#                                     LHK=('Acetoin', 'BDO'),
#                                     is_divided=True,
#                                     product_specification_format='Recovery',
#                                     Lr=0.9995, Hr=0.9995, k=1.2,
#                                     vessel_material = 'Stainless steel 316')
# D403_H = bst.units.HXutility('D403_H', ins=D403-0, V=0, rigorous=True)
# D403_P = units.BDOPump('D403_P', ins=D403-1)



S404 = bst.units.Splitter('S404', ins = D405_P-0, split = 0.92)
S404-0-1-R401
# F402 = bst.units.MultiEffectEvaporator('F402', ins=D404_P-0,
#                                     outs=('H2O', 'BDO_rich'),
#                                     P = (101325, 73581, 50892, 32777, 20000), V = 0.724)
# def adjust_F402_V():
#     instream = F402.ins[0]
#     F402.V = 1.3*instream.imol['H2O']/instream.F_mol
#     F402._run()
# F402.specification = adjust_F402_V
# # F402_H = bst.units.HXutility('F402_H', ins=F402-0, V=0, rigorous=True)
# F402_P = units.BDOPump('F402_P', ins=F402-0)

R402 = units.HydrogenationReactor('R402', ins = (D403_H-0, '', ''),
                                tau = 2,
                                vessel_material = 'Stainless steel 316')


# F404 = bst.units.MultiEffectEvaporator('F404', ins=R402-0,
#                                     outs=('F404_l', 'water'),
#                                     P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
# def adjust_F404_V():
#     instream = F404.ins[0]
#     F404.V = 0.5*instream.imol['Water']/instream.F_mol
#     F404._run()
# F404.specification = adjust_F404_V

# F404_P = units.BDOPump('F404_P', ins=F404-0)

D406 = bst.units.BinaryDistillation('D406', ins=R402-0,
                                    outs=('IBA', 'IBO'),
                                    LHK=('IBA', 'Isobutanol'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.995, Hr=0.995, k=1.2,
                                    vessel_material = 'Stainless steel 316')

D406_H = bst.units.HXutility('D406_H', ins=D406-0, V=0, rigorous=True)
D406_P = units.BDOPump('D406_P', ins=D406-1)

D406_H-0-2-R402

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
M501 = bst.units.Mixer('M501', ins=(S404-1, F301-1, F401_H-0, M204-0, D405_H-0, S403-1))

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

# S401_out1_F_mass = S401.outs[1].F_mass

# if not (S401_out1_F_mass == 0):
#     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M401.ins[3].imass['Ethanol']
#     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M401.ins[3].imass['Dipotassium hydrogen phosphate']
    
# else:
ethanol_fresh = Stream('ethanol_fresh', Ethanol = get_feedstock_dry_mass()*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
DPHP_fresh = Stream('DPHP_fresh', DPHP = get_feedstock_dry_mass()*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
H2_fresh = Stream('H2_fresh', price = price['H2'])
# BDO stream
# BDO = Stream('BDO', units='kg/hr', price=price['BDO'])
# MEK product
MEK = Stream('MEK', units='kg/hr', price=price['MEK'])
# Acetoin product
Acetoin = Stream('Acetoin', units='kg/hr', price=price['Acetoin'])
# Isobutyraldehyde product
Isobutanol = Stream('Isobutanol', units='kg/hr', price=price['Isobutanol'])
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
T604 = units.DPHPStorageTank('T604', ins=DPHP_fresh)
T604.line = 'DPHP storage tank'
T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0)

# 7-day storage time, similar to ethanol's in Humbird et al.
T605 = bst.units.StorageTank('T605', ins=ethanol_fresh,
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T605.line = 'Ethanol storage tank'
T605_P = units.BDOPump('T605_P', ins=T605-0)



# Connections to ATPE Mixer
T604_P-0-1-M401
T605_P-0-0-M401_b

# 7-day storage time, similar to ethanol's in Humbird et al.
T606 = units.BDOStorageTank('T606', ins=D404_H-0, tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Stainless steel')



T606.line = 'MEKStorageTank'
T606_P = units.BDOPump('T606_P', ins=T606-0, outs=MEK)


T607 = units.HydrogenGasStorageTank('T607', ins=H2_fresh)
T607.line = 'H2 storage tank'
T607-0-1-R402

T608 = units.BDOStorageTank('T608', ins=D406_P-0, tau=7*24, V_wf=0.9,
                                      vessel_type='Floating roof',
                                      vessel_material='Stainless steel')



T608.line = 'IsobutanolStorageTank'
T608_P = units.BDOPump('T608_P', ins=T608-0, outs=Isobutanol)


CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                     ratio=get_flow_tpd()/2205)


FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')

#!!! M304_H uses chilled water, thus requiring CWP
CWP = facilities.CWP('CWP', ins='return_chilled_water',
                     outs='process_chilled_water')

# M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
BT = facilities.BT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chems,
                                          baghouse_bag, natural_gas,
                                          'BT_makeup_water'),
                                B_eff=0.8, TG_eff=0.85,
                                combustibles=combustibles,
                                side_streams_to_heat=(water_M201, water_M202, steam_M203),
                                outs=('gas_emission', ash, 'boiler_blowdown_water'))


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
process_water_streams = (water_M201, water_M202, steam_M203, water_M205, 
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

# BDO_sys = System('BDO_sys',
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
BDO_sys = bst.main_flowsheet.create_system(
    'BDO_sys', feeds=[i for i in bst.main_flowsheet.stream
                            if i.sink and not i.source])

BT_sys = System('BT_sys', path=(BT,))


TEA_feeds = set([i for i in BDO_sys.feeds if i.price]+ \
    [i for i in BT_sys.feeds if i.price])

TEA_products = set([i for i in BDO_sys.products if i.price]+ \
    [i for i in BT_sys.products if i.price]+[MEK])
# %%
# =============================================================================
# TEA
# =============================================================================

#!!! Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legistration)
BDO_no_BT_tea = BDOTEA(
        system=BDO_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(U101, WWT_cost,
                    T601, T602, T603, T604, T604_P, T605, T605_P, T606, T606_P, T607,
                    CWP, CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

BDO_no_BT_tea.units.remove(BT)

# # Removed because there is not double counting anyways.
# # Removes feeds/products of BT_sys from BDO_sys to avoid double-counting
# for i in BT_sys.feeds:
#     BDO_sys.feeds.remove(i)
# for i in BT_sys.products:
#     BDO_sys.products.remove(i)

# Boiler turbogenerator potentially has different depreciation schedule
BT_tea = bst.TEA.like(BT_sys, BDO_no_BT_tea)
BT_tea.labor_cost = 0

# Changed to MACRS 20 to be consistent with Humbird
BT_tea.depreciation = 'MACRS20'
BT_tea.OSBL_units = (BT,)

BDO_tea = bst.CombinedTEA([BDO_no_BT_tea, BT_tea], IRR=0.10)
BDO_sys._TEA = BDO_tea

# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

System.converge_method = 'fixed-point' # aitken isn't stable
System.maxiter = 1500
System.molar_tolerance = 0.1

# def get_BDO_MPSP():
#     BDO_sys.simulate()
    
#     for i in range(3):
#         BDO.price = BDO_tea.solve_price(BDO, BDO_no_BT_tea)
#     return BDO.price

def get_MEK_MPSP():
    BDO_sys.simulate()
    
    for i in range(3):
        MEK.price = BDO_tea.solve_price(MEK, BDO_no_BT_tea)
    return MEK.price

# get_MEK_MPSP()

# R301 = F('R301') # Fermentor
# yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = F301,
    mixer = M304,
    reactor=R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('BDO',),
    spec_1=100,
    spec_2=0.909,
    spec_3=18.5,
    path = (M304_H, M304_H_P),
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = R401,
    byproduct_streams = [Isobutanol],
    evaporator_pump = F301_P)

path = (F301, R302)
@np.vectorize
def calculate_titer(V):
    F301.V = V
    for i in path: i._run()
    return spec._calculate_titer()

@np.vectorize   
def calculate_MPSP(V):
    F301.V = V
    BDO_sys.simulate()
    MPSP = MEK.price = BDO_tea.solve_price(MEK, BDO_no_BT_tea)
    return MPSP

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

get_Isobutanol_GWP = lambda: - Isobutanol.F_mass * CFs['GWP_CFs']['Isobutanol']/MEK.F_mass
get_Isobutanol_FEC = lambda: - Isobutanol.F_mass * CFs['FEC_CFs']['Isobutanol']/MEK.F_mass

def get_material_GWP():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
    # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
    return chemical_GWP.sum()/MEK.F_mass

# GWP from combustion of non-biogenic carbons
get_non_bio_GWP = lambda: (natural_gas.get_atomic_flow('C')+ethanol_fresh.get_atomic_flow('C')) \
    * BDO_chemicals.CO2.MW / MEK.F_mass

# GWP from electricity
get_electricity_use = lambda: sum(i.power_utility.rate for i in BDO_sys.units)
get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
    / MEK.F_mass

# CO2 fixed in lactic acid product
get_fixed_GWP = lambda: \
    MEK.get_atomic_flow('C')*BDO_chemicals.CO2.MW/MEK.F_mass

get_GWP = lambda: get_material_GWP()+get_non_bio_GWP()+get_electricity_GWP() - get_Isobutanol_GWP()

# Fossil energy consumption (FEC) from materials
def get_material_FEC():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
    # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
    return chemical_FEC.sum()/MEK.F_mass

# FEC from electricity
get_electricity_FEC = lambda: \
    get_electricity_use()*CFs['FEC_CFs']['Electricity']/MEK.F_mass

# Total FEC
get_FEC = lambda: get_material_FEC()+get_electricity_FEC() - get_Isobutanol_FEC()

get_SPED = lambda: BT.system_heating_demand*0.001/MEK.F_mass
MEK_LHV = 31.45 # MJ/kg MEK

# %% Full analysis
def simulate_and_print():
    get_MEK_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_MEK_MPSP():.3f}/kg')
    print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg MEK')
    print(f'Non-bio GWP is {get_non_bio_GWP():.3f} kg CO2-eq/kg MEK')
    print(f'FEC is {get_FEC():.2f} MJ/kg MEK or {get_FEC()/MEK_LHV:.2f} MJ/MJ MEK')
    print(f'SPED is {get_SPED():.2f} MJ/kg MEK or {get_SPED()/MEK_LHV:.2f} MJ/MJ MEK')
    print('--------------------\n')

simulate_and_print()

# %% 

# =============================================================================
# For Monte Carlo and analyses
# =============================================================================

BDO_sub_sys = {
#     'feedstock_sys': (U101,),
#     'pretreatment_sys': (T201, M201, M202, M203, 
#                          R201, R201_H, T202, T203,
#                          F201, F201_H,
#                          M204, T204, T204_P,
#                          M205, M205_P),
#     'conversion_sys': (H301, M301, M302, R301, R302, T301),
    'separation_sys': (S401, M401, M401_P,
                        S402, 
                        # F401, F401_H, F401_P,
                        D401, D401_H, D401_P, S403,
                        M402_P, S403,
                        D403, D403_H, D403_P,
                        M501,
                        T606, T606_P)
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

# for unit in sum(BDO_sub_sys.values(), ()):
#     if not unit in BDO_sys.units:
#         print(f'{unit.ID} not in BDO_sys.units')

# for unit in BDO_sys.units:
#     if not unit in sum(BDO_sub_sys.values(), ()):
#         print(f'{unit.ID} not in BDO_sub_sys')









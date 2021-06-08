#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:15:23 2019

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

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

@author: yalinli_cabbi
"""


# %% Setup

import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import main_flowsheet as F
from biorefineries import BST222
from biosteam import System
from thermosteam import Stream
from biorefineries.BDO import units, facilities
from biorefineries.BDO._process_specification import ProcessSpecification
from biorefineries.BDO.process_settings import price
from biorefineries.BDO.utils import find_split, splits_df, baseline_feedflow
from biorefineries.BDO.chemicals_data import BDO_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.BDO.tea import BDOTEA

bst.speed_up()
flowsheet = bst.Flowsheet('BDO')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0
# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7
# Baseline cost year is 2016
bst.CE = 541.7
_labor_2007to2016 = 22.71 / 19.55

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
U101.feedstock_flow_rate = feedstock.F_mass # TODO: Check this

# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0


# %% 

# =============================================================================
# Pretreatment streams
# =============================================================================

# To be used for feedstock conditioning, flow is adjusted in PretreatmentMixer
pretreatment_feedstock_water = Stream('pretreatment_feedstock_water',
                                      T=95+273.15, units='kg/hr')

# For pretreatment, baseline is (18+4.1) mg/g dry biomass
# based on P21 in Humbird et al., 93% purity
feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O'] # TODO: Check this
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid', 
                                    H2SO4=feedstock_dry_mass*22.1/1000*0.93,
                                    H2O=feedstock_dry_mass*22.1/1000*0.07,
                                    units='kg/hr')

# Flow adjusted in SulfuricAcidMixer, stream 516 in Humbird et al.
pretreatment_acid_water = Stream('pretreatment_acid_water', T=114+273.15)

# To be added to the feedstock/sulfuric acid mixture,
# will be adjusted by the SteamMixer
pretreatment_steam = Stream('pretreatment_steam', phase='g',
                            T=268+273.15, P=13*101325,
                            Water=(3490+24534)*U101.feedstock_flow_rate/2205,
                            units='kg/hr')

# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia', units='kg/hr', phase='l')
# To be used for ammonia addition, will be updated by AmmoniaMixer
pretreatment_ammonia_water = Stream('pretreatment_ammonia_water', units='kg/hr')


# =============================================================================
# Pretreatment units
# =============================================================================

# Prepare sulfuric acid
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid,
                                      feedstock_dry_mass=feedstock.F_mass - feedstock.imass['Water'])
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))

# Mix sulfuric acid and feedstock, adjust water loading
M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0,
                                            pretreatment_feedstock_water))

# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, pretreatment_steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0,
                                       outs=('R201_g', 'R201_l'))
R201_H = bst.units.HXutility('R201_H', ins=R201-0, V=0, rigorous=True)

# Pump bottom of the pretreatment products to the oligomer conversion tank
T202 = units.BlowdownTank('T202', ins=R201-1)
T203 = units.OligomerConversionTank('T203', ins=T202-0)
F201 = units.PretreatmentFlash('F201', ins=T203-0, outs=('F201_g', 'F201_l'),
                               P=101325, Q=0)
F201_H = bst.units.HXutility('F201_H', ins=F201-0, V=0, rigorous=True)

# Neutralize pretreatment hydrolysate
def update_ammonia_and_mix():
    hydrolysate = F201.outs[1]
    # Load 5% extra
    ammonia.imol['AmmoniumHydroxide'] = (2*hydrolysate.imol['H2SO4']) * 1.05
    M204._run()
    
M204 = units.AmmoniaMixer('M204', ins=(ammonia, pretreatment_ammonia_water))
M204.specification = update_ammonia_and_mix

T204 = units.AmmoniaAdditionTank('T204', ins=(F201-1, M204-0))
T204_P = units.HydrolysatePump('T204_P', ins=T204-0)




M205 = bst.units.Mixer('M205', ins=(R201_H-0, F201_H-0))
M205_P = units.BDOPump('M205_P', ins=M205-0, outs='condensed_pretreatment_waste_vapor')


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
H301 = bst.units.HXutility('H301', ins=T204_P-0, T=50+273.15)

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

# Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
S301_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
S301_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
S301_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
S301 = units.CellMassFilter('S301', ins=R301-0, outs=('solids', ''),
                            moisture_content=0.35,
                            split=find_split(S301_index,
                                              S301_cell_mass_split,
                                              S301_filtrate_split,
                                              chemical_groups))

# S302 = bst.units.Splitter('S302', ins=S301-1, outs=('to_cofermentation', 
#                                                     'to_evaporator'),
#                           split=0.2)

# def adjust_M303_water():
#     M303.ins[1].imol['Water'] = (M303.water_multiplier - 1) * M303.ins[0].imol['Water']
#     M303._run()
    
# M303 = bst.units.Mixer('M303', ins=(S302-1, ''))
# M303.water_multiplier = 6.5
# M303.specification = adjust_M303_water


# F301 = bst.units.MultiEffectEvaporator('F301', ins=M303-0, outs=('F301_l', 'F301_g'),
#                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.9695)
# F301_H = bst.units.HXutility('F301_H', ins=F301-0, T=30+273.15)
# F301_H_P = units.BDOPump('F301_H_P', ins=F301_H-0)

F301 = bst.units.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                       P = (101325, 73581, 50892, 32777, 20000), V = 0.797)
# F301.V = 0.797 for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)


F301_P = units.BDOPump('F301_P', ins=F301-0)


def adjust_M303_water():
    M303.ins[1].imol['Water'] = (M303.water_multiplier - 1) * M303.ins[0].imol['Water']
    M303._run()
    
M303 = bst.units.Mixer('M303', ins=(F301_P-0, dilution_water))
M303.water_multiplier = 1.001
M303.specification = adjust_M303_water
M303_H = bst.units.HXutility('M303_H', ins=M303-0, T=30+273.15)
M303_H_P = units.BDOPump('M303_H_P', ins=M303_H-0)



# Cofermentation
R302 = units.CoFermentation('R302', 
                                ins=('', M303_H_P-0, CSL),
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



M401 = bst.units.Mixer('M401', ins=(S401-1, '', '', '', ''))

def adjust_M401_ethanol_and_DPHP():
    feed_DPHP, feed_ethanol,  recycle_ethanol, recycle_DPHP = M401.ins[1:5]
    # print(S401.outs[0].F_mass)
    # print(recycle_ethanol.imass['Ethanol'])
    M401.ins[2].imass['Ethanol'] = max(0, S401.outs[1].F_mass * 0.24*(1.25*24/15.8) - recycle_ethanol.imass['Ethanol'])
    M401.ins[1].imass['Dipotassium hydrogen phosphate'] = max(0, S401.outs[1].F_mass * 0.25*(1.25*25/16.5) - recycle_DPHP.imass['Dipotassium hydrogen phosphate'])
    # print(M401.ins[2].imass['Ethanol'])
    M401._run()


M401.specification = adjust_M401_ethanol_and_DPHP
M401_P = units.BDOPump('M401_P', ins=M401-0, outs='mixed_stream')

# k_23bdo = 28.34
# k_glucose = 0.046
# k_etoh = 1
# k_h2o = 0

def adjust_S402_split():
    feed = S402.ins[0]
    IDs = ('2,3-Butanediol', 'Glucose', 'Ethanol', 'Water', 'Acetoin', 'Xylose')
    Ks = np.array([28.34, 0.046, 10000, 0, 28.34, 0.046])
    zs = feed.get_normalized_mass(IDs)
    L = tmo.equilibrium.binary_phase_fraction.solve_phase_fraction(zs, Ks, 0.9)
    # print (zs)
    # print(zs[0])
    # print(L)
    # (1/(L*K + 1 - L))*zs = x2
    x2 = ((1)/(L*Ks + 1 - L))*zs
    x1 = Ks * x2
    # print(x1)
    isplit = S402.isplit
    isplit['Water'] = 0.01
    isplit['Ethanol'] = 0.99
    isplit['Dipotassium hydrogen phosphate'] = 0.0001
    
    F_mass_eq = feed.imass[IDs].sum()
    L_mass = L * F_mass_eq
    isplit['2,3-Butanediol'] = L_mass * x1[0] / feed.imass['2,3-Butanediol']
    isplit['Acetoin'] = L_mass * x1[4] / feed.imass['Acetoin']
    isplit['Glucose'] = L_mass * x1[1] / feed.imass['Glucose']
    isplit['Xylose'] = L_mass * x1[5] / feed.imass['Xylose']
    S402._run()
                                    
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
    
    
    

split = [0.001 for i in range(len(BDO_chemicals))]
# split[BDO_chemicals.index('Dipotassium hydrogen phosphate')] = 0

S402 = bst.units.LiquidsSplitSettler('S402', ins = M401_P-0, split=split)
S402.specification = adjust_S402_split


# Separate out the majority of water,
# no need to include agitator thus using biosteam Flash
# F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'),
#                                     # LHK=('AceticAcid', '2,3-Butanediol'),
#                                     # is_divided=True,
#                                     # product_specification_format='Recovery',
#                                     # Lr=0.8, Hr=0.8, k=1.2,
#                                     T = 379, P = 101325,
#                                     vessel_material = 'Stainless steel 316')


F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'),
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
S403 = units.DPHP_Separation('S403', ins = F401_P-0, outs = ('separated_DPHP', 'impurities'))
S403_DPHP_recovery = 0.90 # w/w
S403_power_utility = 136850 # kW

def adjust_S403_streams():
    S403.outs[0].imass['Dipotassium hydrogen phosphate'] = S403_DPHP_recovery*S403.ins[0].imass['Dipotassium hydrogen phosphate']
    S403.outs[1] = S403.ins[0].copy()
    S403.outs[1].imass['Dipotassium hydrogen phosphate'] -= S403.outs[0].imass['Dipotassium hydrogen phosphate']
    S403.power_utility(S403_power_utility)
    
    S403._run()
    
S403.specification = adjust_S403_streams

# DPHP recycle
S403-0-4-M401

# D401 = bst.units.BinaryDistillation('D401', ins=F401_H-0, outs=('D401_g', 'D401_l'),
#                                     LHK=('Ethanol', 'Water'),
#                                     is_divided=True,
#                                     product_specification_format='Recovery',
#                                     Lr=0.99, Hr=0.99, k=1.2,
#                                     vessel_material = 'Stainless steel 316')


# D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0, rigorous=True)
# D401_P = units.BDOPump('D401_P', ins=D401-1)


# Separate out ethanol

# H_D401 = bst.units.HXutility('H_D401', ins=S402-0, V=0.75, rigorous=True)
D401 = bst.units.ShortcutColumn('D401', ins=S402-0,
                                    outs=('D401_g', 'D401_l'),
                                    LHK=('Ethanol', 'BDO'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.999, Hr=0.999, k=1.2,
                                    vessel_material = 'Stainless steel 316')
# def adjust_D401_V():
#     D401_ins_0 = D401.ins[0]
#     D401.V = D401_ins_0.imol['Ethanol']/ (D401_ins_0.F_mol*.95)
#     D401._run()
    
# D401 = bst.units.Flash('D401', ins=S402-0,
#                                     outs=('D401_g', 'D401_l'),
#                                     V=0.79, P = 101325,
#                                     vessel_material = 'Stainless steel 316')
# D401.process_specification = adjust_D401_V

D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0, rigorous=True)
D401_P = units.BDOPump('D401_P', ins=D401-1)





# M402 = bst.units.Mixer('M402', ins=(D401_H-0))
M402_P = units.BDOPump('M402_P', ins=D401_H-0, outs='ethanol_recycle')

# Ethanol recycle
M402_P-0-3-M401

# Separate out Acetoin

def adjust_D402_streams():
    feed, top = D402.ins[0], D402.outs[0]
    
    IDs = ('HMF', 'AceticAcid', 'Furfural')
    flows = feed.imol[IDs]
    
    feed.imol[IDs] = 0
    D402._run()
    top.imol[IDs] = flows
    feed.imol[IDs] = flows
    
    
D402 = bst.units.ShortcutColumn('D402', ins=D401_P-0,
                                    outs=('D402_g', 'D402_l'),
                                    LHK=('Acetoin', 'BDO'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.9995, Hr=0.9995, k=1.2,
                                    vessel_material = 'Stainless steel 316')
D402.specification = adjust_D402_streams
D402_H = bst.units.HXutility('D402_H', ins=D402-0, V=0, rigorous=True)
D402_P = units.BDOPump('D402_P', ins=D402-1)


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
M501 = bst.units.Mixer('M501', ins=(F401_H-0, F301-1, M205_P-0))

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

# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
                              outs=('aerobic_vent', 'aerobic_treated_water'),
                              reactants=soluble_organics,
                              ratio=U101.feedstock_flow_rate/2205)

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
ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
# lime_fresh = Stream('lime_fresh', price=price['Lime'])

# S401_out1_F_mass = S401.outs[1].F_mass

# if not (S401_out1_F_mass == 0):
#     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M401.ins[3].imass['Ethanol']
#     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M401.ins[3].imass['Dipotassium hydrogen phosphate']
    
# else:
ethanol_fresh = Stream('ethanol_fresh', Ethanol = feedstock_dry_mass*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
DPHP_fresh = Stream('DPHP_fresh', DPHP = feedstock_dry_mass*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

# Final product, not pure acid (which should be the case in reality)
BDO = Stream('BDO', units='kg/hr', price=price.get('BDO', 1.0)) # TODO: Check this

# Acetoin product
Acetoin = Stream('Acetoin', units='kg/hr', price=price['Acetoin'])
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
CIP_chems_in = Stream('CIP_chems_in', Water=145*U101.feedstock_flow_rate/2205, 
                      units='kg/hr')
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_like(CIP_chems_in)

# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including enzyme production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                      N2=0.79*1372608*U101.feedstock_flow_rate/2205,
                      O2=0.21*1372608*U101.feedstock_flow_rate/2205)

# 8021 based on stream 713 in Humbird et al.
fire_water_in = Stream('fire_water_in', 
                       Water=8021*U101.feedstock_flow_rate/2205, units='kg/hr')

# =============================================================================
# Facilities units
# =============================================================================

T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh)
T601.line = 'Sulfuric acid storage tank'
S601 = bst.units.ReversedSplitter('S601', ins=T601-0, 
                                  outs=(pretreatment_sulfuric_acid, 
                                        ''))

T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia)
T602.line = 'Ammonia storage tank'

T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
T603.line = 'CSL storage tank'

# DPHP storage
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
T605_P-0-2-M401

# 7-day storage time, similar to ethanol's in Humbird et al.
T606 = units.BDOStorageTank('T606', ins=D402_P-0, tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Stainless steel')



T606.line = 'BDOStorageTank'
T606_P = units.BDOPump('T606_P', ins=T606-0, outs=BDO)

# 7-day storage time, similar to ethanol's in Humbird et al.
T607 = units.BDOStorageTank('T607', ins=D402_H-0, tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Stainless steel')



T607.line = 'AcetoinStorageTank'
T607_P = units.BDOPump('T607_P', ins=T607-0, outs=Acetoin)


CIP = facilities.OrganicAcidsCIP('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = facilities.OrganicAcidsADP('ADP', ins=plant_air_in, outs='plant_air_out')


FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')

# M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
BT = facilities.OrganicAcidsBT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chems,
                                          baghouse_bag, natural_gas,
                                          'BT_makeup_water'),
                                B_eff=0.8, TG_eff=0.85,
                                combustibles=combustibles,
                                side_streams_to_heat=(pretreatment_feedstock_water,
                                                      pretreatment_acid_water,
                                                      pretreatment_steam),
                                outs=('gas_emission', ash, 'boiler_blowdown_water'))


# BT = bst.BDunits.BoilerTurbogenerator('BT',
#                                    ins=(M505-0, R501-0, 'boiler_makeup_water', 'natural_gas', FGD_lime, boiler_chems),
#                                    boiler_efficiency=0.80,
#                                    turbogenerator_efficiency=0.85)

# Blowdown is discharged
CT = facilities.OrganicAcidsCT('CT', 
                                ins=('return_cooling_water',
                                    'CT_makeup_water',
                                    cooling_tower_chems),
                                outs=('process_cooling_water',
                                      'cooling_tower_blowdown'))

# All water used in the system, here only consider water usage,
# if heating needed, then heeating duty required is considered in BT
process_water_streams = (pretreatment_feedstock_water, pretreatment_acid_water,
                         pretreatment_steam, pretreatment_ammonia_water,
                         enzyme_water,
                         aerobic_caustic, 
                         CIP.ins[-1], BT.ins[-1], CT.ins[-1])

PWC = facilities.OrganicAcidsPWC('PWC', ins=system_makeup_water, 
                                 process_water_streams=process_water_streams,
                                 outs='process_water')

# Heat exchange network
HXN = bst.facilities.HeatExchangerNetwork('HXN')


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
#                   D404, D404_H, D404_P, # separate out ethanol for recylcing
#                   F402, F402_H, F402_P], # separate out volatiles
#                 recycle=F402_H-0), # recycle ester
#           ],
#           recycle=D404_H-0), # recycle ethanol
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

BDO_sys = bst.main_flowsheet.create_system(
    'BDO_sys', feeds=[i for i in bst.main_flowsheet.stream
                            if i.sink and not i.source])

BT_sys = System('BT_sys', path=(BT,))


# =============================================================================
# TEA
# =============================================================================

BDO_no_BT_tea = BDOTEA(
        system=BDO_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(U101, WWT_cost,
                    T601, T602, T603, T604, T604_P, T605, T605_P, T606, T606_P,
                    CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=2.5e6*_labor_2007to2016*U101.feedstock_flow_rate/2205,
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


# =============================================================================
# Simulate system and get results
# =============================================================================

def get_BDO_MPSP():
    BDO_sys.simulate()
    
    for i in range(3):
        BDO.price = BDO_tea.solve_price(BDO, BDO_no_BT_tea)
    return BDO.price

get_BDO_MPSP()

# R301 = F('R301') # Fermentor
yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = F301,
    mixer = M303,
    reactor=R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('BDO',),
    yield_=0.909,
    titer=100,
    productivity=18.5,
    path = (M303_H, M303_H_P),
    xylose_utilization_fraction = 0.80)

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
    MPSP = BDO.price = BDO_tea.solve_price(BDO, BDO_no_BT_tea)
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
                        F401, F401_H, F401_P,
                        D401, D401_H, D401_P, S403,
                        M402_P, S403,
                        D402, D402_H, D402_P,
                        M501,
                        T606, T606_P, T607, T607_P)
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









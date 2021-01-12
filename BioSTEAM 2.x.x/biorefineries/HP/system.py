#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Sun Aug 23 12:11:15 2020

@author: sarangbhagwat

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 3-Hydroxypropionic acid instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

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


"""


# %% Setup

import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import main_flowsheet as F
from biosteam.process_tools import BoundedNumericalSpecification
from biosteam import System
from thermosteam import Stream
from biorefineries.HP import units, facilities
from biorefineries.HP._process_specification import ProcessSpecification
from biorefineries.HP.process_settings import price, CFs
from biorefineries.HP.utils import find_split, splits_df, baseline_feedflow
from biorefineries.HP.chemicals_data import HP_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.HP.tea import HPTEA
from biosteam.process_tools import UnitGroup

import matplotlib.pyplot as plt
# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2
bst.speed_up()
flowsheet = bst.Flowsheet('HP')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(HP_chemicals)

System.default_maxiter = 1500
# System.default_converge_method = 'fixed-point'
System.default_converge_method = 'aitken'
System.default_molar_tolerance = 0.5

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
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid', units='kg/hr')
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
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid,
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
fermentation_lime = Stream('fermentation_lime', units='kg/hr')

# For diluting concentrated, inhibitor-reduced hydrolysate
dilution_water = Stream('dilution_water', units='kg/hr')


# =============================================================================
# Conversion units
# =============================================================================

# Cool hydrolysate down to fermentation temperature at 50°C
H301 = bst.units.HXutility('H301', ins=P201-0, T=50+273.15)

# Mix enzyme with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, enzyme, enzyme_water))



# Saccharification and Cofermentation
# R301 = units.SaccharificationAndCoFermentation('R301', 
#                                                ins=(M302-0, CSL),
#                                                outs=('fermentation_effluent', 
#                                                      'sidedraw'))

# Saccharification
R301 = units.Saccharification('R301', 
                                ins=M301-0,
                                outs='saccharification_effluent')

# M303 = bst.units.Mixer('M303', ins=(R301-0, ''))
# M303_P = units.HPPump('M303_P', ins=M303-0)
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



# F301 = bst.units.MultiEffectEvaporator('F301', ins=M304-0, outs=('F301_l', 'F301_g'),
#                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.9695)
# F301_H = bst.units.HXutility('F301_H', ins=F301-0, T=30+273.15)
# F301_H_P = units.HPPump('F301_H_P', ins=F301_H-0)

F301 = bst.units.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.813)
                                        # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)
F301.V = 0.797 #for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)


F301_P = units.HPPump('F301_P', ins=F301-1)
# F301_H = bst.units.HXutility('F301_H', ins=F301-0, V = 0.)

    
M304_H_P = units.HPPump('M304_H_P', ins=F301-0)
M304 = bst.units.Mixer('M304', ins=(M304_H_P-0, dilution_water))
# M304 = bst.units.Mixer('M304', ins=(S301-1, dilution_water, ''))

M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15)

# Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
M302 = bst.units.Mixer('M302', ins=(M304_H-0, ''))


# inoculum_ratio = 0.07

S302 = bst.Splitter('S302', ins=M302-0,
                    outs = ('to_cofermentation', 'to_seedtrain'),
                    split = 0.07) # split = inoculum ratio

# Cofermentationv
# R302 = units.CoFermentation_original('R302', 
#                                 ins=(M304_H_P-0, '', CSL),
#                                 outs=('fermentation_effluent', 'CO2'))


R302 = units.CoFermentation('R302', 
                                ins=(S302-1, '', CSL, fermentation_lime),
                                outs=('fermentation_effluent', 'CO2'),
                                vessel_material='Stainless steel 316',
                                neutralization=True)


# ferm_ratio is the ratio of conversion relative to the fermenter
R303 = units.SeedTrain('R303', ins=S302-0, outs=('seed', 'CO2'), ferm_ratio=0.9)

T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-M302)


# %% 

# =============================================================================
# Separation streams
# =============================================================================

separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')

gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])

separation_decanol = Stream('separation_sulfuric_acid', units='kg/hr')
separation_TOA = Stream('separation_sulfuric_acid', units='kg/hr')
separation_AQ336 = Stream('separation_sulfuric_acid', units='kg/hr')


# separation_water = Stream('water', units='kg/hr')
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

R401 = units.AcidulationReactor('R401', ins = (S401-1, separation_sulfuric_acid),
                                outs = ('acidulated_broth'),
                                vessel_material='Stainless steel 316')
R401_P = bst.units.Pump('R401_P', ins=R401-0)

S402_index = S401_index + ['Gypsum']
S402_gypsum_split = S401_cell_mass_split + [0.995]
S402_filtrate_split = S401_filtrate_split + [0.005]
S402 = units.GypsumFilter('S402', ins=R401_P-0,
                          moisture_content=0.2,
                          split=find_split(S402_index,
                                           S402_gypsum_split,
                                           S402_filtrate_split,
                                           chemical_groups),
                          outs=(gypsum, ''))
def S402_spec():
    if S402.ins[0].imol['CaSO4']>0:
        S402._run()
    else:
        S402.outs[0].mol[:] = 0
        S402.outs[1].mol = S402.ins[0].mol
    
    S402.outs[1].imol['MethylHP'] = S402.outs[1].imol['HP']
    S402.outs[1].imol['HP'] = 0.
    
S402.specification = S402_spec

M401 = bst.units.Mixer('M401', ins=(separation_decanol, ''))


Kds = dict(IDs=('MethylHP',),
           K=np.array([1./7.269950288809821]), # 7. in g/L per g/L
           raffinate_chemicals = ('Water',),
           extract_chemicals = ('Decanol'))
S404 = bst.units.MultiStageMixerSettlers('S404', ins = (S402-1, M401-0),
                                     outs = ('raffinate', 'extract'),
                                     N_stages = 25, partition_data = Kds)
                                     

def adjust_S404_streams():
    feed_decanol, solvent_recycle = M401.ins
    process_stream = S404.ins[0]
    
    existing_decanol = solvent_recycle.ivol['Decanol'] + process_stream.ivol['Decanol']
    
    reqd_decanol = process_stream.F_vol * 0.5 * 0.8 * 2.5
    
    feed_decanol.ivol['Decanol'] = max(0, reqd_decanol - existing_decanol)
    
    M401._run()
    S404._run()

S404.specification = adjust_S404_streams


def partition_coefficients_mass_to_mol(IDs, ys, xs, chemicals):
    IDs = tuple(IDs)
    ys = np.asarray(ys)
    xs = np.asarray(xs)
    MWs = np.array([i.MW for i in chemicals[IDs]])
    ys_mol = ys / MWs
    ys_mol /= ys_mol.sum()
    xs_mol = xs / MWs
    xs_mol /= xs_mol.sum()
    return ys_mol / xs_mol

def composition_from_partition_data(zs, Ks, phi):
    return zs * Ks / (phi * Ks + (1 - phi))


# IDs = ('Water', 'Decanol', '3-Hydroxypropionic acid')
# chemicals = tmo.Chemicals(IDs)
# chemicals['3-Hydroxypropionic acid'].copy_models_from(tmo.Chemical('Lactic acid'), ['V'])
# chemicals.compile(skip_checks=True)
# tmo.settings.set_thermo(chemicals)
# z = 0.5
# K = 7
# phi = 0.5
# total = 1000
# y = composition_from_partition_data(z, K, phi) # For solute
# x = z - y * phi
# ys = np.array([0, 1000 - y, y]) # Extract
# xs = np.array([1000 - x, 1e-6, x]) # Raffinate
# extract = tmo.Stream(flow=ys)
# raffinate = tmo.Stream(flow=xs)
# ys /= extract.get_property('rho', 'g / L')
# xs /= raffinate.get_property('rho', 'g / L')
# Ks = partition_coefficients_mass_to_mol(IDs, ys, xs,
#                                         chemicals=chemicals)
# print(Ks)
    
# sugars recycle (assumes solvent is non-toxic to fermentation microbes)
# S404-0-1-R302
# Split_S404_raffinate = bst.units.Splitter('Split_S404_raffinate', ins = S404-0, split = 0.85)
# Split_S404_raffinate-0-1-R302

# D402 = bst.units.BinaryDistillation('D402', ins=S404-1, outs=('D402_g', 'D402_l'),
#                                     LHK=('HP', 'AQ336'),
#                                     is_divided=True,
#                                     product_specification_format='Recovery',
#                                     Lr=0.9999, Hr=0.9999, k=1.2,
#                                     P = 101325/100,
#                                     vessel_material = 'Stainless steel 316')

# D402-1-3-M401


# M402 = bst.units.Mixer('M402', ins=('separation_water', ''))

# Kds = dict(IDs=('HP',),
#            # K=np.array([7., 100000., 100000., 40132., 0.113]),
#            K = np.array([122.9155075014064]), # 7. in g/L per g/L
#            raffinate_chemicals = ('TOA', 'AQ336', 'Decanol'),
#            extract_chemicals = ('Water',))
# S405 = bst.units.MultiStageMixerSettlers('S405', ins = (S404-1, M402-0),
#                                      outs = ('raffinate', 'extract'),
#                                      N_stages = 25, partition_data = Kds)
          
# S405-0-3-M401
                
# def split_objective_function(split):
#     S406.split = split
#     S406._run()
#     M402._run()
#     S405._run()
#     feed_water, water_recycle = M402.ins
#     process_stream = S405.ins[0]
#     existing_water = water_recycle.ivol['Water'] + process_stream.ivol['Water']
#     reqd_water = process_stream.F_vol * 0.75
#     return existing_water - reqd_water

# def adjust_S405_streams():
#     available_recycle_water, = S406.ins
#     feed_water, water_recycle = M402.ins
#     process_stream = S405.ins[0]
#     available_recycle_water_flow = available_recycle_water.F_vol
#     existing_water = available_recycle_water_flow + process_stream.ivol['Water']
#     reqd_water = process_stream.F_vol * 0.55 * 10
#     waste_water_flow = existing_water - reqd_water
#     feed_water.ivol['Water'] = max(0, -waste_water_flow)
#     if waste_water_flow > 0:
#         S406.split[:] = 1. - (waste_water_flow / available_recycle_water_flow)
#         S406._run()
#     M402._run()
#     S405._run()
#     S405_raffinate = S405.outs[0]
#     S405_raffinate.imol['Water'] = 0
#     S405.outs[1].imol['Water'] = S405.ins[0].imol['Water'] + S405.ins[1].imol['Water']
#     if S405_raffinate.imol['HP'] < 0:
#         S405_raffinate.imol['HP'] = 0

# # def partition_coefficients_mass_to_mol(IDs):
    
# S405.specification = adjust_S405_streams
D401 = bst.units.ShortcutColumn('D401', ins=S404-1, outs=('D401_g', 'D401_l'),
                                    LHK=('MethylHP', 'Decanol'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.995, Hr=0.995, k=1.2,
                                    vessel_material = 'Stainless steel 316')

# def D402_remove_heat_utilities():
#     D402._run()
#     D402.heat_utilities = ()
# D402.specification = D402_remove_heat_utilities

D401_P = units.HPPump('D401_P', ins=D401-1)
D401_P-0-1-M401

#!!! TODO: Make rigorous=True after implementing Esterification and Hydrolysis
D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0., rigorous=False)
def D401_H_spec():
    D401_H._run()
    outstream = D401_H.outs[0]
    outstream.imol['HP'] = outstream.imol['MethylHP']
    outstream.imol['MethylHP'] = 0.
    
    outstream.imol['Water'] = 7.*outstream.imol['HP']
#     instream.vle(T=instream.T, P=instream.P)
#     D401_H._run()
#     # D401_H.outs[0].phases = ('g',)
D401_H.specification = D401_H_spec

R402 = units.DehydrationReactor('R402', ins = (D401_H-0),
                                outs = ('dilute_acryclic_acid'),
                                tau = 57.34/1.5, # Dishisha et al.
                                T = 230 + 273.15,
                                vessel_material='Stainless steel 316')

# def R402_specification():

#     R402._run()

# R402.specification = R402_specification


# R402_H =
# Separate out the majority of water,
# no need to include agitator thus using biosteam Flash
# D401 = bst.units.Flash('D401', ins=S402-1, outs=('D401_g', 'D401_l'),
#                                     # LHK=('AceticAcid', '2,3-Butanediol'),
#                                     # is_divided=True,
#                                     # product_specification_format='Recovery',
#                                     # Lr=0.8, Hr=0.8, k=1.2,
#                                     T = 379, P = 101325,
#                                     vessel_material = 'Stainless steel 316')


# D401 = bst.units.Flash('D401', ins=R401-0, outs=('D401_g', 'D401_l'),
#                                     T = 375, P = 101325,
#                                     vessel_material = 'Stainless steel 316')

# D401 = bst.units.BinaryDistillation('D401', ins=R402-0, outs=('D401_g', 'D401_l'),
#                                     LHK=('Water', 'AcrylicAcid'),
#                                     is_divided=True,
#                                     product_specification_format='Recovery',
#                                     Lr=0.99, Hr=0.99, k=1.2,
#                                     vessel_material = 'Stainless steel 316')

# H401 = bst.units.Flash('F401', ins=R402-0, outs=('F401_l', 'F401_g'),
#                                     P = 101325, V = 0., vessel_material='Stainless steel 316')

# def F401_spec():
#     F401_instream = F401.ins[0]
#     F401.V = F401_instream.imol['Water']/F401_instream.F_mol
#     F401._run()

# F401.specification = F401_spec
# # # Condense waste vapor for recycling
# F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)



# def S406.specification

# F401_P = units.HPPump('F401_P', ins=F401-1)
# D401_P-0-3-M401 # solvent recycle

D402 = bst.units.ShortcutColumn('D402', ins=R402-0, outs=('D402_g', 'D402_l'),
                                    LHK=('Water', 'AcrylicAcid'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.99, k=1.2,
                                    vessel_material = 'Stainless steel 316')

# def D402_remove_heat_utilities():
#     D402._run()
#     D402.heat_utilities = ()
# D402.specification = D402_remove_heat_utilities

D402_P = units.HPPump('D402_P', ins=D402-1)
D402_H = bst.units.HXutility('D402_H', ins=D402-0, T = 308.15, rigorous=True)

# S406 = bst.units.Splitter('S406', ins = D402_H-0, outs = ('recycled_water', 'waste_water'), split = 0.95)


# S406-0-1-M402
# def D402_spec():
#     try:
#         D402._run()
        
#     except:
#         count = 0
#         feasible = False
#         while not feasible:
#             print('Tried ' + str(i))
#             try:
#                 D402.Lr-=0.01
#                 D402._run()
#                 feasible = True
#             except:
#                 feasible = False
#             count+=1
#             if count>9:
#                 break

# def D402_spec():
#     D402._run()
#     D402.outs[0].imol['AQ336']=0
# D402.specification = D402_spec


# D402_H = bst.units.HXutility('D402_H', ins=D402-0, V=0, rigorous=True)
# D402_P = units.HPPump('D402_P', ins=D402-1)





# # # Condense waste vapor for recycling
# F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)
# F401_P = units.HPPump('F401_P', ins=F401-1)


# S403 = bst.units.Splitter('S403', ins=F401_P-0, outs=('to_fermentor', 
#                                                       'to_M501'),
#                                                       split=0.96)

# S403-0-1-R302




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
M501 = bst.units.Mixer('M501', ins=(F301_P-0, S401-0, S404-0)) # without sugars recycle
# M501 = bst.units.Mixer('M501', ins=(F301_P-0, S401-0, S406-1)) # with sugars recycle

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
M505 = bst.units.Mixer('M505', ins=(S503-1, S301-0, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')


# %% 

# =============================================================================
# Facilities streams
# =============================================================================

sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
sulfuric_acid_fresh2 = Stream('sulfuric_acid_fresh2',  price=price['Sulfuric acid'])
# TCP_fresh = Stream('TCP_fresh',  price=price['TCP'])
ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
lime_fresh = Stream('lime_fresh', price=price['Lime'])

decanol_fresh = Stream('decanol_fresh', price=price['Decanol'])
TOA_fresh = Stream('TOA_fresh', price=price['TOA'])
AQ336_fresh = Stream('AQ336_fresh', price=price['AQ336'])

# S401_out1_F_mass = S401.outs[1].F_mass

# if not (S401_out1_F_mass == 0):
#     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M401.ins[3].imass['Ethanol']
#     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M401.ins[3].imass['Dipotassium hydrogen phosphate']
    
# else:
# ethanol_fresh = Stream('ethanol_fresh', Ethanol = get_feedstock_dry_mass()*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
# DPHP_fresh = Stream('DPHP_fresh', DPHP = get_feedstock_dry_mass()*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

# HP stream
# HP = Stream('HP', units='kg/hr', price=price['HP'])
# AA product
AA = Stream('AcrylicAcid', units='kg/hr', price=price['AA'])
# Acetoin product
Acetoin = Stream('Acetoin', units='kg/hr', price=price['Acetoin'])
# Isobutyraldehyde product
IBA = Stream('IBA', units='kg/hr', price=price['IBA'])
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
                                     outs=pretreatment_sulfuric_acid)
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
# T604 = units.DPHPStorageTank('T604', ins=DPHP_fresh)
# T604.line = 'DPHP storage tank'
# T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0)

# 7-day storage time, similar to ethanol's in Humbird et al.
# T605 = bst.units.StorageTank('T605', ins=ethanol_fresh,
#                                      tau=7*24, V_wf=0.9,
#                                      vessel_type='Floating roof',
#                                      vessel_material='Carbon steel')
# T605.line = 'Ethanol storage tank'
# T605_P = units.HPPump('T605_P', ins=T605-0)



# # Connections to ATPE Mixer
# T604_P-0-1-M401
# T605_P-0-2-M401

# 7-day storage time, similar to ethanol's in Humbird et al.
T606 = units.HPStorageTank('T606', ins=D402_P-0, tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Stainless steel')



T606.line = 'AcrylicAcidStorageTank'
T606_P = units.HPPump('T606_P', ins=T606-0, outs=AA)

T607 = units.LimeStorageBin('T607', ins=lime_fresh, outs=fermentation_lime)
T607.line = 'Lime storage tank'


T608 = units.SulfuricAcidStorageTank('T608', ins = sulfuric_acid_fresh2, outs = separation_sulfuric_acid)
T608.line = 'Sulfuric acid storage tank'


T609 = bst.units.StorageTank('T609', ins = decanol_fresh, outs = separation_decanol)
T609.line = 'Decanol storage tank'

T610 = bst.units.StorageTank('T610', ins = TOA_fresh, outs = separation_TOA)
T610.line = 'TOA storage tank'

T611 = bst.units.StorageTank('T611', ins = AQ336_fresh, outs = separation_AQ336)
T611.line = 'AQ336 storage tank'


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
                         dilution_water,
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

# HP_sys = System('HP_sys',
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
#       D401, D401_H, D401_P, # separate water
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
HP_sys = bst.main_flowsheet.create_system(
    'HP_sys', feeds=[i for i in bst.main_flowsheet.stream
                            if i.sink and not i.source])

BT_sys = System('BT_sys', path=(BT,))

# %%
# =============================================================================
# TEA
# =============================================================================

#!!! Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legistration)
HP_no_BT_tea = HPTEA(
        system=HP_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(U101, WWT_cost,
                    T601, T602, T603, T606, T606_P,
                    CWP, CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

HP_no_BT_tea.units.remove(BT)

# # Removed because there is not double counting anyways.
# # Removes feeds/products of BT_sys from HP_sys to avoid double-counting
# for i in BT_sys.feeds:
#     HP_sys.feeds.remove(i)
# for i in BT_sys.products:
#     HP_sys.products.remove(i)

# Boiler turbogenerator potentially has different depreciation schedule
BT_tea = bst.TEA.like(BT_sys, HP_no_BT_tea)
BT_tea.labor_cost = 0

# Changed to MACRS 20 to be consistent with Humbird
BT_tea.depreciation = 'MACRS20'
BT_tea.OSBL_units = (BT,)

HP_tea = bst.CombinedTEA([HP_no_BT_tea, BT_tea], IRR=0.10)
HP_sys._TEA = HP_tea

# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================


# def get_HP_MPSP():
#     HP_sys.simulate()
    
#     for i in range(3):
#         HP.price = HP_tea.solve_price(HP, HP_no_BT_tea)
#     return HP.price

num_sims = 1
num_solve_tea = 3
def get_AA_MPSP():
    for i in range(num_sims):
        HP_sys.simulate()
    for i in range(num_solve_tea):
        AA.price = HP_tea.solve_price(AA, HP_no_BT_tea)
    return AA.price

get_AA_MPSP()

# R301 = F('R301') # Fermentor
# yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = F301,
    pump = M304_H_P,
    mixer = M304,
    heat_exchanger = M304_H,
    reactor= R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('HP','CalciumLactate'),
    spec_1=100,
    spec_2=0.909,
    spec_3=18.5,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = R401,
    byproduct_streams = [Acetoin, IBA],
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = R201)

path = (F301, R302)
@np.vectorize
def calculate_titer(V):
    F301.V = V
    for i in path: i._run()
    return spec._calculate_titer()

@np.vectorize   
def calculate_MPSP(V):
    F301.V = V
    HP_sys.simulate()
    MPSP = AA.price = HP_tea.solve_price(AA, HP_no_BT_tea)
    return MPSP

fermentation_broth = R302.outs[0]

get_titer = lambda: ((fermentation_broth.imass['HP'] + fermentation_broth.imol['CalciumLactate']*2*HP_chemicals.HP.MW)/fermentation_broth.F_vol)


def set_titer(titer):
    curr_titer = get_titer()
    M304.water_multiplier *= curr_titer/titer
    get_AA_MPSP()
    
# Unit groups


pretreatment = UnitGroup('Pretreatment', 
                                          units = [i for i in HP_sys.units if i.ID[1]=='2'])

saccharification_and_fermentation = UnitGroup('Saccharification and Fermentation', 
                                          units = [i for i in HP_sys.units if i.ID[1]=='3'])

separation = UnitGroup('Separation', 
                                          units = [i for i in HP_sys.units if i.ID[1]=='4'])

waste_treatment = UnitGroup('Waste treatment', 
                                          units = [i for i in HP_sys.units if i.ID[1]=='5'])

product_storage_and_pumping = UnitGroup('Product storage and pumping', 
                                          units = [i for i in HP_sys.units if i.ID[1]=='6'])

unit_groups = [pretreatment, saccharification_and_fermentation, separation, waste_treatment,
               product_storage_and_pumping]


titers_to_plot = np.linspace(10, 300, 30)


def plot_heating_duty_contributions_across_titers(titers):
    contributions = list(np.zeros(len(titers)))
    for i in range(len(titers)):
        titer = titers[i]
        set_titer(titer)
        contributions[i] = [ug.get_heating_duty()/AA.F_mass for ug in unit_groups]
    
    contributions = np.array(contributions)
    fig, ax = plt.subplots()
    for j in range(len(contributions[0])):
        ax.plot(titers, contributions[:,j], label = unit_groups[j].name)
    legend = ax.legend()
    plt.show()
    return contributions

def plot_electricity_contributions_across_titers(titers):
    contributions = list(np.zeros(len(titers)))
    for i in range(len(titers)):
        titer = titers[i]
        set_titer(titer)
        contributions[i] = [ug.get_electricity_consumption()/AA.F_mass for ug in unit_groups]
    
    contributions = np.array(contributions)
    fig, ax = plt.subplots()
    for j in range(len(contributions[0])):
        ax.plot(titers, contributions[:,j], label = unit_groups[j].name)
    legend = ax.legend()
    plt.show()
    return contributions

def plot_installed_cost_contributions_across_titers(titers):
    contributions = list(np.zeros(len(titers)))
    for i in range(len(titers)):
        titer = titers[i]
        set_titer(titer)
        contributions[i] = [ug.get_installed_cost()/AA.F_mass for ug in unit_groups]
    
    contributions = np.array(contributions)
    fig, ax = plt.subplots()
    for j in range(len(contributions[0])):

        ax.plot(titers, contributions[:,j], label = unit_groups[j].name)
    legend = ax.legend()
    plt.show()
    return contributions



# plot_heating_duty_contributions_across_titers(titers_to_plot)
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





HP_sub_sys = {
#     'feedstock_sys': (U101,),
#     'pretreatment_sys': (T201, M201, M202, M203, 
#                          R201, R201_H, T202, T203,
#                          F201, F201_H,
#                          M204, T204, T204_P,
#                          M205, M205_P),
#     'conversion_sys': (H301, M301, M302, R301, R302, T301),
    # 'separation_sys': (S401, M401, M401_P,
    #                     S402, 
    #                     D401, D401_H, D401_P,
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

# for unit in sum(HP_sub_sys.values(), ()):
#     if not unit in HP_sys.units:
#         print(f'{unit.ID} not in HP_sys.units')

# for unit in HP_sys.units:
#     if not unit in sum(HP_sub_sys.values(), ()):
#         print(f'{unit.ID} not in HP_sub_sys')




# %%

# =============================================================================
# Life cycle analysis (LCA), waste disposal emission not included
# =============================================================================

TEA_feeds = set([i for i in HP_sys.feeds if i.price]+ \
    [i for i in BT_sys.feeds if i.price])

TEA_products = set([i for i in HP_sys.products if i.price]+ \
    [i for i in BT_sys.products if i.price]+[AA])
# 100-year global warming potential (GWP) from material flows
LCA_streams = TEA_feeds.copy()
LCA_stream = Stream('LCA_stream', units='kg/hr')
    
def get_material_GWP():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
    # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
    return chemical_GWP.sum()/AA.F_mass

# GWP from combustion of non-biogenic carbons
get_non_bio_GWP = lambda: (natural_gas.get_atomic_flow('C'))* HP_chemicals.CO2.MW / AA.F_mass
                           # +ethanol_fresh.get_atomic_flow('C'))* HP_chemicals.CO2.MW / AA.F_mass
                               

# GWP from electricity
get_electricity_use = lambda: sum(i.power_utility.rate for i in HP_sys.units)
get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
    / AA.F_mass

# CO2 fixed in lactic acid product
get_fixed_GWP = lambda: \
    AA.get_atomic_flow('C')*HP_chemicals.CO2.MW/AA.F_mass

get_GWP = lambda: get_material_GWP()+get_non_bio_GWP()+get_electricity_GWP()

# Fossil energy consumption (FEC) from materials
def get_material_FEC():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
    # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
    return chemical_FEC.sum()/AA.F_mass

# FEC from electricity
get_electricity_FEC = lambda: \
    get_electricity_use()*CFs['FEC_CFs']['Electricity']/AA.F_mass

# Total FEC
get_FEC = lambda: get_material_FEC()+get_electricity_FEC()

get_SPED = lambda: BT.system_heating_demand*0.001/AA.F_mass
AA_LHV = 31.45 # MJ/kg AA

# %% Full analysis
def simulate_and_print():
    get_AA_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_AA_MPSP():.3f}/kg')
    print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg AA')
    print(f'Non-bio GWP is {get_non_bio_GWP():.3f} kg CO2-eq/kg AA')
    print(f'FEC is {get_FEC():.2f} MJ/kg AA or {get_FEC()/AA_LHV:.2f} MJ/MJ AA')
    print(f'SPED is {get_SPED():.2f} MJ/kg AA or {get_SPED()/AA_LHV:.2f} MJ/MJ AA')
    print('--------------------\n')

simulate_and_print()





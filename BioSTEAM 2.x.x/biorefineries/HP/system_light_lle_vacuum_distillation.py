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
from biosteam.exceptions import InfeasibleRegion
import matplotlib.pyplot as plt
import copy
from biorefineries.cornstover import CellulosicEthanolTEA
from biosteam import SystemFactory
# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

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

System.default_maxiter = 100
# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'
# System.default_converge_method = 'wegstein'
# System.default_molar_tolerance = 0.5
feedstock_ID = 'Corn stover'

System.default_converge_method = 'wegstein'
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.default_molar_tolerance = 0.1
System.strict_convergence = True # True => throw exception if system does not converge; false => continue with unconverged system


@SystemFactory(ID = 'HP_sys')
def create_HP_sys(ins, outs):
    
    process_groups = []
    # %% 
    
    # =============================================================================
    # Feedstock
    # =============================================================================
    
    feedstock = Stream('feedstock',
                        baseline_feedflow.copy(),
                        units='kg/hr',
                        price=price[feedstock_ID])
    
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock, outs='milled_feedstock')
    # U101 = units.FeedstockSizeReduction('U101', ins=feedstock, outs='milled_feedstock')
    
    # IFF handling costs/utilities are included in feedstock cost and thus not considered here
    U101.cost_items['System'].cost = 0
    U101.cost_items['System'].kW = 0
    
    feedstock_group = UnitGroup('feedstock_group', units=(U101,))
    process_groups.append(feedstock_group)
    
    # %% 
    
    # =============================================================================
    # Pretreatment streams
    # =============================================================================
    
    # For pretreatment, 93% purity
    pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid', units='kg/hr')
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
    # H_M202._cost = lambda: None
    # H_M202._design = lambda: None
    # H_M202.heat_utilities[0].heat_exchanger = None
    # P_M203 =  bst.units.Pump('P_M203', ins=water_M203, P=13.*101325.)
    
    # H_M203 = bst.units.HXutility('H_M203', ins=water_M203,
    #                                   outs='steam_M203',
    #                                   T=268.+273.15, rigorous=True)
    # def H_M203_specification():
    #     # water_M203.show()
    #     H_M202._specification()
    #     M202._run()
    #     M203._run()
    #     H_M203.ins[0].mol[:] = M203.ins[1].mol[:]
    #     H_M203.ins[0].phase = 'l'
    #     # H_M203.ins[0].T = water_M203.T
    #     H_M203._run()
    #     H_M203.outs[0].phase = 'g'
    # H_M203.specification = H_M203_specification
    
    # Prepare sulfuric acid
    get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
    T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid,
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
                                     V=0., rigorous=True)
    H201-0-3-M202
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
    
    
    pretreatment_group = UnitGroup('pretreatment_group', 
                                   units=(H_M201, H_M202, T201, M201, M202, M203,
                                          R201, T202, T203, F201, M204, H201,
                                          M205, T204, P201,))
    process_groups.append(pretreatment_group)
    
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
    H301 = bst.units.HXutility('H301', ins=P201-0, T=50+273.15, rigorous=True)
    
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
                                    outs=('saccharification_effluent',''))
    R301.tau_saccharification = 24. # !!! 24 or 84 h?
    
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'Acetate', 'AceticAcid'] = isplit[ID]
                           
        
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
    fix_split(S301.isplit, 'Glucose')
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
    
        
    M304_P = units.HPPump('M304_P', ins=F301-0)
    M304 = bst.units.Mixer('M304', ins=(M304_P-0, dilution_water))
    # M304 = bst.units.Mixer('M304', ins=(S301-1, dilution_water, ''))
    
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    # M302 = bst.units.Mixer('M302', ins=(M304_H-0, ''))
    
    
    # inoculum_ratio = 0.07
    
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.07) # split = inoculum ratio
    
    # Cofermentation
    # R302 = units.CoFermentation_original('R302', 
    #                                 ins=(M304_H_P-0, '', CSL),
    #                                 outs=('fermentation_effluent', 'CO2'))
    
    
    R302 = units.CoFermentation('R302', 
                                    ins=(S302-1, '', CSL, fermentation_lime),
                                    outs=('fermentation_effluent', 'CO2_fermentation'),
                                    vessel_material='Stainless steel 316',
                                    neutralization=True)
    
    def include_seed_CSL_in_cofermentation(): # note: effluent always has 0 CSL
        R302._run()
        R302.ins[2].F_mass*=1./(1-S302.split[0])
    R302.specification = include_seed_CSL_in_cofermentation
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = units.SeedTrain('R303', ins=S302-0, outs=('seed', 'CO2_seedtrain'), ferm_ratio=0.9)
    
    T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-R302)
    
    
    conversion_group = UnitGroup('conversion_group', 
                                   units=(H301, M301, R301, S301, F301, F301_P,
                                          M304_P, M304, M304_H, S302, R302,
                                          R303, T301))
    process_groups.append(conversion_group)
    
    # %% 
    
    # =============================================================================
    # Separation streams
    # =============================================================================
    
    separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
    
    gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
    
    separation_hexanol = Stream('separation_hexanol', units='kg/hr')
    
    makeup_TiO2_catalyst = Stream('makeup_TiO2_catalyst', units='kg/hr', price=price['TiO2'])
    
    # =============================================================================
    # Separation units
    # =============================================================================
    
             
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=(R302-0,''), outs=('cell_mass', ''),
                                moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    fix_split(S401.isplit, 'Glucose')
    # NOTE: if there is not enough moisture content, it's impossible to pump
    # the fluid into the centrifuge. In fact, the centrifuge would not be able
    # to separate anything.
    # def S401_moisture_capper():
    #     try:
    #         S401._run()
    #     except:
    #         S401.moisture_content *= 0.9
    #         S401_moisture_capper()
    #     S401.moisture_content = 0.4
    # S401.specification = S401_moisture_capper
    
    R401 = units.AcidulationReactor('R401', ins = (S401-1, separation_sulfuric_acid),
                                    outs = ('acidulated_broth'),
                                    vessel_material='Stainless steel 316',
                                    tau = 1.)
    
    R401_H = bst.units.HXutility('R401_H', ins = R401-0, T = 320, rigorous = False)
    R401_P = bst.units.Pump('R401_P', ins=R401_H-0)
    
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
        
        
    S402.specification = S402_spec
    
    
    M401 = bst.units.Mixer('M401', ins=(separation_hexanol,
                                        ''))
    
    M401_H = bst.units.HXutility('M401_H', ins = M401-0, T = 80. + 273.15, rigorous = False)
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=S402-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
    # target_water_x = 0.35
    target_HP_x = 0.10
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol['AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'])
    
    def F401_specification():
        instream = F401.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        ratio = get_x('HP', instream)/target_HP_x
        # no need to check for ratio>1 becasue our target_water_x is consistently lower than the max possible titer
        F401.V = 1. - ratio
        
        F401._run()
        # import pdb
        # pdb.set_trace()
    
    F401.specification = F401_specification
    
    
    def F401_no_run_cost():
        F401.heat_utilities = tuple()
        F401._installed_cost = 0.
    # F401._cost = F401_no_run_cost
        
    F401_H = bst.units.HXutility('F401_H', ins = F401-0, T = 80. + 273.15, rigorous = False)
    F401_P = bst.units.Pump('F401_P', ins=F401_H-0)

    
    Kds = dict(IDs=('HP', 'Water', 'Hexanol', 'AceticAcid'),
               # K=np.array([1./1.941747572815534, 3.606, 0.006]),
               K=np.array([1./1.9379484051844278, 3.690183610720956, 0.0060176892697821486, 1./0.4867537504125923]), # T = 80. + 273.15 K
               phi = 0.5)
    
    S404 = bst.units.MultiStageMixerSettlers('S404', ins = (F401_P-0, M401_H-0),
                                         outs = ('raffinate', 'extract'),
                                         N_stages = 15, partition_data = Kds,) 
                              
                              
    S404.vol_frac = 0.05
    
    
    tolerable_loss_fraction = 0.001
    
    def adjust_S404_Ks_streams_with_comments():
        S404.N_stages = 15 # reset
        S404._setup() # reset
        feed_hexanol, solvent_recycle = M401.ins
        process_stream = S404.ins[0]
        process_stream_F_mol = process_stream.F_mol
        existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
        # N_stages = S404.N_stages
        # S404.T = 340
        K_raffinate = S404.partition_data['K'][0]
        # K_extract = 1./S404.partition_data['K'][0]
        HP_recovery = 1-tolerable_loss_fraction
        reqd_hexanol = HP_recovery * K_raffinate * process_stream_F_mol
        
        feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        M401._run()
        Ks_new = update_Ks(S404)
        # print(Ks_new)
        
        if np.any(np.abs(S404.partition_data['K']/ Ks_new - 1) > 1e-1):
            S404.partition_data['K'] = Ks_new 
        K_raffinate = S404.partition_data['K'][0]
        # K_extract = 1./S404.partition_data['K'][0]
        HP_recovery = 1.-tolerable_loss_fraction
        # prev_reqd_hexanol = reqd_hexanol
        
        reqd_hexanol = HP_recovery * K_raffinate * process_stream_F_mol
        
        feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        M401._run()
        S404_run()
        
    
    def adjust_S404_Ks_streams():
        S404.N_stages = 15 # reset
        S404._setup() # reset
        feed_hexanol, solvent_recycle = M401.ins
        process_stream = S404.ins[0]
        process_stream_F_mol = process_stream.F_mol
        existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
    
        K_raffinate = S404.partition_data['K'][0]
        HP_recovery = 1-tolerable_loss_fraction
        reqd_hexanol = HP_recovery * K_raffinate * process_stream_F_mol
    
            
            
        feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        M401._run()
        Ks_new = update_Ks(S404)
    
        if np.any(np.abs(S404.partition_data['K']/ Ks_new - 1) > 1e-1):
            S404.partition_data['K'] = Ks_new 
        K_raffinate = S404.partition_data['K'][0]
        HP_recovery = 1.-tolerable_loss_fraction
    
        reqd_hexanol = HP_recovery * K_raffinate * process_stream_F_mol
    
        feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        M401._run()
        S404_run()
        
    
    # def adjust_S404_streams():
    #     S404.N_stages = 15 # reset
    #     S404._setup() # reset
    #     feed_hexanol, solvent_recycle = M401.ins
    #     process_stream = S404.ins[0]
    #     process_stream_F_mol = process_stream.F_mol
    #     existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
    
    #     K_raffinate = S404.partition_data['K'][0]
    
    #     HP_recovery = 1-tolerable_loss_fraction
    #     reqd_hexanol = HP_recovery * K_raffinate * process_stream_F_mol
    #     # S404.reqd_hexanol = reqd_hexanol # for access in S404_run
    #     # if existing_hexanol > reqd_hexanol:
    #     #     solvent_recycle.empty()
    #     #     existing_hexanol = 0
    #     feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
    
    #     M401._run()
    #     S404_run()

    def adjust_S404_streams():
        S404.N_stages = 15 # reset
        S404._setup() # reset
        feed_hexanol, solvent_recycle = M401.ins
        process_stream = S404.ins[0]
        process_stream_F_mol = process_stream.F_mol
        existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
    
        K_raffinate = S404.partition_data['K'][0]
    
        HP_recovery = 1-tolerable_loss_fraction
        reqd_hexanol =  HP_recovery * K_raffinate * process_stream_F_mol
        if existing_hexanol > reqd_hexanol:
            feed_hexanol.imol['Hexanol'] = 0.
            solvent_recycle.imol['Hexanol'] = reqd_hexanol
        else:
            feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        
        M401._run()
        M401_H._run()
        S404_run()
        
        if existing_hexanol > reqd_hexanol:
            feed_hexanol.imol['Hexanol'] = S404.outs[0].imol['Hexanol']
            
        
    def update_Ks(lle_unit, solute_indices = (0,), carrier_indices = (1,), solvent_indices = (2,)):
        IDs = lle_unit.partition_data['IDs']
        Ks = lle_unit.partition_data['K']
        solute_chemicals = tuple([IDs[index] for index in solute_indices])
        carrier_chemicals = tuple([IDs[index] for index in carrier_indices])
        solvent_chemicals = tuple([IDs[index] for index in solvent_indices])
        process_stream = lle_unit.ins[0]
        solvent_stream = lle_unit.ins[1]
        
        test_stream = bst.Stream('test_stream')
        test_stream_2 = bst.Stream('test_stream_2')
        test_stream_2.mix_from([process_stream, solvent_stream])
        test_stream.imol[solute_chemicals] = process_stream.imol[solute_chemicals]
        test_stream.imol[carrier_chemicals] = process_stream.imol[carrier_chemicals]
        test_stream.imol[solvent_chemicals] = solvent_stream.imol[solvent_chemicals]
        test_stream.lle(T=process_stream.T, top_chemical = 'Hexanol')
        # test_stream.show()
        Ks_new = (test_stream['L'].imol[IDs]/test_stream['L'].F_mol)/(test_stream['l'].imol[IDs]/test_stream['l'].F_mol)
        
        return Ks_new
    
    def S404_run():
        # mixer = S404.mixer
        # mixer._run()
        # mixed = mixer-0
        # raffinate, extract = S404.outs
        # raffinate.copy_like(mixed)
        # raffinate.imol['H2O', 'Hexanol', 'HP'] *= [0.812, 0.007, 0.001]
        # extract.mol = mixed.mol - raffinate.mol
        # extract.T = raffinate.T = mixed.T
        # extract.P = raffinate.P = mixed.P
        try:
            S404._run()
            
            if has_negative_flows(S404):
                raise InfeasibleRegion('negative flows')
        except:
            # import pdb
            # pdb.set_trace()
            S404.N_stages-=1
            if S404.N_stages == 0:
                S404.N_stages = 15 # reset
                S404._setup() # reset
                raise InfeasibleRegion('number of stages in %s'%(S404.ID))   
            else:
                S404._setup()
                print('\nReduced S404.N_stages to %s\n'%S404.N_stages)
            S404_run()
        
    
    def has_negative_flows(unit):
        for stream in unit.outs + unit.ins:
            if (stream.mol < 0.).any():
                return True
        return False
                
    S404.specification = adjust_S404_streams
    
    ideal_thermo = S404.thermo.ideal()
    
    # S404-0-1-R302 # with sugars recycle
    D401 = bst.units.BinaryDistillation('D401', ins=S404-1, outs=('D401_g', 'D401_l'),
                                        LHK=('Hexanol', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P = 101325/20,
                                        vessel_material = 'Stainless steel 316',
                                        partial_condenser = False,
                                        condenser_thermo = ideal_thermo,
                                        boiler_thermo = ideal_thermo,
                                        thermo=ideal_thermo)
    # D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0., rigorous=True)
    D401_H_P = units.HPPump('D401_H_P', ins=D401-0, P = 101325)
    D401_H_P-0-1-M401
    
    def get_concentration_gpL(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_vol
    
    def get_mass_percent(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_mass
    
    M402 = bst.units.Mixer('M402', ins=(D401-1,
                                        'dilution_water2'))
    
    def M402_objective_fn(Water_imol):
        M402.ins[1].imol['Water'] = Water_imol
        M402._run()
        # return get_concentration_gpL('HP', M402.outs[0]) - 600. # predicted "solubility" of 645 g/L at STP https://hmdb.ca/metabolites/HMDB0000700
        # return get_concentration_gpL('HP', M402.outs[0]) - 270.1 # "Ssolubility "at 25 C # https://www.chemicalbook.com/ChemicalProductProperty_EN_CB6711580.htm
        # return get_mass_percent('HP', M402.outs[0]) - .15 # Dehydration reaction paper
        # return get_mass_percent('HP', M402.outs[0]) - .35 # 30-35 wt% in https://patents.google.com/patent/WO2013192451A1/en
        return get_mass_percent('HP', M402.outs[0]) - .30 # 30wt% with 80% conversion in Dunn et al. 2015 LCA of Bioproducts in GREET
    
    def M402_adjust_water():
        flx.IQ_interpolation(M402_objective_fn, 0., 10000, maxiter=50, ytol=1e-2)
        
    M402.specification = M402_adjust_water
    
    D401_P = units.HPPump('D401_P', ins=M402-0, P=506625*5.4)
    
    R402 = units.DehydrationReactor('R402', ins = (D401_P-0, makeup_TiO2_catalyst, '', ''),
                                    outs = ('dilute_acryclic_acid', 'spent_TiO2_catalyst'),
                                    tau = 57.34/1.5, # Dishisha et al.
                                    T = 230. + 273.15,
                                    vessel_material='Stainless steel 316')
    R402.heat_utilities[0].heat_transfer_efficiency = 1. # this doesn't seem to be changing duty
    
    # spent_TiO2_catalyst is assumed to be sold at 0 $/kg
    
    
    R402_H = bst.units.HXutility('R402_H', ins=R402-0, T = 372.00, rigorous=True)
    
    
    D402 = bst.units.ShortcutColumn('D402', ins=R402_H-0, outs=('D402_g', 'D402_l'),
                                        LHK=('Water', 'AcrylicAcid'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P=101325,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    
    D402_t_H = bst.units.HXutility('D402_H', ins=D402-0, T = 330., rigorous=True)
    
    # recycling water makes system convergence fail
    # D402-0-3-R402
    
    
    D402_P = units.HPPump('D402_P', ins=D402-1)
    
    
    # D402_H = bst.units.HXutility('D402_H', ins=D402_P-0, T = 330., rigorous=True)
    
    
    D403 = bst.units.ShortcutColumn('D403', ins=D402_P-0, outs=('D403_g', 'D403_l'),
                                        LHK=('AcrylicAcid', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.05, P=101325./7.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    
    D403-1-2-R402
    
    D403_H = bst.units.HXutility('D403_H', ins=D403-0, T = 25.+273.15, rigorous=True)
    D403_P = units.HPPump('D403_P', ins=D403_H-0)
    
    separation_group = UnitGroup('separation_group', 
                                   units=(S401, R401, R401_H, R401_P, S402,
                                          F401, F401_P, M401, S404, D401,
                                          D401_H_P, D401_P, M402, 
                                          R402, R402_H, D402, D402_t_H, D402_P, D403, D403_P))
    process_groups.append(separation_group)
    
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
    M501 = bst.units.Mixer('M501', ins=(F301_P-0, D402_t_H-0, F401-1, S404-0)) # without sugars recycle
    # M501 = bst.units.Mixer('M501', ins=(F301_P-0, D402_H-0)) # with sugars recycle
    
    # This represents the total cost of wastewater treatment system
    WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)
    
    R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                    outs=('biogas', 'anaerobic_treated_water', 
                                          'anaerobic_sludge'),
                                    reactants=soluble_organics + ['Glycerol', 'HP', 'Hexanol', 'AcrylicAcid'],
                                    split=find_split(splits_df.index,
                                                     splits_df['stream_611'],
                                                     splits_df['stream_612'],
                                                     chemical_groups),
                                    T=35+273.15)
    fix_split(R501.isplit, 'Glucose')
    
    # !!! Make sure to state and explain this conservative assumption
    # Working explanation: In TRY analysis, we aren't looking at the implications of varying yield of byproducts on 
    # sugars, but solely that of 3-HP. 
    # Just show the implications of the extremes of the assumption to show it doesn't affect results much.
    
    R501.byproducts_combustion_rxns = ParallelRxn([
        Rxn('AceticAcid -> 3 CO2 + H2O + O2', 'AceticAcid', 1.-1e-6, correct_atomic_balance=True),
        Rxn('Glycerol -> 3 CO2 + H2O + O2', 'Glycerol', 1.-1e-6, correct_atomic_balance=True)])
    for i in R501.byproducts_combustion_rxns:
        i.istoichiometry['O2'] = 0.
    def R501_specification():
        R501.byproducts_combustion_rxns(R501.ins[0])
        R501._run()
    # R501.specification = R501_specification # Comment this out for anything other than TRY analysis
    
    get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
    
    # Mix recycled stream and wastewater after R501
    M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
    R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
                                  outs=('aerobic_vent', 'aerobic_treated_water'),
                                  reactants=soluble_organics + ['Glycerol', 'HP', 'Hexanol', 'AcrylicAcid'],
                                  ratio=get_flow_tpd()/2205)
    
    # Membrane bioreactor to split treated wastewater from R502
    S501 = bst.units.Splitter('S501', ins=R502-1, outs=('membrane_treated_water', 
                                                        'membrane_sludge'),
                              split=find_split(splits_df.index,
                                               splits_df['stream_624'],
                                               splits_df['stream_625'],
                                               chemical_groups))
    
    S501.line = 'Membrane bioreactor'
    fix_split(S501.isplit, 'Glucose')
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
    fix_split(S503.isplit, 'Glucose')
    S503.line = 'Sludge centrifuge'
    
    # Reverse osmosis to treat membrane separated water
    S504 = bst.units.Splitter('S504', ins=S501-0, outs=('discharged_water', 'waste_brine'),
                              split=find_split(splits_df.index,
                                               splits_df['stream_626'],
                                               splits_df['stream_627'],
                                               chemical_groups))
    S504.line = 'Reverse osmosis'
    
    # Mix solid wastes to boiler turbogeneration
    
    # Mention results with and without S401-0 in manuscript
    M505 = bst.units.Mixer('M505', ins=(S503-1, S301-0, S401-0), 
                            outs='wastes_to_boiler_turbogenerator')
    
    
    WWT_group = UnitGroup('WWT_group', 
                                   units=(M501, WWT_cost,R501, M502, R502,
                                          S501, S502, M503, M504, S503, S504,
                                          M505,))
    process_groups.append(WWT_group)
    
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
    
    hexanol_fresh = Stream('hexanol_fresh', price=price['Hexanol'])
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

    T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia_M205)
    T602.line = 'Ammonia storage tank'
    
    T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
    T603.line = 'CSL storage tank'
    

    
    T604 = units.LimeStorageBin('T604', ins=lime_fresh, outs=fermentation_lime)
    T604.line = 'Lime storage tank'
    
    
    T605 = units.SulfuricAcidStorageTank('T605', ins = sulfuric_acid_fresh2, outs = separation_sulfuric_acid)
    T605.line = 'Sulfuric acid storage tank'
    

    # 7-day storage time, similar to ethanol's in Humbird et al.
    T606 = units.HPStorageTank('T606', ins=D403_P-0, tau=7*24, V_wf=0.9,
                                         vessel_type='Floating roof',
                                         vessel_material='Stainless steel')
   
    T606.line = 'AcrylicAcidStorageTank'
    T606_P = units.HPPump('T606_P', ins=T606-0, outs=AA)
    
    
    T607 = bst.units.StorageTank('T607', ins = hexanol_fresh, outs = separation_hexanol)
    T607.line = 'Hexanol storage tank'
    
    # T610 = bst.units.StorageTank('T610', ins = TOA_fresh, outs = separation_TOA)
    # T610.line = 'TOA storage tank'
    
    # T611 = bst.units.StorageTank('T611', ins = AQ336_fresh, outs = separation_AQ336)
    # T611.line = 'AQ336 storage tank'
    
    
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_tpd()/2205)
    
    
    FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')
    
    #!!! M304_H uses chilled water, thus requiring CWP
    CWP = facilities.CWP('CWP', ins='return_chilled_water',
                          outs='process_chilled_water')
    
    # CWP = bst.facilities.ChilledWaterPackage('CWP')
    
    # M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
    # BT = facilities.BT('BT', ins=(M505-0, R501-0, 
    #                                           FGD_lime, boiler_chems,
    #                                           baghouse_bag, natural_gas,
    #                                           'BT_makeup_water'),
    #                                 B_eff=0.8, TG_eff=0.85,
    #                                 combustibles=combustibles,
    #                                 side_streams_to_heat=(water_M201, water_M202, steam_M203),
    #                                 outs=('gas_emission', ash, 'boiler_blowdown_water'))
    
    # BT = bst.facilities.Boiler('BT',
    #                                               ins=(M505-0,
    #                                                   R501-0, 
    #                                                   'boiler_makeup_water',
    #                                                   'natural_gas',
    #                                                   'lime',
    #                                                   'boilerchems'), 
    #                                               outs=('gas_emission', 'boiler_blowdown_water', ash, '', ''),)
    #                                               # turbogenerator_efficiency=0.85)
    
    BT = bst.facilities.BoilerTurbogenerator('BT',
                                                  ins=(M505-0,
                                                      R501-0, 
                                                      'boiler_makeup_water',
                                                      'natural_gas',
                                                      'lime',
                                                      'boilerchems'), 
                                                  outs=('gas_emission', 'boiler_blowdown_water', ash,),
                                                  turbogenerator_efficiency=0.85)
    
    # Blowdown is discharged
    CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                                  'CT_makeup_water'),
                       outs=('process_cooling_water', 'cooling_tower_blowdown'))
    
    # All water used in the system, here only consider water usage,
    # if heating needed, then heeating duty required is considered in BT
    process_water_streams = (water_M201, water_M202, water_M203, water_M205, 
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
    def HXN_no_run_cost():
        HXN.heat_utilities = tuple()
        HXN._installed_cost = 0.
        
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []
    
    HXN_group = UnitGroup('HXN_group', 
                                   units=(HXN,))
    process_groups.append(HXN_group)
    

    BT_group = UnitGroup('BT_group',
                                   units=(BT,))
    process_groups.append(BT_group)
    
    CT_group = UnitGroup('CT_group',
                                   units=(CT,))
    process_groups.append(CT_group)
    
    facilities_no_hu_group = UnitGroup('facilities_no_hu_group',
                                   units=(T601, T602, T603, T604, T605, 
                                          T606, T606_P, T607, PWC, ADP, CIP))
    process_groups.append(facilities_no_hu_group)

    globals().update({'process_groups': process_groups})
# %% System setup

# HP_sys = bst.main_flowsheet.create_system(
#     'HP_sys', feeds=feeds)
HP_sys = create_HP_sys()
HP_sys.subsystems[-1].relative_molar_tolerance = 0.005

u = flowsheet.unit
s = flowsheet.stream
feedstock = s.feedstock
AA = s.AcrylicAcid
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185


# feeds = [i for i in flowsheet.stream
#                             if i.sink and not i.source]
feeds = HP_sys.feeds

products = [AA] # Don't include gypsum since we want to add carbon impurities to GWP

emissions = [i for i in flowsheet.stream
                            if i.source and not i.sink and not i in products]
    
BT = flowsheet('BT')
BT_sys = System('BT_sys', path=(BT,))

          
# flowsheet('SYS2').molar_tolerance = 3
# flowsheet('SYS2').maxiter = 100

globals().update(flowsheet.to_dict())

process_groups_dict = {}
for i in range(len(process_groups)):
    group = process_groups[i]
    process_groups_dict[group.name] = group
# %%
# =============================================================================
# TEA
# =============================================================================

HP_no_BT_sys = bst.System('HP_no_BT_sys', path = HP_sys.path, facilities = tuple([i for i in HP_sys.facilities if not i.ID=='BT']))


#!!! Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)
# HP_no_BT_tea = HPTEA(
#         system=HP_no_BT_sys, IRR=0.10, duration=(2016, 2046),
#         depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
#         lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
#         startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
#         startup_VOCfrac=0.75, WC_over_FCI=0.05,
#         finance_interest=0.08, finance_years=10, finance_fraction=0.4,
#         # biosteam Splitters and Mixers have no cost, 
#         # cost of all wastewater treatment units are included in WWT_cost,
#         # BT is not included in this TEA
#         OSBL_units=(u.U101, u.WWT_cost,
#                     u.T601, u.T602, u.T603, u.T606, u.T606_P,
#                     u.CWP, u.CT, u.PWC, u.CIP, u.ADP, u.FWT),
#         warehouse=0.04, site_development=0.09, additional_piping=0.045,
#         proratable_costs=0.10, field_expenses=0.10, construction=0.20,
#         contingency=0.10, other_indirect_costs=0.10, 
#         labor_cost=3212962*get_flow_tpd()/2205,
#         labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

# # HP_no_BT_tea.units.remove(BT)

# # # Removed because there is not double counting anyways.
# # # Removes feeds/products of BT_sys from HP_sys to avoid double-counting
# # for i in BT_sys.feeds:
# #     HP_sys.feeds.remove(i)
# # for i in BT_sys.products:
# #     HP_sys.products.remove(i)

# # Boiler turbogenerator potentially has different depreciation schedule
# BT_tea = bst.TEA.like(BT_sys, HP_no_BT_tea)
# BT_tea.labor_cost = 0

# # Changed to MACRS 20 to be consistent with Humbird
# BT_tea.depreciation = 'MACRS20'
# BT_tea.OSBL_units = (BT,)

# HP_tea = bst.CombinedTEA([HP_no_BT_tea, BT_tea], IRR=0.10)
# HP_sys._TEA = HP_tea


HP_tea = CellulosicEthanolTEA(system=HP_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWT_cost,
                    u.T601, u.T602, u.T603, u.T606, u.T606_P,
                    u.CWP, u.CT, u.PWC, u.CIP, u.ADP, u.FWT, u.BT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)

HP_no_BT_tea = HP_tea
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
        # AA.price = HP_tea.solve_price(AA, HP_no_BT_tea)
        AA.price = HP_tea.solve_price(AA)
    return AA.price

get_AA_MPSP()

seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))
# R301 = F('R301') # Fermentor
# yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = u.F301,
    pump = u.M304_P,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = seed_train_system,
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('HP','CalciumLactate'),
    # spec_1=100,
    # spec_2=0.909,
    # spec_3=18.5,
    spec_1=0.49,
    spec_2=54.8,
    spec_3=0.76,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = u.R401,
    byproduct_streams = [],
    HXN = u.HXN,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301],
    pre_conversion_units = HP_sys.split(u.F301.ins[0])[0],
    baseline_titer = 54.8,
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = u.R201)

# Load baseline specifications

spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

# spec.load_yield(0.49)
# spec.load_titer(54.8)
# spec.load_productivity(0.76)

spec.load_specifications(spec_1 = 0.49, spec_2 = 54.8, spec_3 = 0.76)


path = (u.F301, u.R302)
@np.vectorize
def calculate_titer(V):
    u.F301.V = V
    for i in path: i._run()
    return spec._calculate_titer()

@np.vectorize   
def calculate_MPSP(V):
    u.F301.V = V
    HP_sys.simulate()
    MPSP = AA.price = HP_tea.solve_price(AA, HP_no_BT_tea)
    return MPSP

fermentation_broth = u.R302.outs[0]

get_titer = lambda: ((fermentation_broth.imass['HP'] + fermentation_broth.imol['CalciumLactate']*2*HP_chemicals.HP.MW)/fermentation_broth.F_vol)


def set_titer(titer):
    curr_titer = get_titer()
    u.M304.water_multiplier *= curr_titer/titer
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
# Life cycle analysis (LCA), waste emission C assumed to be entirely CO2
# =============================================================================

TEA_feeds = set([i for i in HP_sys.feeds if i.price]+ \
    [i for i in BT_sys.feeds if i.price])

TEA_products = set([i for i in HP_sys.products if i.price]+ \
    [i for i in BT_sys.products if i.price]+[AA])
# 100-year global warming potential (GWP) from material flows
LCA_streams = TEA_feeds.copy()
LCA_stream = Stream('LCA_stream', units='kg/hr')

feed_chem_IDs = [chem.ID for chem in HP_chemicals]

GWP_CF_stream = CFs['GWP_CF_stream']
FEC_CF_stream = CFs['FEC_CF_stream']


# Carbon balance
total_C_in = sum([feed.get_atomic_flow('C') for feed in feeds])
total_C_out = AA.get_atomic_flow('C') + sum([emission.get_atomic_flow('C') for emission in emissions])
C_bal_error = (total_C_out - total_C_in)/total_C_in

def get_unit_atomic_balance(unit, atom='C'):
    return (sum([i.get_atomic_flow(atom) for i in unit.ins]), 
            sum([i.get_atomic_flow(atom) for i in unit.outs]))




############################# GWP ###################################

def get_material_GWP(): # does not include natural gas as it is an invisible BT stream BT.natural_gas with price BT.natural_gas_price
    chemical_GWP = get_material_GWP_array()
    return sum(chemical_GWP)/AA.F_mass

def get_material_GWP_array():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    # chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
    chemical_GWP = [LCA_stream.imass[ID] * GWP_CF_stream.imass[ID] for ID in feed_chem_IDs]
    # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
    # return chemical_GWP.sum()/AA.F_mass
    return chemical_GWP

def get_material_GWP_breakdown():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP_dict = {ID: LCA_stream.imass[ID] * GWP_CF_stream.imass[ID] / AA.F_mass \
                         for ID in feed_chem_IDs if not LCA_stream.imass[ID] * GWP_CF_stream.imass[ID] == 0.}
    return chemical_GWP_dict

def get_material_GWP_breakdown_fractional():
    chemical_GWP_dict = get_material_GWP_breakdown()
    tot_material_GWP = get_material_GWP()
    for k,v in chemical_GWP_dict.items():
        chemical_GWP_dict[k] /= tot_material_GWP
    return chemical_GWP_dict

def get_material_GWP_breakdown_as_fraction_of_tot_GWP():
    chemical_GWP_dict = get_material_GWP_breakdown()
    tot_material_GWP = get_GWP()
    for k,v in chemical_GWP_dict.items():
        chemical_GWP_dict[k] /= tot_material_GWP
    return chemical_GWP_dict


# GWP from combustion of non-biogenic carbons
get_ng_combustion_GWP = get_onsite_gwp =\
    lambda: (s.natural_gas.get_atomic_flow('C')) * HP_chemicals.CO2.MW / AA.F_mass
                           # +ethanol_fresh.get_atomic_flow('C'))* HP_chemicals.CO2.MW / AA.F_mass

get_ng_GWP = lambda: CFs['GWP_CFs']['CH4']*s.natural_gas.F_mass/AA.F_mass

get_FGHTP_GWP = lambda: (feedstock.F_mass-feedstock.imass['Water']) \
    * CFs['GWP_CFs']['FGHTP %s'%feedstock_ID]/AA.F_mass

get_feedstock_CO2_capture = lambda: feedstock.get_atomic_flow('C')* HP_chemicals.CO2.MW/AA.F_mass
get_feedstock_GWP = lambda: get_FGHTP_GWP() - get_feedstock_CO2_capture()
# get_feedstock_GWP = lambda: get_FGHTP_GWP()
get_emissions_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions]) * HP_chemicals.CO2.MW / AA.F_mass
# GWP from electricity
get_net_electricity = lambda: sum(i.power_utility.rate for i in HP_sys.units)
get_net_electricity_GWP = lambda: get_net_electricity()*CFs['GWP_CFs']['Electricity'] \
    / AA.F_mass


get_total_electricity_demand = get_electricity_use = lambda: -BT.power_utility.rate
get_cooling_electricity_demand = lambda: u.CT.power_utility.rate + CWP.power_utility.rate

get_BT_steam_kJph_heating = lambda: sum([i.duty for i in BT.steam_utilities])
get_BT_steam_kJph_turbogen = lambda: 3600.*BT.electricity_demand/BT.turbogenerator_efficiency
get_BT_steam_kJph_total = lambda: get_BT_steam_kJph_heating() + get_BT_steam_kJph_turbogen()

get_steam_frac_heating = lambda: get_BT_steam_kJph_heating()/get_BT_steam_kJph_total()
get_steam_frac_turbogen = lambda: get_BT_steam_kJph_turbogen()/get_BT_steam_kJph_total()
get_steam_frac_cooling = lambda: get_steam_frac_turbogen() * get_cooling_electricity_demand()/get_total_electricity_demand()
get_steam_frac_electricity_non_cooling = lambda: get_steam_frac_turbogen() * (1-(get_cooling_electricity_demand()/get_total_electricity_demand()))

get_non_cooling_electricity_demand = lambda: get_electricity_demand() - get_cooling_electricity_demand()


get_EOL_GWP = lambda: AA.get_atomic_flow('C') * HP_chemicals.CO2.MW/AA.F_mass

get_direct_emissions_GWP = lambda: get_emissions_GWP() - (get_feedstock_CO2_capture() - get_EOL_GWP())

get_BT_direct_emissions_GWP = lambda: ((sum([i.get_atomic_flow('C') for i in BT.outs])*HP_chemicals['CO2'].MW / AA.F_mass)\
    /get_emissions_GWP()) * get_direct_emissions_GWP()

get_non_BT_direct_emissions_GWP = lambda: get_direct_emissions_GWP() - get_BT_direct_emissions_GWP()
                        # - (get_feedstock_CO2_capture() - get_EOL_GWP())
# get_direct_emissions_GWP = lambda: get_non_BT_direct_emissions_GWP() + get_BT_direct_emissions_GWP()
get_total_steam_GWP = lambda: get_ng_GWP() + get_BT_direct_emissions_GWP()
get_heating_demand_GWP = lambda: get_steam_frac_heating() * get_total_steam_GWP()
get_cooling_demand_GWP = lambda: get_steam_frac_cooling() * get_total_steam_GWP()
get_electricity_demand_non_cooling_GWP = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_GWP() + get_net_electricity_GWP()


# # CO2 fixed in lactic acid product
# get_fixed_GWP = lambda: \
#     AA.get_atomic_flow('C')*HP_chemicals.CO2.MW/AA.F_mass


# get_GWP = lambda: get_feedstock_GWP() + get_material_GWP() + get_ng_GWP() +\
#                   get_electricity_GWP() + get_emissions_GWP()


get_GWP = lambda: get_FGHTP_GWP() + get_material_GWP() + get_ng_GWP() +\
                  get_net_electricity_GWP() + get_direct_emissions_GWP()

get_GWP_alternative = lambda: get_FGHTP_GWP() + get_material_GWP() +\
                    get_non_BT_direct_emissions_GWP() + get_heating_demand_GWP() + get_cooling_demand_GWP() +\
                        get_electricity_demand_non_cooling_GWP()
                        
get_GWP_by_ID = lambda ID: LCA_stream.imass[ID] * GWP_CF_stream.imass[ID]/AA.F_mass


############################## FEC #################################
# Fossil energy consumption (FEC) from materials
def get_material_FEC():
    # chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
    chemical_FEC = get_material_FEC_array()
    # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
    # return chemical_FEC.sum()/AA.F_mass
    return sum(chemical_FEC)/AA.F_mass

def get_material_FEC_array():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    # chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
    chemical_FEC = [LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] for ID in feed_chem_IDs]
    # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
    # return chemical_FEC.sum()/AA.F_mass
    return chemical_FEC

def get_material_FEC_breakdown():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_FEC_dict = {ID: LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] / AA.F_mass \
                         for ID in feed_chem_IDs if not LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] == 0.}
    return chemical_FEC_dict

def get_material_FEC_breakdown_fractional():
    chemical_FEC_dict = get_material_FEC_breakdown()
    tot_material_FEC = get_material_FEC()
    for k,v in chemical_FEC_dict.items():
        chemical_FEC_dict[k] /= tot_material_FEC
    return chemical_FEC_dict

def get_material_FEC_breakdown_as_fraction_of_tot_FEC():
    chemical_FEC_dict = get_material_FEC_breakdown()
    tot_material_FEC = get_FEC()
    for k,v in chemical_FEC_dict.items():
        chemical_FEC_dict[k] /= tot_material_FEC
    return chemical_FEC_dict

get_net_electricity_FEC = lambda: \
    (get_net_electricity()*CFs['FEC_CFs']['Electricity'])/AA.F_mass

get_total_steam_FEC = lambda: get_ng_FEC()
get_heating_demand_FEC = lambda: get_steam_frac_heating() * get_total_steam_FEC()
get_cooling_demand_FEC = lambda: get_steam_frac_cooling() * get_total_steam_FEC()
get_electricity_demand_non_cooling_FEC = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_FEC() + get_net_electricity_FEC()

get_feedstock_FEC = lambda: (feedstock.F_mass-feedstock.imass['Water'])\
    * CFs['FEC_CFs']['FGHTP %s'%feedstock_ID]/AA.F_mass
# FEC from electricity

get_FEC_by_ID = lambda ID: LCA_stream.imass[ID] * FEC_CF_stream.imass[ID]/AA.F_mass

get_ng_FEC = lambda: CFs['FEC_CFs']['CH4']*s.natural_gas.F_mass/AA.F_mass
# Total FEC
get_FEC = lambda: get_material_FEC()+get_net_electricity_FEC()+get_feedstock_FEC()+get_ng_FEC()

get_FEC_alternative = lambda: get_material_FEC() + get_feedstock_FEC() + get_heating_demand_FEC() +\
    get_cooling_demand_FEC() + get_electricity_demand_non_cooling_FEC()
    
    

get_SPED = lambda: (-BT.heat_utilities[0].duty)*0.001/AA.F_mass
AA_LHV = 31.45 # MJ/kg AA

get_material_cost = lambda: sum(get_material_cost_array())


# Demand TEA contributions

get_net_electricity_VOC = lambda: HP_tea.operating_hours*sum(i.power_utility.cost for i in HP_sys.units)

get_total_steam_VOC = get_ng_cost = lambda: HP_tea.operating_hours*BT.natural_gas_price*BT.natural_gas.F_mass

get_heating_demand_VOC = lambda: get_steam_frac_heating() * get_total_steam_VOC()
get_cooling_demand_VOC = lambda: get_steam_frac_cooling() * get_total_steam_VOC()
get_electricity_demand_non_cooling_VOC = lambda: get_steam_frac_electricity_non_cooling() * get_total_steam_VOC() + get_net_electricity_VOC()

def get_VOC():
    for i in range(2):
        # AA.price = HP_tea.solve_price(AA, HP_no_BT_tea)
        HP_tea.solve_price(AA)
    return HP_tea.VOC

def get_material_cost_array():
    material_cost_array = [i.cost for i in feeds]
    material_cost_array.append(BT.natural_gas_price*BT.natural_gas.F_mass)
    return material_cost_array

def get_material_cost_breakdown():
    group_material_costs = {}
    for group in process_groups:
        group_material_costs[group.name] = 0
    counted_feeds =[]
    for feed in feeds:
        for group in process_groups:
            if group.name != 'facilities_no_hu_group':
                for unit in group.units:
                    for instream in unit.ins:
                        if instream.shares_flow_rate_with(feed) and not feed in counted_feeds:
                            group_material_costs[group.name] += feed.price*feed.F_mass/AA.F_mass
                            counted_feeds.append(feed)
    group_material_costs['BT_group'] += BT.natural_gas_price*BT.natural_gas.F_mass/AA.F_mass
    
    return group_material_costs

def get_material_cost_breakdown_fractional():
    mcb_dict = get_material_cost_breakdown()
    sum_all = sum([v for k,v in mcb_dict.items()])
    mcbf_dict = {}
    for k,v in mcb_dict.items():
        mcbf_dict[k] = mcb_dict[k]/sum_all
    return mcbf_dict    


def get_main_chem(feed):
    feed_main_chem = None
    main_chem_flow = 0.
    for chem in HP_chemicals:
        chem_ID = chem.ID
        feed_imol_chem = feed.imol[chem_ID]
        if feed_imol_chem>main_chem_flow and not (chem_ID=='H2O' or chem_ID=='Water'):
            main_chem_flow = feed_imol_chem
            feed_main_chem = chem_ID
    return feed_main_chem

def get_material_cost_breakdown_breakdown():
    group_material_costs = {}
    for group in process_groups:
        group_material_costs[group.name] = {}
    counted_feeds =[]
    for feed in feeds:
        for group in process_groups:
            if group.name != 'facilities_no_hu_group':
                for unit in group.units:
                    for instream in unit.ins:
                        if instream.shares_flow_rate_with(feed) and not feed in counted_feeds:
                            feed_main_chem = get_main_chem(feed)
                            # if feed_main_chem == 'Glucan' and group.name=='pretreatment_group':
                            #     import pdb
                            #     pdb.set_trace()
                            group_material_costs[group.name][feed_main_chem]= feed.price*feed.F_mass/AA.F_mass
                            counted_feeds.append(feed)
    group_material_costs['BT_group']['NG'] = BT.natural_gas_price*BT.natural_gas.F_mass/AA.F_mass
    
    return group_material_costs

def get_material_cost_breakdown_breakdown_fractional():
    mcbb_dict = get_material_cost_breakdown_breakdown()
    
    mcbbf_dict = {}
    for group in process_groups:
        mcbbf_dict[group.name] = {}
    for k1,v1 in mcbb_dict.items():
        v1_items = v1.items()
        sum_all = sum([v for k,v in v1_items])
        for k2, v2 in v1_items:
            mcbbf_dict[k1][k2] = v2/sum_all
    return mcbbf_dict    
# def get_GWP_breakdown():
#     group_GWPs = {}
#     for group in process_groups:
#         group_material_costs[group.name] = 0

  
# %% Full analysis
def simulate_and_print():
    MPSP = get_AA_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg AA')
    # print(f'Non-bio GWP is {get_ng_GWP():.3f} kg CO2-eq/kg AA')
    print(f'FEC is {get_FEC():.2f} MJ/kg AA or {get_FEC()/AA_LHV:.2f} MJ/MJ AA\n')
    print(f'SPED is {get_SPED():.2f} MJ/kg AA or {get_SPED()/AA_LHV:.2f} MJ/MJ AA')
    print('----------------------------------------\n')

simulate_and_print()


# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
HP_sys.diagram('cluster')

# GREET 2020

# fossil-based acrylic acid
# FEC: = 114 MJ eq. / kg AA
# GHG100 = 5.82 kgCO2 eq. / kg AA (7.65 kg CO2 eq. including EOL degradation as CO2)

# bio-based acrylic acid (via 3-HP from glycerol from algal oil)
# FEC: = 40 MJ eq. / kg AA
# GHG100 = 3.14 kgCO2 eq. / kg AA


## Ecoinvent 3.7.1 ##

# acrylic acid production
# CED-f = 49.026 MJ eq. / kg AA
# GWP-100a = 2.245 kg CO2 eq. / kg AA (4.075 including product EOL degradation as CO2)
# price = 1.4508 in 2005$ per kg
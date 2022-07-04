#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-2023, Sarang Bhagwat <sarangb2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""

@author: sarangbhagwat

Created on Sun Aug 23 12:11:15 2020

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

All units are explicitly defined here for transparency and easy reference.
Naming conventions:
    D = Distillation column
    AC = Adsorption column
    F = Flash tank or multiple-effect evaporator
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Storage
    700: Co-heat and power
    800: Cooling utility generation
    900: Miscellaneous facilities
    1000: Heat exchanger network

"""


# %% Setup


import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from math import exp as math_exp
# from biosteam import main_flowsheet as F
# from copy import deepcopy
# from biosteam import System
from thermosteam import Stream
# from biorefineries.cornstover import CellulosicEthanolTEA
from biorefineries.TAL import units, facilities
from biorefineries.TAL._process_specification import ProcessSpecification
from biorefineries.TAL.process_settings import price, CFs
from biorefineries.TAL.utils import find_split, splits_df, baseline_feedflow
from biorefineries.TAL.chemicals_data import TAL_chemicals, chemical_groups, \
                                soluble_organics, combustibles
# from biorefineries.TAL.tea import TALTEA
from biorefineries.cornstover import CellulosicEthanolTEA as TALTEA
from biosteam import SystemFactory
from warnings import filterwarnings
from biosteam.process_tools import BoundedNumericalSpecification
from scipy.interpolate import interp2d

# Based on experimental data from Singh group
ts = [0.166666667,	0.5,	1,	2]
Ts = [303.15, 323.15]
recoveries = [[0.791785714,	0.947,	0.960821429,	0.975035714],
[0.92402381,	0.956595238,	0.96297619,	0.9785]]
capacities = [[0.0739,	0.088386667,	0.089676667,	0.091003333],
[0.086242222,	0.089282222,	0.089877778,	0.091326667]]

# Interpolate
rec_interp = interp2d(ts, Ts, recoveries)
cap_interp = interp2d(ts, Ts, capacities)

filterwarnings('ignore')

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

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

def get_TAL_solubility_in_water_gpL(T):
    return get_mol_TAL_dissolved(T, 1000./18.)*TAL_chemicals['TAL'].MW

def get_K(chem_ID, stream, phase_1, phase_2):
    return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/max(1e-6, (stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol))

def get_TAL_solublity_in_solvent_very_rough(T, solvent_ID='Hexanol', units='g/L'):
    temp_stream =\
        tmo.Stream('temp_stream_get_TAL_solublity_in_solvent_very_rough')
    mol_water = mol_solvent = 1000
    mol_TAL = get_mol_TAL_dissolved(T, mol_water)
    temp_stream.imol['Water'] = mol_water
    temp_stream.imol[solvent_ID] = mol_solvent
    temp_stream.imol['TAL'] = mol_TAL
    temp_stream.lle(T=T, P=temp_stream.P)
    # temp_stream.show(N=100)
    phase_1 = 'l' if temp_stream.imol['l', solvent_ID] > temp_stream.imol['L', solvent_ID] else 'L'
    phase_2 = 'L' if phase_1=='l' else 'l'
    K_TAL_in_extract = get_K('TAL', temp_stream, phase_1, phase_2)
    # print(K_TAL_in_extract)
    if units=='g/L':
        temp_stream_2 = tmo.Stream('temp_stream_2_get_TAL_solublity_in_solvent_very_rough')
        temp_stream_2.imol['TAL'] = K_TAL_in_extract*mol_TAL
        temp_stream_2.imol[solvent_ID] = mol_solvent
        return temp_stream_2.imass['TAL']/temp_stream_2.F_vol
    elif units=='mol/mol':
        return K_TAL_in_extract*mol_TAL/(mol_TAL+mol_solvent) # 

def get_TAL_solubility_in_hexanol():
    return 2.*0.0222/(2.*0.0222+0.951) # mol/mol; 2 * Marco's initial experimental solubility of 2.8 wt% at 21 C

def get_TAL_solubility_in_ethanol_ww():
    return 0.167682 # solubility of 157.425 g-TAL per L-ethanol

# %% 
@SystemFactory(ID = 'TAL_sys')
def create_TAL_sys(ins, outs):

    # %% 
    
    # =============================================================================
    # Feedstock
    # =============================================================================
    
    # feedstock = Stream('feedstock',
    #                     baseline_feedflow.copy(),
    #                     units='kg/hr',
    #                     price=price['Feedstock'])
    
    feedstock = Stream('feedstock')
    feedstock.imass['Glucose'] = 29000.
    feedstock.imass['H2O'] = 500.
    feedstock.price = price['Glucose']*feedstock.imass['Glucose']/feedstock.F_mass
    
    feedstock.F_mass = 25802.9 # at the baseline, the amount of TAL produced would exactly satisfy the US demand for sorbic acid with a hypothetical 100% TAL->sorbic acid conversion.
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # Handling costs/utilities included in feedstock cost thus not considered here
    U101.cost_items['System'].cost = 0
    U101.cost_items['System'].kW = 0
    
    

    
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
    H301 = bst.units.HXutility('H301', ins=U101-0, T=50+273.15)
    


    M304 = bst.units.Mixer('M304', ins=(H301-0, dilution_water))
    M304.water_to_sugar_mol_ratio = 0.5
    
    @M304.add_specification()
    def adjust_M304_water():
        M304_ins_1 = M304.ins[1]
        M304_ins_1.imol['Water'] = M304.water_to_sugar_mol_ratio * M304.ins[0].imol['Glucose', 'Xylose'].sum()
        M304._run()
    # M304.specification = adjust_M304_water()
    
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.07) # split = inoculum ratio
    
    # Cofermentation
    
    R302 = units.CoFermentation('R302', 
                                    ins=(S302-1, '', CSL),
                                    outs=('fermentation_effluent', 'CO2_fermentation'))
    
    def include_seed_CSL_in_cofermentation(): # note: effluent always has 0 CSL
        R302._run()
        R302.ins[2].F_mass*=1./(1-S302.split[0])
    R302.specification = include_seed_CSL_in_cofermentation
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = units.SeedTrain('R303', ins=S302-0, outs=('seed', 'CO2_seedtrain'), ferm_ratio=0.95)
    
    T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-R302)
    
    
    
    # %% 
    
    # =============================================================================
    # Separation streams
    # =============================================================================
    
    Octanol_esterification = Stream('Octanol_esterification', units = 'kg/hr')
    
    H2_hydrogenation = Stream('H2_hydrogenation', units='kg/hr')
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
    
    Hexanol_minimal = Stream('Hexanol_minimal', units = 'kg/hr')
    # Heptane = Stream('Heptane', units = 'kg/hr')
    # Toluene = Stream('Toluene', units = 'kg/hr')
    
    # Hexanol_s = Stream('Hexanol_s', units = 'kg/hr')
    Heptane_s = Stream('Heptane_s', units = 'kg/hr')
    Toluene_s = Stream('Toluene_s', units = 'kg/hr')
    
    Hydrogen = Stream('Hydrogen', units = 'kg/hr')
    KOH = Stream('KOH', units = 'kg/hr')
    HCl = Stream('HCl', units = 'kg/hr')
    
    Ethanol_desorption = Stream('Ethanol_desorption', units='kg/hr')
    # =============================================================================
    # Separation units
    # =============================================================================
    
    # # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=R302-0, outs=('S401_solid_fraction', 'S401_liquid_fraction'),
                                # moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), 
                                solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                      'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    
    H401 = bst.units.HXutility('H401', ins=S401-1, outs = ('broth_to_adsorbtion',), T=30. + 273.15)
    
    M401 = bst.Mixer('M401', ins=(Ethanol_desorption, ''), outs=('mixed_ethanol_for_desorption'))
    S402 = bst.FakeSplitter('S402', ins=M401-0, outs=('ethanol_to_AC401', 'ethanol_to_AC2'))
    def M401_spec():
        makeup_ethanol, recycled_ethanol = M401.ins
        # AC401.run()
        M401._run()
        M401_outs_0 = M401.outs[0]
        M401_outs_0.imol['Ethanol'] = S402.outs[0].imol['Ethanol'] + S402.outs[1].imol['Ethanol']
        makeup_ethanol.imol['Ethanol'] = max(0., M401_outs_0.imol['Ethanol'] - recycled_ethanol.imol['Ethanol'])
        # S402.run()
        # M401._run()
    M401.specification = M401_spec
    
    AC401 = bst.AdsorptionColumnTSA(
        'AC401', 
        # ins=[bst.Stream('feed', TAL=0.014, Water=1, units='kg/hr', T=30 + 273.15), 'ethanol'], 
        ins=[H401-0, S402-0, 'hot_air'],
        outs=['broth_post_adsorption', 'TAL_laden_ethanol', 'ethanol_laden_air'],
        superficial_velocity=7.2, # m/h; typical velocities are 4 to 14.4 m/h for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
        
        regeneration_velocity=14.4, # m/h; default value (updated in unit specification based on titer)
        
        cycle_time=2., # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        
        # This is density of activated carbon packing, including voids.
        # So rho_adsorbent = (1 - epsilon) * rho where epsilon is the void fraction
        # and rho is the density of activated carbon with no voids.
        adsorbent='Activated carbon',
        rho_adsorbent=None, # Bulk density including void fraction; calculated based on void fraction and solid density
        rho_adsorbent_solid=700, # Solid density excluding void fraction (in kg/m3)  # Seader et al. Table 15.2
        
        void_fraction = 0.5, # v/v # Seader et al. Table 15.2
        adsorbent_capacity=0.091, # default value for unsaturated capacity (updated in unit specification); conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
        T_regeneration=30. + 273.15, 
        drying_time = 0.55, # h #!!! This is updated to 0.5 h after the first run
        T_air=TAL_chemicals.Ethanol.Tb + 10, # K
        air_velocity = 2160, # m/h
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        regeneration_fluid=dict(phase='l', Ethanol=1., units='kg/hr'),
        adsorbate_ID='TAL',  
        split=dict(TAL=0, Water=1, VitaminA=1., VitaminD2=1., FermMicrobe=1.),
        length_unused = 1.219, # m; 4 ft based on recommendation by Seader et al. (Separation Process Principles)
        target_recovery=0.99,
        wet_retention=0.5, # conservatively assume half a wash's worth of ethanol is retained in the column before dry air is passed through it
        K = 0.07795, # back-calculated for 1 wash from experimental measurements for 3 washes pooled together; 0.125 for 3-wash # constant desorption partition coefficient; calculated for 1 wash from experimental data for 3 washes pooled together
    )
    AC401._default_equipment_lifetime['Activated carbon'] = 1.
    AC401.adsorbent_cost['Activated carbon'] = price['Activated carbon'] # 41. $/ft^3
    
    @AC401.add_specification
    def AC401_spec(): # update recovery and capacity based on user-input adsorption time and temperature
        
        # T = AC401.ins[0].T
        # t = AC401.cycle_time
        # capacity = cap_interp(t, T)
        # AC401.adsorbent_capacity = capacity[0]
        
        AC401._run()
        
        M401.run()
        
        AC401.ins[1].T = M401.outs[0].T
        
        # some moisture retained
        AC401.outs[1].imol['Water'] = AC401.outs[1].imol['TAL']
        AC401.outs[0].imol['Water'] -= AC401.outs[1].imol['TAL']
    

        
    F401 = bst.units.MultiEffectEvaporator('F401', ins=AC401-1, outs=('F401_b', 'F401_t'), chemical='Ethanol',
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.7)

    F401.flash=False
    F401.TAL_solubility_in_ethanol_ww = get_TAL_solubility_in_ethanol_ww()
    def F401_obj_fn(V):
        F401_b = F401.outs[0]
        # F401_ins_0 = F401.ins[0]
        # TAL_mass = F401_ins_0.imass['TAL']
        # F401_ins_0.imass['TAL'] = 0.
        F401.V = V
        F401._run()
        # F401_ins_0.imass['TAL']  =TAL_mass
        # F401_b.imass['TAL'] = TAL_mass

        return F401.TAL_solubility_in_ethanol_ww - F401_b.imass['TAL']/F401_b.F_mass

    F401.specification = BoundedNumericalSpecification(F401_obj_fn, 1e-4, 1.-1e-4, ytol=1e-4)
    
    P401 = bst.Pump('P401', ins=F401-1, P=101325.)
    # H404 = bst.units.HXutility('H404', ins=P401-0, outs=('cooled_ethanol'), 
    #                            T=30.+273.15, rigorous=True)
    
    # M403 = bst.Mixer('M403', ins=(AC401-1, AC2-1))
    F402 = bst.units.Flash('F402', ins=F401-0, outs=('F402_t', 'F402_b'), P=101325.,
                            V=0.5) # !!! TODO: replace with dryer
    
    F402.product_ethanol_content = 0.05 # g-ethanol per g-TAL, not per g-product
    def F402_spec():
            F402_b = F402.outs[1]
            F402_ins_0 = F402.ins[0]
            water_mol = float(F402_ins_0.imol['Water'])
            F402_ins_0.imol['Water'] = 0.
            TAL_mass = F402_ins_0.imass['TAL']
            TAL_mol = F402_ins_0.imol['TAL']
            F402_ins_0.imol['TAL'] = 0.
            F402.V = 1. - (F402.product_ethanol_content*TAL_mass)/(46.06844*F402_ins_0.imol['Ethanol'])
            
            F402._run()
            F402_ins_0.imol['TAL'] = TAL_mol
            F402_b.imol['TAL'] = TAL_mol      
            F402_b.imol['Water'] = water_mol 
            
    F402.specification = F402_spec
    
    H404 = bst.units.HXutility('H404', ins=F402-0, outs=('cooled_ethanol_F402'), 
                               T=30.+273.15, rigorous=True)
    
    H403 = bst.units.HXutility('H403', ins=F402-1, outs=('cooled_TAL'), 
                               T=30.+273.15, rigorous=True)
    
    
    # M404 = bst.Mixer('M402', ins=(AC401-2, AC2-2))
    H402 = bst.units.HXutility(
        'H402', ins=AC401-2, outs=('cooled_ethanol_laden_air'), 
        T=265.,
        rigorous=True
    )
    
    S403 = bst.units.FakeSplitter('S403', ins=H402-0, outs=('cool_air', 'ethanol_recovered_from_air'))
    def S403_spec():
        S403_ins_0 = S403.ins[0]
        S403.outs[0].mol[:] = S403_ins_0['g'].mol[:]
        S403.outs[1].mol[:] = S403_ins_0['l'].mol[:]
    S403.specification = S403_spec
    # S403-0-2-M404
    M402 = bst.Mixer('M402', ins=(P401-0, H404-0, S403-1), outs=('recycled_ethanol',))
    M402-0-1-M401
    
    
    
    M403 = bst.Mixer('M403', ins=(H403-0, Octanol_esterification),)
    
    @M403.add_specification(run=True)
    def M403_spec():
        M403.ins[1].imol['Octanol'] = 1.2*R401.TAL_to_esters_conversion*M403.ins[0].imol['TAL']
    
    R401 = units.HydrogenationEstersReactor('R401', ins=(M403-0, 'recycled_reactants', H2_hydrogenation), 
                                            outs=('mixed_hydrogenation_products', 'vented_gas'),
                                            tau = 2.333,
                                            P=34.5423*101325.)
    
    F403 = bst.Flash('F403', ins=R401-0, outs = ('volatiles', 'bottom_product_esters'), V = 1.-1e-5, P=101325./100.) 
    #!!! TODO: esters solvent extraction
    
    # S406 = bst.MultiStageMixerSettlers('S406', ins=S403-0, outs=('raffinate', 'extract'),
    #         partition_data={
    #         'K': np.array([30., 30., 1e5, 1e-5]),
    #         'IDs': ('Octyl_5_hydroxyhexanoate', 'Octyl_3_5_dihydroxyhexanoate', 'Ethanol', 'Water',),
    #         'phi' : 0.5,
    #     },
    #     N_stages = 4,
    # )
    
    # S406 = bst.FakeSplitter('S406', ins=M403-0, outs=('extract', 'raffinate'),
    #                         )
    
    # @S406.add_specification
    # def S406_spec():
    #     S406.outs[1].copy_like(S406.ins[0])
    #     S406.outs[0].imol['TAL'] = S406.ins[0].imol['TAL']
    #     S406.outs[1].imol['TAL'] -= S406.outs[0].imol['TAL']
        
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
    M501 = bst.units.Mixer('M501', ins=(
                                        AC401-0,
                                        # S403-1,
                                        # AC401-0,
                                        # S402-1,
                                        F403-0,
                                        # r_S402_s-1, r_S403_s-1, r_S404_s-1,
                                        # X401-1, S408-0,
                                        ))
    
    # This represents the total cost of wastewater treatment system
    WWT_cost = units.WastewaterSystemCost('WWTcost501', ins=M501-0)
    
    R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                    outs=('biogas', 'anaerobic_treated_water', 
                                          'anaerobic_sludge'),
                                    reactants=soluble_organics + ['TAL'],
                                    split=find_split(splits_df.index,
                                                     splits_df['stream_611'],
                                                     splits_df['stream_612'],
                                                     chemical_groups),
                                    T=35+273.15)
    
    get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
    
    # Mix recycled stream and wastewater after R501
    M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
    @M502.add_specification(run=True)
    def M502_spec():
        M503._run()
        
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
    
    # Mix solid wastes to boiler turbogenerator
    M505 = bst.units.Mixer('M505', ins=(S503-1, 
                                        S401-0, 
                                        # F403-0, D401-0,
                                        ), 
                            outs='wastes_to_boiler_turbogenerator')
    
    
    # %% 
    
    # =============================================================================
    # Facilities streams
    # =============================================================================
    
    ethanol_fresh = Stream('ethanol_fresh',  price=price['Ethanol'])
    
    sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
    # TCP_fresh = Stream('TCP_fresh',  price=price['TCP'])
    
    ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
    CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
    # lime_fresh = Stream('lime_fresh', price=price['Lime'])
    
    HCl_fresh = Stream('HCl_fresh', price=price['HCl'])
    
    octanol_fresh = Stream('octanol_fresh', price=price['Hexanol'])
    H2_fresh = Stream('H2_fresh', price=price['Hydrogen'])
    # heptane_fresh = Stream('heptane_fresh', price=price['Heptane'])
    # toluene_fresh = Stream('toluene_fresh', price=price['Toluene'])
    
    # hexanol_fresh_s = Stream('hexanol_fresh_s', price=price['Hexanol'])
    heptane_fresh_s = Stream('heptane_fresh_s', price=price['Heptane'])
    toluene_fresh_s = Stream('toluene_fresh_s', price=price['Toluene'])
    
    hydrogen_fresh = Stream('hydrogen_fresh', price=price['Hydrogen'])
    KOH_fresh = Stream('KOH_fresh', price=price['KOH'])
    # S401_out1_F_mass = S401.outs[1].F_mass
    
    # if not (S401_out1_F_mass == 0):
    #     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M403.ins[3].imass['Ethanol']
    #     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M403.ins[3].imass['Dipotassium hydrogen phosphate']
        
    # else:
    # ethanol_fresh = Stream('ethanol_fresh', Ethanol = get_feedstock_dry_mass()*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
    # DPHP_fresh = Stream('DPHP_fresh', DPHP = get_feedstock_dry_mass()*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    # TAL stream
    # TAL = Stream('TAL', units='kg/hr', price=price['TAL'])
    # SA product
    Mixed_esters = Stream('Mixed_esters', units='kg/hr', price=price['SA'])
    # Acetoin product
    # Acetoin = Stream('Acetoin', units='kg/hr', price=price['Acetoin'])
    # # Isobutyraldehyde product
    # IBA = Stream('IBA', units='kg/hr', price=price['IBA'])
    # Chemicals used/generated in BT
    # FGD_lime = Stream('FGD_lime')
    ash = Stream('ash', price=price['Ash disposal'])
    # boiler_chems = Stream('boiler_chems', price=price['Boiler chems'])
    # baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
    # Supplementary natural gas for BT if produced steam not enough for regenerating
    # all steam streams required by the system
    # natural_gas = Stream('natural_gas', price=price['Natural gas'])
    
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
    
    T601 = bst.units.StorageTank('T601', ins=octanol_fresh,
                                          # outs=Octanol_esterification,
                                           )
    T601.line = 'Octanol storage tank'
    T601_P = units.TALPump('T601_P', ins=T601-0, outs = Octanol_esterification)
    
    T602 = bst.units.StorageTank('T602', ins=H2_fresh,
                                          # outs=H2_esterification,
                                           )
    T602.line = 'H2 storage tank'
    T602_P = units.TALPump('T602_P', ins=T602-0, outs = H2_hydrogenation)
    
    
    # S601 = bst.units.ReversedSplitter('S601', ins=T601-0, 
    #                                   outs=(pretreatment_sulfuric_acid, 
    #                                         ''))
    # T608 = units.TCPStorageTank('T608', ins=TCP_fresh,
    #                                      outs='TCP_catalyst')
    # T608-0-3-R401
    # T608.line = 'Tricalcium diphosphate storage tank'
    #
    # T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia_M205)
    # T602.line = 'Ammonia storage tank'
    
    T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
    T603.line = 'CSL storage tank'
    
    # Ethanol storage
    T604 = bst.units.StorageTank('T604', ins=ethanol_fresh)
    T604.line = 'Ethanol storage tank'
    T604_P = units.TALPump('T604_P', ins=T604-0, outs = Ethanol_desorption)
    # T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0, outs = Hexanol)

    # T607_P = units.TALPump('T607_P', ins=T607-0, outs = Hydrogen)
    
    # Connections to ATPE Mixer
    # T604_P-0-1-M403
    # T605_P-0-2-M403
    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T620 = units.TALStorageTank('T620', ins=F403-1, tau=7*24, V_wf=0.9,
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    
    
    
    T620.line = 'EstersStorageTank'
    
    
    T620_P = units.TALPump('T620_P', ins=T620-0, outs=Mixed_esters)
    
    
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
    
    
    CIP = facilities.CIP('CIP901', ins=CIP_chems_in, outs='CIP_chems_out')
    ADP = facilities.ADP('ADP902', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_tpd()/2205)
    
    
    FWT = units.FireWaterTank('FWT903', ins=fire_water_in, outs='fire_water_out')
    
    CWP = facilities.CWP('CWP802', ins='return_chilled_water',
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
    
    BT = bst.facilities.BoilerTurbogenerator('BT701',
                                                      ins=(M505-0,
                                                          R501-0, 
                                                          'boiler_makeup_water',
                                                          'natural_gas',
                                                          'lime',
                                                          'boilerchems'), 
                                                      outs=('gas_emission', 'boiler_blowdown_water', ash,),
                                                      turbogenerator_efficiency=0.85,
                                                      natural_gas_price=price['Natural gas'])
    
    # BT = bst.BDunits.BoilerTurbogenerator('BT',
    #                                    ins=(M505-0, R501-0, 'boiler_makeup_water', 'natural_gas', FGD_lime, boiler_chems),
    #                                    boiler_efficiency=0.80,
    #                                    turbogenerator_efficiency=0.85)
    
    # Blowdown is discharged
    CT = facilities.CT('CT801', ins=('return_cooling_water', cooling_tower_chems,
                                  'CT_makeup_water'),
                       outs=('process_cooling_water', 'cooling_tower_blowdown'))
    
    # All water used in the system, here only consider water usage,
    # if heating needed, then heeating duty required is considered in BT
    process_water_streams = (enzyme_water,
                             aerobic_caustic, 
                             CIP.ins[-1], BT.ins[-1], CT.ins[-1])
    
    PWC = facilities.PWC('PWC904', ins=(system_makeup_water, S504-0),
                         process_water_streams=process_water_streams,
                         recycled_blowdown_streams=None,
                         outs=('process_water', 'discharged_water'))
    
    # Heat exchange network
    HXN = bst.facilities.HeatExchangerNetwork('HXN1001',
                                              ignored=lambda:[
                                                       H401, 
                                                       H402, 
                                                       H403, 
                                                       H404,
                                                       AC401.heat_exchanger_drying,
                                                       AC401.heat_exchanger_regeneration,
                                                       F401.components['condenser'],
                                                       ],
                                                cache_network=True,
                                                force_ideal_thermo=True,
                                               )

    # HXN = HX_Network('HXN')

# %% 

# =============================================================================
# Complete system
# =============================================================================

TAL_sys = create_TAL_sys()


f = bst.main_flowsheet
u = f.unit
s = f.stream

feedstock = s.feedstock
Mixed_esters = s.Mixed_esters
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185



TEA_feeds = set([i for i in TAL_sys.feeds if i.price]+ \
    [i for i in TAL_sys.feeds if i.price])

TEA_products = set([i for i in TAL_sys.products if i.price]+ \
    [i for i in TAL_sys.products if i.price]+[Mixed_esters])
    
    
for ui in u:
    globals().update({ui.ID: ui})
# %%
# =============================================================================
# TEA
# =============================================================================

# TAL_tea = CellulosicEthanolTEA(system=TAL_sys, IRR=0.10, duration=(2016, 2046),
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
#                     u.CWP, u.CT, u.PWC, u.CIP, u.ADP, u.FWT, u.BT),
#         warehouse=0.04, site_development=0.09, additional_piping=0.045,
#         proratable_costs=0.10, field_expenses=0.10, construction=0.20,
#         contingency=0.10, other_indirect_costs=0.10, 
#         labor_cost=3212962*get_flow_tpd()/2205,
#         labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
#         steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)

# TAL_no_BT_tea = TAL_tea

TAL_tea = TALTEA(system=TAL_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWTcost501,
                    # u.T601, u.T602, 
                    u.T603, u.T604, u.T620,
                    # u.T606, u.T606_P,
                    u.CWP802, u.CT801, u.PWC904, u.CIP901, u.ADP902, u.FWT903, u.BT701),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

TAL_no_BT_tea = TAL_tea

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

#%% Define unit groups

area_names = [
    'feedstock',
    # 'pretreatment',
    'conversion',
    'separation',
    'wastewater',
    'storage',
    'co-heat and power',
    'cooling tower and chilled water package',
    'other facilities',
    'heat exchanger network',
]
# u.CWP901.ID = 'CWP802' 
for ui in u:
    if type(ui) == bst.ChilledWaterPackage:
        ui.ID = 'CWP802' # group with CT for system cooling demand
        break
unit_groups = bst.UnitGroup.group_by_area(TAL_sys.units)
unit_groups.append(bst.UnitGroup('natural gas'))


for i, j in zip(unit_groups, area_names): i.name = j
for i in unit_groups: i.autofill_metrics(shorthand=True, 
                                         electricity_production=True, 
                                         material_cost=True)
for i in unit_groups:
    if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling tower and chilled water package':
        i.metrics[-1].getter = lambda: 0. # Material cost
    if i.name == 'cooling tower and chilled water package':
        i.metrics[1].getter = lambda: 0. # Cooling duty
HXN = None
for HXN_group in unit_groups:
    if HXN_group.name == 'heat exchanger network':
        HXN_group.filter_savings = False
        HXN = HXN_group.units[0]
        assert isinstance(HXN, bst.HeatExchangerNetwork)
        
unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Mat. cost', 
                                                    getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
                                                    units='USD/hr',
                                                    element=None)

unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i
# HXN.force_ideal_thermo = True
CT = u.CT801
BT = u.BT701
CWP = u.CWP802


# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================


def get_Mixed_esters_MPSP():
    for i in range(3):
        TAL_sys.simulate()
    for i in range(3):
        Mixed_esters.price = TAL_tea.solve_price(Mixed_esters)
    return Mixed_esters.price*Mixed_esters.F_mass/sum(Mixed_esters.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])

spec = ProcessSpecification(
    evaporator = None,
    pump = None,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = [],
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('TAL',),
    
    spec_1=0.19,
    spec_2=15.,
    spec_3=0.19,

    
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = [],
    HXN = u.HXN1001,
    maximum_inhibitor_concentration = 1.,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301], # if the line below does not work (depends on BioSTEAM version)
    pre_conversion_units = TAL_sys.split(u.M304.ins[0])[0],
    
    # set baseline fermentation performance here
    baseline_yield = 0.19,
    baseline_titer = 15.,
    baseline_productivity = 0.19,
    
    # baseline_yield = 0.30,
    # baseline_titer = 25.,
    # baseline_productivity = 0.19,
    
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = None)


spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

def M304_titer_obj_fn(water_to_sugar_mol_ratio):
    M304, R302 = u.M304, u.R302
    M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
    M304.specification[0][0]()
    u.M304_H._run()
    u.S302._run()
    u.R303._run()
    u.T301._run()
    R302.specification[0][0]()
    # broth = R302.outs[0]
    # return broth.imass['TAL']/broth.F_vol - R302.titer_to_load
    return R302.effluent_titer - R302.titer_to_load

def load_titer_with_glucose(titer_to_load):
    spec.spec_2 = titer_to_load
    u.R302.titer_to_load = titer_to_load
    flx.IQ_interpolation(M304_titer_obj_fn, 1e-3, 20000.)
    # u.AC401.regeneration_velocity = min(14.4, 3.1158 + ((14.4-3.1158)/(30.-3.))*(titer_to_load-3.)) # heuristic to obtain regeneration velocity at which MPSP is minimum fitted to results from simulations at target_recovery=0.99 
    # u.AC401.regeneration_velocity = 14.4
    
spec.load_spec_2 = load_titer_with_glucose

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


# %% Full analysis
def simulate_and_print():
    get_Mixed_esters_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_Mixed_esters_MPSP():.3f}/kg')
    # print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg SA')
    # print(f'FEC is {get_FEC():.2f} MJ/kg SA or {get_FEC()/SA_LHV:.2f} MJ/MJ SA')
    # print(f'SPED is {get_SPED():.2f} MJ/kg SA or {get_SPED()/SA_LHV:.2f} MJ/MJ SA')
    # print('--------------------\n')

# simulate_and_print()
# TAL_sys.simulate()
get_Mixed_esters_MPSP()
# u.AC401.cycle_time = 4.
# u.AC401.drying_time = 0.5 #!!! Drying time is updated to this value (overwritten the value passed during initialization)
spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
simulate_and_print()

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
    # 'separation_sys': (S401, M403, M403_P,
    #                     S402, 
    #                     # F403, F403_H, X401,
    #                     D401, D401_H, D401_P, S403,
    #                     M402_P, S403,
    #                     D403, D403_H, D403_P,
    #                     M501,
    #                     T606, T606_P, T607, T607_P)
                        # F402, F402_H, F402_P,
                        # D405, D405_H1, D405_H2, D405_P,
                        # M403, M403_P)
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

#%% TEA breakdown

def TEA_breakdown(print_output=False):
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]()
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()
                    
                    
            # if ug.name=='natural gas':
            #     if metric.name=='Mat. cost':
            #         metric_breakdowns[metric.name][ug.name] = BT.natural_gas.F_mass*BT.natural_gas_price
            
            
            # else:
            #     storage_metric_val = metric()
                
    # print and return metric_breakdowns
    if print_output:
        for i in unit_groups[0].metrics:
            print(f"\n\n----- {i.name} ({i.units}) -----")
            metric_breakdowns_i = metric_breakdowns[i.name]
            for j in metric_breakdowns_i.keys():
                print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")
    return metric_breakdowns

TEA_breakdown()
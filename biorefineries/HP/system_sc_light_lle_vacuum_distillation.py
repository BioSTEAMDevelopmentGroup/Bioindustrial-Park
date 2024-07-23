#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# !!! Package name and short description goes here.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
from numba import njit
from biosteam import main_flowsheet as F
from biosteam import System
from thermosteam import Stream
from biorefineries.HP import units, facilities
from biorefineries.HP.lca import LCA
from biorefineries.HP._process_specification import ProcessSpecification
from biorefineries.HP.process_settings import price, CFs, chem_index
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
from biorefineries.cellulosic import create_facilities
# from lactic.hx_network import HX_Network

IQ_interpolation = flx.IQ_interpolation
# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# bst.speed_up()
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
# System.default_converge_method = 'wegstein'
feedstock_ID = 'Corn stover'

# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'
# System.default_converge_method = 'wegstein'
# System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.default_relative_molar_tolerance = 0.001 # supersedes absolute tolerance
System.default_molar_tolerance = 0.1
System.strict_convergence = True # True => throw exception if system does not converge; false => continue with unconverged system


@SystemFactory(ID = 'HP_sys')
def create_HP_sys(ins, outs):
    u, s = flowsheet.unit, flowsheet.stream
    process_groups = []
    # %% Feedstock
    
    # Sugarcane juicing subprocess
    sugarcane_juicing_sys = create_juicing_system_up_to_clarification()
    
    u = sugarcane_juicing_sys.flowsheet.unit
    s = sugarcane_juicing_sys.flowsheet.stream
    
    u.U201.diagram()
    sugarcane_juicing_sys.diagram('cluster')
    # u.U201.ins.append(u.M201-0)
    
    # u.M201-0-1-u.U201
    
    # sugarcane_juicing_sys.simulate(update_configuration=True)
    
    U101 = bst.Unit('U101', ins='', outs='')
    @U101.add_specification(run=False)
    def U101_spec():
        U101.outs[0].copy_like(U101.ins[0])
    
    feedstock = s.sugarcane
    feedstock_sink = feedstock.sink
    U101-0-0-feedstock_sink
    feedstock-0-U101
    feedstock_sink.ins[0].price = 0.
    
    feedstock.F_mass = 554171.74 # initial value; updated by spec.set_production_capacity
    
    # Update all prices to 2019$ using chemical indices
    # sugarcane biorefinery base year is 2019
    for sugarcane_sys_stream in list(s):
        sugarcane_sys_stream.price *= chem_index[2019]/chem_index[2019]
        
    #%% Fermentation
    fermentation_sys = create_HP_fermentation_process(ins=(u.C201-0),
                                                   )
    
    
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
    S401 = bst.units.SolidsCentrifuge('S401', ins=(R302-0,), outs=('cell_mass', ''),
                                moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    fix_split(S401.isplit, 'Glucose')
    # NOTE: if there is not enough moisture content, it is impossible to pump
    # the fluid into the centrifuge; in fact, the centrifuge would not be able
    # to separate anything.
    
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
    
    @S402.add_specification(run=False)
    def S402_spec():
        if S402.ins[0].imol['CaSO4']>0:
            S402._run()
        else:
            S402.outs[0].mol[:] = 0
            S402.outs[1].mol = S402.ins[0].mol
        
        
    # S402.specification = S402_spec
    
    
    M401 = bst.units.Mixer('M401', ins=(separation_hexanol,
                                        ''))
    
    M401_H = bst.units.HXutility('M401_H', ins = M401-0, T = 80. + 273.15, rigorous = False)
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=S402-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
    target_HP_x = 0.10
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol['AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'])
    
    @F401.add_specification(run=False)
    def F401_specification():
        instream = F401.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        ratio = get_x('HP', instream)/target_HP_x
        # no need to check for ratio>1 becasue our target_water_x is consistently lower than the max possible titer
        F401.V = 1. - ratio
        
        F401._run()
    
    # F401.specification = F401_specification
    
    
    def F401_no_run_cost():
        F401.heat_utilities = tuple()
        F401._installed_cost = 0.
    # F401._cost = F401_no_run_cost
    
    F401_P = bst.units.Pump('F401_P', ins=F401-0, P=101325.)    
    F401_H = bst.units.HXutility('F401_H', ins = F401_P-0, T = 80. + 273.15, rigorous = False)
    

    
    Kds = dict(IDs=('HP', 'Water', 'Hexanol', 'AceticAcid'),
                # K=np.array([1./1.9379484051844278, 3.690183610720956, 0.0060176892697821486, 1./0.4867537504125923]), # T = 80. + 273.15 K
                K=np.array([1.9379484051844278, 1/3.690183610720956, 1/0.0060176892697821486, 0.4867537504125923]), # T = 80. + 273.15 K
               phi = 0.5)
    
    
    max_N_stages = 3
    
    S404 = bst.units.MultiStageMixerSettlers('S404', ins = (F401_H-0, M401_H-0),
                                         outs = ('extract', 'raffinate'),
                                         N_stages = max_N_stages, 
                                          partition_data = Kds,
                                         ) 
                              
                              
    S404.vol_frac = 0.05
    
    
    tolerable_loss_fraction = 0.001
    
    # @S404.add_specification(run=False)
    def S404_spec_hexanol():
        feed_hexanol, solvent_recycle = M401.ins
        req_hexanol = S404.ins[0].imol['HP'] * 10.
        feed_hexanol.imol['Hexanol'] = max(0, req_hexanol - solvent_recycle.imol['Hexanol'])
        M401._run()
        M401_H._run()
        S404._run()
        
    @S404.add_specification(run=False)
    def adjust_S404_streams():
        S404.N_stages = max_N_stages # reset
        S404._setup() # reset
        feed_hexanol, solvent_recycle = M401.ins
        process_stream = S404.ins[0]
        process_stream_F_mol = process_stream.F_mol
        existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
    
        # K_raffinate = S404.partition_data['K'][0]
        K_raffinate = 1./S404.partition_data['K'][0]
    
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
            # feed_hexanol.imol['Hexanol'] = S404.outs[0].imol['Hexanol']
            feed_hexanol.imol['Hexanol'] = S404.outs[1].imol['Hexanol']
        
        for i in S404.outs: i.T = M401_H.outs[0].T
        
    def update_Ks(lle_unit, solute_indices = (0,), carrier_indices = (1,), solvent_indices = (2,)):
        IDs = lle_unit.partition_data['IDs']
        # Ks = lle_unit.partition_data['K']
        Ks = 1./lle_unit.partition_data['K']
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
        Ks_new = (test_stream['l'].imol[IDs]/test_stream['l'].F_mol)/(test_stream['L'].imol[IDs]/test_stream['L'].F_mol)
        
        return Ks_new
    
    def S404_run():
        try:
            S404._run()
            
            if has_negative_flows(S404):
                raise InfeasibleRegion('negative flows')
        except:
            # import pdb
            # pdb.set_trace()
            S404.N_stages-=1
            if S404.N_stages == 0:
                S404.N_stages = max_N_stages # reset
                S404._setup() # reset
                raise InfeasibleRegion('number of stages in %s'%(S404.ID))   
            else:
                S404._setup()
                # print('\nReduced S404.N_stages to %s\n'%S404.N_stages)
            S404_run()
        
    
    def has_negative_flows(unit):
        for stream in unit.outs + unit.ins:
            if (stream.mol < 0.).any():
                return True
        return False
                
    # S404.specification = adjust_S404_streams
    
    S404_spec = S404.specifications[0]
    globals().update({'S404_spec': S404_spec})
    
    ideal_thermo = S404.thermo.ideal()
    
    D401 = bst.units.BinaryDistillation('D401', ins=S404-0, outs=('D401_g', 'D401_l'),
                                        LHK=('Hexanol', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P = 101325./20.,
                                        vessel_material = 'Stainless steel 316',
                                        partial_condenser = False,
                                        # condenser_thermo = ideal_thermo,
                                        # boiler_thermo = ideal_thermo,
                                        thermo=ideal_thermo)
    
    @D401.add_specification(run=True)
    def D401_spec():
        D401.ins[0].phase = 'l'
    # D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0., rigorous=True)
    D401_H_P = units.HPPump('D401_H_P', ins=D401-0, P = 101325)
    D401_H_P-0-1-M401
    
    def get_concentration_gpL(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_vol
    
    def get_mass_percent(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_mass
    
    @njit(cache=True)
    def mass_percent_helper(chemmass, totmass):
        return chemmass/totmass
    
    D401_bP = bst.Pump('D401_bP', ins=D401-1, P=101325.)
    
    M402 = bst.units.Mixer('M402', ins=(D401_bP-0,
                                        'dilution_water2'))
    
    M402.HP_wt_frac = 0.3
    
    def M402_objective_fn(Water_imol):
        M402.ins[1].imol['Water'] = Water_imol
        M402._run()
        # return get_concentration_gpL('HP', M402.outs[0]) - 600. # predicted "solubility" of 645 g/L at STP https://hmdb.ca/metabolites/HMDB0000700
        # return get_concentration_gpL('HP', M402.outs[0]) - 270.1 # "Solubility "at 25 C # https://www.chemicalbook.com/ChemicalProductProperty_EN_CB6711580.htm
        # return get_mass_percent('HP', M402.outs[0]) - .15 # Dehydration reaction paper
        # return get_mass_percent('HP', M402.outs[0]) - .35 # 30-35 wt% in https://patents.google.com/patent/WO2013192451A1/en
        # return get_mass_percent('HP', M402.outs[0]) - .30 # 30wt% with 80% conversion in Dunn et al. 2015 LCA of Bioproducts in GREET
        # return get_mass_percent('HP', M402.outs[0]) - .99 # ideal
        
        return get_mass_percent('HP', M402.outs[0]) - M402.HP_wt_frac
    
    @M402.add_specification(run=False)
    def M402_adjust_water():
        flx.IQ_interpolation(M402_objective_fn, 0., 10000, maxiter=50, ytol=1e-2)
        
    # M402.specification = M402_adjust_water
    
    D401_P = units.HPPump('D401_P', ins=M402-0, P=506625.*5.4)
    
    R402 = units.DehydrationReactor('R402', ins = (D401_P-0, makeup_TiO2_catalyst, '',),
                                    outs = ('dilute_acryclic_acid', 'spent_TiO2_catalyst'),
                                    tau = 57.34/1.5, # Dishisha et al.
                                    T = 230. + 273.15,
                                    P = D401_P.P,
                                    vessel_material='Stainless steel 316')
    # R402.heat_utilities[0].heat_transfer_efficiency = 1. 
    
    # spent_TiO2_catalyst is assumed to be sold at 0 $/kg

    R402_H = bst.units.HXutility('R402_H', ins=R402-0, T = 89. + 273.15, rigorous=True)
    
    
    D402 = bst.units.ShortcutColumn('D402', ins=R402_H-0, outs=('D402_g', 'D402_l'),
                                        LHK=('Water', 'AcrylicAcid'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P=101325./10.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    
    D402_dP = bst.Pump('D402_dP', ins=D402-0, outs=('D402_wastewater',), P=101325.)
    # recycling water makes system convergence fail
    # D402-0-3-R402
    
    
    D402_P = units.HPPump('D402_P', ins=D402-1, P=101325.)
    
    
    # D402_H = bst.units.HXutility('D402_H', ins=D402_P-0, T = 330., rigorous=True)
    
    
    D403 = bst.units.ShortcutColumn('D403', ins=D402_P-0, outs=('D403_g', 'D403_l'),
                                        LHK=('AcrylicAcid', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.05, P=101325/20.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    D403_dP = units.HPPump('D403_dP', ins=D403-0, P=101325.)
    D403_bP = units.HPPump('D403_bP', ins=D403-1, P=101325.)
    D403_bP-0-2-R402
    
    D403_H = bst.units.HXutility('D403_H', ins=D403_dP-0, T = 25.+273.15, rigorous=True)
    D403_P = units.HPPump('D403_P', ins=D403_H-0)
    
    separation_group = UnitGroup('separation_group', 
                                   units=(S401, R401, R401_H, R401_P, S402,
                                          F401, F401_P, M401, S404, D401,
                                          D403_dP, D403_bP, D402_dP, D401_H_P, D401_P, M402, 
                                          R402, R402_H, D402, D402_P, D403, D403_H, D403_P))
    process_groups.append(separation_group)
    

    #%%# !!!
    
    
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
    M501 = bst.units.Mixer('M501', ins=(F301_P-0, D402_dP-0, F401-1, S404-1,
                                        H201-0,
                                        ))
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        kind='conventional',
        ins=M501-0,
        mockup=True,
        area=500,
    )
    
    # Mix solid wastes to boiler turbogenerator
    M510 = bst.units.Mixer('M510', ins=(
                                    # S503-1,
                                    S301-0, S401-0),
                            outs='wastes_to_boiler_turbogenerator')
    
    MX = bst.Mixer(900, ['', ''])
    
    M503 = u.M503
    @M503.add_specification(run=False)
    def M503_spec():
        for i in M503.ins: i.phase='l'
        M503._run()
        for j in M503.outs: j.phase='l'
    
    WWT_group = UnitGroup('WWT_group', 
                                   units=wastewater_treatment_sys.units)
    process_groups.append(WWT_group)
    
    # %% 
    
    # =============================================================================
    # Facilities streams
    # =============================================================================
    
    sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
    sulfuric_acid_fresh2 = Stream('sulfuric_acid_fresh2',  price=price['Sulfuric acid'])
    ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
    CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
    lime_fresh = Stream('lime_fresh', price=price['Lime'])
    
    hexanol_fresh = Stream('hexanol_fresh', price=price['Hexanol'])
    TOA_fresh = Stream('TOA_fresh', price=price['TOA'])
    AQ336_fresh = Stream('AQ336_fresh', price=price['AQ336'])
    
    
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
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
    
    
    #%% new
    imbibition_water = Stream('imbibition_water', price=price['Makeup water'])
    rvf_wash_water = Stream('rvf_wash_water', price=price['Makeup water'])
    
    #%%
    # =============================================================================
    # Facilities units
    # =============================================================================
    
    T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh,
                                         outs=sulfuric_acid_T201)
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


    ############################
    
    create_facilities(
        solids_to_boiler=M510-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=[
         imbibition_water,
         rvf_wash_water,
         dilution_water,
         system_makeup_water,
         # s.fire_water,
         # s.boiler_makeup_water,
         # s.CIP,
         # s.recirculated_chilled_water,
         # s.s.3,
         # s.cooling_tower_makeup_water,
         # s.cooling_tower_chemicals,
         ],
        feedstock=feedstock,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX-0,
        BT_area=700,
        area=900,
    )
    
    # CWP803 = bst.ChilledWaterPackage('CWP803', agent=bst.HeatUtility.cooling_agents[-2])
    
    BT = u.BT701
    BT.natural_gas_price = price['Natural gas']
    BT.ins[4].price = price['Lime']
    
    HXN = bst.HeatExchangerNetwork('HXN1001',
                                                ignored=
                                                [
                                                 # D401,
                                                 # D403,
                                                 ],
                                              # cache_network=True,
                                              )
    
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
    
    # # To simulate without HXN, simply uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []
    
    HXN_group = UnitGroup('HXN_group', 
                                   units=(HXN,))
    process_groups.append(HXN_group)
    
    BT = u.BT701
    CT = u.CT901
    # except: breakpoint()
    
    BT_group = UnitGroup('BT_group',
                                   units=(BT,))
    process_groups.append(BT_group)
    
    CT_group = UnitGroup('CT_group',
                                   units=(CT,))
    process_groups.append(CT_group)
    
    facilities_no_hu_group = UnitGroup('facilities_no_hu_group',
                                   units=(T601, T602, T603, T604, T605, 
                                          T606, T606_P, T607, u.CIP901,
                                          # u.ADP902, 
                                          # u.FWT903, 
                                          # u.PWC904,
                                          ))
    process_groups.append(facilities_no_hu_group)

    globals().update({'process_groups': process_groups})
    
# %% System setup

HP_sys = create_HP_sys()
# HP_sys.subsystems[-1].relative_molar_tolerance = 0.005
HP_sys.set_tolerance(mol=1e-2, rmol=1e-3, subsystems=True)

u = flowsheet.unit
s = flowsheet.stream
feedstock = s.feedstock
AA = s.AcrylicAcid
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185

feeds = HP_sys.feeds

products = [AA] # Don't include gypsum since we want to include carbon impurities in GWP calculation

emissions_from_sys = [i for i in flowsheet.stream
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


# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

HP_tea = CellulosicEthanolTEA(system=HP_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, 
                    # u.WWT_cost,
                    u.BT701, u.CT901, u.CWP901,  u.CIP901, u.ADP901, u.FWT901, u.PWC901,
                    ),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)


HP_no_BT_tea = HP_tea
# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================
BT = u.BT701
CT = u.CT901
CWP = u.CWP901
HXN = u.HXN1001


num_sims = 3
num_solve_tea = 3

try:
    HP_sys.simulate()
except:
    try:
        S404.add_specification(S404_spec)
        S404.specifications[0]()
        S404.simulate()
        HP_sys.simulate()
    except:
        S404.add_specification(S404_spec)
        S404.specifications[0]()
        S404.simulate()
        HP_sys.simulate()

def get_AA_MPSP():
    for i in range(num_sims):
        HP_sys.simulate()
    for i in range(num_solve_tea):
        AA.price = HP_tea.solve_price(AA)
    return AA.price

# try:
#     get_AA_MPSP()
# except:
#     # S404.add_specification(S404_spec)
#     # S404.specifications[0]()
#     # S404.simulate()
#     # get_AA_MPSP()
#     pass
# a_1
seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))

# yearly_production = 125000 # ton/yr; baseline
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
    spec_1=0.49,
    spec_2=54.8,
    spec_3=0.76,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = u.R401,
    byproduct_streams = [],
    HXN = HXN,
    maximum_inhibitor_concentration = 1.,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301], # if the line below does not work (depends on BioSTEAM version)
    pre_conversion_units = HP_sys.split(u.F301.ins[0])[0],
    baseline_titer = 54.8,
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = u.R201)

# Load baseline specifications

spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

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
    MPSP = AA.price = HP_tea.solve_price(AA)
    return MPSP

# Unit groups


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
total_C_out = AA.get_atomic_flow('C') + sum([emission.get_atomic_flow('C') for emission in emissions_from_sys])
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
get_emissions_GWP = lambda: sum([stream.get_atomic_flow('C') for stream in emissions_from_sys]) * HP_chemicals.CO2.MW / AA.F_mass
# GWP from electricity
get_net_electricity = lambda: sum(i.power_utility.rate for i in HP_sys.units)
get_net_electricity_GWP = lambda: get_net_electricity()*CFs['GWP_CFs']['Electricity'] \
    / AA.F_mass


get_total_electricity_demand = get_electricity_use = lambda: -BT.power_utility.rate
get_cooling_electricity_demand = lambda: u.CT901.power_utility.rate + CWP.power_utility.rate

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


HP_lca = LCA(HP_sys, HP_chemicals, CFs, feedstock, feedstock_ID, AA, BT, CT, CWP)

# %% Full analysis
# p11, p22, p33 = get_AA_MPSP(), HP_lca.GWP, HP_lca.FEC

def simulate_and_print():
    MPSP, GWP, FEC = get_AA_MPSP(), HP_lca.GWP, HP_lca.FEC
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    # print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg AA')
    print(f'GWP is {GWP:.3f} kg CO2-eq/kg AA')
    # print(f'Non-bio GWP is {get_ng_GWP():.3f} kg CO2-eq/kg AA')
    # print(f'FEC is {get_FEC():.2f} MJ/kg AA or {get_FEC()/AA_LHV:.2f} MJ/MJ AA\n')
    print(f'FEC is {FEC:.2f} MJ/kg AA or {FEC/AA_LHV:.2f} MJ/MJ AA')
    # print(f'SPED is {get_SPED():.2f} MJ/kg AA or {get_SPED()/AA_LHV:.2f} MJ/MJ AA')
    print('----------------------------------------\n')

# simulate_and_print()

def load_simulate_print(HP_yield=0.49, HP_titer=54.8, HP_productivity=0.76):
    spec.load_specifications(spec_1=HP_yield, spec_2=HP_titer, spec_3=HP_productivity)
    simulate_and_print()
    
# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
HP_sys.diagram('cluster')

#%%

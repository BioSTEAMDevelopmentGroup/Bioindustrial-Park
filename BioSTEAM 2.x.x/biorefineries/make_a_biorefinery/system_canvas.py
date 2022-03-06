#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020
@author: sarangbhagwat

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

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
from biosteam.process_tools import BoundedNumericalSpecification
from biosteam import System
from thermosteam import Stream
from biosteam.process_tools import UnitGroup
from biosteam import SystemFactory

from biorefineries.make_a_biorefinery import units, facilities
from biorefineries.make_a_biorefinery.process_settings import price
from biorefineries.make_a_biorefinery.utils import find_split, splits_df, baseline_feedflow
from biorefineries.make_a_biorefinery.chemicals_data import chems, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.make_a_biorefinery.tea import TemplateTEA

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

bst.speed_up()
flowsheet = bst.Flowsheet('canvas')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7

# Set default thermo object for the system
tmo.settings.set_thermo(chems)

feedstock_ID = 'Corn stover'

# System settings
System.default_converge_method = 'wegstein'
# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'

System.default_maxiter = 100
System.default_molar_tolerance = 0.1
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.strict_convergence = True # True => throw exception if system does not converge; False => continue with unconverged system

# %% System
@SystemFactory(ID = 'canvas_sys')
def create_canvas_sys(ins, outs):
    def fix_split(isplit, ID):
        isplit['SuccinicAcid', 'LacticAcid', 'Ethanol'] = isplit[ID]
       
    process_groups = []
    
    # %% 
    
    # =============================================================================
    # Feedstock streams
    # =============================================================================
    
    # !!! Replace the default corn stover feedstock stream 
    # with one that has a suitable composition, total flow, and price for your analysis
    
    # feedstock = Stream('feedstock', Water=250., 
    #                    SuccinicAcid=10., LacticAcid=5., Ethanol=5.,
    #                    FermMicrobe=8.,Lignin = 5., 
    #                    SolubleLignin = 10., GlucoseOligomer = 5.,
    #                    units='kg/hr',
    #                    price=price[feedstock_ID])
    
    feedstock = Stream('feedstock',
                        baseline_feedflow.copy(),
                        units='kg/hr',
                        price=price[feedstock_ID])
    
    # =============================================================================
    # Feedstock units
    # =============================================================================

    U101 = units.FeedstockPreprocessing('U101', ins=feedstock, outs='milled_feedstock')
    
    # IFF handling costs/utilities are included in feedstock cost and thus not considered here
    U101.cost_items['System'].cost = 0.
    U101.cost_items['System'].kW = 0.
    
    
    # !!! Replace or remove these fake HX units
    
    H101 = bst.units.HXutility('H101', ins=U101-0,
                                     outs='heated_feedstock',
                                     T=80.+273.15, rigorous=False)
    H101.line = 'Fake HX'
    
    H102 = bst.units.HXutility('H102', ins=H101-0,
                                     outs='cooled_feedstock',
                                     T=30.+273.15, rigorous=False)
    H102.line = 'Fake HX'
    
    # %% 
    
    # =============================================================================
    # Pretreatment streams
    # =============================================================================
    # !!! Any pretreatment input streams (output from storage) aside from the feedstock go here
    
    
    # =============================================================================
    # Pretreatment units
    # =============================================================================
    # !!! Replace these fake units with pretreatment process units
    R201 = units.FakePretreatment('R201', ins=H102-0, outs=('pretreated_slurry'))
    
    # %% 
    
    # =============================================================================
    # Conversion streams
    # =============================================================================
    # !!! Any conversion input streams (output from storage) aside from the feedstock go here
    saccharification_enzyme = tmo.Stream('saccharification_enzyme')
    saccharification_water = tmo.Stream('saccharification_water')
    
    # =============================================================================
    # Conversion units
    # =============================================================================
    # !!! Replace these fake units with conversion process units
    R301 = units.FakeSaccharification('R301', ins=(R201-0, saccharification_water, saccharification_enzyme),
                                      outs=('saccharified_slurry'))
    
    R302 = units.FakeCoFermentation('R302', ins=R301-0, outs=('fermentation_broth', 'vented_CO2'))
    
     # %% 
    
    # =============================================================================
    # Separation streams
    # =============================================================================
    
    # !!! Any separation input streams (output from storage) aside from the feedstock or fermentation broth go here
    
    
    # =============================================================================
    # Separation units
    # =============================================================================
    
    # !!! Replace this fake unit with separation process units
    
    bst.FakeSplitter._N_outs = 4
    S401 = bst.FakeSplitter('S401', ins=R302-0, outs=('s_product', 's_byproduct_1', 's_byproduct_2', 'S401_to_boiler', 'S401_to_WWT'))
    S401.line = 'Fake separation process'
    S401._graphics = bst.Unit._graphics
    def S401_spec():
        S401_ins_0 = S401.ins[0].copy()

        product_mol = S401_ins_0.imol['SuccinicAcid']
        S401_ins_0.imol['SuccinicAcid'] -= product_mol
        S401.outs[0].imol['SuccinicAcid'] = product_mol
        
        sw_mol = S401_ins_0.imol['FermMicrobe', 'Lignin']
        S401_ins_0.imol['FermMicrobe', 'Lignin'] -= sw_mol
        S401.outs[3].imol['FermMicrobe', 'Lignin'] = sw_mol
        
        retained_water_mol = S401_ins_0.imol['Water']*0.05 # arbitrary
        S401_ins_0.imol['Water'] -= retained_water_mol
        S401.outs[3].imol['Water'] = retained_water_mol
        
        byproduct_1_mol = S401_ins_0.imol['LacticAcid']
        S401_ins_0.imol['LacticAcid'] -= byproduct_1_mol
        S401.outs[1].imol['LacticAcid'] = byproduct_1_mol
        
        byproduct_2_mol = S401_ins_0.imol['Ethanol']
        S401_ins_0.imol['Ethanol'] -= byproduct_2_mol
        S401.outs[2].imol['Ethanol'] = byproduct_2_mol
        
        S401.outs[4].mol[:] = S401_ins_0.mol[:]
    S401.specification = S401_spec
    
    # %% 
    
    # =============================================================================
    # Wastewater treatment streams
    # =============================================================================
    
    # For aerobic digestion, flow will be updated in AerobicDigestion
    air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')
    
    # To neutralize nitric acid formed by nitrification in aerobic digestion
    # flow will be updated in AerobicDigestion
    # The active chemical is modeled as NaOH, but the price is cheaper than that of NaOH
    aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20.+273.15, P=2.*101325.,
                              price=price['Caustics'])
    
    # =============================================================================
    # Wastewater treatment units
    # =============================================================================
    
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=(S401-4)) # !!! here, add all wastewater streams (waste streams that are mostly mostly liquid-phase)
    
    # This represents the total cost of wastewater treatment system
    WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)
    
    R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                    outs=('biogas', 'anaerobic_treated_water', 
                                          'anaerobic_sludge'),
                                    reactants=soluble_organics, # !!! Ensure 'reactants' has every soluble organic chemical present in the biorefinery
                                    split=find_split(splits_df.index,
                                                     splits_df['stream_611'],
                                                     splits_df['stream_612'],
                                                     chemical_groups),
                                    T=35+273.15)
    fix_split(R501.isplit, 'Glucose')
    

    get_flow_dry_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
    
    # Mix recycled stream and wastewater after R501
    M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
    R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
                                  outs=('aerobic_vent', 'aerobic_treated_water'),
                                  reactants=soluble_organics,
                                  ratio=get_flow_dry_tpd()/2205)
    
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
    M505 = bst.units.Mixer('M505', ins=(S401-3), # !!! here, add all solid waste streams
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
    
    # !!! All fresh streams (with prices) go here
    saccharification_enzyme_fresh = tmo.Stream('saccharification_enzyme_fresh', price=price['Enzyme'])
    
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    # !!! All product and byproduct streams (with prices) go here
    # product_stream 
    product_stream = Stream('product_stream', units='kg/hr', price=price['Succinic acid']) # set an arbitrary price as this will be solved for
    # byproduct 1 stream
    byproduct_1_stream = Stream('byproduct_1_stream', units='kg/hr', price=price['Lactic acid'])
    # byproduct 2 stream
    byproduct_2_stream = Stream('byproduct_2_stream', units='kg/hr', price=price['Ethanol'])
    # Acetoin product

    # Cooling tower chemicals
    cooling_tower_chems = Stream('cooling_tower_chems', price=price['Cooling tower chems'])
    
    # Boiler ash
    ash = Stream('ash', price=price['Ash disposal'])
    
    # 145 based on equipment M-910 (clean-in-place system) in Humbird et al.
    CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_dry_tpd()/2205, units='kg/hr')
    
    # 1372608 based on stream 950 in Humbird et al.
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                          N2=0.79*1372608*get_flow_dry_tpd()/2205,
                          O2=0.21*1372608*get_flow_dry_tpd()/2205)
    
    # 8021 based on stream 713 in Humbird et al.
    fire_water_in = Stream('fire_water_in', 
                           Water=8021*get_flow_dry_tpd()/2205, units='kg/hr')
    
    # =============================================================================
    # Facilities units
    # =============================================================================

    # 7-day storage time, similar to ethanol's in Humbird et al.   
    # Product storage
    T601 = bst.StorageTank('T601', ins=S401-0, tau=7.*24., V_wf=0.9, # !!! ins should be the product stream
                                         vessel_type='Floating roof',
                                         vessel_material='Stainless steel')
    T601.line = 'ProductStorageTank'
    T601_P = bst.Pump('T601_P', ins=T601-0, outs=product_stream)
    
    # Byproduct 1 storage
    T602 = bst.StorageTank('T602', ins=S401-1, tau=7.*24., V_wf=0.9, # !!! ins should be a byproduct, if any
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')   
    T602.line = 'ByProduct1StorageTank'
    T602_P = bst.Pump('T602_P', ins=T602-0, outs=byproduct_1_stream)
    
    # Byproduct 1 storage
    T603 = bst.StorageTank('T603', ins=S401-2, tau=7.*24., V_wf=0.9, # !!! ins should be a byproduct, if any
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')   
    T603.line = 'ByProduct2StorageTank'
    T603_P = bst.Pump('T603_P', ins=T603-0, outs=byproduct_2_stream)
    
    # Storage tanks for other process inputs (besides water)
    #!!! Replace these default tau values with better assumptions
    T604 = bst.StorageTank('T604', ins=saccharification_enzyme_fresh, tau=7.*24., V_wf=0.9,
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')   
    T604.line = 'EnzymeStorageTank'
    T604_P = bst.Pump('T604_P', ins=T604-0, outs=saccharification_enzyme)
    
    
    # Misc. facilities
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_dry_tpd()/2205)
    FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')
    
    # Heating/cooling/electrical utility-generating utilities
    CWP = facilities.CWP('CWP', ins='return_chilled_water',
                          outs='process_chilled_water')

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
    
    # All water used in the system other than to regenerate steam for heating (which is handled by BT)
    process_water_streams = (CIP.ins[-1], CT.ins[-1], R301.ins[1])
    
    PWC = facilities.PWC('PWC', ins=(system_makeup_water, S504-0),
                         process_water_streams=process_water_streams,
                         recycled_blowdown_streams=None,
                         outs=('process_water', 'discharged_water'))
    
    # Heat exchanger network
    HXN = bst.facilities.HeatExchangerNetwork('HXN', 
                                              ignored=[])
 
    def HXN_no_run_cost():
        HXN.heat_utilities = tuple()
        HXN._installed_cost = 0.
    
    # # To simulate without HXN, uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []
    
    HXN_group = UnitGroup('HXN_group', 
                                   units=(HXN,))
    process_groups.append(HXN_group)
    
    
    # AWM = AutoWasteManagement('AWM', wastewater_mixer=M501, boiler_solids_mixer=M505,
    #                           to_wastewater_mixer_ID_key=liquid_waste_stream_ID_key,
    #                           to_boiler_solids_mixer_ID_key=solid_waste_stream_ID_key)
    
    BT_group = UnitGroup('BT_group',
                                   units=(BT,))
    process_groups.append(BT_group)
    
    CT_group = UnitGroup('CT_group',
                                   units=(CT,))
    process_groups.append(CT_group)
    
    facilities_no_hu_group = UnitGroup('facilities_no_hu_group',
                                   units=(T601, T601_P, PWC, ADP, CIP))
    process_groups.append(facilities_no_hu_group)

    globals().update({'process_groups': process_groups})
    
    globals().update({'get_flow_dry_tpd': get_flow_dry_tpd})
    
# %% System setup

canvas_sys = create_canvas_sys()

u = flowsheet.unit
s = flowsheet.stream
feedstock = s.feedstock
product_stream = s.product_stream
byproduct_1_stream = s.byproduct_1_stream
byproduct_2_stream = s.byproduct_2_stream


feeds = canvas_sys.feeds

products = [product_stream, byproduct_1_stream, byproduct_2_stream] # Don't include gypsum since we want to include carbon impurities in GWP calculation

emissions = [i for i in flowsheet.stream
                            if i.source and not i.sink and not i in products]
    
BT = flowsheet('BT')
BT_sys = System('BT_sys', path=(BT,))


globals().update(flowsheet.to_dict())

process_groups_dict = {}
for i in range(len(process_groups)):
    group = process_groups[i]
    process_groups_dict[group.name] = group


# %%
# =============================================================================
# TEA
# =============================================================================

# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

template_tea = TemplateTEA(system=canvas_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWT_cost,
                   u.T601, u.T601_P,
                    u.CWP, u.CT, u.PWC, u.CIP, u.ADP, u.FWT, u.BT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_dry_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)

# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

num_sims = 3
num_solve_tea = 3
def get_product_stream_MPSP():
    for i in range(num_sims):
        canvas_sys.simulate()
    for i in range(num_solve_tea):
        product_stream.price = template_tea.solve_price(product_stream)
    return product_stream.price

# get_product_stream_MPSP()

def simulate_and_print():
    MPSP = get_product_stream_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    print('----------------------------------------\n')

simulate_and_print()

# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
canvas_sys.diagram('cluster')



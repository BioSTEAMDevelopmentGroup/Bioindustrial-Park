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
from numba import njit
from biosteam import main_flowsheet as F
from biosteam.process_tools import BoundedNumericalSpecification
from biosteam import System
from thermosteam import Stream
from biorefineries.HP import units, facilities
from biorefineries.HP.lca import LCA
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
# System.default_converge_method = 'wegstein'
feedstock_ID = 'Corn stover'

# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'
System.default_converge_method = 'wegstein'
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.default_molar_tolerance = 0.1
System.strict_convergence = True # True => throw exception if system does not converge; false => continue with unconverged system


@SystemFactory(ID = 'pretreated_hydrolysate_sys')
def create_pretreated_hydrolysate_sys(ins, outs):
    
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
        H_M201._run()
    H_M201.specification = H_M201_specification
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
        H_M202._run()
    H_M202.specification = H_M202_specification
    
    
    # Prepare sulfuric acid
    get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
    T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid,
                                          feedstock_dry_mass=get_feedstock_dry_mass())
    
    M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, H_M201-0))
    
    # Mix sulfuric acid and feedstock, adjust water loading for pretreatment
    M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0, H_M202-0, ''))
    
    # Mix feedstock/sulfuric acid mixture and steam
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
    M501 = bst.units.Mixer('M501', ins=(Stream('fake_stream1', Water=1., Lignin=0.01))) # without sugars recycle
    
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
    # fix_split(R501.isplit, 'Glucose')
    
    # In TRY analysis, we aren't looking at the implications of varying yield of byproducts on 
    # sugars, but solely that of 3-HP. 
    # The extremes of this assumption doesn not significantly affect results.
    R501.byproducts_combustion_rxns = ParallelRxn([
        Rxn('AceticAcid -> 3 CO2 + H2O + O2', 'AceticAcid', 1.-1e-6, correct_atomic_balance=True),
        Rxn('Glycerol -> 3 CO2 + H2O + O2', 'Glycerol', 1.-1e-6, correct_atomic_balance=True)])
    for i in R501.byproducts_combustion_rxns:
        i.istoichiometry['O2'] = 0.
    def R501_specification():
        R501.byproducts_combustion_rxns(R501.ins[0])
        R501._run()
    R501.specification = R501_specification # Comment this out for anything other than TRY analysis
    
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
    # fix_split(S501.isplit, 'Glucose')
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
    # fix_split(S503.isplit, 'Glucose')
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
    M505 = bst.units.Mixer('M505', ins=(S503-1, Stream('fake_stream2', Water=1., Lignin=0.01)), 
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
    ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
    CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
    lime_fresh = Stream('lime_fresh', price=price['Lime'])
    
    hexanol_fresh = Stream('hexanol_fresh', price=price['Hexanol'])
    TOA_fresh = Stream('TOA_fresh', price=price['TOA'])
    AQ336_fresh = Stream('AQ336_fresh', price=price['AQ336'])
    
   
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

    # pretreated_hydrolysate product
    pretreated_hydrolysate = Stream('AcrylicAcid', units='kg/hr', price=price['AA'])
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
    
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_tpd()/2205)
    
    
    FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')
    
    # M304_H uses chilled water, thus requiring CWP
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
    
    # All water used in the system, here only consider water usage,
    # if heating needed, then heeating duty required is considered in BT
    process_water_streams = (water_M201, water_M202, water_M203, water_M205, 
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
    
    # To simulate without HXN, uncomment the following 3 lines:
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
                                   units=(T601, T602, PWC, ADP, CIP))
    process_groups.append(facilities_no_hu_group)

    globals().update({'process_groups': process_groups})
    
# %% System setup

pretreated_hydrolysate_sys = create_pretreated_hydrolysate_sys()
pretreated_hydrolysate_sys.subsystems[-1].relative_molar_tolerance = 0.005

u = flowsheet.unit
s = flowsheet.stream
feedstock = s.feedstock
# pretreated_hydrolysate = s.AcrylicAcid
pretreated_hydrolysate = u.P201.outs[0] 
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185

feeds = pretreated_hydrolysate_sys.feeds

products = [pretreated_hydrolysate] # Don't include gypsum since we want to include carbon impurities in GWP calculation

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

HP_no_BT_sys = bst.System('HP_no_BT_sys', path = pretreated_hydrolysate_sys.path, facilities = tuple([i for i in pretreated_hydrolysate_sys.facilities if not i.ID=='BT']))


#!!! Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

pretreated_hydrolysate_tea = CellulosicEthanolTEA(system=pretreated_hydrolysate_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWT_cost,
                    u.T601, u.T602,
                    u.CWP, u.CT, u.PWC, u.CIP, u.ADP, u.FWT, u.BT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)

HP_no_BT_tea = pretreated_hydrolysate_tea
# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

num_sims = 6
num_solve_tea = 3
def get_pretreated_hydrolysate_MPSP():
    for i in range(num_sims):
        pretreated_hydrolysate_sys.simulate()
    for i in range(num_solve_tea):
        pretreated_hydrolysate.price = pretreated_hydrolysate_tea.solve_price(pretreated_hydrolysate)
    return pretreated_hydrolysate.price


HP_lca = LCA(pretreated_hydrolysate_sys, HP_chemicals, CFs, feedstock, feedstock_ID, pretreated_hydrolysate, [CT, CWP])

# %% Full analysis
def simulate_and_print():
    MPSP, GWP, FEC = get_pretreated_hydrolysate_MPSP(), HP_lca.GWP, HP_lca.FEC
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    # print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg pretreated_hydrolysate')
    print(f'GWP is {GWP:.3f} kg CO2-eq/kg pretreated_hydrolysate')
    # print(f'Non-bio GWP is {get_ng_GWP():.3f} kg CO2-eq/kg pretreated_hydrolysate')
    # print(f'FEC is {get_FEC():.2f} MJ/kg pretreated_hydrolysate or {get_FEC()/pretreated_hydrolysate_LHV:.2f} MJ/MJ pretreated_hydrolysate\n')
    print(f'FEC is {FEC:.2f} MJ/kg pretreated_hydrolysate\n')
    print('----------------------------------------\n')

simulate_and_print()


# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
pretreated_hydrolysate_sys.diagram('cluster')

# GREET 2020

# fossil-based acrylic acid
# FEC: = 114 MJ eq. / kg pretreated_hydrolysate
# GHG100 = 5.82 kgCO2 eq. / kg pretreated_hydrolysate (7.65 kg CO2 eq. including EOL degradation as CO2)

# bio-based acrylic acid (via 3-HP from glycerol from algal oil)
# FEC: = 40 MJ eq. / kg pretreated_hydrolysate
# GHG100 = 3.14 kgCO2 eq. / kg pretreated_hydrolysate


## Ecoinvent 3.7.1 ##

# acrylic acid production
# CED-f = 49.026 MJ eq. / kg pretreated_hydrolysate
# GWP-100a = 2.245 kg CO2 eq. / kg pretreated_hydrolysate (4.075 including product EOL degradation as CO2)
# price = 1.4508 in 2005$ per kg
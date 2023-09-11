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

from biorefineries.cellulosic import create_facilities
from biorefineries.sugarcane import create_juicing_system_up_to_clarification
from biorefineries.TAL.process_areas import create_TAL_fermentation_process,\
                        create_TAL_separation_solubility_exploit_process,\
                        create_TAL_separation_adsorption_process

from scipy.interpolate import interp2d

filterwarnings('ignore')
IQ_interpolation = flx.IQ_interpolation

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# from lactic.hx_network import HX_Network

# # Do this to be able to show more streams in a diagram
# bst.units.Mixer._graphics.edge_in *= 2
# bst.speed_up()
flowsheet = bst.Flowsheet('TAL')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(TAL_chemicals)

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
    
    # # =============================================================================
    # # Feedstock
    # # =============================================================================
    
    # # feedstock = Stream('feedstock',
    # #                     baseline_feedflow.copy(),
    # #                     units='kg/hr',
    # #                     price=price['Feedstock'])
    
    # feedstock = Stream('feedstock')
    # feedstock.imass['Glucose'] = 29000.
    # feedstock.imass['H2O'] = 500.
    # feedstock.price = price['Glucose']*feedstock.imass['Glucose']/feedstock.F_mass
    
    # feedstock.F_mass = 25802.9 # at the baseline, the amount of TAL produced would exactly satisfy the US demand for sorbic acid with a hypothetical 100% TAL->sorbic acid conversion.
    # U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # # Handling costs/utilities included in feedstock cost thus not considered here
    # U101.cost_items['System'].cost = 0
    # U101.cost_items['System'].kW = 0
    
    # Sugarcane juicing subprocess
    sugarcane_juicing_sys = create_juicing_system_up_to_clarification()
    u = sugarcane_juicing_sys.flowsheet.unit
    s = sugarcane_juicing_sys.flowsheet.stream
    feedstock = s.sugarcane
    feedstock.F_mass = 833333.3333 # produces TAL at ~20718.2 metric tonne / year at baseline
    
    sugarcane_juicing_sys.flowsheet.diagram('thorough')

    
    #%% Fermentation
    fermentation_sys = create_TAL_fermentation_process(ins=(u.C201-0),
                                                   )
    
    
    
    # %% Separation
    Ethanol_desorption = Stream('Ethanol_desorption', units='kg/hr')
    separation_sys = create_TAL_separation_adsorption_process(ins=(fermentation_sys-1,
                                                                   Ethanol_desorption))
    # %% 
    
    # =============================================================================
    # Wastewater treatment streams
    # =============================================================================
    

    # =============================================================================
    # Wastewater treatment units
    # =============================================================================
    
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=(
                                        u.F301-1,
                                        separation_sys-2,
                                        '', '', '', '',
                                        ))
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=M501-0,
        mockup=True,
        area=500,
    )
    
    # Mix solid wastes to boiler turbogenerator
    M510 = bst.units.Mixer('M510', ins=(separation_sys-1, 
                                        u.U202-0,
                                        ),
                            outs='wastes_to_boiler_turbogenerator')
    
    MX = bst.Mixer(400, ['', ''])
    
    
    # %% 
    
    # =============================================================================
    # Facilities streams
    # =============================================================================
    
    ethanol_fresh = Stream('ethanol_fresh',  price=price['Ethanol'])
    

    CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
    # lime_fresh = Stream('lime_fresh', price=price['Lime'])
    
    # S401_out1_F_mass = S401.outs[1].F_mass
    
    # if not (S401_out1_F_mass == 0):
    #     ethanol_fresh = Stream('ethanol_fresh', Ethanol = 0.24 * S401_out1_F_mass, units='kg/hr', price=price['Ethanol']) - M401.ins[3].imass['Ethanol']
    #     DPHP_fresh = Stream('DPHP_fresh', DPHP = 0.25 * S401_out1_F_mass, units='kg/hr', price=price['DPHP']) - M401.ins[3].imass['Dipotassium hydrogen phosphate']
        
    # else:
    # ethanol_fresh = Stream('ethanol_fresh', Ethanol = get_feedstock_dry_mass()*48*22.1/1000*0.93, units='kg/hr', price=price['Ethanol'])
    # DPHP_fresh = Stream('DPHP_fresh', DPHP = get_feedstock_dry_mass()*50*22.1/1000*0.93, units='kg/hr', price=price['DPHP'])
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    
    imbibition_water = Stream('imbibition_water', price=price['Makeup water'])
    rvf_wash_water = Stream('rvf_wash_water', price=price['Makeup water'])
    dilution_water = Stream('dilution_water', price=price['Makeup water'])
    # TAL stream
    # TAL = Stream('TAL', units='kg/hr', price=price['TAL'])
    # TAL product
    TAL = Stream('TAL', units='kg/hr', price=price['TAL'])
    


    
    # 145 based on equipment M-910 (clean-in-place system) in Humbird et al.
    # CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_tpd()/2205, units='kg/hr')
    
    # 1372608 based on stream 950 in Humbird et al.
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    # plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
    #                       N2=0.79*1372608*get_flow_tpd()/2205,
    #                       O2=0.21*1372608*get_flow_tpd()/2205)
    
    # 8021 based on stream 713 in Humbird et al.
    # fire_water_in = Stream('fire_water_in', 
    #                        Water=8021*get_flow_tpd()/2205, units='kg/hr')
    
    # =============================================================================
    # Facilities units
    # =============================================================================

    
    T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=1-fermentation_sys)
    T603.line = 'CSL storage tank'
    
    # Ethanol storage
    T604 = bst.units.StorageTank('T604', ins=ethanol_fresh)
    T604.line = 'Ethanol storage tank'
    T604_P = units.TALPump('T604_P', ins=T604-0, outs = Ethanol_desorption)

    
    # 7-day storage time, similar to ethanol's in Humbird et al.
    T620 = units.TALStorageTank('T620', ins=separation_sys-0, tau=7*24, V_wf=0.9,
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    
    
    
    T620.line = 'TALStorageTank'
    
    
    T620_P = units.TALPump('T620_P', ins=T620-0, outs=TAL)
    
    

    
    s = flowsheet.stream
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
        feedstock=s.sugarcane,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX-0,
        BT_area=700,
        area=900,
    )
    
    # Heat exchanger network
    HXN = bst.facilities.HeatExchangerNetwork('HXN1001',
                                               ignored=lambda:[
                                                       u.M304_H,
                                                        # u.H202,
                                                        # H401,
                                                        # H402,
                                                        # H403,
                                                        # H404,
                                                        # A401.heat_exchanger_drying,
                                                        # A401.heat_exchanger_regeneration,
                                                        u.F401,
                                                        ],
                                                # cache_network=True,
                                                force_ideal_thermo=False,
                                               )
    
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
    
    # # To simulate without HXN, simply uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []

# %% 

# =============================================================================
# Complete system
# =============================================================================

TAL_sys = create_TAL_sys()


f = bst.main_flowsheet
u = f.unit
s = f.stream

feedstock = s.sugarcane
TAL = s.TAL
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185



TEA_feeds = set([i for i in TAL_sys.feeds if i.price]+ \
    [i for i in TAL_sys.feeds if i.price])

TEA_products = set([i for i in TAL_sys.products if i.price]+ \
    [i for i in TAL_sys.products if i.price]+[TAL])
    
    
for ui in u:
    globals().update({ui.ID: ui})
    
u.BT701.ID = 'BT701'
u.CT901.ID = 'CT801'
u.CWP901.ID = 'CWP802'
u.CIP901.ID = 'CIP901'
u.ADP901.ID = 'ADP902'
u.FWT901.ID = 'FWT903'
u.PWC901.ID = 'PWC904'


BT = flowsheet('BT')

BT.natural_gas_price = 0.2527
BT.ins[4].price = price['Lime']

globals().update(flowsheet.to_dict())
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
get_flow_dry_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
TAL_tea = TALTEA(system=TAL_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, 
        # operating_days=0.9*365,
        operating_days = 200,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(
                    u.U501,
                    # u.T601, u.T602, 
                    # u.T601, u.T602, u.T603, u.T604,
                    # u.T606, u.T606_P,
                    u.BT701, u.CT801, u.CWP802, u.CIP901, u.ADP902, u.FWT903, u.PWC904,
                    ),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_dry_tpd()/2205,
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

# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
TAL_sys.diagram('cluster')

#%% Define unit groups

area_names = [
    'feedstock',
    'fermentation',
    'separation & upgrading',
    'wastewater',
    'storage',
    'boiler & turbogenerator',
    'cooling utility facilities',
    'other facilities',
    'heat exchanger network',
]


unit_groups = bst.UnitGroup.group_by_area(TAL_sys.units)


unit_groups.append(bst.UnitGroup('natural gas (for steam generation)'))

unit_groups.append(bst.UnitGroup('fixed operating costs'))


for i, j in zip(unit_groups, area_names): i.name = j
for i in unit_groups: i.autofill_metrics(shorthand=False, 
                                         electricity_production=False, 
                                         electricity_consumption=True,
                                         material_cost=True)
for i in unit_groups:
    if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling utility facilities':
        i.metrics[-1].getter = lambda: 0. # Material cost
    if i.name == 'cooling utility facilities':
        i.metrics[1].getter = lambda: 0. # Cooling duty
    if i.name == 'boiler & turbogenerator':
        i.metrics[-1] = bst.evaluation.Metric('Material cost', 
                                            getter=lambda: TAL_tea.utility_cost/TAL_tea.operating_hours, 
                                            units='USD/hr',
                                            element=None) # Material cost
        # i.metrics[-2].getter = lambda: BT.power_utility.rate/1e3 # Electricity consumption [MW]
        
for HXN_group in unit_groups:
    if HXN_group.name == 'heat exchanger network':
        HXN_group.filter_savings = False
        assert isinstance(HXN_group.units[0], bst.HeatExchangerNetwork)


unit_groups[-2].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
                                                    units='USD/hr',
                                                    element=None)

# unit_groups[-2].metrics[-1] = bst.evaluation.Metric('Material cost', 
#                                                     getter=lambda: F404.ins[2].cost, 
#                                                     units='USD/hr',
#                                                     element=None)

unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: TAL_tea.FOC/TAL_tea.operating_hours, 
                                                    units='USD/hr',
                                                    element=None)

# for i in unit_groups:
#     i.metric(i.get_net_electricity_production,
#                 'Net electricity production',
#                 'kW')




unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i

BT = u.BT701
CT = u.CT801
CWP = u.CWP802
HXN = u.HXN1001


# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================


def get_TAL_MPSP():
    for i in range(3):
        TAL_sys.simulate()
    for i in range(3):
        TAL.price = TAL_tea.solve_price(TAL)
    return TAL.price*TAL.F_mass/TAL.imass['TAL']

spec = ProcessSpecification(
    evaporator = u.F301,
    pump = None,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = [],
    seed_train = u.R303,
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('TAL',),
    
    spec_1=0.19,
    spec_2=28.,
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
    baseline_titer = 23.,
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
    # M304, R302 = u.M304, u.R302
    M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
    M304.specifications[0]()
    M304_H._run()
    S302._run()
    R303._run()
    T301._run()
    # R302.simulate()
    R302.specifications[0]()
    return R302.effluent_titer - R302.titer_to_load

def F301_titer_obj_fn(V):
    F301.V = V
    F301._run()
    H301._run()
    M304.specifications[0]()
    M304_H._run()
    S302._run()
    R303._run()
    T301._run()
    # R302.simulate()
    R302.specifications[0]()
    return R302.effluent_titer - R302.titer_to_load

def load_titer_with_glucose(titer_to_load):
    spec.spec_2 = titer_to_load
    R302.titer_to_load = titer_to_load
    F301_titer_obj_fn(1e-4)
    if M304_titer_obj_fn(1e-4)<0: # if the slightest dilution results in too low a conc
        IQ_interpolation(F301_titer_obj_fn, 1e-4, 0.8, ytol=1e-4)
    # elif F301_titer_obj_fn(1e-4)>0: # if the slightest evaporation results in too high a conc
    else:
        F301_titer_obj_fn(1e-4)
        IQ_interpolation(M304_titer_obj_fn, 1e-4, 10000., ytol=1e-4) # for titers lower than 5 g/L with high yields, use an upper bound of 20000
    # else:
    #     raise RuntimeError('Unhandled load_titer case.')
    spec.titer_inhibitor_specification.check_sugar_concentration()
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
#     MPSP = TAL.price = TAL_tea.solve_price(TAL, TAL_no_BT_tea)
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
    get_TAL_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${get_TAL_MPSP():.3f}/kg')
    # print(f'GWP is {get_GWP():.3f} kg CO2-eq/kg TAL')
    # print(f'FEC is {get_FEC():.2f} MJ/kg TAL or {get_FEC()/TAL_LHV:.2f} MJ/MJ TAL')
    # print(f'SPED is {get_SPED():.2f} MJ/kg TAL or {get_SPED()/TAL_LHV:.2f} MJ/MJ TAL')
    # print('--------------------\n')

# simulate_and_print()
# TAL_sys.simulate()
get_TAL_MPSP()
# u.A401.cycle_time = 4.
# u.A401.drying_time = 0.5 #!!! Drying time is updated to this value (overwritten the value passed during initialization)
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
    # 'separation_sys': (S401, M401, M401_P,
    #                     S402, 
    #                     # F401, F401_H, X401,
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
#                        M504, S503, S504, M510),
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

# market price range: $6.51 - $6.74 /kg

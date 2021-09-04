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
from biorefineries import BST222
from biorefineries.cornstover import CellulosicEthanolTEA
from biosteam import System
from thermosteam import Stream
from biorefineries import BDO as bdo
from biorefineries.BDO import units, facilities
from biorefineries.BDO._process_specification import ProcessSpecification
from biorefineries.BDO.process_settings import price, CFs
from biorefineries.BDO.utils import find_split, splits_df, baseline_feedflow
from biorefineries.BDO.chemicals_data import BDO_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.BDO.tea import BDOTEA
from biorefineries.HP.lca import LCA

bst.speed_up()
flowsheet = bst.Flowsheet('BDO')
bst.main_flowsheet.set_flowsheet(flowsheet)
bst.System.default_relative_molar_tolerance = 1e-4
bst.System.default_maxiter = 100
bst.System.default_molar_tolerance = 0.1

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(BDO_chemicals)

BDO_sys = bdo.create_system_oleyl_alcohol()
u = flowsheet.unit
feedstock = BDO_sys.ins[0]
# feedstock.mol[:] = baseline_feedflow[:]
MEK, isobutanol = BDO_sys.outs
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
u = F.unit
s = F.stream
# globals().update(F.to_dict())

# %%
# =============================================================================
# TEA
# =============================================================================

# BDO_tea = bst.CombinedTEA([BDO_no_BT_tea, BT_tea], IRR=0.10)
BDO_tea = CellulosicEthanolTEA(system=BDO_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=bst.get_OSBL(BDO_sys.units),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

# sub_units = BDO_sys.units[:BDO_sys.units.index(u.R401)]
# sub_sys = bst.System(sub_units)
# sub_tea = CellulosicEthanolTEA(system=sub_sys, IRR=0.10, duration=(2016, 2046),
#         depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
#         lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
#         startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
#         startup_VOCfrac=0.75, WC_over_FCI=0.05,
#         finance_interest=0.08, finance_years=10, finance_fraction=0.4,
#         # biosteam Splitters and Mixers have no cost, 
#         # cost of all wastewater treatment units are included in WWT_cost,
#         # BT is not included in this TEA
#         OSBL_units=bst.get_OSBL(sub_sys.units),
#         warehouse=0.04, site_development=0.09, additional_piping=0.045,
#         proratable_costs=0.10, field_expenses=0.10, construction=0.20,
#         contingency=0.10, other_indirect_costs=0.10, 
#         labor_cost=3212962*get_flow_tpd()/2205,
#         labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
#         steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT)

area_names = [
    'feedstock',
    'pretreatment',
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
unit_groups = bst.UnitGroup.group_by_area(BDO_sys.units)
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

TEA_feeds = set([i for i in BDO_sys.feeds if i.price])

TEA_products = set([i for i in BDO_sys.products if i.price])
MEK.price = 0.
def get_MEK_MPSP():
    BDO_sys.simulate()
    
    for i in range(3):
        MEK.price = BDO_tea.solve_price(MEK)
    return MEK.price

# get_MEK_MPSP()

# R301 = F('R301') # Fermentor
seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))

# yearly_production = 125000 # ton/yr
spec = ProcessSpecification(
    evaporator = u.F301,
    evaporator_pump = u.F301_P,
    pump = u.M304_H_P,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = seed_train_system,
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('BDO',),
    spec_1=0.36,
    spec_2=109.9,
    spec_3=1.,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = u.R401,
    byproduct_streams = [],
    HXN = u.HXN1001,
    tolerable_HXN_energy_balance_percent_error=5,
    maximum_inhibitor_concentration = 1.,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301], # if the line below does not work (depends on BioSTEAM version)
    # pre_conversion_units = BDO_sys.split(u.F301.ins[0])[0],
    baseline_titer = 54.8,
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = u.R201)



spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity


# spec = ProcessSpecification(
#     evaporator = u.F301,
#     mixer = u.M304,
#     reactor=u.R302,
#     reaction_name='fermentation_reaction',
#     substrates=('Xylose', 'Glucose'),
#     products=('BDO',),
#     spec_1=0.8,
#     spec_2=109.9,
#     spec_3=1.,
#     path = (u.M304_H, u.M304_H_P),
#     xylose_utilization_fraction = 0.80,
#     feedstock = feedstock,
#     dehydration_reactor = u.R401,
#     byproduct_streams = [s.isobutanol],
#     evaporator_pump = u.F301_P)

path = (u.F301, u.R302)
@np.vectorize
def calculate_titer(V):
    u.F301.V = V
    for i in path: i._run()
    return spec._calculate_titer()

@np.vectorize   
def calculate_MPSP(V):
    u.F301.V = V
    BDO_sys.simulate()
    MPSP = MEK.price = BDO_tea.solve_price(MEK)
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

get_Isobutanol_GWP = lambda: - isobutanol.F_mass * CFs['GWP_CFs']['Isobutanol']/MEK.F_mass
get_Isobutanol_FEC = lambda: - isobutanol.F_mass * CFs['FEC_CFs']['Isobutanol']/MEK.F_mass

def get_material_GWP():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
    # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
    return chemical_GWP.sum()/MEK.F_mass

# GWP from combustion of non-biogenic carbons
get_non_bio_GWP = lambda: (s.natural_gas.get_atomic_flow('C') 
                           + s.oleyl_alcohol.get_atomic_flow('C')) \
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
get_FEC = lambda: get_material_FEC() + get_electricity_FEC() - get_Isobutanol_FEC()

get_SPED = lambda: u.BT.system_heating_demand*0.001/MEK.F_mass
MEK_LHV = 31.45 # MJ/kg MEK

BDO_lca = LCA(BDO_sys, BDO_chemicals, CFs, feedstock, 'Corn stover', MEK, [CT, CWP], BT=BT, CT=CT)

#%% TEA breakdown
def TEA_breakdown():
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = metric() + unit_groups[5].metrics[ug.metrics.index(metric)]()
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()
            # else:
            #     storage_metric_val = metric()
                
    # return metric_breakdowns
    for i in unit_groups[0].metrics:
        print(f"\n\n----- {i.name} ({i.units}) -----")
        metric_breakdowns_i = metric_breakdowns[i.name]
        for j in metric_breakdowns_i.keys():
            print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")

# %% LCA breakdown

def FEC_breakdown():
    return {
        'feedstock': BDO_lca.feedstock_FEC,
        'heating demand': BDO_lca.heating_demand_FEC,
        'cooling demand': BDO_lca.cooling_demand_FEC,
        'electricity demand (non-cooling)': BDO_lca.electricity_demand_non_cooling_FEC,
        'other materials': BDO_lca.material_FEC,
    }

def GWP_breakdown():
    return {
        'feedstock': BDO_lca.FGHTP_GWP,
        'heating demand': BDO_lca.heating_demand_GWP,
        'cooling demand': BDO_lca.cooling_demand_GWP,
        'electricity demand (non-cooling)': BDO_lca.electricity_demand_non_cooling_GWP,
        'other direct non-biogenic emissions': BDO_lca.non_BT_direct_emissions_GWP,
        'other materials': BDO_lca.material_GWP,
    }


# %% Full analysis
get_MEK_MPSP()

# Overall yield assuming yield on xylose is 80% of yield on glucose:
# (0.64*0.36+0.36*0.8*0.36) * 2 = 0.33408 * 2
# saccharified stream mass ratio of glucose:xylose is 0.64:0.36


# to achieve $1.98/kg: spec_1=0.353*2
spec.load_specifications(0.33408*2, 109.9, 1.0)

def print_recycles():
    sys = BDO_sys.copy()
    sys.flatten()
    for i in sys.recycle:
        print(i, i.imass['OleylAlcohol'] / 1e3)

def print_capital_cost():
    print('TCI', BDO_tea.TCI / 1e6, 'MM$')
    print('Fermentation installed cost', u.R302.installed_cost / 1e6, 'MM$')
    print(u.R302.results())

def print_specs():
    print(spec.spec_1, spec.spec_2, spec.spec_3)

def simulate_and_print():
    get_MEK_MPSP()
    print('\n---------- Simulation Results ----------')
    # print_recycles()
    # print_specs()
    # print_capital_cost()
    print(f'MPSP is ${get_MEK_MPSP():.3f}/kg')
    print(f'GWP is {BDO_lca.GWP:.3f} kg CO2-eq/kg MEK')
    # print(f'Non-bio GWP is {():.3f} kg CO2-eq/kg MEK')
    print(f'FEC is {BDO_lca.FEC:.2f} MJ/kg MEK or {get_FEC()/MEK_LHV:.2f} MJ/MJ MEK')
    # print(f'SPED is {get_SPED():.2f} MJ/kg MEK or {get_SPED()/MEK_LHV:.2f} MJ/MJ MEK')
    print('--------------------\n')

simulate_and_print()
# spec.load_specifications(0.36, 13.7, 1.0)
# simulate_and_print()
# spec.load_specifications(0.36, 109.9, 1.0)
# simulate_and_print()

# %% Maximum BDO titer (given by maximum sugar concentration)



# %%% MEK environmental impacts

# GREET 2020:
#     GHG100 = 1.04 kgCO2-eq/kg
#     FEC = 43 MJ/kg
# Ecoinvent 3.7.1:
#     IPCC 2013 GWP100 = 1.9166 kgCO2-eq/kg
#     CED-f = 58.707 MJ/kg

# %%

# =============================================================================
# For Monte Carlo and analyses
# =============================================================================

# BDO_sub_sys = {
#     'feedstock_sys': (U101,),
#     'pretreatment_sys': (T201, M201, M202, M203, 
#                          R201, R201_H, T202, T203,
#                          F201, F201_H,
#                          M204, T204, T204_P,
#                          M205, M205_P),
#     'conversion_sys': (H301, M301, M302, R301, R302, T301),
    # 'separation_sys': (u.S401, u.M401, u.M401_P,
    #                     u.S402, 
    #                     # F401, F401_H, F401_P,
    #                     u.D401, u.D401_H, u.D401_P, u.S403,
    #                     u.M402_P, u.S403,
    #                     u.D403, u.D403_H, u.D403_P,
    #                     u.M501,
    #                     u.T606, u.T606_P)
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
    # }





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
from biorefineries.HP.process_areas import create_HP_fermentation_process,\
                                           create_HP_separation_improved_process,\
                                           create_HP_to_acrylic_acid_upgrading_process
from biorefineries.HP.lca import HPLCA
from biorefineries.HP.models._process_specification import ProcessSpecification
from biorefineries.HP.process_settings import price, CFs, chem_index, _GDP_2007_to_2010
from biorefineries.HP.utils import find_split, splits_df, baseline_feedflow
from biorefineries.HP.chemicals_data import HP_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.HP.tea import HPTEA
from biosteam.process_tools import UnitGroup
import matplotlib.pyplot as plt
import copy
from biorefineries.cornstover import CellulosicEthanolTEA as HPTEA
from biosteam import SystemFactory
from biorefineries.cellulosic import create_facilities
from biorefineries.cellulosic.units import Saccharification

from biorefineries.cornstover import create_dilute_acid_pretreatment_system, create_saccharification_system
# from biorefineries import corn
# from lactic.hx_network import HX_Network

from biorefineries.TAL._general_utils import call_all_specifications_or_run,\
                                                get_more_unit_groups,\
                                                add_metrics_to_unit_groups,\
                                                set_production_capacity,\
                                                TEA_breakdown,\
                                                update_facility_IDs
                                                
import autosynthesis

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
bst.CE = bst.units.design_tools.CEPCI_by_year[2019]
# _labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(HP_chemicals)

System.default_maxiter = 100
# System.default_converge_method = 'wegstein'
feedstock_ID = 'Corn stover'
# feedstock_ID = 'Corn'

# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'
# System.default_converge_method = 'wegstein'
# System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.default_relative_molar_tolerance = 0.001 # supersedes absolute tolerance
System.default_molar_tolerance = 0.1
System.strict_convergence = True # True => throw exception if system does not converge; false => continue with unconverged system

# %% System setup

HP_sys = autosynthesis.process_blocks.get_system_block_based(feedstock='corn stover', 
                                                             product='glacial acrylic acid',
                                                             choice=0)

# HP_sys.subsystems[-1].relative_molar_tolerance = 0.005
HP_sys.set_tolerance(mol=1e-3, rmol=1e-3, subsystems=True)
HP_sys.simulate()

f = bst.main_flowsheet
u = f.unit
s = f.stream

feedstock = s.cornstover
AA = s.product_glacial_AA
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185

update_facility_IDs(HP_sys)
feeds = HP_sys.feeds

products = [AA] # Don't include gypsum since we want to include carbon impurities in GWP calculation

emissions_from_sys = [i for i in flowsheet.stream
                            if i.source and not i.sink and not i in products]
    
BT = u.BT701
CT = u.CT801
CWP = u.CWP802
CWP2 = u.CWP803
HXN = u.HXN1001


globals().update(flowsheet.to_dict())

# %%
# =============================================================================
# TEA
# =============================================================================

# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

HP_tea = HPTEA(system=HP_sys, IRR=0.10, duration=(2019, 2049),
        depreciation='MACRS7', income_tax=0.21, operating_days=330.,
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
                    # u.T620, u.T620_P,
                    u.BT701, u.CT801, u.CWP802, u.CWP803, u.CIP901, u.ADP902, u.FWT903, u.PWC904,
                    ),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

HP_no_BT_tea = HP_tea

#%%
# =============================================================================
# LCA
# =============================================================================

HP_lca = HPLCA(system=HP_sys, 
                 CFs=CFs, 
                 feedstock=feedstock, 
                 feedstock_ID=feedstock_ID,
                 main_product=AA, 
                 main_product_chemical_IDs=['AcrylicAcid',], 
                 by_products=[s.gypsum], 
                 cooling_tower=u.CT801, 
                 chilled_water_processing_units=[u.CWP802, u.CWP803], 
                 boiler=u.BT701, has_turbogenerator=True,
                 add_EOL_GWP=True, # (True for cradle-to-grave LCA; False for cradle-to-gate)
                 input_biogenic_carbon_streams=(feedstock, s.CSL), # any streams, including feedstock, for which CO2 fixing credit is non-zero and has not already been included in GWP impact CFs
                 )

#%% Define unit groups and their metrics

feedstock_acquisition_group = bst.UnitGroup('feedstock acquisition', units=[])
feedstock_juicing_group = bst.UnitGroup('feedstock pretreatment and saccharification', units=f.cornstover_pretreatment_sys.units+f.cornstover_saccharification_sys.units)
fermentation_group = f.HP_fermentation_process.to_unit_group('fermentation')
separation_group = f.HP_separation_improved_process.to_unit_group('separation')
upgrading_group = f.HP_to_acrylic_acid_upgrading_process.to_unit_group('upgrading')


unit_groups = [
    feedstock_acquisition_group,
    feedstock_juicing_group,
    fermentation_group,
    separation_group,
    upgrading_group,
    ]

unit_groups += get_more_unit_groups(system=HP_sys,
                         groups_to_get=['wastewater',
                                        'storage & other facilities',
                                        'boiler & turbogenerator',
                                        'cooling utility facilities',
                                        'other facilities',
                                        'heat exchanger network',
                                        'natural gas (for steam generation)',
                                        # 'natural gas (for product drying)',
                                        # 'chilled brine',
                                        'fixed operating cost',
                                        'electricity consumption',
                                        'heating duty',
                                        'excess electricity',
                                        ]
                         )

add_metrics_to_unit_groups(unit_groups=unit_groups, system=HP_sys, TEA=HP_tea, LCA=HP_lca, hxn_class=bst.HeatExchangerNetwork)

unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i

cooling_facilities_unit_group = unit_groups_dict['cooling utility facilities']

for i in cooling_facilities_unit_group.metrics:
    if i.name.lower() in ('electricity consumption', 'power consumption',):
        i.getter = lambda: sum([ui.power_utility.rate for ui in cooling_facilities_unit_group.units])/1e3

# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

HP_fermentation_process = f.fermentation_HP

def get_AA_MPSP():
    for i in range(3):
        HP_sys.simulate()
    for i in range(3):
        AA.price = HP_tea.solve_price(AA)
    return AA.price*AA.F_mass/AA.imass['AcrylicAcid']


theoretical_max_g_HP_per_g_glucose = 2*HP_chemicals.HP.MW/HP_chemicals.Glucose.MW



desired_annual_production = 134_000 # pure metric ton / y # satisfy 50% of 2019 US demand for acrylic acid

# desired_annual_production = (23_802) * kg_SA_to_kg_KSA # pure metric ton / y # satisfy 50% of 2019 US demand for acrylic acid


spec = ProcessSpecification(
    evaporator = u.F301,
    pump = None,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = [],
    seed_train = u.R303,
    reactor= u.R302,
    reaction_name='fermentation_reaction', # pure metric ton / y
    substrates=('Xylose', 'Glucose'),
    products=('HP',),
    
    desired_annual_production = desired_annual_production, 
    
    spec_1=0.721, 
    spec_2=96., 
    spec_3=0.5714, 

    
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = [],
    HXN = u.HXN1001,
    maximum_inhibitor_concentration = 1.,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301], # if the line below does not work (depends on BioSTEAM version)
    # pre_conversion_units = HP_sys.split(u.M304.ins[0])[0],
    pre_conversion_units = [],
    
    # (ranges from Cao et al. 2022)
    # baseline_yield = 0.0815/theoretical_max_g_HP_per_g_glucose, # mean of 0.074 and 0.089 g/g (range from Cao et al. 2022)
    # baseline_titer = 25.5, # mean of 23 and 28 g/L (range from Cao et al. 2022)
    # baseline_productivity = 0.215, # mean of 0.19 and 0.24 g/L/h (range from Cao et al. 2022)
    
    
    # !!! set baseline fermentation performance here
    baseline_yield = 0.73, 
    baseline_titer = 62.5, 
    baseline_productivity = 0.28, 
    
    
    tolerable_HXN_energy_balance_percent_error = 0.01,
    
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = None)


spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

def clear_units(units_to_clear):
    for i in units_to_clear:
        for j in list(i.ins)+list(i.outs):
            j.empty()
        i.simulate()
        

def M304_titer_obj_fn(water_to_sugar_mol_ratio):
    M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
    call_all_specifications_or_run(HP_fermentation_process.path)
    # HP_fermentation_process.run()
    return R302.effluent_titer - R302.titer_to_load

def F301_titer_obj_fn(V):
    F301.V = V
    call_all_specifications_or_run(HP_fermentation_process.path)
    # HP_fermentation_process.run()
    return R302.effluent_titer - R302.titer_to_load

def load_titer_with_glucose(titer_to_load, set_F301_V=0.8):
    # clear_units([V301, K301])
    F301_ub = 1.-1e-3
    F301_lb = 0. if set_F301_V is None else set_F301_V
    M304_lb, M304_ub = 0., 100_000  # for low-titer high-yield combinations, if infeasible, use a higher upper bound
    
    spec.spec_2 = titer_to_load
    R302.titer_to_load = titer_to_load
    F301_titer_obj_fn(F301_lb)
    
    if M304_titer_obj_fn(M304_lb) < 0.: # if there is too low a conc even with no dilution
        IQ_interpolation(F301_titer_obj_fn, F301_lb, F301_ub, ytol=1e-3)
    # elif F301_titer_obj_fn(1e-4)>0: # if the slightest evaporation results in too high a conc
    elif M304_titer_obj_fn(M304_ub) > 0.:
        IQ_interpolation(M304_titer_obj_fn, 
                         M304_lb,
                         M304_ub, 
                         ytol=1e-3)
    else:
        F301_titer_obj_fn(F301_lb)
        IQ_interpolation(M304_titer_obj_fn, 
                         M304_lb, 
                         M304_ub, 
                         ytol=1e-3)

    spec.titer_inhibitor_specification.check_sugar_concentration()
    
    
spec.load_spec_2 = load_titer_with_glucose

#%% Load specifications
# If, during TRY analysis, you'd like to set production capacity rather than feedstock capacity as constant, 
# use system.simulate rather than spec.set_production_capacity

def spec_set_production_capacity(
                        desired_annual_production=spec.desired_annual_production, # pure metric ton /y
                        method='analytical', # 'IQ_interpolation' or 'analytical'
                        system=HP_sys,
                        TEA=None,
                        spec=spec,
                        product_stream=AA, 
                        product_chemical_IDs=['AcrylicAcid',],
                        feedstock_stream=feedstock,
                        feedstock_F_mass_range=[5000, 2000_000], # wet-kg/h)
                        ):
    set_production_capacity(
                            desired_annual_production=desired_annual_production, # pure metric ton /y
                            method=method, # 'IQ_interpolation' or 'analytical'
                            system=system,
                            TEA=TEA,
                            spec=spec,
                            product_stream=product_stream, 
                            product_chemical_IDs=product_chemical_IDs,
                            feedstock_stream=feedstock_stream,
                            feedstock_F_mass_range=feedstock_F_mass_range, # wet-kg/h
                            )
    
spec.set_production_capacity = spec_set_production_capacity
    
# %% Simulate and evaluate baseline

production_capacity_is_fixed = True
def simulate_and_print():
    if production_capacity_is_fixed: spec_set_production_capacity(spec.desired_annual_production, method='analytical')
    # set_production_capacity(25000, 'analytical')
    print('\n---------- Simulation Results ----------')
    MPSP_AA = get_AA_MPSP()
    print(f'MPSP is ${MPSP_AA:.3f}/kg AA')
    GWP_AA, FEC_AA = HP_lca.GWP, HP_lca.FEC
    print(f'GWP-100a is {GWP_AA:.3f} kg CO2-eq/kg AA')
    print(f'FEC is {FEC_AA:.3f} MJ/kg AA')
    GWP_AA_without_electricity_credit, FEC_AA_without_electricity_credit =\
        GWP_AA - HP_lca.net_electricity_GWP, FEC_AA - HP_lca.net_electricity_FEC
    print(f'GWP-100a without electricity credit is {GWP_AA_without_electricity_credit:.3f} kg CO2-eq/kg AA')
    print(f'FEC without electricity credit is {FEC_AA_without_electricity_credit:.3f} MJ/kg AA')
    # print(f'FEC is {get_FEC():.2f} MJ/kg AA or {get_FEC()/HP_LHV:.2f} MJ/MJ HP')
    # print(f'SPED is {get_SPED():.2f} MJ/kg AA or {get_SPED()/HP_LHV:.2f} MJ/MJ HP')
    # print('--------------------\n')

# %% Initial simulation

for i in range(2): get_AA_MPSP()

spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)

for i in range(2):
    spec_set_production_capacity(
                            desired_annual_production=spec.desired_annual_production, # pure metric ton /y
                            )

HP_tea.labor_cost = 3212962*get_flow_tpd()/2205

#%% Misc.

def get_non_gaseous_waste_carbon_as_fraction_of_HP_GWP100():
    return sum([i.get_atomic_flow('C') for i in HP_sys.products if i.F_mol 
                and ('l' in i.phases or 's' in i.phases or i.phase=='l') 
                and (not i==AA)])/AA.imass['AcrylicAcid']/HP_lca.GWP

#%%

# simulate_and_print()

# %% Diagram

bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
HP_sys.diagram('cluster')

#%% TEA breakdown

# TEA_breakdown(unit_groups_dict=unit_groups_dict,
#               print_output=True,
#               )

#%% TEA breakdown figure
import contourplots

system = HP_sys

###### change operating cost unit labels $/h to MM$/y
for i in unit_groups:
    for j in i.metrics:
        if j.name == 'Operating cost':
            j.units = r"$\mathrm{MM\$}$" + '\u00b7y\u207b\u00b9'
######

df_TEA_breakdown = bst.UnitGroup.df_from_groups(
    unit_groups, fraction=True,
    scale_fractions_to_positive_values=True,
)

# totals=[sum([ui.metrics[i]() for ui in unit_groups])
#         for i in range(len(unit_groups[0].metrics))]

totals=[]
metrics = unit_groups[0].metrics
for i in range(len(metrics)):
    curr_total = 0.
    for ui in unit_groups:
        curr_total += ui.metrics[i]()
    if metrics[i].name=='Operating cost':
        # change total operating cost from $/h to MM$/y
        curr_total *= system.TEA.operating_hours/1e6
    totals.append(curr_total)




plot = False
if plot: 
    contourplots.stacked_bar_plot(dataframe=df_TEA_breakdown, 
                     y_ticks = [-40, -20, 0, 20, 40, 60, 80, 100],
                     y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                     y_units = "%", 
                     colors=['#7BBD84', 
                             '#E58835', 
                             '#F7C652', 
                             '#63C6CE', 
                             '#F8858A', 
                             '#94948C', 
                             '#734A8C', 
                             '#D1C0E1', 
                             '#648496', 
                             # '#B97A57', 
                             '#D1C0E1', 
                             # '#F8858A', 
                               '#b00000', 
                             # '#63C6CE', 
                             '#94948C', 
                             # '#7BBD84', 
                             '#b6fcd5', 
                             '#E58835', 
                             # '#648496',
                             '#b6fcd5',
                             ],
                     hatch_patterns=('\\', '//', '|', 'x',),
                     filename='AA_system_corn_improved_separations' + '_TEA_breakdown_stacked_bar_plot',
                     n_minor_ticks=4,
                     fig_height=5.5*1.1777*0.94*1.0975,
                     fig_width=10,
                     show_totals=True,
                     totals=totals,
                     sig_figs_for_totals=3,
                     units_list=[i.units for i in unit_groups[0].metrics],
                     totals_label_text=r"$\bfsum:$",
                     rotate_xticks=45.,
                     )

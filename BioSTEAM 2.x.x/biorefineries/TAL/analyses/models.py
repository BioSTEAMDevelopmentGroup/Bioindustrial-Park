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

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

# import numpy as np
import biosteam as bst
from chaospy import distributions as shape
# from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
# from biosteam.evaluation.evaluation_tools import Setter
from TAL.system_TAL_adsorption_glucose import TAL_sys, TAL_tea, u, s, unit_groups, spec

# get_annual_factor = lambda: TAL_tea._annual_factor


get_annual_factor = lambda: TAL_tea.operating_days*24 # hours per year

_kg_per_ton = 907.18474


system_feeds = [i for i in TAL_sys.feeds if i.price]
system_products = [i for i in TAL_sys.products if i.price]
    
# gypsum = find.stream.gypsum
# system_products.append(gypsum)

baseline_yield, baseline_titer, baseline_productivity =\
    spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.feedstock
TAL = s.SA
ethanol_fresh = s.ethanol_fresh
CSL_fresh = s.CSL_fresh


R302 = u.R302
R303 = u.R303

AC401 = u.AC401
F401 = u.F401
F402 = u.F402
H403 = u.H403

H402 = u.H402

T620 = u.T620

BT = u.BT701

_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of TAL stream
def get_MSP():
    for i in range(3):
        TAL.price = TAL_tea.solve_price(TAL,
                                    )
    return TAL.price

# Mass flow rate of TAL stream
get_yield = lambda: TAL.F_mass*get_annual_factor()/1e6
# Purity (%) of TAL in the final product
get_purity = lambda: TAL.imass['TAL']/TAL.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: TAL.imol['TAL'] \
    /(R302.outs[0].imol['TAL'])
get_overall_TCI = lambda: TAL_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: TAL_tea.AOC/1e6
get_material_cost = lambda: TAL_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: TAL_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (TAL_sys.power_utility.production-TAL_sys.power_utility.consumption)
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Product yield', get_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product purity', get_purity, '%', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg', 'Biorefinery'),
           Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product recovery', get_recovery, '%', 'Biorefinery'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $', 'Biorefinery'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery')
           ]


# To see if TEA converges well for each simulation
get_NPV = lambda: TAL_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))

# for ug in unit_groups:
#     ug_metrics = ug.metrics
#     metrics.append(Metric(ug.name, ug.metrics[0], '10^6 $', 'Installed cost'))
#     metrics.append(Metric(ug.name, ug.metrics[5], 'USD/h', 'Material cost'))
#     metrics.append(Metric(ug.name, ug.metrics[1], 'GJ/h', 'Cooling duty'))
#     metrics.append(Metric(ug.name, ug.metrics[2], 'GJ/h', 'Heating duty'))
#     metrics.append(Metric(ug.name, ug.metrics[3], 'MW', 'Electricity usage'))

metrics_labels_dict = {
    'Installed cost':(0, '10^6 $'), 
    'Material cost':(5,'USD/h'), 
    'Cooling duty':(1,'GJ/h'), 
    'Heating duty':(2,'GJ/h'), 
    'Electricity usage':(3, 'MW'), 
    }

for m, u_i in metrics_labels_dict.items():
    for ug in unit_groups:
        metrics.append(Metric(ug.name, ug.metrics[u_i[0]], u_i[1], m))
#%%
# =============================================================================
# Construct base model and add parameters
# =============================================================================

model = TAL_model = Model(TAL_sys, metrics)

param = model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

def baseline_triangle(baseline, ratio):
    return shape.Triangle(baseline*(1-ratio), baseline, baseline*(1+ratio))

# A fake parameter serving as a "blank" in sensitivity analysis to capture
# fluctuations due to converging errors
D = baseline_uniform(1, 0.1)
@param(name='Blank parameter', element=feedstock, kind='coupled', units='',
       baseline=1, distribution=D)
def set_blank_parameter(anything):
    # This does nothing
    feedstock.T = feedstock.T

#%% ######################## TEA parameters ########################
# U101 = SSCF.U101
# D = baseline_uniform(2205, 0.1)
# @param(name='Feedstock flow rate', element=feedstock, kind='coupled', units='dry-ton/day',
#        baseline=2205, distribution=D)
# def set_feedstock_flow_rate(rate):
#     feedstock.mass *= rate / U101._cached_flow_rate
#     U101._cached_flow_rate = rate

D = shape.Triangle(0.84, 0.9, 0.96)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.9, distribution=D)
def set_plant_uptime(uptime):
    TAL_tea.operating_days = 365. * uptime


# D = baseline_triangle(1, 0.25)
# @param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
#         baseline=1, distribution=D)
# def set_TCI_ratio(new_ratio):
#     old_ratio = TAL_tea._TCI_ratio_cached
#     # old_ratio = TAL_no_CHP_tea._TCI_ratio_cached
#     for unit in TAL_sys.units:
#         if hasattr(unit, 'cost_items'):
#             for item in unit.cost_items:
#                 unit.cost_items[item].cost /= old_ratio
#                 unit.cost_items[item].cost *= new_ratio
#     TAL_tea._TCI_ratio_cached = new_ratio
#     # TAL_no_CHP_tea._TCI_ratio_cached = new_ratio


# Only include materials that account for >5% of total annual material cost,
# enzyme not included as it's cost is more affected by the loading (considered later)
D = shape.Triangle(0.8*0.2573, 0.2573, 1.2*0.2573)
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-ton',
       baseline=0.2573, distribution=D)
def set_feedstock_price(price):
    feedstock.price = price / _feedstock_factor


# #!!! TODO: update
# D = shape.Triangle(0.198, 0.253, 0.304)
# @param(name='Natural gas unit price', element='TEA', kind='isolated', units='$/kg',
#        baseline=0.253, distribution=D)
# def set_natural_gas_price(price):
#     natural_gas.price = price


D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity unit price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(price):
    bst.PowerUtility.price = price


D = shape.Triangle(0.8*41, 41, 1.2*41)
@param(name='Activated carbon adsorbent price', element='TEA', kind='isolated', units='$/ft^3',
       baseline=41, distribution=D)
def set_adsorbent_price(price):
    AC401.adsorbent_cost['Activated carbon'] = price

#%% ######################## Conversion parameters ########################

# Fermentation
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element=R302, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading):
    R302.CSL_loading = loading

R302 = u.R302
# 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14)
D = shape.Triangle(0.9, 0.95, 1-1e-6)
@param(name='Seed train fermentation ratio', element=R302, kind='coupled', units='%',
       baseline=0.95, distribution=D)
def set_ferm_ratio(ratio):
    R303.ferm_ratio = ratio

D = shape.Triangle(baseline_yield*0.8, baseline_yield, baseline_yield*1.2)
@param(name='TAL yield', element=R302, kind='coupled', units='g/g',
       baseline=baseline_yield, distribution=D)
def set_TAL_yield(TAL_yield):
    # spec.load_specifications(TAL_yield,
    #                          spec.spec_2,
    #                          spec.spec_3)
    spec.spec_1 = TAL_yield

D = shape.Triangle(baseline_titer*0.8, baseline_titer, baseline_titer*1.2)
@param(name='TAL titer', element=R302, kind='coupled', units='g/L',
       baseline=baseline_titer, distribution=D)
def set_TAL_titer(TAL_titer):
    # spec.load_specifications(TAL_yield,
    #                          spec.spec_2,
    #                          spec.spec_3)
    spec.spec_2 = TAL_titer
    
D = shape.Triangle(baseline_productivity*0.8, baseline_productivity, baseline_productivity*1.2)
@param(name='Productivity', element=R302, kind='coupled', units='g/L/hr',
       baseline=baseline_productivity, distribution=D)
def set_TAL_productivity(productivity):
    # spec.load_specifications(TAL_yield,
    #                          spec.spec_2,
    #                          spec.spec_3)
    spec.spec_3 = productivity

D = shape.Triangle(0.8*0.005, 0.005, 1.2*0.005)
@param(name='VitaminA yield', element=R302, kind='coupled', units='g/g',
       baseline=0.005, distribution=D)
def set_vitaminA_yield(yield_):
    R302.glucose_to_VitaminA_rxn.X = yield_
    R303.glucose_to_VitaminA_rxn.X = R303.ferm_ratio*yield_

D = shape.Triangle(0.8*0.005, 0.005, 1.2*0.005)
@param(name='VitaminD2 yield', element=R302, kind='coupled', units='g/g',
       baseline=0.005, distribution=D)
def set_vitaminD2_yield(yield_):
    R302.glucose_to_VitaminD2_rxn.X = yield_
    R303.glucose_to_VitaminD2_rxn.X = R303.ferm_ratio*yield_

D = shape.Triangle(0.8*0.05, 0.05, 1.2*0.05)
@param(name='Microbe yield', element=R302, kind='coupled', units='g/g',
       baseline=0.05, distribution=D)
def set_microbe_yield(yield_):
    R302.glucose_to_microbe_rxn.X = yield_
    R303.glucose_to_microbe_rxn.X = R303.ferm_ratio*yield_
   
    
# D = shape.Triangle(0.05, 0.07, 0.1)
# @param(name='Inoculum ratio', element=R301, kind='coupled', units='%',
#        baseline=0.07, distribution=D)
# def set_inoculum_ratio(ratio):
#     R301.inoculum_ratio = ratio

#%%
######################## Separation parameters ########################

D = shape.Triangle(0.0739, 0.091, 0.091326667)
@param(name='Adsorbent unsaturated capacity', element=AC401, kind='coupled', units='g/g',
       baseline=0.091, distribution=D)
def set_adsorbent_cap(cap):
    AC401.adsorbent_capacity = cap

D = shape.Triangle(0.8*0.35, 0.35, 1.2*0.35)
@param(name='Adsorbent void volume fraction', element=AC401, kind='coupled', units='L/L',
       baseline=0.35, distribution=D)
def set_adsorbent_vvf(frac):
    AC401.void_fraction = frac

D = shape.Triangle(0.8*2050, 2050, 1.2*2050)
@param(name='Adsorbent bulk density', element=AC401, kind='coupled', units='kg/m^3',
       baseline=2050, distribution=D)
def set_adsorbent_bulk_rho(rho):
    AC401.rho_adsorbent = rho

D = shape.Triangle(0.8*0.078, 0.078, 1.2*0.078)
@param(name='Desorption single-wash partition coefficient', element=AC401, kind='coupled', units='(g/L)/(g/L)',
       baseline=0.078, distribution=D)
def set_desorption_K(K):
    AC401.K = K
    
#%%
######################## Facility parameters ########################
D = baseline_uniform(0.8, 0.1)
@param(name='boiler efficiency', element=BT, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_boiler_efficiency(efficiency):
    BT.B_eff = efficiency

D = shape.Triangle(0.8*7*24, 7*24, 1.2*7*24)
@param(name='Product TAL storage time', element=BT, kind='coupled', units='%',
       baseline=7*24, distribution=D)
def set_product_storage_time(storage_time):
    T620.tau = storage_time
    
    
parameters = model.get_parameters()



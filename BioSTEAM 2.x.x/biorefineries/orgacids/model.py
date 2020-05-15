#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:26:46 2020

@author: yalinli_cabbi
"""


'''
TODO:
    Add more (esp. on separation process) uncertainty parameters,
        consider using the bst.evaluation.State.load_default_parameters function
    
    Double-check the small errors in some sum calculations    
'''


# %% Setup

import numpy as np
import pandas as pd
from chaospy import distributions as shape
from biosteam.evaluation import Model, Metric
from biosteam.evaluation import evaluation_tools as tools
from orgacids.process_settings import price
from orgacids.system import lactic_acid, gypsum, M301, R301, R302, BT, CT, \
    orgacids_sys, boiler_sys, orgacids_sub_sys, orgacids_tea, orgacids_sys_no_boiler_tea,  \
    boiler_sys_tea
# from orgacids.system import *

#!!! Want to change to orgacids_sys_no_BT_tea
_hr_per_yr = orgacids_sys_no_boiler_tea._annual_factor
real_units = sum([i for i in orgacids_sub_sys.values()], ())


# %% Define metric functions

# Minimum selling price of lactic_acid stream
def get_MSP():
    for i in range(5):
        lactic_acid.price = orgacids_tea.solve_price(lactic_acid,
                                                      orgacids_sys_no_boiler_tea)
    return lactic_acid.price
# Mass flow rate of lactic_acid stream
get_yield = lambda: lactic_acid.F_mass
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: lactic_acid.imass['LacticAcid'] \
    /(R301.outs[1].imass['LacticAcid']+R301.outs[1].imass['CalciumLactate'])
get_overall_TCI = lambda: orgacids_tea.TCI
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: orgacids_tea.AOC
get_material_cost = lambda: orgacids_tea.material_cost
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: orgacids_tea.sales
# System power usage, individual unit power usage should be positive
unit_power_usage = [i.power_utility.rate for i in orgacids_sys.units if i.power_utility]
excess_power = lambda: BT.electricity_generated
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*price['Electricity']*_hr_per_yr)

metrics = [Metric('Minimum selling price', get_MSP, '$/kg'),
           Metric('Product yield', get_yield, 'kg/hr'),
           Metric('Product purity', get_purity, '%'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '$'),
           Metric('Annual operating cost', get_overall_AOC, '$/yr'),
           Metric('Annual material cost', get_material_cost, '$/yr'),
           Metric('Annual product sale', get_annual_sale, '$/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '$/yr')
           ]

# Get capital cost of sub-systems
def get_installation_cost(system):
    return lambda: sum(i.installation_cost for i in orgacids_sub_sys[system])
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend(
        (Metric(system, get_installation_cost(system), '$', 'Installation cost'),))

# All checks should be ~0
check_installation_cost = \
    lambda: sum(get_installation_cost(system)() 
                for system in orgacids_sub_sys.keys()) - orgacids_tea.installation_cost
metrics.extend((Metric('Check', check_installation_cost, '$', 'Installation cost'),))

# Get material cost
system_feeds = [i for i in orgacids_sys.feeds if i.price] + \
    [i for i in boiler_sys.feeds if i.price]
def get_material_cost(feed):
    return lambda: feed.price*feed.F_mass*_hr_per_yr
for feed in system_feeds:
    metrics.extend((Metric(feed.ID, get_material_cost(feed), '$/yr', 'Material cost'),))
check_material_cost = lambda: sum(get_material_cost(feed)()
                                  for feed in system_feeds) - orgacids_tea.material_cost

metrics.extend((Metric('Check', check_material_cost, '$/yr', 'Material cost'),))

# Get product sale
system_products = [i for i in orgacids_sys.products if i.price] + \
    [i for i in boiler_sys.products if i.price]
def get_product_sale(stream):
    return lambda: stream.price * stream.F_mass * _hr_per_yr
for product in system_products:
    metrics.extend((Metric(product.ID, get_product_sale(product), '$/yr', 'Product sale'),))
# Baseline disposal cost of gypsum is 0, thus not caputured
#!!! Need to add other baseline-0 wastes if considering their disposal costs in
# sensitivity and uncertainty analyses
get_gypsum_sale = lambda: gypsum.price * gypsum.F_mass * _hr_per_yr
check_product_sale= \
    lambda: sum(get_product_sale(product)() for product in system_products) \
        + get_gypsum_sale() - orgacids_tea.sales
metrics.extend((
    Metric(gypsum.ID, get_gypsum_sale, '$/yr', 'Product sale'),
    Metric('Check', check_product_sale, '$/yr', 'Product sale')))

# Get heating demand, positive if needs heating
get_system_heating_demand = lambda: BT.system_heating_demand
get_pretreatment_steam_heating_demand = lambda: BT.side_streams_lps.duty
get_BT_heating_demand = lambda: sum(i.duty for i in BT.heat_utilities if i.duty*i.cost>0)

def get_heating_demand(system):
    heat_utilities = sum((i.heat_utilities for i in orgacids_sub_sys[system]), ())
    heating_utilities = [i for i in heat_utilities if i.duty*i.cost>0]
    return lambda: sum(i.duty for i in heating_utilities)

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    # The only heating demand for the pretreatment system is the heat needed to
    # generate the side steam
    if system == 'pretreatment_sys':
        metrics.extend((Metric(system, get_pretreatment_steam_heating_demand, 'kJ/hr', 
                               'Heating demand'),))
    elif system == 'BT':
        metrics.extend((Metric(system, get_BT_heating_demand, 'kJ/hr', 
                               'Heating demand'),))
    else: metrics.extend((Metric(system, get_heating_demand(system), 'kJ/hr', 
                                 'Heating demand'),))

check_heating_demand = \
    lambda: sum((get_heating_demand(system)() 
                 for system in orgacids_sub_sys.keys() if system != 'BT'), 
                get_pretreatment_steam_heating_demand()) - get_system_heating_demand()
                                            # (get_side_streams_heating_demand(),))[0]
                                            
metrics.extend((
    Metric('Overall', get_system_heating_demand, 'kJ/hr', 'Heating demand'),
    Metric('Check', check_heating_demand, 'MJ/hr', 'Heating demand')
    ))

# Get cooling demand, negative if needs cooling
get_system_cooling_demand = lambda: CT.system_cooling_demand
get_CT_cooling_demand = lambda: sum(i.duty for i in CT.heat_utilities if i.duty*i.cost<0)

def get_cooling_demand(system):
    heat_utilities = sum((i.heat_utilities for i in orgacids_sub_sys[system]), ())
    cooling_utilities = [i for i in heat_utilities if i.duty*i.cost<0]
    # return lambda: sum(i.duty for i in cooling_utilities)/1e3/get_system_cooling_demand()
    return lambda: sum(i.duty for i in cooling_utilities)

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    elif system == 'CT':
        metrics.extend((Metric(system, get_CT_cooling_demand, 'kJ/hr', 'Cooling demand'),))
    else:
        metrics.extend((Metric(system, get_cooling_demand(system), 'kJ/hr', 'Cooling demand'),))

check_cooling_demand = \
    lambda: sum(get_cooling_demand(system)()
                for system in orgacids_sub_sys.keys() if system != 'CT') \
        - get_system_cooling_demand()

metrics.extend((
    Metric('Overall', get_system_cooling_demand, 'kJ/hr', 'Cooling demand'),
    Metric('Check', check_cooling_demand, 'kJ/hr', 'Cooling demand')
    ))

# Get power demand, positive if using power
get_system_power_demand = lambda: sum(i.power_utility.rate for i in orgacids_sys.units
                                      if i.power_utility)
def get_power_demand(system):
    power_utilities = [i.power_utility for i in orgacids_sub_sys[system]]
    return lambda: sum(i.rate for i in power_utilities)

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric(system, get_power_demand(system), 'kW', 'Power demand'),))

check_power_demand = lambda: sum(get_power_demand(system)()
                                 for system in orgacids_sub_sys.keys()) - get_system_power_demand()
metrics.extend((
    Metric('Overall', get_system_power_demand, 'kW', 'Power demand'),
    Metric('Check', check_power_demand, 'kW', 'Power demand')
    ))

# Get utility cost including heating, cooling, and power
get_system_utility_cost = lambda: sum(i.utility_cost for i in real_units)
def get_utility_cost(system):
    return lambda: sum(i.utility_cost for i in orgacids_sub_sys[system])

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric(system, get_utility_cost(system), '$/hr', 'Utility cost'),))

check_utility_cost = \
    lambda: sum(get_utility_cost(system)() for system in orgacids_sub_sys.keys()) \
        - get_system_utility_cost()

metrics.extend((
    Metric('Sum', get_system_utility_cost, '$/hr', 'Utility cost'),
    Metric('Overall', check_utility_cost, '$/hr', 'Utility cost')
    ))

# To see if TEA converges well for each simulation
get_NPV = lambda: orgacids_tea.NPV
metrics.extend((Metric('NPV', get_NPV, '$', 'TEA'), ))


# %% Add evaluating parameters
#!!! Parameters incomplete now

orgacids_model = Model(orgacids_sys, metrics)
param = orgacids_model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

# Costs
D = shape.Uniform(-0.0288, 0.00776)
@param(name='Gypsum price', element='TEA', kind='isolated', units='$/kg', 
       baseline=0, distribution=D)
def set_gypsum_price(cost): gypsum.price = cost

D = shape.Triangle(0.9, 1, 1.1)
@param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
       baseline=1, distribution=D)
def set_system_base_CAPEX(ratio): 
    for unit in orgacids_sys.units:
        if hasattr(unit, 'cost_items'):
            for item in unit.cost_items:
                unit.cost_items[item].cost *= ratio

D = shape.Uniform(0.01, 0.03)
@param(name='Enzyme loading', element=M301, kind='coupled', units='g/g cellulose',
        baseline=0.02, distribution=D)
def set_enzyme_loading(loading): M301.enzyme_loading = loading

# Fermentation
D = baseline_uniform(24, 0.1)
@param(name='Prehydrolysis time', element=R301, kind='coupled', units='hr',
        baseline=24, distribution=D)
def set_R301_prehydrolysis_time(tau): R301.tau_prehydrolysis = tau

D = shape.Triangle(76, 120, 145)
@param(name='Fermentation time', element=R301, kind='coupled', units='hr',
        baseline=120, distribution=D)
def set_R301_fermentation_time(tau): R301.tau_cofermentation = tau

D = shape.Uniform(5, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
        baseline=10, distribution=D)
def set_CSL_loading(loading): R301.CSL_loading = loading

D = shape.Triangle(0.55, 0.76, 0.93)
@param(name='Lactic acid yield', element=R301, kind='coupled', units='g/g substrate',
        baseline=0.76, distribution=D)
def set_R301_lactic_acid_yield(X):
    R301.cofermentation_rxns.X[0] = R301.cofermentation_rxns.X[3] = X

D = shape.Triangle(0.004, 0.069, 0.321)
@param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g substrate',
        baseline=0.069, distribution=D)
def set_R301_acetic_acid_yield(X): 
    # 1e6 is to avoid generating tiny negative flow (e.g., 1e-14) in R301
    X = min(X, 1-1e-6-R301.cofermentation_rxns.X[0]-R301.cofermentation_rxns.X[2])
    R301.cofermentation_rxns.X[1] = R301.cofermentation_rxns.X[4] = X

D = shape.Triangle(0.05, 0.1, 0.15)
@param(name='Innoculum ratio', element=R301, kind='coupled', units='%',
        baseline=0.1, distribution=D)
def set_innoculum_ratio(ratio): R301.inoculum_ratio = ratio

# Seed train fermentation yield as a ratio of the main fermenter
D = baseline_uniform(36, 0.1)
@param(name='Seed train time', element=R302, kind='coupled', units='hr',
        baseline=36, distribution=D)
def set_R302_fermentation_time(tau): R302.tau_batch = tau

D = shape.Triangle(0.8, 0.9, 1)
@param(name='Seed train yield', element=R302, kind='coupled', units='% of R301',
        baseline=0.9, distribution=D)
def set_R302_ratio(ratio):
    R302.ferm_ratio = ratio


# %% Metrics to evalute system across IRR coordinate

# Note: not using evaluate_across_coordinate function as the function will 
# simulate the system at each coordinate point, which is not necessary since
# changing IRR only affects the cashflow, not the system

def create_IRR_metrics(IRR):
    def get_IRR_based_MSP():
        # orgacids_tea.IRR = orgacids_sys_no_boiler_tea.IRR = boiler_sys_tea.IRR = IRR
        orgacids_tea.IRR = IRR
        return get_MSP()
    return [Metric('MSP', get_IRR_based_MSP, '$/kg', f'IRR={IRR:.0%}'),
            Metric('NPV', get_NPV, '$', f'IRR={IRR:.0%}')]

IRRs = np.linspace(0, 0.4, 41)
IRR_metrics = sum([create_IRR_metrics(IRR) for IRR in IRRs],[])

orgacids_model_IRR = orgacids_model.copy()
orgacids_model_IRR.metrics = IRR_metrics




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
    orgacids_tea, orgacids_sys_no_boiler_tea, orgacids_sys, orgacids_sub_sys, \
    boiler_sys_tea
# from orgacids.system import *

__all__ = ('orgacids_model', 'orgacids_model_IRR')

#!!! Want to change to orgacids_sys_no_BT_tea
_hr_per_yr = orgacids_sys_no_boiler_tea._annual_factor
real_units = sum([i for i in orgacids_sub_sys.values()], ())


# %% Define metric functions

# Minimum selling price of lactic_acid stream
get_MSP = lambda: [orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
                   for i in range(10)][-1]
# Mass flow rate of lactic_acid stream
get_yield = lambda: lactic_acid.F_mass
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: lactic_acid.imass['LacticAcid'] \
    /(R301.outs[1].imass['LacticAcid']+R301.outs[1].imass['CalciumLactate'])
get_overall_TCI = lambda: orgacids_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: orgacids_tea.AOC/1e6
get_material_costs = lambda: orgacids_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sales = lambda: orgacids_tea.sales/1e6
# System power usage, individual unit power usage should be positive
unit_power_usage = [i.power_utility.rate for i in orgacids_sys.units if i.power_utility]
excess_power = lambda: BT.generated_electricity
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*price['Electricity']*_hr_per_yr/1e6)

metrics = [Metric('Minimum selling price', get_MSP, '$/kg'),
           Metric('Product yield', get_yield, 'kg/hr'),
           Metric('Product purity', get_purity, '%'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr'),
           Metric('Annual material costs', get_material_costs, '10^6 $/yr'),
           Metric('Annual product sales', get_annual_sales, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr')
           ]

# Get capital cost of sub-systems (as ratios of overall capital cost)
def get_CAPEX_ratios(system):
    return lambda: sum([i.installation_cost 
                        for i in orgacids_sub_sys[system]])/orgacids_tea.installation_cost
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend(
        (Metric(system, get_CAPEX_ratios(system), '%', 'CAPEX ratio'),))
get_CAPEX_ratios_sum = lambda: sum(get_CAPEX_ratios(system)() 
                                   for system in orgacids_sub_sys.keys())
metrics.extend((Metric('Sum', get_CAPEX_ratios_sum, '%', 'CAPEX ratio'),))

# Get material costs (as ratios of overall costs)
system_feeds = tuple(i for i in orgacids_sys.feeds if i.price)
def get_material_cost_ratios(feed):
    return lambda: feed.price*feed.F_mass*_hr_per_yr/orgacids_tea.material_cost
for feed in system_feeds:
    metrics.extend((Metric(feed.ID, get_material_cost_ratios(feed), '%', 'Cost ratio'),))
get_get_material_cost_ratios_sum = lambda: sum(get_material_cost_ratios(feed)() 
                                               for feed in system_feeds)
metrics.extend((Metric('Sum', get_get_material_cost_ratios_sum, '%', 'Cost ratio'),))

# Get product sales (as ratios of overall sales)
system_products = tuple(i for i in orgacids_sys.products if i.price)
def get_product_sale_ratios(stream):
    return lambda: stream.price*stream.F_mass*_hr_per_yr/orgacids_tea.sales
for product in system_products:
    metrics.extend((Metric(product.ID, get_product_sale_ratios(product), '%', 'Sale ratio'),))
# Baseline disposal cost of gypsum is 0, thus not caputured
#!!! Need to add other baseline-0 wastes if considering their disposal costs in
# sensitivity and uncertainty analyses
get_gypsum_sale_ratio = lambda: gypsum.price*gypsum.F_mass*_hr_per_yr/orgacids_tea.sales
get_product_sale_ratios_sum = lambda: sum(get_product_sale_ratios(product)()
                                          for product in system_products)
metrics.extend((
    Metric(gypsum.ID, get_gypsum_sale_ratio, '%', 'Sale ratio'),
    Metric('Sum', get_product_sale_ratios_sum, '%', 'Sale ratio')))

# Get heating demand (kJ/hr), positive if needs heating
get_system_heating_demand = lambda: BT.system_heating_demand/1e3
def get_heating_demand_ratios(system):
    heat_utilities = sum([i.heat_utilities for i in orgacids_sub_sys[system]], ())
    heating_utilities = [i for i in heat_utilities if i.duty*i.cost>0]
    return lambda: sum([i.duty for i in heating_utilities])/1e3/get_system_heating_demand()
get_side_streams_heating_demand = lambda: BT.side_streams_lps.duty/1e3/get_system_heating_demand()
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    if system == 'pretreatment_sys':
        metrics.extend((Metric(system, get_side_streams_heating_demand, '%', 
                               'Heating demand ratio'),))
    else: metrics.extend((Metric(system, get_heating_demand_ratios(system), '%', 
                                 'Heating demand ratio'),))
get_heating_demand_ratios_sum = lambda: sum((get_heating_demand_ratios(system)()
                                            for system in orgacids_sub_sys.keys()
                                            if system != 'BT'),
                                            (get_side_streams_heating_demand(),))
metrics.extend((
    Metric('Sum', get_heating_demand_ratios_sum, '%', 'Heating demand ratio'),
    Metric('Overall', get_system_heating_demand, 'MJ/hr', 'Heating demand')
    ))

# Get cooling demand (kJ/hr), negative if needs cooling
get_system_cooling_demand = lambda: CT.system_cooling_demand/1e3
def get_cooling_demand_ratios(system):
    heat_utilities = sum([i.heat_utilities for i in orgacids_sub_sys[system]], ())
    cooling_utilities = [i for i in heat_utilities if i.duty*i.cost<0]
    return lambda: sum([i.duty for i in cooling_utilities])/1e3/get_system_cooling_demand()
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric(system, get_cooling_demand_ratios(system), '%', 'Cooling demand ratio'),))
get_cooling_demand_ratios_sum = lambda: sum(get_cooling_demand_ratios(system)()
                                            for system in orgacids_sub_sys.keys()
                                            if system != 'CT')
metrics.extend((
    Metric('Sum', get_cooling_demand_ratios_sum, '%', 'Cooling demand ratio'),
    Metric('Overall', get_system_cooling_demand, 'MJ/hr', 'Cooling demand')
    ))

# Get power demand (kW), positive if using power
get_system_power_demand = lambda: sum(i.power_utility.rate for i in orgacids_sys.units
                                      if i.power_utility)/1e3
def get_power_demand_ratios(system):
    power_utilities = [i.power_utility for i in orgacids_sub_sys[system]]
    return lambda: sum([i.rate for i in power_utilities])/1e3/get_system_power_demand()
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric(system, get_power_demand_ratios(system), '%', 'Power demand ratio'),))
get_power_demand_ratios_sum = lambda: sum(get_power_demand_ratios(system)()
                                          for system in orgacids_sub_sys.keys())
metrics.extend((
    Metric('Sum', get_power_demand_ratios_sum, '%', 'Power demand ratio'),
    Metric('Overall', get_system_power_demand, 'MW', 'Power demand')
    ))

# Get utility cost including heating, cooling, and power
get_system_utility_cost = lambda: sum(i.utility_cost for i in real_units)*_hr_per_yr/1e6
def get_utility_cost_ratios(system):
    return lambda: sum([i.utility_cost for i in orgacids_sub_sys[system]]) \
        *_hr_per_yr/1e6/get_system_utility_cost()
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric('Utility cost ratio', get_utility_cost_ratios(system), '%', system),))
get_utility_cost_ratios_sum = lambda: sum(get_utility_cost_ratios(system)()
                                          for system in orgacids_sub_sys.keys())
metrics.extend((
    Metric('Utility cost ratio', get_utility_cost_ratios_sum, '%', 'Sum'),
    Metric('Utility cost', get_system_utility_cost, '10^6 $/yr', 'Overall')
    ))

# To see if TEA converges well for each simulation
get_NPV = lambda: orgacids_tea.NPV
metrics.extend((Metric('TEA', get_NPV, '$', 'NPV'), ))


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

def create_IRR_metric(IRR):
    def get_IRR_based_MSP():
        orgacids_tea.IRR = orgacids_sys_no_boiler_tea.IRR = boiler_sys_tea.IRR = IRR
        return [orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea) 
                for i in range(10)][-1]
    return Metric(f'MSP at IRR={IRR}', get_IRR_based_MSP, '$/kg', 'TEA')
#     return Metric(f'MSP at IRR={int(IRR*100)}%', get_IRR_based_MSP, '$/kg', 'TEA')

IRRs = np.linspace(0, 0.4, 41)
IRR_metrics = [create_IRR_metric(IRR) for IRR in IRRs] + [metrics[-1]]

orgacids_model_IRR = orgacids_model.copy()
orgacids_model_IRR.metrics = IRR_metrics




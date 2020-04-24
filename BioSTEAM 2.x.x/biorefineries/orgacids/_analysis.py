#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:24:42 2020

@author: yalinli_cabbi
"""


# %% Setup

from pandas import ExcelWriter
from chaospy import distributions as shape
from biosteam.evaluation import Model, Metric
from biosteam.evaluation import evaluation_tools as tools
from orgacids.process_settings import price
from orgacids.system import *

_hr_per_yr = orgacids_sys_no_BT_tea._annual_factor


# %% Define metric functions

# Minimum selling price of lactic_acid stream ($/kg)
get_MSP = lambda: orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_BT_tea)
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: lactic_acid.imass['LacticAcid'] \
                       /(R301.outs[1].imass['LacticAcid']+R301.outs[1].imass['CalciumLactate'])
get_overall_TCI = lambda: orgacids_tea.TCI/1e6
# Note that operating cost excludes waste disposal costs
get_overall_AOC = lambda: orgacids_tea.AOC/1e6
get_material_costs = lambda: orgacids_tea.material_cost/1e6
# Annual sale revenue from products, note that AOC excludes electricity credit
# but includes negative sales from waste disposal
# (i.e., wastes are products of negative selling price)
get_annual_sales = lambda: orgacids_tea.sales/1e6
# Sysem power usage excluding BT, individual unit power usage should be positive
unit_power_usage = [i.power_utility.rate for i in orgacids_sys.units
                    if i.power_utility and i is not BT]
# Excess power is positive if generating more electricity than needed,
# BT power usage is negative if generating more power than used
excess_power = lambda: (-BT.power_utility.rate-sum(unit_power_usage))
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*price['Electricity']*_hr_per_yr/1e6)

metrics = [Metric('Minimum selling price', get_MSP, '$/kg'),
           Metric('Product purity', get_purity, '%'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr'),
           Metric('Annual material costs', get_material_costs, '10^6 $/yr'),
           Metric('Annual product sales', get_annual_sales, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr')
           ]

# Get total capital investment of sub-systems (as ratios of overall TCI)
def get_TCI_ratios(sub_system_tea): 
    return lambda: sub_system_tea.TCI/orgacids_tea.TCI
# Get annual operating cost (as ratios of overall AOC)
def get_AOC_ratios(sub_system_tea): 
    return lambda: sub_system_tea.AOC/orgacids_tea.AOC

for tea in sub_teas:
    metrics.extend(
        (Metric('TCI', get_TCI_ratios(tea), '%', tea.system.ID),
         Metric('AOC', get_AOC_ratios(tea), '%', tea.system.ID)))

# Get material costs (as ratios of overall costs)
system_feeds = tuple(i for i in orgacids_sys.feeds if i.price)
def get_cost_ratios(feed):
    return lambda: feed.price*feed.F_mass*_hr_per_yr/orgacids_tea.material_cost
for feed in system_feeds:
    metrics.extend((Metric('Cost ratio', get_cost_ratios(feed), '%', feed.ID),))

# Get product sales (as ratios of overall sales)
system_products = tuple(i for i in orgacids_sys.products if i.price)
def get_sale_ratios(stream):
    return lambda: stream.price*stream.F_mass*_hr_per_yr/orgacids_tea.sales
for product in system_products:
    metrics.extend((Metric('Sale ratio', get_sale_ratios(product), '%', product.ID),))
# Baseline disposal cost of gypsum is 0, thus not caputured
#!!! Need to add other baseline-0 wastes if considering their disposal costs in
# sensitivity and uncertainty analyses
get_gypsum_sale_ratio = lambda: gypsum.price*gypsum.F_mass*_hr_per_yr/orgacids_tea.sales
metrics.extend((Metric('Sale ratio', get_gypsum_sale_ratio, '%', gypsum.ID),))

# Get heating demand (kJ/hr), positive if needs heating
def get_heating_demand(system):
    heat_utilities = sum([i.heat_utilities for i in system.units], ())
    heating_utilities = [i for i in heat_utilities if i.ID=='low_pressure_steam']
    return lambda: sum([i.duty for i in heating_utilities])/1e3

# Get cooling demand (kJ/hr), negative if needs cooling
def get_cooling_demand(system):
    heat_utilities = sum([i.heat_utilities for i in system.units], ())
    cooling_utilities = [i for i in heat_utilities if i.ID=='cooling_water']
    return lambda: sum([i.duty for i in cooling_utilities])/1e3

# Get power demand (kW), positive if using power
def get_power_demand(system):
    power_utilities = [i.power_utility for i in system.units]
    return lambda: sum([i.rate for i in power_utilities])/1e3

# Get utility cost, note that as heating/cooling agents are regenerated onset,
# utility cost is acutally power (electricity) cost
def get_utility_cost(system):
    return lambda: sum([i.utility_cost for i in system.units])*_hr_per_yr/1e6

for system in sub_systems:
    metrics.extend(
        (Metric('Heating demand', get_heating_demand(system), 'MJ/hr', system.ID),
         Metric('Cooling demand', get_cooling_demand(system), 'MJ/hr', system.ID),
         Metric('Power demand', get_power_demand(system), 'MW', system.ID),
         Metric('Utility cost', get_utility_cost(system), '10^6 $/yr', system.ID)
        ))

metric_number = len(metrics)


# %% Add evaluating parameters
#!!! Parameters incomplete now

orgacids_model = Model(orgacids_sys, metrics)
param = orgacids_model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

def baseline_triangle(baseline, ratio):
    return shape.Triangle(baseline*(1-ratio), baseline, baseline*(1+ratio))

# Costs
D = shape.Uniform(-0.00776, 0.042)
@param(name='Gypsum disposal', element='TEA', kind='isolated', units='$/kg', 
       baseline=0, distribution=D)
def set_gypsum_cost(cost): gypsum.price = cost

D = baseline_triangle(1, 0.1)
@param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
       baseline=1, distribution=D)
def set_system_base_TCI(ratio): 
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
@param(name='Batch fermentation time', element=R301, kind='coupled', units='hr',
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
    R301.cofermentation_rxns.X[0] = R301.cofermentation_rxns.X[2] = X

D = shape.Triangle(0.004, 0.069, 0.321)
@param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g substrate',
       baseline=0.069, distribution=D)
def set_R301_acetic_acid_yield(X): 
    X = min(X, 1-R301.cofermentation_rxns.X[0])
    R301.cofermentation_rxns.X[1] = R301.cofermentation_rxns.X[3] = X

D = shape.Triangle(0.05, 0.1, 0.15)
@param(name='Innoculum ratio', element=R301, kind='coupled', units='%',
       baseline=0.1, distribution=D)
def set_innoculum_ratio(ratio): R301.inoculum_ratio = ratio

# Seed train
D = baseline_uniform(36, 0.1)
@param(name='Residence time', element=R302, kind='coupled', units='hr',
       baseline=36, distribution=D)
def set_R302_fermentation_time(tau): R302.tau_batch = tau

D = baseline_triangle(0.9, 0.1)
@param(name='Lactic acid yield', element=R302, kind='coupled', units='g/g substrate',
       baseline=0.9, distribution=D)
def set_R302_lactic_acid_yield(ratio):
    X = R301.cofermentation_rxns[0].X * ratio
    R302.cofermentation_rxns.X[0] = R302.cofermentation_rxns.X[2] = X

D = baseline_triangle(0.9, 0.1)
@param(name='Acetic acid yield', element=R302, kind='coupled', units='g/g substrate',
       baseline=0.9, distribution=D)
def set_R302_acetic_acid_yield(ratio):
    X = R301.cofermentation_rxns[1].X * ratio
    X = min(X, 1-R302.cofermentation_rxns.X[0])
    R302.cofermentation_rxns.X[1] = R302.cofermentation_rxns.X[3] = X

parameter_number = len(orgacids_model.get_baseline_sample())


# %% Run evaluation and output results

#!!! Maybe can use `orgacids_model.evaluate_across_coordinate` to look at
# the impact of heavy acid content (maybe also latic acid purity?) vs. MSP

samples = orgacids_model.sample(N=25, rule='L')
orgacids_model.load_samples(samples)
orgacids_model.evaluate()

parameter_table = orgacids_model.table.iloc[:, :parameter_number]
Monte_Carlo_results = orgacids_model.table.iloc[:, parameter_number::]
percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]
Monte_Carlo_percentiles = Monte_Carlo_results.quantile(q=percentiles)

# Note that if only one metric is used, then need to make sure it's a tuple with
# spearman_metrics = (metrics[0],)
spearman_metrics = metrics[0:3]
spearman_results = orgacids_model.spearman(spearman_metrics)

with ExcelWriter('Analysis results.xlsx') as writer:
    parameter_table.to_excel(writer, sheet_name='Parameters')
    Monte_Carlo_results.to_excel(writer, sheet_name='Monte Carlo')
    Monte_Carlo_percentiles.to_excel(writer, sheet_name='Monte Carlo Percentiles')
    spearman_results.to_excel(writer, sheet_name='Spearman')
    












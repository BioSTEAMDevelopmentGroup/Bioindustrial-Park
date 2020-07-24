#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon May 11 16:26:46 2020

Modified from the biorefineries constructed in [1] and [2] for the production of
lactic acid from lignocellulosic feedstocks

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted.
    July, 2020.

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import biosteam as bst
import lactic.system as system
from chaospy import distributions as shape
from biosteam.evaluation import Model, Metric

lactic_no_CHP_tea = system.lactic_no_CHP_tea
get_annual_factor = lambda: lactic_no_CHP_tea._annual_factor
_kg_per_ton = 907.18474

lactic_sys = system.lactic_sys
CHP_sys = system.CHP_sys
system_feeds = [i for i in lactic_sys.feeds if i.price] + \
    [i for i in CHP_sys.feeds if i.price]
system_products = [i for i in lactic_sys.products if i.price] + \
    [i for i in CHP_sys.products if i.price]

gypsum = system.gypsum
system_products.append(gypsum)


# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

# Minimum product selling price of lactic_acid stream
lactic_acid = system.lactic_acid
lactic_tea = system.lactic_tea
def get_MPSP():
    for i in range(3):
        lactic_acid.price = lactic_tea.solve_price(lactic_acid)
    return lactic_acid.price

# Mass flow rate of lactic_acid stream
get_yield = lambda: lactic_acid.F_mass*get_annual_factor()/1e6
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
R301 = system.R301
get_recovery = lambda: lactic_acid.imol['LacticAcid'] \
    /(R301.outs[0].imol['LacticAcid']+2*R301.outs[0].imol['CalciumLactate'])
get_overall_TCI = lambda: lactic_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: lactic_tea.AOC/1e6
get_material_cost = lambda: lactic_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: lactic_tea.sales/1e6
# System power usage, individual unit power usage should be positive
CHP = system.CHP
excess_power = lambda: CHP.electricity_generated
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('Minimum product selling price', get_MPSP, '$/kg'),
           Metric('Product yield', get_yield, '10^6 kg/yr'),
           Metric('Product purity', get_purity, '%'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr')
           ]

# =============================================================================
# Capital cost breakdown
# =============================================================================

process_groups = system.process_groups
def get_installed_cost(group):
    return lambda: group.get_installed_cost()
for group in process_groups:
    if group.name == 'feedstock_group': continue
    metrics.extend(
        (Metric(group.name, get_installed_cost(group), '10^6 $', 'Installed cost'),))

# All checks should be ~0
check_installed_cost = \
    lambda: sum(get_installed_cost(group)() 
                for group in process_groups) - lactic_tea.installed_equipment_cost/1e6
metrics.extend((Metric('Check', check_installed_cost, '10^6 $', 'Installed cost'),))


# =============================================================================
# Material cost breakdown
# =============================================================================

def get_material_cost(feed):
    return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
for feed in system_feeds:
    metrics.extend((Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'),))
lime_R301 = system.lime_R301
lime_CHP = system.lime_CHP
get_fermentation_lime_ratio = lambda: lime_R301.imol['Lime'] \
    / (lime_R301.imol['Lime']+lime_CHP.imol['Lime']) 
T602_S = system.T602_S
get_separation_sulfuric_acid_ratio = lambda: T602_S.outs[1].F_mass/T602_S.F_mass_out
check_material_cost = lambda: sum(get_material_cost(feed)()
                                  for feed in system_feeds) - lactic_tea.material_cost/1e6

metrics.extend((
    Metric('Fermentation lime ratio', get_fermentation_lime_ratio, 
           '%', 'Material cost'),
    Metric('Separation sulfuric acid ratio', get_separation_sulfuric_acid_ratio, 
           '%', 'Material cost'),
    Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')))

def get_product_sale(stream):
    return lambda: stream.price*stream.F_mass*get_annual_factor()/1e6
for product in system_products:
    metrics.extend((Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'),))
check_product_sale= \
    lambda: sum(get_product_sale(product)() for product in system_products) \
        - lactic_tea.sales/1e6
metrics.extend((Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'),))


# =============================================================================
# Heating demand breakdown (positive if needs heating)
# =============================================================================

get_system_heating_demand = lambda: CHP.system_heating_demand*get_annual_factor()/1e9
get_pretreatment_steam_heating_demand = lambda: CHP.side_streams_lps.duty*get_annual_factor()/1e9
HXN = system.HXN
get_HXN_heating_demand = lambda: sum(i.duty for i in HXN.heat_utilities 
                                    if i.duty*i.cost>0)*get_annual_factor()/1e9
get_CHP_heating_demand = lambda: sum(i.duty for i in CHP.heat_utilities 
                                    if i.duty*i.cost>0)*get_annual_factor()/1e9

def get_heating_demand(group):
    return lambda: group.get_heating_duty()*get_annual_factor()/1e3

for group in process_groups:
    if group.name in ('feedstock_group', 'HXN_group', 'CHP_group'): continue
    # The only heating demand for the pretreatment system is the heat needed to
    # generate the side steam
    if group.name == 'pretreatment_group':
        metrics.extend((Metric(group.name, get_pretreatment_steam_heating_demand, '10^9 kJ/yr', 
                                'Heating demand'),))
    else: metrics.extend((Metric(group.name, get_heating_demand(group), '10^9 kJ/yr', 
                                  'Heating demand'),))

check_heating_demand = \
    lambda: sum((get_heating_demand(group)() for group in process_groups), 
                get_pretreatment_steam_heating_demand())
                                            
metrics.extend((
    Metric('HXN', get_HXN_heating_demand, '10^9 kJ/yr', 'Heating demand'),
    Metric('CHP', get_CHP_heating_demand, '10^9 kJ/yr', 'Heating demand'),    
    Metric('Sum', get_system_heating_demand, '10^9 kJ/yr', 'Heating demand'),
    Metric('Check', check_heating_demand, '10^9 kJ/yr', 'Heating demand')
    ))

# =============================================================================
# Cooling demand breakdown (negative if needs cooling)
# =============================================================================

CT = system.CT
get_system_cooling_water_duty = lambda: CT.system_cooling_water_duty*get_annual_factor()/1e9
get_HXN_cooling_demand = lambda: sum(i.duty for i in HXN.heat_utilities 
                                    if i.duty*i.cost<0)*get_annual_factor()/1e9
get_CT_cooling_water_duty = lambda: sum(i.duty for i in CT.heat_utilities 
                                    if i.duty*i.cost<0)*get_annual_factor()/1e9

def get_cooling_demand(group):
    return lambda: group.get_heating_duty()*get_annual_factor()/1e3

for group in process_groups:
    if group.name in ('feedstock_group', 'HXN_group', 'CT_group'): continue
    else: metrics.extend((Metric(group.name, get_cooling_demand(group),
                                  '10^9 kJ/yr', 'Cooling demand'),))

check_cooling_demand = \
    lambda: sum(get_cooling_demand(group)() for group in process_groups)

metrics.extend((
    Metric('HXN', get_HXN_cooling_demand, '10^9 kJ/yr', 'Cooling demand'),
    Metric('CT', get_CT_cooling_water_duty, '10^9 kJ/yr', 'Cooling demand'),    
    Metric('Sum', get_system_cooling_water_duty, '10^9 kJ/yr', 'Cooling demand'),
    Metric('Check', check_cooling_demand, '10^9 kJ/yr', 'Cooling demand')
    ))

# =============================================================================
# Power demand breakdown (positive if using power)
# =============================================================================

get_system_power_demand = lambda: sum(i.power_utility.rate for i in lactic_sys.units
                                      if i.power_utility)

def get_power_demand(group):
    return lambda: sum(i.rate for i in group.power_utilities)

for group in process_groups:
    if group.name == 'feedstock_group': continue
    metrics.extend((Metric(group.name, get_power_demand(group), 'kW', 'Power demand'),))

check_power_demand = lambda: sum(get_power_demand(group)()
                                 for group in process_groups) - get_system_power_demand()
metrics.extend((
    Metric('Sum', get_system_power_demand, 'kW', 'Power demand'),
    Metric('Check', check_power_demand, 'kW', 'Power demand')
    ))

# =============================================================================
# Utility cost breakdown (including heating, cooling, and power)
# =============================================================================

get_system_utility_cost = lambda: lactic_tea.utility_cost/1e6

def get_utility_cost(group):
    return lambda: sum(i.utility_cost for i in group.units)*get_annual_factor()/1e6

for group in process_groups:
    if group.name == 'feedstock_group': continue
    metrics.extend((Metric(group.name, get_utility_cost(group), '10^6 $/yr', 'Utility cost'),))

check_utility_cost = \
    lambda: sum(get_utility_cost(group)() for group in process_groups) \
        - get_system_utility_cost()

metrics.extend((
    Metric('Sum', get_system_utility_cost, '10^6 $/yr', 'Utility cost'),
    Metric('Check', check_utility_cost, '10^6 $/yr', 'Utility cost')
    ))

# To see if TEA converges well for each simulation
get_NPV = lambda: lactic_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))


# %% 

# =============================================================================
# Construct base model
# =============================================================================

model_full = Model(lactic_sys, metrics)
param = model_full.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))


# =============================================================================
# TEA parameters
# =============================================================================

U101 = system.U101
D = baseline_uniform(2205, 0.1)
@param(name='Flow rate', element=U101, kind='coupled', units='dry U.S. ton/day',
       baseline=2205, distribution=D)
def set_feedstock_flow_rate(rate):
    U101.feedstock_flow_rate = rate

D = shape.Uniform(0.9, 1)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.96, distribution=D)
def set_operating_days(uptime):
    lactic_tea.operating_days = 365 * uptime

special_price = {
#     stream           distribution     min        max
    'feedstock':        ('Uniform',   (0.0529,    0.117)),
    'CSL':              ('Uniform',   (0.0673,    0.112)),
    'lime':             ('Uniform',   (0.160,     0.288)),
    'ethanol':          ('Triangle',  (0.460,     0.978)),
    'natural_gas':      ('Triangle',  (0.198,     0.304)),
    'gypsum':           ('Uniform',   (-0.0288,   0.00776))
    }

# Prices for boiler_chems, baghouse_bag, and cooling_tower_chems are not included
# as they are tied to CHP/CT duties
default_price_streams = ('sulfuric_acid', 'ammonia', 'enzyme_R301', 
                         'caustic_R502', 'ash')

Setter = bst.evaluation.evaluation_tools.parameter.Setter
def add_stream_price_param(stream, D):
    param(setter=Setter(stream, 'price'),
          name=f'{stream.ID} price', element='TEA', kind='isolated', units='$/kg',
          baseline=stream.price, distribution=D)

for stream_ID in special_price.keys():
    stream = getattr(system, stream_ID)
    lower = special_price[stream_ID][1][0]
    mid = stream.price
    upper = special_price[stream_ID][1][-1]
    if special_price[stream_ID][0] == 'Triangle':
        D = shape.Triangle(lower, mid, upper)
    elif special_price[stream.ID][0] == 'Uniform':
        D = shape.Uniform(lower, upper)
    add_stream_price_param(stream, D)

for stream_ID in default_price_streams:
    stream = getattr(system, stream_ID)
    baseline = stream.price
    D = baseline_uniform(baseline, 0.1)
    add_stream_price_param(stream, D)

D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(price): 
    bst.PowerUtility.price = price

D = baseline_uniform(1, 0.25)
@param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
        baseline=1, distribution=D)
def set_TCI_ratio(new_ratio):
    old_ratio = lactic_no_CHP_tea._TCI_ratio_cached
    for unit in lactic_sys.units:
        if hasattr(unit, 'cost_items'):
            for item in unit.cost_items:
                unit.cost_items[item].cost /= old_ratio
                unit.cost_items[item].cost *= new_ratio
    lactic_no_CHP_tea._TCI_ratio_cached = new_ratio


# =============================================================================
# Pretreatment parameters
# =============================================================================

M202 = system.M202
D = shape.Uniform(0.25, 0.4)
@param(name='Pretreatment solid loading', element=M202, kind='coupled', units='%', 
       baseline=0.3, distribution=D)
def set_pretreatment_solid_loading(loading): 
    M202.solid_loading = loading
    
T201 = system.T201
D = shape.Uniform(10, 35)
@param(name='Pretreatment sulfuric acid loading', element=T201,
       kind='coupled', units='mg/g-dry feedstock', baseline=22.1, distribution=D)
def set_pretreatment_sulfuric_acid_loading(loading): 
    T201.feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    T201.acid_loading = loading

R201 = system.R201
D = shape.Uniform(0.06, 0.12)
@param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='%',
       baseline=0.099, distribution=D)
def set_R201_glucan_conversion(X): R201.pretreatment_rxns[0].X = X    

D = shape.Uniform(0.8, 0.92)
@param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R201_xylan_conversion(X): R201.pretreatment_rxns[4].X = X        

# =============================================================================
# Conversion parameters
# =============================================================================

M301 = system.M301
D = shape.Uniform(0.175, 0.25)
@param(name='Enzymatic hydrolysis solid loading', element=M301, kind='coupled', units='%',
       baseline=0.2, distribution=D)
def set_R301_hydrolysis_solid_loading(loading): M301.solid_loading = loading

D = shape.Uniform(10, 30)
@param(name='Enzyme loading', element=M301, kind='coupled', units='mg/g glucan',
       baseline=20, distribution=D)
def set_R301_enzyme_loading(loading): M301.enzyme_loading = loading

# Enzymatic hydrolysis
D = baseline_uniform(24, 0.1)
@param(name='Enzymatic hydrolysis time', element=R301, kind='coupled', units='hr',
       baseline=24, distribution=D)
def set_R301_hydrolysis_time(tau): R301.tau_saccharification = tau

D = shape.Uniform(0.75, 0.9)
@param(name='Enzymatic hydrolysis glucan-to-glucose', element=R301, kind='coupled', units='%',
       baseline=0.85, distribution=D)
def set_R301_glucan_conversion(X): R301.saccharification_rxns[2].X = X

# Fermentation
D = shape.Uniform(5, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading): R301.CSL_loading = loading

R302 = system.R302
D = baseline_uniform(0.9, 0.1)
@param(name='Seed train fermentation ratio', element=R302, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_ferm_ratio(ratio):
    R302.ferm_ratio = ratio

D = shape.Triangle(0.18, 0.89, 1.92)
@param(name='Productivity', element=R301, kind='coupled', units='hr',
       baseline=0.89, distribution=D)
def set_productivity(productivity):
    R301.productivity = productivity
    R302.productivity = productivity * R302.ferm_ratio

D = shape.Triangle(0.55, 0.76, 0.93)
@param(name='Lactic acid yield', element=R301, kind='coupled', units='g/g sugar',
       baseline=0.76, distribution=D)
def set_lactic_acid_yield(lactic_yield):
    R301_X = R301.cofermentation_rxns.X
    R302_X = R302.cofermentation_rxns.X
    R301_X[0] = R301_X[3] = lactic_yield
    # 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14) in R301
    R301_X[1] = R301_X[4] = min(R301_X[1], 1-1e-6-R301_X[0]-R301_X[2])
    R302_X[0] = R302_X[3] = lactic_yield * R302.ferm_ratio
    R302_X[1] = R302_X[4] = R301_X[1] * R302.ferm_ratio

D = shape.Triangle(0.004, 0.07, 0.32)
@param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g sugar',
       baseline=0.07, distribution=D)
def set_acetic_acid_yield(acetic_yield): 
    R301_X = R301.cofermentation_rxns.X
    R302_X = R302.cofermentation_rxns.X
    R301_acetic_yield = min(acetic_yield, 1-1e-6-R301_X[0]-R301_X[2])
    R301_X[1] = R301_X[4] = R301_acetic_yield
    R302_acetic_yield = min(R301_acetic_yield, 1-1e-6-R302_X[0]-R302_X[2])
    R302_X[1] = R302_X[4] = R302_acetic_yield * R302.ferm_ratio

D = shape.Uniform(0.05, 0.1)
@param(name='Inoculum ratio', element=R301, kind='coupled', units='%',
       baseline=0.07, distribution=D)
def set_innoculum_ratio(ratio): R301.inoculum_ratio = ratio


# =============================================================================
# Separation parameters
# =============================================================================

S402 = system.S402
D = shape.Triangle(0.95, 0.995, 1)
@param(name='Gypsum split', element=S402, kind='coupled', units='',
       baseline=0.995, distribution=D)
def set_S402_gypsum_split(split):
    gypsum_index = S402.chemicals.index('Gypsum')
    S402.split[gypsum_index] = split

R401 = system.R401
D = baseline_uniform(0.5, 0.1)
@param(name='Acidulation time', element=R401, kind='coupled', units='hr',
       baseline=0.5, distribution=D)
def set_R401_tau(tau):
    R401.tau = tau

R402 = system.R402
D = shape.Triangle(0.95, 1, 1.05)
@param(name='Esterification conversion factor', element=R402, kind='coupled', units='',
       baseline=1, distribution=D)
def set_R402_conversion_factor(factor):
    R402.X_factor = factor
    
R403 = system.R403
D = shape.Triangle(0.72, 0.8, 0.88)
@param(name='Hydrolysis conversion', element=R403, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_R403_conversion_factor(X):
    R403.hydrolysis_rxns.X[:] = X
    

# =============================================================================
# Facilities parameters
# =============================================================================

# Facilities are currently not in unit path. thus not set the element here
D = baseline_uniform(0.8, 0.1)
@param(name='CHP combustion efficiency', element=CHP, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_CHP_combustion_efficiency(efficiency):
    CHP.B_eff = efficiency

# All parameters
parameters = model_full.get_parameters()


# %%

# =============================================================================
# Model to evalute system across internal rate of return
# =============================================================================

def create_IRR_metrics(IRR):
    def get_IRR_based_MPSP():
        lactic_tea.IRR = IRR
        return get_MPSP()
    return [Metric('Minimum product selling price', get_IRR_based_MPSP, '$/kg', f'IRR={IRR:.0%}'),
            Metric('Net present value', get_NPV, '$', f'IRR={IRR:.0%}')]

IRRs = np.linspace(0, 0.4, 41)
IRR_metrics = sum([create_IRR_metrics(IRR) for IRR in IRRs],[])

model_IRR = Model(lactic_sys, IRR_metrics)
model_IRR.set_parameters(parameters)


# %%

# =============================================================================
# Model to evalute system across feedstock price and carbohydate content
# =============================================================================

def create_price_metris(price):
    def get_price_based_MPSP():
        price_per_kg = price / _kg_per_ton * 0.8
        feedstock.price = price_per_kg
        return get_MPSP()
    return [Metric('Minimum product selling price', get_price_based_MPSP,
                   '$/kg', f'Price={price:.0f} [$/dry-ton]'),
            Metric('Net present value', get_NPV,
                   '$', f'Price={price:.0f} [$/dry-ton]')]

prices = np.linspace(50, 300, 26)
# 71.3 is the baseline
prices = np.concatenate((prices, np.array([71.3])))

price_metrics = sum([create_price_metris(price) for price in prices],[])

def set_carbs(carbs_content):
    carbs = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')
    dry_mass = feedstock.F_mass.copy() - feedstock.imass[('H2O',)].copy()
    old_carbs_mass_total = feedstock.imass[carbs].sum().copy()
    ratio = feedstock.get_normalized_mass(carbs)
    new_carbs_mass = dry_mass * carbs_content * ratio
    feedstock.set_flow(new_carbs_mass, 'kg/hr', carbs)
    mass_diff = new_carbs_mass.sum() - old_carbs_mass_total
    feedstock.imass['Extractives'] -= mass_diff
    if any(feedstock.mass < 0):
        raise ValueError(f'Carbohydrate content of {carbs_content*100:.0f}% dry weight is infeasible')

model_feedstock = Model(lactic_sys, price_metrics)

param = model_feedstock.parameter

# Set a fake parameter to enable evaluation across internal rate of return
feedstock = system.feedstock
D = shape.Uniform(0.9, 1.1)
@param(name='Fake parameter', element=feedstock, kind='coupled', units='',
       baseline=1, distribution=D)
def set_fake_parameter(anything):
    # This is to make sure all system is registered in the flowsheet
    feedstock.price = feedstock.price


# %%

# =============================================================================
# Model to evalute system across feedstock succinic acid content
# =============================================================================

def set_succinic(content):
    dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    feedstock.imass['SuccinicAcid'] = content * dry_mass
    # Use Extractives to close mass balance
    feedstock.imass['Extractives'] -= (feedstock.F_mass-feedstock.imass['H2O']) - dry_mass
    if any(feedstock.mass<0):
        raise ValueError(f'Succinic acid content of {content*100:.0f}% dry weight is infeasible')

model_succinic = Model(lactic_sys, metrics)
model_succinic.set_parameters(parameters)







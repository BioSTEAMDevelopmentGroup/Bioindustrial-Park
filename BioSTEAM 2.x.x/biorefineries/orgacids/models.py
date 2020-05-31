#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 16:26:46 2020

@author: yalinli_cabbi
"""


'''
TODO:
    Add separation process uncertainty parameters  
'''


# %% Setup

import numpy as np
import biosteam as bst
from chaospy import distributions as shape
from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import Setter
from orgacids.system import orgacids_sub_sys, orgacids_tea, orgacids_no_BT_tea

get_annual_factor = lambda: orgacids_no_BT_tea._annual_factor
_kg_per_ton = 907.18474
real_units = sum([i for i in orgacids_sub_sys.values()], ())

orgacids_sys = find.system.orgacids_sys
BT_sys = find.system.BT_sys

system_feeds = [i for i in orgacids_sys.feeds if i.price] + \
    [i for i in BT_sys.feeds if i.price]
system_products = [i for i in orgacids_sys.products if i.price] + \
    [i for i in BT_sys.products if i.price]
    
gypsum = find.stream.gypsum
system_products.append(gypsum)


# %% Define metric functions

# =============================================================================
# Overall biorefinery
# =============================================================================

# Minimum selling price of lactic_acid stream
def get_MSP():
    for i in range(5):
        lactic_acid.price = orgacids_tea.solve_price(lactic_acid,
                                                      orgacids_no_BT_tea)
    return lactic_acid.price
# Mass flow rate of lactic_acid stream
lactic_acid = find.stream.lactic_acid
get_yield = lambda: lactic_acid.F_mass*get_annual_factor()/1e6
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
R301 = find.unit.R301
get_recovery = lambda: lactic_acid.imol['LacticAcid'] \
    /(R301.outs[0].imol['LacticAcid']+2*R301.outs[0].imol['CalciumLactate'])
get_overall_TCI = lambda: orgacids_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: orgacids_tea.AOC/1e6
get_material_cost = lambda: orgacids_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: orgacids_tea.sales/1e6
# System power usage, individual unit power usage should be positive
BT = find.unit.BT
excess_power = lambda: BT.electricity_generated
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('Minimum selling price', get_MSP, '$/kg'),
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

def get_installation_cost(system):
    return lambda: sum(i.installation_cost for i in orgacids_sub_sys[system])/1e6
for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend(
        (Metric(system, get_installation_cost(system), '10^6 $', 'Installed cost'),))

# All checks should be ~0
check_installation_cost = \
    lambda: sum(get_installation_cost(system)() 
                for system in orgacids_sub_sys.keys()) - orgacids_tea.installation_cost/1e6
metrics.extend((Metric('Check', check_installation_cost, '10^6 $', 'Installed cost'),))

# =============================================================================
# Material cost breakdown
# =============================================================================

def get_material_cost(feed):
    return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
for feed in system_feeds:
    metrics.extend((Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'),))
S602 = find.unit.S602
get_fermentation_lime_ratio = lambda: S602.outs[0].imol['Lime']/S602.ins[0].imol['Lime']
S601 = find.unit.S601
get_separation_sulfuric_acid_ratio = lambda: S601.outs[1].imol['H2SO4']/S601.ins[0].imol['H2SO4']
check_material_cost = lambda: sum(get_material_cost(feed)()
                                  for feed in system_feeds) - orgacids_tea.material_cost/1e6

metrics.extend((
    Metric('Fermentation lime ratio', get_fermentation_lime_ratio, 
           '%', 'Material cost'),
    Metric('Separation sulfuric acid ratio', get_separation_sulfuric_acid_ratio, 
           '%', 'Material cost'),
    Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')))

def get_product_sale(stream):
    return lambda: stream.price*stream.F_mass *get_annual_factor()/1e6
for product in system_products:
    metrics.extend((Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'),))
check_product_sale= \
    lambda: sum(get_product_sale(product)() for product in system_products) \
        - orgacids_tea.sales/1e6
metrics.extend((Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'),))

# =============================================================================
# Heating demand breakdown (positive if needs heating)
# =============================================================================

get_system_heating_demand = lambda: BT.system_heating_demand*get_annual_factor()/1e9
get_pretreatment_steam_heating_demand = lambda: BT.side_streams_lps.duty*get_annual_factor()/1e9
get_BT_heating_demand = lambda: sum(i.duty for i in BT.heat_utilities 
                                    if i.duty*i.cost>0)*get_annual_factor()/1e9

def get_heating_demand(system):
    heat_utilities = sum((i.heat_utilities for i in orgacids_sub_sys[system]), ())
    heating_utilities = [i for i in heat_utilities if i.duty*i.cost>0]
    return lambda: sum(i.duty for i in heating_utilities)*get_annual_factor()/1e9

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    # The only heating demand for the pretreatment system is the heat needed to
    # generate the side steam
    if system == 'pretreatment_sys':
        metrics.extend((Metric(system, get_pretreatment_steam_heating_demand, '10^9 kJ/yr', 
                               'Heating demand'),))
    elif system == 'BT':
        metrics.extend((Metric(system, get_BT_heating_demand, '10^9 kJ/yr', 
                               'Heating demand'),))
    else: metrics.extend((Metric(system, get_heating_demand(system), '10^9 kJ/yr', 
                                 'Heating demand'),))

check_heating_demand = \
    lambda: sum((get_heating_demand(system)() 
                 for system in orgacids_sub_sys.keys() if system != 'BT'), 
                get_pretreatment_steam_heating_demand()) - get_system_heating_demand()
                                            
metrics.extend((
    Metric('Sum', get_system_heating_demand, '10^9 kJ/yr', 'Heating demand'),
    Metric('Check', check_heating_demand, '10^9 kJ/yr', 'Heating demand')
    ))

# =============================================================================
# Cooling demand breakdown (negative if needs cooling)
# =============================================================================

CT = find.unit.CT
get_system_cooling_demand = lambda: CT.system_cooling_demand*get_annual_factor()/1e9
get_CT_cooling_demand = lambda: sum(i.duty for i in CT.heat_utilities 
                                    if i.duty*i.cost<0)*get_annual_factor()/1e9

def get_cooling_demand(system):
    heat_utilities = sum((i.heat_utilities for i in orgacids_sub_sys[system]), ())
    cooling_utilities = [i for i in heat_utilities if i.duty*i.cost<0]
    # return lambda: sum(i.duty for i in cooling_utilities)/1e3/get_system_cooling_demand()
    return lambda: sum(i.duty for i in cooling_utilities)*get_annual_factor()/1e9

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    elif system == 'CT':
        metrics.extend((Metric(system, get_CT_cooling_demand, '10^9 kJ/yr', 'Cooling demand'),))
    else:
        metrics.extend((Metric(system, get_cooling_demand(system), '10^9 kJ/yr', 'Cooling demand'),))

check_cooling_demand = \
    lambda: sum(get_cooling_demand(system)()
                for system in orgacids_sub_sys.keys() if system != 'CT') \
        - get_system_cooling_demand()

metrics.extend((
    Metric('Sum', get_system_cooling_demand, '10^9 kJ/yr', 'Cooling demand'),
    Metric('Check', check_cooling_demand, '10^9 kJ/yr', 'Cooling demand')
    ))

# =============================================================================
# Power demand breakdown (positive if using power)
# =============================================================================

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
    Metric('Sum', get_system_power_demand, 'kW', 'Power demand'),
    Metric('Check', check_power_demand, 'kW', 'Power demand')
    ))

# =============================================================================
# Utility cost breakdown, (including heating, cooling, and power)
# =============================================================================

get_system_utility_cost = lambda: sum(i.utility_cost for i in real_units)*get_annual_factor()/1e6
def get_utility_cost(system):
    return lambda: sum(i.utility_cost for i in orgacids_sub_sys[system])*get_annual_factor()/1e6

for system in orgacids_sub_sys.keys():
    if system == 'feedstock_sys': continue
    metrics.extend((Metric(system, get_utility_cost(system), '10^6 $/yr', 'Utility cost'),))

check_utility_cost = \
    lambda: sum(get_utility_cost(system)() for system in orgacids_sub_sys.keys()) \
        - get_system_utility_cost()

metrics.extend((
    Metric('Sum', get_system_utility_cost, '10^6 $/yr', 'Utility cost'),
    Metric('Check', check_utility_cost, '10^6 $/yr', 'Utility cost')
    ))

# To see if TEA converges well for each simulation
get_NPV = lambda: orgacids_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))


# %% Add evaluating parameters
#!!! Need to add separation parameters

# This is the base model
orgacids_model = Model(orgacids_sys, metrics)
param = orgacids_model.parameter

def baseline_uniform(baseline, ratio):
    return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

# =============================================================================
# Feedstock
# =============================================================================

U101 = find.unit.U101
D = baseline_uniform(2205, 0.1)
@param(name='Flow rate', element=U101, kind='coupled', units='U.S. ton/day',
       baseline=2205, distribution=D)
def set_feedstock_flow_rate(rate):
    U101.feedstock_flow_rate = rate

feedstock = find.stream.feedstock
# Use Extract to close mass balance
D = shape.Triangle(0.265, 0.319, 0.376)
@param(name='Glucan content', element=feedstock, kind='coupled', 
        units='dry mass %', baseline=0.3505, distribution=D)
def set_feedstock_glucan_content(content):
    dry_mass = feedstock.F_mass-feedstock.imass['H2O']
    old_content = feedstock.imass['Glucan'] / dry_mass
    mass_diff = (content-old_content) * dry_mass
    feedstock.imass['Glucan'] += mass_diff
    feedstock.imass['Extract'] -= mass_diff

D = shape.Triangle(0.148, 0.189, 0.227)
@param(name='Xylan content', element=feedstock, kind='coupled', 
        units='dry mass %', baseline=0.1953, distribution=D)
def set_feedstock_xylan_content(content):
    dry_mass = feedstock.F_mass-feedstock.imass['H2O']
    old_content = feedstock.imass['Xylan'] / dry_mass
    mass_diff = (content-old_content) * dry_mass
    feedstock.imass['Xylan'] += mass_diff
    feedstock.imass['Extract'] -= mass_diff
    
D = shape.Uniform(0, 0.05)
@param(name='Succinic acid content', element=feedstock, kind='coupled', 
        units='dry mass %', baseline=0, distribution=D)
def set_feedstock_succinic_acid_content(content):
    dry_mass = feedstock.F_mass-feedstock.imass['H2O']
    old_content = feedstock.imass['SuccinicAcid'] / dry_mass
    mass_diff = (content-old_content) * dry_mass
    feedstock.imass['SuccinicAcid'] += mass_diff
    feedstock.imass['Extract'] -= mass_diff

# =============================================================================
# TEA
# =============================================================================

D = baseline_uniform(0.9, 0.1)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.9, distribution=D)
def set_operating_days(uptime):
    orgacids_tea.operating_days = 365 * uptime

# Impactful parameters are set to triangular distribution based on literature,
# less important ones are set to ±10% of baseline value
special_price = {
#     stream           distribution     min        max
    'feedstock':        ('Uniform',   (0.0529,    0.117)),
    'CSL_fresh':        ('Uniform',   (0.0673,    0.112)),
    'lime_fresh':       ('Uniform',   (0.160,     0.288)),
    'ethanol_fresh':    ('Triangle',  (0.460,     0.978)),
    'natural_gas':      ('Triangle',  (0.198,     0.304)),
    'gypsum':           ('Uniform',   (-0.0288,   0.00776))
    }
default_price_streams = (
    'sulfuric_acid_fresh', 'ammonia_fresh', 'enzyme', 'system_makeup_water',
    'aerobic_caustic', 'cooling_tower_chems', 'boiler_chems', 'baghouse_bag', 'ash',
    )

def add_stream_price_param(stream, D):
    param(setter=Setter(stream, 'price'),
          name=f'{stream.ID} price', element='TEA', kind='isolated', units='$/kg',
          baseline=stream.price, distribution=D)

for stream_ID in special_price.keys():
    stream = getattr(find.stream, stream_ID)
    lower = special_price[stream_ID][1][0]
    mid = stream.price
    upper = special_price[stream_ID][1][-1]
    if special_price[stream_ID][0] == 'Triangle':
        D = shape.Triangle(lower, mid, upper)
    elif special_price[stream.ID][0] == 'Uniform':
        D = shape.Uniform(lower, upper)
    add_stream_price_param(stream, D)

for stream_ID in default_price_streams:
    stream = getattr(find.stream, stream_ID)
    baseline = stream.price
    D = baseline_uniform(baseline, 0.1)
    add_stream_price_param(stream, D)

D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity price', element='TEA', kind='isolated', units='$/kWh',
        baseline=0.070, distribution=D)
def set_electricity_price(price): 
    bst.PowerUtility.price = price

D = baseline_uniform(1, 0.1)
@param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
       baseline=1, distribution=D)
def set_TCI_ratio(ratio): 
    for unit in orgacids_sys.units:
        if hasattr(unit, 'cost_items'):
            for item in unit.cost_items:
                unit.cost_items[item].cost *= ratio

# =============================================================================
# Pretreatment
# =============================================================================

M202 = find.unit.M202
D = shape.Uniform(0.25, 0.4)
@param(name='Pretreatment solid loading', element=M202, kind='coupled', units='%', 
       baseline=0.3, distribution=D)
def set_pretreatment_solid_loading(loading): 
    M202.solid_loading = loading
    
pretreatment_sulfuric_acid = find.stream.pretreatment_sulfuric_acid
D = shape.Uniform(10, 35)
@param(name='Pretreatment sulfuric acid loading', element=pretreatment_sulfuric_acid,
       kind='coupled', units='mg/g-dry feedstock', baseline=22.1, distribution=D)
def set_pretreatment_sulfuric_acid_loading(loading): 
    feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    pretreatment_sulfuric_acid.imass['H2SO4'] = feedstock_dry_mass*loading/1000*0.93
    pretreatment_sulfuric_acid.imass['H2O'] = feedstock_dry_mass*loading/1000*0.07

R201 = find.unit.R201
D = shape.Uniform(0.06, 0.12)
@param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='%',
       baseline=0.099, distribution=D)
def set_R201_glucan_conversion(X): R201.pretreatment_rxns[0].X = X    

D = shape.Uniform(0.8, 0.92)
@param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R201_xylan_conversion(X): R201.pretreatment_rxns[4].X = X        

# =============================================================================
# Conversion
# =============================================================================

M301 = find.unit.M301
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
D = shape.Triangle(76, 120, 145)
@param(name='Fermentation time', element=R301, kind='coupled', units='hr',
        baseline=120, distribution=D)
def set_R301_fermentation_time(tau): R301.tau_cofermentation = tau

D = shape.Uniform(5, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
        baseline=10, distribution=D)
def set_CSL_loading(loading): R301.CSL_loading = loading

R302 = find.unit.R302
D = shape.Triangle(0.55, 0.76, 0.93)
@param(name='Lactic acid yield', element=R301, kind='coupled', units='g/g substrate',
        baseline=0.76, distribution=D)
def set_R301_lactic_acid_yield(X):
    R301_X = R301.cofermentation_rxns.X
    R301_X[0] = R301_X[3] = X
    R302_X = R302.cofermentation_rxns.X
    R302_X[0] = R302_X[3] = X * R302.ferm_ratio

D = shape.Triangle(0.004, 0.07, 0.32)
@param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g substrate',
        baseline=0.07, distribution=D)
def set_R301_acetic_acid_yield(X): 
    # 1e6 is to avoid generating tiny negative flow (e.g., 1e-14) in R301
    R301_X = R301.cofermentation_rxns.X
    X = min(X, 1-1e-6-R301_X[0]-R301_X[2])
    R301_X[1] = R301_X[4] = X
    R302_X = R302.cofermentation_rxns.X
    X = min(X*R302.ferm_ratio, 1-1e-6-R302_X[0]-R302_X[2])
    R302_X[1] = R302_X[4] = X

D = shape.Uniform(0.05, 0.1)
@param(name='Innoculum ratio', element=R301, kind='coupled', units='%',
        baseline=0.07, distribution=D)
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
    R301_X = R301.cofermentation_rxns.X
    R302_X = R302.cofermentation_rxns.X
    ratio = min(ratio, (1-1e-6-R302_X[2])/(R301_X[0]+R301_X[1]))
    R302.ferm_ratio = ratio

# =============================================================================
# Separation
# =============================================================================


parameters = orgacids_model.get_parameters()


# %% Construct other models

# =============================================================================
# Model to evalute system across lactic acid yield
# =============================================================================

orgacids_model_LA_yield = Model(orgacids_sys, metrics)

def set_LA_yield(LA_yield):
    R301_X = R301.cofermentation_rxns.X
    R301_X[0] = R301_X[3] = LA_yield
    R301_X[1] = R301_X[4] = min(1-1e-6-R301_X[0]-R301_X[2], R301_X[1])
    R302_X = R302.cofermentation_rxns.X
    R302_X[0] = R302_X[3] = R301_X[0] * R302.ferm_ratio
    R302_X[1] = R302_X[4] = R301_X[1] * R302.ferm_ratio

LA_yield_parameters = tuple([i for i in parameters if not i.name=='Lactic acid yield'])
orgacids_model_LA_yield.set_parameters(LA_yield_parameters)

# =============================================================================
# Model to evalute system across internal rate of return
# =============================================================================

def create_IRR_metrics(IRR):
    def get_IRR_based_MSP():
        orgacids_tea.IRR = IRR
        return get_MSP()
    return [Metric('Minimum selling price', get_IRR_based_MSP, '$/kg', f'IRR={IRR:.0%}'),
            Metric('Net present value', get_NPV, '$', f'IRR={IRR:.0%}')]

IRRs = np.linspace(0, 0.4, 41)
IRR_metrics = sum([create_IRR_metrics(IRR) for IRR in IRRs],[])

orgacids_model_IRR = Model(orgacids_sys, IRR_metrics)
orgacids_model_IRR.set_parameters(parameters)

# =============================================================================
# Model to evalute system across feedstock price and carbohydate content
# =============================================================================

def create_feedstock_price_metris(price):
    def get_price_based_MSP():
        price_per_kg = price / _kg_per_ton * 0.8
        feedstock.price = price_per_kg
        return get_MSP()
    return [Metric('Minimum selling price', get_price_based_MSP, '$/kg', f'Price={price:.0f} [$/dry-ton]'),
            Metric('Net present value', get_NPV, '$', f'Price={price:.0f} [$/dry-ton]')]

prices = np.linspace(50, 300, 26)
feedstock_price_metrics = sum([create_feedstock_price_metris(price) 
                               for price in prices],[])

def set_feedstock_carbs(content):
    ratio = feedstock.imass['Glucan'] / \
        (feedstock.imass['Glucan']+feedstock.imass['Xylan'])
    dry_mass = feedstock.F_mass-feedstock.imass['H2O']
    old_content = sum(feedstock.imass[i] 
                      for i in ('Glucan', 'Xylan', 'Arabinan', 
                                'Galactan', 'Mannan')) / dry_mass
    mass_diff = (content - old_content) * dry_mass
    feedstock.imass['Glucan'] += mass_diff * ratio
    feedstock.imass['Xylan'] += mass_diff * (1-ratio)
    feedstock.imass['Extract'] -= mass_diff

orgacids_model_feedstock = Model(orgacids_sys, feedstock_price_metrics)

param = orgacids_model_feedstock.parameter

# Set a fake parameter to enable evaluation
D = shape.Uniform(0.9, 1.1)
@param(name='Fake parameter', element='TEA', kind='isolated', units='',
       baseline=1, distribution=D)
def set_fake_parameter(anything): pass





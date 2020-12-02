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

References:
[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted,
    2020.
    
@author: yalinli_cabbi
"""


# %%

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import biosteam as bst
from biosteam.evaluation import Model, Metric
from chaospy import distributions as shape
from biorefineries.lactic._process_settings import CFs
from biorefineries.lactic._utils import set_yield
from biorefineries.lactic import system

_kg_per_ton = 907.18474
_feedstock_factor = _kg_per_ton / 0.8


# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

# Minimum product selling price of lactic_acid stream
lactic_acid = system.lactic_acid
lactic_tea = system.lactic_tea
def get_MPSP():
    lactic_acid.price = 0
    for i in range(3):
        MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
    return MPSP

feedstock = system.feedstock
# Yield in 10^6 kg/yr
lactic_no_CHP_tea = system.lactic_no_CHP_tea
get_annual_factor = lambda: lactic_tea.operating_days*24
get_total_yield = lambda: lactic_acid.F_mass*get_annual_factor()/1e6
# Yield in % of dry feedstock
get_mass_yield = lambda: lactic_acid.F_mass/(feedstock.F_mass-feedstock.imass['Water'])
R301 = system.R301
get_titer = lambda: R301.effluent_titer
# Purity (%) of LacticAcid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: lactic_acid.imol['LacticAcid'] \
    /(R301.outs[0].imol['LacticAcid']+2*R301.outs[0].imol['CalciumLactate'])
get_overall_TCI = lambda: lactic_tea.TCI/1e6
get_lactic_sale = lambda: get_total_yield()*lactic_acid.price
# Including negative product sales (ash/gypsum disposal) but excluding electricity credit
gypsum = system.gypsum
get_gypsum_sale = lambda: gypsum.F_mass*gypsum.price*get_annual_factor()/1e6
ash = system.ash
get_ash_sale = lambda: ash.F_mass*ash.price*get_annual_factor()/1e6
get_operating_cost = lambda: lactic_tea.AOC/1e6-get_gypsum_sale()-get_ash_sale()
# Including negative product sales (ash/gypsum disposal) but excluding electricity credit
get_material_cost = lambda: lactic_tea.material_cost/1e6-get_gypsum_sale()-get_ash_sale()
# System power usage, individual unit power usage should be positive
CHP = system.CHP
excess_power = lambda: CHP.electricity_generated
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('MPSP', get_MPSP, '$/kg'),
           Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
           Metric('Product mass yield', get_mass_yield, '%'),
           Metric('Fermentation titer', get_titer, 'g/L'),
           Metric('Product purity', get_purity, '%'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $'),
           Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
           Metric('Annual product sale', get_lactic_sale, '10^6 $/yr'),
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
    metrics.append(
        Metric(group.name, get_installed_cost(group), '10^6 $', 'Installed cost'))

# All checks should be ~0
check_installed_cost = \
    lambda: sum(get_installed_cost(group)() 
                for group in process_groups) - lactic_tea.installed_equipment_cost/1e6
metrics.append(Metric('Check', check_installed_cost, '10^6 $', 'Installed cost'))


# =============================================================================
# Material cost and product sale breakdown
# =============================================================================

TEA_feeds = system.TEA_feeds
TEA_products = system.TEA_products

def get_material_cost(feed):
    return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
for feed in TEA_feeds:
    metrics.append(Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'))

# Sulfuric acid used in pretreatment and separation processes
T602_S = system.T602_S
get_pretreatment_sulfuric_acid_ratio = lambda: T602_S.outs[0].F_mass/T602_S.F_mass_out
# Ammonia used in pretreatment and CHP process
T603_S = system.T603_S
get_pretreatment_ammonia_ratio = lambda: T603_S.outs[0].F_mass/T602_S.F_mass_out

check_material_cost = lambda: sum(get_material_cost(feed)()
                                  for feed in TEA_feeds) - lactic_tea.material_cost/1e6

metrics.extend((
    Metric('Pretreatment sulfuric acid ratio', get_pretreatment_sulfuric_acid_ratio, 
           '%', 'Material cost'),
    Metric('Pretreatment ammonia ratio', get_pretreatment_ammonia_ratio, 
           '%', 'Material cost'),
    Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')
    ))

def get_product_sale(stream):
    return lambda: stream.price*stream.F_mass*get_annual_factor()/1e6
for product in TEA_products:
    metrics.append(Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'))
check_product_sale= \
    lambda: sum(get_product_sale(product)() for product in TEA_products) \
        - lactic_tea.sales/1e6
metrics.append(Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'))


# =============================================================================
# Heating demand breakdown (positive if needs heating)
# =============================================================================

get_system_heating_demand = lambda: CHP.system_heating_demand/1e6
get_pretreatment_steam_heating_demand = lambda: CHP.side_streams_lps.duty/1e6
HXN = system.HXN
get_HXN_heating_demand = lambda: sum(i.duty for i in HXN.heat_utilities
                                    if i.flow*i.duty>0)/1e6
get_CHP_heating_demand = lambda: sum(i.duty for i in CHP.heat_utilities 
                                    if i.flow*i.duty>0)/1e6

def get_heating_demand(group):
    return lambda: group.get_heating_duty()

for group in process_groups:
    if group.name in ('feedstock_group', 'HXN_group', 'CHP_group'): continue
    # The only heating demand for the pretreatment system is the heat needed to
    # generate the side steam
    if group.name == 'pretreatment_group':
        metrics.append(Metric(group.name, get_pretreatment_steam_heating_demand, '10^6 kJ/hr', 
                                'Heating demand'))
    else: metrics.append(Metric(group.name, get_heating_demand(group), '10^6 kJ/hr', 
                                  'Heating demand'))

check_heating_demand = \
    lambda: sum((get_heating_demand(group)() for group in process_groups), 
                get_pretreatment_steam_heating_demand()) + \
                get_HXN_heating_demand() + get_CHP_heating_demand()
                                            
metrics.extend((
    Metric('HXN', get_HXN_heating_demand, '10^6 kJ/hr', 'Heating demand'),
    Metric('CHP', get_CHP_heating_demand, '10^6 kJ/hr', 'Heating demand'),    
    Metric('Sum', get_system_heating_demand, '10^6 kJ/hr', 'Heating demand'),
    Metric('Check', check_heating_demand, '10^6 kJ/hr', 'Heating demand')
    ))

# =============================================================================
# Cooling demand breakdown (negative if needs cooling)
# =============================================================================

CT = system.CT
get_system_cooling_water_duty = lambda: CT.system_cooling_water_duty/1e6
get_HXN_cooling_demand = lambda: sum(i.duty for i in HXN.heat_utilities 
                                    if i.flow*i.duty<0)/1e6
get_CT_cooling_water_duty = lambda: sum(i.duty for i in CT.heat_utilities 
                                    if i.flow*i.duty<0)/1e6

def get_cooling_demand(group):
    return lambda: -group.get_cooling_duty()

for group in process_groups:
    if group.name in ('feedstock_group', 'HXN_group', 'CT_group'): continue
    else: metrics.append(Metric(group.name, get_cooling_demand(group),
                                  '10^9 kJ/yr', 'Cooling demand'))

check_cooling_demand = \
    lambda: sum(get_cooling_demand(group)() for group in process_groups) + \
                get_HXN_cooling_demand() + get_CT_cooling_water_duty()

metrics.extend((
    Metric('HXN', get_HXN_cooling_demand, '10^6 kJ/hr', 'Cooling demand'),
    Metric('CT', get_CT_cooling_water_duty, '10^6 kJ/hr', 'Cooling demand'),    
    Metric('Sum', get_system_cooling_water_duty, '10^6 kJ/hr', 'Cooling demand'),
    Metric('Check', check_cooling_demand, '10^6 kJ/hr', 'Cooling demand')
    ))

# =============================================================================
# Power demand breakdown (positive if using power)
# =============================================================================

lactic_sys = system.lactic_sys
get_electricity_use = system.get_electricity_use

def get_power_demand(group):
    return lambda: sum(i.rate for i in group.power_utilities)

for group in process_groups:
    if group.name == 'feedstock_group': continue
    metrics.append(Metric(group.name, get_power_demand(group), 'kW', 'Power demand'))

check_power_demand = lambda: sum(get_power_demand(group)()
                                 for group in process_groups) - get_electricity_use()
metrics.extend((
    Metric('Sum', get_electricity_use, 'kW', 'Power demand'),
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
    metrics.append(Metric(group.name, get_utility_cost(group), '10^6 $/yr', 'Utility cost'))

check_utility_cost = \
    lambda: sum(get_utility_cost(group)() for group in process_groups) \
        - get_system_utility_cost()

metrics.extend((
    Metric('Sum', get_system_utility_cost, '10^6 $/yr', 'Utility cost'),
    Metric('Check', check_utility_cost, '10^6 $/yr', 'Utility cost')
    ))

# To see if TEA converges well for each simulation
get_NPV = lambda: lactic_tea.NPV
metrics.append(Metric('NPV', get_NPV, '$', 'TEA'))
metrics_no_IRRs = metrics.copy()

index_TEA = len(metrics)


# %%

# =============================================================================
# Metrics to evalute system at different internal rate of return metrics
# =============================================================================

def create_IRR_metrics(IRR):
    def get_IRR_based_MPSP():
        lactic_tea.IRR = IRR
        return get_MPSP()
    return [Metric('MPSP', get_IRR_based_MPSP, '$/kg', f'IRR={IRR:.0%}'),
            Metric('NPV', get_NPV, '$', f'IRR={IRR:.0%}')]

# This is to ensure Monte Carlo results will be at 10% IRR
IRR1 = np.arange(0, 0.10, 0.01)
IRR2 = np.arange(0.11, 0.41, 0.01)
IRRs = IRR1.tolist() + IRR2.tolist() + [0.1]
for IRR in IRRs:
    metrics.extend((i for i in create_IRR_metrics(IRR)))

index_IRR = len(metrics)


# %%

# =============================================================================
# Metrics for global warming potential
# =============================================================================

get_GWP = system.get_GWP
get_FEC = system.get_FEC
metrics.extend((
    Metric('Total GWP', get_GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Total FEC', get_FEC, 'MJ/kg', 'LCA')
    ))

get_material_GWP = system.get_material_GWP
get_electricity_GWP = system.get_electricity_GWP

natural_gas = system.natural_gas
get_CH4_production_GWP = lambda: \
    CFs['GWP_CFs']['CH4']*natural_gas.F_mass/lactic_acid.F_mass
get_CH4_onsite_GWP = lambda: \
     natural_gas.get_atomic_flow('C')*natural_gas.chemicals.CO2.MW/lactic_acid.F_mass

lime = system.lime
lime_CHP = system.lime_CHP
get_lime_GWP  = lambda: \
    CFs['GWP_CFs']['Lime']*(lime.F_mass+lime_CHP.F_mass)/lactic_acid.F_mass

ethanol = system.ethanol
get_ethanol_onsite_GWP = lambda: \
     ethanol.get_atomic_flow('C')*ethanol.chemicals.CO2.MW/lactic_acid.F_mass

# For all other materials, including both production and combustion
get_other_materials_GWP = lambda: get_material_GWP()+get_ethanol_onsite_GWP()-\
    get_CH4_production_GWP()-get_lime_GWP()

check_GWP = lambda: get_GWP()-get_material_GWP()-get_electricity_GWP()- \
    get_CH4_onsite_GWP()-get_ethanol_onsite_GWP()
    
metrics.extend((
    Metric('Electricity', get_electricity_GWP, 'kg CO2-eq/kg', 'GWP'),
    Metric('Natural gas production', get_CH4_production_GWP, 'kg CO2-eq/kg', 'GWP'),
    Metric('Natural gas combustion', get_CH4_onsite_GWP, 'kg CO2-eq/kg', 'GWP'), 
    Metric('Lime', get_lime_GWP, 'kg CO2-eq/kg', 'GWP'),
    Metric('Other materials', get_other_materials_GWP, 'kg CO2-eq/kg', 'GWP'),
    Metric('GWP check', check_GWP, 'kg CO2-eq/kg', 'GWP')
    ))

get_electricity_FEC = system.get_electricity_FEC
get_CH4_FEC = lambda: CFs['FEC_CFs']['CH4']*natural_gas.F_mass/lactic_acid.F_mass
get_lime_FEC = lambda: \
    CFs['FEC_CFs']['Lime']*(lime.F_mass+lime_CHP.F_mass)/lactic_acid.F_mass

get_material_FEC = system.get_material_FEC
get_other_materials_FEC = lambda: get_material_FEC()-get_CH4_FEC()-get_lime_FEC()



check_FEC = lambda: get_FEC() - get_material_FEC() - \
    get_system_utility_cost()*1e6/get_annual_factor()/bst.PowerUtility.price * \
    CFs['FEC_CFs']['Electricity']/lactic_acid.F_mass

metrics.extend((
    Metric('Electricity', get_electricity_FEC, 'MJ/kg', 'FEC'),
    Metric('Natural gas', get_CH4_FEC, 'MJ/kg', 'FEC'),
    Metric('Lime', get_lime_FEC, 'MJ/kg', 'FEC'),
    Metric('Other materials', get_other_materials_FEC, 'MJ/kg', 'FEC'),
    Metric('FEC check', check_FEC, 'MJ/kg', 'FEC')
    ))


# %% 

# =============================================================================
# Construct base model
# =============================================================================

model_full = Model(lactic_sys, metrics)
param = model_full.parameter

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

# =============================================================================
# TEA parameters
# =============================================================================

# U101 = system.U101
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
    lactic_tea.operating_days = 365 * uptime

D = baseline_triangle(1, 0.25)
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

# Only include materials that account for >5% of total annual material cost,
# enzyme not included as it's cost is more affected by the loading (considered later)
D = shape.Triangle(60, 71.3, 83.7)
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-ton',
       baseline=71.3, distribution=D)
def set_feedstock_price(price):
    feedstock.price = price / _feedstock_factor

sulfuric_acid = system.sulfuric_acid
D = shape.Triangle(0.0910, 0.0948, 0.1046)
@param(name='Sulfuric acid unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0.0948, distribution=D)
def set_sulfuric_acid_price(price):
    sulfuric_acid.price = price

D = shape.Triangle(0.160, 0.262, 0.288)
@param(name='Lime unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0.262, distribution=D)
def set_lime_price(price):
    lime.price = price
    lime_CHP.price = price

D = shape.Triangle(0.198, 0.253, 0.304)
@param(name='Natural gas unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0.253, distribution=D)
def set_natural_gas_price(price):
    natural_gas.price = price

D = shape.Uniform(-0.0288, 0.00776)
@param(name='Gypsum unit price', element='TEA', kind='isolated', units='$/kg',
       baseline=0, distribution=D)
def set_gypsum_price(price):
    gypsum.price = price

D = shape.Triangle(0.067, 0.070, 0.074)
@param(name='Electricity unit price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(price): 
    bst.PowerUtility.price = price


# =============================================================================
# Pretreatment parameters
# =============================================================================

M202 = system.M202
D = shape.Triangle(0.25, 0.3, 0.4)
@param(name='Pretreatment solid loading', element=M202, kind='coupled', units='%', 
       baseline=0.3, distribution=D)
def set_pretreatment_solid_loading(loading): 
    M202.solid_loading = loading
    
T201 = system.T201
D = shape.Triangle(10, 22.1, 35)
@param(name='Pretreatment sulfuric acid loading', element=T201,
       kind='coupled', units='mg/g', baseline=22.1, distribution=D)
def set_pretreatment_sulfuric_acid_loading(loading): 
    T201.feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    T201.acid_loading = loading

R201 = system.R201
D = shape.Triangle(0.06, 0.099, 0.12)
@param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='%',
       baseline=0.099, distribution=D)
def set_R201_glucan_conversion(X):
    R201.pretreatment_rxns[0].X = X    

D = shape.Triangle(0.8, 0.9, 0.92)
@param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R201_xylan_conversion(X):
    R201.pretreatment_rxns[4].X = X        

# =============================================================================
# Conversion parameters
# =============================================================================

M301 = system.M301
D = shape.Triangle(0.175, 0.2, 0.25)
@param(name='Enzymatic hydrolysis solid loading', element=M301, kind='coupled', units='%',
       baseline=0.2, distribution=D)
def set_R301_solid_loading(loading):
    M301.solid_loading = loading

D = shape.Triangle(10, 20, 30)
@param(name='Enzyme loading', element=M301, kind='coupled', units='mg/g',
       baseline=20, distribution=D)
def set_R301_enzyme_loading(loading):
    M301.enzyme_loading = loading

# Enzymatic hydrolysis
D = shape.Triangle(0, 24, 56)
@param(name='Enzymatic hydrolysis time', element=R301, kind='coupled', units='hr',
       baseline=24, distribution=D)
def set_R301_saccharification_time(tau):
    R301.tau_saccharification = tau

D = shape.Triangle(0.75, 0.9, 0.948-1e-6)
@param(name='Enzymatic hydrolysis glucan-to-glucose', element=R301, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R301_glucan_conversion(X):
    R301.saccharification_rxns[2].X = X

# Fermentation
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading):
    R301.CSL_loading = loading

R302 = system.R302
# 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14)
D = shape.Triangle(0.9, 0.95, 1-1e-6)
@param(name='Seed train fermentation ratio', element=R302, kind='coupled', units='%',
       baseline=0.95, distribution=D)
def set_ferm_ratio(ratio):
    R302.ferm_ratio = ratio

D = shape.Triangle(0.55, 0.76, 0.93)
@param(name='Lactic acid yield', element=R301, kind='coupled', units='g/g',
       baseline=0.76, distribution=D)
def set_lactic_yield(lactic_yield):
    R301.yield_limit = lactic_yield
    set_yield(lactic_yield, R301, R302)
    # R301_X = R301.cofermentation_rxns.X
    # R302_X = R302.cofermentation_rxns.X
    # R301_X[0] = R301_X[3] = R301.yield_limit = lactic_yield
    # R301_X[1] = R301_X[4] = min(R301_X[1], 1-1e-6-R301_X[0]-R301_X[2])
    # R302_X[0] = R302_X[3] = lactic_yield * R302.ferm_ratio
    # R302_X[1] = R302_X[4] = min(R301_X[1]*R302.ferm_ratio, 1-1e-6-R302_X[0]-R302_X[2])    

D = shape.Triangle(0.33, 0.89, 1.66)
@param(name='Productivity', element=R301, kind='coupled', units='g/L/hr',
       baseline=0.89, distribution=D)
def set_lactic_productivity(productivity):
    R301.productivity = productivity
    R302.productivity = productivity * R302.ferm_ratio

D = shape.Triangle(0.004, 0.07, 0.32)
@param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g',
       baseline=0.07, distribution=D)
def set_acetic_yield(acetic_yield): 
    R301_X = R301.cofermentation_rxns.X
    R302_X = R302.cofermentation_rxns.X
    R301_X[1] = R301_X[4] = min(acetic_yield, 1-1e-6-R301_X[0]-R301_X[2])
    R302_X[1] = R302_X[4] = min(R301_X[1]*R302.ferm_ratio, 1-1e-6-R302_X[0]-R302_X[2])

D = shape.Triangle(0.05, 0.07, 0.1)
@param(name='Inoculum ratio', element=R301, kind='coupled', units='%',
       baseline=0.07, distribution=D)
def set_inoculum_ratio(ratio):
    R301.inoculum_ratio = ratio


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
D = baseline_triangle(1, 0.1)
@param(name='Acidulation time', element=R401, kind='coupled', units='hr',
       baseline=1, distribution=D)
def set_R401_tau(tau):
    R401.tau = tau

R402 = system.R402
D = baseline_triangle(1, 0.1)
@param(name='Esterification conversion factor', element=R402, kind='coupled', units='',
       baseline=1, distribution=D)
def set_R402_conversion_factor(factor):
    R402.X_factor = factor
    
R403 = system.R403
D = baseline_triangle(0.8, 0.1)
@param(name='Hydrolysis conversion', element=R403, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_R403_conversion_factor(X):
    R403.hydrolysis_rxns.X[:] = X
    

# =============================================================================
# Facilities parameters
# =============================================================================

D = baseline_uniform(0.8, 0.1)
@param(name='boiler efficiency', element=CHP, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_boiler_efficiency(efficiency):
    CHP.B_eff = efficiency

# All parameters
parameters = model_full.get_parameters()




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

# %% 

# =============================================================================
# Setup
# =============================================================================

import biosteam as bst
from biosteam.evaluation import Model, Metric
from biorefineries.isobutanol.system import corn_EtOH_IBO_sys, corn_EtOH_IBO_sys_tea, model_specification, unit_groups, unit_groups_dict

IBO_sys = corn_EtOH_IBO_sys
IBO_tea = corn_EtOH_IBO_sys_tea

get_annual_factor = lambda: IBO_tea.operating_hours # hours per year

_kg_per_ton = 907.18474

system_feeds = [i for i in IBO_sys.feeds if i.price]
system_products = [i for i in IBO_sys.products if i.price]

f = IBO_sys.flowsheet

u, s = f.unit, f.stream

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.corn
product_stream = s.ethanol

V406 = u.V406


_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of TAL stream
def get_MSP():
    for i in range(3):
        product_stream.price = IBO_tea.solve_price(product_stream)
    # return product_stream.price*product_stream.F_mass/sum(product_stream.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])
    return product_stream.price 

# Mass flow rate of TAL stream
get_yield = lambda: product_stream.F_mass*get_annual_factor()/1e6
# Purity (%) of TAL in the final product
get_purity = lambda: product_stream.imass['Ethanol']/product_stream.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: product_stream.imol['Ethanol']\
    /(V406.outs[1].imol['Ethanol'])
get_overall_TCI = lambda: IBO_tea.TCI/1e6

get_overall_installed_cost = lambda: IBO_tea.installed_equipment_cost/1e6

# Annual operating cost, note that AOC includes electricity credit
get_overall_AOC = lambda: IBO_tea.AOC/1e6
get_material_cost = lambda: (IBO_tea.material_cost)/1e6
get_overall_FOC = lambda: IBO_tea.FOC/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: IBO_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (IBO_sys.power_utility.production-IBO_sys.power_utility.consumption)
get_electricity_price = lambda: bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*get_electricity_price()*get_annual_factor())/1e6

metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Production rate', get_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product purity', get_purity, '%', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg IBO', 'Biorefinery'),
           Metric('Adjusted production rate', get_adjusted_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product recovery', get_recovery, '%', 'Biorefinery'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $', 'Biorefinery'),
           Metric('Total installed equipment cost', get_overall_installed_cost, '10^6 $', 'Biorefinery'),
           Metric('Annual material cost (incl. boiler ash disposal)', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual operating cost (incl. electricity credit)', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale (excl. electricity)', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Fixed operating cost', get_overall_FOC, '10^6 $/yr', 'Biorefinery'),
           ]

# To see if TEA converges well for each simulation
get_NPV = lambda: IBO_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))

#%% Unit group metrics - absolute
ug_metrics = unit_groups[0].metrics
for mi in range(len(ug_metrics)):
    m = ug_metrics[mi]
    for ug in unit_groups:
        metrics.append(Metric(ug.name, ug.metrics[mi], m.units, m.name))

#%% Metric totals
def metric_total_func_generator(metric_index):
    mname =  unit_groups[0].metrics[metric_index].name
    
    if not mname == 'Operating cost':
        return lambda: sum([ugr.metrics[metric_index]() 
                            for ugr in unit_groups])
    else:
        return lambda: sum([ugr.metrics[metric_index]() 
                            for ugr in unit_groups])\
                        * IBO_tea.operating_hours / 1e6
                        
def metric_total_without_offset_func_generator(metric_index):
    mname =  unit_groups[0].metrics[metric_index].name
    if mname == 'Operating cost':
        return lambda: (sum([ugr.metrics[metric_index]() 
                        for ugr in unit_groups]) + excess_power()*get_electricity_price())\
                        * IBO_tea.operating_hours / 1e6
    elif mname in ('Heating duty', 'Cooling duty'):
        return lambda: sum([ugr.metrics[metric_index]() 
                        for ugr in unit_groups]) - unit_groups_dict['heat exchanger network'].metrics[metric_index]()

for mi in range(len(ug_metrics)):
    m = ug_metrics[mi]
    metrics.append(Metric('Total', 
                          metric_total_func_generator(mi),
                          m.units,
                          m.name))
    if m.name in ('Operating cost', 'Heating duty', 'Cooling duty'):
        metrics.append(Metric('Total before offset',
                              metric_total_without_offset_func_generator(mi), 
                              m.units, 
                              m.name))
metrics.append(Metric('Total', 
                      lambda: IBO_sys.power_utility.rate/1e3,
                      'MW',
                      'Net electricity'))

#%% Unit group metrics - contributions
def metric_fraction_func_generator(u_group, metric_index):
  mname =  u_group.metrics[metric_index].name
  if not mname in ('Operating cost', 'Heating duty', 'Cooling duty'):
      return lambda: u_group.metrics[metric_index]() / max(1e-5, sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups]))
  elif mname == 'Operating cost':
      return lambda: u_group.metrics[metric_index]() / max(1e-5, (sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups]) + excess_power()*get_electricity_price()))
  elif mname in ('Heating duty', 'Cooling duty'):
      return lambda: u_group.metrics[metric_index]() / max(1e-5, (sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups]) - unit_groups_dict['heat exchanger network'].metrics[metric_index]()))

for mi in range(len(ug_metrics)):
    m = ug_metrics[mi]
    for ug in unit_groups:
        metrics.append(Metric(ug.name, 
                              metric_fraction_func_generator(ug, mi),
                              m.units, 
                              'Contributions [%] - ' + m.name,
                              ))

#%% Material cost contributions to total operating cost
def mat_cost_frac_of_op_cost_func_generator(si):
  return lambda: si.cost*IBO_tea.operating_hours/(IBO_tea.AOC +
                                                  excess_power()*get_electricity_price()*IBO_tea.operating_hours)

for i in IBO_sys.feeds:
    if i.price:
        metrics.append(Metric(i.ID, 
                              mat_cost_frac_of_op_cost_func_generator(i),
                              'MM$/y', 
                              'Contributions to total operating cost [%]',
                              ))

#%% Material cost contributions to total material cost
def mat_cost_frac_of_mat_cost_func_generator(si):
  return lambda: si.cost*IBO_tea.operating_hours/IBO_tea.material_cost

for i in IBO_sys.feeds:
    if i.price:
        metrics.append(Metric(i.ID, 
                              mat_cost_frac_of_mat_cost_func_generator(i),
                              'MM$/y', 
                              'Contributions to total material cost [%]',
                              ))


#%% Generate the required namespace
namespace_dict = {}
exclude_from_globals = [
    'search',
    'register',
    'register_safely',
    'discard',
    'clear',
    'mark_safe_to_replace',
    'unmark_safe_to_replace']

namespace_dict.update({k:s.__getitem__(k) for k in s.__dir__() if not k in exclude_from_globals})
namespace_dict.update({k:u.__getitem__(k) for k in u.__dir__() if not k in exclude_from_globals})
namespace_dict['feedstock'] = feedstock
namespace_dict['product_stream'] = product_stream
namespace_dict['IBO_tea'] = namespace_dict['tea'] = IBO_tea
# namespace_dict['spec'] = spec
PowerUtility = bst.PowerUtility
namespace_dict['PowerUtility'] = PowerUtility
# namespace_dict['PD'] = s.PD

#%% 

model = IBO_model = Model(IBO_sys, metrics)

model.specification = model_specification

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

# import numpy as np
import biosteam as bst
from chaospy import distributions as shape
# from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
# from biosteam.evaluation.evaluation_tools import Setter
from biorefineries.oxalic.systems.sugarcane.system_sc_broth import oxalic_sys, oxalic_tea, oxalic_lca, u, s, unit_groups, unit_groups_dict, spec, price, TEA_breakdown, simulate_and_print, theoretical_max_g_oxalic_per_g_glucose, oxalic_chemicals

from biorefineries.oxalic.models.model_utils import EasyInputModel, codify
# get_annual_factor = lambda: oxalic_tea._annual_factor

def create_function(code, namespace_dict, var_to_return='y'):
    def wrapper_fn(statement):
        def f():
            exec(codify(statement), namespace_dict)
            return namespace_dict[var_to_return]
        return f
    function = wrapper_fn(code)
    return function

get_annual_factor = lambda: oxalic_tea.operating_hours # hours per year

_kg_per_ton = 907.18474

system_feeds = [i for i in oxalic_sys.feeds if i.price]
system_products = [i for i in oxalic_sys.products if i.price]

# gypsum = find.stream.gypsum
# system_products.append(gypsum)

baseline_yield, baseline_titer, baseline_productivity =\
    spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity

spec.reactor.neutralization = True

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.sugarcane
product_stream = s.OxalicAcid
# CSL = s.CSL_fresh


R302 = u.R302
R303 = u.R303

BT = u.BT701
# F404 = u.F404

_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of oxalic stream
def get_MSP():
    for i in range(3):
        product_stream.price = oxalic_tea.solve_price(product_stream)
    # return product_stream.price*product_stream.F_mass/sum(product_stream.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])
    return product_stream.price 

# Mass flow rate of oxalic stream
get_yield = lambda: product_stream.F_mass*get_annual_factor()/1e6
# Purity (%) of oxalic in the final product
get_purity = lambda: product_stream.imass['OxalicAcid']/product_stream.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()


get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: product_stream.imol['OxalicAcid']\
    /(R302.outs[1].imol['OxalicAcid'])
get_overall_TCI = lambda: oxalic_tea.TCI/1e6

get_overall_installed_cost = lambda: oxalic_tea.installed_equipment_cost/1e6

# Annual operating cost, note that AOC includes electricity credit
get_overall_AOC = lambda: oxalic_tea.AOC/1e6
get_material_cost = lambda: (oxalic_tea.material_cost +
                             abs(BT.ash_disposal_price*
                                 BT.ash_disposal.F_mass*
                                 oxalic_tea.operating_hours))/1e6
get_overall_FOC = lambda: oxalic_tea.FOC/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: oxalic_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (oxalic_sys.power_utility.production-oxalic_sys.power_utility.consumption)
get_electricity_price = lambda: bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*get_electricity_price()*get_annual_factor())/1e6

get_T620_cost = lambda: oxalic_sys.flowsheet.T620.installed_cost/1e6
metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Production rate', get_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product purity', get_purity, '%', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg AA', 'Biorefinery'),
           Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product recovery', get_recovery, '%', 'Biorefinery'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $', 'Biorefinery'),
           Metric('Total installed equipment cost', get_overall_installed_cost, '10^6 $', 'Biorefinery'),
           Metric('Annual material cost (incl. boiler ash disposal)', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual operating cost (incl. electricity credit)', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale (excl. electricity)', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Fixed operating cost', get_overall_FOC, '10^6 $/yr', 'Biorefinery'),
           Metric('OA storage cost', get_T620_cost, '10^6', 'Biorefinery')
           ]

# To see if TEA converges well for each simulation
get_NPV = lambda: oxalic_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))


# metrics_labels_dict = {
#     'Installed cost':(0, '10^6 $'), 
#     'Material cost':(4,'USD/h'), 
#     'Cooling duty':(1,'GJ/h'), 
#     'Heating duty':(2,'GJ/h'), 
#     'Electricity usage':(3, 'MW'), 
#     }

# for m, u_i in metrics_labels_dict.items():
#     for ug in unit_groups:
#         metrics.append(Metric(ug.name, ug.metrics[u_i[0]], u_i[1], m))

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
                        * oxalic_tea.operating_hours / 1e6
                        
def metric_total_without_offset_func_generator(metric_index):
    mname =  unit_groups[0].metrics[metric_index].name
    if mname == 'Operating cost':
        return lambda: (sum([ugr.metrics[metric_index]() 
                        for ugr in unit_groups]) + excess_power()*get_electricity_price())\
                        * oxalic_tea.operating_hours / 1e6
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
                      lambda: oxalic_sys.power_utility.rate/1e3,
                      'MW',
                      'Net electricity'))

# metrics.append(Metric('Total', 
#                       lambda: unit_groups_dict['natural gas (for product drying)'].metrics[1]() * oxalic_tea.operating_hours/1e6,
#                       'MM$/y',
#                       'Natural gas (product drying) material cost'))

#%% Unit group metrics - contributions
def metric_fraction_func_generator(u_group, metric_index):
  mname =  u_group.metrics[metric_index].name
  if not mname in ('Operating cost', 'Heating duty', 'Cooling duty'):
      return lambda: u_group.metrics[metric_index]() / sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups])
  elif mname == 'Operating cost':
      return lambda: u_group.metrics[metric_index]() / (sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups]) + excess_power()*get_electricity_price())
  elif mname in ('Heating duty', 'Cooling duty'):
      return lambda: u_group.metrics[metric_index]() / (sum([ugr.metrics[metric_index]() 
                                                            for ugr in unit_groups]) - unit_groups_dict['heat exchanger network'].metrics[metric_index]())

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
  return lambda: si.cost*oxalic_tea.operating_hours/(oxalic_tea.AOC +
                                                  excess_power()*get_electricity_price()*oxalic_tea.operating_hours)

for i in oxalic_sys.feeds:
    if i.price:
        metrics.append(Metric(i.ID, 
                              mat_cost_frac_of_op_cost_func_generator(i),
                              'MM$/y', 
                              'Contributions to total operating cost [%]',
                              ))

#%% Material cost contributions to total material cost
def mat_cost_frac_of_mat_cost_func_generator(si):
  return lambda: si.cost*oxalic_tea.operating_hours/oxalic_tea.material_cost

for i in oxalic_sys.feeds:
    if i.price:
        metrics.append(Metric(i.ID, 
                              mat_cost_frac_of_mat_cost_func_generator(i),
                              'MM$/y', 
                              'Contributions to total material cost [%]',
                              ))

#%% LCA - absolute impacts
get_GWP_before_offset = lambda: oxalic_lca.GWP - oxalic_lca.net_electricity_GWP
get_other_materials_GWP = lambda: oxalic_lca.material_GWP - oxalic_lca.material_GWP_breakdown['AceticAcid'] - oxalic_lca.material_GWP_breakdown['CSL'] - oxalic_lca.material_GWP_breakdown['DAP'] - oxalic_lca.material_GWP_breakdown['CH4'] # all materials other than feedstock, CSL, acetate, and CH4

# IPCC 2013 GWP100a
metrics.append(Metric('Total GWP100a', lambda: oxalic_lca.GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('Total GWP100a excl. net electricity', get_GWP_before_offset, 'kg-CO2-eq/kg', 'Biorefinery'))
# metrics.append(Metric('GWP100a - Heating demand', lambda: oxalic_lca.heating_demand_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
# metrics.append(Metric('GWP100a - Cooling demand', lambda: oxalic_lca.cooling_demand_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Net electricity', lambda: oxalic_lca.net_electricity_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Direct non-biogenic emissions', lambda: oxalic_lca.direct_emissions_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))

metrics.append(Metric('GWP100a - Feedstock (FGHTP) ', lambda: oxalic_lca.FGHTP_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Materials (except feedstock) ', lambda: oxalic_lca.material_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))

for k in oxalic_lca.material_GWP_breakdown.keys():
    code = f"y = oxalic_lca.material_GWP_breakdown['{k}']"
    f = create_function(code, {'oxalic_lca': oxalic_lca})
    metrics.append(Metric(f'GWP100a - Materials breakdown - {k}', f, 'kg-CO2-eq/kg', 'Biorefinery'))
    

# FEC
get_FEC_before_offset = lambda: oxalic_lca.FEC - oxalic_lca.net_electricity_FEC
get_other_materials_FEC = lambda: oxalic_lca.material_FEC - oxalic_lca.material_FEC_breakdown['AceticAcid'] - oxalic_lca.material_FEC_breakdown['CH4'] - oxalic_lca.material_FEC_breakdown['CSL'] - oxalic_lca.material_FEC_breakdown['DAP']

metrics.append(Metric('Total FEC', lambda: oxalic_lca.FEC, 'MJ/kg', 'Biorefinery'))
metrics.append(Metric('Total FEC excl. net electricity', get_FEC_before_offset, 'MJ/kg', 'Biorefinery'))
# metrics.append(Metric('FEC - Heating demand', lambda: oxalic_lca.heating_demand_FEC, 'MJ/kg', 'Biorefinery'))
# metrics.append(Metric('FEC - Cooling demand', lambda: oxalic_lca.cooling_demand_FEC, 'MJ/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Net electricity', lambda: oxalic_lca.net_electricity_FEC, 'MJ/kg', 'Biorefinery'))

metrics.append(Metric('FEC - Feedstock (FGHTP) ', lambda: oxalic_lca.feedstock_FEC, 'MJ/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Materials (except feedstock) ', lambda: oxalic_lca.material_FEC, 'MJ/kg', 'Biorefinery'))

for k in oxalic_lca.material_FEC_breakdown.keys():
    code = f"y = oxalic_lca.material_FEC_breakdown['{k}']"
    f = create_function(code, {'oxalic_lca': oxalic_lca})
    metrics.append(Metric(f'FEC - Materials breakdown - {k}', f, 'MJ/kg', 'Biorefinery'))
    
#%% LCA - % contributions relative to total positive impacts

# IPCC 2013 GWP100a

metrics.append(Metric('GWP100a % - Net electricity', lambda: oxalic_lca.net_electricity_GWP/get_GWP_before_offset(), '%', 'Biorefinery'))
metrics.append(Metric('GWP100a %- Direct non-biogenic emissions', lambda: oxalic_lca.direct_emissions_GWP/get_GWP_before_offset(), '%', 'Biorefinery'))

metrics.append(Metric('GWP100a % - Feedstock (FGHTP) ', lambda: oxalic_lca.FGHTP_GWP/get_GWP_before_offset(), '%', 'Biorefinery'))
metrics.append(Metric('GWP100a % - Materials (except feedstock and BT.natural_gas) ', lambda: oxalic_lca.material_GWP/get_GWP_before_offset(), '%', 'Biorefinery'))

for k in oxalic_lca.material_GWP_breakdown.keys():
    code = f"y = oxalic_lca.material_GWP_breakdown['{k}'] / get_GWP_before_offset()"
    f = create_function(code, {'oxalic_lca': oxalic_lca, 'get_GWP_before_offset': get_GWP_before_offset})
    metrics.append(Metric(f'GWP100a % - Materials breakdown - {k}', f, '%', 'Biorefinery'))
    
# FEC

metrics.append(Metric('FEC % - Net electricity', lambda: oxalic_lca.net_electricity_FEC/get_FEC_before_offset(), '%', 'Biorefinery'))

metrics.append(Metric('FEC % - Feedstock (FGHTP) ', lambda: oxalic_lca.feedstock_FEC/get_FEC_before_offset(), '%', 'Biorefinery'))
metrics.append(Metric('FEC % - Materials (except feedstock) ', lambda: oxalic_lca.material_FEC/get_FEC_before_offset(), 'MJ/kg', 'Biorefinery'))

for k in oxalic_lca.material_FEC_breakdown.keys():
    code = f"y = oxalic_lca.material_FEC_breakdown['{k}'] / get_FEC_before_offset()"
    f = create_function(code, {'oxalic_lca': oxalic_lca, 'get_FEC_before_offset': get_GWP_before_offset})
    metrics.append(Metric(f'FEC % - Materials breakdown - {k}', f, '%', 'Biorefinery'))
    
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
namespace_dict['oxalic_tea'] = namespace_dict['tea'] = oxalic_tea
namespace_dict['spec'] = spec
PowerUtility = bst.PowerUtility
namespace_dict['PowerUtility'] = PowerUtility
# namespace_dict['PD'] = s.PD
namespace_dict['theoretical_max_g_oxalic_per_g_glucose'] = theoretical_max_g_oxalic_per_g_glucose

#%% 
model = oxalic_model = EasyInputModel(oxalic_sys, metrics, namespace_dict=namespace_dict)
    
#%% Bugfix barrage
baseline_spec = {'spec_1': spec.baseline_yield,
                 'spec_2': spec.baseline_titer,
                 'spec_3': spec.baseline_productivity,}

system=model._system
def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_specifications(**baseline_spec)
    spec.set_production_capacity(spec.desired_annual_production)
    # system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    spec.set_production_capacity(spec.desired_annual_production)
    # system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    spec.set_production_capacity(spec.desired_annual_production)
    # system.simulate()
    
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                # breakpoint()
                raise e
###############################

#%% Model specification
pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
pre_fermenter_units_path.reverse()
def model_specification():
    try:
        for i in pre_fermenter_units_path: i.simulate()
        # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        spec.set_production_capacity(spec.desired_annual_production)
        # system.simulate()
        # model._system.simulate()
    

    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # breakpoint()
        # raise e
        if 'sugar concentration' in str_e:
            # flowsheet('OxalicAcid').F_mass /= 1000.
            raise e
        else:
            run_bugfix_barrage()
            
model.specification = model_specification



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

import numpy as np
import biosteam as bst
from chaospy import distributions as shape
from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools.parameter import Setter
from biorefineries.HP.system_light_lle_vacuum_distillation import process_groups_dict, HP_tea, HP_no_BT_tea,\
flowsheet, get_GWP, get_FEC, CFs, spec, process_groups, process_groups_dict,\
get_material_cost_breakdown, get_GWP_by_ID, get_FEC_by_ID, get_ng_GWP, get_ng_FEC,\
get_FGHTP_GWP, get_feedstock_FEC, get_electricity_demand_non_cooling_GWP, get_electricity_demand_non_cooling_FEC,\
get_direct_emissions_GWP, get_net_electricity_GWP, get_net_electricity_FEC,\
get_heating_demand_GWP, get_heating_demand_FEC,\
get_cooling_demand_GWP, get_cooling_demand_FEC, get_heating_demand_VOC, get_cooling_demand_VOC, get_electricity_demand_non_cooling_VOC, get_VOC
from warnings import warn

find.set_flowsheet(flowsheet)
# get_annual_factor = lambda: HP_no_BT_tea._annual_factor
get_annual_factor = lambda: 1.
_kg_per_ton = 907.18474
_feedstock_factor = _kg_per_ton / 0.8

HP_sys = find.system.HP_sys
BT_sys = find.system.BT_sys
HXN = find.unit.HXN

system_feeds = [i for i in HP_sys.feeds if i.price] + \
    [i for i in BT_sys.feeds if i.price]
system_products = [i for i in HP_sys.products if i.price] + \
    [i for i in BT_sys.products if i.price]
    
gypsum = find.stream.gypsum
system_products.append(gypsum)


# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

# Minimum selling price of AA stream
def get_MSP():
    for i in range(3):
        AA.price = HP_tea.solve_price(AA)
    return AA.price

# Mass flow rate of HP stream
AA = find.unit.T606_P-0
feedstock = find.stream.feedstock
get_yield = lambda: AA.F_mass*get_annual_factor()/1e6
# Purity (%) of HP in the final product
get_purity = lambda: AA.imass['AA']/AA.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
R301 = find.unit.R301
get_recovery = lambda: AA.imol['HP'] \
    /(R302.outs[0].imol['HP']+2*R302.outs[0].imol['CalciumLactate'])
get_overall_TCI = lambda: HP_tea.TCI/1e6
# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: HP_tea.AOC/1e6
get_material_cost = lambda: HP_tea.material_cost/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: HP_tea.sales/1e6
# System power usage, individual unit power usage should be positive
BT = find.unit.BT
# excess_power = lambda: BT.electricity_generated
excess_power = lambda: BT.power_utility.production
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6


#%% 
class Metrics(list):
    def __init__(self, metrics):
        self.extend(metrics)
        
    def append(self, metric):
        # self_copy = copy.deepcopy(self)
        repeated = False
        for i in self:
            if i.index == metric.index:
                repeated = True
                break
        if repeated:
            warn("Metric {metric.index} already exists in Metrics object. Second instance has been deleted.")
        else:
            super().append(metric)
        
    def extend(self,metrics):
        for i in metrics:
            self.append(i)
#%%

metrics = Metrics([Metric('Minimum selling price', get_MSP, '$/kg'),
           Metric('Product yield', get_yield, '10^6 kg/yr'),
           Metric('Product purity', get_purity, '%'),
            Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg'),
           Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr'),
           Metric('Product recovery', get_recovery, '%'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr')
           ])

# %% Breakdowns by process groups
def get_group_heating_demand(group):
    return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow>0.]) for unit in group.units])

def get_group_cooling_demand(group):
    return sum([sum([hu.duty for hu in unit.heat_utilities if hu.duty*hu.flow<0.]) for unit in group.units])

# Heating duty

metrics.extend((Metric('feedstock_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['feedstock_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('pretreatment_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['pretreatment_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('conversion_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['conversion_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('separation_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['separation_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('WWT_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['WWT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('HXN_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['HXN_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('BT_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['BT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('CT_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['CT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('facilities_no_hu_group - heating demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow>0.]) for unit in \
                                          process_groups_dict['facilities_no_hu_group'].units])/AA.F_mass,
                       'MJ/kg'),))

# Cooling duty

metrics.extend((Metric('feedstock_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['feedstock_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('pretreatment_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['pretreatment_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('conversion_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['conversion_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('separation_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['separation_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('WWT_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['WWT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('HXN_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['HXN_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('BT_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['BT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('CT_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['CT_group'].units])/AA.F_mass,
                       'MJ/kg'),))

metrics.extend((Metric('facilities_no_hu_group - cooling demand',
                       lambda: 0.001*sum([sum([hu.duty for hu in unit.heat_utilities \
                                               if hu.duty*hu.flow<0.]) for unit in \
                                          process_groups_dict['facilities_no_hu_group'].units])/AA.F_mass,
                       'MJ/kg'),))  


# Installed equipment cost

metrics.extend((Metric('feedstock_group - installed equipment cost',
                       lambda:process_groups_dict['feedstock_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('pretreatment_group - installed equipment cost',
                       lambda:process_groups_dict['pretreatment_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('conversion_group - installed equipment cost',
                       lambda:process_groups_dict['conversion_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('separation_group - installed equipment cost',
                       lambda:process_groups_dict['separation_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('WWT_group - installed equipment cost',
                       lambda:process_groups_dict['WWT_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('HXN_group - installed equipment cost',
                       lambda:process_groups_dict['HXN_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('BT_group - installed equipment cost',
                       lambda:process_groups_dict['BT_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('CT_group - installed equipment cost',
                       lambda:process_groups_dict['CT_group'].get_installed_cost(),
                       '10^6 $'),))

metrics.extend((Metric('facilities_no_hu_group - installed equipment cost',
                       lambda:process_groups_dict['facilities_no_hu_group'].get_installed_cost(),
                       '10^6 $'),))  

# Power utility demand

metrics.extend((Metric('feedstock_group - power utility demand',
                       lambda:process_groups_dict['feedstock_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('pretreatment_group - power utility demand',
                       lambda:process_groups_dict['pretreatment_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('conversion_group - power utility demand',
                       lambda:process_groups_dict['conversion_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('separation_group - power utility demand',
                       lambda:process_groups_dict['separation_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('WWT_group - power utility demand',
                       lambda:process_groups_dict['WWT_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('HXN_group - power utility demand',
                       lambda:process_groups_dict['HXN_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('BT_group - power utility demand',
                       lambda:process_groups_dict['BT_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('CT_group - power utility demand',
                       lambda:process_groups_dict['CT_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))

metrics.extend((Metric('facilities_no_hu_group - power utility demand',
                       lambda:process_groups_dict['facilities_no_hu_group'].get_electricity_consumption()/AA.F_mass,
                       'MW/kg'),))  

# Material cost

metrics.extend((Metric('feedstock_group - material cost',
                       lambda:get_material_cost_breakdown()['feedstock_group'],
                       '$/kg'),))

metrics.extend((Metric('pretreatment_group - material cost',
                       lambda:get_material_cost_breakdown()['pretreatment_group'],
                       '$/kg'),))

metrics.extend((Metric('conversion_group - material cost',
                       lambda:get_material_cost_breakdown()['conversion_group'],
                       '$/kg'),))

metrics.extend((Metric('separation_group - material cost',
                       lambda:get_material_cost_breakdown()['separation_group'],
                       '$/kg'),))

metrics.extend((Metric('WWT_group - material cost',
                       lambda:get_material_cost_breakdown()['WWT_group'],
                       '$/kg'),))

metrics.extend((Metric('HXN_group - material cost',
                       lambda:get_material_cost_breakdown()['HXN_group'],
                       '$/kg'),))

metrics.extend((Metric('BT_group - material cost',
                       lambda:get_material_cost_breakdown()['BT_group'],
                       '$/kg'),))

metrics.extend((Metric('CT_group - material cost',
                       lambda:get_material_cost_breakdown()['CT_group'],
                       '$/kg'),))

metrics.extend((Metric('facilities_no_hu_group - material cost',
                       lambda:get_material_cost_breakdown()['facilities_no_hu_group'],
                       '$/kg'),))
# # =============================================================================
# # Capital cost breakdown
# # =============================================================================

# def get_installed_equipment_cost(group):
#     # return lambda: sum(i.installed_equipment_cost for i in process_groups_dict[system])/1e6
#     return process_groups_dict[group].get_installed_cost()/AA.F_mass
# for system in process_groups_dict.keys():
#     # if system == 'feedstock_sys': continue
#     metrics.extend(
#         (Metric(system, get_installed_equipment_cost(system), '10^6 $', 'Installed cost'),))

# # All checks should be ~0
# check_installed_equipment_cost = \
#     lambda: sum(get_installed_equipment_cost(system)() 
#                 for system in process_groups_dict.keys()) - HP_tea.installed_equipment_cost/1e6
# metrics.extend((Metric('Check', check_installed_equipment_cost, '10^6 $', 'Installed cost'),))


# =============================================================================
# Material cost breakdown
# =============================================================================

def get_material_cost(feed):
    return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
for feed in system_feeds:
    metrics.extend((Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'),))
fermentation_lime = find.stream.fermentation_lime
FGD_lime = find.stream.FGD_lime
# get_fermentation_lime_ratio = lambda: fermentation_lime.imol['Lime'] \
#     / (fermentation_lime.imol['Lime']+FGD_lime.imol['Lime']) 
# S601 = find.unit.S601
# get_separation_sulfuric_acid_ratio = lambda: S601.outs[1].imol['H2SO4']/S601.ins[0].imol['H2SO4']
# get_separation_sulfuric_acid_ratio = 
check_material_cost = lambda: sum(get_material_cost(feed)()
                                  for feed in system_feeds) - HP_tea.material_cost/1e6

# metrics.extend((
#     Metric('Fermentation lime ratio', get_fermentation_lime_ratio, 
#            '%', 'Material cost'),
#     # Metric('Separation sulfuric acid ratio', get_separation_sulfuric_acid_ratio, 
#     #        '%', 'Material cost'),
#     Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')))

def get_product_sale(stream):
    return lambda: stream.price*stream.F_mass*get_annual_factor()/1e6
for product in system_products:
    metrics.extend((Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'),))
check_product_sale= \
    lambda: sum(get_product_sale(product)() for product in system_products) \
        - HP_tea.sales/1e6
metrics.extend((Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'),))


# # =============================================================================
# # Heating demand breakdown (positive if needs heating)
# # =============================================================================

# # get_system_heating_demand = lambda: BT.system_heating_demand*get_annual_factor()/1e9
# get_system_heating_demand = lambda: sum([i.duty for i in BT.steam_utilities])*get_annual_factor()/1e9
# # get_pretreatment_steam_heating_demand = lambda: BT.side_streams_lps.duty*get_annual_factor()/1e9
# get_pretreatment_steam_heating_demand = lambda: unit_groups[0].get_heating_duty()*get_annual_factor()/1e3

# get_HXN_heating_demand = lambda: sum(i.duty for i in HXN.heat_utilities 
#                                     if i.duty*i.cost>0)*get_annual_factor()/1e9
# get_BT_heating_demand = lambda: sum(i.duty for i in BT.heat_utilities 
#                                     if i.duty*i.cost>0)*get_annual_factor()/1e9

# def get_heating_demand(system):
#     heat_utilities = sum([i.heat_utilities for i in process_groups_dict[system]], ())
#     return lambda: sum([i.duty for i in heat_utilities
#                         if i.duty*i.cost>0])*get_annual_factor()/1e9

# for system in process_groups_dict.keys():
#     if system in ('feedstock_sys', 'HXN', 'BT'): continue
#     # The only heating demand for the pretreatment system is the heat needed to
#     # generate the side steam
#     if system == 'pretreatment_sys':
#         metrics.extend((Metric(system, get_pretreatment_steam_heating_demand, '10^9 kJ/yr', 
#                                 'Heating demand'),))
#     else: metrics.extend((Metric(system, get_heating_demand(system), '10^9 kJ/yr', 
#                                   'Heating demand'),))

# check_heating_demand = \
#     lambda: sum((get_heating_demand(system)() for system in process_groups_dict.keys()), 
#                 get_pretreatment_steam_heating_demand())
                                            
# metrics.extend((
#     Metric('HXN', get_HXN_heating_demand, '10^9 kJ/yr', 'Heating demand'),
#     Metric('BT', get_BT_heating_demand, '10^9 kJ/yr', 'Heating demand'),    
#     Metric('Sum', get_system_heating_demand, '10^9 kJ/yr', 'Heating demand'),
#     Metric('Check', check_heating_demand, '10^9 kJ/yr', 'Heating demand')
#     ))

# # =============================================================================
# # Cooling demand breakdown (negative if needs cooling)
# # =============================================================================

# CT = find.unit.CT
# get_system_cooling_water_duty = lambda: CT.system_cooling_water_duty*get_annual_factor()/1e9
# get_HXN_cooling_demand = lambda: sum(i.duty for i in HXN.heat_utilities 
#                                     if i.duty*i.cost<0)*get_annual_factor()/1e9
# get_CT_cooling_demand = lambda: sum(i.duty for i in CT.heat_utilities 
#                                     if i.duty*i.cost<0)*get_annual_factor()/1e9

# def get_cooling_demand(system):
#     heat_utilities = sum((i.heat_utilities for i in process_groups_dict[system]), ())
#     return lambda: sum([i.duty for i in heat_utilities
#                         if i.duty*i.cost<0])*get_annual_factor()/1e9

# for system in process_groups_dict.keys():
#     if system in ('feedstock_sys', 'HXN', 'CT'): continue
#     else: metrics.extend((Metric(system, get_cooling_demand(system),
#                                   '10^9 kJ/yr', 'Cooling demand'),))

# check_cooling_demand = \
#     lambda: sum(get_cooling_demand(system)() for system in process_groups_dict.keys())

# metrics.extend((
#     Metric('HXN', get_HXN_cooling_demand, '10^9 kJ/yr', 'Cooling demand'),
#     Metric('CT', get_CT_cooling_demand, '10^9 kJ/yr', 'Cooling demand'),    
#     Metric('Sum', get_system_cooling_water_duty, '10^9 kJ/yr', 'Cooling demand'),
#     Metric('Check', check_cooling_demand, '10^9 kJ/yr', 'Cooling demand')
#     ))

# # =============================================================================
# # Power demand breakdown (positive if using power)
# # =============================================================================

# get_system_power_demand = lambda: sum(i.power_utility.rate for i in HP_sys.units
#                                       if i.power_utility)

# def get_power_demand(system):
#     power_utilities = [i.power_utility for i in process_groups_dict[system]]
#     return lambda: sum(i.rate for i in power_utilities)

# for system in process_groups_dict.keys():
#     if system == 'feedstock_sys': continue
#     metrics.extend((Metric(system, get_power_demand(system), 'kW', 'Power demand'),))

# check_power_demand = lambda: sum(get_power_demand(system)()
#                                  for system in process_groups_dict.keys()) - get_system_power_demand()
# metrics.extend((
#     Metric('Sum', get_system_power_demand, 'kW', 'Power demand'),
#     Metric('Check', check_power_demand, 'kW', 'Power demand')
#     ))

# # =============================================================================
# # Utility cost breakdown (including heating, cooling, and power)
# # =============================================================================

# get_system_utility_cost = lambda: HP_tea.utility_cost/1e6

# def get_utility_cost(system):
#     return lambda: sum(i.utility_cost for i in process_groups_dict[system])*get_annual_factor()/1e6

# for system in process_groups_dict.keys():
#     if system == 'feedstock_sys': continue
#     metrics.extend((Metric(system, get_utility_cost(system), '10^6 $/yr', 'Utility cost'),))

# check_utility_cost = \
#     lambda: sum(get_utility_cost(system)() for system in process_groups_dict.keys()) \
#         - get_system_utility_cost()

# metrics.extend((
#     Metric('Sum', get_system_utility_cost, '10^6 $/yr', 'Utility cost'),
#     Metric('Check', check_utility_cost, '10^6 $/yr', 'Utility cost')
#     ))


# To see if TEA converges well for each simulation
get_NPV = lambda: HP_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))

# To check HXN energy balance error
# metrics.append(Metric('HXN energy balance error', lambda: HXN.energy_balance_percent_error))

metrics.extend((Metric('HXN energy balance error', lambda: HXN.energy_balance_percent_error, '%', 'TEA'), ))

##### LCA #####
metrics.extend((
    Metric('Total GWP', get_GWP, 'kg CO2-eq/kg', 'LCA'),
    Metric('Total FEC', get_FEC, 'MJ/kg', 'LCA')
    ))

# Material GWP
metrics.extend((Metric('GWP - H2SO4',
                       lambda:get_GWP_by_ID('H2SO4'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - NaOH',
                       lambda:get_GWP_by_ID('NaOH'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - AmmoniumHydroxide',
                       lambda:get_GWP_by_ID('AmmoniumHydroxide'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - CalciumDihydroxide',
                       lambda:get_GWP_by_ID('CalciumDihydroxide'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - Hexanol',
                       lambda:get_GWP_by_ID('Hexanol'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - Enzyme',
                       lambda:get_GWP_by_ID('Enzyme'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - TiO2',
                       lambda:get_GWP_by_ID('TiO2'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('GWP - CSL',
                       lambda:get_GWP_by_ID('CSL'),
                       'kg CO2/kg', 'LCA'),))

# Material FEC
metrics.extend((Metric('FEC - H2SO4',
                       lambda:get_FEC_by_ID('H2SO4'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - NaOH',
                       lambda:get_FEC_by_ID('NaOH'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - AmmoniumHydroxide',
                       lambda:get_FEC_by_ID('AmmoniumHydroxide'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - CalciumDihydroxide',
                       lambda:get_FEC_by_ID('CalciumDihydroxide'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - Hexanol',
                       lambda:get_FEC_by_ID('Hexanol'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - Enzyme',
                       lambda:get_FEC_by_ID('Enzyme'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - TiO2',
                       lambda:get_FEC_by_ID('TiO2'),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - CSL',
                       lambda:get_FEC_by_ID('CSL'),
                       'kg CO2/kg', 'LCA'),))
# Natural gas
metrics.extend((Metric('GWP - natural gas',
                       lambda:get_ng_GWP(),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - natural gas',
                       lambda:get_ng_FEC(),
                       'MJ/kg', 'LCA'),))

# Natural gas
metrics.extend((Metric('GWP - natural gas',
                       lambda:get_ng_GWP(),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - natural gas',
                       lambda:get_ng_FEC(),
                       'MJ/kg', 'LCA'),))

# Electricity
metrics.extend((Metric('GWP - electricity, net',
                       lambda:get_net_electricity_GWP(),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - electricity, net',
                       lambda:get_net_electricity_FEC(),
                       'MJ/kg', 'LCA'),))

# Feedstock growth, harvesting, transportation, and preprocessing
metrics.extend((Metric('GWP - Feedstock GHTP',
                       lambda:get_FGHTP_GWP(),
                       'kg CO2/kg', 'LCA'),))
metrics.extend((Metric('FEC - Feedstock GHTP',
                       lambda:get_feedstock_FEC(),
                       'MJ/kg', 'LCA'),))

# Direct non-biogenic emissions GWP
metrics.extend((Metric('GWP - Other direct non-bio emmissions',
                       lambda:get_direct_emissions_GWP(),
                       'kg CO2/kg', 'LCA'),))

metrics.extend((Metric('GWP - Other direct non-bio emmissions',
                       lambda:get_direct_emissions_GWP(),
                       'kg CO2/kg', 'LCA'),))

# Demand LCA contributions
metrics.extend((Metric('cGWP - System heating demand',
                       lambda:get_heating_demand_GWP()/get_GWP(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cGWP - System cooling demand',
                       lambda:get_cooling_demand_GWP()/get_GWP(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cGWP - System non-cooling electricity demand',
                       lambda:get_electricity_demand_non_cooling_GWP()/get_GWP(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cFEC - System heating demand',
                       lambda:get_heating_demand_FEC()/get_FEC(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cFEC - System cooling demand',
                       lambda:get_cooling_demand_FEC()/get_FEC(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cFEC - System non-cooling electricity demand',
                       lambda:get_electricity_demand_non_cooling_FEC()/get_FEC(),
                       'frac', 'LCA'),))

# Demand TEA contributions
metrics.extend((Metric('cVOC - System heating demand',
                       lambda:get_heating_demand_VOC()/get_VOC(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cVOC - System cooling demand',
                       lambda:get_cooling_demand_VOC()/get_VOC(),
                       'frac', 'LCA'),))

metrics.extend((Metric('cVOC - System non-cooling electricity demand',
                       lambda:get_electricity_demand_non_cooling_VOC()/get_VOC(),
                       'frac', 'LCA'),))
# %% 

# =============================================================================
# Construction base model
# =============================================================================

HP_model = Model(HP_sys, metrics)
param = HP_model.parameter

def baseline_uniform(baseline, ratio):
    lb, ub = baseline*(1-ratio), baseline*(1+ratio)
    if lb > ub: ub, lb = lb, ub
    return shape.Uniform(lb, ub)

def baseline_triangle(baseline, ratio):
    lb, mid, ub = baseline*(1-ratio), baseline, baseline*(1+ratio)
    if lb > ub: ub, lb = lb, ub
    return shape.Triangle(lb, mid, ub)


D = baseline_uniform(1, 0.1)
@param(name='Blank parameter', element=feedstock, kind='coupled', units='',
       baseline=1, distribution=D)
def set_blank_parameter(anything):
    # This does nothing
    feedstock.T = feedstock.T
    
# =============================================================================
# TEA parameters
# =============================================================================

# U101 = find.unit.U101
# D = baseline_uniform(2205, 0.1)
# @param(name='Flow rate', element=U101, kind='coupled', units='U.S. ton/day',
#        baseline=2205, distribution=D)
# def set_feedstock_flow_rate(rate):
#     U101.feedstock_flow_rate = rate

D = shape.Triangle(0.84, 0.9, 0.96)
@param(name='Plant uptime', element='TEA', kind='isolated', units='%',
       baseline=0.9, distribution=D)
def set_operating_days(uptime):
    HP_tea.operating_days = 365. * uptime


# D = baseline_triangle(1, 0.25)
# @param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
#         baseline=1, distribution=D)
# def set_TCI_ratio(new_ratio):
#     old_ratio = HP_no_BT_tea._TCI_ratio_cached
#     for unit in HP_sys.units:
#         if hasattr(unit, 'cost_items'):
#             for item in unit.cost_items:
#                 unit.cost_items[item].cost /= old_ratio
#                 unit.cost_items[item].cost *= new_ratio
#     HP_no_BT_tea._TCI_ratio_cached = new_ratio


feedstock = find('feedstock')
D = shape.Triangle(60, 71.3, 83.7)
@param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-ton',
       baseline=71.3, distribution=D)
def set_feedstock_price(price):
    feedstock.price = price / _feedstock_factor


# Impactful parameters are set to triangular distribution based on literature,
# less important ones are set to ±10% of baseline value
special_price = {
#     stream           distribution     min  mid   max
    # 'feedstock':        ('Triangle',   (60, 71.3, 83.7)),
    # 'CSL_fresh':        ('Uniform',   (0.0673,    0.112)),
    'lime_fresh':       ('Triangle',   (0.160, 0.262, 0.288)),
    # 'ethanol_fresh':    ('Triangle',  (0.460,     0.978)),
    # 'hexanol_fresh':    ('Uniform',   (0.551*0.9,     0.551*1.1)),
    # 'natural_gas':      ('Triangle',  (0.198,     0.304)),
    'gypsum':           ('Uniform',   (-0.0288,   0.00776))
    }

# Prices for boiler_chems, baghouse_bag, and cooling_tower_chems are not included
# as they are tied to BT/CT duties
default_price_streams = ('sulfuric_acid_fresh',
                         'makeup_TiO2_catalyst', 'ammonia_fresh', 'enzyme', 
                         'system_makeup_water', 'aerobic_caustic', 'ash', 'hexanol_fresh')

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
    D = baseline_triangle(baseline, 0.1)
    add_stream_price_param(stream, D)

D = shape.Triangle(0.067, 0.070, 0.073)
@param(name='Electricity price', element='TEA', kind='isolated', units='$/kWh',
       baseline=0.070, distribution=D)
def set_electricity_price(price): 
    bst.PowerUtility.price = price

BT = flowsheet('BT')
natural_gas_price = BT.natural_gas_price
D = shape.Triangle(natural_gas_price*0.9, natural_gas_price, natural_gas_price*1.1)
@param(name='Natural gas price', element='TEA', kind='isolated', units='$/kWh',
       baseline=natural_gas_price, distribution=D)
def set_natural_gas_price(price): 
    BT.natural_gas_price = price
    

# HP_sys._cached_TCI_ratio = 1.
# new_ratio = 1.
D = baseline_triangle(1., 0.25)
@param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
        baseline=1., distribution=D)
def set_TCI_ratio(ratio): 
    for unit in HP_sys.units:
        if hasattr(unit, 'cost_items'):
            for item in unit.cost_items:
                unit.cost_items[item].cost *= ratio/HP_no_BT_tea._TCI_ratio_cached
                HP_no_BT_tea._TCI_ratio_cached = ratio

# =============================================================================
# LCA parameters
# =============================================================================

D = shape.Uniform(0.09646, 0.12894) # see Feedstock_impacts_YL
@param(name='Feedstock GHTP GWP100 CF', element='LCA', kind='isolated', units='kg-CO2eq. / dry-kg',
       baseline=0.10945, distribution=D)
def set_feedstock_GWP_CF(CF):
    CFs['GWP_CFs']['FGHTP Corn stover'] = CF
    
D = shape.Uniform(1.32576, 1.75741) # see Feedstock_impacts_YL
@param(name='Feedstock GHTP FEC CF', element='LCA', kind='isolated', units='MJeq. / dry-kg',
       baseline=1.68, distribution=D)
def set_feedstock_FEC_CF(CF):
    CFs['FEC_CFs']['FGHTP Corn stover'] = CF
    
# =============================================================================
# Pretreatment parameters
# =============================================================================

M202 = find.unit.M202
D = shape.Triangle(0.25, 0.3, 0.4)
@param(name='Pretreatment solid loading', element=M202, kind='coupled', units='%', 
       baseline=0.3, distribution=D)
def set_pretreatment_solid_loading(loading): 
    M202.solid_loading = loading


# baseline imass discrepancy
pretreatment_sulfuric_acid = find.stream.pretreatment_sulfuric_acid
D = shape.Triangle(10, 22.1, 35)
@param(name='Pretreatment sulfuric acid loading', element=pretreatment_sulfuric_acid,
       kind='coupled', units='mg/g-dry feedstock', baseline=22.1, distribution=D)
def set_pretreatment_sulfuric_acid_loading(loading): 
    feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    pretreatment_sulfuric_acid.imass['H2SO4'] = feedstock_dry_mass*loading/1000*0.93
    pretreatment_sulfuric_acid.imass['H2O'] = feedstock_dry_mass*loading/1000*0.07

R201 = find.unit.R201
D = shape.Triangle(0.06, 0.099, 0.12)
@param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='%',
       baseline=0.099, distribution=D)
def set_R201_glucan_conversion(X): R201.pretreatment_rxns[0].X = X    

D = shape.Triangle(0.8, 0.9, 0.92)
@param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R201_xylan_conversion(X): R201.pretreatment_rxns[4].X = X        


# =============================================================================
# Conversion parameters
# =============================================================================

M301 = find.unit.M301
R302 = find.unit.R302
R303 = find.unit.R303

D = shape.Triangle(0.175, 0.2, 0.25)
@param(name='Enzymatic hydrolysis solid loading', element=M301, kind='coupled', units='%',
       baseline=0.2, distribution=D)
def set_R301_hydrolysis_solid_loading(loading): M301.solid_loading = loading

D = shape.Triangle(10, 20, 30)
@param(name='Enzyme loading', element=M301, kind='coupled', units='mg/g glucan',
       baseline=20, distribution=D)
def set_R301_enzyme_loading(loading): M301.enzyme_loading = loading

# Enzymatic hydrolysis
D = shape.Triangle(0, 24, 56)
@param(name='Enzymatic hydrolysis time', element=R301, kind='coupled', units='hr',
       baseline=24, distribution=D)
def set_R301_hydrolysis_time(tau): R301.tau_saccharification = tau

D = shape.Triangle(0.75, 0.9, 0.948-1e-6)
@param(name='Enzymatic hydrolysis glucan-to-glucose', element=R301, kind='coupled', units='%',
       baseline=0.9, distribution=D)
def set_R301_glucan_conversion(X): R301.saccharification_rxns[2].X = X

# Fermentation
# D = shape.Triangle(76, 120, 145)
# @param(name='Fermentation time', element=R302, kind='coupled', units='hr',
#        baseline=120, distribution=D)
# def set_R302_fermentation_time(tau): R302.tau_cofermentation = tau


D = shape.Triangle(0.684, 0.76, 0.8360000000000001)
@param(name='Productivity', element=R302, kind='coupled', units='g/L/hr',
       baseline=0.76, distribution=D)
def set_HP_productivity(productivity):
    # R301.productivity = productivity
    # R302.productivity = productivity * R302.ferm_ratio
    # spec.load_productivity(productivity)
    # spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=productivity)
    
    spec.spec_3 = productivity
    
D = shape.Triangle(5, 10, 15)
@param(name='CSL loading', element=R301, kind='coupled', units='g/L',
       baseline=10, distribution=D)
def set_CSL_loading(loading): R302.CSL_loading = loading



########################### COMMENT Y,T OUT FOR TRY P_MC ANALYSIS ############################

D = shape.Triangle(0.49*0.8, 0.49, 0.49*1.2) # +/- 20% of baseline
@param(name='3-Hydroxypropionic acid yield', element=R302, kind='coupled', units='% theoretical',
        baseline=0.49, distribution=D)
def set_R302_HP_yield(X):
    # R302_X = R302.cofermentation_rxns.X
    # R302_X[0] = R302_X[2] = X
    # R303_X = R303.cofermentation_rxns.X
    # R303_X[0] = R303_X[3] = X * R303.ferm_ratio
    # spec.load_specifications(spec_1=X, spec_2=spec.spec_2, spec_3=spec.spec_3)
    # spec.load_yield(X)
    spec.spec_1 = X

D = shape.Triangle(54.8*0.8, 54.8, 54.8*1.2) # +/- 20% of baseline
@param(name='3-Hydroxypropionic acid titer', element=R302, kind='coupled', units='g/L',
        baseline=54.8, distribution=D)
def set_R302_HP_titer(X):
    # R302_X = R302.cofermentation_rxns.X
    # R302_X[0] = R302_X[2] = X
    # R303_X = R303.cofermentation_rxns.X
    # R303_X[0] = R303_X[3] = X * R303.ferm_ratio
    # spec.load_specifications(spec_1=spec.spec_1, spec_2=X, spec_3=spec.spec_3)
    # spec.load_titer(X)
    spec.spec_2 = X
    
##############################################################################################

D = shape.Triangle(0.032, 0.040, 0.048)
@param(name='Acetic acid and glycerol yield', element=R303, kind='coupled', units='% theoretical',
        baseline=0.040, distribution=D)
def set_R301_acetic_acid_yield(X): 
    # 1e6 is to avoid generating tiny negative flow (e.g., 1e-14) in R301
    # R302_X = R302.cofermentation_rxns.X
    ferm_ratio = R303.ferm_ratio
    
    X1 = min(X, 1-1e-6-R302.glucose_to_HP_rxn.X-R302.glucose_to_biomass_rxn.X)
    X2 = min(X, 1-1e-6-R302.xylose_to_HP_rxn.X-R302.xylose_to_biomass_rxn.X)
    
    R302.glucose_to_acetic_acid_rxn.X = X1
    R303.glucose_to_acetic_acid_rxn.X = X1 * ferm_ratio
    
    R302.xylose_to_acetic_acid_rxn.X = X2
    R303.xylose_to_acetic_acid_rxn.X = X2 * ferm_ratio
    
    X1_glycerol = X1 if X1==X else 0.
    X2_glycerol = X2 if X2==X else 0.  
    
    R302.glucose_to_glycerol_rxn.X = X1_glycerol
    R303.glucose_to_glycerol_rxn.X = X1_glycerol * ferm_ratio
    
    R302.xylose_to_glycerol_rxn.X = X2_glycerol
    R303.xylose_to_glycerol_rxn.X = X2_glycerol * ferm_ratio
    
    
    
    
    
    
    
    
    
    
    # R302_X[1] = R302_X[3] = X
    # R303_X = R303.cofermentation_rxns.X
    # X = min(X*R303.ferm_ratio, 1-1e-6-R303_X[0]-R303_X[2])
    # R303_X[1] = R303_X[3] = X




# D = shape.Uniform(0.4, 0.6)
# @param(name='Unfermented sugars routed to CO2 generation (balance is routed to cell mass generation)', element=R303, kind='coupled', units='% theoretical',
#        baseline=0.5, distribution=D)
# def set_sugars_to_CO2_gen(X): 
#    R302.CO2_generation_rxns.X[0] = R302.CO2_generation_rxns.X[1] = X
#    R303.CO2_generation_rxns.X[0] = R303.CO2_generation_rxns.X[1] = X * R303.ferm_ratio
   
   
   
# D = shape.Triangle(0.05, 0.07, 0.1)
# @param(name='Innoculum ratio', element=R301, kind='coupled', units='%',
#        baseline=0.07, distribution=D)
# def set_innoculum_ratio(ratio): R301.inoculum_ratio = ratio



# Seed train fermentation yield as a ratio of the main fermenter
# D = baseline_uniform(36, 0.1)
# @param(name='Seed train time', element=R303, kind='coupled', units='hr',
#        baseline=36, distribution=D)
# def set_R303_fermentation_time(tau): R303.tau_batch = tau





# D = shape.Triangle(0.8, 0.9, 1)
# @param(name='Seed train yield', element=R302, kind='coupled', units='% of R301',
#        baseline=0.9, distribution=D)
# def set_R303_ratio(ratio):
#     R302_X = R302.cofermentation_rxns.X
#     R303_X = R303.cofermentation_rxns.X
#     ratio = min(ratio, (1-1e-6-R303_X[2])/(R302_X[0]+R302_X[1]))
#     R303.ferm_ratio = ratio

# =============================================================================
# Separation parameters
# =============================================================================

S402 = find.unit.S402
D = shape.Triangle(0.95, 0.995, 1.)
@param(name='Gypsum split', element=S402, kind='coupled', units='',
       baseline=0.995, distribution=D)
def set_S402_gypsum_split(split):
    gypsum_index = S402.chemicals.index('Gypsum')
    S402.split[gypsum_index] = split

R401 = find.unit.R401
D = baseline_triangle(1., 0.1)
@param(name='Acidulation time', element=R401, kind='coupled', units='hr',
       baseline=1., distribution=D)
def set_R401_tau(tau):
    R401.tau = tau

R402 = find.unit.R402
# D = baseline_triangle(0.95, 0.05)
D = shape.Triangle(0.64, 0.8, 0.96)
@param(name='Dehydration conversion', element=R402, kind='coupled', units='',
       baseline=0.8, distribution=D)
def set_R402_conversion(X):
    R402.dehydration_reactions[0].X = X
    
# R403 = find.unit.R403
# D = shape.Triangle(0.72, 0.8, 0.88)
# @param(name='Hydrolysis conversion', element=R403, kind='coupled', units='%',
#        baseline=0.8, distribution=D)
# def set_R403_conversion_factor(X):
#     R403.hydrolysis_rxns.X[:] = X
    

# =============================================================================
# Facilities parameters
# =============================================================================

# Facilities are currently not in unit path. thus not set the element here
D = baseline_uniform(0.8, 0.1)
@param(name='BT combustion efficiency', element=BT, kind='coupled', units='%',
       baseline=0.8, distribution=D)
def set_BT_combustion_efficiency(efficiency):
    BT.boiler_efficiency = efficiency

# All parameters
parameters = HP_model.get_parameters()

index_TEA = len(metrics)


# %%

# =============================================================================
# Model to evalute system across internal rate of return
# =============================================================================

def create_IRR_metrics(IRR):
    def get_IRR_based_MSP():
        HP_tea.IRR = IRR
        return get_MSP()
    return [Metric('Minimum selling price', get_IRR_based_MSP, '$/kg', f'IRR={IRR:.0%}'),
            Metric('Net present value', get_NPV, '$', f'IRR={IRR:.0%}')]

IRRs = np.linspace(0, 0.4, 41)
IRR_metrics = sum([create_IRR_metrics(IRR) for IRR in IRRs],[])

HP_model_IRR = Model(HP_sys, IRR_metrics)
HP_model_IRR.set_parameters(parameters)


# %% 

# =============================================================================
# Model to evalute system across 3-Hydroxypropionic acid yield
# =============================================================================

HP_model_LA_yield = Model(HP_sys, metrics)

def set_LA_yield(LA_yield):
    R301_X = R301.cofermentation_rxns.X
    R301_X[0] = R301_X[3] = LA_yield
    R301_X[1] = R301_X[4] = min(1-1e-6-R301_X[0]-R301_X[2], R301_X[1])
    R302_X = R302.cofermentation_rxns.X
    R302_X[0] = R302_X[3] = R301_X[0] * R302.ferm_ratio
    R302_X[1] = R302_X[4] = R301_X[1] * R302.ferm_ratio

LA_yield_parameters = tuple([i for i in parameters if not i.name=='3-Hydroxypropionic acid yield'])
HP_model_LA_yield.set_parameters(LA_yield_parameters)


# %%

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
prices = np.concatenate((prices, np.array([71.26])))
feedstock_price_metrics = sum([create_feedstock_price_metris(price) 
                                for price in prices],[])

def set_feedstock_carbs(carbs_content):
    carbs = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')
    dry_mass = feedstock.F_mass.copy() - feedstock.imass[('H2O',)].copy()
    old_carbs_mass_total = feedstock.imass[carbs].sum().copy()
    ratio = feedstock.get_normalized_mass(carbs)
    new_carbs_mass = dry_mass * carbs_content * ratio
    feedstock.set_flow(new_carbs_mass, 'kg/hr', carbs)
    mass_diff = new_carbs_mass.sum() - old_carbs_mass_total
    feedstock.imass['Extract'] -= mass_diff
    if any(feedstock.mass < 0):
        raise ValueError(f'Carbohydrate content of {carbs_content*100:.0f}% dry weight is infeasible')

HP_model_feedstock = Model(HP_sys, feedstock_price_metrics)

param = HP_model_feedstock.parameter

# Set a fake parameter to enable evaluation across internal rate of return
feedstock = find.stream.feedstock
D = shape.Uniform(0.9, 1.1)
@param(name='Fake parameter', element=feedstock, kind='coupled', units='',
       baseline=1, distribution=D)
def set_fake_parameter(anything): pass

index_IRR = len(metrics)
# %%

# =============================================================================
# Model to evalute system across feedstock succinic acid content
# =============================================================================

def set_feedstock_succinic_acid_content(SA_content):
    dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    feedstock.imass['SuccinicAcid'] = SA_content * dry_mass
    # Use Extract to close mass balance
    feedstock.imass['Extract'] -= (feedstock.F_mass-feedstock.imass['H2O']) - dry_mass
    if any(feedstock.mass<0):
        raise ValueError(f'Succinic acid content of {SA_content*100:.0f}% dry weight is infeasible')

HP_model_SA_content = Model(HP_sys, metrics)
HP_model_SA_content.set_parameters(parameters)


# %% 

# =============================================================================
# Model to evalute system across HXN minimum approach temperature
# =============================================================================
get_HXN_util_cost = lambda: HXN.utility_cost
get_HXN_util_savings = lambda: - 100*(HXN.utility_cost*350*24) / (HP_sys.utility_cost - HXN.utility_cost*350*24)
def create_HXN_T_min_app_metrics(T_min_app):
    return [Metric('HXN utility savings', get_HXN_util_savings, '%', f'T_min_app={T_min_app:.0f} [K]')]

T_min_apps = np.linspace(1, 20, 21)
HXN_T_min_app_metrics = sum([create_HXN_T_min_app_metrics(T_min_app) for T_min_app in T_min_apps],[])


HP_model_HXN_T_min_app = Model(HP_sys, metrics)

def set_HXN_T_min_app(T_min_app):
    HXN.T_min_app = T_min_app

HP_model_HXN_T_min_app.set_parameters(parameters)



# %% Evaluate
# N_samples = 100
# rule = 'L' # For Latin-Hypercube sampling
# model=HP_model
# samples = model.sample(N_samples, rule)
# model.load_samples(samples)
# model.evaluate()
# model.table # All evaluations are stored as a pandas DataFrame



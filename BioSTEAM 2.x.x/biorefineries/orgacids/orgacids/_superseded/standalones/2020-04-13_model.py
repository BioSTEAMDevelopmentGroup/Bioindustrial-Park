#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:24:42 2020

@author: yalinli_cabbi
"""


# %% Setup

from biosteam.evaluation import evaluation_tools as tools
from biosteam.evaluation import Model, Metric
from orgacids.system import lactic_acid, R301, orgacids_tea, orgacids_sys_no_BT_tea, \
    sub_systems, sub_teas
            

# %% Define metric functions

# Minimum selling price of lactic acid ($/kg)
get_MSP = lambda: orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_BT_tea)
# Purity (%) of lactic acid in the final product
get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
# Recovery (%) = lactic acid recovered/total produced lactic acid
get_recovery = lambda: lactic_acid.imass['LacticAcid']/R301.outs[1].imass['LacticAcid']

metrics = [Metric('Minimum selling price', get_MSP, '$/kg'),
           Metric('Product purity', get_purity, '%'),
           Metric('Separation recovery', get_recovery, '%')]

# Get utility usage including heating (kJ/hr), cooling (kJ/hr), and power (kW)
def get_utility_demands(systems):
    heating_demands = cooling_demands = power_demands = {}
    for system in systems:
        if system.ID == 'BT_sys':
            heating_demands[system.ID] = -system.path[0].heat_generated
            cooling_demands[system.ID] = system.path[0].cooling_need
            power_demands['BT_sys_equipment'] = system.path[0].equipment_electricity
            power_demands['BT_sys_steam'] = -system.path[0].generated_electricity
        else:                
            system_heating_demand = system_cooling_demand = system_power_demand = 0
            for unit in system.units:
                system_heating_demand += sum(unit.heat_utilities[i].duty 
                                          for i in range(len(unit.heat_utilities))
                                          if unit.heat_utilities[i].ID=='low_pressure_steam')
                system_cooling_demand += sum(unit.heat_utilities[i].duty 
                                          for i in range(len(unit.heat_utilities))
                                          if unit.heat_utilities[i].ID=='cooling_water')
                if unit.power_utility:
                    system_power_demand += unit.power_utility.rate

            heating_demands[system.ID] = system_heating_demand
            cooling_demands[system.ID] = system_cooling_demand
            power_demands[system.ID] = system_power_demand
    
    heating_demands['overall'] = sum(heating_demands[i] for i in heating_demands)
    cooling_demands['overall'] = sum(cooling_demands[i] for i in cooling_demands)
    power_demands['overall'] = sum(power_demands[i] for i in power_demands)
    return (heating_demands, cooling_demands, power_demands)
heating_demands, cooling_demands, power_demands = get_utility_demands(sub_systems)


# Get system costs
def get_costs(teas):
    TCI = AOC = Depreciation = {}
    for tea in teas:
        TCI[tea.system.ID] = tea.TCI
        AOC[tea.system.ID] = tea.AOC
        Depreciation[tea.system.ID] = tea.annual_depreciation
    return (TCI, AOC, Depreciation)
TCI, AOC, Depreciation = get_costs(sub_teas)


for system in sub_systems:
    metrics.extend(
        Metric('Heating demand', get_utility_demands())
        )


# %% 
'''
sub_systems = (feedstock_sys,
               pretreatment_sys,
               fermentation_sys,
               separation_sys,
               wastewater_sys,
               facilities_sys_no_utilities,
               BT_sys, CT_sys, PWC_sys)



sub_teas = (feedstock_sys_tea,
            pretreatment_sys_tea,
            fermentation_sys_tea,
            separation_sys_tea,
            wastewater_sys_tea,
            BT_sys_tea,
            facilities_sys_no_BT_tea)

'''











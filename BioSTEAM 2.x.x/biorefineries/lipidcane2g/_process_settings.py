# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 10:02:05 2019

@author: yoelr
"""
import biosteam as bst

__all__ = ('load_process_settings',)

# %% Process settings

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 607.5 # 2019
    bst.PowerUtility.price = 0.065
    HeatUtility = bst.HeatUtility
    steam_utility = HeatUtility.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    HeatUtility.get_agent('cooling_water').regeneration_price = 0
    HeatUtility.get_agent('chilled_water').heat_transfer_price = 0

ton2kg = 1./907.18474 # ton/hr to kg/hr
ethanol_density_kgL = 0.789 # kg/L
liter_per_gallon = 3.78541
ethanol_cost = 2.15 # USD/gal
ethanol_density_kggal = liter_per_gallon * ethanol_density_kgL # kg/gal
enzyme_price = 4.24 * 50/1000 # USD/kg

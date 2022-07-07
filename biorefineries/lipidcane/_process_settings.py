# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 10:02:05 2019

@author: yoelr
"""
import biosteam as bst

__all__ = ('load_process_settings', 'price',)

# %% Process settings

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 567 # 2013
    bst.PowerUtility.price = 0.065
    HeatUtility = bst.HeatUtility
    steam_utility = HeatUtility.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    HeatUtility.get_agent('cooling_water').regeneration_price = 0
    HeatUtility.get_agent('chilled_water').heat_transfer_price = 0

# Raw material price (USD/kg)
price = {'Lipid cane': 0.03455, # 70% m.c
         'Methanol': 0.547,
         'Water': 0.000353,
         'HCl': 0.205,
         'Lime': 0.077,
         'H3PO4': 0, # 0.700, # TODO: find price
         'NaOCH3': 2.93,
         'NaOH':0.41,
         'Protease': 0.5,
         'Polymer': 0, # TODO: find price
         'Crude glycerol': 0.21,
         'Biodiesel': 1.38, # assuming density = 870 kg/m3
         'Ethanol': 0.789,
         'Waste': -0.33,
         'Gasoline': 0.756,
         'Crude oil': 0.661} # 30 cts / lb (vegetable); 5yr average https://www.ams.usda.gov/mnreports/lswagenergy.pdf

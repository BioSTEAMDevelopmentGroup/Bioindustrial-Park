# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:15 2019

@author: yoelr
"""
import biosteam as bst

bst.CE = 525


factor = 1/907.18474 # ton/hr to kg/hr
ethanol_density_kgL = 0.789 # kg/L
liter_per_gallon = 3.78541
ethanol_cost = 2.15 # USD/gal
ethanol_density_kggal = liter_per_gallon * ethanol_density_kgL # kg/gal
enzyme_price = 0.66 # USD/kg
biogas_price = 1.03 #It includes the cost to purify the biogas

price = {'Ethanol':   0.8031,
         'Feedstock': 0.06, #checked in Elizabeth's article
         'Sulfuric acid': 81.39 * factor,
         'Ammonia': 406.96 * factor,
         'Caustic': 135.65 * factor,
         'Boiler chems': 4532.17 * factor,
         'FGD lime': 180.87 * factor,
         'Cooling tower chems': 2716.1 * factor,
         'Makeup water': 0.1e-3, # checked with Ola
         'Ash disposal': -28.86 * factor,
         'Electricity': 0.0572, # USD/kWh
         'Denaturant': 0.756,
         'Enzyme': enzyme_price, #checked in Elizabeth's article 
         'Pure biogas': biogas_price} #checked in Elizabeth's article
bst.PowerUtility.price = price['Electricity']
_ha = bst.HeatUtility.get_heating_agent('low_pressure_steam')
_ha.heat_transfer_efficiency = 0.85
_ha.T = 529.2
_ha.P = 44e5
_ha.regeneration_price = 0.30626
_CW = bst.HeatUtility.get_cooling_agent('cooling_water')
_CW.T = 28 + 273.15
_CW.T_limit = _CW.T + 9
_CW.regeneration_price = 0
bst.HeatUtility.get_cooling_agent('chilled_water').heat_transfer_price = 0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:30:10 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

@author: yalinli_cabbi
"""

import biosteam as bst

#TODO: need to consider labor and chemical costs across different years
bst.CE = 541.7 # Year 2016

kg_per_ton = 907.18474
liter_per_gallon = 3.78541

feedstock_cost = 0.0546 # USD/kg
enzyme_cost = 0.507 # USD/kg
# Placeholder for the assumed chemicals used in separation system
methanol_cost = 1 #!!! USD/kg need to find a methanol cost

HPA_price = 50 # USD/kg, 3-hydroxypropionic acid
# Not including adipic acid in fermentation modeling now
AA_price = 1.639 # USD/kg, adipic acid
BA_price = 2 # USD/kg, butyric acid
CA_price = 0.7 # USD/, citric acid
LA_price = 1.1 # USD/kg, lactic acid
# Not including muconic acid in fermentation modeling now
MA_price = AA_price # USD/kg, cis, cis-muconic acid; no data availalbe
PA_price = 1.5 # USD/kg, propionic acid
SA_price = 2 # USD/kg

# Prices need double-checking
price = {'3-hydroxypropionic acid': HPA_price,
         'Adipic acid': AA_price,
         'Butyric acid': BA_price,
         'Citric acid': CA_price,
         'Lactic acid': LA_price,
         'cis, cis-Muconic acid': MA_price,
         'Propionic aicd': PA_price,
         'Succinic acid': SA_price,
         'Feedstock': feedstock_cost, # changed from Humbird et al., 2011 to Davis et al., 2018
         'Sulfuric acid': 81.39 / kg_per_ton,
         'Ammonia': 406.96 / kg_per_ton,
         'CSL': 51.55 / kg_per_ton,
         'DAP': 895.32 / kg_per_ton,
         'Sorbitol': 1021.93 / kg_per_ton,
         'Glucose': 526.52 / kg_per_ton,
         'Caustic': 135.65 / kg_per_ton,
         'Boiler chems': 4532.17 / kg_per_ton,
         'Lime': 180.87 / kg_per_ton,
         'Cooling tower chems': 2716.1 / kg_per_ton,
         'Makeup water': 0.23 / kg_per_ton,
         # Price of ash is negative because it's a product stream
         'Ash disposal': -28.86 / kg_per_ton,
         'Electricity': 0.0572, # USD/kWh
         'Denaturant': 0.756,
         'Enzyme': enzyme_cost,
         'Methanol': methanol_cost} 
bst.PowerUtility.price = price['Electricity']

_ha = bst.HeatUtility.heating_agents['Low pressure steam']
_ha.efficiency = 0.85
_ha.T = 529.2
_ha.P = 44e5
_ha.price_kmol = 0.30626
_ha.Hvap = 30235.86
_CW = bst.HeatUtility.cooling_agents['Cooling water']
_CW.T = 28 + 273.15
_CW.T_limit = _CW.T + 9
_CW.price_kmol = 0
bst.HeatUtility.cooling_agents['Chilled water'].price_kJ = 0

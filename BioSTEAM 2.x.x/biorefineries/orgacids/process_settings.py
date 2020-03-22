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
# Should convert those costs using respective indices, not GDP indices
bst.CE = 541.7 # Year 2016

_kg_per_ton = 907.18474
_liter_per_gallon = 3.78541
_GDP_2007to2016 = 1.114 / 0.961 #!!! Should be changed to chemical index

# USD/kg, 20% moisture content
feedstock_cost = 58.5 /_kg_per_ton*0.8 

# $4.24/kg protein in P51 of Humbird et al., 2007$
enzyme_cost = 4.24 # USD/kg

# 2.21 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 
# based on Annual Energy Outlook from Energy Information Adiministration
# (https://www.eia.gov/outlooks/aeo/)
# 0.789 is density of ethanol in kg/L
ethanol_price = 2.21 / (_liter_per_gallon*0.789)

# $2.10/gal from Humbird et al. in 2007$, density of gasoline is 2.819 kg/gal
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels
# from H2 Tools maintained by Pacific Northwest National Laboratory
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)
denaturant_price = 2.10 / 2.819 * _GDP_2007to2016

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_cost = 466833/5/(24*365*0.96)

HPA_price = 50 # USD/kg, 3-hydroxypropionic acid
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
         'Sulfuric acid': 81.39 / _kg_per_ton * _GDP_2007to2016,
         'Ammonia': 406.96 / _kg_per_ton * _GDP_2007to2016,
         'CSL': 51.55 / _kg_per_ton * _GDP_2007to2016,
         'DAP': 895.32 / _kg_per_ton * _GDP_2007to2016,
         'Sorbitol': 1021.93 / _kg_per_ton * _GDP_2007to2016,
         'Glucose': 526.52 / _kg_per_ton * _GDP_2007to2016,
         'Caustic': 135.65 / _kg_per_ton * _GDP_2007to2016,
         'Boiler chemicals': 4532.17 / _kg_per_ton * _GDP_2007to2016,
         'Lime': 180.87 / _kg_per_ton * _GDP_2007to2016,
         'Cooling tower chems': 2716.1 / _kg_per_ton * _GDP_2007to2016,
         'Makeup water': 0.23 / _kg_per_ton * _GDP_2007to2016,
         # Cost of ash is negative because it's a product stream
         'Ash disposal': -28.86 / _kg_per_ton * _GDP_2007to2016,
         'Electricity': 0.0572 * _GDP_2007to2016, # USD/kWh
         'Enzyme': enzyme_cost * _GDP_2007to2016,
         'Ethanol': ethanol_price,
         'Denaturant': denaturant_price,
         'Baghouse bag': baghouse_bag_cost
         } 
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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:30:10 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

@author: yalinli_cabbi
"""

import biosteam as bst

bst.CE = 541.7 # Year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_chemical_2011to2016 = 102.5 / 91.7
# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
# changed from Humbird et al., 2011 to Davis et al., 2018
feedstock_cost = 71.26 /_kg_per_ton * 0.8 

# 2.18 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 
# based on Annual Energy Outlook from Energy Information Adiministration
# (https://www.eia.gov/outlooks/aeo/)
# 0.789 is density of ethanol in kg/L
ethanol_price = 2.18 / (_liter_per_gallon*0.789)

# 2.86 is the average motor gasoline price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook from Energy Information Adiministration
# (https://www.eia.gov/outlooks/aeo/)
# density of gasoline is 2.819 kg/gal
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels
# from H2 Tools maintained by Pacific Northwest National Laboratory
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)
denaturant_price = 2.86 / 2.819

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of Davis et al., 2018 (BDO scenario)
# 7880 is operating hours/yr on Page 10 of Davis et al., 2018,
# cost is negative because it's a product stream
ash_disposal_cost = -1.41e6 / (4279*7880)

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

# All in 2016$/kg
price = {'3-hydroxypropionic acid': HPA_price,
         'Adipic acid': AA_price,
         'Butyric acid': BA_price,
         'Citric acid': CA_price,
         'Lactic acid': LA_price,
         'cis, cis-Muconic acid': MA_price,
         'Propionic aicd': PA_price,
         'Succinic acid': SA_price,
         'Feedstock': feedstock_cost, 
         'Sulfuric acid': 0.0430 * _lb_per_kg,
         'Ammonia': 0.1900 * _lb_per_kg,
         'CSL': 0.0339 * _lb_per_kg,
         'DAP': 0.1645 * _lb_per_kg,
         'Caustic': 0.0832 * _lb_per_kg * _chemical_2011to2016, # Davis et al. 2013
         'Boiler chemicals': 2.9772 * _lb_per_kg,
         'Lime': 0.1189 * _lb_per_kg,
         'Cooling tower chems': 1.7842 * _lb_per_kg,
         'Makeup water': 0.0002 * _lb_per_kg,
         # Cost of ash is negative because it's a product stream
         'Ash disposal': ash_disposal_cost,
         'Electricity': 0.068, # USD/kWh from Davis et al., 2018
         # $6.16/kg protein in 2016$, P25 of Davis et al., 2018$
         'Enzyme': 6.16,
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

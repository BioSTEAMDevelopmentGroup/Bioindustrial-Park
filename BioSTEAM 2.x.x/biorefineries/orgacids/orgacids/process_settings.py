#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:30:10 2019

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

@author: yalinli_cabbi
"""

import biosteam as bst
from orgacids.chemicals import chems

bst.CE = 541.7 # Year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_chemical_2011to2016 = 102.5 / 91.7
_ft3_per_m3 = 35.3147

# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
# changed from Humbird et al., 2011 to Davis et al., 2018
feedstock_cost = 71.26 / _kg_per_ton * 0.8 

# Methanol price from Goellner et al., production from natural gas (Case 3),
# average of two load structures in 2011$,
# crude methanol with ~1% CO2 and 1% H2O
methanol_cost = (311.17+345.39)/2/_kg_per_ton * _chemical_2011to2016

# 2.18 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 
# based on Annual Energy Outlook from Energy Information Adiministration
# (https://www.eia.gov/outlooks/aeo/)
# 0.789 is density of ethanol in kg/L
ethanol_cost = 2.18 / (_liter_per_gallon*0.789)

# 2.86 is the average motor gasoline price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook from Energy Information Adiministration
# (https://www.eia.gov/outlooks/aeo/)
# density of gasoline is 2.819 kg/gal
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels
# from H2 Tools maintained by Pacific Northwest National Laboratory
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)
denaturant_cost = 2.86 / 2.819

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of Davis et al., 2018 (BDO scenario)
# 7880 is operating hours/yr on Page 10 of Davis et al., 2018,
# cost is negative because it's a product stream
ash_disposal_cost = -1.41e6 / (4279*7880)

# Assums no cost/credit for baseline, the same as ash disposal for lower end,
# for the higher end (i.e., positive selling price indicating profit), use 
# USGS 2015-2019 average  free on bard price in $/metric ton for crude gypsum. 
# National Minerals Information Center. Mineral Commodity Summaries 2020; 
# U.S. Geological Survey, 2020.
# Assuming all prices were in their nominal year (e.g., 2015 price in 2015$)
# and adjusted to 2016$
# 2015: 7.80 * 1.114 / 1.100 = 7.90
# 2016: 8.00
# 2017: 7.50 * 1.114 / 1.134 = 7.37
# 2018: 8.30 * 1.114 / 1.157 = 7.99
# 2019: 8.00 * 1.114 / 1.185 = 7.52
# (7.90+8.00+7.37+7.99+7.52) / 5 = 7.76 (in metric tonne)
# For the lower end (i.e., negative selling price indicating cost), use price from
# Aden et al., 2002: $0.0094/lb in 2000$ = 0.0094*1.114/0.802*2.20462 = $0.0288/kg
# in 2016$
gypsum_cost = 0

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_cost = 466833 / 5 / (24*365*0.96)

# 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# Using the same conversion, lower and upper bounds should be 
# $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chems.CH4.MW
natural_gas_price = 4.70/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW)

# All in 2016$/kg
price = {# IHS markit report, 2016 price, US market, 
         # technical grade, average of lower and higher range
         'Lactic acid': 1.88, 
         'Feedstock': feedstock_cost, 
         'Sulfuric acid': 0.0430 * _lb_per_kg,
         # 0.1900 is for NH3
         'AmmoniumHydroxide': 0.1900 * _lb_per_kg * 17.031/35.046,
         'CSL': 0.0339 * _lb_per_kg,
         'DAP': 0.1645 * _lb_per_kg,
         'Caustic': 0.0832 * _lb_per_kg * _chemical_2011to2016, # Davis et al. 2013
         'Boiler chemicals': 2.9772 * _lb_per_kg,
         'Lime': 0.1189 * _lb_per_kg,
         # Cooling tower chemis in Humbird et al. (CIP is a part of cooling tower)
         'CIP chems': 1.7842 * _lb_per_kg,
         'Makeup water': 0.0002 * _lb_per_kg,
         # Cost of ash is negative because it's a product stream
         'Ash disposal': ash_disposal_cost,
         'Gypsum': gypsum_cost,
         'Electricity': 0.068, # USD/kWh from Davis et al., 2018
         # $6.16/kg protein in 2016$, P25 of Davis et al., 2018
         'Enzyme': 6.16,
         'Methanol': methanol_cost,
         'Ethanol': ethanol_cost,
         'Denaturant': denaturant_cost,
         'Baghouse bag': baghouse_bag_cost,
         'Natural gas': natural_gas_price
         } 
bst.PowerUtility.price = price['Electricity']

_mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
_hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
# Adjusted to Humbird LP/HP steam temperature
# biosteam native mps T is lower than LP steam T in Humbird, thus adjust mps.T
_mps.T = 233 + 273.15
_hps.T = 266 + 273.15
_cw = bst.HeatUtility.get_cooling_agent('cooling_water')
_cw.T = 28 + 273.15
_cw.T_limit = _cw.T + 9
_cw.regeneration_price = 0
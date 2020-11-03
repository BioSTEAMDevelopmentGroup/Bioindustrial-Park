#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Sat Jun 27 13:51:49 2020

References:
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
@author: yalinli_cabbi
"""


# %%

import biosteam as bst
from ethanol_adipic.chemicals import chems

bst.CE = 541.7 # year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147
_J_per_BTU = 1055.06
_GDP_2007to2016 = 1.160

# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
_feedstock_factor = _kg_per_ton / 0.8
feedstock_price = 71.3 / _feedstock_factor

# 2.86 is the average motor gasoline price between 2010-2019 in 2016 $/gal	
# based on AEO from EIA, density of gasoline is 2.819 kg/gal	
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels	
# from H2 Tools maintained by Pacific Northwest National Laboratory	
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)	
denaturant_price = 2.86 / 2.819

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of ref [2] (BDO scenario)
# 7880 is operating hours/yr on Page 10 of ref [2],
# cost is negative because it's a product stream
ash_disposal_price = -1.41e6 / (4279*7880)

# Baseline from ref [2], lower bound is 2015-2019 average of 	
# hydrate lime in $/ton at plant from Mineral Commodity Summaries 2020.	
# 2015: 146.40 * (1.114/1.100) / 907.18474 = 0.163	
# 2016: 145.50 / 907.18474 = 0.160	
# 2017: 147.10 * (1.114/1.134) / 907.18474 = 0.159	
# 2018: 151.50 * (1.114/1.157) / 907.18474 = 0.161	
# 2019: 151.00 * (1.114/1.185) / 907.18474 = 0.156	
# (0.163+0.160+0.159+0.161+0.156) / 5 = 0.160	
# Upper bound is +10% from baseline = 0.1189 * _lb_per_kg * 1.1 = 0.288
lime_price = 0.1189 * _lb_per_kg

# The original cost is $466,183 every 5 years in ref [1], converted to per hour using
# the assumed 96% uptime
baghouse_bag_price = 466183/5/(24*365*0.96) * _GDP_2007to2016

# $5/MM BTU
CH4_LHV = chems.CH4.LHV
CH4_MW = chems.CH4.MW
CH4_cost_per_J = 5/(1e6*_J_per_BTU)
CH4_cost_per_mol = CH4_cost_per_J * -CH4_LHV
natural_gas_price = CH4_cost_per_mol * (1000/CH4_MW)


# All in 2016$/kg
price = {'Feedstock': feedstock_price, 	
         'H2SO4': 0.0430 * _lb_per_kg,	
         # 0.1900 is for NH3	
         'NH4OH': 0.1900 * _lb_per_kg * 17.031/35.046,
         'NaOH': 0.2384 * _lb_per_kg,
         'CSL': 0.0339 * _lb_per_kg,	
         'DAP': 0.1645 * _lb_per_kg,
         # $6.16/kg protein in 2016$, P25 of ref [2]
         'Enzyme': 6.16,
         'H2': 0.7306 * _lb_per_kg,
         'Hydrogenation catalyst': 528 * _lb_per_kg,
         'WWT polymer': 2.6282 * _lb_per_kg,
         'Natural gas': natural_gas_price,
         'Lime': lime_price,         
         'Boiler chems': 2.9772 * _lb_per_kg,	
         'Baghouse bag': baghouse_bag_price,
         'Cooling tower chems': 1.7842 * _lb_per_kg,	
         'Makeup water': 0.0002 * _lb_per_kg,
         # Cost of ash is negative because it's a product stream	
         'Ash disposal': ash_disposal_price,
         'Electricity': 0.068,
         'Denaturant': denaturant_price,
         'Ethanol': 0.3370 * _lb_per_kg,
         'Adipic acid': 0.8554 * _lb_per_kg,
         'Sodium sulfate': 0.0706 * _lb_per_kg
         }
    
bst.PowerUtility.price = price['Electricity']

_lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
_mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
_hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')

# Adjusted to LP/HP steam temperature as in ref [1]
# biosteam native mps T is lower than LP steam T in ref [1], thus adjust mps.T
_mps.T = 233 + 273.15
_hps.T = 266 + 273.15
_cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
_chilled = bst.HeatUtility.get_cooling_agent('chilled_water')
_cooling.T = 28 + 273.15
_cooling.T_limit = _cooling.T + 9

# Side steam in CHP not a heat utility, thus will cause problem in TEA utility
# cost calculation if price not set to 0 here, costs for regeneration of heating
# and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
for i in (_lps, _mps, _hps, _cooling, _chilled):
    i.heat_transfer_price = i.regeneration_price = 0
    
    

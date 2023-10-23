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

@author: sarangbhagwat
"""

import biosteam as bst
import thermosteam as tmo
from biorefineries.TAL.chemicals_data import chems

bst.CE = 541.7 # year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147

_GDP_2007_to_2016 = 1.114 / 0.961
_GDP_2008_to_2016 = 1.114 / 0.990
_chemical_2011to2016 = 102.5 / 91.7
_chemical_2013to2016 = 102.5 / 101.3
_chemical_2014to2016 = 102.5 / 105.3

_chemical_2017to2016 = 102.5 / 106.9
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb
_chemical_2022to2016 = 102.5 / 145.3


# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
# changed from Humbird et al., 2011 to Davis et al., 2018
feedstock_price = 71.26 / _kg_per_ton * 0.8 

# 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)
# (https://www.eia.gov/outlooks/aeo/), which is $0.7328/gal and similar to the
# 2.2/(2988/1e3) = $0.736/gal based on a density of 2988 g/gal from H2 Tools
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg
_ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol
_ethanol_MW = chems.Ethanol.MW
_ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
ethanol_price = 2.2 / _ethanol_kg_2_gal
	

# Dipotassium hydrogen phosphate (DPHP)
# https://www.alibaba.com/product-detail/Food-Grade-Dipotassium-Dydrogen-Phosphate-Trihydrate_60842047866.html?spm=a2700.7724857.normalList.2.4ef2457e3gPbfv&s=p
# DISREGARD: https://www.sigmaaldrich.com/catalog/product/mm/105104?lang=en&region=US

DPHP_price = 1.15

# 2.86 is the average motor gasoline price between 2010-2019 in 2016 $/gal	
# based on AEO from EIA, density of gasoline is 2.819 kg/gal	
# based on Lower and Higher Heating Values of Hydrogen and Other Fuels	
# from H2 Tools maintained by Pacific Northwest National Laboratory	
# (https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels)	
denaturant_price = 2.86 / 2.819

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of Davis et al., 2018 (TAL scenario)
# 7880 is operating hours/yr on Page 10 of Davis et al., 2018,
# cost is negative because it's a product stream
ash_disposal_price = -1.41e6 / (4279*7880)

# Assums no cost/credit for baseline, the same as ash disposal for the lower bound,	
# for the upper bound (i.e., positive selling price indicating profit), use 	
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
# For the lower bound (i.e., negative selling price indicating cost), use price from
# Aden et al., 2002: $0.0094/lb in 2000$ = 0.0094*1.114/0.802*2.20462 = $0.0288/kg
# in 2016$
gypsum_price = 0

# Baseline from Davis et al., 2018, lower bound is 2015-2019 average of 	
# hydrate lime in $/ton at plant from Mineral Commodity Summaries 2020.	
# 2015: 146.40 * (1.114/1.100) / 907.18474 = 0.163	
# 2016: 145.50 / 907.18474 = 0.160	
# 2017: 147.10 * (1.114/1.134) / 907.18474 = 0.159	
# 2018: 151.50 * (1.114/1.157) / 907.18474 = 0.161	
# 2019: 151.00 * (1.114/1.185) / 907.18474 = 0.156	
# (0.163+0.160+0.159+0.161+0.156) / 5 = 0.160	
# Upper bound is +10% from baseline = 0.1189 * _lb_per_kg * 1.1 = 0.288
lime_price = 0.1189 * _lb_per_kg

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_price = 466833 / 5 / (24*365*0.96)

# 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# Using the same conversion, lower and upper bounds should be 
# $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chems.CH4.MW
natural_gas_price = 4.70/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW)

# https://www.rightpricechemicals.com/buy-amberlyst-15-ion-exchange-resin.html	
# USD 383.13 for 2.5kg (largest available size order), accessed 06/11/2020
amberlyst_15_price = 153.252 * _chemical_2020to2016

# https://www.alibaba.com/product-detail/Tricalcium-Phosphate-Tricalcium-Phosphate-TCP-Tricalcium_60744013678.html?spm=a2700.galleryofferlist.0.0.42f16684C9iJhz&s=p
TCP_price = 850 / _kg_per_ton # tricalcium (di)phosphate


# TAL_price = 1.88 # initial value
SA_price = 1.88 # initial value

# Currently not in use
# Methanol price from Goellner et al., production from natural gas (Case 3),
# average of two load structures in 2011$,
# crude methanol with ~1% CO2 and 1% H2O
methanol_price = (311.17+345.39)/2/_kg_per_ton * _chemical_2011to2016


# Acetoin product selling price
# https://www.alibaba.com/product-detail/Acetoin-CAS-NO-513-86-0_60673118759.html?spm=a2700.galleryofferlist.0.0.4d906f82dIhSkn
# acetoin_price = 5.925
acetoin_price = 3. # assumed
# Isobutyraldehyde product selling price
# https://www.alibaba.com/product-detail/China-manufacture-Isobutyraldehyde-price_60837108075.html?spm=a2700.galleryofferlist.0.0.753369fcbZcNhe
# IBA_price = 1.2
IBA_price = 0. # assumed

TAL_price = 4. # assumed; when solving for MPSP, this is merely the initial value and has no effect on results

glucose_price = 236. /_kg_per_ton # refer to email: "sugar price from maravelias group"

# https://www.alibaba.com/product-detail/Manufacturer-of-Hexyl-Alcohol-Hexanol-n_60403061175.html?spm=a2700.galleryofferlist.0.0.1021992cx1VYY8
hexanol_price = 6.

# https://www.alibaba.com/product-detail/N-heptane-Heptane-Heptane-Heptane-Supply_62451341262.html?spm=a2700.galleryofferlist.0.0.2b407553zLGZvf&s=p
heptane_price = 1900./_kg_per_ton

# Viswanathan et al. 2020
toluene_price = 1.1

# assumed
AuPd_price = 60000.

# https://www.alibaba.com/product-detail/Gas-Hydrogen-45kg-Lpg-Gas-Cylinder_62018626105.html?spm=a2700.galleryofferlist.0.0.e9ba7ce2O0TyvK
hydrogen_price = 1.

# https://www.energy.gov/eere/fuelcells/hydrogen-shot
hydrogen_renewable_price = 5.
# https://www.alibaba.com/product-detail/catalyst-raney-nickel-stainless-steel-square_60731133248.html?spm=a2700.galleryofferlist.0.0.2a583524IfakZZ
RaneyNi_price = 30.

NiSiO2_price = RaneyNi_price

# https://www.alibaba.com/product-detail/Reagent-Grade-90-caustic-potash-potassium_62118969650.html?spm=a2700.galleryofferlist.0.0.28555ed4pKlEVC&s=p
KOH_price = 1.6

# https://www.alibaba.com/product-detail/Hydrochloric-acid-HCl-7647-01-0_60439085052.html?spm=a2700.galleryofferlist.0.0.4fa42c515nP2GH
HCl_price = 0.3

activated_carbon_price = 41. # $/ft^3 # Seader et al.


PdC_price = 0.075*(2045./0.0311035) \
      + (1-0.075)*0.45 # Pd : 2045 EIB (USD/troy-ounce) # https://apps.catalysts.basf.com/apps/eibprices/mp/ (accessed 6/3/2022)
                       # activated carbon: $0.45/kg # https://www-sciencedirect-com.proxy2.library.illinois.edu/science/article/pii/S2590174522000411?via%3Dihub
spent_PdC_price = 1. # assumed

acetone_price = 0.63 * _GDP_2008_to_2016 * _lb_per_kg # average of range ($0.44 - $0.82 /lb) from https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/

acetic_acid_price = 1.135 * _GDP_2008_to_2016 # average of range ($ 0.772 - 1.499 /kg) from # https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
sodium_acetate_price = acetic_acid_price # unused

amberlyst70_price = 30. #!!!

acetylacetone_price = 1.5 # 2,4-pentanedione or acetylacetone

price = {'SA': SA_price,
         'PD': acetylacetone_price, # 2,4-pentanedione or acetylacetone
         'TCP': TCP_price,
         'AuPd': AuPd_price,
         'IBA': IBA_price,
         'TAL': TAL_price,
         'KOH': KOH_price,
         'HCl': HCl_price,
         'Hydrogen': hydrogen_price,
         'Renewable hydrogen': hydrogen_renewable_price,
         'Acetoin': acetoin_price,
         'RaneyNi': RaneyNi_price,
         'Ni-SiO2': NiSiO2_price,
         'Amberlyst-70': amberlyst70_price,
         'PdC': PdC_price,
         'Spent PdC': spent_PdC_price,
         'Feedstock': feedstock_price,
         'Glucose': glucose_price,
         'Hexanol': hexanol_price,
         'Heptane': heptane_price,
         'Toluene': toluene_price,
         'Acetone': acetone_price,
         'Sulfuric acid': 0.0430 * _lb_per_kg,	
         # 0.1900 is for NH3	
         'AmmoniumHydroxide': 0.1900 * _lb_per_kg * 17.031/35.046,	
         'CSL': 0.0339 * _lb_per_kg,
         'Caustics': 0.2384 * _lb_per_kg * 0.5, # 50 wt% NaOH/water mixture	
         'Boiler chems': 2.9772 * _lb_per_kg,	
         'Lime': lime_price,
         'Cooling tower chems': 1.7842 * _lb_per_kg,	
         'Makeup water': 0.0002 * _lb_per_kg,	
         # Cost of ash is negative because it's a product stream	
         'Ash disposal': ash_disposal_price,
         'Electricity': 0.070, # AEO from EIA, 2010-2019 average (0.067-0.074 range)	
         # $6.16/kg protein in 2016$, P25 of Davis et al., 2018	
         'Enzyme': 6.16,
         'DPHP': DPHP_price,
         'Baghouse bag': baghouse_bag_price,	
         'Natural gas': natural_gas_price,
         'Methanol': methanol_price,
         'Ethanol': ethanol_price,
         # Below currently not in use
         'Gypsum': gypsum_price,
         'Denaturant': denaturant_price,
         'Amberlyst15': amberlyst_15_price,
         'DAP': 0.1645 * _lb_per_kg,
         'Activated carbon': activated_carbon_price,
         'Sodium acetate': sodium_acetate_price,
         'Acetic acid': acetic_acid_price,
         }
    
#!!! Round all prices to 4 *decimal places*
for k in price.keys():
    price[k] = round(price[k], 4)

bst.PowerUtility.price = price['Electricity']

_lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
_mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
_hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
_mps.T = 233 + 273.15
_hps.T = 266 + 273.15

_cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
_chilled = bst.HeatUtility.get_cooling_agent('chilled_water')
_cooling.regeneration_price = 0
_cooling.T = 28 + 273.15
_cooling.T_limit = _cooling.T + 9

# Side steam in CHP not a heat utility, thus will cause problem in TEA utility
# cost calculation if price not set to 0 here, costs for regeneration of heating
# and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
for i in (_lps, _mps, _hps, _cooling, _chilled):
    i.heat_transfer_price = i.regeneration_price = 0
    # if i == _cooling: continue
    # i.heat_transfer_efficiency = 0.85


# %%

# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA), all from ref [5] if not noted otherwise
# =============================================================================

CFs = {}

# =============================================================================
# 100-year global warming potential (GWP) in kg CO2-eq/kg
# =============================================================================
GWP_CFs = {
    'CH4': 0.40, # NA NG from shale and conventional recovery
    'CSL': 1.55,
    
    # 'Enzyme': 2.24, 
    'Ethanol': 1.44,
    'Acetone': 2.5435, #  ecoinvent 3.8 market for acetone, liquid, RoW
    
    'H2SO4': 44.47/1e3,   
    'Lime': 1.29 * 56.0774/74.093, # CaO to Ca(OH)2
    'CaO': 1.29,
    'NaOH': 2.11,
    # 'NH4OH': 2.64 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,   
    # 'MEA': 3.4062, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 1.3598, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 0.87104, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
    'H2': 2.3716, # ecoinvent 3.8 market for hydrogen, liquid, RoW
    
    'SodiumAcetate': 1., # !!! update
    'AceticAcid': 1.6198, # market for acetic acid, without water, in 98% solution state, GLO
    
    'PD': (2*(58.080)*3.5917 + 102.089*2.5435)/(2*100.117), # Acetylacetone; based on GLO/RoW IPCC 2013 CFs of markets for precursors acetone and acetyl anhydride
    # 'DiammoniumSulfate': 1.2901, # ecoinvent 3.8 market for ammonium sulfate, RoW
    # 'MagnesiumSulfate': 1.0411, # ecoinvent 3.8 market for magnesium sulfate, GLO
    
    'NiSiO2':10., # !!! update
    'Amberlyst70_':10., # !!! update
    }



GWP_CF_array = chems.kwarray(GWP_CFs)
# In kg CO2-eq/kg of material
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')
CFs['GWP_CF_stream'] = GWP_CF_stream

GWP_CFs['Electricity'] = 0.4490 # kg CO2-eq/kWh GREET 2022 US Mix  # assume production==consumption, both in kg CO2-eq/kWh


GWP_CFs['Sugarcane'] = 0.12043 * 0.3/0.286 # ecoinvent 3.8 market for sugarcane, RoW
# # adjusted from dry wt content of 28.6% (their assumption) to 30% (our assumption)


# GWP_CFs['Sugarcane'] = 0.02931 # GREET 2022


CFs['GWP_CFs'] = GWP_CFs

# =============================================================================
# Fossil energy consumption (FEC), in MJ/kg of material
# =============================================================================


FEC_CFs = {
    'CH4': 50, # NA NG from shale and conventional recovery
    'CSL': 12,
    
    'Ethanol': 16,
    'Acetone': 66.852, #  ecoinvent 3.8 market for acetone, liquid, RoW
    # 'Enzyme': 26,
    
    'H2SO4': 568.98/1e3,
    'Lime': 4.896 * 56.0774/74.093, # CaO to Ca(OH)2
    'CaO': 4.896, 
    'NaOH': 29,
    # 'NH4OH': 42 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,
    # 'MEA': 67.898, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 16.538, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 7.4243, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
    'H2': 75.747, # ecoinvent 3.8 market for hydrogen, liquid, RoW

    'SodiumAcetate': 1., # !!! update
    'AceticAcid': 45.611, # market for acetic acid, without water, in 98% solution state, GLO
    
    'PD': (2*(58.080)*66.852 + 102.089*70.817)/(2*100.117), # Acetylacetone; based on GLO/RoW cumulative energy demand CFs of markets for precursors acetone and acetyl anhydride
    
    # 'DiammoniumSulfate': 15.166, # ecoinvent 3.8 market for ammonium sulfate, RoW
    # 'MagnesiumSulfate': 13.805, # ecoinvent 3.8 market for magnesium sulfate, GLO
    
    'NiSiO2':10., # !!! update
    'Amberlyst70_':10., # !!! update
    
    }

FEC_CF_array = chems.kwarray(FEC_CFs)
# In MJ/kg of material
FEC_CF_stream = tmo.Stream('FEC_CF_stream', FEC_CF_array, units='kg/hr')

CFs['FEC_CF_stream'] = FEC_CF_stream

FEC_CFs['Electricity'] = 5.724 # MJ/kWh # GREET 2022 US Mix #assume production==consumption, both in MJ/kWh

FEC_CFs['Sugarcane'] = 	0.40192 * 0.3/0.286 # ecoinvent 3.8 market for sugarcane, RoW
# # adjusted from dry wt content of 28.6% (their assumption) to 30% (our assumption)

# FEC_CFs['Sugarcane'] = 	0.2265 # GREET 2022

CFs['FEC_CFs'] = FEC_CFs
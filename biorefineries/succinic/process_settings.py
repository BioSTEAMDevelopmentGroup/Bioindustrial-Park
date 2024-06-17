#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat

References
----------
[1] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234
[2] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.
[3] Argonne National Laboratory. The Greenhouse gases, Regulated Emissions,
    and Energy use in Transportation (GREET) Model https://greet.es.anl.gov/
    (accessed Aug 25, 2020).
[4] Roni et al., Herbaceous Feedstock 2018 State of Technology Report;
    INL/EXT-18-51654-Rev000; Idaho National Lab. (INL), 2020.
    https://doi.org/10.2172/1615147.
[5] ecoinvent 3.6 https://www.ecoinvent.org/home.html (accessed Aug 26, 2020).

"""

import biosteam as bst
import thermosteam as tmo
from biorefineries.succinic.chemicals_data import chems

bst.CE = 541.7 # year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147

_GDP_2007_to_2016 = 1.114 / 0.961
_chemical_2011to2016 = 102.5 / 91.7
_chemical_2013to2016 = 102.5 / 101.3
_chemical_2014to2016 = 102.5 / 105.3

_chemical_2017to2016 = 102.5 / 106.9
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb
_chemical_2022to2016 = 102.5 / 145.3

# From sugarcane biorefinery
feedstock_price = 0.03455

# 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)
# (https://www.eia.gov/outlooks/aeo/), which is $0.7328/gal and similar to the
# 2.2/(2988/1e3) = $0.736/gal based on a density of 2988 g/gal from H2 Tools
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg
_ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol
_ethanol_MW = chems.Ethanol.MW
_ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
ethanol_price = 2.2 / _ethanol_kg_2_gal
	

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

# Assumes no cost/credit for baseline, the same as ash disposal for the lower bound,	
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

# 365-380 2022$/MT according to https://www.chemanalyst.com/Pricing-data/liquid-carbon-dioxide-1090
# Baseline: 0.263 2016$/kg
# Uncertainty range: 0.257 - 0.268 2016$/kg
liquid_CO2_price = 0.3725 * _chemical_2022to2016


succinic_acid_price = ((2.86*38e3 + 2.50*40e3)/(38e3+40e3))*_chemical_2013to2016 # $2.75/kg in 2016$ # initial value to solve MPSP; when solving NPV, use this as the selling price

# Monoethanol amine
# 1021.69 $/MT in Nov 2014 - 1855.65 $/MT in Oct 2017 range from https://www.intratec.us/chemical-markets/monoethanolamine-price
# Baseline: 1.426 2016$/kg
# Uncertainty range: 1.032 - 1.819 2016$/kg
MEA_price = (1021.69*_chemical_2014to2016 + 1855.65*_chemical_2017to2016)/(2*1000)


# Diammonium sulfate (ammonium sulfate)
# mean of 0.161 (range 0.154-0.167)	$/kg in 2007 $	
# mean of 0.187 (range of 0.178-0.194) $/kg when converted to 2016$
# https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/
diammonium_sulfate_price = 0.161 * _GDP_2007_to_2016 


# Magnesium sulfate
# mean of 0.436 (range 0.397-0.474)	$/kg in 2007 $	
# mean of 0.505 (range of 0.460-0.549) $/kg when converted to 2016$
# https://web.archive.org/web/20161125084558/http://www.icis.com:80/chemicals/channel-info-chemicals-a-z/                        
magnesium_sulfate_price = 0.436 * _GDP_2007_to_2016

# All in 2016$/kg
price = {
         'Feedstock': feedstock_price,
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
         'Enzyme': 6.16, # $6.16/kg protein in 2016$, P25 of Davis et al., 2018	
         'Baghouse bag': baghouse_bag_price,	
         'Natural gas': natural_gas_price,
         'Ethanol': ethanol_price,
         'Liquid carbon dioxide': liquid_CO2_price,
         'Gypsum': gypsum_price,
         'Monoethanolamine': MEA_price,
         'Succinic acid': succinic_acid_price, # initial value to solve MPSP
         # Below currently not in use
         'Denaturant': denaturant_price,
         'Diammonium sulfate': diammonium_sulfate_price,
         'Magnesium sulfate': magnesium_sulfate_price,
         
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
    # 'Ethanol': 1.44,
    
    'H2SO4': 44.47/1e3,   
    'Lime': 1.29 * 56.0774/74.093, # CaO to Ca(OH)2
    'NaOH': 2.11,
    'NH4OH': 2.64 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,   
    'MEA': 3.4062, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 1.3598, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 0.87104, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
    'DiammoniumSulfate': 1.2901, # ecoinvent 3.8 market for ammonium sulfate, RoW
    'MagnesiumSulfate': 1.0411, # ecoinvent 3.8 market for magnesium sulfate, GLO
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
    
    # 'Ethanol': 16,
    # 'Enzyme': 26,
    
    'H2SO4': 568.98/1e3,
    'Lime': 4.896 * 56.0774/74.093, # CaO to Ca(OH)2
    'NaOH': 29,
    'NH4OH': 42 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,
    'MEA': 67.898, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 16.538, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 7.4243, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
    'DiammoniumSulfate': 15.166, # ecoinvent 3.8 market for ammonium sulfate, RoW
    'MagnesiumSulfate': 13.805, # ecoinvent 3.8 market for magnesium sulfate, GLO
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





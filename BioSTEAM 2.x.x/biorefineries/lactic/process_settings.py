#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon Dec 30 09:30:10 2019

References:
[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted,
    2020.

[3] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[4] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.
    
[5] Neupane et al., Life-Cycle Greenhouse Gas and Water Intensity of Cellulosic
    Biofuel Production Using Cholinium Lysinate Ionic Liquid Pretreatment.
    ACS Sustainable Chem. Eng. 2017, 5 (11), 10176–10185.
    https://doi.org/10.1021/acssuschemeng.7b02116.

[6] Nordahl et al., Life-Cycle Greenhouse Gas Emissions and Human Health
    Trade-Offs of Organic Waste Management Strategies. Environ. Sci. Technol. 2020.
    https://doi.org/10.1021/acs.est.0c00364.

[7] Jiménez Rivero et al., Life Cycle Energy and Material Flow Implications
    of Gypsum Plasterboard Recycling in the European Union.
    Resources, Conservation and Recycling 2016, 108, 171–181.
    https://doi.org/10.1016/j.resconrec.2016.01.014.

[8] U.S. Energy Information Administration. U.S. electricity generation by energy source.
    https://www.eia.gov/tools/faqs/faq.php (accessed Aug 11, 2020).

[9] Roni et al., Herbaceous Feedstock 2018 State of Technology Report;
    INL/EXT-18-51654-Rev000; Idaho National Lab. (INL), 2020.
    https://doi.org/10.2172/1615147.

@author: yalinli_cabbi
"""


# %%

# =============================================================================
# Setup
# =============================================================================

import biosteam as bst
import thermosteam as tmo
from lactic.chemicals import chems


# %%

# =============================================================================
# Energy balances
# =============================================================================

_lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
_mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
_hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
_mps.T = 233 + 273.15
_hps.T = 266 + 273.15

_cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
_cooling.regeneration_price = 0
_cooling.T = 28 + 273.15
_cooling.T_limit = _cooling.T + 9

# Side steam in CHP not a heat utility, thus will cause problem in TEA utility
# cost calculation if price not set to 0 here, costs for regeneration of heating
# and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
for i in (_lps, _mps, _hps, _cooling):
    i.heat_transfer_price = i.regeneration_price = 0
    # if i == _cooling: continue
    # i.heat_transfer_efficiency = 0.85


# %%

# =============================================================================
# Prices for techno-economic analysis (TEA), all in $/kg (electricity in $/kWh)
# and from ref [3] if not noted
# =============================================================================

bst.CE = 541.7 # year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb
_GDP_2007to2016 = 1.160

# From USD/dry-ton to USD/kg in 2016$, 20% moisture content,
# baseline and lower bound (60) from ref [3], upper bound (146.4) from
# Hartley et al., ACS Sustainable Chem. Eng. 2020, 8 (19), 7267–7277.
# https://doi.org/10.1021/acssuschemeng.9b06551.
_feedstock_factor = _kg_per_ton / 0.8
feedstock_price = 71.3 / _feedstock_factor
		
# Baseline from ref [3], lower and upper bounds (96% and 110% of baseline)
# calculated using the index of sulfuric acid from U.S. Bureau of Labor Statistics
# (Producer Price Index by Commodity for Chemicals and Allied Products)
# https://fred.stlouisfed.org/series/WPU0613020T1 (accessed Jul 31, 2020).
sulfuric_acid_price = 0.0430 * _lb_per_kg

# 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 	
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)	
# (https://www.eia.gov/outlooks/aeo/), which is $0.732/gal and similar to the 	
# 2.2/(2988/1e3) = $0.736/gal based on a density of 2988 g/gal from H2 Tools
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg
_ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol	
_ethanol_MW = chems.Ethanol.MW
_ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
ethanol_price = 2.2 / _ethanol_kg_2_gal

# Cost is negative because it's a product stream
ash_disposal_price = -1.41e6 / (4279*7880)

# Assums no cost/credit for baseline, the same as ash disposal for the lower bound,	
# for the upper bound (i.e., positive selling price indicating profit), use 	
# USGS 2015-2019 average  free on bard price in $/metric tonne for crude gypsum. 
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
# ref [4]: $0.0094/lb in 2000$ = 0.0094*1.114/0.802*2.20462 = $0.0288/kg in 2016$
gypsum_price = 0

# Baseline from ref [3], lower bound is 2015-2019 average of 	
# hydrate lime in $/ton at plant from Mineral Commodity Summaries 2020.	
# 2015: 146.40 * (1.114/1.100) / 907.18474 = 0.163	
# 2016: 145.50 / 907.18474 = 0.160	
# 2017: 147.10 * (1.114/1.134) / 907.18474 = 0.159	
# 2018: 151.50 * (1.114/1.157) / 907.18474 = 0.161	
# 2019: 151.00 * (1.114/1.185) / 907.18474 = 0.156	
# (0.163+0.160+0.159+0.161+0.156) / 5 = 0.160	
# Upper bound is +10% from baseline = 0.1189 * _lb_per_kg * 1.1 = 0.288
lime_price = 0.1189 * _lb_per_kg

baghouse_bag_price = 466183/5/(24*365*0.96) * _GDP_2007to2016

# 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# Using the same conversion, lower and upper bounds (min/max of 2010-2019) should be 
# $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
_CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
_CH4_MW = chems.CH4.MW
natural_gas_price = 4.70/1e3*_ft3_per_m3*_CH4_V * (1e3/_CH4_MW)

# https://www.rightpricechemicals.com/buy-amberlyst-15-ion-exchange-resin.html	
# USD 383.13 for 2.5kg (largest available size order), accessed 06/11/2020
amberlyst_15_price = 153.252 * _chemical_2020to2016

# All in 2016$/kg
price = {'Feedstock': feedstock_price, 	
         'H2SO4': 0.0430 * _lb_per_kg,	
         # 0.1900 is for NH3	
         'NH4OH': 0.1900 * _lb_per_kg * chems.NH3.MW/chems.NH4OH.MW,
         'CSL': 0.0339 * _lb_per_kg,
         'Enzyme': 6.16,
         'Lime': lime_price,
         'Amberlyst15': amberlyst_15_price,
         'NaOH': 0.2384 * _lb_per_kg,
         'WWT polymer': 2.6282 * _lb_per_kg,
         'Boiler chems': 2.9772 * _lb_per_kg,	
         'Cooling tower chems': 1.7842 * _lb_per_kg,	
         'Makeup water': 0.0002 * _lb_per_kg,	
         # Cost of ash is negative because it's a product stream	
         'Ash disposal': ash_disposal_price,	
         'Gypsum': gypsum_price,
         'Ethanol': ethanol_price,	
         'Baghouse bag': baghouse_bag_price,	
         'Natural gas': natural_gas_price}

# $/kWh, from EIA AEO, 2010-2019 average (0.067-0.074 range)
bst.PowerUtility.price = 0.070


# %%

# =============================================================================
# Global warming potential (GWP) for life cycle analysis (LCA), all from ref [5] 
# if not noted
# =============================================================================

# _CH4_MJ_per_kg = - chems.CH4.HHV/1e6/(chems.CH4.MW/1e3)

# In kg CO2-eq/kg
# GWP_CFs = {
#     'NH4OH': 1.987 * chems.NH3.MW/chems.NH4OH.MW,
#     'CSL': 0.1105,
#     'CH4': 0.00101 * _CH4_MJ_per_kg,
#     'Enzyme': 0.1128,
#     # From ref [6], 0.94 kg CO2 and 6e-7 kg N2O per kg of lime
#     'Lime': 0.94+6e-7*298,
#     'NaOH': 0.2186,
#     'H2SO4': 1.36e-5,
#     'Gypsum': 0
#     # 'Gypsum': -2.03/1e3 # ref[7]
#     }


_CH4_m3_per_kg = _CH4_V/_CH4_MW * 1e3

# In kg CO2-eq/kg
GWP_CFs = {
    'NH4OH': 2.0462 * chems.NH3.MW/chems.NH4OH.MW,
    'CSL': 2.7797, # fodder yeast
    'CH4': 0.40519 * _CH4_m3_per_kg,
    'Enzyme': 12.237,
    'Lime': 0.96303,
    'NaOH': 1.2909,
    'H2SO4': 0.13803,
    'Ethanol': 0.79329
    }

GWP_CF_array = chems.kwarray(GWP_CFs)

# Actually in kg CO2-eq/kg
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')

# In kg CO2-eq/kWh
GWP_CF_electricity = 0.58019

# In kg CO2-eq/kg from ref [9]
GWP_CF_feedstock = 66.2 / _feedstock_factor


# =============================================================================
# Cumulative energy demand - fossil (CEDf)
# =============================================================================


CEDf_CFs = {
    'NH4OH': 37.89 * chems.NH3.MW/chems.NH4OH.MW,
    'CSL': 19.058, # fodder yeast
    'CH4': 42.378 * _CH4_m3_per_kg,
    'Enzyme': 136.28,
    'Lime': 4.7698,
    'NaOH': 15.128,
    'H2SO4': 4.1051,
    'Ethanol': 5.0354
    }

CEDf_CF_array = chems.kwarray(CEDf_CFs)

# Actually in kg CO2-eq/kg
CEDf_CF_stream = tmo.Stream('CEDf_CF_stream', CEDf_CF_array, units='kg/hr')

# In kg CO2-eq/kWh
CEDf_CF_electricity = 6.9438















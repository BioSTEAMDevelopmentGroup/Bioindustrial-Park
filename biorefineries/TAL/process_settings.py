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
_chemical_2011to2016 = 102.5 / 91.7
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb

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

# https://www.alibaba.com/product-detail/Reagent-Grade-90-caustic-potash-potassium_62118969650.html?spm=a2700.galleryofferlist.0.0.28555ed4pKlEVC&s=p
KOH_price = 1.6

# https://www.alibaba.com/product-detail/Hydrochloric-acid-HCl-7647-01-0_60439085052.html?spm=a2700.galleryofferlist.0.0.4fa42c515nP2GH
HCl_price = 0.3

activated_carbon_price = 41. # $/ft^3 # Seader et al.
# All in 2016$/kg

PdC_price = 0.075*(2045./0.0311035) \
      + (1-0.075)*0.45 # Pd : 2045 EIB (USD/troy-ounce) per https://apps.catalysts.basf.com/apps/eibprices/mp/ (accessed 6/3/2022)
                       # activated carbon: $0.45/kg # https://www-sciencedirect-com.proxy2.library.illinois.edu/science/article/pii/S2590174522000411?via%3Dihub
spent_PdC_price = 1. # assumed

price = {'SA': SA_price,
         'PD': 1.5, # 2,4-pentanedione or acetylacetone
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
         'PdC': PdC_price,
         'Spent PdC': spent_PdC_price,
         'Feedstock': feedstock_price,
         'Glucose': glucose_price,
         'Hexanol': hexanol_price,
         'Heptane': heptane_price,
         'Toluene': toluene_price,
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
         'Activated carbon': activated_carbon_price}
    
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
# Characterization factors (CFs) for life cycle analysis (LCA), all from ref [5] 
# if not noted, note that it is unclear if in-plant receiving and preprocessing
# (~50% of the total impact per ref [6]) of feedstock is included in ref [5]
# =============================================================================

CFs = {}

# =============================================================================
# 100-year global warming potential (GWP) in kg CO2-eq/kg
# =============================================================================
GWP_CFs = {
    'NH4OH': 2.64 * chems.NH3.MW/chems.NH4OH.MW,
    'CSL': 1.55,
    'CH4': 0.40, # NA NG from shale and conventional recovery
    'Enzyme': 2.24,
    'Lime': 1.29,
    'NaOH': 2.11,
    'H2SO4': 44.47/1e3,
    'Ethanol': 1.44,
    }

GWP_CF_array = chems.kwarray(GWP_CFs)
# In kg CO2-eq/kg of material
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')

GWP_CFs['Corn stover'] = 44.70/1e3 * 0.8
GWP_CFs['Switchgrass'] = 87.81/1e3 * 0.8
GWP_CFs['Miscanthus'] = 78.28/1e3 * 0.8
GWP_CFs['CaCO3'] = 10.30/1e3
GWP_CFs['Gypsum'] = -4.20/1e3
# In kg CO2-eq/kWh
GWP_CFs['Electricity'] = 0.48
# From corn stover
GWP_CFs['LacticAcid_GREET'] = 1.80
# From ref [7], lactic acid production, RoW, TRACI global warming
GWP_CFs['LacticAcid_fossil'] = 4.1787

CFs['GWP_CFs'] = GWP_CFs
CFs['GWP_CF_stream'] = GWP_CF_stream

# =============================================================================
# Fossil energy consumption (FEC), in MJ/kg of material
# =============================================================================

FEC_CFs = {
    'NH4OH': 42 * chems.NH3.MW/chems.NH4OH.MW,
    'CSL': 12,
    'CH4': 50, # NA NG from shale and conventional recovery
    'Enzyme': 26,
    'Lime': 4.896,
    'NaOH': 29,
    'H2SO4': 568.98/1e3,
    'Ethanol': 16
    }

FEC_CF_array = chems.kwarray(FEC_CFs)
# In MJ/kg of material
FEC_CF_stream = tmo.Stream('FEC_CF_stream', FEC_CF_array, units='kg/hr')

FEC_CFs['Corn stover'] = 688.60/1e3 * 0.8
FEC_CFs['Switchgrass'] = 892.41/1e3 * 0.8
FEC_CFs['Miscanthus'] = 569.05/1e3 * 0.8
FEC_CFs['CaCO3'] = 133.19/1e3
FEC_CFs['Gypsum'] = -44.19/1e3
# In MJ/kWh
FEC_CFs['Electricity'] = 5.926
# From corn stover
FEC_CFs['LacticAcid'] = 29
# From ref [7], lactic acid production, RoW, cumulative energy demand, fossil
FEC_CFs['LacticAcid_fossil'] = 79.524

CFs['FEC_CFs'] = FEC_CFs
CFs['FEC_CF_stream'] = FEC_CF_stream






# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# # BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# # Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# # Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# # Copyright (C) 2020, Yalin Li <mailto.yalin.li@gmail.com>,
# # Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# # 
# # This module is under the UIUC open-source license. See 
# # github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# # for license details.

# """
# Created on Mon Dec 30 09:30:10 2019
# Modified from the biorefineries constructed in [1] and [2] for the production of
# lactic acid from lignocellulosic feedstocks
# [1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
#     Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
#     ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
#     https://doi.org/10.1021/acssuschemeng.9b07040
    
# [2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
#     Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted.
#     July, 2020.
# [3] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
#     Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
#     NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
#     https://doi.org/10.2172/1483234
# [4] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
#     Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
#     Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
#     National Renewable Energy Lab (NREL), 2002.
#     https://doi.org/10.2172/1218326.
# @author: yalinli_cabbi
# """


# # %%

# import biosteam as bst
# from lactic.chemicals import chems

# bst.CE = 541.7 # year 2016
# _kg_per_ton = 907.18474
# _lb_per_kg = 2.20462
# _liter_per_gallon = 3.78541
# _ft3_per_m3 = 35.3147
# _chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb

# # From USD/dry-ton to USD/kg in 2016$, 20% moisture content,
# # baseline and lower bound (60) from ref [3], upper bound (146.4) from
# # Hartley et al., ACS Sustainable Chem. Eng. 2020, 8 (19), 7267–7277.
# # https://doi.org/10.1021/acssuschemeng.9b06551.
# _feedstock_factor = _kg_per_ton / 0.8
# feedstock_price = 71.3 / _feedstock_factor
# 		
# # Baseline from ref [3], lower and upper bounds (96% and 110% of baseline)
# # calculated using the index of sulfuric acid from U.S. Bureau of Labor Statistics
# # (Producer Price Index by Commodity for Chemicals and Allied Products)
# # https://fred.stlouisfed.org/series/WPU0613020T1 (accessed Jul 31, 2020).
# sulfuric_acid_price = 0.0430 * _lb_per_kg

# # 2.2 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 	
# # based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)	
# # (https://www.eia.gov/outlooks/aeo/), which is $0.732/gal and similar to the 	
# # 2.2/(2988/1e3) = $0.736/gal based on a density of 2988 g/gal from H2 Tools
# # Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg
# _ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol	
# _ethanol_MW = chems.Ethanol.MW
# _ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
# ethanol_price = 2.2 / _ethanol_kg_2_gal

# # Cost is negative because it's a product stream
# ash_disposal_price = -1.41e6 / (4279*7880)

# # Assums no cost/credit for baseline, the same as ash disposal for the lower bound,	
# # for the upper bound (i.e., positive selling price indicating profit), use 	
# # USGS 2015-2019 average  free on bard price in $/metric ton for crude gypsum. 
# # National Minerals Information Center. Mineral Commodity Summaries 2020; 
# # U.S. Geological Survey, 2020.
# # Assuming all prices were in their nominal year (e.g., 2015 price in 2015$)
# # and adjusted to 2016$
# # 2015: 7.80 * 1.114 / 1.100 = 7.90
# # 2016: 8.00
# # 2017: 7.50 * 1.114 / 1.134 = 7.37
# # 2018: 8.30 * 1.114 / 1.157 = 7.99
# # 2019: 8.00 * 1.114 / 1.185 = 7.52
# # (7.90+8.00+7.37+7.99+7.52) / 5 = 7.76 (in metric tonne)
# # For the lower bound (i.e., negative selling price indicating cost), use price from
# # ref [4]: $0.0094/lb in 2000$ = 0.0094*1.114/0.802*2.20462 = $0.0288/kg in 2016$
# gypsum_price = 0

# # Baseline from ref [3], lower bound is 2015-2019 average of 	
# # hydrate lime in $/ton at plant from Mineral Commodity Summaries 2020.	
# # 2015: 146.40 * (1.114/1.100) / 907.18474 = 0.163	
# # 2016: 145.50 / 907.18474 = 0.160	
# # 2017: 147.10 * (1.114/1.134) / 907.18474 = 0.159	
# # 2018: 151.50 * (1.114/1.157) / 907.18474 = 0.161	
# # 2019: 151.00 * (1.114/1.185) / 907.18474 = 0.156	
# # (0.163+0.160+0.159+0.161+0.156) / 5 = 0.160	
# # Upper bound is +10% from baseline = 0.1189 * _lb_per_kg * 1.1 = 0.288
# lime_price = 0.1189 * _lb_per_kg

# baghouse_bag_price = 466833 / 5 / (24*365*0.96)

# # 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# # which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# # similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# # Using the same conversion, lower and upper bounds (min/max of 2010-2019) should be 
# # $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
# _CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
# _CH4_MW = chems.CH4.MW
# natural_gas_price = 4.70/1e3*_ft3_per_m3*_CH4_V * (1e3/_CH4_MW)

# # https://www.rightpricechemicals.com/buy-amberlyst-15-ion-exchange-resin.html	
# # USD 383.13 for 2.5kg (largest available size order), accessed 06/11/2020
# amberlyst_15_price = 153.252 * _chemical_2020to2016

# # All in 2016$/kg
# price = {'Feedstock': feedstock_price, 	
#          'Sulfuric acid': 0.0430 * _lb_per_kg,	
#          # 0.1900 is for NH3	
#          'NH4OH': 0.1900 * _lb_per_kg * 17.031/35.046,	
#          'CSL': 0.0339 * _lb_per_kg,
#          'Enzyme': 6.16,
#          'Lime': lime_price,
#          'Amberlyst15': amberlyst_15_price,
#          'NaOH': 0.2384 * _lb_per_kg,
#          'WWT polymer': 2.6282 * _lb_per_kg,
#          'Boiler chems': 2.9772 * _lb_per_kg,	
#          'Cooling tower chems': 1.7842 * _lb_per_kg,	
#          'Makeup water': 0.0002 * _lb_per_kg,	
#          # Cost of ash is negative because it's a product stream	
#          'Ash disposal': ash_disposal_price,	
#          'Gypsum': gypsum_price,	
#          'Electricity': 0.070, # AEO from EIA, 2010-2019 average (0.067-0.074 range)	
#          'Ethanol': ethanol_price,	
#          'Baghouse bag': baghouse_bag_price,	
#          'Natural gas': natural_gas_price}
    
# bst.PowerUtility.price = price['Electricity']

# _lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
# _mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
# _hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
# _mps.T = 233 + 273.15
# _hps.T = 266 + 273.15

# _cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
# _cooling.regeneration_price = 0
# _cooling.T = 28 + 273.15
# _cooling.T_limit = _cooling.T + 9

# # Side steam in CHP not a heat utility, thus will cause problem in TEA utility
# # cost calculation if price not set to 0 here, costs for regeneration of heating
# # and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
# for i in (_lps, _mps, _hps, _cooling):
#     i.heat_transfer_price = i.regeneration_price = 0
#     # if i == _cooling: continue
#     # i.heat_transfer_efficiency = 0.85

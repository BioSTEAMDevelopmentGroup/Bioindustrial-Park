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
import thermosteam as tmo
import biosteam as bst
from biorefineries.HP.chemicals_data import HP_chemicals as chems
tmo.settings.set_thermo(chems)
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147
_chemical_2011to2016 = 102.5 / 91.7
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb

# From USD/dry-ton to USD/kg in 2016$, 20% moisture content
# changed from Humbird et al., 2011 to Davis et al., 2018
cornstover_price = 71.26 / _kg_per_ton * 0.8 

# https://www.alibaba.com/product-detail/Fresh-sugarcane-Organic-HIGH-QUALITY-LOW_62501985626.html?spm=a2700.galleryofferlist.normal_offer.d_title.58f539c887hNgz
sugarcane_price = 0.3 # !!! temporary

# 2.18 is the average whole-sale ethanol price between 2010-2019 in 2016 $/gal 	
# based on Annual Energy Outlook (AEO) from Energy Information Adiministration (EIA)	
# (https://www.eia.gov/outlooks/aeo/), which is $0.732/gal and similar to the 	
# 2.18/(2988/1e3) = $0.730/gal based on a density of 2988 g/gal from H2 Tools	
# Lower and upper bounds are $1.37/gal and $2.79/gal, or $0.460/kg and $0.978/kg	
ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol	
ethanol_MW = chems.Ethanol.MW	
ethanol_price = 2.18 / (_liter_per_gallon/chems.Ethanol.V('l', 298.15, 101325)*ethanol_MW/1e6)
	

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

# 1.41e6 is $/yr and 4279 in kg/hr from Table 33 of Davis et al., 2018 (HP scenario)
# 7880 is operating hours/yr on Page 10 of Davis et al., 2018,
# price is negative because it's a product stream
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
gypsum_price = 0.

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

# #https://www.alibaba.com/product-detail/Latest-promotion-price-Hot-selling-high_62503274885.html?spm=a2700.7724857.normalList.7.54351dad5oDzYH&s=p&fullFirstScreen=true
# TiO2_price = 1.784

#https://www.alibaba.com/product-detail/Titanium-Dioxide-Chinese-Tio2-Producer-Supply_1600082019592.html?spm=a2700.galleryofferlist.normal_offer.d_title.a3415562TmEgSy&s=p
TiO2_price = 2130./_kg_per_ton

# HP_price = 1.88 # initial value
AA_price = 1.88 # initial value

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


# https://www.alibaba.com/product-detail/Decyl-Alcohol-98-min_62002993466.html?spm=a2700.galleryofferlist.normal_offer.d_title.67b2ac2fnUEnUY
Decanol_price = 1.25

# https://www.alibaba.com/product-detail/Trioctylamine-CAS-NO-1116-76-3_60139027874.html?spm=a2700.galleryofferlist.normal_offer.d_title.2b9b3d2fWfzxb0
TOA_price = 7.5

# https://www.alibaba.com/product-detail/Methyl-Trioctyl-Ammonium-Chloride-Aliquat-336-_50039302076.html?spm=a2700.galleryofferlist.normal_offer.d_title.6f3340baJbm4l2
AQ336_price = 1.115

# https://www.alibaba.com/product-detail/Factory-Supply-High-Purity-Industrial-Grade_1600155716170.html?spm=a2700.galleryofferlist.normal_offer.d_title.68217a2fOdbzfi
Octanol_price = 1.05

# https://www.alibaba.com/product-detail/Best-price-high-quality-Butyl-Acetate_60659284944.html?spm=a2700.galleryofferlist.normal_offer.d_title.8a6f2e61vk4mhP
Butyl_acetate_price = 0.75

# https://www.alibaba.com/product-detail/1-hexanol-cas-111-27-3_1600063760860.html?spm=a2700.galleryofferlist.normal_offer.d_title.16b44add1yPF8c
Hexanol_price = 500. /_kg_per_ton

# https://www.alibaba.com/product-detail/1-6-Hexanediol-99-5-1_1600163052618.html?spm=a2700.galleryofferlist.normal_offer.d_title.4f596babnJcGBc
Hexanediol_price = 2.9

# https://www.alibaba.com/product-detail/Supply-Best-Price-98-1-8_1600152253896.html?spm=a2700.galleryofferlist.normal_offer.d_title.4ff91a18MzPO0F
Octanediol_price = 6.5

# All in 2016$/kg
price = {'AA': AA_price,
         'TiO2': TiO2_price,
         'IBA': IBA_price,
         'Acetoin': acetoin_price,
         'Corn stover': cornstover_price,
         'Sugarcane': sugarcane_price,
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
         'Decanol': Decanol_price,
         'TOA': TOA_price,
         'AQ336': AQ336_price,
         'Octanol': Octanol_price,
         'Butyl acetate': Butyl_acetate_price,
         'Hexanol': Hexanol_price,
         'Hexanediol': Hexanediol_price,
         'Octanediol': Octanediol_price}
    

def load_process_settings():
    # from biorefineries import cornstover as cs
    # cs.load_process_settings()
    # tmo.settings.set_thermo(chems)
    bst.CE = 541.7 # year 2016
    bst.PowerUtility.price = price['Electricity']
    
    _lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    _mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
    _hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    _mps.T = 233 + 273.15
    _hps.T = 266 + 273.15
    
    # _lps.heat_transfer_efficiency = 0.900
    # _lps.heat_transfer_efficiency = 0.950
    # # Do this OR bst.default() since both do the same thing:
    # _lps.T = 412.19
    # _lps.P = 344738.
    # 
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
    
load_process_settings()

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
    'CH4': 0.33, # NA NG from shale and conventional recovery
    'Enzyme': 2.24,
    'Lime': 1.29,
    'NaOH': 2.11,
    'H2SO4': 0.04344,
    'Ethanol': 1.44
    }
H3PO4_GWP_CF = 2.5426
KOH_GWP_CF = 2.299
GWP_CFs['DPHP'] = 174.2*(H3PO4_GWP_CF/97.944 + 2*KOH_GWP_CF/56.106)



# GWP_CFs['CaCO3'] = 10.30/1e3
GWP_CFs['Gypsum'] = -4.20/1e3


# from ecoinvent 3.6 IPCC 2013:
# GWP_CFs['CalciumDihydroxide'] = 1.2105 * 56.0774 / 74.093 # /kg-quicklime converted to kg-slaked_lime assuming CF of 0 for water

# from GREET 2020 GHG-100:
GWP_CFs['CalciumDihydroxide'] = 1.29 * 56.0774 / 74.093 # /kg-quicklime converted to kg-slaked_lime assuming CF of 0 for water

# from ecoinvent 3.7 IPCC 2013
GWP_CFs['Hexanol'] = 0.87409 # currently set to CF of 1,1-dimethylcyclopentane to generic market for solvent, organic (per kg)

GWP_CFs['TiO2'] = 7.8029 # ecoinvent 3.7.1, market for titanium dioxide [RoW] - IPCC 2013 GWP100a

GWP_CF_array = chems.kwarray(GWP_CFs)


# In kg CO2-eq/kg of material
# print(tmo.settings.get_thermo())
# import pdb
# pdb.set_trace()
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')


# GWP_CFs['FGHTP Corn stover'] = 68.82/1000. # Wendt et al. 2018: Techno-Economic Assessment of a Chopped Feedstock Logistics Supply Chain for Corn Stover

GWP_CFs['FGHTP Corn stover'] = 0.10945 # see Feedstock_impacts_YL
GWP_CFs['FGHTP Sugarcane'] = 0.10945 # placeholder

# GWP_CFs['Corn stover'] = 44.70/1e3 * 0.8
# GWP_CFs['Switchgrass'] = 87.81/1e3 * 0.8
# GWP_CFs['Miscanthus'] = 78.28/1e3 * 0.8
# In kg CO2-eq/kWh
GWP_CFs['Electricity'] = 0.48
# GWP_CFs['Electricity'] = 0.
# From corn stover
GWP_CFs['LacticAcid_GREET'] = 1.80
# From ref [7], lactic acid production, RoW, TRACI global warming
GWP_CFs['LacticAcid_fossil'] = 4.1787

GWP_CFs['Sugarcane'] = 0.12158 # ecoinvent 3.6, sugarcane production, RoW, IPCC 2013 GWP-100a

CFs['GWP_CFs'] = GWP_CFs
CFs['GWP_CF_stream'] = GWP_CF_stream
# GWP_CFs['']
# =============================================================================
# Fossil energy consumption (FEC), in MJ/kg of material
# =============================================================================

FEC_CFs = {
    'NH4OH': 42. * chems.NH3.MW/chems.NH4OH.MW,
    'CSL': 12.,
    'CH4': 50., # NA NG from shale and conventional recovery
    'Enzyme': 26.,
    'Lime': 4.896,
    'NaOH': 29.,
    'H2SO4': 0.56898,
    'Ethanol': 16.
    }
H3PO4_FEC_CF = 39.542
KOH_FEC_CF = 30.421
FEC_CFs['DPHP'] = 174.2*(H3PO4_FEC_CF/97.944 + 2*KOH_FEC_CF/56.106)



# FEC_CFs['CaCO3'] = 133.19/1e3
FEC_CFs['Gypsum'] = -44.19/1e3


# from ecoinvent 3.6 IPCC 2013:
FEC_CFs['CalciumDihydroxide'] = 5.2339  * 56.0774 / 74.093 # /kg-quicklime converted to kg-slaked_lime assuming CF of 0 for water
FEC_CFs['Hexanol'] = 64.652 # currently set to CF of 1,1-dimethylcyclopentane to generic market for solvent, organic (per kg)

FEC_CFs['TiO2'] = 82.361 # ecoinvent 3.7.1, market for titanium dioxide [RoW] - CED fossil

FEC_CF_array = chems.kwarray(FEC_CFs)
# In MJ/kg of material
FEC_CF_stream = tmo.Stream('FEC_CF_stream', FEC_CF_array, units='kg/hr')

# FEC_CFs['FGHTP Corn stover'] = 767.3/1000. # Wendt et al. 2018: Techno-Economic Assessment of a Chopped Feedstock Logistics Supply Chain for Corn Stover


FEC_CFs['FGHTP Corn stover'] = 1.68000 # see Feedstock_impacts_YL
FEC_CFs['FGHTP Sugarcane'] = 1.68000 # placeholder


CFs['FEC_CFs'] = FEC_CFs
CFs['FEC_CF_stream'] = FEC_CF_stream
# In MJ/kWh
FEC_CFs['Electricity'] = 5.926
# FEC_CFs['Electricity'] = 0.
# From corn stover
FEC_CFs['LacticAcid'] = 29.
# From ref [7], lactic acid production, RoW, cumulative energy demand, fossil
FEC_CFs['LacticAcid_fossil'] = 79.524

FEC_CFs['Sugarcane'] = 0.37338 # ecoinvent 3.6, sugarcane production, RoW, IPCC 2013 GWP-100a




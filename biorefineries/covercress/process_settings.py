#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 26 09:34:40 2026

@author: princyk2
"""
"""
#[1] ecoinvent 3.6 https://www.ecoinvent.org/home.html (accessed Aug 26, 2020).
[2] :Life cycle energy and environmental impacts of hydroprocessed renewable
jet fuel production from pennycress
"""
import biosteam as bst
import thermosteam as tmo
from biorefineries.covercress.chemicals_cc import create_covercress_chemicals
chemicals=create_covercress_chemicals()

bst.CE = 541.7 # year 2016
_kg_per_ton = 907.18474
_lb_per_kg = 2.20462
_liter_per_gallon = 3.78541
_ft3_per_m3 = 35.3147
_cents_dollars=0.01

_corn_bushel_to_kg = 25.402 # https://www.ers.usda.gov/webdocs/publications/41880/33132_ah697_002.pdf
#PPI 2014 ≈ WPU06_2014 × (PPI2016 / WPU06_2016)
#WP06_2019 avg(jan-dec)-avg(293+288)-291
#PPI for 2016 is 102.5
# WPU06_2016-260
#WPU06_2014 from https://fred.stlouisfed.org/series/WPU06   
_chemical_2011to2016 = 102.5 / 91.7
_chemical_2013to2016 = 102.5 / 101.3
_chemical_2014to2016 = 102.5 / 105.3
_chemical_2017to2016 = 102.5 / 106.9
_chemical_2020to2016 = 102.5 / 113.8 # average of Jan and Feb
_chemical_2022to2016 = 102.5 / 145.3
_chemical_2019to2016 = 102.5 / 114.72 # refer eq above
_chemical_index_2025= 126.6
_chemical_index_2025to2016= 102.5/126.6
_GDP_2007_to_2016 = 1.114 / 0.961
_GDP_2008_to_2016 = 1.114 / 0.990
_GDP_2008_to_2010 = 1.012 / 0.990
_GDP_2007_to_2010 = 1.012 / 0.961

chem_index = { # Dictionary of chemical indices
                    2010: 82.2,
                    2011: 79.5,
                    2012: 83.7,
                    2013: 87.9,
                    2014: 91.3,
                    2015: 93.1,
                    2016: 88.8,
                    2017: 92.7,
                    2018: 93.3,
                    2019: 97.0, # the TEA year for the 2023 TAL production study
                    2020: 100.2,
                    2021: 112.0,
                    2022: 125.829,
                    }
#%% Feedstocks

#TODO need to update
covercress_price = 0.20 * _chemical_2019to2016 #refer[4] table2 assumed penycress price
#crude oil price - 59cts/lb for soybean oil chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.ams.usda.gov/mnreports/lswagenergy.pdf
#crudeoil from canola-https://dir.tridge.com/prices/canola-oil
#highest price 3.85 $/Kg and 2.09 $/Kg-avg=2.99 for 2022
_crudeoil_price = 2.99*_chemical_2022to2016 #as per 2016

ash_disposal_price = -1.41e6 / (4279*7880)

gypsum_price = 0

lime_price = 0.1189 * _lb_per_kg

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_price = 466833 / 5 / (24*365*0.96)
enzyme_price=4.5*102.5/126.6  #https://www.chemanalyst.com/Pricing-data/protease-enzyme-2277

# 4.70 is the average natural gas price in 2016$/Mcf based on AEO from EIA,
# which is $0.231/kg at 273.15 K or $0.253/kg at 298.15 K using BioSTEAM,   
# similar to the 4.7/1000/22*1000 = $0.214/kg at 273.15 K using 22 g/ft3 from H2 Tools
# Using the same conversion, lower and upper bounds should be 
# $3.68/Mcf and $5.65/Mcf, or $0.198/kg and $0.304/kg
CH4_V = chemicals.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chemicals.CH4.MW
natural_gas_price = 4.70/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW)
Makeup_water_price=0.0002 * _lb_per_kg
NaOH_price=0.85
NaCl_price=0.029 #https://businessanalytiq.com/procurementanalytics/index/sodium-chloride-price-index/
hexane_price=2.06*_chemical_2019to2016 #refer [4] *0.655 #https://catcost.chemcatbio.org/materials-library
methanol_price =0.39*_chemical_2022to2016 #https://www.alibaba.com/showroom/price-industrial-methanol.html?src=sem_bing&field=UG&from=sem_bing&cmpgn=487347983&adgrp=1230354838541717&tgt=kwd-2328696910231351:loc-190&KwdID=2328696910231351&mtchtyp=b&bdmtchtyp+=bb&ntwrk=o&device=c&creative=76897320341689&p1=default&p2=default&p3=default&Query=methanol%20price&msclkid=90348bc31aa710d432a4e263cc5a3a3b
ethanol_price =1*_chemical_index_2025to2016 #https://www.alibaba.com/product-detail/Ethanol-for-Industrial-Use-99-96_11000015420804.html?spm=a2700.galleryofferlist.normal_offer.d_image.145113a0tOayNM&priceId=c297325c6f1f44a0815ea1f5057ee641
crudeoil_price= 0.995*_chemical_index_2025to2016 #https://www.imarcgroup.com/crude-soybean-oil-pricing-report?utm
protein_price=  1.77*_chemical_index_2025to2016# 1.77 https://www.chemanalyst.com/Pricing-data/soy-protein-isolate-1578
citricacid_price=1.07*_chemical_2019to2016 #*-----refer[4]......https://www.alibaba.com/product-detail/Industrial-Grade-Anhydrous-Citric-Acid-Wholesale_1601574006767.html?spm=a2700.prosearch.normal_offer.d_title.568667afYAhbyD&priceId=ef4ba8fed82641b8abc6a270972fe145
# SAF_price=2.305*_chemical_index_2025to2016 
catalyst_price=337.03*_chemical_2019to2016 #[refer[4]]
Hydrogen_price=2.96*_chemical_2019to2016 #[refer[4]]
SAF_price= 0.6*_chemical_2019to2016 #[refer[4]] USD/L-0.48 density -(0.76-0.84) avg 0.8 so 0.48/0.8
Propane_price= 0.23*_chemical_2019to2016 #[refer[4]C3-C4 #0.96 #refer [34]  or 
Naptha_price=  0.87*_chemical_2019to2016 #[refer[4]C5-C7 # 0.96 #[34]avg refer table2  or
Green_diesel_price= 0.975*_chemical_2019to2016 #[refer[4] C17-C22  #1.4 #refer[34] or density of GD is (0.77-0.83) avg is 0.8 g/ml and for 0.78$/L to $/kg  0.78/0.8
meal_price=0.345*102.5/126.6 #later update to 2025 , this one is for 2026
# All in 2016$/kg
price = {
          'Feedstock': covercress_price,
         # 'Sulfuric acid': 0.0430 * _lb_per_kg,	
         # 0.1900 is for NH3	
         # 'AmmoniumHydroxide': 0.1900 * _lb_per_kg * 17.031/35.046,	
         'Methanol': methanol_price,
          'Ethanol':ethanol_price,
         'Caustics': 0.2384 * _lb_per_kg * 0.5, # 50 wt% NaOH/water mixture	
          'Hexane':hexane_price,
          'NaOH':NaOH_price,
          'NaCl': NaCl_price,
         'Lime': lime_price,
         'Makeup water': Makeup_water_price,	
         # Cost of ash is negative because it's a product stream	
         'Ash disposal': ash_disposal_price,
         'Electricity': 0.070, # AEO from EIA, 2010-2019 average (0.067-0.074 range)	
         'Baghouse bag': baghouse_bag_price,	
         'Natural gas': natural_gas_price,
         'crude oil':_crudeoil_price,
         'Gypsum': gypsum_price,
         'crude_oil': crudeoil_price,
         'protein_isolate':protein_price,
         'citric_acid':citricacid_price,
         'SAF':SAF_price,
         'Pt/Al2O3':catalyst_price,
         'Hydrogen':Hydrogen_price,
         'Propane':Propane_price,
         'Naptha':Naptha_price,
         'Green_diesel':Green_diesel_price,
         'enzyme':enzyme_price,
         'meal':meal_price,
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

for i in (_lps, _mps, _hps, _cooling, _chilled):
    i.heat_transfer_price = i.regeneration_price = 0
# heat transfer price na dregenraion price for _lps, _mps, _hps, _cooling, _chilled=0 because as boiler produces the steam required fro entire plant and as natural gas is 


# %%

# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA), all from ref [1] if not noted otherwise
# =============================================================================

CFs = {}

# =============================================================================
# 100-year global warming potential (GWP) in kg CO2-eq/kg
# =============================================================================
GWP_CFs = {
    'CH4': 0.40, # NA NG from shale and conventional recovery
    # 'H2SO4': 44.47/1e3,   
    'Lime': 1.29 * 56.0774/74.093, # CaO to Ca(OH)2
    'NaOH': 2.11,
    # 'Ethanol': 
    # 'NH4OH': 2.64 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,   
    # 'MEA': 3.4062, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 1.3598, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 0.87104, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
     'Hexane':  0.55002, # [kg*CO2*eq / kg] Ecoinvent 3.6 Cut off; IPCC 2013; Global market for hexane
     'Methanol':0.66082,	#ecoinvent 3.8, market for methanol,from oxalic process settings
     'NaOCH3':1, #TODO update 
     'Hydrogen': 1.63, #ecoinvent https://ecoquery.ecoinvent.org/3.8/cutoff/dataset/20555/impact_assessment
     'Pt/Al2O3':6.71e+4,	#kg CO2-Eq #TODO update assumed as only platinum 
     'Citric_acid': 5.9048,  # kg CO2-eq/kg, ecoinvent 3.8
     'NaCl': 2.11, #TODO update
     'Enzyme':1 #TODO update place holder check for protease enzyme 
    }


GWP_CF_array = chemicals.kwarray(GWP_CFs)
#GWP_CF_array will do the (some_mass_flow_array * GWP_CF_array).sum() (Kg,CO2/kg*[kg CO₂-eq/kg)
# will give In kg CO2-eq/kg of material
GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')
#make this array as a single stream in Kg/hr to calculate total GWP
# CFs['GWP_CF_stream'] = GWP_CF_stream

GWP_CFs['Electricity'] = 0.4490 # kg CO2-eq/kWh GREET 2022 US Mix  # assume production==consumption, both in kg CO2-eq/kWh
GWP_CFs['Feedstock'] = 0.115  # refer excel cal in the [ref]
CFs['GWP_100'] = GWP_CFs


# =============================================================================
# Fossil energy consumption (FEC), in MJ/kg of material
# =============================================================================


FEC_CFs = {
    'CH4': 50, # NA NG from shale and conventional recovery
    'Methanol': 34.713, #Oxalic and HP process_settings.py (ecoinvent 3.8, GLO
    # 'H2SO4': 568.98/1e3,
    'Lime': 4.896 * 56.0774/74.093, # CaO to Ca(OH)2
    'NaOH': 29,
    # 'NH4OH': 42 * 0.4860, # multiplied by chemicals.NH3.MW/chemicals.NH4OH.MW,
    # 'MEA': 67.898, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'H3PO4': 16.538, # ecoinvent 3.8 purification of wet-process phosphoric acid to industrial grade, product in 85% solution state, RoW # cradle-to-gate
    'CO2': 7.4243, # ecoinvent 3.8 carbon dioxide production, liquid, RoW
    'NaOCH3':1,  #TODO update 
    'Hexane':44,   #Refer[2] 
    'Hydrogen':183.2, #Refer[2]
    # 'Platinum on Aluminum (Pt/ɣ-Al2O3) ': 5.4 #MJ/kg  #Trefer [2]
    #'Citric_acid'
    # 'NaCl'
    }

FEC_CF_array = chemicals.kwarray(FEC_CFs)
# In MJ/kg of material
FEC_CF_stream = tmo.Stream('FEC_CF_stream', FEC_CF_array, units='kg/hr')

# CFs['FEC_CF_stream'] = FEC_CF_stream

FEC_CFs['Electricity'] = 5.724 # MJ/kWh # GREET 2022 US Mix #assume production==consumption, both in MJ/kWh
#TODO update for covercress
FEC_CFs['Feedstock'] = 1
CFs['FEC'] = FEC_CFs





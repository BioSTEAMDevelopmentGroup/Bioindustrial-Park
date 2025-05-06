#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 10:45:27 2024

@author: wenjun
"""

import numpy as np
import biosteam as bst
import qsdsan as qs
import biorefineries.SAF._chemicals as chems


__all__=('load_preferences_and_process_settings',
         'add_utility_agent',
         'price',
         'GWP_CFs')


# %%
 

# =============================================================================
# System preferences and settings
# =============================================================================
def load_preferences_and_process_settings(T,flow_units,N,P_units,CE,
                                          indicator,electricity_price,
                                          electricity_EI):
    bst.preferences.T = T
    bst.preferences.flow = flow_units
    bst.preferences.N = N
    bst.preferences.P = P_units
    bst.preferences.composition = True
    bst.preferences.light_mode()
    bst.preferences.save()
    bst.settings.CEPCI = CE  
    bst.settings.define_impact_indicator(key=indicator, units='kg*CO2e')
    bst.settings.electricity_price = electricity_price
    bst.settings.set_electricity_CF(indicator,electricity_EI, basis='kWhr', units='kg*CO2e')
   
    lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
    hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    
    # hps.T = 266 + 273.15
    # hps.P = 44e5

    # Side steam in CHP not a heat utility, thus will cause problem in TEA utility
    # cost calculation if price not set to 0 here, costs for regeneration of heating
    # and cooling utilities will be considered as CAPEX and OPEX of CHP and CT, respectively
    for i in (lps, mps, hps):
        i.heat_transfer_price = i.regeneration_price = 0
    
    bst.settings.get_agent('cooling_water').regeneration_price = 0
    bst.settings.get_agent('chilled_water').heat_transfer_price = 0
    bst.settings.get_agent('chilled_brine').heat_transfer_price = 0
    bst.settings.get_agent('propane').heat_transfer_price = 0
    



# =============================================================================
# Add utility agent: heating + cooling(for cryogenic distillation)
# =============================================================================
def add_utility_agent():
    # Add a heating agent
    DPO_chem = qs.Chemical('DPO_chem', search_ID='101-84-8')
    BIP_chem = qs.Chemical('BIP_chem', search_ID='92-52-4')
    DPO = qs.Component.from_chemical('DPO', chemical=DPO_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    BIP = qs.Component.from_chemical('BIP', chemical=BIP_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    HTF_thermo = bst.Thermo((DPO, BIP,))
    HTF = bst.UtilityAgent('HTF', DPO=0.735, BIP=0.265, T=673.15, P=10.6*101325, phase='g',
                           # 400 C (673.15 K) and 138 psig (951477 pa) are max temp and pressure for HTF
                           thermo=HTF_thermo,
                           # T_limit = 495 F (530.372 K) is the highest temp that vapor can exist
                           regeneration_price=1) 
    # from Lang, C.; Lee, B. Heat Transfer Fluid Life Time Analysis of Diphenyl Oxide/Biphenyl Grades for Concentrated Solar Power Plants.
                          
    # Temperature and pressure: https://www.dow.com/content/dam/dcc/documents/\
    # en-us/app-tech-guide/176/176-01334-01-dowtherm-heat-transfer-fluids-\
    # engineering-manual.pdf?iframe=true (accessed on 11-16-2022)
    bst.HeatUtility.heating_agents.append(HTF)

    # Add a cooling agent
    Decamethyltetrasiloxane_chem = qs.Chemical('Decamethyltetrasiloxane_chem',search_ID='141-62-8')
    Octamethyltrisiloxane_chem = qs.Chemical('Octamethyltrisiloxane',search_ID='107-51-7')
    DEC = qs.Component.from_chemical('DEC', chemical=Decamethyltetrasiloxane_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    OCT = qs.Component.from_chemical('OCT', chemical=Octamethyltrisiloxane_chem, particle_size='Soluble',
                                     degradability='Slowly', organic=True)
    LTF_thermo = bst.Thermo((DEC,OCT,))
    LTF = bst.UtilityAgent('LTF', DEC=0.4, OCT=0.6, T=173.15, P=5.2*101325, phase='l',
                           thermo = LTF_thermo,
                           heat_transfer_price=0.00001317) # Based on propane heat transfer price
    bst.HeatUtility.cooling_agents.append(LTF)

# =============================================================================
# Prices for techno-economic analysis (TEA)
# Adjust price to 2023 to match up with jet fuel tax credit for 2023-2024
# =============================================================================
_lb2kg = 0.453592
_ft3_per_m3 = 35.3147
_kgal_to_L = 3.785412 * 1000
_gal_to_L = 3.785412
gasoline_density = 0.7489 # in kg/L
gasoline_price_conversion_index = 1 / _gal_to_L / gasoline_density # from $/gal to $/kg
diesel_density = 0.838 # in kg/L
diesel_price_conversion_index = 1 / _gal_to_L / diesel_density # from $/gal to $/kg

# Producer price index (PPI) by Commodity: Chemicals and Allied Products, 
# from https://fred.stlouisfed.org/series/WPU06
PPI_2016 = 265.108 # Average 2016
PPI_2017 = 280.825 # Average 2017
PPI_2018 = 295.150 # Average 2018
PPI_2019 = 289.125 # Average 2019
PPI_2021 = 331.413 # Average 2021
PPI_2022 = 366.788 # Average 2022
PPI_2023 = 346.236 # Aug 2023

_chemical_2016to2023 = PPI_2023/PPI_2016
_chemical_2021to2023 = PPI_2023/PPI_2021
_chemical_2019to2023 = PPI_2023/PPI_2019

# Producer Price Index by Industry: Utilities,
# from https://fred.stlouisfed.org/series/PCU221221
PPI_utility_2016 = 134.842 # Average 2016
PPI_utility_2017 = 141.042 # Average 2017
PPI_utility_2018 = 145.808 # Average 2018
PPI_utility_2019 = 146.767 # Average 2019
PPI_utility_2020 = 143.283 # Average 2020
PPI_utility_2021 = 174.764 # Average 2021
PPI_utility_2022 = 201.339 # Average 2022
PPI_utility_2016_2022 = [PPI_utility_2016,
                         PPI_utility_2017,
                         PPI_utility_2018,
                         PPI_utility_2019,
                         PPI_utility_2020,
                         PPI_utility_2021,
                         PPI_utility_2022]

PPI_utility_2023 = 187.389 # Average 2023
_utility_2016_2022_period_to2023 = PPI_utility_2023/np.mean(PPI_utility_2016_2022)
_utility_2021_to2023 = PPI_utility_2023/PPI_utility_2021

_GDP_2007to2023 = 122.254/86.346 # https://fred.stlouisfed.org/series/GDPDEF

energycane_price = 0.035

H3PO4_price = 1.25 # https://catcost.chemcatbio.org/materials-library accessed Jan 7 2024

flocculant_price = 1 # https://jucheng88888.en.made-in-china.com/product/yTBrhovAZnWK/China-Flocculation-and-Sedimentation-Anionic-Polyacrylamide.html accessed Jan 8 2024

lime_price = 0.1 # https://catcost.chemcatbio.org/materials-library accessed Jan 7 2024

enzyme_price = 6.16 * _chemical_2016to2023 # From lactic in $2016

CSL_price = 0.0339 / _lb2kg *_chemical_2016to2023# From lactic in $2016

DAP_price = 0.1645 /_lb2kg *_chemical_2016to2023 # From ethanol_adiapic

NaOH_price = 0.86 * _chemical_2021to2023 # https://catcost.chemcatbio.org/materials-library accessed Mar 19

Syndol_catalyst_price = 20.48 * _chemical_2021to2023 # From Advanced fuels from ethanol – a superstructure optimization approach(2021)

Ni_loaded_aluminosilicate_catalyst_price = 12.3 # Bulk price accessed in Jan 7 2024 from https://www.sigmaaldrich.com/US/en/product/aldrich/208779 accessed Jan 7 2024

Aluminosilicate_catalyst_price = 1.34 # $1340/MT for Sodium Aluminosilicate https://www.chemanalyst.com/Pricing-data/aluminosilicate-1518 accessed Feb 2 2024 

Como_catalyst_price = 5 # Cobalt molybdenum catalyst, for 50 kg above, https://biz.alibaba.com/contract/poBuy.htm?productData=[{%22productId%22:1600703116039,%22quantity%22:50}]&channelType=PRODUCT_DETAIL accessed Jan 7 2024

H2_price = 1.07 * _chemical_2021to2023 # From 'Advanced fuels from ethanol – a superstructure optimization approach (2021)' 

ash_disposal_price = -1.41e6/(4279*7880) * _chemical_2016to2023

boiler_chems_price = 809/55/9.1/_lb2kg # 55 gal cost $809 https://www.chemworld.com/ProductDetails.asp?ProductCode=CHEMWORLD%2D1394
                                       # density is 9.1 lb/gal, from pdf.on website accessed Jan 10 2024

cooling_tower_price = 519/55/9.7/_lb2kg # 55 gal cost $519 https://www.coolingtowerchemicals.com/ProductDetails.asp?ProductCode=CTC1334NM 
                                        # density is 9.70 lb/gal, from https://www.chemworld.com/v/vspfiles/assets/images/sds-ChemWorld1334.pdf accessed Jan 10 2024

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_price = 466833 / 5 / (24*365*0.96) *_GDP_2007to2023

natural_gas_price_2016_to_2022 = [3.51,4.08,4.19,3.90,3.32,5.44,7.66] # https://www.eia.gov/dnav/ng/hist/n3035us3a.htm, $/Mcf(thousand cubic feet)

electricity_price_2016_to_2022 = [6.76,6.88,6.92,6.81,6.67,7.18,8.45]  # https://www.eia.gov/electricity/data/browser/#/topic/7?agg=2,0,1&geo=g&freq=M cents/Kwh

water_price_2021 = 3.86 # $2021/kGal https://www.osti.gov/servlets/purl/1975260

CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chems.CH4.MW
natural_gas_price = np.mean(natural_gas_price_2016_to_2022)/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW) * _utility_2016_2022_period_to2023
     
electricity_price = np.mean(electricity_price_2016_to_2022) * _utility_2016_2022_period_to2023/100 # cents to $

water_price = water_price_2021 / _kgal_to_L * _utility_2021_to2023

# Co-product credits
gasoline_price_2021_to_2023 = [3.1,4.059,3.635] # $/gallon https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_nus_a.htm

diesel_price_2021_to_2023 = [3.287,4.989,4.214] # $/gallon https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_nus_a.htm

jet_fuel_price_2023 = 2.437 # $/kg https://www.iata.org/en/iata-repository/publications/economic-reports/sustainable-aviation-fuel-output-increases-but-volumes-still-low/

gasoline_price = np.mean(gasoline_price_2021_to_2023) * gasoline_price_conversion_index # s/gallon to $/kg

diesel_price = np.mean(diesel_price_2021_to_2023) * diesel_price_conversion_index # s/gallon to $/kg

jet_fuel_price = jet_fuel_price_2023

MEA_price = 1.44 # from https://businessanalytiq.com/procurementanalytics/index/monoethanolamine-price-index/

compressed_CO2_price = -1*(0.022 + 0.03) # Transportation and storage price from https://netl.doe.gov/projects/files/ComparativeAnalysisofTransportandStorageOptionsfromaCO2SourcePerspective_042018.pdf

# All in 2023 $/kg
price = {
    'feedstock': 0.0875*0.8,
    'water': water_price,
    'H3PO4': H3PO4_price * _chemical_2021to2023,
    'flocculant': flocculant_price, 
    'lime': lime_price * _chemical_2021to2023,  
    # 'enzyme': enzyme_price,
    'CSL': CSL_price, 
    'DAP': DAP_price,
    'NaOH': NaOH_price * 0.5, # 50% NaOH
    'Syndol catalyst': Syndol_catalyst_price, 
    'Ni-loaded aluminosilicate catalyst': Ni_loaded_aluminosilicate_catalyst_price, 
    'Aluminosilicate catalyst': Aluminosilicate_catalyst_price,
    'Como catalyst': Como_catalyst_price, 
    'h2': H2_price,
    'ash disposal': ash_disposal_price,
    'boiler chems': boiler_chems_price,
    'cooling tower chems': cooling_tower_price,
    'baghouse bag': baghouse_bag_price,
    'natural gas': natural_gas_price,
    'electricity': electricity_price,
    'gasoline': gasoline_price,
    'jet fuel': jet_fuel_price,
    'diesel': diesel_price,
    'MEA': MEA_price * 0.3, # 30% MEA
    'compressed CO2': compressed_CO2_price,
    'oc': 0.4/3.78/845*1000,
    'ethanol': 0.41/0.789, # $0.41/L is from https://www.icao.int/environmental-protection/Pages/SAF_RULESOFTHUMB.aspx
        }



# %%

# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA)
# =============================================================================

# 100-year global warming potential (GWP) in kg CO2-eq/kg

GWP_CFs = {
    'switchgrass': -0.14315*0.8, # -154.61 gCO2 is dry miscanthus average CI of PADD2 from Xinxin's data
    'sulfuric_acid': 1.3527/1000/0.03, # GREET
    'ammonia': 3.1996, # ecoinvent 3.8 market for polyacrylamide-based anionic flocculants, GLO
    'denaturant': 0.8415,
    # 'ethanol':-1.28*0.915,  # ethanol is a co-product from https://ecoquery.ecoinvent.org/3.8/cutoff/dataset/5353/impact_assessment 91.5% is wt%
    'cellulase': 2.24*0.05, # 5% concentration
    'CSL': 1.56,
    'DAP': 1.6445, # ecoinvent 3.8 market for diammonium phosphate, RoW
    'NaOH': 2.01, # Greet
    'caustic':2.01,
    'lime': 1.29,
    'lime_boiler_1':1.29,
    'FGD_lime':1.29/2, # 50% concentration
    'boiler_chems':1.5568,
    'boiler_chems_1':1.5568, # # Based on the production of sulfite which is commonly used as an oxygen scavenger in the boiler
    'hydrogen': 1.5624, # ecoinvent 3.8 market for hydrogen, gaseous, GLO
    'natural_gas': 0.40, # NA NG from shale and conventional recovery
    'natural_gas_1': 0.40,
    'electricity': 0.48, # in kg CO2-eq/kWh https://www.eia.gov/tools/faqs/faq.php?id=74&t=11
    'ash':0.01, 
    'ash_1':0.01,# from Production of Sustainable Aviation Fuels in Petroleum Refineries: Evaluation of New Bio-Refinery Concepts
    'Syndol_catalyst':10+1, # production+disposal, assume prodution=10 based on similar chemicalss
    'first_catalyst':10+1, # production+disposal, assume prodution=10 based on similar chemicalss
    'second_catalyst':10+1, # production+disposal, assume prodution=10 based on similar chemicalss
    'Como_catalyst':9.52+1, # production+disposal, from Production of Sustainable Aviation Fuels in Petroleum Refineries: Evaluation of New Bio-Refinery Concepts
    'gasoline':-0.8415, # gasoline is a co-product; from GREET
    'diesel':-0.6535 # diesel is a co-product; from GREET
    }







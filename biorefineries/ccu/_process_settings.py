# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 14:27:18 2025

@author: IGB
"""

import numpy as np
import biosteam as bst
from EtOH._chemicals import create_MeOH_chemicals
import biorefineries.cellulosic.streams as cellulosic_stream

chems = create_MeOH_chemicals()
__all__=('load_preferences_and_process_settings',
         'price',
         'CFs')


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
# Prices for techno-economic analysis (TEA)
# Adjust price to 2023
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
PPI_2023 = 351.821 # Average 2023
PPI_2024 = 346.560 # Average 2024

_chemical_2016to2023 = PPI_2023/PPI_2016
_chemical_2021to2023 = PPI_2023/PPI_2021
_chemical_2019to2023 = PPI_2023/PPI_2019
_chemical_2024to2023 = PPI_2023/PPI_2024

# Producer Price Index by Industry: Utilities,
# from https://fred.stlouisfed.org/series/PCU221221
PPI_utility_2016 = 134.842 # Average 2016
PPI_utility_2017 = 141.042 # Average 2017
PPI_utility_2018 = 145.808 # Average 2018
PPI_utility_2019 = 146.767 # Average 2019
PPI_utility_2020 = 143.283 # Average 2020
PPI_utility_2021 = 174.764 # Average 2021
PPI_utility_2022 = 201.339 # Average 2022
PPI_utility_2023 = 187.203 # Average 2023
PPI_utility_2016_2023 = [PPI_utility_2016,
                         PPI_utility_2017,
                         PPI_utility_2018,
                         PPI_utility_2019,
                         PPI_utility_2020,
                         PPI_utility_2021,
                         PPI_utility_2022,
                         PPI_utility_2023]
                         
_utility_2016_2023_period_to2023 = PPI_utility_2023/np.mean(PPI_utility_2016_2023)
_utility_2021_to2023 = PPI_utility_2023/PPI_utility_2021

_GDP_2007to2023 = 122.254/86.346 # https://fred.stlouisfed.org/series/GDPDEF

cornstover_price = 0.068682 * _chemical_2016to2023 # from Sarang's 3HP

sulfuric_acid_price = 0.1375 * _chemical_2021to2023 # from catcost; average of US except east coast

water_price_2021 = 3.86 # $2021/kGal https://www.osti.gov/servlets/purl/1975260

water_price = water_price_2021 / _kgal_to_L * _utility_2021_to2023

ammonia_price = 0.825 * _chemical_2021to2023 # from catcost; average

cellulase_price = 6.16 * _chemical_2016to2023 # From lactic in $2016; enzyme

CSL_price = 0.0747 * _chemical_2016to2023 # from Sarang's 3HP

DAP_price = cellulosic_stream.DAP['price'] * _chemical_2019to2023

NaOH_price = 0.86 * _chemical_2021to2023 # https://catcost.chemcatbio.org/materials-library accessed Mar 19

MEA_price = 1.44 # from https://businessanalytiq.com/procurementanalytics/index/monoethanolamine-price-index/

catalyst_MeOH_price = 32 * _chemical_2024to2023 # https://pubs.acs.org/doi/full/10.1021/acssuschemeng.4c08820 and https://www.sciencedirect.com/science/article/pii/S1385894723066056?pes=vor&utm_source=acs&getft_integrator=acs

denaturant_price = cellulosic_stream.denaturant['price'] * _chemical_2019to2023

lime_price = 0.13514 # https://www.usgs.gov/media/files/lime-historical-statistics-data-series-140

boiler_chems_price = 809/55/9.1/_lb2kg # 55 gal cost $809 https://www.chemworld.com/ProductDetails.asp?ProductCode=CHEMWORLD%2D1394
                                       # density is 9.1 lb/gal, from pdf.on website accessed Jan 10 2024

CH4_V = chems.CH4.V(298.15, 101325) # molar volume in m3/mol
CH4_MW = chems.CH4.MW
natural_gas_price_2016_to_2023 = [3.51,4.08,4.19,3.90,3.32,5.44,7.69,4.53] # https://www.eia.gov/dnav/ng/hist/n3035us3a.htm $/Mcf(thousand cubic feet)
natural_gas_price = np.mean(natural_gas_price_2016_to_2023)/1e3*_ft3_per_m3*CH4_V * (1e3/CH4_MW) * _utility_2016_2023_period_to2023

electricity_price_2016_to_2023 = [6.76,6.88,6.92,6.81,6.67,7.18,8.32,8.04]  # https://www.eia.gov/electricity/data/browser/#/topic/7?agg=2,0,1&geo=g&freq=M cents/Kwh
electricity_price = np.mean(electricity_price_2016_to_2023) * _utility_2016_2023_period_to2023/100 # cents to $
electricity_price_low = 6.32 * _utility_2016_2023_period_to2023/100
electricity_price_high = 9.38 * _utility_2016_2023_period_to2023/100


electricity_wind_price_2016_to_2022 = [0.0276,0.0213,0.0195,0.0269,0.0295,0.025,0.0255,] # https://view.officeapps.live.com/op/view.aspx?src=https%3A%2F%2Femp.lbl.gov%2Fsites%2Fdefault%2Ffiles%2F2024-08%2FLand-Based%2520Wind%2520Market%2520Report_2024%2520Edition_Data_File.xlsx&wdOrigin=BROWSELINK
electricity_wind_price = np.mean(electricity_wind_price_2016_to_2022)

ash_disposal_price = -1.41e6/(4279*7880) * _chemical_2016to2023

ethanol_price = 0.866/0.747561512219763 # average value

H3PO4_price = 1.25 # https://catcost.chemcatbio.org/materials-library accessed Jan 7 2024

flocculant_price = 1 # https://jucheng88888.en.made-in-china.com/product/yTBrhovAZnWK/China-Flocculation-and-Sedimentation-Anionic-Polyacrylamide.html accessed Jan 8 2024

H2_price = 1.07 * _chemical_2021to2023 # From 'Advanced fuels from ethanol – a superstructure optimization approach (2021)' 

cooling_tower_price = 519/55/9.7/_lb2kg # 55 gal cost $519 https://www.coolingtowerchemicals.com/ProductDetails.asp?ProductCode=CTC1334NM 
                                        # density is 9.70 lb/gal, from https://www.chemworld.com/v/vspfiles/assets/images/sds-ChemWorld1334.pdf accessed Jan 10 2024

# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# The original cost is $466,183 every 5 years, converted to per hour assuming 96% uptime
baghouse_bag_price = 466833 / 5 / (24*365*0.96) *_GDP_2007to2023


# Co-product credits
oxygen_price = 0.23 # bulk oxygen price; https://www.imarcgroup.com/oxygen-pricing-report#:~:text=Oxygen%20Prices%20Q4%202023&text=The%20United%20States%20oxygen%20price,2023%20was%20348%20USD/MT.

methanol_price = 0.35 # https://www.methanol.org/methanol-price-supply-demand/
# All in 2023 $/kg
price = {
    'feedstock': cornstover_price,
    'sulfuric_acid': sulfuric_acid_price,
    'makeup_process_water': water_price, # PWC makeup water
    'makeup_RO_water': water_price/27*56, # PWC makeup RO water
    'ammonia': ammonia_price,
    'cellulase': cellulase_price * 0.05, # Cellulase concentration is 5%
    'CSL': CSL_price, 
    'DAP': DAP_price,
    'caustic': NaOH_price * 0.5, # 50% NaOH
    'makeup_MEA': MEA_price * 0.3, # 30% MEA
    'denaturant': denaturant_price,
    'catalyst_MeOH': catalyst_MeOH_price,
    'cooling_tower_chemicals': cooling_tower_price,
    'natural_gas': natural_gas_price,
    'FGD_lime': lime_price,
    'boiler_chemicals': boiler_chems_price,
    'ethanol': ethanol_price,
    # 'H3PO4': H3PO4_price * _chemical_2021to2023,
    # 'flocculant': flocculant_price,
    'h2': H2_price,
    'ash': ash_disposal_price,
    'baghouse_bag': baghouse_bag_price,
    'oxygen': oxygen_price,
    'methanol': methanol_price,
    # 'electricity': electricity_price,
    'Electricity': electricity_wind_price
        }



# %%

# =============================================================================
# Characterization factors (CFs) for life cycle analysis (LCA)
# =============================================================================
CFs = {}

# 100-year global warming potential (GWP) in kg CO2-eq/kg
GWP_CFs = {
    'cornstover': 0.10945, # Table S4 of the SI of Bhagwat et al. 2021
    'sugarcane': 0.044535, # GREET 2023, Sugarcane Production for Brazil Plant
    'corn': 0.2610, # GREET 2023, Corn Production for Biofuel Refinery
    'H2SO4': 0.159, # ecoinvent 3.8 market for sulfuric acid (RoW); https://ecoquery.ecoinvent.org/3.8/cutoff/dataset/18422/impact_assessment
    'NH4OH': 2.53, # ecoinvent 3.8 market for liquid ammonia (North America); https://ecoquery.ecoinvent.org/3.8/cutoff/dataset/22603/impact_assessment
    'Cellulase': 2.24,
    'CSL': 1.55,
    'DAP': 1.6445, # ecoinvent 3.8 market for diammonium phosphate, RoW
    'MEA': 3.4062, # ecoinvent 3.8 ethanolamine production, RoW [monoethanolamine]
    'NaOH': 2.11,
    'Lime': 1.29,
    # 'makeup_RO_water': 0.00587, # ecoinvent 3.8 market for water, ultrapure (RoW)
    # 'makeup_process_water': 0.00102, # ecoinvent 3.8 market for tap water (RoW, industrial or household), RoW
    'CaO': 10+1, # catalyst_MeOH; production+disposal, assume prodution=10, disposal=1, based on similar chemicalss
    'CH4': 0.40, # NA NG from shale and conventional recovery
    # 'electricity': 0.435, # in kg CO2-eq/kWh https://www.eia.gov/tools/faqs/faq.php?id=74&t=11
    'Electricity': 0.011, # wind; https://www.energy.gov/eere/wind/articles/how-wind-can-help-us-breathe-easier#:~:text=In%20general%2C%20lifecycle%20greenhouse%20gas,2%2FkWh%20for%20natural%20gas.
    'H2': 3.53,
    'MeOH': 0.58, # ecoinvent 3.8 for methanol production from natural gas
    'oxygen': 0.1726, # GREET 2022; at gate; production of oxygen
    # 'H3PO4': 0.86829, # ecoinvent 3.8 market for phosphoric acid, RoW
    # 'polymer': 3.1996, # ecoinvent 3.8 market for polyacrylamide-based anionic flocculants, GLO
    }
CFs['GWP_100'] = GWP_CFs
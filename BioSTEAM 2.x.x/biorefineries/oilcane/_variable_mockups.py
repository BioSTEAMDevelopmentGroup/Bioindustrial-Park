# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:44:17 2021

@author: yrc2
"""
from biosteam import MockVariable

(set_cane_oil_content, set_relative_sorghum_oil_content, set_bagasse_oil_retention, 
 set_bagasse_oil_extraction_efficiency, 
 set_plant_capacity, set_ethanol_price,
 set_biodiesel_price, set_natural_gas_price, 
 set_electricity_price, set_operating_days, 
 set_IRR, set_crude_glycerol_price, set_pure_glycerol_price,
 set_saccharification_reaction_time, set_cellulase_price, 
 set_cellulase_loading, set_reactor_base_cost,
 set_cane_glucose_yield, set_cane_xylose_yield, 
 set_sorghum_glucose_yield, set_sorghum_xylose_yield, 
 set_glucose_to_ethanol_yield, set_xylose_to_ethanol_yield,
 set_cofermentation_titer, set_cofermentation_productivity,
 set_cane_PL_content, set_sorghum_PL_content, set_cane_FFA_content,
 set_sorghum_FFA_content, set_TAG_to_FFA_conversion, set_oilcane_GWP,
 set_methanol_GWP, set_pure_glycerine_GWP, set_cellulase_GWP, set_natural_gas_GWP
 ) = all_parameter_mockups = (
    MockVariable('Oil retention', '%', 'Stream-sugarcane'),
    MockVariable('Bagasse oil extraction efficiency', '%', 'Stream-sugarcane'),
    MockVariable('Capacity', 'ton/hr', 'Stream-sugarcane'),
    MockVariable('Price', 'USD/gal', 'Stream-ethanol'),
    MockVariable('Price', 'USD/gal', 'Stream-biodiesel'),
    MockVariable('Price', 'USD/cf', 'Stream-natural gas'),
    MockVariable('Electricity price', 'USD/kWhr', 'biorefinery'),
    MockVariable('Operating days', 'day/yr', 'biorefinery'),
    MockVariable('IRR', '%', 'biorefinery'),
    MockVariable('Price', 'USD/kg', 'Stream-crude glycerol'),
    MockVariable('Price', 'USD/kg', 'Stream-pure glycerol'),
    MockVariable('Reaction time', 'hr', 'Saccharification'),
    MockVariable('Price', 'USD/kg', 'Stream-cellulase'),
    MockVariable('Cellulase loading', 'wt. % cellulose', 'Stream-cellulase'),
    MockVariable('Base cost', 'million USD', 'Pretreatment reactor system'),
    MockVariable('Cane glucose yield', '%', 'Pretreatment and saccharification'),
    MockVariable('Sorghum glucose yield', '%', 'Pretreatment and saccharification'),
    MockVariable('Cane xylose yield', '%', 'Pretreatment and saccharification'),
    MockVariable('Sorghum xylose yield', '%', 'Pretreatment and saccharification'),
    MockVariable('Glucose to ethanol yield', '%', 'Cofermentation'),
    MockVariable('Xylose to ethanol yield', '%', 'Cofermentation'),
    MockVariable('Titer', 'g/L', 'Cofermentation'),
    MockVariable('Productivity', 'g/L/hr', 'Cofermentation'),
    MockVariable('Cane PL content', '% oil', 'oilcane'),
    MockVariable('Sorghum PL content', '% oil', 'oilsorghum'),
    MockVariable('Cane FFA content', '% oil', 'oilcane'),
    MockVariable('Sorghum FFA content', '% oil', 'oilsorghum'),    
    MockVariable('Cane oil content', 'dry wt. %', 'Stream-sugarcane'),
    MockVariable('Relative sorghum oil content', 'dry wt. %', 'Stream-sugarcane'),
    MockVariable('TAG to FFA conversion', '% theoretical', 'Biorefinery'), 
    MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-oilcane'),
    MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-methanol'),
    MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-pure glycerine'),
    MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-cellulase'),
    MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-natural gas'),
)
     
(MFPP, feedstock_consumption, biodiesel_production, ethanol_production, 
 electricity_production, natural_gas_consumption, TCI, 
 heat_exchanger_network_error, GWP_economic, GWP_ethanol, GWP_biodiesel, 
 GWP_crude_glycerol, GWP_electricity, GWP_ethanol_displacement,
 GWP_biofuel_allocation, GWP_ethanol_allocation,
 GWP_biodiesel_allocation, GWP_crude_glycerol_allocation,
 MFPP_derivative, 
 biodiesel_production_derivative, ethanol_production_derivative, 
 electricity_production_derivative, natural_gas_consumption_derivative, 
 TCI_derivative, GWP_economic_derivative, 
 GWP_ethanol_derivative, GWP_biodiesel_derivative,
 GWP_crude_glycerol_derivative, GWP_electricity_derivative
 ) = all_metric_mockups = (
    MockVariable('MFPP', 'USD/ton', 'Biorefinery'),
    MockVariable('Feedstock consumption', 'ton/yr', 'Biorefinery'),
    MockVariable('Biodiesel production', 'Gal/ton', 'Biorefinery'),
    MockVariable('Ethanol production', 'Gal/ton', 'Biorefinery'),
    MockVariable('Electricity production', 'kWhr/ton', 'Biorefinery'),
    MockVariable('Natural gas consumption', 'cf/ton', 'Biorefinery'),
    MockVariable('TCI', '10^6*USD', 'Biorefinery'),
    MockVariable('Heat exchanger network error', '%', 'Biorefinery'),
    MockVariable('GWP', 'kg*CO2*eq / USD', 'Economic allocation'),
    MockVariable('Ethanol GWP', 'kg*CO2*eq / gal', 'Economic allocation'),
    MockVariable('Biodiesel GWP', 'kg*CO2*eq / gal', 'Economic allocation'),
    MockVariable('Crude glycerol GWP', 'kg*CO2*eq / kg', 'Economic allocation'),
    MockVariable('Electricity GWP', 'kg*CO2*eq / MWhr', 'Economic allocation'),
    MockVariable('Ethanol GWP', 'kg*CO2*eq / gal', 'Displacement allocation'),
    MockVariable('Biofuel GWP', 'kg*CO2*eq / GGE', 'Energy allocation'),
    MockVariable('Ethanol GWP', 'kg*CO2*eq / gal', 'Energy allocation'),
    MockVariable('Biodiesel GWP', 'kg*CO2*eq / gal', 'Energy allocation'),
    MockVariable('Crude-glycerol GWP', 'kg*CO2*eq / kg', 'Energy allocation'),
    MockVariable('MFPP derivative', 'USD/ton', 'Biorefinery'),
    MockVariable('Biodiesel production derivative', 'Gal/ton', 'Biorefinery'),
    MockVariable('Ethanol production derivative', 'Gal/ton', 'Biorefinery'),
    MockVariable('Electricity production derivative', 'kWhr/ton', 'Biorefinery'),
    MockVariable('Natural gas consumption derivative', 'cf/ton', 'Biorefinery'),
    MockVariable('TCI derivative', '10^6*USD', 'Biorefinery'),
    MockVariable('GWP derivative', 'kg*CO2*eq / USD', 'Economic allocation'),
    MockVariable('Ethanol GWP derivative', 'kg*CO2*eq / gal', 'Ethanol'),
    MockVariable('Biodiesel GWP derivative', 'kg*CO2*eq / gal', 'Biodiesel'),
    MockVariable('Crude glycerol GWP derivative', 'kg*CO2*eq / kg', 'Crude glycerol'),
    MockVariable('Electricity GWP derivative', 'kg*CO2*eq / MWhr', 'Electricity'),
)

tea_monte_carlo_metric_mockups = (
    MFPP, 
    TCI,
    ethanol_production,
    biodiesel_production,
    electricity_production,
    natural_gas_consumption
)

tea_monte_carlo_derivative_metric_mockups = (
    MFPP_derivative, 
    TCI_derivative,
    ethanol_production_derivative,
    biodiesel_production_derivative,
    electricity_production_derivative,
    natural_gas_consumption_derivative,
)

lca_monte_carlo_metric_mockups = (
    GWP_economic,
    GWP_ethanol,
    GWP_biodiesel,
    GWP_electricity,
    GWP_crude_glycerol,
)

lca_monte_carlo_derivative_metric_mockups = (
    GWP_economic_derivative,
    GWP_ethanol_derivative,
    GWP_biodiesel_derivative, 
    GWP_electricity_derivative,
    GWP_crude_glycerol_derivative,
)
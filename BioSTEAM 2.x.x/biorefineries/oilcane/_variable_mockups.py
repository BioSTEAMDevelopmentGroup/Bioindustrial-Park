# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:44:17 2021

@author: yrc2
"""
import biosteam as bst

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
 set_methanol_GWP, set_pure_glycerine_GWP, set_cellulase_GWP,
 ) = all_parameter_mockups = (
    bst.MockVariable('Oil retention', '%', 'Stream-sugarcane'),
    bst.MockVariable('Bagasse oil extraction efficiency', '%', 'Stream-sugarcane'),
    bst.MockVariable('Capacity', 'ton/hr', 'Stream-sugarcane'),
    bst.MockVariable('Price', 'USD/gal', 'Stream-ethanol'),
    bst.MockVariable('Price', 'USD/gal', 'Stream-biodiesel'),
    bst.MockVariable('Price', 'USD/cf', 'Stream-natural gas'),
    bst.MockVariable('Electricity price', 'USD/kWhr', 'biorefinery'),
    bst.MockVariable('Operating days', 'day/yr', 'biorefinery'),
    bst.MockVariable('IRR', '%', 'biorefinery'),
    bst.MockVariable('Price', 'USD/kg', 'Stream-crude glycerol'),
    bst.MockVariable('Price', 'USD/kg', 'Stream-pure glycerol'),
    bst.MockVariable('Reaction time', 'hr', 'Saccharification'),
    bst.MockVariable('Price', 'USD/kg', 'Stream-cellulase'),
    bst.MockVariable('Cellulase loading', 'wt. % cellulose', 'Stream-cellulase'),
    bst.MockVariable('Base cost', 'million USD', 'Pretreatment reactor system'),
    bst.MockVariable('Cane glucose yield', '%', 'Pretreatment and saccharification'),
    bst.MockVariable('Sorghum glucose yield', '%', 'Pretreatment and saccharification'),
    bst.MockVariable('Cane xylose yield', '%', 'Pretreatment and saccharification'),
    bst.MockVariable('Sorghum xylose yield', '%', 'Pretreatment and saccharification'),
    bst.MockVariable('Glucose to ethanol yield', '%', 'Cofermentation'),
    bst.MockVariable('Xylose to ethanol yield', '%', 'Cofermentation'),
    bst.MockVariable('Titer', 'g/L', 'Cofermentation'),
    bst.MockVariable('Productivity', 'g/L/hr', 'Cofermentation'),
    bst.MockVariable('Cane PL content', '% oil', 'oilcane'),
    bst.MockVariable('Sorghum PL content', '% oil', 'oilsorghum'),
    bst.MockVariable('Cane FFA content', '% oil', 'oilcane'),
    bst.MockVariable('Sorghum FFA content', '% oil', 'oilsorghum'),    
    bst.MockVariable('Cane oil content', 'dry wt. %', 'Stream-sugarcane'),
    bst.MockVariable('Relative sorghum oil content', 'dry wt. %', 'Stream-sugarcane'),
    bst.MockVariable('TAG to FFA conversion', '% theoretical', 'Biorefinery'), 
    bst.MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-oilcane'),
    bst.MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-methanol'),
    bst.MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-pure glycerine'),
    bst.MockVariable('GWP', 'kg*CO2-eq/kg', 'Stream-cellulase'),
)
     
(MFPP, biodiesel_production, ethanol_production, electricity_production, 
 natural_gas_consumption, TCI, feedstock_consumption, 
 heat_exchanger_network_error, GWP_economic, GWP_ethanol, GWP_biodiesel, 
 GWP_crude_glycerol, GWP_electricity, MFPP_derivative,
 biodiesel_production_derivative, ethanol_production_derivative, 
 electricity_production_derivative, natural_gas_consumption_derivative, 
 TCI_derivative, GWP_economic_derivative, 
 GWP_biodiesel_derivative, GWP_ethanol_derivative, GWP_electricity_derivative, 
 GWP_crude_glycerol_derivative,
 ) = all_metric_mockups = (
    bst.MockVariable('MFPP', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('Biodiesel production', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Ethanol production', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Electricity production', 'kWhr/ton', 'Biorefinery'),
    bst.MockVariable('Natural gas consumption', 'cf/ton', 'Biorefinery'),
    bst.MockVariable('TCI', '10^6*USD', 'Biorefinery'),
    bst.MockVariable('Feedstock consumption', 'ton/yr', 'Biorefinery'),
    bst.MockVariable('Heat exchanger network error', '%', 'Biorefinery'),
    bst.MockVariable('GWP', 'kg*CO2*eq / USD', 'Economic allocation'),
    bst.MockVariable('Ethanol GWP', 'kg*CO2*eq / gal', 'Ethanol'),
    bst.MockVariable('Biodiesel GWP', 'kg*CO2*eq / gal', 'Biodiesel'),
    bst.MockVariable('Crude glycerol GWP', 'kg*CO2*eq / kg', 'Crude glycerol'),
    bst.MockVariable('Electricity GWP', 'kg*CO2*eq / MWhr', 'Electricity'),
    bst.MockVariable('MFPP derivative', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('Biodiesel production derivative', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Ethanol production derivative', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Electricity production derivative', 'kWhr/ton', 'Biorefinery'),
    bst.MockVariable('Natural gas consumption derivative', 'cf/ton', 'Biorefinery'),
    bst.MockVariable('TCI derivative', '10^6*USD', 'Biorefinery'),
    bst.MockVariable('GWP derivative', 'kg*CO2*eq / USD', 'Economic allocation'),
    bst.MockVariable('Ethanol GWP derivative', 'kg*CO2*eq / gal', 'Ethanol'),
    bst.MockVariable('Biodiesel GWP derivative', 'kg*CO2*eq / gal', 'Biodiesel'),
    bst.MockVariable('Crude glycerol GWP derivative', 'kg*CO2*eq / kg', 'Crude glycerol'),
    bst.MockVariable('Electricity GWP derivative', 'kg*CO2*eq / MWhr', 'Electricity'),
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
    GWP_biodiesel, 
    GWP_ethanol,
    GWP_electricity,
    GWP_crude_glycerol,
)

lca_monte_carlo_derivative_metric_mockups = (
    GWP_biodiesel_derivative, 
    GWP_ethanol_derivative,
    GWP_electricity_derivative,
    GWP_crude_glycerol_derivative,
)
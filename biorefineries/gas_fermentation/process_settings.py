# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst

__all__ = (
    'load_process_settings',
)
GWP = 'GWP'
characterization_factors = {
    # 'Electricity': 0.38, # [kg*CO2*eq / kWhr] From GREET 2020; NG-Fired Simple-Cycle Gas Turbine CHP Plant
    # https://www2.gov.bc.ca/assets/gov/business/natural-resource-industries/reports/carbon_intensity_of_hydrogen_production_methods.pdf
    'H2': 1.8023, # [kg*CO2*eq / kg] From GREET 2022; Compressed G.H2 produced from PEM electrolysis (solar or wind)
    'Cornstover': 0.03819, # [kg*CO2*eq / kg] From GREET 2023; Miscanthus for ethanol production
    'Miscanthus': 0.07606, # [kg*CO2*eq / kg] From GREET 2023; Miscanthus for ethanol production
    'Flue gas': -0.07134147183700712, # [kg*CO2*eq / kg] Assume C is displaced from flue gas
    'Natural gas': 0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
    'Ethyl acetate': 4.2677, # [kg*CO2*eq / kg] From GREET 2022; Ethyl acetate production pathway
    'Hexane': 0.55002, # [kg*CO2*eq / kg] Ecoinvent 3.6 Cut off; IPCC 2013; Global market for hexane
}

def load_process_settings():
    settings = bst.settings
    settings.define_impact_indicator(GWP, 'kg*CO2e')
    # for i in settings.heating_agents[:-1]:
    #     settings.set_utility_agent_CF(
    #         i.ID, GWP, 0.2958e-3, 'kJ' # [kg*CO2*eq / kJ] From GREET 2020; Natural gas, on-site
    #     )
    settings.CEPCI = 816.0 # 2022
    bst.settings.electricity_price = 0.060 # Maryland solar REC (renewable energy certificates) # https://escholarship.org/uc/item/80n4q8xc
    bst.PowerUtility.set_CF('GWP', 0)
    cw = settings.get_cooling_agent("cooling_water")
    cw.set_property('T', 70, 'degF') # https://www.sciencedirect.com/topics/earth-and-planetary-sciences/cooling-tower
    hps = settings.get_heating_agent("high_pressure_steam")
    hps.heat_transfer_efficiency = 0.85
    # hps.regeneration_price = 0.08064
    hps.T = 529.2
    hps.P = 44e5
    mps = settings.get_heating_agent("medium_pressure_steam")
    mps.heat_transfer_efficiency = 0.90
    # mps.regeneration_price = 0.07974
    mps.T = 480.3
    mps.P = 18e5
    lps = settings.get_heating_agent("low_pressure_steam")
    lps.heat_transfer_efficiency = 0.95
    # lps.regeneration_price = 0.06768
    lps.T = 428.6
    lps.P = 55e4
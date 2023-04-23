# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                     Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

__all__ = (
    'BiorefinerySettings',
    'default_feedstock_composition',
    'default_load_utility_settings',
    'default_stream_GWP_CFs',
    'default_stream_prices',
    )


class BiorefinerySettings:
    '''
    Containing the general parameters needed in biorefinery design/TEA/LCA
    that can be modified to update the biorefinery.
    '''
    
    def __init__(
            self,
            system_ID='corn_sys',
            CEPCI=567, # 2013 chemical engineering plant capital index
            load_utility_settings=None,
            feedstock_composition={},
            feedstock_hourly_mass_flow=46375, # kg / hr
            process_parameters={},
            stream_prices={},
            stream_GWP_CFs={},
            other_unit_parameters={},
            tea_parameters={},
            **kwargs,
            ):
        self.system_ID = system_ID
        self.CEPCI = CEPCI
        self.load_utility_settings = load_utility_settings or default_load_utility_settings
        self.feedstock_composition = feedstock_composition or default_feedstock_composition
        self.feedstock_hourly_mass_flow = feedstock_hourly_mass_flow
        self.process_parameters = process_parameters or default_process_parameters
        self.stream_prices = stream_prices or default_stream_prices
        self.stream_GWP_CFs = stream_GWP_CFs or default_stream_GWP_CFs
        self.other_unit_parameters = other_unit_parameters
        self.tea_parameters = tea_parameters
        for k, v in kwargs: setattr(self, k, v)
        
    def load_process_settings(self):
        '''Set capital index and load utility settings'''
        bst.CE = self.CEPCI
        self.load_utility_settings()
        

def default_load_utility_settings(
        power_price=0.07, # Value from Chinmay's report; 0.07, in original
        power_GWP_CFs=(1., 1.,), # consumption, production #!!! needs to be updated
        ):
    bst.process_tools.default_utilities()
    bst.PowerUtility.price = power_price
    bst.PowerUtility.characterization_factors['GWP'] = power_GWP_CFs
    HeatUtility = bst.HeatUtility
    lps = HeatUtility.get_agent('low_pressure_steam')
    Water = lps.chemicals.Water
    lps.T = 152 + 273.15
    lps.P = Water.Psat(lps.T)
    lps.heat_transfer_efficiency = 1.0
    lps.regeneration_price = default_stream_prices['Steam'] * Water.MW
    for agent in HeatUtility.heating_agents:
        agent.regeneration_price = default_stream_prices['Steam'] * Water.MW
    cw = HeatUtility.get_agent('cooling_water')
    cw.regeneration_price = 0.073e-3 * Water.MW
    cw.T = 25. + 273.15


# Feedstock weight composition
default_feedstock_composition = {
    'Starch': 0.612,
    'Water': 0.15,
    'Fiber': 0.1067,
    'SolubleProtein': 0.034,
    'InsolubleProtein': 0.0493,
    'Oil': 0.034,
    'Ash': 0.014
    }


# Process parameters
default_process_parameters = {
    'slurry_solids_content': 0.311,   # g Corn / g Total
    'slurry_ammonia_loading': 0.002,  # g Ammonia / g dry Corn
    'slurry_lime_loading': 0.00012,   # g Lime / g dry Corn
    'liquefaction_alpha_amylase_loading': 0.00082,    # g Enzyme / g dry Corn
    'saccharification_sulfuric_acid_loading': 0.001, # g H2SO4 / g dry Corn
    'saccharification_gluco_amylase_loading': 0.0011, # g Enzyme / g dry Mash
    'scrubber_wash_water_over_vent': 1.21, # g Water / g Vent
    'yeast_loading': 8e-5, # g Yeast (5 wt %) / g Mash
    }

# Characterization factors for 100-year global warming potential
#!!! Numbers are all fake now
default_stream_GWP_CFs = {
    'Ethanol': 1,
    'Corn': 1,
    'DDGS': 1,
    'Yeast': 1,
    'Enzyme': 1,
    'Denaturant': 1,
    'Sulfuric acid': 1,
    'Ammonia': 1,
    'Caustic': 1,
    'Lime': 1,
    'Steam': 1,
    'Crude oil': 1,
}


# Raw material price (USD/kg)
default_stream_prices = {
    'Ethanol': 0.48547915353569393, 
    'Corn': 0.13227735731092652, # Value from Chinmay's report; 0.08476585075177462 in original
    'DDGS': 0.12026, # Value from Chinmay's report; 0.09687821462905594, in original
    'Yeast': 1.86, # 1.86 in Chinmay's report; 0.907310526171629 in original
    'Enzyme': 2.25, # Value from Chinmay's report; 3.20739717428309 in original
    'Denaturant': 0.4355069727002459,
    'Sulfuric acid': 0.08971711759613593,
    'Ammonia': 0.4485966110937889,
    'Caustic': 0.14952852932689323,
    'Lime': 0.19937504680689405,
    'Steam': 12.86e-3, # Value from Chinmay's report; 17.08e-3 in original lipidcane report 
    'Crude oil': 0.56,
}



# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

__all__ = ('load_process_settings', 'price',)

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 567 # 2013
    bst.PowerUtility.price = 0.05
    HeatUtility = bst.HeatUtility
    lps = HeatUtility.get_agent('low_pressure_steam')
    Water = lps.chemicals.Water
    lps.T = 152 + 273.15
    lps.P = Water.Psat(lps.T)
    lps.heat_transfer_efficiency = 0.95 
    lps.regeneration_price = 17.08e-3 * Water.MW
    cw = HeatUtility.get_agent('cooling_water')
    cw.regeneration_price = 0.073e-3 * Water.MW
    cw.T = 25. + 273.15


# Raw material price (USD/kg)
price = {
    'Ethanol': 0.48547915353569393, # 0.789,
    'Corn': 0.08476585075177462,
    'DDGS': 0.09687821462905594,
    'Yeast': 0.907310526171629,
    'Enzyme': 3.20739717428309,
    'Denaturant': 0.4355069727002459,
    'Sulfuric acid': 0.08971711759613593,
    'Ammonia': 0.4485966110937889,
    'Caustic': 0.14952852932689323,
    'Lime': 0.19937504680689405,
    'Steam': 17.08e-3,
}
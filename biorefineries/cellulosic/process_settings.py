# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst

__all__ = (
    'load_process_settings',
)

def load_process_settings():
    bst.process_tools.default_utilities()
    settings = bst.settings 
    settings.CEPCI = 525.4
    settings.electricity_price = 0.0572 # USD/kWh
    lps = settings.get_heating_agent('low_pressure_steam')
    lps.heat_transfer_efficiency = 0.90
    lps.T = 529.2
    lps.P = 44e5
    lps.regeneration_price = 0.
    cw = settings.get_cooling_agent('cooling_water')
    cw.T = 28 + 273.15
    cw.T_limit = cw.T + 9
    cw.regeneration_price = 0.
    settings.get_cooling_agent('chilled_water').heat_transfer_price = 0
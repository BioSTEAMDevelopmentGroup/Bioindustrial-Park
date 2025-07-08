# -*- coding: utf-8 -*-
"""
Created on 2025-07-07 14:55:35

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst

__all__ = ('load_process_settings',)

# %% Process settings

def load_process_settings():
    settings = bst.settings
    bst.process_tools.default_utilities()
    settings.CEPCI = 607.5 # 2019
    settings.electricity_price = 0.065
    steam_utility = settings.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    settings.get_agent('cooling_water').regeneration_price = 0
    settings.get_agent('chilled_water').heat_transfer_price = 0

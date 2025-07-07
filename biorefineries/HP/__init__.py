#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
All units are explicitly defined here for transparency and easy reference.
Naming conventions:
    D = Distillation column
    C = Crystallization
    AC = Adsorption column
    F = Flash tank or multiple-effect evaporator
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
Processes:
    100: Feedstock preprocessing
    200: Feedstock pretreatment and juicing
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Storage
    700: Co-heat and power
    800: Cooling utility generation
    900: Miscellaneous facilities
    1000: Heat exchanger network

"""

#%%
from .models import load

__all__ = ['load_model']

def load_model(feedstock, product, fermentation_performance):
    model, system, spec, tea, lca, get_adjusted_MSP, simulate_and_print, TEA_breakdown, unit_groups_dict =\
        load.load_HP_model(feedstock, product, fermentation_performance)
    globals().update({
        'system':system, 'spec':spec, 'tea':tea, 'lca':lca, 'get_adjusted_MSP':get_adjusted_MSP,
        'simulate_and_print':simulate_and_print, 'TEA_breakdown':TEA_breakdown,
        'unit_groups_dict': unit_groups_dict,
        })

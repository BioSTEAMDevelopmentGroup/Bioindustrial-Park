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

def load_TAL_model(mode='A'):
    implemented_modes = ['A', 'B', 'C', 'D']
    if not mode in implemented_modes: raise ValueError(f"Mode must be one of {implemented_modes}, not {mode}.")
    system, spec, tea, lca, get_adjusted_MSP, simulate_and_print, TEA_breakdown, run_TAL_uncertainty_analysis = load.load_TAL_model(mode=mode)
    globals().update({
        'system':system, 'spec':spec, 'tea':tea, 'lca':lca, 'get_adjusted_MSP':get_adjusted_MSP,
        'simulate_and_print':simulate_and_print, 'TEA_breakdown':TEA_breakdown,
        'run_TAL_uncertainty_analysis':run_TAL_uncertainty_analysis,
        })
    
def load_KS_model(mode='A'):
    implemented_modes = ['THF_Ethanol_A', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    if not mode in implemented_modes: raise ValueError(f"Mode must be one of {implemented_modes}, not {mode}.")
    system, spec, tea, lca, get_adjusted_MSP, simulate_and_print, TEA_breakdown = load.load_KS_model(mode=mode)
    globals().update({
        'system':system, 'spec':spec, 'tea':tea, 'lca':lca, 'get_adjusted_MSP':get_adjusted_MSP,
        'simulate_and_print':simulate_and_print, 'TEA_breakdown':TEA_breakdown,
        })
    
__all__ = ['load_TAL_model', 'load_KS_model']

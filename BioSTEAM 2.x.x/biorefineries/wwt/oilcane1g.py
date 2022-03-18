#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/oilcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

info = {
    'abbr': 'oc1g',
    'WWT_ID': '11',
    'is2G': False,
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_oc1g_comparison_systems():
    from biorefineries.wwt import create_comparison_systems
    from biorefineries.oilcane import (
        create_chemicals,
        create_oilcane_to_biodiesel_and_ethanol_1g as create_system,
        create_tea,
        load_process_settings,
        )
    functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    sys_dct = {
        'create_system': {'operating_hours': 24*200},
        'rename_storage_to': 1000,
        'create_wastewater_process': {'skip_AeF': True},
        # `vinasse`, `wastewater`, `fiber_fines`,
        # COD of `evaporator_condensate` is only ~20 mg/L
        'ww_streams': (('M403', 0), ('M603', 0), ('U206', 1)),
        'solids_streams': (('M701', 0), ('U205', 0)), # the second one is `filter_cake`
        'BT': 'BT701',
        'new_wwt_connections': {'solids': ('BT701', 0), 'biogas': ('BT701', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)
    return exist_sys, new_sys


def simulate_oc1g_systems():
    from biorefineries.wwt import simulate_systems
    exist_sys, new_sys = create_oc1g_comparison_systems()
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys



# %%

# =============================================================================
# Models
# =============================================================================


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    exist_sys, new_sys = simulate_oc1g_systems()
    # exist_model, new_model = create_oc1g_comparison_models()
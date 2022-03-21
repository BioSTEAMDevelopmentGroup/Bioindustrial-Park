#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the lactic acid biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lactic
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

info = {
    'abbr': 'la',
    'WWT_ID': '5',
    'is2G': True,
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_la_comparison_systems():
    # Create from scratch
    from biorefineries.wwt import create_comparison_systems, add_wwt_chemicals
    from biorefineries.lactic import (
        create_chemicals,
        create_system,
        create_tea,
        load_process_settings,
        get_splits,
        )
    # Add WWT chemicals to the existing splits array,
    # splits of chemicals that do now exist in the original chemicals obj
    # will be copied from the splits of the corresponding group
    la_chems = add_wwt_chemicals(create_chemicals())
    def create_new_splits(original_splits):
        new_splits = la_chems.zeros()
        new_splits[la_chems.indices(('Bisulfite', 'CitricAcid', 'HCl', 'NaOCl'))] = \
            original_splits[la_chems.index('NaOH')]
        return new_splits
    cell_mass_split, gypsum_split, AD_split, MB_split = get_splits(la_chems)
    new_cell_mass_split = create_new_splits(cell_mass_split)
    new_gypsum_split = create_new_splits(gypsum_split)

    functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    sys_dct = {
        'create_system': {'cell_mass_split': new_cell_mass_split, 'gypsum_split': new_gypsum_split},
        'BT': 'CHP',
        'new_wwt_connections': {'sludge': ('M601', 0), 'biogas': ('CHP', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)

    # # IRR doesn't match with direct loading as closely as creating from scratch
    # from biorefineries.wwt import create_comparison_systems
    # from biorefineries import lactic as la
    # sys_dct = {
    #     'load': {'print_results': False},
    #     'system_name': 'lactic_sys',
    #     'BT': 'CHP',
    #     'new_wwt_connections': {'sludge': ('M601', 0), 'biogas': ('CHP', 1)},
    #     }
    # exist_sys, new_sys = create_comparison_systems(info, la, sys_dct, from_load=True)

    return exist_sys, new_sys


def simulate_la_systems():
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_la_comparison_systems()
    # ~504 mg/L COD, soluble lignin, arabinose, and galactose all >10%,
    # lactic acid, extract, xylose, and mannose ~5-10%
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
    exist_sys, new_sys = simulate_la_systems()
    # exist_model, new_model = create_la_comparison_models()
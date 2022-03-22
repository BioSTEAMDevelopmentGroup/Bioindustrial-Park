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
    'abbr': 'oc2g',
    'WWT_ID': '5',
    'is2G': True,
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_oc2g_comparison_systems():
    # # Does not work for oilcane biorefineries due to the many settings
    # # not included in the system creation function
    # from biorefineries.wwt import create_comparison_systems
    # from biorefineries.oilcane import (
    #     create_chemicals,
    #     create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation as create_system,
    #     create_tea,
    #     load_process_settings,
    #     )
    # functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    # sys_dct = {
    #     'create_system': {'operating_hours': 24*200},
    #     'rename_storage_to': 1100,
    #     'BT': 'BT701',
    #     'new_wwt_connections': {'sludge': ('BT701', 0), 'biogas': ('BT701', 1)},
    #     }
    # exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)

    from biorefineries.wwt import create_comparison_systems
    from biorefineries import oilcane as oc
    sys_dct = {
        'load': {'name': 'O2', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'BT': 'BT701',
        'new_wwt_connections': {'sludge': ('BT701', 0), 'biogas': ('BT701', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct, from_load=True)

    return exist_sys, new_sys


def simulate_oc2g_systems():
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_oc2g_comparison_systems()
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_oc2g_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_oc2g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'oilcane',
        'FERM_product': 'ethanol',
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_solids_mixer': 'M301',
        'PT_rx': 'R301',
        'EH_mixer': 'M401',
        'EH_rx': 'U401',
        'fermentor': 'R401',
        'TE_rx': 'U802',
        'isplit_efficiency_is_reversed': True,
        'bagasse_oil_extraction': 'U402',
        'bagasse_oil_retention': 'U201',
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification', 2),
            'FERM glucan-to-product': ('cofermentation', 0),
            'FERM xylan-to-product': ('cofermentation', 1),
            'FERM oil-to-FFA': ('oil_reaction', 0), # not fermentation, but happens in the fermentor
            'TE oil-to-product': ('transesterification', (0, 1, 2)),
            },
        'BT': 'BT701',
        'BT_eff': ('boiler_efficiency', 'turbogenerator_efficiency'),
        'wwt_system': 'exist_sys_wwt',
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['sludge'] = 'sludge'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model_dct['new_wwt_ID'] = info['WWT_ID']
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_oc2g_models(**kwargs):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_oc2g_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **kwargs)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_oc2g_systems()
    # exist_model, new_model = create_oc2g_comparison_models()
    exist_model, new_model = evaluate_oc2g_models(N=10)
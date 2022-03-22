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
    'abbr': 'sc1g',
    'WWT_ID': '8',
    'is2G': False,
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_sc1g_comparison_systems():
    # # Does not work for oilcane biorefineries due to the many settings
    # # not included in the system creation function
    # from biorefineries.wwt import create_comparison_systems
    # from biorefineries.oilcane import (
    #     create_chemicals,
    #     create_sugarcane_to_ethanol_system as create_system,
    #     create_tea,
    #     load_process_settings,
    #     )
    # functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    # sys_dct = {
    #     'create_system': {'operating_hours': 24*200, 'use_area_convention': True, 'pellet_bagasse': True},
    #     'rename_storage_to': 700,
    #     'create_wastewater_process': {'skip_AeF': True},
    #     # `vinasse`, `fiber_fines`,
    #     # not using `wastewater` as it contains `evaporator_condensate` (all water)
    #     'ww_streams': (('H302', 1), ('U211', 0)),
    #     'solids_streams': (('U207', 0), ('U210', 0)), # `bagasse`, `filter_cake`
    #     'BT': 'BT401',
    #     'new_wwt_connections': {'solids': ('BT401', 0), 'biogas': ('BT401', 1)},
    #     }
    # exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)

    from biorefineries.wwt import create_comparison_systems
    from biorefineries import oilcane as oc
    sys_dct = {
        'load': {'name': 'S1', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'create_wastewater_process': {'skip_AeF': True},
        # `vinasse`, `fiber_fines`,
        # not using `wastewater` as it contains `evaporator_condensate` (all water)
        'ww_streams': (('H302', 1), ('U211', 1)),
        'solids_streams': (('U207', 0), ('U210', 0)), # `bagasse`, `filter_cake`
        'BT': 'BT401',
        'new_wwt_connections': {'solids': ('BT401', 0), 'biogas': ('BT401', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct, from_load=True)

    return exist_sys, new_sys


def simulate_sc1g_systems():
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_sc1g_comparison_systems()
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_sc1g_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_sc1g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'sugarcane',
        'FERM_product': 'ethanol',
        'PT_rx': 'R301',
        'fermentor': 'R301',
        'reactions': {
            'PT glucan-to-glucose': ('hydrolysis_reaction', ),
            'FERM glucan-to-product': ('fermentation_reaction', ),
            },
        'BT': 'BT401',
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


def evaluate_sc1g_models(**kwargs):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_sc1g_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **kwargs)



# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_sc1g_systems()
    # exist_model, new_model = create_sc1g_comparison_models()
    exist_model, new_model = evaluate_sc1g_models(N=10)
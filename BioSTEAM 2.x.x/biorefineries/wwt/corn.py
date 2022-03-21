#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the corn biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/corn
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


info = {
    'abbr': 'cn',
    'WWT_ID': '7',
    'is2G': False,
    'add_CHP': True,
    'ww_price': -0.003, # the default -0.03 leads to two solutions
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_cn_comparison_systems():
    # # Create from scratch
    # import biosteam as bst
    # from biorefineries.wwt import create_comparison_systems
    # from biorefineries.corn import (
    #     create_chemicals,
    #     create_system,
    #     create_tea,
    #     load_process_settings
    #     )
    # functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    # sys_dct = {
    #     'create_system': {'flowsheet': bst.main_flowsheet},
    #     'create_wastewater_process': {'skip_AeF': True},
    #     'ww_streams': (('MH103', 1), ('MX5', 0)),
    #     }
    # exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)

    from biorefineries.wwt import create_comparison_systems
    from biorefineries import corn as cn
    sys_dct = {
        'system_name': 'corn_sys',
        'create_wastewater_process': {'skip_AeF': True},
        'ww_streams': (('MH103', 1), ('MX5', 0)),
        }
    exist_sys, new_sys = create_comparison_systems(info, cn, sys_dct, from_load=True)

    return exist_sys, new_sys


def simulate_cn_systems():
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_cn_comparison_systems()
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_cn_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_cn_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'feedstock': 'corn',
        'PT_rx': 'V310',
        'fermentor': 'V405',
        'reactions': {
            'PT glucan-to-glucose': ('reaction',),
            'FERM glucan-to-product': ('reaction',),
            },
        'BT': None,
        'wwt_system': 'exist_sys_wwt',
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['BT'] = 'CHP'
    new_model_dct['sludge'] = 'sludge'
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model_dct['WWT_ID'] = info['WWT_ID']
    new_model = create_comparison_models(new_sys, new_model_dct)

    return exist_model, new_model


def evaluate_cn_models(**kwargs):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_cn_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **kwargs)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_cn_systems()
    # exist_model, new_model = create_cn_comparison_models()
    exist_model, new_model = evaluate_cn_models(N=10, notify=100)
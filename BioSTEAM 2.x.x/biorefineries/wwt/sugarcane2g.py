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
    'abbr': 'sc2g',
    'WWT_ID': '5',
    'is2G': True,
    'FERM_product': 'ethanol',
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_sc2g_comparison_systems(default_BD=True):
    from biorefineries.wwt import create_comparison_systems
    from biorefineries import oilcane as oc
    BD = {} if not default_BD else 1.
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': BD,})
    sys_dct = {
        'load': {'name': 'S2', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'create_wastewater_process': wwt_kwdct,
        'BT': 'BT701',
        'new_wwt_connections': {'sludge': ('M701', 0), 'biogas': ('BT701', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct, from_load=True)
    return exist_sys, new_sys


def simulate_sc2g_systems(**sys_kwdct):
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_sc2g_comparison_systems(**sys_kwdct)
    # If using conservative biodegradability,
    # ~123 mg/L COD, mostly (~100/>80%) due to soluble lignin and arabinose
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_sc2g_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_sc2g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'sugarcane',
        'FERM_product': info['FERM_product'],
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_solids_mixer': 'M301',
        'PT_rx': 'R301',
        'EH_mixer': 'M401',
        'EH_rx': 'U401',
        'fermentor': 'R401',
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification', 2),
            'FERM glucan-to-product': ('cofermentation', 0),
            'FERM xylan-to-product': ('cofermentation', 1),
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


def evaluate_sc2g_models(**eval_kwdct):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_sc2g_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **eval_kwdct)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_sc2g_systems(default_BD=True)
    # exist_model, new_model = create_sc2g_comparison_models()
    exist_model, new_model = evaluate_sc2g_models(N=1000)
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
    'FERM_product': 'lactic_acid',
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

#
def create_la_comparison_systems(default_BD=True):
    from biorefineries.wwt import create_comparison_systems
    from biorefineries import lactic as la
    BD = {} if not default_BD else 1.
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': BD,})
    sys_dct = {
        # 'load': {'print_results': False}, # need to run `simulate_and_print` for results to match
        'system_name': 'lactic_sys',
        'create_wastewater_process': wwt_kwdct,
        'BT': 'CHP',
        'new_wwt_connections': {'sludge': ('M601', 0), 'biogas': ('CHP', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, la, sys_dct, from_load=True)
    return exist_sys, new_sys


def simulate_la_systems(**sys_kwdct):
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_la_comparison_systems(**sys_kwdct)
    # If using conservative biodegradability,
    # ~504 mg/L COD, soluble lignin, arabinose, and galactose all >10%,
    # lactic acid, extract, xylose, and mannose ~5-10%
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_la_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_la_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'feedstock',
        'FERM_product': info['FERM_product'],
        'sludge': 'wastes_to_CHP',
        'biogas': 'biogas',
        'PT_solids_mixer': 'M202',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R301',
        'reactions': {
            'PT glucan-to-glucose': ('pretreatment_rxns', 0),
            'PT xylan-to-xylose': ('pretreatment_rxns', 4),
            'EH glucan-to-glucose': ('saccharification_rxns', 2),
            'FERM glucan-to-product': ('cofermentation_rxns', 0),
            'FERM xylan-to-product': ('cofermentation_rxns', 3),
            },
        'BT': 'CHP',
        'BT_eff': ('B_eff', 'TG_eff'),
        'wwt_system': 'exist_sys_wwt',
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model_dct['new_wwt_ID'] = info['WWT_ID']
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_la_models(**eval_kwdct):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_la_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **eval_kwdct)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_la_systems(default_BD=True)
    # exist_model, new_model = create_la_comparison_models()
    exist_model, new_model = evaluate_la_models(N=1000)
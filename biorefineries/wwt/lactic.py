#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the lactic acid biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lactic
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries import lactic as la
from biorefineries.wwt import (
    create_comparison_systems, simulate_systems,
    create_comparison_models,evaluate_models,
    )

info = {
    'abbr': 'la',
    'WWT_ID': '5',
    'is2G': True,
    'FERM_product': ['lactic_acid',],
    'add_BT': False,
    'ww_price': None,
    }

CF_dct = {
    ##### Feeds #####
    'ammonia': ('NH4OH',), # NH4OH
    'caustic_R502': ('NaOH',),
    'CSL': ('CSL',),
    'enzyme_M301': ('Cellulase',), # this is pure enzyme
    'ethanol': ('Ethanol',),
    'feedstock': ('CornStover',),
    'lime': ('Lime',),
    'lime_CHP': ('Lime',),
    'natural_gas': ('CH4',),
    'polymer_R502': ('Polymer',),
    'sulfuric_acid': ('H2SO4',),
    ##### Co-products #####
    # following the lactic paper baseline, gypsum not included
    }


# %%

# =============================================================================
# Systems
# =============================================================================

#
def create_la_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    sys_dct = {
        # 'load': {'print_results': False}, # need to run `simulate_and_print` for results to match
        'system_name': 'lactic_sys',
        'create_wastewater_process': wwt_kwdct,
        'BT': 'BT',
        'new_wwt_connections': {'sludge': ('M601', 0), 'biogas': ('BT', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, la, sys_dct)
    return exist_sys, new_sys


def simulate_la_systems(**sys_kwdct):
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
    exist_sys, new_sys = create_la_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'feedstock',
        'FERM_product': info['FERM_product'],
        'sludge': 'solids_to_boiler',
        'biogas': 'biogas',
        'PT_solids_mixer': 'M203',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R301',
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification_rxns', 2),
            'FERM glucose-to-product': ('cofermentation_rxns', 0),
            'FERM xylose-to-product': ('cofermentation_rxns', 3),
            },
        'BT': 'BT',
        'BT_eff': ('boiler_efficiency', 'turbogenerator_efficiency'),
        'wwt_system': 'exist_sys_wwt',
        'wwt_ID': info['WWT_ID'],
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_la_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_la_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_la_systems(biodegradability=1)
    # exist_model, new_model = create_la_comparison_models()
    exist_model, new_model = evaluate_la_models(
        # include_baseline=True,
        # include_uncertainty=True,
        include_BMP=True,
        # N_uncertainty=1000,
        # uncertainty_skip_exist=True, # for testing
        N_BMP=100,
        # BMPs=(0.5, 0.9499,), # for testing, 0.9499 allows for minor error
        )
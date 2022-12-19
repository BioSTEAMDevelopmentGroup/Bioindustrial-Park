#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery (configuration S1/-1):
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/oilcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries import oilcane as oc
from biorefineries.wwt import (
    create_comparison_systems, simulate_systems,
    create_comparison_models,evaluate_models,
    )

info = {
    'abbr': 'sc1g',
    'WWT_ID': '8',
    'is2G': False,
    'FERM_product': ['ethanol',],
    'add_BT': False,
    'ww_price': None,
    }

CF_dct = {
    ##### Feeds #####
    'denaturant': ('Denaturant'),
    # 'dryer_natural_gas': ('CH4',), # used to be in the registry, but not any more
    'H3PO4': ('H3PO4',), # moisture content already adjusted
    'lime': ('CaO', 0.046), # CaO and water
    'natural_gas': ('CH4',), # actually empty
    'polymer': ('Polymer',),
    'sugarcane': ('Sugarcane',), # moisture content already adjusted
    ##### Co-products #####
    # `fiber_fines`, `wastewater`, `vinasse` taken care of by WWT
    # `filter_cake` taken care of by BT
    # 'Yeast': ('Yeast',), # no price considered, no GWP considered (probably used in fermentation)
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_sc1g_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    wwt_kwdct['skip_AeF'] = True
    sys_dct = {
        'load': {'name': 'S1', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'sugarcane_sys',
        'create_wastewater_process': wwt_kwdct,
        # `wastewater` is mixed from `fiber_fines` (taken care of),
        # `stripper_bottoms_product` (~20 mg/L COD), and `evaporator_condensate` (only water)
        'ww_streams': ('fiber_fines', 'vinasse'),
        'solids_streams': ('bagasse', 'filter_cake'), # `bagasse`, `filter_cake`
        'BT': 'BT401',
        'new_wwt_connections': {'solids': ('BT401', 0), 'biogas': ('BT401', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct)
    return exist_sys, new_sys


def simulate_sc1g_systems(**sys_kwdct):
    global exist_sys, new_sys
    exist_sys, new_sys = create_sc1g_comparison_systems(**sys_kwdct)
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_sc1g_comparison_models():
    exist_sys, new_sys = create_sc1g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'sugarcane',
        'FERM_product': info['FERM_product'],
        'PT_rx': 'R301',
        'fermentor': 'R301',
        'reactions': {
            'PT glucan-to-glucose': ('hydrolysis_reaction', ),
            'FERM glucose-to-product': ('fermentation_reaction', ),
            },
        'BT': 'BT401',
        'BT_eff': ('boiler_efficiency', 'turbogenerator_efficiency'),
        'wwt_system': 'exist_sys_wwt',
        'wwt_ID': info['WWT_ID'],
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['sludge'] = 'sludge'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_sc1g_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_sc1g_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_sc1g_systems(biodegradability=1)
    # exist_model, new_model = create_sc1g_comparison_models()
    exist_model, new_model = evaluate_sc1g_models( # 1G BMP should be high
        include_baseline=True,
        include_uncertainty=True,
        N_uncertainty=1000,
        # uncertainty_skip_exist=True,
        )
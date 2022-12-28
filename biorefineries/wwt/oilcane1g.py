#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery (configuration O1/1):
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
    'abbr': 'oc1g',
    'WWT_ID': '11',
    'is2G': False,
    'FERM_product': ['ethanol',],
    'add_BT': False,
    'ww_price': None,
    }

CF_dct = {
    ##### Feeds #####
    'catalyst': ('TEcatalyst',), # methanol and NaOCH3
    'denaturant': ('Denaturant'),
    'dryer_natural_gas': ('CH4',),
    'H3PO4': ('H3PO4',), # moisture content already adjusted
    'HCl': ('HCl',), # moisture content already adjusted
    'lime': ('CaO', 0.046), # CaO and water
    'methanol': ('Methanol',),
    'NaOH': ('NaOH',),
    'natural_gas': ('CH4',), # usually not needed
    'oilcane': ('Oilcane',), # moisture content already adjusted
    'polymer': ('Polymer'),
    'pure_glycerine': ('GlycerinPure',),
    ##### Co-products #####
    'biodiesel': ('Biodiesel',), # has <0.01 wt% impurities
    'crude_glycerol': ('GlycerinCrude',),
    # `fiber_fines`, `wastewater`, `vinasse` taken care of by WWT
    # `filter_cake` taken care of by BT
    # 'Yeast': ('Yeast',), # no price considered, no GWP considered (probably used in fermentation)
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_oc1g_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    wwt_kwdct['skip_AeF'] = True
    sys_dct = {
        'load': {'name': 'O1', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'BT': 'BT601',
        'create_wastewater_process': wwt_kwdct,
        'ww_streams': ('fiber_fines', 'wastewater', 'vinasse',),
        'solids_streams': ('filter_cake', ('M601', 0),),
        'new_wwt_connections': {'solids': ('BT601', 0), 'biogas': ('BT601', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct)
    return exist_sys, new_sys


def simulate_oc1g_systems(**sys_kwdct):
    global exist_sys, new_sys
    exist_sys, new_sys = create_oc1g_comparison_systems(**sys_kwdct)
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_oc1g_comparison_models():
    exist_sys, new_sys = create_oc1g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'oilcane',
        'FERM_product': info['FERM_product'],
        'PT_rx': 'R301',
        'fermentor': 'R301',
        'TE_rx': 'U501',
        'isplit_efficiency_is_reversed': False,
        'reactions': {
            'PT glucan-to-glucose': ('hydrolysis_reaction', ),
            'FERM glucose-to-product': ('fermentation_reaction', ),
            'FERM oil-to-FFA': ('oil_reaction', 0), # not fermentation, but happens in the fermentor
            'TE oil-to-product': ('transesterification', (0, 1, 2)),
            },
        'BT': 'BT601',
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


def evaluate_oc1g_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_oc1g_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_oc1g_systems(biodegradability=1)
    # exist_model, new_model = create_oc1g_comparison_models()
    exist_model, new_model = evaluate_oc1g_models( # 1G BMP should be high
        include_baseline=True,
        include_uncertainty=True,
        N_uncertainty=1000,
        # uncertainty_skip_exist=True,
        )
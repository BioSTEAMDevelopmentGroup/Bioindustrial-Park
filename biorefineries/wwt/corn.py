#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the corn biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/corn
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries import corn as cn
from biorefineries.wwt import (
    create_comparison_systems, simulate_systems,
    create_comparison_models, evaluate_models,
    )

info = {
    'abbr': 'cn',
    'WWT_ID': '7',
    'is2G': False,
    'FERM_product': ['ethanol',],
    'add_BT': True,
    'ww_price': None,
    }

CF_dct = {
    ##### Feeds #####
    'alpha_amylase': ('AlphaAmylase',), # 0.00082 soluble protein and water
    'ammonia': ('NH3',),
    'corn': ('Corn',), # adjust for the moisture content
    'denaturant': ('Denaturant',),
    'gluco_amylase': ('GlucoAmylase',), # 0.0011 soluble protein and water
    'lime': ('CaO',),
    'natural_gas': ('CH4',),
    'steam':('Steam',),
    'sulfuric_acid': ('H2SO4',),
    ('X611', 'ins', 2): ('CH4',),
    'yeast': ('Yeast',),
    ##### Co-products #####
    'crude_oil': ('CornOil',), # triolein
    'DDGS': ('DDGS',),
     # `s4` (from MH103) taken care of by WWT
}


# %%

# =============================================================================
# Systems
# =============================================================================

def create_cn_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    wwt_kwdct['skip_AeF'] = True
    sys_dct = {
        'system_name': 'corn_sys',
        'BT': 'CHP',
        # 'BT': 'BT',
        'create_wastewater_process': wwt_kwdct,
        'ww_streams': (('MH103', 1), ('MX5', 0)),
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, cn, sys_dct)
    return exist_sys, new_sys


def simulate_cn_systems(**sys_kwdct):
    global exist_sys, new_sys
    exist_sys, new_sys = create_cn_comparison_systems(**sys_kwdct)
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_cn_comparison_models():
    exist_sys, new_sys = create_cn_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'corn',
        'FERM_product': info['FERM_product'],
        'PT_rx': 'V310',
        'fermentor': 'V405',
        'reactions': {
            'PT glucan-to-glucose': ('reaction',),
            'FERM glucose-to-product': ('reaction',),
            },
        'wwt_system': 'exist_sys_wwt',
        'wwt_ID': info['WWT_ID'],
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['BT'] = 'CHP'
    # new_model_dct['BT'] = 'BT'
    new_model_dct['BT_eff'] = ('eff',) # need to be an Iterable
    # new_model_dct['BT_eff'] = ('boiler_efficiency', 'turbogenerator_efficiency')
    new_model_dct['sludge'] = 'sludge'
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = create_comparison_models(new_sys, new_model_dct)

    return exist_model, new_model


def evaluate_cn_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_cn_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_cn_systems(biodegradability=1)
    # exist_model, new_model = create_cn_comparison_models()
    exist_model, new_model = evaluate_cn_models( # 1G BMP should be high
        include_baseline=True,
        include_uncertainty=True,
        N_uncertainty=1000,
        # uncertainty_skip_exist=True,
        )
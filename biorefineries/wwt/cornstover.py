#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from biorefineries import cornstover as cs
from biorefineries.wwt import (
    create_comparison_systems, simulate_systems,
    create_comparison_models,evaluate_models,
    )

info = {
    'abbr': 'cs',
    'WWT_ID': '6',
    'is2G': True,
    'FERM_product': ['ethanol',],
    'add_BT': False,
    'ww_price': None,
    }

CF_dct = { # all streams are feeds
    'ammonia': ('NH4OH',), # NH4OH
    'caustic': ('NaOH', 0.5), # NaOH and water
    'cellulase': ('Cellulase', 0.05), # cellulase and water
    'cornstover': ('CornStover',), # adjust for the moisture content
    'CSL': ('CSL',),
    'DAP': ('DAP',),
    'denaturant': ('Denaturant'),
    'FGD_lime': ('Lime', 0.4513), # lime and water
    'natural_gas': ('CH4',), # CH4 actually not used due to heat surplus
    'sulfuric_acid': ('H2SO4',),
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_cs_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    sys_dct = {
        'system_name': 'cornstover_sys',
        'create_wastewater_process': wwt_kwdct,
        'BT': 'BT',
        'new_wwt_connections': {'sludge': ('slurry_mixer', 0), 'biogas': ('gas_mixer', 0)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, cs, sys_dct)
    return exist_sys, new_sys


def simulate_cs_systems(**sys_kwdct):
    global exist_sys, new_sys
    exist_sys, new_sys = create_cs_comparison_systems(**sys_kwdct)
    # If using conservative biodegradability,
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_cs_comparison_models():
    exist_sys, new_sys = create_cs_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'cornstover',
        'FERM_product': info['FERM_product'],
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_solids_mixer': 'M203',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R303',
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification', 2),
            'FERM glucose-to-product': ('cofermentation', 0),
            'FERM xylose-to-product': ('cofermentation', 4),
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
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_cs_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_cs_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_cs_systems(biodegradability=1)
    # exist_model, new_model = create_cs_comparison_models()
    exist_model, new_model = evaluate_cs_models(
        # include_baseline=True,
        # include_uncertainty=True,
        include_BMP=True,
        # N_uncertainty=1000,
        # uncertainty_skip_exist=True, # for testing
        N_BMP=100,
        # BMPs=(0.5, 0.9499,), # for testing, 0.9499 allows for minor error
        )
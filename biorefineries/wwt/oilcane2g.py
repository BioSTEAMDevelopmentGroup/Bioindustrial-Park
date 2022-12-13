#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery (configuration O2/2):
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
    'abbr': 'oc2g',
    'WWT_ID': '5',
    'is2G': True,
    'FERM_product': ['advanced_ethanol', 'cellulosic_ethanol',],
    'add_BT': False,
    'ww_price': None,
    }

CF_dct = {
    ##### Feeds #####
    'catalyst': ('TEcatalyst',), # methanol and NaOCH3
    'caustic': ('NaOH', 0.5), # NaOH and water
    'cellulase': ('Cellulase', 0.05), # cellulase and water
    'CSL': ('CSL',),
    'DAP': ('DAP',),
    'denaturant': ('Denaturant'),
    'dryer_natural_gas': ('CH4',),
    'FGD_lime': ('Lime', 0.4513), # lime and water
    'H3PO4': ('H3PO4',), # moisture content already adjusted
    'HCl': ('HCl',), # moisture content already adjusted
    'lime': ('CaO', 0.046), # CaO and water
    'methanol': ('Methanol',),
    'NaOH': ('NaOH',),
    'natural_gas': ('CH4',),
    'oilcane': ('Oilcane',), # moisture content already adjusted
    'polymer': ('Polymer',),
    'pure_glycerine': ('GlycerinPure',),
    'urea': ('Urea',),
    ##### Co-products #####
    'biodiesel': ('Biodiesel',), # has <0.01 wt% impurities
    'crude_glycerol': ('GlycerinCrude',),
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_oc2g_comparison_systems(biodegradability=1): # will be multiplied by 0.86/0.05 for biogas/cell mass
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': biodegradability,})
    sys_dct = {
        'load': {'name': 'O2', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'create_wastewater_process': wwt_kwdct,
        'BT': 'BT701',
        'new_wwt_connections': {'sludge': ('M701', 0), 'biogas': ('BT701', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct)
    return exist_sys, new_sys


def simulate_oc2g_systems(**sys_kwdct):
    global exist_sys, new_sys
    exist_sys, new_sys = create_oc2g_comparison_systems(**sys_kwdct)
    # If using conservative biodegradability,
    # ~184 mg/L COD, mostly (~150/>80%) due to soluble lignin and arabinose
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_oc2g_comparison_models():
    exist_sys, new_sys = create_oc2g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'CF_dct': CF_dct,
        'feedstock': 'oilcane',
        'FERM_product': info['FERM_product'],
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_solids_mixer': 'M301',
        'PT_rx': 'R301',
        'EH_mixer': 'M401',
        'EH_rx': 'U401',
        'fermentor': 'R401',
        'TE_rx': 'U802',
        'isplit_efficiency_is_reversed': True,
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification', 2),
            'FERM glucose-to-product': ('cofermentation', 0),
            'FERM xylose-to-product': ('cofermentation', 1),
            'FERM oil-to-FFA': ('oil_reaction', 0), # not fermentation, but happens in the fermentor
            'TE oil-to-product': ('transesterification', (0, 1, 2)),
            },
        'BT': 'BT701',
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


def evaluate_oc2g_models(**eval_kwdct):
    global exist_model, new_model
    exist_model, new_model = create_oc2g_comparison_models()
    evaluate_models(exist_model, new_model, info['abbr'], **eval_kwdct)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_oc2g_systems(biodegradability=1)
    # exist_model, new_model = create_oc2g_comparison_models()
    exist_model, new_model = evaluate_oc2g_models(
        # include_baseline=True,
        # include_uncertainty=True,
        include_BMP=True,
        # N_uncertainty=1000,
        # uncertainty_skip_exist=True, # for testing
        N_BMP=100,
        # BMPs=(0.5, 0.9499,), # for testing, 0.9499 allows for minor error
        )
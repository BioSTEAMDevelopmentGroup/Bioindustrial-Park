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
    'FERM_product': 'ethanol',
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_sc1g_comparison_systems(default_BD=True):
    from biorefineries.wwt import create_comparison_systems
    from biorefineries import oilcane as oc
    BD = {} if not default_BD else 1.
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': BD,})
    wwt_kwdct['skip_AeF'] = True
    CF_dct = {
        ##### Feeds #####
        'denaturant': ('Denaturant'),
        'dryer_natural_gas': ('CH4',),
        'H3PO4': ('H3PO4',),
        'lime': ('CaO', 0.046), # CaO and water
        'polymer': ('Polymer'),
        'sugarcane': ('Sugarcane',), # moisture content already adjusted
        ##### Co-products #####
        # 'Yeast': ('Yeast',), # no price considered, no GWP considered (probably used in fermentation)
        # `fiber_fines`, `wastewater`, `vinasse` taken care of by WWT
        # `filter_cake` taken care of by BT
        # `s41` (from `T302`) is empty
        }
    sys_dct = {
        'load': {'name': 'S1', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'create_wastewater_process': wwt_kwdct,
        # `fiber_fines`, `vinasse`
        # `wastewater` is mixed from `fiber_fines` (taken care of),
        # `stripper_bottoms_product` (~20 mg/L COD), and `evaporator_condensate` (only water)
        'ww_streams': (('U211', 1), ('H302', 1)),
        'solids_streams': (('U207', 0), ('U210', 0)), # `bagasse`, `filter_cake`
        'BT': 'BT401',
        'new_wwt_connections': {'solids': ('BT401', 0), 'biogas': ('BT401', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct)
    return exist_sys, new_sys


def simulate_sc1g_systems(**sys_kwdct):
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_sc1g_comparison_systems(**sys_kwdct)
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
        'FERM_product': info['FERM_product'],
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


def evaluate_sc1g_models(**eval_kwdct):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_sc1g_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **eval_kwdct)



# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    exist_sys, new_sys = simulate_sc1g_systems(default_BD=True)
    # exist_model, new_model = create_sc1g_comparison_models()
    # exist_model, new_model = evaluate_sc1g_models(N=10)
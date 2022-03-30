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
    'abbr': 'oc1g',
    'WWT_ID': '11',
    'is2G': False,
    'FERM_product': 'ethanol',
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_oc1g_comparison_systems(default_BD=True):
    from biorefineries.wwt import create_comparison_systems
    from biorefineries import oilcane as oc
    BD = {} if not default_BD else 1.
    wwt_kwdct = dict.fromkeys(('IC_kwargs', 'AnMBR_kwargs',), {'biodegradability': BD,})
    wwt_kwdct['skip_AeF'] = True
    CF_dct = {
        ##### Feeds #####
        'catalyst': ('TEcatalyst',), # methanol and NaOCH3
        'denaturant': ('Denaturant'),
        'dryer_natural_gas': ('CH4',),
        'H3PO4': ('H3PO4', 0.5), # H3PO4 and water
        'HCl': ('HCl', 0.3498), # HCl and water
        'lime': ('CaO', 0.046), # CaO and water
        'methanol': ('Methanol',),
        'NaOH': ('NaOH',),
        'natural_gas': ('CH4',), # probably not needed
        'oilcane': ('Oilcane', (1-0.7)), # adjust for the moisture content
        'polymer': ('Flocculant'),
        'pure_glycerine': ('Glycerol',),
        ##### Co-products #####
        'biodiesel': ('Biodiesel',), # has <0.01 wt% impurities
        'crude_glycerol': ('Glycerol', 0.8), # glycerol and water #!!! need to see if there's a CF for crude glycerol
        'Yeast': ('Yeast',), #!!! not pure yeast, need to consider
        # `fiber_fines`, `wastewater`, `vinasse` taken care of by WWT
        # `filter_cake` taken care of by BT
        # `s46` (from `T302`) is empty
        }
    sys_dct = {
        'load': {'name': 'O1', 'cache': None, 'reduce_chemicals': False},
        'system_name': 'oilcane_sys',
        'BT': 'BT701',
        'create_wastewater_process': wwt_kwdct,
        # `fiber_fines`, `wastewater`, `vinasse`
        'ww_streams': (('U206', 1), ('M603', 0), ('M403', 0),),
        'solids_streams': (('U205', 0), ('M701', 0),), # the first one is `filter_cake`
        'BT': 'BT701',
        'new_wwt_connections': {'solids': ('BT701', 0), 'biogas': ('BT701', 1)},
        'CF_dct': CF_dct,
        }
    exist_sys, new_sys = create_comparison_systems(info, oc, sys_dct)
    return exist_sys, new_sys


def simulate_oc1g_systems(**sys_kwdct):
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_oc1g_comparison_systems(**sys_kwdct)
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys



# %%

# =============================================================================
# Models
# =============================================================================

def create_oc1g_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_oc1g_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'oilcane',
        'FERM_product': info['FERM_product'],
        'PT_rx': 'R301',
        'fermentor': 'R301',
        'TE_rx': 'U601',
        'isplit_efficiency_is_reversed': False,
        'bagasse_oil_extraction': 'U406',
        'bagasse_oil_retention': 'U201',
        'reactions': {
            'PT glucan-to-glucose': ('hydrolysis_reaction', ),
            'FERM glucan-to-product': ('fermentation_reaction', ),
            'FERM oil-to-FFA': ('oil_reaction', 0), # not fermentation, but happens in the fermentor
            'TE oil-to-product': ('transesterification', (0, 1, 2)),
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


def evaluate_oc1g_models(**eval_kwdct):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_oc1g_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **eval_kwdct)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_oc1g_systems(default_BD=True)
    # exist_model, new_model = create_oc1g_comparison_models()
    exist_model, new_model = evaluate_oc1g_models(N=1000)
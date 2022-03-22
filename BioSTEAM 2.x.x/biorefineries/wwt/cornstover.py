#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

info = {
    'abbr': 'cs',
    'WWT_ID': '6',
    'is2G': True,
    'add_CHP': False,
    'ww_price': None,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_cs_comparison_systems():
    # # Create from scratch, IRR doesn't match as closely as the method below
    # from biorefineries.wwt import create_comparison_systems
    # from biorefineries import cornstover as cs
    # from biorefineries.cornstover import (
    #     create_chemicals,
    #     create_system,
    #     create_tea,
    #     load_process_settings,
    #     )
    # OSBL_IDs = [u.ID for u in cs.cornstover_tea.OSBL_units]
    # functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    # sys_dct = {
    #     'create_system': {'include_blowdown_recycle': True},
    #     'BT': 'BT',
    #     'new_wwt_connections': {'sludge': ('M501', 0), 'biogas': ('BT', 1)},
    #     }
    # exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)

    # exist_f, new_f = exist_sys.flowsheet, new_sys.flowsheet
    # exist_sys.TEA.OSBL_units = [getattr(exist_f.unit, ID) for ID in OSBL_IDs]
    # OSBL_IDs.remove('WWTC')
    # OSBL_IDs.extend([u.ID for u in new_f.system.new_sys_wwt.units])
    # new_sys.TEA.OSBL_units = [getattr(new_f.unit, ID) for ID in OSBL_IDs]

    from biorefineries.wwt import create_comparison_systems
    from biorefineries import cornstover as cs
    sys_dct = {
        'system_name': 'cornstover_sys',
        'BT': 'BT',
        'new_wwt_connections': {'sludge': ('M501', 0), 'biogas': ('BT', 1)},
        }
    exist_sys, new_sys = create_comparison_systems(info, cs, sys_dct, from_load=True)

    return exist_sys, new_sys


def simulate_cs_systems():
    from biorefineries.wwt import simulate_systems
    global exist_sys, new_sys
    exist_sys, new_sys = create_cs_comparison_systems()
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_cs_comparison_models():
    from biorefineries.wwt import create_comparison_models
    exist_sys, new_sys = create_cs_comparison_systems()

    ##### Existing system #####
    exist_model_dct = {
        'abbr': info['abbr'],
        'feedstock': 'cornstover',
        'primary_product': 'ethanol',
        'sulfuric_acid': 'sulfuric_acid',
        'acid_dilution_water': 'warm_process_water_1',
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_acid_mixer': 'M201',
        'PT_solids_mixer': 'M203',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R303',
        'reactions': {
            'PT glucan-to-glucose': ('reactions', 0),
            'PT xylan-to-xylose': ('reactions', 8),
            'EH glucan-to-glucose': ('saccharification', 2),
            'FERM glucan-to-product': ('cofermentation', 0),
            'FERM xylan-to-product': ('cofermentation', 4),
            },
        'BT': 'BT',
        'BT_eff': ('boiler_efficiency', 'turbogenerator_efficiency'),
        'wwt_system': 'exist_sys_wwt',
        'is2G': info['is2G'],
        }
    exist_model = create_comparison_models(exist_sys, exist_model_dct)

    ##### With the new wastewater treatment process #####
    new_model_dct = exist_model_dct.copy()
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model_dct['new_wwt_ID'] = info['WWT_ID']
    new_model = create_comparison_models(new_sys, new_model_dct)
    return exist_model, new_model


def evaluate_cs_models(**kwargs):
    from biorefineries.wwt import evaluate_models
    global exist_model, new_model
    exist_model, new_model = create_cs_comparison_models()
    return evaluate_models(exist_model, new_model, abbr=info['abbr'], **kwargs)


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    # exist_sys, new_sys = simulate_cs_systems()
    # exist_model, new_model = create_cs_comparison_models()
    exist_model, new_model = evaluate_cs_models(N=10)
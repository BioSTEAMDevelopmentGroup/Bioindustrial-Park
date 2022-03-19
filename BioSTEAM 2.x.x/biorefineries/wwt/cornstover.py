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
    exist_sys, new_sys = create_cs_comparison_systems()
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================

def create_cs_comparison_models():
    from biosteam import Model
    from biorefineries.wwt import add_2G_parameters, add_new_wwt_parameters, add_metrics
    exist_sys, new_sys = create_cs_comparison_systems()

    ##### Existing system #####
    exist_model = Model(exist_sys)
    exist_model_dct = {
        'feedstock': 'cornstover',
        'sulfuric_acid': 'sulfuric_acid',
        'acid_dilution_water': 'warm_process_water_1',
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_acid_mixer': 'M201',
        'PT_solids_mixer': 'M203',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R303',
        'BT': 'BT',
        'wwt_system': 'exist_sys_wwt',
        }

    exist_model = add_2G_parameters(exist_model, exist_model_dct)
    BT = exist_sys.flowsheet.unit.BT
    eff = BT.boiler_efficiency * BT.turbogenerator_efficiency
    exist_model = add_metrics(exist_model, exist_model_dct, eff=eff)

    ##### With the new wastewater treatment process #####
    new_model = Model(new_sys)
    new_model_dct = exist_model_dct.copy()
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = add_2G_parameters(new_model, new_model_dct)
    new_model = add_new_wwt_parameters(new_model, process_ID=info['WWT_ID'])
    new_model = add_metrics(new_model, new_model_dct, eff=eff)
    return exist_model, new_model


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    exist_sys, new_sys = simulate_cs_systems()
    # exist_model, new_model = create_cs_comparison_models()
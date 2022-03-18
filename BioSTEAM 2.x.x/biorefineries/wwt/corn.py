#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the corn biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/corn
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


info = {
    'abbr': 'cn',
    'WWT_ID': '7',
    'is2G': False,
    'add_CHP': True,
    'ww_price': -0.003,
    }


# %%

# =============================================================================
# Systems
# =============================================================================

def create_cn_comparison_systems():
    import biosteam as bst
    from biorefineries.wwt import create_comparison_systems
    from biorefineries.corn import (
        create_chemicals,
        create_system,
        create_tea,
        load_process_settings
        )
    functions = (create_chemicals, create_system, create_tea, load_process_settings,)
    sys_dct = {
        'create_system': {'flowsheet': bst.main_flowsheet},
        'create_wastewater_process': {'skip_AeF': True},
        'ww_streams': (('MH103', 1), ('MX5', 0)),
        }
    exist_sys, new_sys = create_comparison_systems(info, functions, sys_dct)
    return exist_sys, new_sys


def simulate_cn_systems():
    from biorefineries.wwt import simulate_systems
    exist_sys, new_sys = create_cn_comparison_systems()
    simulate_systems(exist_sys, new_sys, info)
    return exist_sys, new_sys


# %%

# =============================================================================
# Models
# =============================================================================


# %%

# =============================================================================
# Run
# =============================================================================

if __name__ == '__main__':
    exist_sys, new_sys = simulate_cn_systems()
    # exist_model, new_model = create_cn_comparison_models()
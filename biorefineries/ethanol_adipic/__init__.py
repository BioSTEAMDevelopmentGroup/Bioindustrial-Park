#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%


import biosteam as bst
bst.speed_up()

from . import (
    _chemicals,
    _utils,
    _settings,
    _units,
    _facilities,
    _tea,
    _processes,
    systems
    )


from ._chemicals import chems
from ._utils import auom, CEPCI
from ._settings import price, CFs

getattr = getattr

def load_system(system_kind='acid', depot_kind='HMPP'):
    if not system_kind in ('acid', 'base', 'AFEX'): \
        raise ValueError('system_kind can only be "acid", "base", or "AFEX", '
                         f'not {system_kind}.')
    if not depot_kind in ('CPP', 'CPP_AFEX', 'HMPP', 'HMPP_AFEX'): \
        raise ValueError('depot_kind can only be "CPP", "CPP_AFEX", "HMPP", '
                         f'or "HMPP_AFEX", not {depot_kind}.')

    global flowsheet, groups, teas, funcs, biorefinery, tea

    depot_dct = systems.depot_dct[depot_kind]
    create_sys = getattr(systems, f'create_{system_kind}_biorefinery')

    flowsheet, groups, teas, funcs = create_sys(depot_dct['preprocessed'])
    biorefinery = flowsheet.system.biorefinery
    tea = teas['tea']

    bst.settings.set_thermo(chems)
    bst.main_flowsheet.set_flowsheet(flowsheet)


load_system('acid', 'HMPP')


# %%

# =============================================================================
# Simulate system and get results
# =============================================================================

feedstock_GWPs = systems.feedstock_GWPs

def simulate_and_print(depot_for_GWP=None):
    system_kind = flowsheet.ID
    s = flowsheet.stream
    line = system_kind if system_kind.isupper() else system_kind.capitalize()

    print(f'\n---------- {line} Biorefinery ----------')
    print(f'MESP: ${funcs["simulate_get_MESP"]()*systems._ethanol_kg_2_gal:.2f}/gal')
    print(f'GWP: {funcs["get_GWP"]():.3f} kg CO2-eq/gal ethanol without feedstock')
    if depot_for_GWP:
        GWP = funcs['get_GWP']()
        GWP = funcs['get_GWP']() + (feedstock_GWPs[depot_for_GWP]*s.feedstock.F_mass) \
            / (s.ethanol.F_mass/systems._ethanol_kg_2_gal)
        print(f'GWP: {GWP:.3f} kg CO2-eq/gal ethanol with feedstock')
    print('--------------------------------------')
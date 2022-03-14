#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <zoe.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

from biosteam import main_flowsheet
from . import (
    _chemicals,
    utils,
    _settings,
    _units,
    _facilities,
    _tea,
    _processes,
    systems,
    )

from ._chemicals import chems
from .utils import auom, CEPCI
from ._settings import price, CFs
from ._processes import update_settings
from .systems import *
from .models import *

def load_system(kind='SSCF'):
    if not kind in ('SSCF', 'SHF'):
        raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')
    global flowsheet, groups, teas, funcs, \
        lactic_sys, lactic_tea, feedstock, lactic_acid, \
        simulate_and_print, simulate_fermentation_improvement, \
        simulate_separation_improvement, simulate_operating_improvement

    flowsheet = getattr(systems, f'{kind}_flowsheet')
    groups = getattr(systems, f'{kind}_groups')
    teas = getattr(systems, f'{kind}_teas')
    funcs = getattr(systems, f'{kind}_funcs')

    update_settings(chems)
    main_flowsheet.set_flowsheet(flowsheet)

    lactic_sys = flowsheet.system.lactic_sys
    lactic_tea = teas['lactic_tea']
    feedstock = flowsheet.stream.feedstock
    lactic_acid = flowsheet.stream.lactic_acid

    simulate_and_print = lambda : systems.simulate_and_print(kind)
    simulate_fermentation_improvement = \
        lambda: systems.simulate_fermentation_improvement(kind)
    simulate_separation_improvement = \
        lambda: systems.simulate_separation_improvement(kind)
    simulate_operating_improvement = \
        lambda: systems.simulate_operating_improvement(kind)
    return

load_system('SSCF')

# __all__ = (
#     flowsheet, groups, teas, funcs, \
#     lactic_sys, lactic_tea, feedstock, lactic_acid, \
#     simulate_and_print, simulate_fermentation_improvement, \
#     simulate_separation_improvement, simulate_operating_improvement,
#     *lactic_sys.units, # not seem to be working
#     )
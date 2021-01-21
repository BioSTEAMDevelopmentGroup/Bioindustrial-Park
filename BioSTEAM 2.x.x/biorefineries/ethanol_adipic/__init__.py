#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

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

# getattr = getattr

# def load_system(kind='SSCF'):
#     global flowsheet, groups, teas, funcs, lactic_sys, lactic_tea, \
#         simulate_and_print, simulate_fermentation_improvement, \
#         simulate_separation_improvement, simulate_operating_improvement
    
#     flowsheet = getattr(systems, f'{kind}_flowsheet')
#     groups = getattr(systems, f'{kind}_groups')
#     teas = getattr(systems, f'{kind}_teas')
#     funcs = getattr(systems, f'{kind}_funcs')
    
#     update_settings(chems)
#     bst.main_flowsheet.set_flowsheet(flowsheet)
        
#     lactic_sys = flowsheet.system.lactic_sys
#     lactic_tea = teas['lactic_tea']
#     simulate_and_print = lambda : systems.simulate_and_print(kind)
#     simulate_fermentation_improvement = \
#         lambda: systems.simulate_fermentation_improvement(kind)
#     simulate_separation_improvement = \
#         lambda: systems.simulate_separation_improvement(kind)
#     simulate_operating_improvement = \
#         lambda: systems.simulate_operating_improvement(kind)


# load_system('SSCF')
# simulate_and_print()
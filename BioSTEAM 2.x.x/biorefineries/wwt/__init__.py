#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os
path = os.path.dirname(__file__)
results_path = os.path.join(path, 'results')
del os

import biosteam as bst
from biorefineries import(
    sugarcane as sc,
    lipidcane as lc,
    cornstover as cs,
    )
ethanol_density_kggal = cs.ethanol_density_kggal

from biorefineries.sugarcane import chemicals as sc_chems
from biorefineries.lipidcane import chemicals as lc_chems
from biorefineries.cornstover import chemicals as cs_chems

bst.speed_up()

# from . import (
#     _chemicals,
#     utils,
#     _settings,
#     _internal_circulation_rx,
#     _wwt_pump,
#     _polishing_filter,
#     _membrane_bioreactor,
#     _sludge_handling,
#     _wwt_sys,
#     _lca,
#     )

# from ._chemicals import *
# from ._settings import *
# from ._internal_circulation_rx import *
# from ._wwt_pump import *
# from ._polishing_filter import *
# from ._membrane_bioreactor import *
# from ._sludge_handling import *
# from ._wwt_sys import *
# from ._lca import *

import \
    _chemicals, \
    utils, \
    _settings, \
    _internal_circulation_rx, \
    _wwt_pump, \
    _polishing_filter, \
    _membrane_bioreactor, \
    _sludge_handling, \
    _wwt_sys, \
    _lca


from _chemicals import *
from _settings import *
from _internal_circulation_rx import *
from _wwt_pump import *
from _polishing_filter import *
from _membrane_bioreactor import *
from _sludge_handling import *
from _wwt_sys import *
from _lca import *

__all__ = (
    'results_path', 'ethanol_density_kggal',
    'sc', 'lc', 'cs', 'sc_chems', 'lc_chems', 'cs_chems',
    *_chemicals.__all__,
    *_settings.__all__,
    *_internal_circulation_rx.__all__,
    *_wwt_pump.__all__,
    *_polishing_filter.__all__,
    *_membrane_bioreactor.__all__,
    *_sludge_handling.__all__,
    *_wwt_sys.__all__,
    *_lca.__all__,
    )
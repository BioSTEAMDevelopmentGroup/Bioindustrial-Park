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

import thermosteam as tmo, biosteam as bst
auom = tmo.units_of_measure.AbsoluteUnitsOfMeasure

# Using BioSTEAM's default values
# (about the same as the previous dict, only as integers)
CEPCI = bst.design_tools.CEPCI_by_year
# # Legact Chemical Engineering Plant Cost Index from Chemical Engineering Magazine
# # (https://www.chemengonline.com/the-magazine/)
# CEPCI = {1997: 386.5,
#          1998: 389.5,
#          2007: 525.4,
#          2009: 521.9,
#          2010: 550.8,
#          2011: 585.7,
#          2012: 584.6,
#          2013: 567.3,
#          2014: 576.1,
#          2016: 541.7}


# %%

from ._chemicals import *
global chemicals
chemicals = create_chemicals()
from .utils import *
from ._process_settings import *
from ._tea import *
from ._processes import *
from .systems import *
from .models import *

from . import (
    _chemicals,
    utils,
    _process_settings,
    _units,
    _facilities,
    _tea,
    _processes,
    systems,
    models,
    )


_system_loaded = False
def load(kind='SSCF'):
    _load_system(kind)
    dct = globals()
    dct.update(flowsheet.to_dict())

def _load_system(kind='SSCF'):
    if not kind in ('SSCF', 'SHF'):
        raise ValueError(f'kind can only be "SSCF" or "SHF", not "{kind}".')
    global flowsheet, groups, teas, funcs, lactic_sys, lactic_tea
    load_process_settings(chems)
    flowsheet, groups, teas, funcs = create_system(kind, return_all=True)
    lactic_sys = flowsheet.system.lactic_sys
    lactic_tea = teas['lactic_tea']
    return

from .. import PY37
if PY37:
    def __getattr__(name):
        if name == 'chemicals': return chemicals
        if not _system_loaded:
            _load_system()
            dct = globals()
            dct.update(flowsheet.to_dict())
            if name in dct: return dct[name]
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
else:
    load()
del PY37

__all__ = (
    'auom', 'CEPCI',
    'flowsheet', 'groups', 'teas', 'funcs',
    'lactic_sys', 'lactic_tea',
    *_chemicals.__all__,
    *_process_settings.__all__,
    *systems.__all__,
    *models.__all__,
    )
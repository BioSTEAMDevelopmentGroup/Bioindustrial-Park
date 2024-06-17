# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import (units,
               _process_settings,
               _chemicals,
               systems,
               _tea,
)

__all__ = [*_process_settings.__all__,
           *_chemicals.__all__,
           *systems.__all__,
           *_tea.__all__,
           'cornstover_sys',
           'cornstover_tea', 
           'flowsheet',
           'Area100',
           'Area200',
           'Area300',
           'Area400',
           'Area500',
           'Area600',
           'Area700',
           'Area800',
           'AllAreas',
           'areas',
]

from .units import *
from ._process_settings import *
from ._chemicals import *
from .systems import *
from ._tea import *
from biorefineries.cellulosic import Biorefinery

_biorefinery_loaded = False

def load():
    global _biorefinery_loaded
    br = Biorefinery('corn stover ethanol')
    dct = globals()
    dct.update(br.__dict__)
    dct['chemicals'] = br.chemicals
    _biorefinery_loaded = True

def __getattr__(name):
    if not _biorefinery_loaded:
        load()
        dct = globals()
        if name in dct: return dct[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
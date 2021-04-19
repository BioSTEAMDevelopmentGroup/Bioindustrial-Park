#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

@author: Sarang Bhagwat and Yoel Cortes-Pena
"""
from biorefineries import PY37
from . import (
    process_settings, 
    chemicals_data, 
    tea, 
    units, 
    facilities
)
from .chemicals_data import *
from .process_settings import *
from .tea import *
from .units import *
from .facilities import *

__all__ = [
    'process_settings', 
    'chemicals_data', 
    'tea', 
    'units', 
    'facilities',
]

_system_loaded = False
_chemicals_loaded = False

default_configuration = 'cellulosic'

def load(configuration=None):
    if not _chemicals_loaded: _load_chemicals()
    _load_system(configuration)
    dct = globals()
    dct.update(flowsheet.system.__dict__)
    dct.update(flowsheet.stream.__dict__)
    dct.update(flowsheet.unit.__dict__)

def _load_system(configuration=None):
    load_process_settings()
    if not configuration: configuration = default_configuration
    if configuration == 'cellulosic':
        _load_celluloic_system()
    elif configuration == 'sugarcane':
        _load_sugarcane_system()
    else:
        raise ValueError("configuration must be either 'cellulosic' or 'sugarcane'; "
                        f"not '{configuration}'")

def _load_chemicals():
    global chemicals
    from .chemicals_data import HP_chemicals
    chemicals = HP_chemicals
    _chemicals_loaded = True

def _load_celluloic_system():
    global HP_sys, HP_tea, flowsheet, _system_loaded
    from .system import HP_sys, HP_tea, flowsheet
    _system_loaded = True

def _load_sugarcane_system():
    global HP_sys, HP_tea, flowsheet, _system_loaded
    from .system_sugarcane import HP_sys, HP_tea, flowsheet
    _system_loaded = True

if PY37:    
    def __getattr__(name):
        if not _chemicals_loaded:
            _load_chemicals()
            if name == 'chemicals': return chemicals
        if not _system_loaded: 
            try:
                _load_system()
            finally:
                dct = globals()
                dct.update(flowsheet.system.__dict__)
                dct.update(flowsheet.stream.__dict__)
                dct.update(flowsheet.unit.__dict__)
            if name in dct: return dct[name]
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
else:
    try:
        from lazypkg import LazyPkg
    except:
        from warnings import warn
        warn('Python 3.7 or newer is required to lazy load biorefinery; import '
             'and run the load function to load a biorefinery')
    else:
        LazyPkg(__name__, ['system'])
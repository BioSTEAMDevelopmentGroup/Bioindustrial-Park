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
               _system,
               _tea,
)

__all__ = [*units.__all__,
           *_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           'wheatstraw_sys',
           'wheatstraw_tea', 
           'flowsheet',
]

from .units import *
from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *

_system_loaded = False
_chemicals_loaded = False

def load():
    if not _chemicals_loaded: _load_chemicals()
    _load_system()
    dct = globals()
    dct.update(flowsheet.system.__dict__)
    dct.update(flowsheet.stream.__dict__)
    dct.update(flowsheet.unit.__dict__)

def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def _load_system():
    import biosteam as bst
    from biosteam import main_flowsheet as F
    global wheatstraw_sys, wheatstraw_tea, specs, flowsheet, _system_loaded
    flowsheet = bst.Flowsheet('wheatstraw')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    wheatstraw_sys = create_system()
    wheatstraw_sys.simulate()
    u = F.unit
    OSBL_units = (u.WWTC, u.CWP, u.CT, u.PWC, u.ADP, u.BT)
    wheatstraw_tea = create_tea(wheatstraw_sys, OSBL_units)
    ethanol = F.stream.ethanol
    ethanol.price = wheatstraw_tea.solve_price(ethanol)
    _system_loaded = True
    
def __getattr__(name):
    if not _chemicals_loaded:
        _load_chemicals()
        if name == 'chemicals': return chemicals
    if not _system_loaded: 
        _load_system()
        dct = globals()
        dct.update(flowsheet.system.__dict__)
        dct.update(flowsheet.stream.__dict__)
        dct.update(flowsheet.unit.__dict__)
        if name in dct: return dct[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from .. import PY37
from . import (utils,
               units,
               _process_settings,
               _chemicals,
               _system,
               _tea,
)

__all__ = [*utils.__all__,
           *units.__all__,
           *_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           'corn_sys',
           'corn_tea', 
           'flowsheet',
]

from .abm import *
from .utils import *
from .units import *
from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *

_system_loaded = False
_chemicals_loaded = False

def load():
    if not _chemicals_loaded: _load_chemicals()
    try:
        _load_system()
    finally:
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
    global corn_sys, corn_tea, flowsheet, all_areas, _system_loaded
    flowsheet = bst.Flowsheet('corn')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    corn_sys = create_system()
    corn_sys.simulate()
    corn_tea = create_tea(corn_sys)
    corn_tea.IRR = corn_tea.solve_IRR()
    all_areas = bst.process_tools.UnitGroup('All Areas', corn_sys.units)
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
    load()

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
from . import (_process_settings,
               _chemicals,
               _system,
               _tea,
)

__all__ = [*_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           'sugarcane_sys',
           'sugarcane_tea', 
           'flowsheet',
]

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
    global sugarcane_sys, sugarcane_tea, specs, flowsheet, _system_loaded
    flowsheet = bst.Flowsheet('sugarcane')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    sugarcane_sys = create_system()
    sugarcane_sys.simulate()
    sugarcane_tea = create_tea(sugarcane_sys)
    sugarcane_tea.IRR = sugarcane_tea.solve_IRR()
    _system_loaded = True
  
if PY37:
    def __getattr__(name):
        if not _chemicals_loaded:
            _load_chemicals()
            if name == 'chemicals': return chemicals
        if not _system_loaded: 
            try:
                _load_system()
            except Exception as Error:
                dct = globals()
                dct.update(flowsheet.system.__dict__)
                dct.update(flowsheet.stream.__dict__)
                dct.update(flowsheet.unit.__dict__)
                if name in dct: return dct[name]
                raise Error
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
else: 
    load()

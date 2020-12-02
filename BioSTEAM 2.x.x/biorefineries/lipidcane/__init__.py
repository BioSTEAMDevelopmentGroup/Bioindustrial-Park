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
               _process_settings,
               _chemicals,
               _system,
               _tea,
)

__all__ = [*utils.__all__,
           *_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           'lipidcane_sys',
           'lipidcane_tea', 
           'flowsheet',
]

from .utils import *
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
    global lipidcane_sys, lipidcane_tea, specs, flowsheet, _system_loaded
    flowsheet = bst.Flowsheet('lipidcane')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    lipidcane_sys = create_system()
    lipidcane_sys.simulate()
    lipidcane_tea = create_tea(lipidcane_sys)
    lipidcane_tea.IRR = lipidcane_tea.solve_IRR()
    _system_loaded = True

if PY37:    
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
else:
    load()

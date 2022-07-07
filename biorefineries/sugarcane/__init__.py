# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import (_process_settings,
               _chemicals,
               systems,
               _tea,
)

__all__ = [*_process_settings.__all__,
           *_chemicals.__all__,
           *systems.__all__,
           *_tea.__all__,
           'sugarcane_sys',
           'sugarcane_tea', 
           'flowsheet',
]

from ._process_settings import *
from ._chemicals import *
from .systems import *
from ._tea import *

_system_loaded = False
_chemicals_loaded = False

def load(pellet_bagasse=None):
    if not _chemicals_loaded: _load_chemicals()
    try: 
        _load_system(pellet_bagasse)
    finally: 
        dct = globals()
        dct.update(flowsheet.to_dict())

def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def _load_system(pellet_bagasse=None):
    import biosteam as bst
    from biosteam import main_flowsheet as F
    global sugarcane_sys, sugarcane_tea, flowsheet, tea, sys, _system_loaded
    flowsheet = bst.Flowsheet('sugarcane')
    F.set_flowsheet(flowsheet)
    bst.Stream.ticket_name = 's'
    bst.Stream.ticket_number = 0
    bst.Stream.unregistered_ticket_number = 0
    bst.System.ticket_name = 'SYS'
    bst.System.ticket_number = 0
    bst.System.unregistered_ticket_number = 0
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    sys = sugarcane_sys = create_sugarcane_to_ethanol_system(pellet_bagasse=pellet_bagasse)
    sugarcane_sys.simulate()
    tea = sugarcane_tea = create_tea(sugarcane_sys)
    sugarcane_tea.IRR = sugarcane_tea.solve_IRR()
    _system_loaded = True
  
def __getattr__(name):
    if not _chemicals_loaded:
        _load_chemicals()
        if name == 'chemicals': return chemicals
    if not _system_loaded: 
        try:
            _load_system()
        finally:
            dct = globals()
            dct.update(flowsheet.to_dict())
        if name in dct: return dct[name]
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

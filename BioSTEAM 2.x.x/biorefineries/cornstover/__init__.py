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
    global cornstover_sys, cornstover_tea, specs, flowsheet, _system_loaded
    global Area100, Area200, Area300, Area400, Area500, Area600, Area700, Area800
    global AllAreas, areas, ethanol_price_gal
    flowsheet = bst.Flowsheet('cornstover')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    cornstover_sys = create_system()
    cornstover_sys.simulate()
    u = F.unit
    OSBL_units = (u.WWTC, u.CWP, u.CT, u.PWC, u.ADP,
                  u.T701, u.T702, u.P701, u.P702, u.M701, u.FT,
                  u.CSL_storage, u.DAP_storage, u.BT)
    cornstover_tea = create_tea(cornstover_sys, OSBL_units, [u.U101])
    ethanol = F.stream.ethanol
    ethanol.price = cornstover_tea.solve_price(ethanol)
    ethanol_price_gal = ethanol.price * ethanol_density_kggal
    UnitGroup = bst.process_tools.UnitGroup
    Area100 = UnitGroup('Area 100', (u.U101,))
    Area200 = UnitGroup('Area 200', (u.T201, u.M201, u.R201, u.P201,
                                    u.T202, u.F201, u.H201, u.M210, u.T203))
    Area300 = UnitGroup('Area 300', (u.H301, u.M301, u.R301,
                                    u.R302, u.T301, u.T302))                 
    Area400 = UnitGroup('Area 400', (u.D401, u.H401, u.D402, u.P401,
                                    u.M402, u.D403, u.P402, u.H402,
                                    u.U401, u.H403, u.M701, u.S401))
    Area500 = UnitGroup('Area 500', (u.WWTC,))
    Area600 = UnitGroup('Area 600', (u.T701, u.T702, u.P701, u.P702, u.M701, u.FT,
                                    u.CSL_storage, u.DAP_storage))
    Area700 = UnitGroup('Area 700', (u.BT,))
    Area800 = UnitGroup('Area 800', (u.CWP, u.CT, u.PWC, u.ADP, u.CIP_package))
    areas = (Area100, Area200, Area300, Area400,
             Area500, Area600, Area700, Area800)
    AllAreas = UnitGroup('All Areas', cornstover_sys.units)
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
del PY37
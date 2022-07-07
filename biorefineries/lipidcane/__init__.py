# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
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
    _load_system(globals())

def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def _load_system(dct):
    import biosteam as bst
    from biosteam import main_flowsheet as F, UnitGroup
    global lipidcane_sys, lipidcane_tea, sys, tea, specs, flowsheet, _system_loaded
    global biodiesel_production_units, ethanol_production_units
    flowsheet = bst.Flowsheet('lipidcane')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    sys = lipidcane_sys = create_lipidcane_to_biodiesel_and_conventional_ethanol_system()
    lipidcane_sys.simulate()
    tea = lipidcane_tea = create_tea(lipidcane_sys)
    lipidcane_tea.IRR = lipidcane_tea.solve_IRR()
    dct.update(flowsheet.to_dict())
    biodiesel_production_units = UnitGroup('Biodisel production units', 
                                           [T401, P401, T402, P402, T403,
                                            P403, T404, P404, S401, R401,
                                            C401, P405, R402, R402, C402,
                                            T405, P406, C403, F401, P407,
                                            H401, P408, T406, P409, C404,
                                            T407, P410, D401, H402, D402,
                                            H403, P411, H404, P412, T408,
                                            T409])
    ethanol_production_units = UnitGroup('Ethanol production units',
                                         [S301, F301, P306, M301, H301,
                                          T305, R301, T301, C301, D301,
                                          M302, P301, H302, D302, P302,
                                          M303, D303, H303, U301, H304,
                                          T302, P304, T303, P305, M304,
                                          T304, P303])
    _system_loaded = True

def __getattr__(name):
    if not _chemicals_loaded:
        _load_chemicals()
        if name == 'chemicals': return chemicals
    if not _system_loaded: 
        dct = globals()
        _load_system(dct)
        if name in dct: return dct[name]
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

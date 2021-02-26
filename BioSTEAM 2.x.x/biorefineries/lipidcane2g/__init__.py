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
           'lipidcane2g_sys',
           'lipidcane2g_tea', 
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
    dct = globals()
    _load_system(dct)

def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def _load_system(dct):
    import biosteam as bst
    from biosteam import main_flowsheet as F, UnitGroup
    global lipidcane2g_sys, lipidcane2g_tea, specs, flowsheet, _system_loaded
    global unit_groups
    flowsheet = bst.Flowsheet('lipidcane2g')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    lipidcane2g_sys = create_lipidcane_to_biodiesel_and_both_cellulosic_and_conventional_ethanol_system()
    lipidcane2g_sys.simulate()
    u = flowsheet.unit
    bst.rename_unit(u.BT, 1000)
    bst.rename_units([u.FT, u.CWP, u.CIP_package, u.ADP, u.CT, u.PWC], 1100)
    bst.rename_units([i for i in lipidcane2g_sys.units if bst.is_storage_unit(i)], 1200)
    unit_groups = UnitGroup.group_by_area(lipidcane2g_sys.units)
    area_names = ['Feedstock handling', 
                  'Juicing', 
                  'Biod. prod.',
                  'Conv. ferm.', 
                  'pretreatment',
                  'Cofementation',
                  'Ethanol sep.', 
                  'Wastewater treatment', 
                  'Boiler turbogenerator',
                  'Utilities',
                  'Storage']
    for i, j in zip(unit_groups, area_names): i.name = j
    lipidcane2g_tea = create_tea(lipidcane2g_sys)
    lipidcane2g_tea.IRR = lipidcane2g_tea.solve_IRR()
    dct.update(flowsheet.to_dict())
    _system_loaded = True

def cellulosic_ethanol():
    return flowsheet.unit.M801.ins[0].imass['Ethanol']

def conventional_ethanol():
    return flowsheet.unit.M801.ins[1].imass['Ethanol']

if PY37:    
    def __getattr__(name):
        if not _chemicals_loaded:
            _load_chemicals()
            if name == 'chemicals': return chemicals
        if not _system_loaded:
            dct = globals()
            _load_system(dct)
            if name in dct: return dct[name]
        raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
else:
    load()

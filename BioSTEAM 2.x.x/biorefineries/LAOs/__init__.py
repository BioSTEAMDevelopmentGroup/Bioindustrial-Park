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
               utils,
               _chemicals,
               _process_specifications,
               _system,
               _tea,
)

__all__ = [*units.__all__,
           *utils.__all__,
           *_chemicals.__all__,
           *_process_specifications.__all__,
           *_system.__all__,
           *_tea.__all__,
           'LAOs_sys',
           'LAOs_tea', 
           'flowsheet', 
           'unit_groups', 
           'OSBL_unit_group',
           'specs',
           'utils',
]

from .utils import *
from .units import *
from ._chemicals import *
from ._process_specifications import *
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
    global chemicals
    chemicals = create_chemicals()    

def _load_system():
    import biosteam as bst
    from biosteam import main_flowsheet as F
    global LAOs_sys, LAOs_tea, specs, flowsheet, unit_groups, OSBL_unit_group
    global _system_loaded, products
    flowsheet = bst.Flowsheet('LAOs')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    LAOs_sys = create_system()
    OSBL_units = (F.unit.CWP, F.unit.BT, F.unit.CT, F.unit.T101, F.unit.T102, 
                  F.unit.T103, F.unit.T104, F.unit.T107, F.unit.T108, 
                  F.unit.T109, F.unit.T110)
    LAOs_tea = create_tea(LAOs_sys, OSBL_units)
    for i in LAOs_tea.TEAs: i.duration = (2017, 2047)
    for i in LAOs_tea.TEAs: i.contingency = 0.3
    specs = LAOsProcessSpecifications(LAOs_sys, LAOs_tea)
    UnitGroup = bst.process_tools.UnitGroup
    name_types = {
        'Tanks': bst.Tank,
        'Pumps': bst.Pump,
        'Facilities':bst.Facility,
        'Fermentation': bst.BatchBioreactor,
        'Heat exchangers': bst.HX,
        'Dehydration reactor': units.AdiabaticFixedbedGasReactor,
        'Mixer-settler': (bst.MixerSettler),
        'Centrifuges': (bst.LiquidsCentrifuge,
                        SolidLiquidsSplitCentrifuge),
        'Distillation': bst.BinaryDistillation
    }
    unit_groups = UnitGroup.group_by_types(LAOs_tea.units, name_types)
    OSBL_unit_group = UnitGroup('OSBL', OSBL_units)
    bst.System.molar_tolerance = 0.1
    bst.System.converge_method = 'aitken'
    specs.run_specifications() # Sets process specifications and simulates system
    products = (F('hexene'), F('octene'), F('decene'))
    for i in range(2): set_LAOs_MPSP(get_LAOs_MPSP())
    _system_loaded = True

if PY37:
    def __getattr__(name):
        if not _chemicals_loaded:
            _load_chemicals()
            if name == 'chemicals': return chemicals
        if not _system_loaded: _load_system()
        dct = globals()
        dct.update(flowsheet.system.__dict__)
        dct.update(flowsheet.stream.__dict__)
        dct.update(flowsheet.unit.__dict__)
        if name in dct:
            return dct[name]
        else:
            raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
else:
    load()
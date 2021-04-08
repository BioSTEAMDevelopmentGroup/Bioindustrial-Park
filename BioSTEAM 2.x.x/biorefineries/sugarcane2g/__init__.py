# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
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


def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True


def load(name):
    import biosteam as bst
    from biosteam import main_flowsheet as F, UnitGroup
    global sugarcane_sys, sugarcane_tea, specs, flowsheet, _system_loaded
    global unit_groups
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('sugarcane2g')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    dct = globals()
    u = flowsheet.unit
    s = flowsheet.stream
    sugarcane_sys = create_sugarcane_to_ethanol_2g()
    if name == 'sugarcane2g':
        # sugarcane_tea = create_tea(sugarcane_sys)
        # sugarcane_tea.operating_days = 200
        # sugarcane_sys.simulate()
        sugarcane_sys.simulate()
        scenario_a = sugarcane_sys.get_scenario_costs(200 * 24.)
        cornstover_sys = trim_to_cornstover_hot_water_cellulosic_ethanol(sugarcane_sys)
        cornstover_sys.simulate()
        scenario_b = cornstover_sys.get_scenario_costs(130 * 24.)
        sugarcane_tea = create_agile_tea([scenario_a, scenario_b])
        # sugarcane_tea.IRR = 0.1
        dct.update(flowsheet.to_dict())
        # s.ethanol.price = sugarcane_tea.solve_price(s.ethanol)

def ethanol_price():
    return sugarcane_tea.solve_price(flowsheet.stream.ethanol) * 2.98668849
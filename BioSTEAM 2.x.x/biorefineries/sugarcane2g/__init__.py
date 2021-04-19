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
    global sugarcane_sys, sugarcane_tea, specs, flowsheet, _system_loaded, agile_sugarcane_sys
    global unit_groups
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('sugarcane2g')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    dct = globals()
    u = flowsheet.unit
    s = flowsheet.stream
    sugarcane_sys = create_sugarcane_to_ethanol_2g(operating_hours=200 * 24)
    if name == 'sugarcane2g':
        # sugarcane_tea = create_tea(sugarcane_sys)
        # sugarcane_tea.operating_days = 200
        # sugarcane_sys.simulate()
        sugarcane_tea = create_agile_tea(sugarcane_sys.units)
        agile_sugarcane_sys = AgileSugarcaneSystem(
            sugarcane_sys, [3.3e5, 3.3e5], [2400, 2400], s.sugarcane
        )
        agile_sugarcane_sys.simulate()
        scenario_a = sugarcane_tea.create_scenario(agile_sugarcane_sys)
        cornstover_sys = trim_to_cornstover_hot_water_cellulosic_ethanol(
            sugarcane_sys,
            operating_hours=2400.,
        )
        cornstover_sys.simulate()
        scenario_b = sugarcane_tea.create_scenario(cornstover_sys)
        sugarcane_tea.compile_scenarios([scenario_a, scenario_b])
        # sugarcane_tea.IRR = 0.1
        dct.update(flowsheet.to_dict())
        # s.ethanol.price = sugarcane_tea.solve_price(s.ethanol)

def ethanol_price():
    return sugarcane_tea.solve_price(flowsheet.stream.ethanol) * 2.98668849
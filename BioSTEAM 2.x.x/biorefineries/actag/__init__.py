# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from ._chemicals import create_conventional_chemicals, create_cellulosic_chemicals
from ._system import (
    create_cellulosic_acTAG_system,
    create_conventional_acTAG_system,
)
from ._units import *

_conventional_chemicals_loaded = False
_cellulosic_chemicals_loaded = False

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 607.5 # 2019
    bst.PowerUtility.price = 0.065
    HeatUtility = bst.HeatUtility
    steam_utility = HeatUtility.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    HeatUtility.get_agent('cooling_water').regeneration_price = 0
    HeatUtility.get_agent('chilled_water').heat_transfer_price = 0

def load_conventional_chemicals():
    global conventional_chemicals, _chemicals_loaded
    conventional_chemicals = create_conventional_chemicals()
    _conventional_chemicals_loaded = True

def load_cellulosic_chemicals():
    global cellulosic_chemicals, _chemicals_loaded
    cellulosic_chemicals = create_cellulosic_chemicals()
    _cellulosic_chemicals_loaded = True

def load(configuration):
    configuration = int(configuration)
    if configuration == 1:
        if not _conventional_chemicals_loaded: 
            load_conventional_chemicals()
            chemicals = conventional_chemicals
        flowsheet = bst.Flowsheet('conventional_acTAG')
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(chemicals)
        sys = create_conventional_acTAG_system()
    elif configuration == 2:
        if not _cellulosic_chemicals_loaded: 
            load_cellulosic_chemicals()
            chemicals = cellulosic_chemicals
        flowsheet = bst.Flowsheet('cellulosic_acTAG')
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(chemicals)
        sys = create_cellulosic_acTAG_system()
    else:
        raise ValueError(f"invalid configuration '{configuration}'; only 1 and 2 are valid")
    u = flowsheet.unit
    s = flowsheet.stream
    load_process_settings()
    return sys
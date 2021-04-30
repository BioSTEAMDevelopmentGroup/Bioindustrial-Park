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
from . import (utils,
               _process_settings,
               _chemicals,
               _system,
               _tea,
               _lipid_extraction_specification,
)

__all__ = [*utils.__all__,
           *_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           *_lipid_extraction_specification.__all__,
           'lipidcane_sys',
           'lipidcane_tea', 
           'flowsheet',
]

from .utils import *
from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *
from ._lipid_extraction_specification import *
from biorefineries.sugarcane import create_tea as create_conventional_ethanol_tea

_system_loaded = False
_chemicals_loaded = False

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def load(name, agile=False):
    import biosteam as bst
    from biosteam import main_flowsheet as F, UnitGroup
    global lipidcane_sys, lipidcane_tea, specs, flowsheet, _system_loaded
    global lipid_extraction_specification
    global unit_groups
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('lipidcane2g')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    dct = globals()
    u = flowsheet.unit
    s = flowsheet.stream
    operating_hours = 24 * 200
    area_names = None
    def rename_storage_units(storage):
        bst.rename_units([i for i in lipidcane_sys.units if bst.is_storage_unit(i)], storage)
    
    if name == 5 or name == 'divided 1 and 2g front end oil separation':
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,    
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Biod. prod.',
            'Conv. ferm.', 
            'Pretreatment',
            'Cofementation',
            'Ethanol sep.', 
            'Wastewater treatment', 
            'Boiler turbogenerator',
            'Utilities',
            'Storage'
        ]
        rename_storage_units()
    elif name == 4 or name == 'divided 1 and 2g hydrolyzate oil separation':
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Biod. prod.', 
            'Conv. ferm.', 
            'Pretreatment', 
            'Cofementation',
            'Ethanol sep.', 
            'Wastewater treatment', 
            'Boiler turbogenerator',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_hydrolyzate_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units()
    elif name == 3 or name == 'combined 1 and 2g post fermentation oil separation':
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'Oil ext.',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1200)
    elif name == 2 or name == 'divided 1 and 2g post fermentation oil separation':
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Conv. ferm.',
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'Oil ext.',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1300)
    elif name == 1 or name == 'divided 1 and 2g bagasse expression':
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Conv. ferm.',
            'Oil ext.',
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_1_and_2g_bagasse_expression(
            operating_hours=operating_hours,
        )
    elif name == 0 or name == '1g':
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_1g(
            operating_hours=operating_hours,
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'EtOH prod.', 
            'Oil ext.',
            'Biod. prod.', 
            'CH&P',
            'Utilities',
            'Storage',
        ]
        rename_storage_units(900)
    else:
        raise NotImplementedError(name)
    unit_groups = UnitGroup.group_by_area(lipidcane_sys.units)
    if area_names:
            for i, j in zip(unit_groups, area_names): i.name = j
    if name in (1, 2, 3, 4) if isinstance(name, int) else '2g' in name:
        if agile:
            lipidcane_tea = create_agile_tea(lipidcane_sys.units)
            lipidcane_sys.simulate()
            scenario_a = sugarcane_tea.create_scenario(lipidcane_sys)
            cornstover_sys = trim_to_cornstover_hot_water_cellulosic_ethanol(
                lipidcane_sys,
                operating_hours=24 * 100,
            )
            cornstover_sys.simulate()
            scenario_b = lipidcane_tea.create_scenario(cornstover_sys)
            sugarcane_tea.compile_scenarios([scenario_a, scenario_b])
            dct.update(flowsheet.to_dict())
            lipidcane_tea.IRR = 0.10
            s.ethanol.price = lipidcane_tea.solve_price(s.ethanol)
            return
    lipidcane_tea = create_tea(lipidcane_sys)
    lipidcane_tea.operating_days = 200
    try: 
        lipidcane_sys.simulate()
    except Exception as e:
        raise e
    else:
        lipidcane_tea.IRR = 0.10
    finally:
        dct.update(flowsheet.to_dict())
    
    for i in lipidcane_sys.units:
        if getattr(i, 'tag', None) == 'lipid extraction efficiency unit a':
            isplit_a = i.isplit
            break
    
    for i in lipidcane_sys.units:
        if getattr(i, 'tag', None) == 'lipid extraction efficiency unit b':
            isplit_b = i.isplit
            break
    
    lipid_extraction_specification = LipidExtractionSpecification(
            lipidcane_sys, lipidcane, isplit_a, isplit_b,
        )
        
        

def ethanol_price():
    return lipidcane_tea.solve_price(flowsheet.stream.ethanol) * 2.98668849
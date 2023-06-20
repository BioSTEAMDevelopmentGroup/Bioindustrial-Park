# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import os
os.environ["NUMBA_DISABLE_JIT"] = '1' # In case numba or numba cache not working properly
import numpy as np
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import pytest
from biosteam.process_tools import UnitGroup
from importlib import import_module

__all__ = (
    'test_sugarcane',
    'test_lipidcane',
    'test_cornstover',
    'test_LAOs',
    'test_lactic',
    'test_ethanol_adipic',
    'generate_all_code',
    'generate_code',
    'print_results',
)

feedstocks_by_module = {
    'corn': 'corn',
    'lipidcane': 'lipidcane',
    'cornstover': 'cornstover',
    'sugarcane': 'sugarcane',
    'oilcane': 'feedstock',
    'LAOs': 'glucose',
    'lactic': 'feedstock',
    'ethanol_adipic': 'feedstock',
    'wheatstraw': 'wheatstraw',
    'animal_bedding': 'bedding',
    # 'HP': 'feedstock',
}
products_by_module = {
    'corn': 'ethanol',
    'lipidcane': 'ethanol',
    'cornstover': 'ethanol',
    'sugarcane': 'ethanol',
    'LAOs': 'octene',
    'wheatstraw': 'ethanol',
    'animal_bedding': 'ethanol',
    'lactic': 'lactic_acid',
    'ethanol_adipic': 'ethanol',
    'HP': 'AcrylicAcid',
    'oilcane': 'ethanol',
}

must_load = {
    'oilcane', 'corn', 'sugarcane', 'cornstover', 'lipidcane'
}

marked_slow = {'wheatstraw', 'animal_bedding'}

configurations = {
    'HP': ('cellulosic', 'sugarcane'),
    'oilcane': ('S1', 'S2', 'O1', 'O2', 'O3', 'O4', 'S1*', 'S2*', 'O1*', 'O2*'),
}

def default_settings(f):
    def wrapper(*args, **kwargs):
        bst.process_tools.default()
        try:
            return f(*args, **kwargs)
        finally:
            bst.process_tools.default()
    wrapper.__name__ = f.__name__
    return wrapper
        

def generate_code(module_name, feedstock_name=None, product_name=None, configuration=None):
    if not feedstock_name:
        feedstock_name = feedstocks_by_module[module_name]
    if not product_name:
        product_name = products_by_module[module_name]
    bst.process_tools.default()
    module = import_module('biorefineries.' + module_name)
    if module_name in must_load:
        if configuration is None:
            module.load()
        else:
            module.load(configuration)
    else:
        try: 
            if configuration is None:
                pass
            else:
                module.load(configuration)
        except: pass
    feedstock = getattr(module, feedstock_name)
    product = getattr(module, product_name)
    for tea_name in ('tea', f'{module_name}_tea',  f'{feedstock_name}_tea', f'{product_name}_tea'):
        try:
            tea = getattr(module, tea_name)
        except AttributeError:
            continue
        else:
            break
    try: units = UnitGroup('Biorefinery', tea.units)
    except:
        breakpoint()
    IRR = tea.IRR
    sales = tea.sales
    material_cost = tea.material_cost
    installed_equipment_cost = tea.installed_equipment_cost
    utility_cost = tea.utility_cost
    heating_duty = units.get_heating_duty()
    cooling_duty = units.get_cooling_duty()
    electricity_consumption = units.get_electricity_consumption()
    electricity_production = units.get_electricity_production()
    configuration_tag = f"_{configuration}".replace('*', '_agile') if configuration else ''
    configuration_name = f"'{configuration}'" if configuration else ''
    print(
    ("@pytest.mark.slow\n" if module_name in marked_slow else "") +
    f"""@default_settings
def test_{module_name}{configuration_tag}():
    from biorefineries import {module_name} as module
    if '{module_name}' in must_load:
        module.load({configuration_name})
    else:
        try: module.load({configuration_name})
        except: pass
    feedstock = module.{feedstock_name}
    product = module.{product_name}
    tea = module.{tea_name}
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, {IRR}, rtol=5e-2)
    assert np.allclose(feedstock.price, {feedstock.price}, rtol=5e-2)
    assert np.allclose(product.price, {product.price}, rtol=5e-2)
    assert np.allclose(tea.sales, {sales}, rtol=5e-2)
    assert np.allclose(tea.material_cost, {material_cost}, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, {installed_equipment_cost}, rtol=5e-2)
    assert np.allclose(tea.utility_cost, {utility_cost}, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), {heating_duty}, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), {cooling_duty}, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), {electricity_consumption}, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), {electricity_production}, rtol=5e-2)
    """
    )

def generate_all_code():
    from warnings import filterwarnings
    filterwarnings('ignore')
    for i in feedstocks_by_module: 
        if i in configurations:
            for j in configurations[i]:
                generate_code(i, configuration=j)
        else:
            generate_code(i)
        
def print_results(tea):
    units = UnitGroup('Biorefinery', tea.units)
    print('Sales:', tea.sales)
    print('Material cost:', tea.material_cost)
    print('Installed equipment cost:', tea.installed_equipment_cost)
    print('Utility cost:', tea.utility_cost)
    print('Heating duty:', units.get_heating_duty())
    print('Cooling duty:', units.get_cooling_duty())
    print('Electricity consumption:', units.get_electricity_consumption())
    print('Electricity production:', units.get_electricity_production())

@default_settings
def test_corn():
    from biorefineries import corn as module
    if 'corn' in must_load:
        module.load()
    else:
        try: module.load()
        except: pass
    feedstock = module.corn
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.04769406769940579, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74973086.86827628, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52917743.6806844, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 65987734.083792426, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 10274134.399572838, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 115.31849520368542, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 105.52628418496067, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.0669963081288474, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=5e-2)
    
@default_settings
def test_lipidcane():
    from biorefineries import lipidcane as module
    if 'lipidcane' in must_load:
        module.load()
    else:
        try: module.load()
        except: pass
    feedstock = module.lipidcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.20865355832165103, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693809.0099629, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798770.286153086, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 224150640.60415673, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34735955.640232466, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 201.57736544521418, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 249.89918148814343, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.039308145969065, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 124.04147980101291, rtol=5e-2)
    
@default_settings
def test_cornstover():
    from biorefineries import cornstover as module
    if 'cornstover' in must_load:
        module.load()
    else:
        try: module.load()
        except: pass
    feedstock = module.cornstover
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.05158816935126135, rtol=5e-2)
    assert np.allclose(product.price, 0.7154530694282512, rtol=5e-2)
    assert np.allclose(tea.sales, 132154005.82497102, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82365283.00148688, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207924897.31670856, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -10953746.157227168, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 357.42602131639455, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 374.0397216764084, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.504608092973328, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.04230437729655, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    if 'sugarcane' in must_load:
        module.load()
    else:
        try: module.load()
        except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.13417361425007043, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301963.96150482, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165736.478367165, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 198936569.1451582, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -17067514.575046588, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 248.0565662637677, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 289.52405057724155, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.993365934157058, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 68.44784811552151, rtol=5e-2)
    
@default_settings
def test_oilcane_S1():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('S1')
    else:
        try: module.load('S1')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.038419151419652874, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81216602.83525659, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63362758.2628729, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 133982834.75687276, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19381517.804545257, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 306.04784820793185, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 295.7889692925312, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.993317735219735, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 77.8229743904066, rtol=5e-2)
    
@default_settings
def test_oilcane_S2():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('S2')
    else:
        try: module.load('S2')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.038125778058066086, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 142066649.40761417, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74131304.72697282, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 242904199.5263379, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1798657.9002356618, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 382.4812949202671, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 403.21631938552963, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.927915521607964, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 24.927915521607964, rtol=5e-2)
    
@default_settings
def test_oilcane_O1():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O1')
    else:
        try: module.load('O1')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.04407758588479899, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 73143973.7292505, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74590576.60482024, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 170217577.13064548, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46687809.76778478, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 284.86503068488025, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.7954819103185, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.457554468398685, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 166.29101959757548, rtol=5e-2)
    
@default_settings
def test_oilcane_O2():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O2')
    else:
        try: module.load('O2')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.05077816261667996, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166847026.91905588, rtol=5e-2)
    assert np.allclose(tea.material_cost, 100526990.29969029, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 267103924.97999448, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6370983.569182227, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 407.5161995687945, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 406.69515229531396, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.036553507876388, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 49.31435445122779, rtol=5e-2)
    
@default_settings
def test_oilcane_O3():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O3')
    else:
        try: module.load('O3')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03877900972715943, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 60656253.07557489, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63596223.773526296, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 162548796.0438482, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46123794.55702143, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 271.4745501463212, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 259.72335801538134, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.150398899240725, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 164.13856910788985, rtol=5e-2)
    
@default_settings
def test_oilcane_O4():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O4')
    else:
        try: module.load('O4')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.04199211588322421, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148592757.57557717, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82809352.97257052, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 255936544.69088632, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4004922.5124022355, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 387.83630319009046, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 382.997868492699, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.727407059775523, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 41.31516106332716, rtol=5e-2)
    
@default_settings
def test_oilcane_S1_agile():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('S1*')
    else:
        try: module.load('S1*')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.04031531995389104, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 79579640.10023846, rtol=5e-2)
    assert np.allclose(tea.material_cost, 66357343.61215794, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 120097913.5272843, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -20611222.96492123, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 140.15227969846998, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 165.21529788883913, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.5016800601044, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 72.60376471164439, rtol=5e-2)
    
@default_settings
def test_oilcane_S2_agile():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('S2*')
    else:
        try: module.load('S2*')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.041898042468137255, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 140030748.4514662, rtol=5e-2)
    assert np.allclose(tea.material_cost, 80334075.7114395, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 214443929.9573282, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 780612.1419791012, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 136.89152832842447, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 218.9592193351276, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.195104598662823, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 28.645246107788267, rtol=5e-2)
    
@default_settings
def test_oilcane_O1_agile():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O1*')
    else:
        try: module.load('O1*')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.046152230303086246, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71543985.23665722, rtol=5e-2)
    assert np.allclose(tea.material_cost, 77193974.57287961, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 147306960.4218294, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45540323.84092915, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 192.8021481616401, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 187.41525341217994, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.319365795289885, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 113.2227750973087, rtol=5e-2)
    
@default_settings
def test_oilcane_O2_agile():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O2*')
    else:
        try: module.load('O2*')
        except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.054975598100391716, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 162273402.16250658, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106406849.11189632, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 228546356.60882708, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6983748.004429089, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 179.1217218313402, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 231.64685009631594, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.390858975864738, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 43.28862817101215, rtol=5e-2)
    
    
@default_settings
def test_LAOs():
    from biorefineries import LAOs as module
    try: module.load()
    except: pass
    feedstock = module.glucose
    product = module.octene
    tea = module.LAOs_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.265, rtol=5e-2)
    assert np.allclose(product.price, 1.289861475509996, rtol=5e-2)
    assert np.allclose(tea.sales, 161175673.52938324, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.45874006, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 54699221.91419212, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2531383.5854514386, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.310783441715046, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.69592718862404, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.5278713460948774, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.527871346094877, rtol=5e-2)
    
### DO NOT DELETE:
@pytest.mark.slow
def test_lactic():
    from biorefineries import lactic
    bst.process_tools.default_utilities()
    system = lactic.system
    MPSP = system.lactic_acid.price
    assert np.allclose(MPSP, 1.5285811040727408, atol=0.01)
    tea = system.lactic_tea
    assert np.allclose(tea.sales, 332247972.3322271, rtol=0.01)
    assert np.allclose(tea.material_cost, 218911499.70700422, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 321744540.84435904, rtol=0.01)
    assert np.allclose(tea.utility_cost, 24995492.58140952, rtol=0.01)
    units = UnitGroup('Biorefinery', system.lactic_tea.units)
    assert np.allclose(system.CHP.system_heating_demand/1e6, 1763.4942302374052, rtol=0.01)
    assert np.allclose(-system.CT.system_cooling_water_duty/1e6, 1629.0964343101657, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 45.291535445041546, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 0.0)

@pytest.mark.slow
def test_ethanol_adipic():
    bst.process_tools.default_utilities()
    from biorefineries import ethanol_adipic
    acid = ethanol_adipic.system_acid
    MESP = acid.ethanol.price * acid._ethanol_kg_2_gal
    assert np.allclose(MESP, 2.486866594017598, atol=0.01)
    tea = acid.ethanol_tea
    assert np.allclose(tea.sales, 152002785.9542236, rtol=0.01)
    assert np.allclose(tea.material_cost, 112973059.93909653, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 202880204.7890376, rtol=0.01)
    assert np.allclose(tea.utility_cost, -18225095.0315181, rtol=0.01)
    units = UnitGroup('Biorefinery', acid.ethanol_tea.units)
    assert np.allclose(acid.CHP.system_heating_demand/1e6, 333.8802418517945, rtol=0.01)
    assert np.allclose(units.get_cooling_duty(), 327.05990128304114, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 23.84999102548279, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 55.72024685271333, rtol=0.01)
    
    base = ethanol_adipic.system_base
    MESP = base.ethanol.price * base._ethanol_kg_2_gal
    assert np.allclose(MESP, 2.703699730342356, atol=0.01)
    tea = base.ethanol_adipic_tea
    assert np.allclose(tea.sales, 219850046.2308626, rtol=0.01)
    assert np.allclose(tea.material_cost, 120551649.9082640, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 283345546.0556234, rtol=0.01)
    assert np.allclose(tea.utility_cost, 14241964.663629433, rtol=0.01)
    units = UnitGroup('Biorefinery', base.ethanol_adipic_tea.units)
    assert np.allclose(base.CHP.system_heating_demand/1e6, 296.8322623939439, rtol=0.01)
    assert np.allclose(units.get_cooling_duty(), 247.79573940245342, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 26.565278642577344, rtol=0.001)
    assert np.allclose(units.get_electricity_production(), 0.0)

if __name__ == '__main__':
    generate_all_code()
    # test_corn()
    # test_sugarcane()
    # test_lipidcane()
    # test_cornstover()
    # test_oilcane_S1()
    # test_oilcane_S2()
    # test_oilcane_O1()
    # test_oilcane_O2()
    # test_oilcane_O3()
    # test_oilcane_O4()
    # test_oilcane_S1_agile()
    # test_oilcane_S2_agile()
    # test_oilcane_O1_agile()
    # test_oilcane_O2_agile()
    # test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





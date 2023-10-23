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
    assert np.allclose(tea.IRR, 0.053222010734474134, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74973086.8078639, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52917743.62293877, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 66066236.03106276, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 9924897.194763457, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 108.48791237650839, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 99.98403111215264, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.0612291376300735, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.20974921103082406, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693808.6378837, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798689.823035285, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 223273551.4686045, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34849320.92391035, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 200.20930252115042, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 249.40864010555856, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.998971553763419, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 124.36449465808356, rtol=5e-2)
    
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
    assert np.allclose(tea.utility_cost, -10953746.157227173, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 357.42602131639455, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 374.03972167640836, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.504608092973335, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.04230437729657, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.13755002253674653, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301963.95368072, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165736.47820579, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 199128913.5350583, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18004656.781389583, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 236.74383547570676, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.26921703638527, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.668643379652474, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.11890955159149, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03901749327797711, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81216602.88499977, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64320172.25829746, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 133966028.3905485, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -20310170.196726065, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 294.67418727093514, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 281.4917276514267, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.66856129694358, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 80.5084221116419, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038125784357298784, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 142066647.74009016, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74131313.97250621, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 242904216.0705572, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1798643.1704908144, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 382.4813120260001, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 403.21633696874926, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.927918146885965, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 24.927918146885983, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.043995511741750766, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 73143973.72688098, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74458166.31612569, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169634385.9646353, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46392020.27493236, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 285.5614762751667, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 274.86143471215365, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.416471256578337, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 165.29235054170647, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05084011135645966, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166847026.91884053, rtol=5e-2)
    assert np.allclose(tea.material_cost, 100626932.36037634, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 263585001.22436345, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5520716.121781402, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 393.22365951257257, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 392.3914039071641, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.03824228012872, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 46.552563480957716, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038879522027389125, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 60656253.09414633, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63758380.39277257, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 161989353.65385896, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46125094.93885466, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 271.45664327793725, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 259.7260272962395, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.150399246222575, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 164.1427990117023, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04304175741255468, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148592757.58234242, rtol=5e-2)
    assert np.allclose(tea.material_cost, 84502740.88709588, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 254145611.657441, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5118624.683903244, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 372.4957125467686, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 371.8065124343648, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.729813681301188, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 44.93724453002472, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0408771523550341, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 79579640.13537754, rtol=5e-2)
    assert np.allclose(tea.material_cost, 67256295.46154399, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 119363338.20539875, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -21311176.975464605, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 140.16399013887926, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 165.21350884324144, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.5016814313403, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 72.60099876198963, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04189805056190648, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 140030747.069449, rtol=5e-2)
    assert np.allclose(tea.material_cost, 80334087.48191133, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 214443957.4431826, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 780592.7722834187, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 136.89152832482412, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 218.9592212050165, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.195104732416453, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 28.64524528880774, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0461615275423034, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71543985.23513563, rtol=5e-2)
    assert np.allclose(tea.material_cost, 77208850.16573285, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 146332802.00683832, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45287206.77045706, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 192.85742030861974, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 185.98904678664437, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.2506194650854425, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 112.92701121292627, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05487234681971159, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 162273402.1621236, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106241647.1315949, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 226272638.3933427, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6206778.157946691, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 167.42265795960347, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 220.73238151625182, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.39227250117013, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 41.364516816155934, rtol=5e-2)
    
@default_settings
def test_LAOs():
    from biorefineries import LAOs as module
    if 'LAOs' in must_load:
        module.load()
    else:
        try: module.load()
        except: pass
    feedstock = module.glucose
    product = module.octene
    tea = module.LAOs_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.265, rtol=5e-2)
    assert np.allclose(product.price, 1.2927068503051573, rtol=5e-2)
    assert np.allclose(tea.sales, 161541697.4292124, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135665307.92992145, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 55078641.5870012, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2763304.51643572, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.28531181149635, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.70084024745563, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.577370692927678, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.577370692927677, rtol=5e-2)
    
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
    # generate_all_code()
    test_corn()
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_oilcane_S1()
    test_oilcane_S2()
    test_oilcane_O1()
    test_oilcane_O2()
    test_oilcane_O3()
    test_oilcane_O4()
    test_oilcane_S1_agile()
    test_oilcane_S2_agile()
    test_oilcane_O1_agile()
    test_oilcane_O2_agile()
    test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





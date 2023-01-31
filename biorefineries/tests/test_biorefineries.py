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
    try: module.load(configuration)
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
    units = UnitGroup('Biorefinery', tea.units)
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
    try: module.load()
    except: pass
    feedstock = module.corn
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.05524120063646625, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723596.88174264, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525176.92772423, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 61891581.13350075, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7484262.229481168, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.0630956399827, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.83385306281221, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9369534236256791, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=5e-2)
    
@default_settings
def test_lipidcane():
    from biorefineries import lipidcane as module
    try: module.load()
    except: pass
    feedstock = module.lipidcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.21479807710155088, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.72636937, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58863950.04494265, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 224182525.149431, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -37093725.944828905, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 207.1049068884513, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 234.97872550762915, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.556381520720417, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 130.87263434050848, rtol=5e-2)
    
@default_settings
def test_cornstover():
    from biorefineries import cornstover as module
    try: module.load()
    except: pass
    feedstock = module.cornstover
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.05158816935126135, rtol=5e-2)
    assert np.allclose(product.price, 0.7079324847924882, rtol=5e-2)
    assert np.allclose(tea.sales, 130764818.4240719, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82181784.09989835, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 205211850.35825574, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -11505007.455217296, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 363.19777385364335, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 274.8385532341637, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.475394720502877, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 45.588907540821644, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.14081684743443543, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606105, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57282540.46485529, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 197050337.91986912, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18519108.798572544, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.72055792287915, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 236.41034434059142, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.725492410409335, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.49649020165383, rtol=5e-2)
    
@default_settings
def test_oilcane_S1():
    from biorefineries import oilcane as module
    try: module.load('S1')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.0389376887664696, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73735543, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64310168.74309213, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 132288800.07336637, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19883435.487783715, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 304.4518210497561, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 252.85735309138698, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.004170830859518, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 78.10484461452931, rtol=5e-2)
    
@default_settings
def test_oilcane_S2():
    from biorefineries import oilcane as module
    try: module.load('S2')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03480828971097124, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 131703149.88672155, rtol=5e-2)
    assert np.allclose(tea.material_cost, 67093269.68508662, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 232140169.4900541, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1376539.0221294055, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 349.38850511223006, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 328.1499248645642, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.060977086401, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.109449934047024, rtol=5e-2)
    
@default_settings
def test_oilcane_O1():
    from biorefineries import oilcane as module
    try: module.load('O1')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.04460813819226493, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71993360.47700194, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74937999.62896879, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169899934.9570202, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -48027307.71977071, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 252.196300240475, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 234.53114267956863, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.15620338632552, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 169.72816837060586, rtol=5e-2)
    
@default_settings
def test_oilcane_O2():
    from biorefineries import oilcane as module
    try: module.load('O2')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.05003897939324008, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166664496.71806905, rtol=5e-2)
    assert np.allclose(tea.material_cost, 98766407.3714235, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 259187639.068608, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2826493.7674901537, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 390.90675326426367, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 348.2320224340548, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.569728835728675, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 37.27454253081435, rtol=5e-2)
    
@default_settings
def test_oilcane_O3():
    from biorefineries import oilcane as module
    try: module.load('O3')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03917453784823141, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.7419263, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63796309.072414234, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 162106474.48227784, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47284373.07499168, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.03602210470234, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 227.6050914317077, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.923507783370574, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 167.0808291510632, rtol=5e-2)
    
@default_settings
def test_oilcane_O4():
    from biorefineries import oilcane as module
    try: module.load('O4')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.042059112633287626, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148409512.76797798, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82455517.24802905, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 249198456.31523743, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2119022.0128182517, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 379.93351816273235, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 335.7556017296543, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.240069544380084, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 34.645466260545234, rtol=5e-2)
    
@default_settings
def test_oilcane_S1_agile():
    from biorefineries import oilcane as module
    try: module.load('S1*')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.043188533528732025, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.41629672, rtol=5e-2)
    assert np.allclose(tea.material_cost, 92398314.84855245, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 139064433.32294938, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27382442.381128136, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 178.7763836857405, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 202.808452984611, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.30771674889143, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 93.74907148845497, rtol=5e-2)
    
@default_settings
def test_oilcane_S2_agile():
    from biorefineries import oilcane as module
    try: module.load('S2*')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.04234974307131911, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 171945869.47831404, rtol=5e-2)
    assert np.allclose(tea.material_cost, 103698757.96446578, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 243204080.48905522, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1686533.3713901546, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 157.72723170910083, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 319.96378472348573, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.316734581667276, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 25.509420384284013, rtol=5e-2)
    
@default_settings
def test_oilcane_O1_agile():
    from biorefineries import oilcane as module
    try: module.load('O1*')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.049141125413655624, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 92435850.58665347, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106653220.31743592, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 171256107.83772278, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -59821585.121718176, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 227.03389380488375, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 203.9329712577954, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.049211133094778, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 149.90405126461286, rtol=5e-2)
    
@default_settings
def test_oilcane_O2_agile():
    from biorefineries import oilcane as module
    try: module.load('O2*')
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.057289328022565034, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 210830375.18554768, rtol=5e-2)
    assert np.allclose(tea.material_cost, 143105486.91172585, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 263215511.72079462, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4268223.9268260375, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 311.076222133392, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.4730338465815, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.40997715651808, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 43.34032027386667, rtol=5e-2)
    
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
    # generate_all_code()
    test_corn()
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_LAOs()
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
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





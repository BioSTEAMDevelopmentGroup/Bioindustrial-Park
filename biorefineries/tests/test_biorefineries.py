# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
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
    tea = module.corn_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.054793135761420905, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723596.88174264, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525176.94552828, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 62093937.0411999, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7484262.234701365, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.0630956270271, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.8338531988122, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.936953423748958, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.21462041059411274, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.72636937, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58863950.04494265, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 224377181.991511, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -37093725.94482891, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 207.10490688845127, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 234.9787255076292, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.556381520720414, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 130.8726343405085, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.710598833148354, rtol=5e-2)
    assert np.allclose(tea.sales, 131327312.04763526, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82213161.58336543, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 206707328.22897246, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -11353990.27079505, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 363.3753814500929, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 274.940915842381, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.733008173236037, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 45.53260351027676, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.14081684760783747, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606105, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57282540.464877196, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 197050338.00378352, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18519108.863740746, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.72055703822804, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 236.41034447328917, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.725492410412995, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.49649041052992, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03893768879824197, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73735546, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64310168.79396416, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 132288800.12735988, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19883435.55277177, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 304.4518201551533, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 252.85735331009988, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.00417083086553, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 78.10484482575495, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.034808289549146804, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 131703149.27372651, rtol=5e-2)
    assert np.allclose(tea.material_cost, 67093269.32139692, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 232140169.73738986, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1376538.7481811198, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 349.3885049227761, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 328.1499293423478, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.06097711577237, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.109450854667372, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.044613328857648606, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71993360.47700194, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74946304.69359091, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169871042.53970197, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -48027307.69518688, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 252.19630057888367, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 234.53114262880734, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.15620338632412, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 169.72816829070376, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05004519235574291, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166664496.71806905, rtol=5e-2)
    assert np.allclose(tea.material_cost, 98776348.1114479, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 259153056.42877504, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2826493.7674901467, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 390.9067532642638, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 348.2320224340548, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.569728835728668, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 37.27454253081432, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03917453783330245, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.7419263, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63796309.04851951, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 162106474.47010747, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47284373.05040786, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.036022443111, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 227.6050913809463, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.923507783369175, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 167.08082907116116, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.042059112630315754, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148409512.76797798, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82455517.24327405, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 249198456.31523758, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2119022.012818315, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 379.93351816273133, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 335.7556017296544, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.24006954438008, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 34.645466260545454, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04318853355872526, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.41629674, rtol=5e-2)
    assert np.allclose(tea.material_cost, 92398314.91097715, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 139064433.3189067, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27382442.44225976, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 178.77638368574102, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 202.8084529704335, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.307716748893183, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 93.74907146613862, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04234959602200342, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 171945869.2723481, rtol=5e-2)
    assert np.allclose(tea.material_cost, 103698452.06622273, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 243205476.73650748, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1686533.387474242, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 157.72723103253057, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 319.96378415968206, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.316734581073774, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 25.509419497171976, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04914511824497407, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 92435850.58665347, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106661525.4065949, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 171227215.43249747, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -59821585.12171817, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 227.03389380488372, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 203.93297125779543, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.049211133094778, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 149.90405126461283, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0572941072852774, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 210830375.1855477, rtol=5e-2)
    assert np.allclose(tea.material_cost, 143115427.77818748, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 263180929.0799267, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4268224.05807124, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 311.07622213339306, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.4730338465808, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.409977156518096, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 43.34032027386641, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.29078979468983, rtol=5e-2)
    assert np.allclose(tea.sales, 161296136.29097626, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.45874006, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 55049531.41964978, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2531383.585451439, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.31078344171506, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.69592718862403, rtol=5e-2)
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





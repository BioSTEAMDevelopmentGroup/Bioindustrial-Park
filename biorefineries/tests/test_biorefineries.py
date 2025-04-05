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
    # 'test_LAOs',
    # 'test_lactic',
    # 'test_ethanol_adipic',
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
    # 'LAOs': 'glucose',
    # 'lactic': 'feedstock',
    # 'ethanol_adipic': 'feedstock',
    # 'wheatstraw': 'wheatstraw',
    # 'animal_bedding': 'bedding',
    # 'HP': 'feedstock',
}
products_by_module = {
    'corn': 'ethanol',
    'lipidcane': 'ethanol',
    'cornstover': 'ethanol',
    'sugarcane': 'ethanol',
    # 'LAOs': 'octene',
    # 'wheatstraw': 'ethanol',
    # 'animal_bedding': 'ethanol',
    # 'lactic': 'lactic_acid',
    # 'ethanol_adipic': 'ethanol',
    # 'HP': 'AcrylicAcid',
    'oilcane': 'ethanol',
    'O6': 'biodiesel',
    'O7': 'biodiesel', 
    'O8': 'biodiesel',
    'O9': 'biodiesel',
}

must_load = {
    'oilcane', 'corn', 'sugarcane', 'cornstover', 'lipidcane'
}

marked_slow = {'wheatstraw', 'animal_bedding'}

configurations = {
    'HP': ('cellulosic', 'sugarcane'),
    'oilcane': ('S1', 'S2', 'O1', 'O2', 'O3', 'O4', 'S1*', 'S2*', 'O1*', 'O2*',
                'O6', 'O8', 'O9'),
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
        if configuration in products_by_module:
            product_name = products_by_module[configuration] 
        else:
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
    assert np.allclose(tea.IRR, 0.06760145804238414, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74976201.15513338, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52917786.93817723, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 66139140.44031943, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 9396522.56718782, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 108.19366840623299, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 99.99118048363374, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.92300294915781, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.21010076225105967, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693810.9072454, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58785512.91504841, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 223366592.47328365, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34996794.732216835, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 198.6815712355573, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 248.78635422182805, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.887013584818798, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 124.7252091516593, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.6845271894041846, rtol=5e-2)
    assert np.allclose(tea.sales, 126441590.25659768, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82341081.26750697, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207662643.7741959, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -10808374.834058467, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 357.32861206090973, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 373.70358659633587, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.42419116160912, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 46.650867404106165, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.1393406198726412, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301968.5027519, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165736.57206202, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 199177662.8915299, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18104610.221218303, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 236.16185047299098, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.3576079652774, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.538372498597424, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.25632300276929, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.46631835806534216, rtol=5e-2)
    assert np.allclose(tea.sales, 81216703.38995753, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57891735.256732605, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 142809178.71277824, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -20434820.483626157, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 326.68221545803823, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 312.87993682715984, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.533091759502069, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 89.62700162366329, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 142107917.71520165, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69136354.07706703, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 260896517.1648265, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1789246.798722768, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 425.04400181599163, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 447.9257234543171, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 27.431081294067887, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 27.431081294067866, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.46631835806534216, rtol=5e-2)
    assert np.allclose(tea.sales, 73144071.17496133, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59935860.31495893, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 180885656.26972273, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46639564.38836228, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 315.2628578639308, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 304.552918793788, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.198841665581883, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 184.13707727663416, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166847018.64939293, rtol=5e-2)
    assert np.allclose(tea.material_cost, 75062560.00817204, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 282707539.67539394, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5775640.955512564, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 434.15019602019834, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 434.77804385923054, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.36626843678665, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 52.300773536297726, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.46631835806534216, rtol=5e-2)
    assert np.allclose(tea.sales, 60655502.14182142, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57499545.52708925, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 172793812.18859252, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46275594.171999834, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 300.8740870661126, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 288.69707275173494, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.94828130728928, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 182.556742735425, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148592752.4063514, rtol=5e-2)
    assert np.allclose(tea.material_cost, 71529302.36710705, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 272675985.6884632, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5233269.984793013, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 412.99755732620787, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 413.3155522741972, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.086644876688943, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 50.062464358456324, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.46631835806534216, rtol=5e-2)
    assert np.allclose(tea.sales, 51350651.15547196, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57860463.10031674, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 130678243.30031638, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -21243196.374852374, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 161.39909953088167, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 191.01709670499895, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.49431007955507, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 84.05613664897392, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 75607969.19947377, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69279559.74404275, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 235270542.61037537, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 963277.4683855053, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 158.24664866227837, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 253.11204911829913, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.996455453589206, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 33.056035817904565, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.46631835806534216, rtol=5e-2)
    assert np.allclose(tea.sales, 47691059.515063815, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59448621.04083655, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 162771827.82858622, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47902363.12590465, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 191.99609067360586, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 209.50092558012324, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.240031443890345, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 160.5410569819557, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 89586731.34953079, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74780043.05217299, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 254949979.03274456, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6793681.631487491, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 216.01960832043798, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 268.14114349552824, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.716720856403903, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 54.083622604833955, rtol=5e-2)
    
@default_settings
def test_oilcane_O6():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O6')
    else:
        try: module.load('O6')
        except: pass
    feedstock = module.feedstock
    product = module.biodiesel
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 1.38, rtol=5e-2)
    assert np.allclose(tea.sales, 94628961.54871106, rtol=5e-2)
    assert np.allclose(tea.material_cost, 70890745.90760629, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 496884487.4959105, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1399367.8051586035, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 694.942870585679, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 933.9772836749571, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 54.51411990797513, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 54.53645321235478, rtol=5e-2)
    
@default_settings
def test_oilcane_O8():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O8')
    else:
        try: module.load('O8')
        except: pass
    feedstock = module.feedstock
    product = module.biodiesel
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 1.38, rtol=5e-2)
    assert np.allclose(tea.sales, 87099699.09221663, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74843010.11414064, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 514460248.8393247, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -32968590.870251596, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 206.2126745565225, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 704.3541527018425, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 57.55140299117374, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 181.68763253446355, rtol=5e-2)
    
@default_settings
def test_oilcane_O9():
    from biorefineries import oilcane as module
    if 'oilcane' in must_load:
        module.load('O9')
    else:
        try: module.load('O9')
        except: pass
    feedstock = module.feedstock
    product = module.biodiesel
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03500000000000001, rtol=5e-2)
    assert np.allclose(product.price, 1.38, rtol=5e-2)
    assert np.allclose(tea.sales, 97231141.2305434, rtol=5e-2)
    assert np.allclose(tea.material_cost, 75361133.95928662, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 498371274.83804065, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -30226850.602626987, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 224.07782259173135, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 678.5283881633434, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 55.480703569916635, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 169.71567485527137, rtol=5e-2)
  
# TODO: Potentially deprecate biorefinery
# @default_settings
# def test_LAOs():
#     from biorefineries import LAOs as module
#     if 'LAOs' in must_load:
#         module.load()
#     else:
#         try: module.load()
#         except: pass
#     feedstock = module.glucose
#     product = module.octene
#     tea = module.LAOs_tea
#     units = UnitGroup('Biorefinery', tea.units)
#     assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
#     assert np.allclose(feedstock.price, 0.265, rtol=5e-2)
#     assert np.allclose(product.price, 1.2759748416953078, rtol=5e-2)
#     assert np.allclose(tea.sales, 159370478.7059679, rtol=5e-2)
#     assert np.allclose(tea.material_cost, 135665282.30809242, rtol=5e-2)
#     assert np.allclose(tea.installed_equipment_cost, 54579433.65324043, rtol=5e-2)
#     assert np.allclose(tea.utility_cost, 2728013.015676101, rtol=5e-2)
#     assert np.allclose(units.get_heating_duty(), 38.28531194772121, rtol=5e-2)
#     assert np.allclose(units.get_cooling_duty(), 123.59847080242945, rtol=5e-2)
#     assert np.allclose(units.get_electricity_consumption(), 3.416233374810637, rtol=5e-2)
#     assert np.allclose(units.get_electricity_production(), 3.4162333748106346, rtol=5e-2)

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
    # test_oilcane_O6()
    # test_oilcane_O8()
    # test_oilcane_O9()
    # test_oilcane_S1_agile()
    # test_oilcane_S2_agile()
    # test_oilcane_O1_agile()
    # test_oilcane_O2_agile()
    # test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





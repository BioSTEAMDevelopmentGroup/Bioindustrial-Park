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
marked_slow = {
    'wheatstraw', 'animal_bedding'
}
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
    assert np.allclose(tea.IRR, 0.05381572438333119, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74992038.8556337, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52925382.42511212, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 65387501.9390586, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 10442066.327145828, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 111.3588815958425, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 105.08735576937576, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.6248530112233648, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.2096184874236113, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693810.90724541, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58791909.359143384, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 224993001.08352256, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -35051662.1479298, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 196.03543411890655, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 248.74360483108703, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.335938857939432, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 125.34999152642408, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.6865466977048817, rtol=5e-2)
    assert np.allclose(tea.sales, 126814634.63623773, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82270811.86953875, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 209112840.52373964, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -10675099.091695543, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 351.9303289817331, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 347.548288866479, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.05578511575179, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.93524983359892, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.13927430031288052, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301968.50275187, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57173796.168331705, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 201382960.99778348, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18258007.571078993, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 231.71520460357308, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.1487459328047, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.09763804588573, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 72.3062254997152, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 81216703.38995746, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57899806.22775801, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 141266073.40874788, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -16698257.182554739, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 381.2012571371377, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 303.72883725099325, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 11.155375956039892, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 76.75445011609807, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 142237746.46385777, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69156446.04847895, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 262214925.83913735, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1752477.5030571576, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 419.21135806535466, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 440.1365785915077, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 28.23560017347859, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 28.235600173478584, rtol=5e-2)
    
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
    assert np.allclose(tea.material_cost, 59941590.45130447, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 179058936.52354053, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -41571339.24648653, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 390.90054156403824, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 292.68962161516185, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.640667230808956, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 166.27564789336768, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 166847018.57179162, rtol=5e-2)
    assert np.allclose(tea.material_cost, 75055649.87958117, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 284079794.73630774, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5953751.651936062, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 428.02725600887777, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 425.1743322067999, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 27.169572389185724, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 53.74711851318012, rtol=5e-2)
    
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
    assert np.allclose(tea.material_cost, 57504827.471605845, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 170898761.4138952, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -41203776.35488143, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 376.72383507765727, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 276.8132251809804, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.355542043729026, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 164.64777445492118, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 148592752.84767812, rtol=5e-2)
    assert np.allclose(tea.material_cost, 71521798.12757736, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 273980124.0171497, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5403704.210478699, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 407.1857538126716, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 403.6792764194469, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.844126113125732, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 51.43526492864981, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 51350651.15547192, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57868463.28708232, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 129289997.63716191, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -17446725.727028634, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 212.40574416786785, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 183.9112307198191, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.970399437290968, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 72.87097722785816, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 75648326.35464674, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69294551.60661857, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 236558724.27895167, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 892623.9924263542, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 158.24810943502771, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 247.74161418869107, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.60122291182192, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 34.23013795039522, rtol=5e-2)
    
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
    assert np.allclose(tea.material_cost, 59454284.44353091, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 161079591.31872085, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -42745171.73383796, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 258.4931903315889, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 199.80673709467587, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.577142107885178, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 145.28168621782208, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 89586731.3110852, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74772936.9182138, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 256084574.36616516, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6961298.235638824, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 215.9550286633267, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 260.96213963819093, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.363879199889762, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 55.11560474413598, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 91006175.06693201, rtol=5e-2)
    assert np.allclose(tea.material_cost, 70070199.61688498, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 489186430.8591811, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1396644.3983708634, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 711.9518530355249, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 896.1831596698129, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 54.47379412534113, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 54.50548041056594, rtol=5e-2)
    
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
    assert np.allclose(tea.material_cost, 74828290.96493468, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 515526381.2209381, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -32770585.514114555, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 205.74408141645563, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 690.9484367882108, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 58.37797218308237, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 181.79891812083665, rtol=5e-2)
    
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
    assert np.allclose(tea.material_cost, 75354255.64536168, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 498073100.4135921, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27440503.318503812, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 263.19023480751235, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 663.2754246494304, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 56.30705460772081, rtol=5e-2)
    # assert np.allclose(units.get_electricity_production(), 160.47964930153773, rtol=5e-2)
  
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
    test_oilcane_O6()
    test_oilcane_O8()
    test_oilcane_O9()
    test_oilcane_S1_agile()
    test_oilcane_S2_agile()
    test_oilcane_O1_agile()
    test_oilcane_O2_agile()
    # test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





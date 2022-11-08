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
    'test_wheatstraw',
    'test_animal_bedding',
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
    configuration_tag = f"_{configuration}".replace('*', 'agile') if configuration else ''
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
    assert np.allclose(tea.IRR, 0.05484137947443737, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723599.41753717, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525177.04024876, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 62101158.6567123, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7480503.5096371, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.08307155629652, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.69704367515685, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9369502187665468, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.2113160935470353, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 103632376.79166704, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58858510.67200561, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 206548409.19150496, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -28942675.979259465, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 206.98209150873876, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 217.9662298884939, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.253708655344464, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 104.4447376459022, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7382534422812541, rtol=5e-2)
    assert np.allclose(tea.sales, 136365493.22811958, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82338339.40833226, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 202809309.1641161, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5343615.419122992, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 363.2978269464661, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 310.18429398661016, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.96420317229188, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 33.26974052126595, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.14295807786505277, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606541, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57277034.59309775, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 167497255.24592966, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -12214721.859830968, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.47182596212357, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 222.97316413184618, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.485040401747591, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 51.049620456611784, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 38.72786405515464, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73736511, rtol=5e-2)
    assert np.allclose(tea.material_cost, 61966747928.798935, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 120393463.29938526, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -16619199.456099978, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 246.43319824149506, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 245.01879601587152, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.969452183294071, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 66.46087263831664, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 30.298329358930594, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 121082189.94588801, rtol=5e-2)
    assert np.allclose(tea.material_cost, 48487151275.51035, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 220113258.27978802, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 3040845.7517896923, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 335.77985196845947, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 309.3471136143858, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.605729734873215, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 21.605729734873208, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 41.451503228944574, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71993360.47703359, rtol=5e-2)
    assert np.allclose(tea.material_cost, 66325967353.69603, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 163349892.86799583, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -41321070.27424532, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 252.1405067883693, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 223.61256956220166, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.96236348055846, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 152.65978522998645, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 34.79517738296873, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 136218442.09948644, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55686300242.64717, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 232885996.50894728, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 3073281.8327416694, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 510.09808214059336, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 379.07947530399883, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.60470348782333, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.605524292953255, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 36.28377617348664, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.74167342, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58055156156.09348, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 155595234.6461211, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -41026646.85053678, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 246.9800900477171, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 216.77598585710086, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.731256311182769, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 150.1534601155833, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 27.786609251098497, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 118273243.15192656, rtol=5e-2)
    assert np.allclose(tea.material_cost, 44469169460.17138, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 222699099.13798767, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2539262.603052804, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 498.81963333089385, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 367.0390220362551, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.138812868268865, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.139640400379495, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 42.55755678950908, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.4163145, rtol=5e-2)
    assert np.allclose(tea.material_cost, 88522507090.33224, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 126827608.96577129, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -23052567.399197616, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 116.58223220714964, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 194.4210315724837, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.216905546072425, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 81.18528634869719, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 37.48391110023801, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 158417310.84075725, rtol=5e-2)
    assert np.allclose(tea.material_cost, 77980105097.48112, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 230984689.52941796, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 3945925.036118205, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 124.8325866661483, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 321.28351077655685, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.952207878688043, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.05474196088139, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 45.94366514811367, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 92435868.22021642, rtol=5e-2)
    assert np.allclose(tea.material_cost, 95567254388.9696, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 164805803.40080687, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -51539355.6350529, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 226.99597514720185, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 199.49601258509207, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.840388930547745, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 131.171307045442, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 41.303884431766335, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 171616246.20992085, rtol=5e-2)
    assert np.allclose(tea.material_cost, 85929799051.34785, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 234217725.98683402, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 3940733.851130674, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 442.8708176234132, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 209.20490569748029, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.11534555288025, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.20590538501947, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.2918667227938558, rtol=5e-2)
    assert np.allclose(tea.sales, 161435883.1972842, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.45874006, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 55049531.41964977, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2702388.184876362, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.31078344171506, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.69592718862404, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.527871346094877, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.5278713460948787, rtol=5e-2)
    

### DO NOT DELETE:
### Code commented for legacy purposes
# @pytest.mark.slow
# def test_lactic():
#     from biorefineries import lactic
#     bst.process_tools.default_utilities()
#     system = lactic.system
#     MPSP = system.lactic_acid.price
#     assert np.allclose(MPSP, 1.5285811040727408, atol=0.01)
#     tea = system.lactic_tea
#     assert np.allclose(tea.sales, 332247972.3322271, rtol=0.01)
#     assert np.allclose(tea.material_cost, 218911499.70700422, rtol=0.01)
#     assert np.allclose(tea.installed_equipment_cost, 321744540.84435904, rtol=0.01)
#     assert np.allclose(tea.utility_cost, 24995492.58140952, rtol=0.01)
#     units = UnitGroup('Biorefinery', system.lactic_tea.units)
#     assert np.allclose(system.CHP.system_heating_demand/1e6, 1763.4942302374052, rtol=0.01)
#     assert np.allclose(-system.CT.system_cooling_water_duty/1e6, 1629.0964343101657, rtol=0.01)
#     assert np.allclose(units.get_electricity_consumption(), 45.291535445041546, rtol=0.01)
#     assert np.allclose(units.get_electricity_production(), 0.0)

# @pytest.mark.slow
# def test_ethanol_adipic():
#     bst.process_tools.default_utilities()
#     from biorefineries import ethanol_adipic
#     acid = ethanol_adipic.system_acid
#     MESP = acid.ethanol.price * acid._ethanol_kg_2_gal
#     assert np.allclose(MESP, 2.486866594017598, atol=0.01)
#     tea = acid.ethanol_tea
#     assert np.allclose(tea.sales, 152002785.9542236, rtol=0.01)
#     assert np.allclose(tea.material_cost, 112973059.93909653, rtol=0.01)
#     assert np.allclose(tea.installed_equipment_cost, 202880204.7890376, rtol=0.01)
#     assert np.allclose(tea.utility_cost, -18225095.0315181, rtol=0.01)
#     units = UnitGroup('Biorefinery', acid.ethanol_tea.units)
#     assert np.allclose(acid.CHP.system_heating_demand/1e6, 333.8802418517945, rtol=0.01)
#     assert np.allclose(units.get_cooling_duty(), 327.05990128304114, rtol=0.01)
#     assert np.allclose(units.get_electricity_consumption(), 23.84999102548279, rtol=0.01)
#     assert np.allclose(units.get_electricity_production(), 55.72024685271333, rtol=0.01)
    
#     base = ethanol_adipic.system_base
#     MESP = base.ethanol.price * base._ethanol_kg_2_gal
#     assert np.allclose(MESP, 2.703699730342356, atol=0.01)
#     tea = base.ethanol_adipic_tea
#     assert np.allclose(tea.sales, 219850046.2308626, rtol=0.01)
#     assert np.allclose(tea.material_cost, 120551649.9082640, rtol=0.01)
#     assert np.allclose(tea.installed_equipment_cost, 283345546.0556234, rtol=0.01)
#     assert np.allclose(tea.utility_cost, 14241964.663629433, rtol=0.01)
#     units = UnitGroup('Biorefinery', base.ethanol_adipic_tea.units)
#     assert np.allclose(base.CHP.system_heating_demand/1e6, 296.8322623939439, rtol=0.01)
#     assert np.allclose(units.get_cooling_duty(), 247.79573940245342, rtol=0.01)
#     assert np.allclose(units.get_electricity_consumption(), 26.565278642577344, rtol=0.001)
#     assert np.allclose(units.get_electricity_production(), 0.0)
    
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





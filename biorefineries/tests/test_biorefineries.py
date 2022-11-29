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
    assert np.allclose(tea.IRR, 0.05484137947443766, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723599.41753717, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525177.04024876, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 62101158.65671231, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7480503.509637097, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.08307155629652, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.69704367515685, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9369502187665466, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.2113160935470352, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 103632376.79166704, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58858510.67200561, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 206548409.19150504, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -28942675.979259465, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 206.98209150873885, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 217.96622988849387, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.253708655344462, rtol=5e-2)
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
    assert np.allclose(product.price, 0.7382534422848163, rtol=5e-2)
    assert np.allclose(tea.sales, 136365493.22877756, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82338339.40833226, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 202809309.16411608, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5343615.419122993, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.14295807786505282, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606541, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57277034.59309775, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 167497255.24592963, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -12214721.859830972, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.47182596212357, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 222.97316413184618, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.48504040174759, rtol=5e-2)
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
    assert np.allclose(feedstock.price, 0.038727864055154645, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73736511, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63969076.23320489, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 120393463.2993852, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -16619199.45609998, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 246.433198241495, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 245.01879601587152, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.969452183294074, rtol=5e-2)
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
    assert np.allclose(feedstock.price, 0.03132168467589462, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 121082189.94588801, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59813085.95090494, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 220113258.27978796, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1403472.9847080866, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 335.77985196845947, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 309.34711361438576, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.60572973487322, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 21.60572973487318, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03132168467589462, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 121082189.94588801, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59813085.95090494, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 220113258.27978796, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1403472.9847080866, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 335.77985196845947, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 309.34711361438576, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.60572973487322, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 21.60572973487318, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0358466818852938, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 136218442.09948644, rtol=5e-2)
    assert np.allclose(tea.material_cost, 71371009.68382457, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 232885996.46071145, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1390874.637524818, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 510.0980824685457, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 379.0794751805062, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.60470348737718, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.60552421635094, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.036976650892765514, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.74167342, rtol=5e-2)
    assert np.allclose(tea.material_cost, 60276803.95352954, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 155595234.64612114, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -42135324.82550437, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 246.98009004771708, rtol=5e-2)
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
    assert np.allclose(feedstock.price, 0.028503844547222545, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 118273243.2145202, rtol=5e-2)
    assert np.allclose(tea.material_cost, 56200720.874732584, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 222700495.9676085, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1391380.372901436, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 498.81963500075017, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 367.0390240764973, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.1388128814724, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.139640553259618, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04255755678950908, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.4163145, rtol=5e-2)
    assert np.allclose(tea.material_cost, 91078765.15968868, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 126827608.96577129, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -23052567.399197616, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 116.58223220714966, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 194.4210315724837, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.216905546072427, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 81.1852863486972, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03850594196847115, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 158417310.84075725, rtol=5e-2)
    assert np.allclose(tea.material_cost, 93459863.52950847, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 230984688.76621056, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1820018.8583263834, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 124.83259313914101, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 321.2835104244336, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.952207880946357, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.054739676397364, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03850594196847115, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 158417310.84075725, rtol=5e-2)
    assert np.allclose(tea.material_cost, 93459863.52950847, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 230984688.76621056, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1820018.8583263834, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 124.83259313914101, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 321.2835104244336, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.952207880946357, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.054739676397364, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04232874996832941, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 171615833.18740612, rtol=5e-2)
    assert np.allclose(tea.material_cost, 105762996.58664085, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 234216693.8860991, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1808931.7663218838, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 442.8693228566128, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 209.2022798714919, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.114981136097715, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 19.20553567294961, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.291866722793068, rtol=5e-2)
    assert np.allclose(tea.sales, 161435883.19718194, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.45874006, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 55049531.419649765, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2702388.1848763623, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.31078344171506, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.69592718862403, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.527871346094877, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.5278713460948787, rtol=5e-2)

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





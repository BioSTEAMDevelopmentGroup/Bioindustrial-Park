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
    assert np.allclose(tea.IRR, 0.06139011178456174, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74973086.74735086, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52917743.675847545, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 66083566.37324064, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 9832127.768854192, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 108.17555454164444, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 99.98403221363675, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9224686061832343, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.2100514791090884, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693808.63788584, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798587.08599357, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 223381724.61662412, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34996723.807722345, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 198.6815430225962, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 248.7863469967366, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.887247406054514, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 124.72521565079774, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.6934706653247219, rtol=5e-2)
    assert np.allclose(tea.sales, 128093570.98172233, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82341315.08800615, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 208126596.35525888, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -9263997.915747166, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 357.3426519926122, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 373.71015228705033, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.55632362881821, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 46.64586595970775, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.1355824716933896, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301963.95367965, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165736.47820577, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 200225766.39200342, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -17297876.692693803, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 236.16184775838, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.35760281712936, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.66383246735897, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.25632331818251, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 81216699.20927736, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57891735.16294234, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 143508135.6921627, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19627657.380752247, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 326.6822677645666, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 312.87994069793905, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.672278059659584, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 89.62698887108415, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 142059925.14741462, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69128486.24747512, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 261982005.79193866, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 3209403.528556343, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 424.9785995593406, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 447.9872698058816, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 27.688761907631022, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 27.68876190763101, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 73144066.99332994, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59935860.21922503, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 181244807.48321533, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46268517.36951378, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 315.262837698138, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 304.5529019506643, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.257724887119764, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 184.13708577468518, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 166847023.86146083, rtol=5e-2)
    assert np.allclose(tea.material_cost, 75062553.54596199, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 283349599.34683865, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4874170.121213364, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 434.1563584047804, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 434.7837866029941, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.525151992910974, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 52.29949682506389, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 60655497.95723747, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57499545.43135535, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 173126792.27919227, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45924176.34956955, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 300.8740843655788, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 288.6970694933892, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 10.003610257116053, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 182.5567469543418, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 148592757.58212915, rtol=5e-2)
    assert np.allclose(tea.material_cost, 71529295.90505801, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 273303528.5354287, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4337795.373403454, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 413.003745238578, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 413.32131418871904, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.244259236506807, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 50.061181500690715, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 51350648.51161736, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57860463.00801965, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 131267580.97376747, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -20441929.16990094, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 161.39909332406498, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 191.01709389268333, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.596679694677176, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 84.05613891215202, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 75594412.32655069, rtol=5e-2)
    assert np.allclose(tea.material_cost, 69298501.16793385, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 236089599.175185, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2389859.6875048475, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 158.22772476614065, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 252.97866825252206, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.198114044486562, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 33.06590886969057, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 47691056.89102804, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59448620.95333576, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 163062525.36518234, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47542308.99393952, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 191.99607945214262, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 209.50091536252026, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.279767833445963, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 160.54106385563222, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 89586729.29048933, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74620252.77562289, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 255455299.66366962, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5901813.933428848, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 215.93822771295146, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 268.04801439597844, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.83101691944999, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 54.020066981015546, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 101050580.35400365, rtol=5e-2)
    assert np.allclose(tea.material_cost, 72376570.1729823, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 498944271.89657515, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2287427.4494785187, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 720.4745232533116, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 769.1747182006673, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 41.38615239011329, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 41.386003430143575, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 87101702.1865109, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74842849.5636859, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 505194260.6006373, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -35834591.49779868, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 206.21271908650357, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 528.1363391392526, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 44.0826728281144, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 181.687618110022, rtol=5e-2)
    
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
    assert np.allclose(tea.sales, 97231171.8543844, rtol=5e-2)
    assert np.allclose(tea.material_cost, 75361130.82022867, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 489053572.2963267, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -33138903.9061295, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 224.0778059274724, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 502.3155362496767, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 42.00421898178156, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 169.7158558567133, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.2759748416953078, rtol=5e-2)
    assert np.allclose(tea.sales, 159370478.7059679, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135665282.30809242, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 54579433.65324043, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2728013.015676101, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.28531194772121, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.59847080242945, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.416233374810637, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.4162333748106346, rtol=5e-2)
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





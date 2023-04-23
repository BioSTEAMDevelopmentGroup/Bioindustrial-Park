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
    assert np.allclose(tea.IRR, 0.04922733123672491, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74995091.73070191, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52916786.73434014, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 65442127.70623944, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 10271412.679106986, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 115.83227253969748, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 105.40843363176639, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.0345530068091895, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.2087969269572819, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.2272782, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798824.46333956, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 222683054.95395762, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34302930.87741131, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 208.07987555959306, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 247.6140635537415, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.8918931386888564, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 122.50616491289568, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7237629470132995, rtol=5e-2)
    assert np.allclose(tea.sales, 133688954.56797652, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82366246.60839348, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207023090.2182626, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -9204629.753290413, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 373.0930414985088, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 382.0337968578203, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.430868103935907, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 43.343583217122216, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.13595611785517797, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302441.90067281, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165746.3360821, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 197174697.00375932, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -17147611.842259977, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 249.55812694712847, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 270.15238349973663, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.390086345637947, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 68.09331295417243, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038668493418044764, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211568.14673989, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63757445.55934654, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 132738518.64029773, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19453007.22274715, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 307.61461856303066, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 276.334892978732, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.389404466621892, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 77.44292647858074, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03734261878205841, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 140855458.63127875, rtol=5e-2)
    assert np.allclose(tea.material_cost, 72668315.57603355, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 243061542.15145144, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2064625.4580368823, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 400.039204198183, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 412.28629382918893, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.731581625448484, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 24.731581625448484, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04432004121577088, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 72469014.54514585, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74375341.8766146, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169361146.0206496, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46882143.62274073, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 290.89365338183285, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 272.30767767459804, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.20842824244405, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 166.6487062719999, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05030374891288528, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166772454.0665997, rtol=5e-2)
    assert np.allclose(tea.material_cost, 99174148.6776965, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 261746956.4511415, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -3738175.4723080383, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 423.16465571245186, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 410.84020153967833, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.921463959420592, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 40.60483252744789, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03886752601674389, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59978850.91960837, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63205163.394760385, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 161638861.81801897, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46132782.10481969, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 277.39099435186506, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 257.86433292207556, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.990220688788664, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 163.98471863466187, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.042152024053291576, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148517751.4527836, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82588282.31889394, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 251929035.68585032, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2798272.972854699, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 403.3528103503903, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 389.9418387747109, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.61282390834223, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 37.24134519005745, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04292627955792694, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103452130.03145395, rtol=5e-2)
    assert np.allclose(tea.material_cost, 91693675.39573096, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 139126741.86755073, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -26713297.16668028, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 188.16448699669445, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 208.69991826034314, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.363263125746057, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 91.66935352198716, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.044866356484318796, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 180935384.3825863, rtol=5e-2)
    assert np.allclose(tea.material_cost, 110566361.79525942, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 252738901.5434872, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1550914.6086441134, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 177.52177297304314, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 282.81134637074655, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.859119756908314, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 34.31171959228939, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04920561573311784, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 93007865.23945124, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106662762.23774314, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 170791961.56928682, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -59158639.98942079, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 241.52128798369225, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 231.3522851931485, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.11083720474494, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 148.94304450369347, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.057627746676254525, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 210974460.47876817, rtol=5e-2)
    assert np.allclose(tea.material_cost, 143805834.36405486, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 265869645.3561265, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5457709.305558443, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 232.80763313548175, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 295.26877881935224, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.785148148293988, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 46.6937885512565, rtol=5e-2)
    
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
    # test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





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
    assert np.allclose(tea.IRR, 0.054626557115021235, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74973086.74735086, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52917743.675847545, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 66010546.39191148, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 9844797.713497443, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 108.4879126317954, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 99.98403221363675, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9167496359672997, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.210033140138807, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102693808.63788584, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798573.08879176, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 223284755.50565106, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34957643.56858203, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 199.24708293383873, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 248.7015160100502, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.878974325965104, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 124.59168539397659, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7134129336617239, rtol=5e-2)
    assert np.allclose(tea.sales, 131777182.26117879, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82361657.3298363, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207640710.5095279, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -10805703.676522434, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 357.3426519926122, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 373.71015228705033, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.423770836398187, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.13778334033156747, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301963.95367965, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165736.47820577, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 199044415.61415282, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18047819.17345761, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 236.74571056830356, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.27002339564086, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.529859167636321, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.118466821395, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03909901727316123, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81216602.88501228, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64450619.78219299, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 133887008.04192507, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -20353183.679144043, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 294.67014595879596, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 281.49342160864705, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.529714547983003, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 80.50937631034155, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03825450437104777, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 142100832.9900629, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74343245.50635117, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 242667293.76864335, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1791158.6037448407, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 382.530952123583, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 403.14308296698505, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.725525568593074, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 24.725525568593078, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04422371475194139, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 73143973.81859528, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74826205.21284597, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169303860.25186098, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46607710.98645762, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 284.4099271585844, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 262.7160690133169, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.995172211986278, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 165.5643830195706, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05104468966718037, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166847026.9192766, rtol=5e-2)
    assert np.allclose(tea.material_cost, 100957912.13914822, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 263523879.2391125, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5706788.952956749, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 391.56554900277195, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 391.1816036293526, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.757353248618163, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 46.874670755096254, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.039054528365420914, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 60656252.42314774, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64040597.04032007, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 161592015.69687632, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -46244696.21561435, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 271.4574054713763, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 248.4455710627821, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.768248235748286, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 164.1426989840477, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04317128413931401, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148592757.5821292, rtol=5e-2)
    assert np.allclose(tea.material_cost, 84712639.73620406, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 253967485.56708425, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5164898.424652384, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 372.52819769506664, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 371.8653892457309, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.504141277823138, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 44.860203229794976, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.040953135166192464, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 79579640.13538662, rtol=5e-2)
    assert np.allclose(tea.material_cost, 67377870.66513433, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 119298990.81935692, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -21355194.3325402, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 140.16490588853017, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 165.2141996364127, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.41801746589942, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 72.60076874325792, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04201665038225913, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 140087618.21769792, rtol=5e-2)
    assert np.allclose(tea.material_cost, 80541942.36328456, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 214336296.34728137, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 760888.9257047343, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 136.9442313222078, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 219.0534211812585, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 19.08693426416488, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 28.628947505856978, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04633259624749116, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71543984.5715841, rtol=5e-2)
    assert np.allclose(tea.material_cost, 77482467.27363421, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 146272868.74803126, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45481725.3566978, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 192.8876553959246, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 183.80909903551674, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.093248077100177, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 113.06769072621711, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05505235072315981, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 162237111.92611295, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106517176.26915126, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 226216492.7800231, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -6388793.126066338, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 167.16303440846175, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 220.74654715880496, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.234308730123544, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 41.57029379972555, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.288805979721651, rtol=5e-2)
    assert np.allclose(tea.sales, 161035503.37609148, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135665282.30809242, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 54579433.653240435, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2728013.015676101, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.28531194772121, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.59847080242946, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.4162333748106373, rtol=5e-2)
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





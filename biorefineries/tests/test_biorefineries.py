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
    assert np.allclose(tea.IRR, 0.04922733123672505, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74995091.73070191, rtol=5e-2)
    assert np.allclose(tea.material_cost, 52916786.73434014, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 65442127.706239425, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 10271412.67910698, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 115.83227253969739, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 105.40843363176639, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.034553006809189, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.20932919399877384, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.2272782, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58798737.910276994, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 222740280.64491123, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -34517330.76179105, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 205.17342749647366, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 247.08985731746512, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.890958703657048, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 123.19240959446556, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7145111617876865, rtol=5e-2)
    assert np.allclose(tea.sales, 131980022.79714778, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82366246.60839348, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207640203.212064, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -11041041.106571833, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 356.99137288069835, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 373.13690860863034, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.4272660859673, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.145366085216345, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.1372113082983206, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302441.90067281, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57165746.3360821, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 197624315.82163516, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -17578015.941603705, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 243.72708389992727, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 269.1006983378128, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.38858713981407, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 69.47008700698386, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038906765965980375, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211568.14673989, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64138682.62725732, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 132957617.30555125, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19877457.122838948, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 301.7835763258438, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 275.28320762945066, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.387905260791424, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 78.81970034013877, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038091083249825616, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 141782972.35415757, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74019579.54032566, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 242058951.96357125, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1856458.8437069105, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 381.73423549500194, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 401.4358937347683, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.79328178494203, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 24.793281784942018, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04381903933083459, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 72541162.23399265, rtol=5e-2)
    assert np.allclose(tea.material_cost, 73561169.3065002, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 167904676.72549188, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45667863.91251517, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 286.7340545557082, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 270.24956697588044, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.240299491405967, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 162.72194363407567, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05100049937427724, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 165471446.99165088, rtol=5e-2)
    assert np.allclose(tea.material_cost, 100053205.3827274, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 262103651.71402207, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5920428.319076779, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 404.00963416799794, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 401.9419701212242, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.74726211930081, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.52300168041126, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03870430795144252, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 60154361.83981891, rtol=5e-2)
    assert np.allclose(tea.material_cost, 62950801.382851556, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 160424575.5450812, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -45428439.93708473, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 273.3424982315873, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 255.0361356306883, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.942650297243977, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 161.63014291518465, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04227783727378679, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 147367587.0345435, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82582239.42326018, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 250967289.93424806, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -3658839.0076199747, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 384.3608304399581, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 378.47041113003166, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.4408065563754, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 39.86604520866369, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.043169694393018876, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103452130.03146541, rtol=5e-2)
    assert np.allclose(tea.material_cost, 92199979.56910989, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 139305529.37205353, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27253989.037229188, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 188.48021056474133, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 209.5042241589957, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.36186124488584, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 92.93539035660477, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04550810337704206, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 181788090.91717443, rtol=5e-2)
    assert np.allclose(tea.material_cost, 111925631.58243252, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 252483405.54725802, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1125039.47345165, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 177.93178804843694, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 284.144190858698, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.830707961559863, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 36.4696954915187, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04908419340472239, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 93007865.23945123, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106410203.79417554, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 170453639.44925314, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -58827900.327230155, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 241.7199349035914, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 233.27185997786694, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.20796433247521, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 148.39641891205272, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0589118930242581, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 210974460.47876817, rtol=5e-2)
    assert np.allclose(tea.material_cost, 146476858.77324435, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 267760597.61769038, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -8520343.481336148, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 233.19182278523107, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 300.35939091803255, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.782223104716127, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 54.70678907695694, rtol=5e-2)
    
    
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





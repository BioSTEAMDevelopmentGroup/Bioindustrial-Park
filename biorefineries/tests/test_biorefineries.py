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
    assert np.allclose(tea.IRR, 0.055241200299056194, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723596.88174266, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525176.94552829, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 61891581.16887585, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7484262.226815703, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.0630956454741, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.83385198664736, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.936953423657254, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.21478565883250303, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.74923621, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58864075.39522756, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 224182651.6582225, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -37089258.15332504, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 207.16520420513572, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 235.01806569112503, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.556464498844547, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 130.858397474069, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7050642372630093, rtol=5e-2)
    assert np.allclose(tea.sales, 130235000.4471838, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82263697.41365018, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 205816655.19519845, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -12215300.052522901, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 356.2743693558788, rtol=5e-2)
    # TODO: Possibly CI fails in this test because of changes to chemicals and thermo.
    # Double check updates to these packages and update thermosteam dependencies.
    # assert np.allclose(units.get_cooling_duty(), 298.3267172215581, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.623806970125774, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.214345655327804, rtol=5e-2)
    
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
    assert np.allclose(tea.IRR, 0.14081844241721223, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606188, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57282511.78382494, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 197049762.1894435, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -18519371.790406622, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.70632887653446, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 236.39928468248982, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.725475811227108, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 71.49731652501673, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.0389383348413644, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73735796, rtol=5e-2)
    assert np.allclose(tea.material_cost, 64311172.42803233, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 132288786.67392789, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -19884394.428716373, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 304.43869244437116, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 252.84773285811895, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.004153957765865, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 78.10794442413409, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03573458535189133, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 134003547.72568287, rtol=5e-2)
    assert np.allclose(tea.material_cost, 68973047.60336359, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 233483285.20677984, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1403797.1886558784, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 340.1625680779562, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 330.32183815957416, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.349094123165134, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.349094123165134, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04449688209003297, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71993360.48574364, rtol=5e-2)
    assert np.allclose(tea.material_cost, 74764522.34697373, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 169925319.8530473, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47865574.496539906, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 254.375961600393, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 236.33352439267088, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.16721512244948, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 169.2135261050696, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.050600185897238026, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 166664496.79685155, rtol=5e-2)
    assert np.allclose(tea.material_cost, 99692586.37760958, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 260577732.39089537, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -4057102.602630536, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 373.79658082208766, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 351.76465645262596, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.60835447217956, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 41.31363544392717, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.039063223808054215, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.761615224, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63622739.119250365, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 162132018.40067413, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -47122581.99247755, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 249.21647989277696, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 229.40765207586395, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.934519524150376, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 166.56599883999004, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04261805307805807, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 148409512.41547424, rtol=5e-2)
    assert np.allclose(tea.material_cost, 83378068.60324614, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 250625185.26240283, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -3354024.0923702535, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 362.7628389923203, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 339.24224657010865, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.278695051899295, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 38.69883752539266, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.043188920449417585, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.4161051, rtol=5e-2)
    assert np.allclose(tea.material_cost, 92399107.23331147, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 139062816.75508913, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27382839.066928346, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 178.7990581961187, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 202.8250635377964, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.307745965637384, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 93.74371767378848, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.043201700396196294, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 174202861.4871691, rtol=5e-2)
    assert np.allclose(tea.material_cost, 105867997.51442836, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 245081643.52989417, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1341452.264877621, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 141.9291222568288, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 320.1968643406256, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 24.344523910275065, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 29.247420693280652, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04898209744023252, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 92435850.59753108, rtol=5e-2)
    assert np.allclose(tea.material_cost, 106327889.90350282, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 171836794.78330606, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -59628498.16921598, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 232.52044044596025, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 205.00969618992272, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.057779566048495, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 149.55887714829697, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05791391024738966, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 210830375.2851004, rtol=5e-2)
    assert np.allclose(tea.material_cost, 144436222.17026004, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 264545042.7121095, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5887944.203758875, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 291.62713375951375, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 276.84268935860524, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.452776854391285, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 47.5303540096168, rtol=5e-2)
    
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





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
    'test_wheatstraw',
    'test_animal_bedding',
    'test_lactic',
    'test_ethanol_adipic',
    'generate_code',
    'print_results',
)

# Block speed-up for consistent testing. Note that speed-up wrapps functions,
# thus preventing code coverage to be analyzed
bst.speed_up = tmo.speed_up = flx.speed_up = lambda: None 

feedstocks_by_module = {
    'lipidcane': 'lipidcane',
    'cornstover': 'cornstover',
    'sugarcane': 'sugarcane',
    'LAOs': 'glucose',
    'wheatstraw': 'wheatstraw',
    'animal_bedding': 'bedding',
    'lactic': 'feedstock',
    'ethanol_adipic': 'feedstock',
}
products_by_module = {
    'lipidcane': 'ethanol',
    'cornstover': 'ethanol',
    'sugarcane': 'ethanol',
    'LAOs': 'octene',
    'wheatstraw': 'ethanol',
    'animal_bedding': 'ethanol',
    'lactic': 'lactic_acid',
    'ethanol_adipic': 'ethanol',
}

def generate_code(module_name, feedstock_name=None, product_name=None):
    if not feedstock_name:
        feedstock_name = feedstocks_by_module[module_name]
    if not product_name:
        product_name = products_by_module[module_name]
    bst.process_tools.default_utilities()
    module = import_module('biorefineries.' + module_name)
    try: module.load()
    except: pass
    feedstock = getattr(module, feedstock_name)
    product = getattr(module, product_name)
    tea_name = f'{module_name}_tea'
    try:
        tea = getattr(module, tea_name)
    except AttributeError:
        try:
            tea_name = f'{feedstock_name}_tea'
            tea = getattr(module, tea_name)
        except AttributeError:
            tea_name = f'{product_name}_tea'
            tea = getattr(module, tea_name)
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
    print(f"""def test_{module_name}():
    bst.process_tools.default_utilities()
    from biorefineries import {module_name} as module
    try: module.load()
    except: pass
    feedstock = module.{feedstock_name}
    product = module.{product_name}
    tea = module.{tea_name}
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, {IRR}, rtol=1e-2)
    assert np.allclose(feedstock.price, {feedstock.price}, rtol=1e-2)
    assert np.allclose(product.price, {product.price}, rtol=1e-2)
    assert np.allclose(tea.sales, {sales}, rtol=1e-2)
    assert np.allclose(tea.material_cost, {material_cost}, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, {installed_equipment_cost}, rtol=1e-2)
    assert np.allclose(tea.utility_cost, {utility_cost}, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), {heating_duty}, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), {cooling_duty}, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), {electricity_consumption}, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), {electricity_production}, rtol=1e-2)
    """
    )

def generate_all_code():
    from warnings import filterwarnings
    filterwarnings('ignore')
    for i in feedstocks_by_module: generate_code(i)
        
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
    
def test_lipidcane():
    bst.process_tools.default_utilities()
    from biorefineries import lipidcane as module
    try: module.load()
    except: pass
    feedstock = module.lipidcane
    product = module.ethanol
    tea = module.lipidcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.18000945341182453, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 102493171.42530276, rtol=1e-2)
    assert np.allclose(tea.material_cost, 61819200.176152, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 202849847.56774512, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -21065582.11058501, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 276.94917232588114, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 245.98384752626487, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.581267484994973, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 77.09915886507521, rtol=1e-2)
    
def test_cornstover():
    bst.process_tools.default_utilities()
    from biorefineries import cornstover as module
    try: module.load()
    except: pass
    feedstock = module.cornstover
    product = module.ethanol
    tea = module.cornstover_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.05158816935126135, rtol=1e-2)
    assert np.allclose(product.price, 0.7140906601434479, rtol=1e-2)
    assert np.allclose(tea.sales, 134376292.22291023, rtol=1e-2)
    assert np.allclose(tea.material_cost, 80366427.42525397, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 218887698.8347628, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -8369523.722273722, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 341.6190758403897, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 384.3630438280598, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.854606165361538, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 41.25381002701907, rtol=1e-2)
    
def test_sugarcane():
    bst.process_tools.default_utilities()
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.sugarcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.10713543101350805, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 87955077.0538273, rtol=1e-2)
    assert np.allclose(tea.material_cost, 60192842.50809233, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 163303973.8430249, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -6446162.4286760185, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 312.98253582785566, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 250.2056638160375, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.630688343934322, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 30.291465358921567, rtol=1e-2)
    
def test_LAOs():
    bst.process_tools.default_utilities()
    from biorefineries import LAOs as module
    try: module.load()
    except: pass
    feedstock = module.glucose
    product = module.octene
    tea = module.LAOs_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.265, rtol=1e-2)
    assert np.allclose(product.price, 1.3099409534281277, rtol=1e-2)
    assert np.allclose(tea.sales, 163427716.3748871, rtol=1e-2)
    assert np.allclose(tea.material_cost, 135661583.58974993, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 62338869.2109012, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 2788634.0182542275, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 40.00054620871384, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 125.31693495458933, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.596591914799831, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 3.5965919147998306, rtol=1e-2)
    
def test_wheatstraw():
    bst.process_tools.default_utilities()
    from biorefineries import wheatstraw as module
    try: module.load()
    except: pass
    feedstock = module.wheatstraw
    product = module.ethanol
    tea = module.wheatstraw_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.056999999999999995, rtol=1e-2)
    assert np.allclose(product.price, 0.9111601864149698, rtol=1e-2)
    assert np.allclose(tea.sales, 128280893.19815566, rtol=1e-2)
    assert np.allclose(tea.material_cost, 63430743.5116114, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 237017812.01782107, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -5358091.609591275, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 214.26465184757853, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 261.6818186290857, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.954763387282593, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 34.09357263357368, rtol=1e-2)
    
def test_animal_bedding():
    bst.process_tools.default_utilities()
    from biorefineries import animal_bedding as module
    try: module.load()
    except: pass
    feedstock = module.bedding
    product = module.ethanol
    tea = module.bedding_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.003, rtol=1e-2)
    assert np.allclose(product.price, 0.786931542947226, rtol=1e-2)
    assert np.allclose(tea.sales, 113067734.31494412, rtol=1e-2)
    assert np.allclose(tea.material_cost, 33147675.58326452, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 279518315.8576579, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -2040299.6574320383, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 157.8854022174697, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 198.23736472092472, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 32.038638068681244, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 36.280168513723815, rtol=1e-2)
    
def test_lactic():
    bst.process_tools.default_utilities()
    from biorefineries import lactic as module
    try: module.load()
    except: pass
    feedstock = module.feedstock
    product = module.lactic_acid
    tea = module.lactic_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.06287583717512708, rtol=1e-2)
    assert np.allclose(product.price, 1.5290315668201484, rtol=1e-2)
    assert np.allclose(tea.sales, 332346881.2497965, rtol=1e-2)
    assert np.allclose(tea.material_cost, 219074528.09864402, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 321333077.7058257, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 25050313.99229334, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 1829.09043477846, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 1895.6439958358949, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 45.390871189920496, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    
def test_ethanol_adipic():
    bst.process_tools.default_utilities()
    from biorefineries import ethanol_adipic as module
    try: module.load()
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.ethanol_adipic_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.06287581915485817, rtol=1e-2)
    assert np.allclose(product.price, 0.8360250946567223, rtol=1e-2)
    assert np.allclose(tea.sales, 219891872.86398175, rtol=1e-2)
    assert np.allclose(tea.material_cost, 120581699.4153729, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 283387972.8581399, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 14242922.86112985, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 296.834875396348, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 247.78343626334382, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.56706595101369, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    
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
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_wheatstraw()
    test_LAOs()
    test_animal_bedding()
    test_lactic()
    test_ethanol_adipic()





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
    'generate_all_code',
    'generate_code',
    'print_results',
)

# Block speed-up for consistent testing. Note that speed-up wrapps functions,
# thus preventing code coverage to be analyzed
bst.speed_up = tmo.speed_up = flx.speed_up = lambda: None 

feedstocks_by_module = {
    'corn': 'corn',
    'lipidcane': 'lipidcane',
    'cornstover': 'cornstover',
    'sugarcane': 'sugarcane',
    'LAOs': 'glucose',
    'lactic': 'feedstock',
    'ethanol_adipic': 'feedstock',
    'HP': 'feedstock',
    'wheatstraw': 'wheatstraw',
    'animal_bedding': 'bedding',
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
}

marked_slow = {'wheatstraw', 'animal_bedding'}

configurations = {
    'HP': ('cellulosic', 'sugarcane')
}

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
    configuration_tag = f"_{configuration}" if configuration else ''
    configuration_name = f"'{configuration}'" if configuration else ''
    print(
    ("@pytest.mark.slow\n" if module_name in marked_slow else "") +
    f"""def test_{module_name}{configuration_tag}():
    bst.process_tools.default()
    from biorefineries import {module_name} as module
    try: module.load({configuration_name})
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
    bst.process_tools.default()
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
    
def test_corn():
    bst.process_tools.default()
    from biorefineries import corn as module
    try: module.load()
    except: pass
    feedstock = module.corn
    product = module.ethanol
    tea = module.corn_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.06385040370556491, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=1e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=1e-2)
    assert np.allclose(tea.sales, 73911903.8787555, rtol=1e-2)
    assert np.allclose(tea.material_cost, 54229612.27968329, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 67765794.19850408, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 6646814.904066811, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 93.75005124166793, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 96.87211085963307, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.338927565823793, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    bst.process_tools.default()
    
def test_lipidcane():
    bst.process_tools.default()
    from biorefineries import lipidcane as module
    try: module.load()
    except: pass
    feedstock = module.lipidcane
    product = module.ethanol
    tea = module.lipidcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.18049812456892045, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 102493171.42530274, rtol=1e-2)
    assert np.allclose(tea.material_cost, 61819200.176152, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 202305311.36327308, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -21065582.110587113, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 276.94917232585243, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 245.98384752627595, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.581267484995006, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 77.09915886508198, rtol=1e-2)
    bst.process_tools.default()
    
def test_cornstover():
    bst.process_tools.default()
    from biorefineries import cornstover as module
    try: module.load()
    except: pass
    feedstock = module.cornstover
    product = module.ethanol
    tea = module.cornstover_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.05158816935126135, rtol=1e-2)
    assert np.allclose(product.price, 0.7130236347531335, rtol=1e-2)
    assert np.allclose(tea.sales, 134171543.2799207, rtol=1e-2)
    assert np.allclose(tea.material_cost, 80366427.42505018, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 218206214.35534477, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -8369523.761661947, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 341.61907550660806, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 384.36304359465396, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.854606162607357, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 41.25381010614813, rtol=1e-2)
    bst.process_tools.default()
    
def test_sugarcane():
    bst.process_tools.default()
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.sugarcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.10748625239125097, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 87955077.0538273, rtol=1e-2)
    assert np.allclose(tea.material_cost, 60192842.50809233, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 162936659.34374583, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -6446162.428681537, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 312.98253582778057, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 250.2056638160482, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.630688343934354, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 30.2914653589393, rtol=1e-2)
    bst.process_tools.default()
    
def test_LAOs():
    bst.process_tools.default()
    from biorefineries import LAOs as module
    try: module.load()
    except: pass
    feedstock = module.glucose
    product = module.octene
    tea = module.LAOs_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.265, rtol=1e-2)
    assert np.allclose(product.price, 1.3040294044240837, rtol=1e-2)
    assert np.allclose(tea.sales, 163011798.28645134, rtol=1e-2)
    assert np.allclose(tea.material_cost, 135661590.95377028, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 61120398.09733076, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 2790455.914728851, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 40.0359632161767, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 125.33344794500302, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.596749741993358, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 3.5967497419933556, rtol=1e-2)
    bst.process_tools.default()

def test_lactic():
    bst.process_tools.default()
    from biorefineries import lactic as module
    try: module.load()
    except: pass
    feedstock = module.feedstock
    product = module.lactic_acid
    tea = module.lactic_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.06287583717512708, rtol=1e-2)
    assert np.allclose(product.price, 1.527148003894687, rtol=1e-2)
    assert np.allclose(tea.sales, 331933301.22110283, rtol=1e-2)
    assert np.allclose(tea.material_cost, 219074369.5100679, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 319884092.3488297, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 25053751.788215984, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 1827.5889438443508, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 1895.6439958372046, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 45.3971004352685, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    bst.process_tools.default()
    
def test_ethanol_adipic():
    bst.process_tools.default()
    from biorefineries import ethanol_adipic as module
    try: module.load()
    except: pass
    feedstock = module.feedstock
    product = module.ethanol
    tea = module.ethanol_adipic_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.06287581915485817, rtol=1e-2)
    assert np.allclose(product.price, 0.8352254290065602, rtol=1e-2)
    assert np.allclose(tea.sales, 219751260.20215675, rtol=1e-2)
    assert np.allclose(tea.material_cost, 120581699.41536081, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 282909147.08106965, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 14242922.861128656, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 296.8348753960908, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 247.7834362631533, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 26.567065951011458, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    bst.process_tools.default()
    
def test_HP_cellulosic():
    bst.process_tools.default()
    from biorefineries import HP as module
    try: module.load('cellulosic')
    except: pass
    feedstock = module.feedstock
    product = module.AcrylicAcid
    tea = module.HP_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.0628405632131775, rtol=1e-2)
    assert np.allclose(product.price, 1.2840111058021453, rtol=1e-2)
    assert np.allclose(tea.sales, 277692011.45432574, rtol=1e-2)
    assert np.allclose(tea.material_cost, 157608833.90096268, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 370125335.8581012, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 13917753.041283982, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 518.2382739594598, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 673.8593732993994, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 25.218803075458396, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    bst.process_tools.default()
    
def test_HP_sugarcane():
    bst.process_tools.default()
    from biorefineries import HP as module
    try: module.load('sugarcane')
    except: pass
    feedstock = module.feedstock
    product = module.AcrylicAcid
    tea = module.HP_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 1.0605606284638849, rtol=1e-2)
    assert np.allclose(tea.sales, 247364134.69047028, rtol=1e-2)
    assert np.allclose(tea.material_cost, 136306398.34921095, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 339171532.9914042, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 9618608.817279955, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 677.6593151558566, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 598.6913720328237, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 17.521718931637015, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.09291408768551465, rtol=1e-2)
    bst.process_tools.default()
    
@pytest.mark.slow
def test_wheatstraw():
    bst.process_tools.default()
    from biorefineries import wheatstraw as module
    try: module.load()
    except: pass
    feedstock = module.wheatstraw
    product = module.ethanol
    tea = module.wheatstraw_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.056999999999999995, rtol=1e-2)
    assert np.allclose(product.price, 0.9129911493074223, rtol=1e-2)
    assert np.allclose(tea.sales, 128434121.84797268, rtol=1e-2)
    assert np.allclose(tea.material_cost, 63431712.368926525, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 236647372.93625015, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -5101942.525598399, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 216.35913004517622, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 265.61952180127986, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.992736050523217, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 33.59904305663211, rtol=1e-2)
    bst.process_tools.default()
    
@pytest.mark.slow
def test_animal_bedding():
    bst.process_tools.default()
    from biorefineries import animal_bedding as module
    try: module.load()
    except: pass
    feedstock = module.bedding
    product = module.ethanol
    tea = module.bedding_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.003, rtol=1e-2)
    assert np.allclose(product.price, 0.7993436451648969, rtol=1e-2)
    assert np.allclose(tea.sales, 113884982.64989139, rtol=1e-2)
    assert np.allclose(tea.material_cost, 33147675.635541618, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 278658315.9940021, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -1000318.3282334469, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 166.72574074820272, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 207.67028281325636, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 32.11332834719463, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 34.19286633490201, rtol=1e-2)
    bst.process_tools.default()
    
    

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
    test_corn()
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_HP_cellulosic()
    test_HP_sugarcane()
    test_LAOs()
    test_lactic()
    test_ethanol_adipic()





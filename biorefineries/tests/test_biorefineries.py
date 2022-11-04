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
}

marked_slow = {'wheatstraw', 'animal_bedding'}

configurations = {
    'HP': ('cellulosic', 'sugarcane')
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
    assert np.allclose(tea.IRR, 0.054805285497855516, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723599.4175357, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525177.040248714, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 62117503.15997671, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7480503.478283916, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.08307189333905, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.69704380331558, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9369502188176067, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=5e-2)
    
@default_settings
def test_lipidcane():
    from biorefineries import lipidcane as module
    try: module.load()
    except: pass
    feedstock = module.lipidcane
    product = module.ethanol
    tea = module.lipidcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.2073843960349924, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 103632141.11986372, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58858162.489711784, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207001382.66925383, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -27751718.421849728, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 228.09384951510418, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 217.59755805859746, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.383142834194489, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 100.75700328840807, rtol=5e-2)

@default_settings
def test_cornstover():
    from biorefineries import cornstover as module
    try: module.load()
    except: pass
    feedstock = module.cornstover
    product = module.ethanol
    tea = module.cornstover_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.05158816935126135, rtol=5e-2)
    assert np.allclose(product.price, 0.7186468160404061, rtol=5e-2)
    assert np.allclose(tea.sales, 132743931.48928116, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82331609.7091591, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 207085777.19071975, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -9970210.2916593, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 313.8006121492292, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 383.8838046404844, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 20.75825813167742, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 44.681720569971326, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.sugarcane_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.14142881268514967, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88301969.9470638, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57278006.49173075, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 168767839.36499676, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -12161514.230486756, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 240.5558533310666, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 227.80222188805848, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.728459445437558, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 51.122502226762606, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.28830046834297, rtol=5e-2)
    assert np.allclose(tea.sales, 160973110.40423048, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.4604146, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 53669195.26697574, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2604134.508858938, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 38.31082754293686, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 123.69592726380715, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.5278714644386904, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 3.527871464438691, rtol=5e-2)
    
def test_lactic():
    bst.process_tools.default()
    from biorefineries import lactic as module
    module.load_system('SSCF')
    feedstock = module.feedstock
    product = module.lactic_acid
    tea = module.lactic_tea
    units = UnitGroup('Biorefinery', tea.units)
    module.lactic_sys.simulate()
    product.price = product_price = tea.solve_price(product)
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.06287583717512708, rtol=5e-2)
    assert np.allclose(product.price, 1.365160205800196, rtol=5e-2)
    assert np.allclose(tea.sales, 297352245.8641948, rtol=5e-2)
    assert np.allclose(tea.material_cost, 202146016.56312773, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 282073346.9023893, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 20669375.58387342, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 1842.5273042407728, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 1786.5220868424758, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 37.452662868510195, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=5e-2)
    bst.process_tools.default()

# Work is ongoing, at this stage, it is fine as long as these modules can load
# and systems can be simulated
def test_ethanol_adipic():
    bst.process_tools.default()
    from biorefineries import ethanol_adipic as module
    module.load_system('acid', 'HMPP')
    module.biorefinery.simulate()
    module.load_system('AFEX', 'CPP_AFEX')
    module.biorefinery.simulate()
    module.load_system('base', 'HMPP')
    module.biorefinery.simulate()
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
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.056999999999999995, rtol=5e-2)
    assert np.allclose(product.price, 0.9033321094343063, rtol=5e-2)
    assert np.allclose(tea.sales, 128751220.49360518, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63431132.07841703, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 237441843.41517347, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5034648.66311628, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 216.35813625516792, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 265.7465556328792, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.125593758589847, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 33.592005154050234, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.1, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.003, rtol=5e-2)
    assert np.allclose(product.price, 0.7890118048567955, rtol=5e-2)
    assert np.allclose(tea.sales, 107275778.85699502, rtol=5e-2)
    assert np.allclose(tea.material_cost, 32460366.631773323, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 262912439.97063047, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -2560757.2115471847, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 165.1624907961539, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 211.85054691475685, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 29.19232139864, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 34.515818677863315, rtol=5e-2)
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
    generate_all_code()
    # test_corn()
    # test_sugarcane()
    # test_lipidcane()
    # test_cornstover()
    # test_LAOs()
    # test_HP_cellulosic()
    # test_HP_sugarcane()
    # test_lactic()
    # test_ethanol_adipic()





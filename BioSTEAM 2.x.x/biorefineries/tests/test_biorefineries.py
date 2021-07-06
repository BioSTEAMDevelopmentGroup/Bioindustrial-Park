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
    

    bst.process_tools.default()
    from biorefineries import corn as module
    try: module.load()
    except: pass
    feedstock = module.corn
    product = module.ethanol
    tea = module.corn_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.037843194745526644, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=1e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=1e-2)
    assert np.allclose(tea.sales, 73911402.04639979, rtol=1e-2)
    assert np.allclose(tea.material_cost, 54229609.23083669, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 67931892.15209496, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 8267419.421774661, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 135.60034920381048, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 141.00001509562642, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.340125651483439, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
    bst.process_tools.default()
    
def test_corn():
    bst.process_tools.default()
    from biorefineries import corn as module
    try: module.load()
    except: pass
    feedstock = module.corn
    product = module.ethanol
    tea = module.corn_tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.03989508733609, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=1e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=1e-2)
    assert np.allclose(tea.sales, 73909879.42964743, rtol=1e-2)
    assert np.allclose(tea.material_cost, 54229305.139442496, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 67980658.30941775, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 8136250.075012307, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 135.8019116579886, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 139.05181894956826, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 2.33924346153039, rtol=1e-2)
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
    assert np.allclose(tea.IRR, 0.21047391367462487, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 108715833.87257311, rtol=1e-2)
    assert np.allclose(tea.material_cost, 59023796.119482465, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 239038533.43349013, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -34582229.25333484, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 245.92023657075794, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 280.75554419221567, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 16.42156294204761, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 126.51505024946395, rtol=1e-2)
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
    assert np.allclose(product.price, 0.7141047799147656, rtol=1e-2)
    assert np.allclose(tea.sales, 131519594.87222457, rtol=1e-2)
    assert np.allclose(tea.material_cost, 82696666.93499637, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 213729761.48123008, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -12216823.378096422, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 308.3841985204103, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 298.4404395355398, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.456532820335028, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 46.85379645812982, rtol=1e-2)
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
    assert np.allclose(tea.IRR, 0.12659387829384772, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=1e-2)
    assert np.allclose(product.price, 0.789, rtol=1e-2)
    assert np.allclose(tea.sales, 88129976.99654803, rtol=1e-2)
    assert np.allclose(tea.material_cost, 59842659.397153035, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 205172809.39941, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -19201338.64548274, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 284.23468528204995, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 304.1271932338852, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 15.82128374128548, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 76.76409793701727, rtol=1e-2)
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
    assert np.allclose(product.price, 1.3269049845098067, rtol=1e-2)
    assert np.allclose(tea.sales, 161335906.29907194, rtol=1e-2)
    assert np.allclose(tea.material_cost, 135661164.88724267, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 57347781.65516735, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 1554844.4698447809, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 39.587331392856626, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 127.39469106251019, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.5306112844050603, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 7.271296144805024, rtol=1e-2)
    bst.process_tools.default()
    
def test_lactic():
    bst.process_tools.default()
    from biorefineries import lactic as module
    module.load_system('SSCF')
    feedstock = module.feedstock
    product = module.lactic_acid
    tea = module.lactic_tea
    units = UnitGroup('Biorefinery', tea.units)
    module.lactic_sys.simulate()
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.06287583717512708, rtol=1e-2)
    assert np.allclose(tea.solve_price(product), 1.237548148885074, rtol=1e-2)
    assert np.allclose(tea.sales, 54143215.78374311, rtol=1e-2)
    assert np.allclose(tea.material_cost, 219104068.95373908, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 320073822.05670184, rtol=1e-2)
    assert np.allclose(tea.utility_cost, 22275814.89016047, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 1827.5926437830065, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 1895.6520852552394, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 40.363511796333356, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 0.0, rtol=1e-2)
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
    assert np.allclose(tea.IRR, 0.1, rtol=1e-2)
    assert np.allclose(feedstock.price, 0.056999999999999995, rtol=1e-2)
    assert np.allclose(product.price, 0.9033321094343063, rtol=1e-2)
    assert np.allclose(tea.sales, 128751220.49360518, rtol=1e-2)
    assert np.allclose(tea.material_cost, 63431132.07841703, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 237441843.41517347, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -5034648.66311628, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 216.35813625516792, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 265.7465556328792, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.125593758589847, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 33.592005154050234, rtol=1e-2)
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
    assert np.allclose(product.price, 0.7890118048567955, rtol=1e-2)
    assert np.allclose(tea.sales, 107275778.85699502, rtol=1e-2)
    assert np.allclose(tea.material_cost, 32460366.631773323, rtol=1e-2)
    assert np.allclose(tea.installed_equipment_cost, 262912439.97063047, rtol=1e-2)
    assert np.allclose(tea.utility_cost, -2560757.2115471847, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 165.1624907961539, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 211.85054691475685, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 29.19232139864, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 34.515818677863315, rtol=1e-2)
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
    
# if __name__ == '__main__':
    # test_corn()
    # test_sugarcane()
    # test_lipidcane()
    # test_cornstover()
    # # test_HP_cellulosic()
    # # test_HP_sugarcane()
    # test_LAOs()
    # test_lactic()
    # test_ethanol_adipic()





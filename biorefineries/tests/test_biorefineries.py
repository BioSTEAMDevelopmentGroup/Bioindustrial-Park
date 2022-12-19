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
    assert np.allclose(tea.IRR, 0.054793135761420995, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.13227735731092652, rtol=5e-2)
    assert np.allclose(product.price, 0.48547915353569393, rtol=5e-2)
    assert np.allclose(tea.sales, 74723596.88174264, rtol=5e-2)
    assert np.allclose(tea.material_cost, 55525176.94552828, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 62093937.04119991, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 7484262.234701366, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 95.0630956270271, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 114.83385319881219, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 1.9369534237489578, rtol=5e-2)
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
    assert np.allclose(tea.IRR, 0.20854069159126307, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 102694047.72636937, rtol=5e-2)
    assert np.allclose(tea.material_cost, 58856880.18876951, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 206545752.4568234, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -28931835.470725324, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 207.10490688845118, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 218.17006801938587, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.258278600210271, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 104.41456236827366, rtol=5e-2)
    
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
    assert np.allclose(product.price, 0.7333719888983871, rtol=5e-2)
    assert np.allclose(tea.sales, 135536068.36974084, rtol=5e-2)
    assert np.allclose(tea.material_cost, 82209367.24874295, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 200715685.00584513, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -5526978.776637238, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 363.3753814500929, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 267.1471947922442, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 18.578785100558147, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 33.26470926469097, rtol=5e-2)
    
@default_settings
def test_sugarcane():
    from biorefineries import sugarcane as module
    try: module.load()
    except: pass
    feedstock = module.sugarcane
    product = module.ethanol
    tea = module.tea
    units = UnitGroup('Biorefinery', tea.units)
    assert np.allclose(tea.IRR, 0.1428763497389623, rtol=5e-2)
    assert np.allclose(feedstock.price, 0.03455, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 88302442.78606105, rtol=5e-2)
    assert np.allclose(tea.material_cost, 57277062.32223333, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 167490968.14053738, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -12194793.465249699, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.720557038228, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 223.38600265957055, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.494504639075067, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 50.99521162967652, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.038715056359816404, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 81211491.73735546, rtol=5e-2)
    assert np.allclose(tea.material_cost, 63948614.355520636, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 120388026.76088502, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -16597694.443256002, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 246.6805648809778, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 245.42152922165081, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 7.978780546507185, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 66.40030690875176, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03132217680573546, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 121030074.16768698, rtol=5e-2)
    assert np.allclose(tea.material_cost, 59771930.64433406, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 220078597.56212905, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1402757.249752447, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 335.9393189712925, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 309.64371458799303, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 21.60763866514533, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 21.607638665145313, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.042396337907574855, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 71993360.47700194, rtol=5e-2)
    assert np.allclose(tea.material_cost, 71394555.13860849, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 163350004.57130826, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -42831260.10448998, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 252.19630057888364, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 223.6801185520896, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.963759583476811, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 152.64785261438897, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04631588328993245, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 163071076.93074393, rtol=5e-2)
    assert np.allclose(tea.material_cost, 92970839.40374137, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 250788391.9581871, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1391860.2888903257, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 390.459823210292, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 338.59112999347474, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.386696911774806, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.386696911774795, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03697502362337886, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 59505464.7419263, rtol=5e-2)
    assert np.allclose(tea.material_cost, 60272559.94946644, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 155595541.53343153, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -42131213.804979876, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 247.03602244311102, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 216.84363200874705, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.732652417909573, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 150.14149487418058, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04017535907850782, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 146687116.13255346, rtol=5e-2)
    assert np.allclose(tea.material_cost, 79161090.83887991, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 241177664.53088266, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1397939.8456321142, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 377.19554892676285, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 325.14475343375193, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.924499079714153, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 22.92449907971412, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.042544735854539106, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 103451891.41629674, rtol=5e-2)
    assert np.allclose(tea.material_cost, 91052136.37564363, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 126823497.00086042, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -23025179.19294905, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 116.80745940094549, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 194.79140976769625, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.22548156833173, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 81.12965160997803, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.03850019156182423, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 158349820.86452833, rtol=5e-2)
    assert np.allclose(tea.material_cost, 93403868.68606952, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 230951296.42368263, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1818646.9865419918, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 151.82709847956076, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 302.30563988986455, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.954406254392953, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 23.063573067248704, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.04671756133087959, rtol=5e-2)
    assert np.allclose(product.price, 0.7256416231373818, rtol=5e-2)
    assert np.allclose(tea.sales, 92435850.58665347, rtol=5e-2)
    assert np.allclose(tea.material_cost, 101606324.27962075, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 164707855.938459, rtol=5e-2)
    assert np.allclose(tea.utility_cost, -53124217.03632314, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 227.03389380488375, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 193.43642943115623, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 8.86305410205487, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 133.38171690786595, rtol=5e-2)
    
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
    assert np.allclose(feedstock.price, 0.05379967252808456, rtol=5e-2)
    assert np.allclose(product.price, 0.789, rtol=5e-2)
    assert np.allclose(tea.sales, 207226491.06851304, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135996457.0549305, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 254562018.17076746, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 1312895.6652508373, rtol=5e-2)
    assert np.allclose(units.get_heating_duty(), 311.20847327831564, rtol=5e-2)
    assert np.allclose(units.get_cooling_duty(), 265.88313200248007, rtol=5e-2)
    assert np.allclose(units.get_electricity_consumption(), 23.229887222503045, rtol=5e-2)
    assert np.allclose(units.get_electricity_production(), 28.321981978621604, rtol=5e-2)
    
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
    assert np.allclose(product.price, 1.2918667227938558, rtol=5e-2)
    assert np.allclose(tea.sales, 161435883.1972842, rtol=5e-2)
    assert np.allclose(tea.material_cost, 135658829.45874006, rtol=5e-2)
    assert np.allclose(tea.installed_equipment_cost, 55049531.41964978, rtol=5e-2)
    assert np.allclose(tea.utility_cost, 2702388.184876362, rtol=5e-2)
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





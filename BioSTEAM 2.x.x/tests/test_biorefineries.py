# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Each function in this module tests a biorefinery.If all tests passed, no error 
is raised.

"""
import numpy as np
import biosteam as bst
from biosteam.process_tools import UnitGroup

__all__ = ('test_sugarcane', 
           'test_lipidcane', 
           'test_cornstover',
           'test_LAOs',
           'test_ethanol_adipic',
           'test_lactic',)

def test_sugarcane():
    from biorefineries import sugarcane as sc
    sc.load()
    units = UnitGroup('Biorefinery', sc.sugarcane_tea.units)
    assert np.allclose(sc.sugarcane_tea.IRR, 0.10326317303317518)
    assert np.allclose(sc.sugarcane_tea.sales, 87641108.77580193)
    assert np.allclose(sc.sugarcane_tea.material_cost, 60096414.780067384)
    assert np.allclose(sc.sugarcane_tea.installed_equipment_cost, 111490520.90754643)
    assert np.allclose(sc.sugarcane_tea.utility_cost, -8595640.622981431)
    assert np.allclose(units.get_heating_duty(), 272.488831440168)
    assert np.allclose(units.get_cooling_duty(), 343.762814732534)
    assert np.allclose(units.get_electricity_consumption(), 9.827812210952846)
    assert np.allclose(units.get_electricity_production(), 37.37794241281641)

def test_lipidcane():
    from biorefineries import lipidcane as lc
    lc.load()
    units = UnitGroup('Biorefinery', lc.lipidcane_tea.units)
    assert np.allclose(lc.lipidcane_tea.IRR, 0.17829951498180963)
    assert np.allclose(lc.lipidcane_tea.sales, 102410583.81934823)
    assert np.allclose(lc.lipidcane_tea.material_cost, 61762909.6376681)
    assert np.allclose(lc.lipidcane_tea.installed_equipment_cost, 139554763.73812002)
    assert np.allclose(lc.lipidcane_tea.utility_cost, -25856786.959638454)
    assert np.allclose(units.get_heating_duty(), 198.1915651157178)
    assert np.allclose(units.get_cooling_duty(), 304.5184189296353)
    assert np.allclose(units.get_electricity_consumption(), 9.645038160426235)
    assert np.allclose(units.get_electricity_production(), 92.5193553387546)

def test_cornstover():
    from biorefineries import cornstover as cs
    cs.load()
    MESP = cs.cornstover_tea.solve_price(cs.ethanol)
    units = UnitGroup('Biorefinery', cs.cornstover_tea.units)
    assert np.allclose(MESP, 0.7382850997534475)
    assert np.allclose(cs.cornstover_tea.sales, 134041032.28829786)
    assert np.allclose(cs.cornstover_tea.material_cost, 82900270.63843909)
    assert np.allclose(cs.cornstover_tea.installed_equipment_cost, 219990786.55181444)
    assert np.allclose(cs.cornstover_tea.utility_cost, -11054921.042512089)
    assert np.allclose(units.get_heating_duty(), 318.0490631730075)
    assert np.allclose(units.get_cooling_duty(), 364.9938841327776)
    assert np.allclose(units.get_electricity_consumption(), 22.36637203703606)
    assert np.allclose(units.get_electricity_production(), 45.3481845362712)
 
def test_LAOs():
    from biorefineries import LAOs as laos
    laos.load()
    MPSP = laos.get_LAOs_MPSP()
    units = UnitGroup('Biorefinery', laos.LAOs_tea.units)
    assert np.allclose(MPSP, 1388.5255959156789, rtol=0.01)
    assert np.allclose(laos.LAOs_tea.sales, 194084400.8831837, rtol=0.01)
    assert np.allclose(laos.LAOs_tea.material_cost, 163075600.3024821, rtol=0.01)
    assert np.allclose(laos.LAOs_tea.installed_equipment_cost, 70477488.66232976, rtol=0.03)
    assert np.allclose(laos.LAOs_tea.utility_cost, 3760758.8485403117, rtol=0.01)
    assert np.allclose(units.get_heating_duty(), 60.46838019012832, rtol=0.01)
    assert np.allclose(units.get_cooling_duty(), 136.49148628223284, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 3.155372464756524, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 3.155372464756526, rtol=0.01)   
    
def test_ethanol_adipic():
    bst.process_tools.default_utilities()
    from biorefineries import ethanol_adipic
    acid = ethanol_adipic.system_acid
    MESP = acid.simulate_get_MESP()
    assert np.allclose(MESP, 2.3383974799850873, atol=0.001)
    tea = acid.ethanol_tea
    assert np.allclose(tea.sales, 142814170.6871417, rtol=0.001)
    assert np.allclose(tea.material_cost, 98030872.485875, rtol=0.001)
    assert np.allclose(tea.installed_equipment_cost, 201617183.59793562, rtol=0.001)
    assert np.allclose(tea.utility_cost, -12104425.737786738, rtol=0.001)
    units = UnitGroup('Biorefinery', acid.ethanol_tea.units)
    assert np.allclose(acid.CHP.system_heating_demand/1e6, 359.9811664916131, rtol=0.001)
    assert np.allclose(units.get_cooling_duty(), 330.92129172498926, rtol=0.001)
    assert np.allclose(units.get_electricity_consumption(), 23.900435032178247, rtol=0.001)
    assert np.allclose(units.get_electricity_production(), 45.06746566975969, rtol=0.001)
    
    base = ethanol_adipic.system_base
    MESP = base.simulate_get_MESP()
    assert np.allclose(MESP, 2.206842775418839, atol=0.001)
    tea = base.ethanol_adipic_tea
    assert np.allclose(tea.sales, 235176081.26473045, rtol=0.001)
    assert np.allclose(tea.material_cost, 131288487.04028216, rtol=0.001)
    assert np.allclose(tea.installed_equipment_cost, 298528854.16088444, rtol=0.001)
    assert np.allclose(tea.utility_cost, 14579802.123597547, rtol=0.001)
    units = UnitGroup('Biorefinery', base.ethanol_adipic_tea.units)
    assert np.allclose(base.CHP.system_heating_demand/1e6, 336.5735818785919, rtol=0.001)
    assert np.allclose(units.get_cooling_duty(), 359.1181296316559, rtol=0.001)
    assert np.allclose(units.get_electricity_consumption(), 27.19544073551338, rtol=0.001)
    assert np.allclose(units.get_electricity_production(), 0.0)
    
def test_lactic():
    bst.process_tools.default_utilities()
    from biorefineries import lactic
    system = lactic.system
    MPSP = system.simulate_get_MPSP()
    assert np.allclose(MPSP, 1.47295738396085, atol=0.001)
    tea = system.lactic_tea
    assert np.allclose(tea.sales, 349336905.47623944, rtol=0.001)
    assert np.allclose(tea.material_cost, 231659140.1072927, rtol=0.001)
    assert np.allclose(tea.installed_equipment_cost, 337250395.9673476, rtol=0.001)
    assert np.allclose(tea.utility_cost, 25533406.595646832, rtol=0.001)
    units = UnitGroup('Biorefinery', system.lactic_tea.units)
    assert np.allclose(system.CHP.system_heating_demand/1e6, 1982.4854859150275, rtol=0.001)
    assert np.allclose(-system.CT.system_cooling_water_duty/1e6, 1850.3994601177017, rtol=0.001)
    assert np.allclose(units.get_electricity_consumption(), 46.26622924484822, rtol=0.001)
    assert np.allclose(units.get_electricity_production(), 0.0)

    
if __name__ == '__main__':
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_LAOs()
    test_lactic()
















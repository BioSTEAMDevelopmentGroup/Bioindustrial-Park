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
from biosteam.process_tools import UnitGroup

__all__ = ('test_sugarcane', 
           'test_lipidcane', 
           'test_cornstover',
           'test_LAOs')

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
    assert np.allclose(MPSP, 1388.5255959156789)
    assert np.allclose(laos.LAOs_tea.sales, 194084400.8831837)
    assert np.allclose(laos.LAOs_tea.material_cost, 163075600.3024821)
    assert np.allclose(laos.LAOs_tea.installed_equipment_cost, 70477488.66232976)
    assert np.allclose(laos.LAOs_tea.utility_cost, 3760758.8485403117)
    assert np.allclose(units.get_heating_duty(), 60.46838019012832)
    assert np.allclose(units.get_cooling_duty(), 136.49148628223284)
    assert np.allclose(units.get_electricity_consumption(), 3.155372464756524)
    assert np.allclose(units.get_electricity_production(), 3.155372464756526)   
    
if __name__ == '__main__':
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_LAOs()
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
    assert np.allclose(sc.sugarcane_tea.IRR, 0.09952201383654127)
    assert np.allclose(sc.sugarcane_tea.sales, 87641104.25548652)
    assert np.allclose(sc.sugarcane_tea.material_cost, 60096508.050463505)
    assert np.allclose(sc.sugarcane_tea.installed_equipment_cost, 117725174.12298785)
    assert np.allclose(sc.sugarcane_tea.utility_cost, -8589365.387334356)
    assert np.allclose(units.get_heating_duty(), 272.487235278769)
    assert np.allclose(units.get_cooling_duty(), 344.88392082847355)
    assert np.allclose(units.get_electricity_consumption(), 9.853875285162502)
    assert np.allclose(units.get_electricity_production(), 37.38389255225982)

def test_lipidcane():
    from biorefineries import lipidcane as lc
    lc.load()
    units = UnitGroup('Biorefinery', lc.lipidcane_tea.units)
    assert np.allclose(lc.lipidcane_tea.IRR, 0.1733032594846777)
    assert np.allclose(lc.lipidcane_tea.sales, 102410581.55979738)
    assert np.allclose(lc.lipidcane_tea.material_cost, 61762971.048900105)
    assert np.allclose(lc.lipidcane_tea.installed_equipment_cost, 148647752.99924695)
    assert np.allclose(lc.lipidcane_tea.utility_cost, -25886609.91411801)
    assert np.allclose(units.get_heating_duty(), 197.7426409332782)
    assert np.allclose(units.get_cooling_duty(), 305.16703491311733)
    assert np.allclose(units.get_electricity_consumption(), 9.65815356491132)
    assert np.allclose(units.get_electricity_production(), 92.62805713580244)

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
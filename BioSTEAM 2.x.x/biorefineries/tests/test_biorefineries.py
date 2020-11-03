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
from biosteam.process_tools import UnitGroup
import pytest

def test_sugarcane():
    from biorefineries import sugarcane as sc
    sc.load()
    units = UnitGroup('Biorefinery', sc.sugarcane_tea.units)
    assert np.allclose(sc.sugarcane_tea.IRR, 0.10276414072162979, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.sales, 87558012.26683149, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.material_cost, 60096508.050463505, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.installed_equipment_cost, 178821241.29052082, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.utility_cost, -8691762.69137734, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 271.1415397273802, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 341.9636663890628, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.843412469700933, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 37.70162622411549, rtol=1e-2)

def test_lipidcane():
    from biorefineries import lipidcane as lc
    lc.load()
    units = UnitGroup('Biorefinery', lc.lipidcane_tea.units)
    assert np.allclose(lc.lipidcane_tea.IRR, 0.17804673529883677, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.sales, 102369101.03579773, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.material_cost, 61761912.9045123, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.installed_equipment_cost, 222127644.25250342, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.utility_cost, -25937821.45942718, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 197.06976435736067, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 303.70937200373174, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.65288763259321, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 92.78693077178295, rtol=1e-2)

def test_cornstover():
    from biorefineries import cornstover as cs
    cs.load()
    MESP = cs.cornstover_tea.solve_price(cs.ethanol)
    units = UnitGroup('Biorefinery', cs.cornstover_tea.units)
    assert np.allclose(MESP, 0.7383610895500932, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.sales, 134055638.87837549, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.material_cost, 82902092.80482696, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.installed_equipment_cost, 220012587.48251432, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.utility_cost, -11047774.698866017, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 318.05635812416523, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 365.67806026618115, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.371322764496814, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 45.33827889984683, rtol=1e-2)
    
def test_LAOs():
    from biorefineries import LAOs as laos
    laos.load()
    MPSP = laos.get_LAOs_MPSP()
    units = UnitGroup('Biorefinery', laos.LAOs_tea.units)
    assert np.allclose(MPSP, 1226.0016718824597, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.sales, 168789574.412942, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.material_cost, 135661583.58974993, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.installed_equipment_cost, 76661769.51734665, rtol=0.05)
    assert np.allclose(laos.LAOs_tea.utility_cost, 3811246.803774001, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 59.19501587183199, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 134.66774238880174, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.689496361470118, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 3.6894963614701215, rtol=1e-2)   

@pytest.mark.slow
def test_wheatstraw():
    from biorefineries import wheatstraw as ws
    ws.load()
    MESP = ws.wheatstraw_tea.solve_price(ws.ethanol)
    units = UnitGroup('Biorefinery', ws.wheatstraw_tea.units)
    assert np.allclose(MESP, 0.9384667423768497, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.sales, 130275585.21251984, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.material_cost, 63298497.065546915, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.installed_equipment_cost, 242292530.10266995, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.utility_cost, -4373067.719624106, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 219.0919211247056, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 265.2733205865057, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.88974837441627, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 31.98081487705149, rtol=1e-2)

@pytest.mark.slow
def test_annimal_bedding():
    from biorefineries import animal_bedding as ab
    ab.load()
    MESP = ab.bedding_tea.solve_price(ab.ethanol)
    units = UnitGroup('Biorefinery', ab.bedding_tea.units)
    assert np.allclose(MESP, 0.8188585523640879, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.sales, 115097527.36582053, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.material_cost, 33081968.8237324, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.installed_equipment_cost, 285002340.76542723, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.utility_cost, -1178035.514879562, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 161.56120885595166, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 201.26648298810298, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 31.999853945401934, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 34.44884396679516, rtol=1e-2)
    
@pytest.mark.slow
def test_lactic():
    bst.process_tools.default_utilities()
    from biorefineries import lactic
    system = lactic.system
    MPSP = system.lactic_acid.price
    assert np.allclose(MPSP, 1.5704836997441853, atol=0.01)
    tea = system.lactic_tea
    assert np.allclose(tea.sales, 340805307.1489925, rtol=0.01)
    assert np.allclose(tea.material_cost, 225035229.6439402, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 329160559.8556272, rtol=0.01)
    assert np.allclose(tea.utility_cost, 25551416.34519198, rtol=0.01)
    units = UnitGroup('Biorefinery', system.lactic_tea.units)
    assert np.allclose(system.CHP.system_heating_demand/1e6, 1881.5136539040284, rtol=0.01)
    assert np.allclose(-system.CT.system_cooling_water_duty/1e6, 1714.9228013906093, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 46.29886269694857, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 0.0)
    
if __name__ == '__main__':
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_wheatstraw()
    test_LAOs()
    test_annimal_bedding()
    test_lactic()





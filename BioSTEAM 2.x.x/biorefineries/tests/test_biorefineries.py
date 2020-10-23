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
    assert np.allclose(sc.sugarcane_tea.IRR, 0.10319729375184486, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.sales, 87641104.25548652, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.material_cost, 60096508.050463505, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.installed_equipment_cost, 111512561.90122335, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.utility_cost, -8589365.38733437, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 272.4872352787688, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 344.88392082847315, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.853875285162502, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 37.383892552259844, rtol=1e-2)

def test_lipidcane():
    from biorefineries import lipidcane as lc
    lc.load()
    units = UnitGroup('Biorefinery', lc.lipidcane_tea.units)
    assert np.allclose(lc.lipidcane_tea.IRR, 0.1783634347376759, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.sales, 102410581.5597974, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.material_cost, 61762971.048900105, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.installed_equipment_cost, 139559824.39299855, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.utility_cost, -25886609.91411801, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 197.74264093327824, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 305.1670349131174, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.65815356491132, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 92.62805713580244, rtol=1e-2)

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
    
def test_ethanol_adipic():
    bst.process_tools.default_utilities()
    from biorefineries import ethanol_adipic
    acid = ethanol_adipic.system_acid
    MESP = acid.ethanol.price
    assert np.allclose(MESP, 2.512942530141015, atol=0.01)
    tea = acid.ethanol_tea
    assert np.allclose(tea.sales, 153621809.87439528, rtol=0.01)
    assert np.allclose(tea.material_cost, 112976158.84156583, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 202999332.38051096, rtol=0.01)
    assert np.allclose(tea.utility_cost, -16654751.636141604, rtol=0.01)
    units = UnitGroup('Biorefinery', acid.ethanol_tea.units)
    assert np.allclose(acid.CHP.system_heating_demand/1e6, 342.19255164427744, rtol=0.01)
    assert np.allclose(units.get_cooling_duty(), 335.3383652559594, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 23.921945234219077, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 53.04613879616649, rtol=0.01)
    
    base = ethanol_adipic.system_base
    MESP = base.ethanol.price
    assert np.allclose(MESP, 2.723886457801242, atol=0.01)
    tea = base.ethanol_adipic_tea
    assert np.allclose(tea.sales, 220968148.48076087, rtol=0.01)
    assert np.allclose(tea.material_cost, 121378733.07275604, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 284333240.3730375, rtol=0.01)
    assert np.allclose(tea.utility_cost, 14299236.824111573, rtol=0.01)
    units = UnitGroup('Biorefinery', base.ethanol_adipic_tea.units)
    assert np.allclose(base.CHP.system_heating_demand/1e6, 307.55991674861224, rtol=0.01)
    assert np.allclose(units.get_cooling_duty(), 259.37733506080116, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 26.67210736583321, rtol=0.001)
    assert np.allclose(units.get_electricity_production(), 0.0)
    
def test_lactic():
    bst.process_tools.default_utilities()
    from biorefineries import lactic
    system = lactic.system
    MPSP = system.lactic_acid.price
    assert np.allclose(MPSP, 1.548015650016793, atol=0.01)
    tea = system.lactic_tea
    assert np.allclose(tea.sales, 335881119.8206509, rtol=0.01)
    assert np.allclose(tea.material_cost, 221588718.1301599, rtol=0.01)
    assert np.allclose(tea.installed_equipment_cost, 324175792.74265414, rtol=0.01)
    assert np.allclose(tea.utility_cost, 25350111.57914885, rtol=0.01)
    units = UnitGroup('Biorefinery', system.lactic_tea.units)
    assert np.allclose(system.CHP.system_heating_demand/1e6, 1812.5943347629961, rtol=0.01)
    assert np.allclose(-system.CT.system_cooling_water_duty/1e6, 1715.9264327014253, rtol=0.01)
    assert np.allclose(units.get_electricity_consumption(), 45.934100853716096, rtol=0.01)
    assert np.allclose(units.get_electricity_production(), 0.0)

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
    
if __name__ == '__main__':
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_wheatstraw()
    test_LAOs()
    test_annimal_bedding()
#    test_ethanol_adipic()
#    test_lactic()





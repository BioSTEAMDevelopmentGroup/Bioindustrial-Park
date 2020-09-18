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
















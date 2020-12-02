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
from biosteam.process_tools import UnitGroup
import pytest

# Block speed-up for consistent testing. Note that speed-up wrapps functions,
# thus preventing code coverage to be analyzed
bst.speed_up = tmo.speed_up = flx.speed_up = lambda: None 

def print_results(tea):
    units = UnitGroup('Biorefinery', tea.units)
    print(tea.IRR)
    print(tea.sales)
    print(tea.material_cost)
    print(tea.installed_equipment_cost)
    print(tea.utility_cost)
    print(units.get_heating_duty())
    print(units.get_cooling_duty())
    print(units.get_electricity_consumption())
    print(units.get_electricity_production())
    
def test_sugarcane():
    from biorefineries import sugarcane as sc
    sc.load()
    units = UnitGroup('Biorefinery', sc.sugarcane_tea.units)
    assert np.allclose(sc.sugarcane_tea.IRR, 0.1064405821621119, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.sales, 87558012.26683149, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.material_cost, 60096508.050463505, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.installed_equipment_cost, 181877578.2850932, rtol=1e-2)
    assert np.allclose(sc.sugarcane_tea.utility_cost, -10028790.139157953, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 263.50013407609373, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 329.94966728221334, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.831303573342776, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 41.974861711669554, rtol=1e-2)

def test_lipidcane():
    from biorefineries import lipidcane as lc
    lc.load()
    units = UnitGroup('Biorefinery', lc.lipidcane_tea.units)
    assert np.allclose(lc.lipidcane_tea.IRR, 0.17980836639978406, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.sales, 102369101.03579773, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.material_cost, 61760752.25323551, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.installed_equipment_cost, 224875041.54906386, rtol=1e-2)
    assert np.allclose(lc.lipidcane_tea.utility_cost, -27329649.666423105, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 191.59071424606782, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 294.9121447943192, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 9.658208262905742, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 97.25323924503112, rtol=1e-2)

def test_cornstover():
    from biorefineries import cornstover as cs
    cs.load()
    MESP = cs.cornstover_tea.solve_price(cs.ethanol)
    units = UnitGroup('Biorefinery', cs.cornstover_tea.units)
    assert np.allclose(MESP, 0.7308776819169308, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.sales, 132680811.62702852, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.material_cost, 82899286.3495096, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.installed_equipment_cost, 221139367.31134215, rtol=1e-2)
    assert np.allclose(cs.cornstover_tea.utility_cost, -12676528.501955263, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 310.33962099605156, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 359.055441898979, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.34435621558365, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 48.69728990482738, rtol=1e-2)
    
def test_LAOs():
    from biorefineries import LAOs as laos
    laos.load()
    MPSP = laos.get_LAOs_MPSP()
    units = UnitGroup('Biorefinery', laos.LAOs_tea.units)
    assert np.allclose(MPSP, 1206.123978150559, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.sales, 165952296.1656475, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.material_cost, 135661583.5897499, rtol=1e-2)
    assert np.allclose(laos.LAOs_tea.installed_equipment_cost, 70837450.03577328, rtol=0.05)
    assert np.allclose(laos.LAOs_tea.utility_cost, 2795235.378962637, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 40.00054620871382, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 125.33654746457061, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 3.6274634582889065, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 3.627463458288911, rtol=1e-2)   

@pytest.mark.slow
def test_wheatstraw():
    from biorefineries import wheatstraw as ws
    ws.load()
    MESP = ws.wheatstraw_tea.solve_price(ws.ethanol)
    units = UnitGroup('Biorefinery', ws.wheatstraw_tea.units)
    assert np.allclose(MESP, 0.9298747302315753, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.sales, 129763770.04047716, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.material_cost, 63431955.06015967, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.installed_equipment_cost, 243695707.37332124, rtol=1e-2)
    assert np.allclose(ws.wheatstraw_tea.utility_cost, -5355626.786941081, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 214.26464305084303, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 261.681819819361, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 22.959889541896427, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 34.09357472698683, rtol=1e-2)

@pytest.mark.slow
def test_annimal_bedding():
    from biorefineries import animal_bedding as ab
    ab.load()
    MESP = ab.bedding_tea.solve_price(ab.ethanol)
    units = UnitGroup('Biorefinery', ab.bedding_tea.units)
    assert np.allclose(MESP, 0.8075091100223122, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.sales, 114420341.00709279, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.material_cost, 33147675.847261157, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.installed_equipment_cost, 285644893.62946343, rtol=1e-2)
    assert np.allclose(ab.bedding_tea.utility_cost, -2044647.6771812998, rtol=1e-2)
    assert np.allclose(units.get_heating_duty(), 157.88539960199708, rtol=1e-2)
    assert np.allclose(units.get_cooling_duty(), 198.23736376471848, rtol=1e-2)
    assert np.allclose(units.get_electricity_consumption(), 32.02959950750465, rtol=1e-2)
    assert np.allclose(units.get_electricity_production(), 36.280168947419845, rtol=1e-2)
    
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
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_wheatstraw()
    test_LAOs()
    test_annimal_bedding()
    test_lactic()
    test_ethanol_adipic()





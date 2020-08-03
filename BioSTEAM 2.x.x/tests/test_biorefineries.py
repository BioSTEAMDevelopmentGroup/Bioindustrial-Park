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
from biosteam.process_tools import UnitGroup

__all__ = ('test_lipidcane', 'test_cornstover')

def test_lipidcane():
    """
    Test lipidcane biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_lipidcane()
    
    """
    from biorefineries.lipidcane import lipidcane_tea
    units = UnitGroup('Biorefinery', lipidcane_tea.units)
    assert np.allclose(lipidcane_tea.IRR, 0.17829951498180963)
    assert np.allclose(lipidcane_tea.sales, 102410583.81934823)
    assert np.allclose(lipidcane_tea.material_cost, 61762909.6376681)
    assert np.allclose(lipidcane_tea.installed_equipment_cost, 139554763.73812002)
    assert np.allclose(lipidcane_tea.utility_cost, -25856786.959638454)
    assert np.allclose(units.get_heating_duty(), 198.1915651157178)
    assert np.allclose(units.get_cooling_duty(), 304.5184189296353)
    assert np.allclose(units.get_electricity_consumption(), 9.645038160426235)
    assert np.allclose(units.get_electricity_production(), 92.5193553387546)
    
def test_cornstover():
    """
    Test cornstover biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_cornstover()
    
    """
    from biorefineries.cornstover import cornstover_tea, ethanol
    MESP = cornstover_tea.solve_price(ethanol)
    units = UnitGroup('Biorefinery', cornstover_tea.units)
    assert np.allclose(MESP, 0.7387411012675781)
    assert np.allclose(cornstover_tea.sales, 134124825.10761952)
    assert np.allclose(cornstover_tea.material_cost, 82900270.63843909)
    assert np.allclose(cornstover_tea.installed_equipment_cost, 220353539.85009322)
    assert np.allclose(cornstover_tea.utility_cost, -11054921.042512074)
    assert np.allclose(units.get_heating_duty(), 318.0490631730075)
    assert np.allclose(units.get_cooling_duty(), 364.9938841327776)
    assert np.allclose(units.get_electricity_consumption(), 22.36637203703606)
    assert np.allclose(units.get_electricity_production(), 45.34818453627118)
    
if __name__ == '__main__':
    test_lipidcane()
    test_cornstover()
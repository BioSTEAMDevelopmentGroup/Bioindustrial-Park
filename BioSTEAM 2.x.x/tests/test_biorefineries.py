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
    assert np.allclose(lipidcane_tea.IRR, 0.17824969755468204)
    assert np.allclose(lipidcane_tea.sales, 102394506.53266351)
    assert np.allclose(lipidcane_tea.material_cost, 61762932.60602826)
    assert np.allclose(lipidcane_tea.installed_equipment_cost, 139562256.39522457)
    assert np.allclose(lipidcane_tea.utility_cost, -25857650.125692848)
    assert np.allclose(units.get_heating_duty(), 198.17981776084886)
    assert np.allclose(units.get_cooling_duty(), 304.52025279499816)
    assert np.allclose(units.get_electricity_consumption(), 9.645045283570326)
    assert np.allclose(units.get_electricity_production(), 92.52212901976533)
    
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
    assert np.allclose(MESP, 0.779430272130384)
    assert np.allclose(cornstover_tea.sales, 141622910.26714197)
    assert np.allclose(cornstover_tea.material_cost, 82900270.63843912)
    assert np.allclose(cornstover_tea.installed_equipment_cost, 220353539.85009325)
    assert np.allclose(cornstover_tea.utility_cost, -11054921.042512089)
    assert np.allclose(units.get_heating_duty(), 318.0490631730075)
    assert np.allclose(units.get_cooling_duty(), 364.9938841327776)
    assert np.allclose(units.get_electricity_consumption(), 22.36637203703606)
    assert np.allclose(units.get_electricity_production(), 45.3481845362712)
    
if __name__ == '__main__':
    test_lipidcane()
    test_cornstover()
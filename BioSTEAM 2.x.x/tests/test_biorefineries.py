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
    from biorefineries.lipidcane.system import lipidcane_tea
    units = UnitGroup('Biorefinery', lipidcane_tea.units)
    assert np.allclose(lipidcane_tea.IRR, 0.1777740209105262)
    assert np.allclose(lipidcane_tea.sales, 102751539.14512987)
    assert np.allclose(lipidcane_tea.material_cost, 62073249.769924864)
    assert np.allclose(lipidcane_tea.installed_cost, 139776480.90494874)
    assert np.allclose(lipidcane_tea.utility_cost, -25737178.078972895)
    assert np.allclose(units.get_heating_duty(), 199.5814871059657)
    assert np.allclose(units.get_cooling_duty(), 306.64848667499194)
    assert np.allclose(units.get_electricity_consumption(), 9.708269935954073)
    assert np.allclose(units.get_electricity_production(), 92.19922531727745)
    
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
    assert np.allclose(MESP, 0.7127573892429486)
    assert np.allclose(cornstover_tea.sales, 130722340.25312951)
    assert np.allclose(cornstover_tea.material_cost, 80542496.40653823)
    assert np.allclose(cornstover_tea.installed_cost, 204437262.02031183)
    assert np.allclose(cornstover_tea.utility_cost, -15306742.335879013)
    assert np.allclose(units.get_heating_duty(), 305.00593360849587)
    assert np.allclose(units.get_cooling_duty(), 318.85012928383793)
    assert np.allclose(units.get_electricity_consumption(), 18.40075443477195)
    assert np.allclose(units.get_electricity_production(), 50.22157712400753)
    
if __name__ == '__main__':
    test_lipidcane()
    test_cornstover()
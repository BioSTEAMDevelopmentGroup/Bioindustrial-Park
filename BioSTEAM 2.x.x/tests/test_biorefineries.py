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
from biosteam import default_utilities

__all__ = ('test_sugarcane', 
           'test_lipidcane', 
           'test_cornstover',
           'test_LAOs')

def test_sugarcane():
    """
    Test sugarcane biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_sugarcane()
    
    """
    default_utilities()
    from biorefineries.sugarcane import sugarcane_tea
    units = UnitGroup('Biorefinery', sugarcane_tea.units)
    assert np.allclose(sugarcane_tea.IRR, 0.10326317303317518)
    assert np.allclose(sugarcane_tea.sales, 87641108.77580193)
    assert np.allclose(sugarcane_tea.material_cost, 60096414.780067384)
    assert np.allclose(sugarcane_tea.installed_equipment_cost, 111490520.90754643)
    assert np.allclose(sugarcane_tea.utility_cost, -8595640.622981431)
    assert np.allclose(units.get_heating_duty(), 272.488831440168)
    assert np.allclose(units.get_cooling_duty(), 343.762814732534)
    assert np.allclose(units.get_electricity_consumption(), 9.827812210952846)
    assert np.allclose(units.get_electricity_production(), 37.37794241281641)

def test_lipidcane():
    """
    Test lipidcane biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_lipidcane()
    
    """
    default_utilities()
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
    default_utilities()
    from biorefineries.cornstover import cornstover_tea, ethanol
    MESP = cornstover_tea.solve_price(ethanol)
    units = UnitGroup('Biorefinery', cornstover_tea.units)
    assert np.allclose(MESP, 0.7382850997534475)
    assert np.allclose(cornstover_tea.sales, 134041032.28829786)
    assert np.allclose(cornstover_tea.material_cost, 82900270.63843909)
    assert np.allclose(cornstover_tea.installed_equipment_cost, 219990786.55181444)
    assert np.allclose(cornstover_tea.utility_cost, -11054921.042512089)
    assert np.allclose(units.get_heating_duty(), 318.0490631730075)
    assert np.allclose(units.get_cooling_duty(), 364.9938841327776)
    assert np.allclose(units.get_electricity_consumption(), 22.36637203703606)
    assert np.allclose(units.get_electricity_production(), 45.3481845362712)
 
def test_LAOs():
    """
    Test linear-alpha olefins biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_LAOs()
    
    """
    default_utilities()
    from biorefineries.LAOs import LAOs_tea, get_LAOs_MPSP
    MPSP = get_LAOs_MPSP()
    units = UnitGroup('Biorefinery', LAOs_tea.units)
    assert np.allclose(MPSP, 1386.9140838394856)
    assert np.allclose(LAOs_tea.sales, 194137353.07884187)
    assert np.allclose(LAOs_tea.material_cost, 163028024.89425507)
    assert np.allclose(LAOs_tea.installed_equipment_cost, 70249155.97098511)
    assert np.allclose(LAOs_tea.utility_cost, 3808618.491773257)
    assert np.allclose(units.get_heating_duty(), 61.35846036766088)
    assert np.allclose(units.get_cooling_duty(), 139.64639353406835)
    assert np.allclose(units.get_electricity_consumption(), 3.1670096267106196)
    assert np.allclose(units.get_electricity_production(), 3.167009626710621)   

if __name__ == '__main__':
    test_sugarcane()
    test_lipidcane()
    test_cornstover()
    test_LAOs()
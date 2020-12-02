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
from biosteam import main_flowsheet as find

__all__ = ('test_wheatstraw','test_bedding')

 
def test_wheatstraw():
    """
    Test wheatstraw biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_wheatstraw()
    
    """
    from biorefineries.wheatstraw import wheatstraw_tea, ethanol,wheatstraw_sys
    sys=wheatstraw_sys
    heating_set = set()
    cooling_set = set()
    power_set = set()
    
    for u in sys.units:
        for hu in u.heat_utilities:
            if hu.duty < 0:
                cooling_set.add(hu)
            else:
                heating_set.add(hu) 
        if u.ID == 'BT': continue   
        power_set.add(u.power_utility)        
        
        
    heating_duty = sum([i.duty/3600 for i in heating_set])
    cooling_duty = sum([-i.duty/3600 for i in cooling_set])
    power_consumption = sum([i.rate for i in power_set])
    power_production = find.unit.BT.get_design_result('Work','kW')
    
    MESP = wheatstraw_tea.solve_price(ethanol)
    
    assert np.allclose(MESP, 0.7709848574557352)
    assert np.allclose(wheatstraw_tea.sales, 110002799.85052818)
    assert np.allclose(wheatstraw_tea.material_cost, 62781438.64308582)
    assert np.allclose(wheatstraw_tea.installation_cost, 232107199.63202378)
    assert np.allclose(wheatstraw_tea.utility_cost, -10772692.43243473)
    assert np.allclose(heating_duty, 57768.68128136157)
    assert np.allclose(cooling_duty, 76946.63233675627)
    assert np.allclose(power_consumption, 21304.3438418358)
    assert np.allclose(power_production, 44340.18684321244)

def test_bedding():
    """
    Test bedding biorefinery. If all tests passed, no error is raised.
    
    Examples
    --------
    >>> test_bedding()
    
    """
    from biorefineries.bedding import bedding_tea, ethanol,bedding_sys
    sys=bedding_sys
    heating_set = set()
    cooling_set = set()
    power_set = set()
    
    for u in sys.units:
        
        for hu in u.heat_utilities:
            if hu.duty < 0:
                cooling_set.add(hu)
            else:
                heating_set.add(hu) 
        if u.ID == 'BT': continue        
        power_set.add(u.power_utility)        
        
        
    heating_duty = sum([i.duty/3600 for i in heating_set])
    cooling_duty = sum([-i.duty/3600 for i in cooling_set])
    power_consumption = sum([i.rate for i in power_set])
    power_production = find.unit.BT.get_design_result('Work','kW')
    
    MESP = bedding_tea.solve_price(ethanol)
    
    assert np.allclose(MESP, 0.4781603383903969)
    assert np.allclose(bedding_tea.sales, 100260507.63708565)
    assert np.allclose(bedding_tea.material_cost, 29991699.505678747)
    assert np.allclose(bedding_tea.installation_cost, 269786738.1692209)
    assert np.allclose(bedding_tea.utility_cost, 3060882.9816840775)
    assert np.allclose(heating_duty, 49755.30892975407)
    assert np.allclose(cooling_duty, 51140.490836284385)
    assert np.allclose(power_consumption, 31088.87688252176)
    assert np.allclose(power_production, 25178.04613625532)
    
if __name__ == '__main__':
    test_wheatstraw()
    test_bedding()
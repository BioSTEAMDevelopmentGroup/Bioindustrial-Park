# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:47:43 2025

@author: IGB
"""

import biosteam as bst
import biorefineries.nitric._units as _units

thermo = bst.Thermo([bst.Chemical('N2'), bst.Chemical('O2'), bst.Chemical('H2O'), bst.Chemical('HNO3'),])
bst.settings.set_thermo(thermo)

#%%

air_in = bst.Stream('air_in', N2=0.79, O2=0.21, total_flow=100, units='kg/hr', phase='g')

water_in = bst.Stream('water_in', H2O=1, total_flow=100, units='kg/hr', phase='l')

C101 = bst.IsothermalCompressor('C101', ins=air_in, outs='', P=9.3*101325,
                                eta=0.7, compressor_type='reciprocating')

R101 = _units.PlasmaReactor('R101', ins=(C101-0, water_in), outs=('vent', ''))

CWP = bst.ChilledWaterPackage('CWP')

U101 = _units.PowerUnit('U101', ins=R101-1, outs='product', power=10)
@U101.add_specification(run=True)
def adjust_U101_power():
    R101.run()
    U101.power = R101.power
    
sys_plasma = bst.main_flowsheet.create_system('sys_plasma')


                                

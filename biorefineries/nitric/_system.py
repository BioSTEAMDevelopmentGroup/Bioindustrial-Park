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

air_in = bst.Stream('air_in', N2=0.79, O2=0.21, units='m3/hr', phase='g')

water_in = bst.Stream('water_in', H2O=1, units='kg/hr', phase='l')

C101 = bst.IsothermalCompressor('C101', ins=air_in, outs='', P=9.3*101325,
                                eta=0.7, compressor_type='reciprocating')
    
# R101 = _units.PlasmaReactor('R101', ins=(C101-0, water_in), outs=('vent', ''), 
#                             tau=24, HNO3_scale=10000, concentration=4.778/62*63,
#                             electricity_consumption=0.2, electricity_to_heat_ratio=0.4,
#                             air_scale_ratio=1)

R101 = _units.PlasmaReactor('R101', ins=(C101-0, water_in), outs=('vent', ''),
                            HNO3_scale=10000, concentration=0.005, # in weight
                            electricity_consumption=15, cost_per_power=0.9, 
                            heat_by_electricity=0.4)
                            
@R101.add_specification(run=True)
def adjust_water_flow():
    water_in = R101.ins[1]
    water_in.imass['H2O'] = R101.HNO3_scale / 24 * (1 / R101.concentration - 1)


U101 = _units.PowerUnit('U101', ins=R101-1, outs='product', cost_per_power=1, power=10)
@U101.add_specification(run=True)
def adjust_U101_power():
    R101.run()
    U101.power = R101.power

@C101.add_specification(run=True)
def adjust_air_flow():
    air_in = C101.ins[0]
    U101.simulate()
    air_in.ivol['N2'] = 0.79 * U101.power / 10 * 0.025 * 60  # scaled by power; 25 lpm for 10 kW
    air_in.ivol['O2'] = 0.21 * U101.power / 10 * 0.025 * 60  # scaled by power; 25 lpm for 10 kW
    
sys_plasma = bst.main_flowsheet.create_system('sys_plasma')


                                

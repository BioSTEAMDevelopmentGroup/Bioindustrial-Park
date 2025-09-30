# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:47:43 2025

@author: IGB
"""

import biosteam as bst
import thermosteam as tmo
# import biorefineries.nitric._units as _units

thermo = bst.Thermo([bst.Chemical('N2'), bst.Chemical('O2'), bst.Chemical('H2O'),])
bst.settings.set_thermo(thermo)

#%%

air_in = bst.Stream('air_in', N2=0.79, O2=0.21, total_flow=100, units='kg/hr', phase='g')
water_in = bst.Stream('water_in', H2O=1, total_flow=100, units='kg/hr', phase='l')

C101 = bst.IsothermalCompressor('C101', ins=air_in, outs='', P=9.3*101325,
                                eta=0.7, compressor_type='reciprocating')


                                

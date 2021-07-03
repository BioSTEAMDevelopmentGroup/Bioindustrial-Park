# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 03:55:18 2021

@author: yrc2
"""

import thermosteam as tmo

OleylAlcohol = tmo.Chemical('OleylAlcohol', search_ID='143-28-2', default=True)
BDO = tmo.Chemical('BDO', search_ID='2,3-Butanediol')
Glucose = tmo.Chemical('Glucose')
Water = tmo.Chemical('Water')
Acetoin = tmo.Chemical('Acetoin', search_ID='3-Hydroxybutanone')
tmo.settings.set_thermo([OleylAlcohol, Acetoin, BDO, Glucose, Water], skip_checks=True)
print('BDO', 0.1)
s_partition = tmo.Stream(None, OleylAlcohol=1, Water=1, BDO=0.0001, Acetoin=0.10,
                         Glucose=0.01, units='kg/hr')
s_partition.lle(T=300.15)
K = tmo.separations.lle_partition_coefficients(s_partition['L'], s_partition['l'])
print(K)
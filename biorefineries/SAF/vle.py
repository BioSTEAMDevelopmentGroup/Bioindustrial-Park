#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 09:14:25 2023

@author: wenjun
"""

from thermosteam import indexer, equilibrium, settings

settings.set_thermo(['Ethylene'], cache=True)
imol = indexer.MolarFlowIndexer(
            l=[('Ethylene',0)],
            g=[('Ethylene', 40)])
vle = equilibrium.VLE(imol)
vle(T=298.15,P=199438)
vle


settings.set_thermo(['Ethylene','Ethane','Propene','Butene','H2','CH4','H2O'], cache=True)
imol = indexer.MolarFlowIndexer(
            l=[('Ethylene',0),('Ethane', 0),('Propene', 0),('Butene', 0),('H2', 0), ('CH4', 0),('H2O', 0)],
            g=[('Ethylene', 6.7),('Ethane', 0),('Propene', 0),('Butene', 0),('H2', 0), ('CH4', 0), ('H2O', 657)])
               
vle = equilibrium.VLE(imol)
vle(T=308.75,P=101325)
vle


settings.set_thermo(['Ethylene','Butadiene','C6H12','C8H16','1-decene','CH4','H2O'], cache=True)
imol = indexer.MolarFlowIndexer(
            l=[('Ethylene',0),('Butadiene', 0),('C6H12', 0),('C8H16', 0),('1-decene', 0), ('CH4', 0),('H2O', 0)],
            g=[('Ethylene', 4.42),('Butadiene', 337),('C6H12', 0.118),('C8H16', 0.341),('1-decene', 0.409), ('CH4', 0), ('H2O', 0)])
               
vle = equilibrium.VLE(imol)
vle(T=358.15,P=101325)
vle
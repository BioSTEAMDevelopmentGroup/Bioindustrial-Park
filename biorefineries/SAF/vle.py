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


settings.set_thermo(['Ethylene','Ethane','Propene','Butene'], cache=True)
imol = indexer.MolarFlowIndexer(
            l=[('Ethylene',0),('Ethane', 0),('Propene', 0),('Butene', 0)],
            g=[('Ethylene', 682),('Ethane', 1.86),('Propene', 0.138),('Butene', 1.73)])
vle = equilibrium.VLE(imol)
vle(T=248.15,P=27*101325)
vle


settings.set_thermo(['Ethylene','1-Octene','1-hexene','Butene'], cache=True)
imol = indexer.MolarFlowIndexer(
            l=[('Ethylene',0),('Butene', 0),('1-Octene', 0),('1-hexene', 0)],
            g=[('Ethylene', 20.5),('Butene', 265),('1-Octene', 3.31),('1-hexene', 39.7)])
vle = equilibrium.VLE(imol)
vle(T=358.15,P=3.5e6)
vle
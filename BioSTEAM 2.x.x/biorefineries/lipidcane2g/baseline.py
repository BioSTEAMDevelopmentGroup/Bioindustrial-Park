# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 02:14:56 2021

@author: yrc2
"""

from biorefineries import lipidcane as lc
lc.lipidcane_sys.diagram('cluster')

F_mass = lc.lipidcane.F_mass
lc.lipidcane.empty()
lc.lipidcane.imass['']
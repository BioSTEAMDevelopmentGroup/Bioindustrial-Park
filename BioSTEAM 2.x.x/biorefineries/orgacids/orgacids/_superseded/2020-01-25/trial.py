#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""

from biosteam import *
species = Species('Water')
species.Sugar = compounds.Substance('Sugar')
Stream.species = species
s1 = Stream('s1', Water=20)
Stream.species = Species('Water', 'Ethanol')
s2 = Stream('s2') # Note that s2 and s1 have different Species objects
J1 = units.Junction('J1', s1, s2)
#I added the ID since Junction has ID now
J1.simulate()
J1.show()
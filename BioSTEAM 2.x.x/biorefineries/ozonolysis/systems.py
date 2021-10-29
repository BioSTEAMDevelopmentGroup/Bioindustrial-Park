# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 08:18:19 2021

@author: yrc2
"""
from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals import chemicals
import biosteam as bst

bst.settings.set_thermo(chemicals)
feed = bst.Stream('feed', Water=100)
R1 = units.OzonolysisReactor('R1', feed)


ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram()
ozonolysis_sys.simulate()
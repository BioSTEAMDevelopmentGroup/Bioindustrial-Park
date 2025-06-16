# -*- coding: utf-8 -*-
"""
Created on 2025-06-05 18:57:00

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import biosteam as bst
from biorefineries.prefers.systems.LegH import create_LegH_system
from biorefineries.prefers import _chemicals as c
# %%
# Set up the global thermosteam environment
bst.settings.set_thermo(c.create_chemicals_LegH())

# Create the system using the factory
legh_sys = create_LegH_system()
# Simulate the system
legh_sys.simulate()

# View results
legh_sys.show()
legh_sys.diagram(format='html')
# print(s.LegH_Ingredients.F_mass)
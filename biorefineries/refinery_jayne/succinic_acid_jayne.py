#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 18:45:24 2022

@author: jayneallen
"""

#%% Setup

import biosteam as bst
from biosteam import settings, units
from biorefineries.refinery_jayne.crystallization_class_test_module_jayne import SuccinicAcidCrystallizer

#%% Streams

settings.set_thermo(['Glucose','Quicklime','Water','H2SO4'])

pure_glucose = bst.Stream(Glucose=1000.)
pure_glucose.show()

quicklime = bst.Stream(Quicklime=10000.,Water=100.)
quicklime.show()

sulfuric_acid = bst.Stream(H2SO4=30.)
sulfuric_acid.show()

elution_water = bst.Stream(Water=40.)
elution_water.show()

#%% Units

# Fermentation
R1 = units.BatchBioreactor('R1', ins=('pure_glucose','quicklime'), outs=('vent','fermentation_broth'),
                           tau = 0.95, N=8.)
# R1.simulate()     #abstract class and simulate batch bioreactor class (plug in def set up and def run)

# Solid-liquid separation
S1 = units.SolidsSeparator('S1', ins='fermentation_broth', outs=('crude_calcium_succinate','insoluble_cellular_materials'),
                           split=0.3)

# Acidulation
S2 = units.SolidsSeparator('S2', ins=('crude_calcium_succinate','sulfuric_acid'), outs=('crude_succinic_acid','gypsum'),
                           split=0.3)

# Mixer before ion exchange column
M1 = units.Mixer('M1', ins=('crude_succinic_acid','elution_water'),outs='')

# # Ion exchange column
# S3 = units.LiquidsSplitSettler('S3', ins=M1.outs[0], outs=('dilute_water_to_WWT','liquid_succinic_acid'),
#                                split=0.3)

# Crystallization
C1 = units.BatchCrystallizer('C1', ins=M1.outs[0], outs='wet_succinic_acid_crystals')

# Drying
S4 = units.RotaryVacuumFilter('S4', ins='wet_succinic_acid_crystals',outs=('dry_succinic_acid_crystals',
                                                                           'water_and_other_insolubles'),
                              split=0.5)

#%% Diagram

bst.main_flowsheet.diagram()













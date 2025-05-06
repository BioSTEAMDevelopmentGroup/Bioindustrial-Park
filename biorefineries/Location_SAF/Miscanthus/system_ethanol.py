#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:21:43 2024

@author: wenjun
"""

import biosteam as bst
import thermosteam
from biorefineries.cellulosic.chemicals import create_cellulosic_ethanol_chemicals
from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system


F = bst.Flowsheet('switchgrass_flor')
bst.main_flowsheet.set_flowsheet(F)

# #create chemicals
chem = create_cellulosic_ethanol_chemicals()
# #TODO: change extractives to extract
chem.set_synonym('Extract','Extractives')
bst.settings.set_thermo(chem, cache= True)


#create streams

swg = bst.Stream('switchgrass', # This is miscanthus, but I don't want to change it and get errors for now
Water=0.2,
Ash=0.0301,
Glucan=0.3472,
Xylan=0.2048,
Galactan=0.0112,
Arabinan=0.0136,
Lignin=0.1616,
Extract=0.0223,
Acetate=0.0034,
Protein=0.0058,
total_flow=104229.16,
units='kg/hr',
price=0.08)

#create system
sys_ethanol = create_cellulosic_ethanol_system('sys_switchgrass',ins = swg)

bst.rename_units([i for i in F.unit if i.ID[1] == '6'], 700)
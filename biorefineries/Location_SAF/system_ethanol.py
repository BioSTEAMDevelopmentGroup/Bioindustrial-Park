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
# SWITCHGRASS
swg = bst.Stream('switchgrass',
Arabinan=0.02789023841655421,
Galactan=0.010436347278452543,
Glucan=0.2717049032838507,
Xylan=0.21214574898785432,
Mannan=0.005937921727395412,
Lignin=0.17112010796221322,
Ash=0.016194331983805668,
Extractives=0.08457040035987407,
Water=0.2,
total_flow=104229.16,
units='kg/hr',
#price = 0.1)
price=0.08) #0.08 $/kg

#MISCANTHUS
# swg = bst.Stream('switchgrass',
# Water=0.2,
# Ash=0.0301,
# Glucan=0.3472,
# Xylan=0.2048,
# Galactan=0.0112,
# Arabinan=0.0136,
# Lignin=0.1616,
# Extract=0.0223,
# Acetate=0.0034,
# Protein=0.0058,
# total_flow=104229.16,
# units='kg/hr',
# price=0.08)

#create system
sys_ethanol = create_cellulosic_ethanol_system('sys_switchgrass',ins = swg)

bst.rename_units([i for i in F.unit if i.ID[1] == '6'], 700)

# F.switchgrass # to access stream
# F.ethanol # to access stream

# ethanol = F.ethanol

# _gal_per_m3 = 1000/3.785412
# get_ethanol_conversion = lambda: ethanol.F_vol * _gal_per_m3 / (F.switchgrass.F_mass * 0.8/1000)

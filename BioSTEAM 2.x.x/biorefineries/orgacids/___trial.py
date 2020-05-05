#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.16.7
thermosteam v0.16.5
'''

import biosteam as bst
# from Sarang.system import *
from orgacids.system import *

find = bst.main_flowsheet
broth_for_separation = find.stream.broth_for_separation
# MPSP = 0
System = bst.System
System.maxiter = 1000
System.molar_tolerance = 5

for i in range(1):
    # try:
    System.converge_method = 'Fixed-point'
    orgacids_sys.simulate()
for i in range(10):
    MPSP = lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    # except: 
    #     System.converge_method = 'Wegstein'
    #     orgacids_sys.simulate()
    #     print(F_pre_S404.results())
    #     MPSP = lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    # print(F401.outs[1].F_mass)
    # Right now MPSP changes a bit each time the system simulates due to separation_sys
    # starting F401, the outs change each simulation
    
print('MPSP = %s $/kg'%MPSP)
# orgacids_sys.diagram()

# S405 = bst.units.ShortcutColumn('S404', ins=S404-1,
#                                 LHK=('LacticAcid', 'H2SO4'),
#                                 product_specification_format='Recovery',
#                                 Lr=0.95, Hr=0.85, k=1.2)
# S405.simulate()
# S405.show(N=100)

# S404.show(N=100)
# tot_mass = 0

# for species in S404.outs[1].chemicals:
#     if species.phase_ref=='l' and S404.outs[1].imass[species.ID]>0:
        
#         print(species.ID)
#         print(species.Tb)
#         if not species.ID == 'H2SO4':
#             tot_mass +=S404.outs[1].imass[species.ID]
# LA_purity = round(100*S404.outs[1].imass['LacticAcid']/tot_mass)

LA_purity = round(100*S404.outs[1].imass['LacticAcid']/S404.outs[1].F_mass)
LA_yield = round(100*S404.outs[1].imol['LacticAcid']/(M401.outs[0].imol['LacticAcid'] + 2*M401.outs[0].imol['CalciumLactate']))
print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(LA_purity, LA_yield, round(MPSP,3)))

#!!! M203.heat_utilities raises error
#!!! save_report not saving diagrams
# orgacids_sys.save_report('orgacids_sys.xlsx')



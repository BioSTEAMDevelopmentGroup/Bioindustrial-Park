#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.17.1
thermosteam v0.17.1
'''

import biosteam as bst
# from Sarang.system import *
from orgacids.system import *


# %%

# find = bst.main_flowsheet
# broth_for_separation = find.stream.broth_for_separation
# # MPSP = 0
# System = bst.System
# System.maxiter = 1000
# System.molar_tolerance = 5

System.converge_method = 'Fixed-point'

for i in range(10):
    orgacids_sys.simulate()
    for i in range(10):
        MPSP = lactic_acid.price = \
            orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)

MPSP = lactic_acid.price
print('MPSP = %s $/kg'%MPSP)
# orgacids_sys.diagram()


LA_purity = round(100*S404.outs[1].imass['LacticAcid']/S404.outs[1].F_mass)
LA_yield = round(100*S404.outs[1].imol['LacticAcid']/(M401.outs[0].imol['LacticAcid'] + 2*M401.outs[0].imol['CalciumLactate']))
print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(LA_purity, LA_yield, round(MPSP,3)))

# orgacids_sys.save_report('orgacids_sys.xlsx')






# %%























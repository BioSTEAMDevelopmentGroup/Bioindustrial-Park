#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.18.0
thermosteam v0.18.0
'''

import biosteam as bst
from orgacids.system import *

System = bst.System
System.converge_method = 'Fixed-point'
System.maxiter = 1000
System.molar_tolerance = 1

for i in range(5):
    orgacids_sys.simulate()
for i in range(5):
    MSP = lactic_acid.price = orgacids_tea.solve_price(
        lactic_acid, orgacids_no_BT_tea)

purity = lactic_acid.imass['LacticAcid'] / lactic_acid.F_mass
print(f'\nPurity is {purity:.1%}')

broth = R301.outs[0]
produced_lactic_acid = broth.imol['LacticAcid'] + 2*broth.imol['CalciumLactate']
separation_yield = lactic_acid.imol['LacticAcid'] / produced_lactic_acid
print(f'\nSeparation yield is {separation_yield:.1%}')
print(f'\nMinimum selling price is {MSP:.3f} $/kg')


from orgacids.utils import get_sugar_conc

R301 = bst.main_flowsheet.unit.R301
sugar_conc = get_sugar_conc(R301.saccharified_stream, ('Glucose','Xylose'))
print(f'\nSugar concentration of R301 after saccharification is {sugar_conc:.1f} g/L')



















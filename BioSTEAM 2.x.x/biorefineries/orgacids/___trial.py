#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.17.7
thermosteam v0.17.15
'''

import biosteam as bst
# from Sarang.system import *
from orgacids.system import *


# %% MPSP

MPSP = lactic_acid.price
# orgacids_sys.diagram()


LA_purity = round(100*S404.outs[1].imass['LacticAcid']/S404.outs[1].F_mass)
LA_yield = round(100*S404.outs[1].imol['LacticAcid']/(M401.outs[0].imol['LacticAcid'] + 2*M401.outs[0].imol['CalciumLactate']))
print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(LA_purity, LA_yield, round(MPSP,3)))

# orgacids_sys.save_report('orgacids_sys.xlsx')


# %% Sugar concentration

from orgacids.utils import get_sugar_conc
sugar_conc = get_sugar_conc(R301.saccharified_stream, ('Glucose','Xylose'))
print(f'Sugar concentration of R301 after saccharification is {sugar_conc}')






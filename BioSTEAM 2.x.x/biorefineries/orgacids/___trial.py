#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.12.6
thermosteam v0.12.4
'''

import biosteam as bst
from orgacids.system import *

find = bst.main_flowsheet
broth_for_separation = find.stream.broth_for_separation

# # Right now minimum selling price changes a bit each time the system simulates 
# # due to separation_sys starting F401, the outs change each simulation
# print(F401.outs[1].F_mass)

# for i in range(3): orgacids_sys.simulate()
# for i in range(3): 
#     lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_BT_tea)
    
print(lactic_acid.price)

# orgacids_sys.save_report('orgacids_sys.xlsx')





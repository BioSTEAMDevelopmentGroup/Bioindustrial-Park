#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.11.3
thermosteam v0.11.2
'''

import biosteam as bst
from orgacids.system import *

find = bst.main_flowsheet
broth_for_separation = find.stream.broth_for_separation

for i in range(5):
    orgacids_sys.simulate()
    # print(F401.outs[1].F_mass)
    # Right now MPSP changes a bit each time the system simulates due to separation_sys
    # starting F401, the outs change each simulation
    MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    print(MPSP)

# orgacids_sys.save_report('orgacids_sys.xlsx')





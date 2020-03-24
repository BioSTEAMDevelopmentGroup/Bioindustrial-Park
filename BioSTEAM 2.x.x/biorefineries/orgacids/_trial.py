#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Set up

from orgacids.system import *

for i in range(0, 10):
    orgacids_sys.simulate()
    MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    print(MPSP)
    

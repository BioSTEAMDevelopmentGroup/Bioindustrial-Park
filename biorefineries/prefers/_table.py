# -*- coding: utf-8 -*-
"""
Created on 2025-09-19 14:51:38

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import numpy as np
import pandas as pd
import thermosteam as tmo
import biosteam as bst
from biosteam._heat_utility import HeatUtility

DataFrame = pd.DataFrame
ExcelWriter = pd.ExcelWriter

def all_cost_table(tea):
    
    columns = ('Unit operation',
               'Purchase cost (10^6 USD)',
               'Utility cost (10^6 USD/yr)',
               'Installed cost (10^6 USD)',
               'Installing cost (10^6 USD)',)
    units = sorted([i for i in tea.system.units if i._design or i._cost], key=lambda x: x.line)
    operating_days = tea.operating_days
    N_units = len(units)
    array = np.empty((N_units, 5), dtype=object)
    IDs = []
    types = array[0:, 0]
    C_cap = array[0:, 1]
    C_op = array[0:, 2]
    C_inst = array[0:, 3]
    C_in = array[0:, 4]
    # Get data
    for i in range(N_units):
        unit = units[i]
        types[i] = unit.line
        C_cap[i] = unit.purchase_cost / 1e6
        C_op[i] = unit.utility_cost * operating_days * 24  / 1e6
        C_inst[i] = unit.installed_cost / 1e6
        C_in[i] = C_inst[i] - C_cap[i]
        IDs.append(unit.ID)
    
    df = DataFrame(array, columns=columns, index=IDs)    
    if not tea.lang_factor:
        df['Installed cost (10^6 USD)'] = [u.installed_cost / 1e6 for u in units]
    
    return df
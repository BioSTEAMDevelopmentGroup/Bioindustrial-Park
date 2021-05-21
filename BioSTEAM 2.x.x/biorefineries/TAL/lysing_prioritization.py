# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:11:59 2021

@author: sarangbhagwat
"""

from biorefineries.TAL.system_split import *
import numpy as np
import pandas as pd
# from matplotlib import pyplot as plt
from datetime import datetime

spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

steps = 8
cellmass_TAL_retention_gpg_array = np.linspace(0., 0.5, steps)

S404_s.cellmass_TAL_retention_gpg = 0.
no_retention_cost = get_SA_MPSP()

additional_costs = []
for retention in cellmass_TAL_retention_gpg_array:
    S404_s.cellmass_TAL_retention_gpg = retention
    additional_costs.append(get_SA_MPSP() - no_retention_cost)

#%% Plot
df = pd.DataFrame({'waste all retained TAL':additional_costs}, index = cellmass_TAL_retention_gpg_array)
df.plot.line(style={'waste all retained TAL':'r--'})

# %% Save as excel file
dateTimeObj = datetime.now()
file_to_save = 'lysing_prioritization_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

with pd.ExcelWriter(file_to_save+'.xlsx') as writer:
    
    df.transpose().to_excel(writer, sheet_name='Lysing priorotization')
    
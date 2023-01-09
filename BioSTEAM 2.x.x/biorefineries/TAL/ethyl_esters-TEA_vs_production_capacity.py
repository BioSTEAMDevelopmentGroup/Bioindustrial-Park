# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 15:38:35 2023

@author: sarangbhagwat
"""

import numpy as np
from matplotlib import pyplot as plt
from pandas import DataFrame
print('\nImporting and running biorefinery (single-point) ...\n')
from biorefineries.TAL.system_ethyl_esters import get_Mixed_esters_MPSP, TAL_tea, s

#%% Setup
print('\nSetting up ...\n')
feedstock = s.feedstock
Mixed_esters = s.Mixed_esters
feedstock_Fmasses = np.linspace(feedstock.F_mass*0.01, feedstock.F_mass*2, 100)
prod_Fmasses = []
MPSPs = []
TCIs = []
AOCs = []

#%% Simulate
print('\nSimulating ...\n')
ct = 0
for ffm in feedstock_Fmasses:
    ct +=1
    feedstock.F_mass = ffm
    MPSPs.append(get_Mixed_esters_MPSP())
    TCIs.append(TAL_tea.TCI/1e6)
    AOCs.append(TAL_tea.AOC/1e6)
    prod_Fmasses.append(Mixed_esters.F_mass*TAL_tea.operating_hours/1e3)
    if ct%10 == 0:
        print(f'{ct}/{len(feedstock_Fmasses)} ...')
if not ct%10 == 0:
    print(f'{ct}/{len(feedstock_Fmasses)} ...')

#%% Plot
print('\nPlotting ...\n')
plt.plot(prod_Fmasses, MPSPs)
plt.show()

plt.plot(prod_Fmasses, TCIs)
plt.show()

plt.plot(prod_Fmasses, AOCs)
plt.show()

#%% Create dataframe
print('\nCreating dataframe ...\n')
df_cap = DataFrame(data={'Production capacity [MT,y]': prod_Fmasses,
                         'MPSP [$/kg-mixed_esters]': MPSPs,
                         'TCI [M$]': TCIs,
                         'AOC [M$/y]': AOCs,})

#%% Save as .xlsx
filename = 'MPSP_mixed_ethyl_esters_vs_capacity.xlsx'
print('\nSaving as excel file ...\n')
df_cap.to_excel('MPSP_mixed_ethyl_esters_vs_capacity.xlsx')

print(f'\nSaved as {filename}.')
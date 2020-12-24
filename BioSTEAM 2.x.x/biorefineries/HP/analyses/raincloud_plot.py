# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 10:30:18 2020

@author: saran
"""
# %% Imports

from matplotlib import pyplot as plt
import seaborn
import pandas as pd
from ptitprince import PtitPrince as pt

#%% Design the plot
# fig, ax = plt.subplots(figsize=(7, 5))

# # Data
# # array_data = []
# PANDAS_DATAFRAME = pd.read_csv("just_MPSP.csv", sep= ",")
# CATEGORICAL_INFO="MPSP"
# ACTUAL_VALUES= "Biorefinery"

# # # Labels
# # XLABEL=
# # YLABEL=

# #%% Plot the plot
# pt.RainCloud(data=PANDAS_DATAFRAME, x=CATEGORICAL_INFO, y=ACTUAL_VALUES,  ax=ax, orient='h')
#              # palette=[LIST_OF_COLORS])
 
# # Change some aspects of the plot to make it look prettier
# # ax.set(ylabel=YLABEL, xlabel=XLABEL)
# ax.tick_params(axis='x', which='major', bottom=True, direction='in')
# ax.tick_params(axis='y', which='major', left=True, direction='in')
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)

# # Save the plot
# plt.savefig('raincloud.png', transparent=True, dpi=300)


# %% Alternative plot
dx = "Biorefinery"; dy = "MPSP"; ort = "h"; pal = "Set3"; sigma = .2
df = pd.read_csv("just_MPSP.csv", sep= ",")
f, ax = plt.subplots(figsize=(7, 5))

ax=pt.RainCloud(x = dx, y = dy, data = df, palette = pal, bw = sigma,
                 width_viol = .3, ax = ax, orient = ort, move = .2)

# plt.title("Figure P8\n Rainclouds with Shifted Rain")

plt.savefig('raincloud.png', transparent=True, dpi=300)
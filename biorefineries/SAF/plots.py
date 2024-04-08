#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 20:41:38 2024

@author: wenjun
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import biosteam as bst
from biosteam import preferences
from biosteam import report
from biosteam.plots.utils import CABBI_green_colormap
from thermosteam.utils.colors import *
#import contourplots
from biorefineries.SAF.systems import sys, F, process_groups_dict, process_groups
#from biorefineries.SAF.models import *

#%%
# =============================================================================
# Sensitivity
# =============================================================================
# 1. MPSP
def MSP_sensitivity_plot_no_CCS():
    bst.plots.plot_spearman_1d(rhos=(-0.145221138,
                                     -0.173385778,
                                     -0.177101627,
                                     -0.177816292,
                                     -0.190147753,
                                     0.19169614,
                                     -0.215473368,
                                     0.30783141,
                                     0.440609761,
                                     0.623233098), 
                               index = ['Enzymatic hydrolysis solids loading',
                                        'Enzymatic hydrolysis glucan-to-glucose',
                                        'Boiler efficiency',
                                        'Gasoline price',
                                        'Bagasse split for ethanol',
                                        'Enzyme price',
                                        'Plant uptime',
                                        'Feedstock price',
                                        'Enzyme loading',
                                        'TCI ratio'],
                               xlabel= 'MSP',
                               color = GG_colors.blue.RGBn,edgecolors = 'black', sort = False, w = 1./2)
    



def MSP_sensitivity_plot_CCS():

    bst.plots.plot_spearman_1d(rhos=(-0.122046053,
                                     -0.130803244,
                                     -0.14547189,
                                     0.154736631,
                                     -0.174742754,
                                     0.204332532,
                                     0.305287507,
                                     -0.445919381,
                                     -0.454413025,
                                     0.50072789,),
                               index=['fermentation glucose-to-ethanol',
                                        'plant uptime',
                                        'boiler efficiency',
                                        'feedstock unit price',
                                        'feedstock flow rate',
                                        'enzyme loading',
                                        'TCI ratio',
                                        'CO2 transportation and storage unit cost'
                                        'bagasse split raio for ethanol',
                                        'CO2 capture ratio'],
                               xlabel= 'MSP',
                               color = GG_colors.blue.RGBn,edgecolors = 'black', sort = False, w = 1./2)
# 2. GWP
def GWP_sensitivity_plot():
    bst.plots.plot_spearman_1d(rhos=(0.045471445,
                                     0.046671414,
                                     0.055375035,
                                     0.059740621,
                                     -0.069694718,
                                     -0.071422701,
                                     -0.087769537,
                                     -0.285025321,
                                     0.568716744,
                                     0.74068337),    
                               index = ['Enzyme GWP',
                                        'Natural gas GWP',
                                        'Fermentation glucose-to-ethanol',
                                        'Enzymatic hydrolysis glucan-to-glucose',
                                        '1st oligomerization residence time',
                                        'Enzyme loading',
                                        'Pretreatment solids loading',
                                        'Boiler efficiency',
                                        'Feedstock GWP',
                                        'Bagasse split ratio for ethanol'],
                               xlabel= 'GWP Jet fuel',
                               color = GG_colors.blue.RGBn,edgecolors = 'black', sort = False, w = 1./2)




import pandas as pd, numpy as np, matplotlib.pyplot as plt
from qsdsan.utils import palettes
from colorpalette import Color
from matplotlib.mathtext import _mathtext as mathtext

# color palette
Guest = palettes['Guest']
b = Guest.blue.HEX
g = Guest.green.HEX
r = Guest.red.HEX
o = Guest.orange.HEX
y = Guest.yellow.HEX
a = Guest.gray.HEX
p = Guest.purple.HEX

db = Color('dark_blue', (53, 118, 127)).HEX
dg = Color('dark_green', (77, 126, 83)).HEX
dr = Color('dark_red', (156, 75, 80)).HEX
do = Color('dark_orange', (167, 95, 62)).HEX
dy = Color('dark_yellow', (171, 137, 55)).HEX
da = Color('dark_gray', (78, 78, 78)).HEX
dp = Color('dark_purple', (76, 56, 90)).HEX

path = 'path'

centralized = pd.read_excel(path + 'centralized')

Al_Mazarih = pd.read_excel(path + 'Al_Mazarih')

Kisumu = pd.read_excel(path + 'Kisumu')

Malindi = pd.read_excel(path + 'Malindi')

Sweto = pd.read_excel(path + 'Sweto')

Dakar = pd.read_excel(path + 'Dakar')

Kampala = pd.read_excel(path + 'Kampala')

Alibagh = pd.read_excel(path + 'Alibagh')

Sharjah = pd.read_excel(path + 'Sharjah')

#%% visualization (7 cities)

fig, ax = plt.subplots(figsize=(21, 8))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 42])
ax.set_ylim([0, 0.7])

ax.tick_params(direction='inout', length=15, width=1.5, bottom=False, top=False, left=True, right=False, pad=0)

ax_right = ax.twinx()
ax_right.set_xlim([0, 42])
ax_right.set_ylim([0, 0.7])
ax_right.tick_params(direction='in', length=7.5, width=1.5, bottom=False, top=True, left=False, right=True, labelcolor='none')

mathtext.FontConstantsBase.sup1 = 0.35

def plot_city(city, ID, start_position):
    centralized_data = centralized[ID]

    box1 = ax_right.boxplot(centralized_data,
                            positions=[start_position,],
                            widths=[0.8,],
                            patch_artist=True,
                            showfliers=False)
                            # flierprops={'marker':'o',
                            #             'markersize':3,
                            #             'markerfacecolor':g}
    
    for box in box1['boxes']:
        box.set(color='k', facecolor=g, linewidth=1.5)
        
    for whisker in box1['whiskers']:
        whisker.set(color='k', linewidth=1.5)
    
    for median in box1['medians']:
        median.set(color='k', linewidth=1.5)
        
    for cap in box1['caps']:
        cap.set(color='k', linewidth=1.5)
    
    box2 = ax_right.boxplot(city[f'{ID}_br_A'],
                            positions=[start_position+1,],
                            widths=[0.8,],
                            patch_artist=True,
                            showfliers=False)
    
    for box in box2['boxes']:
        box.set(color='k', facecolor=r, linewidth=1.5)
        
    for whisker in box2['whiskers']:
        whisker.set(color='k', linewidth=1.5)
    
    for median in box2['medians']:
        median.set(color='k', linewidth=1.5)
        
    for cap in box2['caps']:
        cap.set(color='k', linewidth=1.5)
        
    box3 = ax_right.boxplot(city[f'{ID}_br_B'],
                            positions=[start_position+2,],
                            widths=[0.8,],
                            patch_artist=True,
                            showfliers=False)
    
    for box in box3['boxes']:
        box.set(color='k', facecolor=p, linewidth=1.5)
        
    for whisker in box3['whiskers']:
        whisker.set(color='k', linewidth=1.5)
    
    for median in box3['medians']:
        median.set(color='k', linewidth=1.5)
        
    for cap in box3['caps']:
        cap.set(color='k', linewidth=1.5)
        
    box4 = ax_right.boxplot(city[f'{ID}_ng_B'],
                            positions=[start_position+3,],
                            widths=[0.8,],
                            patch_artist=True,
                            showfliers=False)
    
    for box in box4['boxes']:
        box.set(color='k', facecolor=a, linewidth=1.5)
        
    for whisker in box4['whiskers']:
        whisker.set(color='k', linewidth=1.5)
    
    for median in box4['medians']:
        median.set(color='k', linewidth=1.5)
        
    for cap in box4['caps']:
        cap.set(color='k', linewidth=1.5)
    
    box5 = ax_right.boxplot(city[f'{ID}_re_B'],
                            positions=[start_position+4,],
                            widths=[0.8,],
                            patch_artist=True,
                            showfliers=False)
    
    for box in box5['boxes']:
        box.set(color='k', facecolor=b, linewidth=1.5)
        
    for whisker in box5['whiskers']:
        whisker.set(color='k', linewidth=1.5)
    
    for median in box5['medians']:
        median.set(color='k', linewidth=1.5)
        
    for cap in box5['caps']:
        cap.set(color='k', linewidth=1.5)

plot_city(Al_Mazarih, 'Al_Mazarih', 1)
plot_city(Kisumu, 'Kisumu', 7)
plot_city(Malindi, 'Malindi', 13)
plot_city(Sweto, 'Sweto', 19)
plot_city(Dakar, 'Dakar', 25)
plot_city(Alibagh, 'Alibagh', 31)
plot_city(Sharjah, 'Sharjah', 37)

for x_position in np.arange(6, 42, 6):
    plt.plot([x_position, x_position],
             [0, 1.75],
             lw=1.5, color='k', zorder=1)

plt.xticks(np.arange(3, 41, 6),
           ['Al Mazarih','Kisumu','Malindi','Sweto','Dakar','Alibagh','Sharjah'])

ax.set_ylabel(r'$\mathbf{Daily\ user\ fee}$ [\$·cap${^{-1}}$·day${^{-1}}$]', fontname='Arial', fontsize=30, labelpad=0)

#%% visualization (Kampala)

fig, ax = plt.subplots(figsize=(3, 8))

plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['hatch.linewidth'] = 1.5
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25

plt.xticks(fontname='Arial')
plt.yticks(fontname='Arial')

plt.rcParams.update({'mathtext.fontset': 'custom'})
plt.rcParams.update({'mathtext.default': 'regular'})
plt.rcParams.update({'mathtext.bf': 'Arial: bold'})

ax = plt.gca()
ax.set_xlim([0, 6])
ax.set_ylim([0, 1.75])

ax.tick_params(direction='inout', length=15, width=1.5, bottom=False, top=False, left=True, right=False, pad=0)

ax_right = ax.twinx()
ax_right.set_xlim([0, 6])
ax_right.set_ylim([0, 1.75])
ax_right.tick_params(direction='in', length=7.5, width=1.5, bottom=False, top=True, left=False, right=True, labelcolor='none')

mathtext.FontConstantsBase.sup1 = 0.35

plot_city(Kampala, 'Kampala', 1)

plt.xticks([3,], ['Kampala',])

ax.set_ylabel(r'$\mathbf{Daily\ user\ fee}$ [\$·cap${^{-1}}$·day${^{-1}}$]', fontname='Arial', fontsize=30, labelpad=0)
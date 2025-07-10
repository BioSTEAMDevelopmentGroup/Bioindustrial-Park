#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

# Run this cell first
from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors
import numpy as np

from biorefineries import HP
# from biorefineries.HP.systems.system_HP_solubility_exploit_ethanol_sugarcane import HP_tea, lca, R302, spec, HP_product, simulate_and_print, theoretical_max_g_HP_per_g_SA

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os

# from biorefineries.HP.models import models_HP_solubility_exploit as models

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')
# bst.speed_up()

#%%
HP.load_model('glucose', 'acrylic', '300 L')

system = HP.system
f = system.flowsheet
s,u = f.stream, f.unit

simulate_and_print = HP.simulate_and_print
get_adjusted_MSP = HP.get_adjusted_MSP
tea = HP.tea
lca = HP.lca

#%% Filepaths
HP_filepath = HP.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
HP_results_filepath = HP_filepath + '\\analyses\\results\\'

#%%
simulate_and_print()

#%%\

#%%
steps = 100
V_wfs = np.linspace(0.1, 0.99, steps)
MPSPs, TCIs, AOCs, GWPs, FECs = [], [], [], [], []

for i, V_wf in zip(range(steps), V_wfs):
    print(f'{i+1}/{steps}')
    u.R302.V_wf = V_wf
    # simulate_and_print()
    system.simulate()
    MPSPs.append(get_adjusted_MSP())
    TCIs.append(tea.TCI)
    AOCs.append(tea.AOC)
    GWPs.append(lca.GWP)
    FECs.append(lca.FEC)

#%%

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
def plot_multiple_metrics(x_axis_list, y_axis_list_of_lists, 
                          xlim=(0,100), 
                          ylims_list=[],
                          xticks=[0 + 10*i for i in range(0,11)],
                          yticks_list=[],
                          save_fig=True,
                          filename='supernatant_recycling_fixed_decarboxylation_plot.png',
                          ylabels_list=[r"$\bfHP$"+" " +r"$\bfRecovery$",
                                        r"$\bfMPSP$",
                                        r"$\mathrm{\bfGWP}_{\bf100}$",
                                        r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$",
                                        r"$\bfFEC$" ,
                                        ],
                          xlabel=r"$\bfSupernatant$"+" " +r"$\bfRecycling$",
                          y_axis_units_list= [r"$\mathrm{\%}$", 
                                              r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                                              r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                                              r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$",
                                              ],
                          x_axis_units=r"$\mathrm{\%}$",
                          units_on_newline=False, # (x, y)
                          fontname="Arial Unicode",
                          fontsize="12",
                          y_plotline_colors_list=[ 
                                                   "black",
                                                   "#A97802",
                                                   "blue",
                                                  '#607429',
                                                  '#A100A1',
                                                  ],
                          dpi=600,
                          bbox_inches='tight',
                          ):

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=1)
    plt.rcParams['font.sans-serif'] = fontname
    plt.rcParams['font.size'] = fontsize
    plt.xlim(xlim[0], xlim[1])
    
    p1, = ax.plot(x_axis_list, y_axis_list_of_lists[0],
                  color=y_plotline_colors_list[0],
                  zorder=1)
    
    if ylims_list:
        ax.set_ylim(ylims_list[0])
    ax.yaxis.label.set_color(p1.get_color())
    
    units_opening_brackets = ["[" , "[" ]
    for j in range(len(units_opening_brackets)):
        if units_on_newline:
            units_opening_brackets[j] = "\n[" 
            
    ax.set_ylabel(ylabels_list[0] +" " + units_opening_brackets[1] + y_axis_units_list[0] + "]")
    ax.tick_params(axis='y', colors=p1.get_color(), which='both', direction='inout')
    ax.spines['right'].set_color(p1.get_color())
    
    ax.set_yticks(yticks_list[0])
    ax.tick_params(axis='x', direction='inout')
    ax.tick_params(axis='x', which='minor', direction='inout')
    ax.tick_params(axis='y', which='minor', direction='inout')
    
    # Add minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    for i in range(1, len(y_axis_list_of_lists)):
        twini = ax.twinx()
        pi, = twini.plot(x_axis_list, y_axis_list_of_lists[i],
                         color=y_plotline_colors_list[i], 
                         zorder=i+1)
        
        twini.yaxis.label.set_color(pi.get_color())
        twini.spines['right'].set_color(pi.get_color())
        twini.spines.right.set_position(("axes", 1 + 0.2*(i-1)))
        twini.tick_params(axis='y', colors=pi.get_color(), which='both', direction='inout')
        
        if units_on_newline:
            units_opening_brackets = ["\n[", "\n["]
        else:
            units_opening_brackets = ["[", "["]

        twini.set_ylabel(ylabels_list[i] +" " + units_opening_brackets[1] + y_axis_units_list[i] + "]")
        
        if ylims_list:
            twini.set_ylim(ylims_list[i][0], ylims_list[i][1])
        
        twini.set_yticks(yticks_list[i])
        twini.yaxis.set_minor_locator(AutoMinorLocator())  # <-- Minor ticks for twin axes

    ax.set_xlabel(xlabel +" " + units_opening_brackets[0] + x_axis_units + "]")
    ax.set_xticks(xticks)
    
    if save_fig:
        plt.savefig(filename, dpi=dpi, bbox_inches=bbox_inches)
    
    plt.show()
    
#%%
plot_multiple_metrics(100*np.array(V_wfs),
                      [
                       MPSPs,
                       np.array(TCIs)/1e6,
                       np.array(AOCs)/1e6,
                       GWPs, 
                       FECs,
                       ],
                      
                      xlim=(0,100), 
                      ylims_list=[
                                  (0,2.5),
                                  (0,1200),
                                  (0,200),
                                  (0,5),
                                  (0,50)
                                  ],
                      xticks=[0 + 10*i for i in range(0,11)],
                      yticks_list = [
                                 np.arange(0, 3.51, 0.5),
                                 np.arange(0, 2000.1, 250),
                                 np.arange(0, 250.1, 50),
                                 np.arange(0, 6.1, 1),
                                 np.arange(0, 60.1, 10),
                                 ],
                      
                      xlabel=r"$\bfFermentation$"+" " +r"$\bfWorking$"+" " +r"$\bfVolume$"+" " +r"$\bfFraction$",
                      
                      ylabels_list=[
                                    r"$\bfMPSP$",
                                    r"$\bfTCI$",
                                    r"$\bfAOC$",
                                    r"$\bfCI$",
                                    # r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$",
                                    r"$\bfFEC$",
                                    ],
                      
                      # y_axis_units_list= [
                      #                     r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                      #                     r"$\mathrm{MM\$}$", 
                      #                     r"$\mathrm{MM\$}\cdot\mathrm{y}^{-1}$",
                      #                     r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                      #                     r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$",
                      #                     ]
                      
                      y_axis_units_list= [
                                          r"$\mathrm{\$/kg}$",
                                          r"$\mathrm{MM\$}$", 
                                          r"$\mathrm{MM\$/y}$",
                                          r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq./kg}$",
                                          r"$\mathrm{MJ/kg}$",
                                          ]
                      
                      )

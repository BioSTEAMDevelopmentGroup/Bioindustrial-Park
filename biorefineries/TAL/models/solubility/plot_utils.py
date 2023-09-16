#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2023-2024, Sarang Bhagwat <sarangb2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Created on Fri Aug 18 17:32:13 2023

@author: sarangbhagwat
"""

from matplotlib import pyplot as plt
import numpy as np

def plot_solubility_model(experimental_Ts, experimental_solubilities, 
                          get_TAL_solubility_in_water_gpL,
                          model_Ts=273.15 + np.linspace(0., 99., 50),
                          R_squared=None,
                          xlim=(0,100), ylim=(0,200),
                          xticks=[0 + 10*i for i in range(0,11)],
                          save_fig=True,
                          filename='solubility_plot.png',
                          ylabel=r"$\bfTAL$" +" " +r"$\bfSolubility$" +" " + r"$\bfin$" +" " + r"$\bfWater$",
                          xlabel=r"$\bfTemperature$",
                          solubility_units= r"$\mathrm{g-TAL}\cdot\mathrm{L-water}^{-1}$",
                          temperature_units="\u00b0" + r"$\mathrm{C}$",
                          units_on_newline=(False, True), # (x, y)
                          legend_loc='center left',
                          fontname="Arial Unicode",
                          fontsize="12",
                          predicted_plot_color='blue',
                          experimental_plot_color='orange',
                          experimental_plot_marker='x',
                          dpi=600,
                          bbox_inches='tight',
                          ):


    plt.rcParams['font.sans-serif'] = fontname
    plt.rcParams['font.size'] = fontsize
    
    model_solubilities = [get_TAL_solubility_in_water_gpL(T) for T in model_Ts]
    
    plt.plot(model_Ts-273.15, model_solubilities, color=predicted_plot_color, 
             label='predicted')
    
    plt.scatter(experimental_Ts-273.15, experimental_solubilities, marker=experimental_plot_marker, color=experimental_plot_color,
                label='experimental')
    
    plt.text(5, 180, f'R\u00b2 = {round(R_squared, 3)}')
    
    units_opening_brackets = ["[" , "[" ]
    for i in range(len(units_opening_brackets)):
        if units_on_newline[i]:
            units_opening_brackets[i] = "\n[" 
            
    plt.xlabel(xlabel +" " + units_opening_brackets[0] + temperature_units + "]")
    plt.ylabel(ylabel +" " + units_opening_brackets[1] + solubility_units + "]")
    
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    plt.xticks(xticks)
    plt.legend(loc=legend_loc)
    
    if save_fig:
        plt.savefig(filename, dpi=dpi, bbox_inches=bbox_inches,
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
    
    plt.show()

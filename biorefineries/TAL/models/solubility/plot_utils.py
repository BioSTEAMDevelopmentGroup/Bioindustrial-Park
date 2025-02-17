#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

def plot_solubility_model(experimental_Ts, experimental_solubilities, 
                          get_TAL_solubility_in_water_gpL,
                          model_Ts=273.15 + np.linspace(0., 99., 50),
                          R_squared=None,
                          xlim=(0,100), ylim=(0,200),
                          xticks=[0 + 10*i for i in range(0,11)],
                          save_fig=False,
                          filename='solubility_plot.png',
                          # ylabel=r"$\bfTAL$" +" " +r"$\bfSolubility$" +" " + r"$\bfin$" +" " + r"$\bfWater$",
                          # xlabel=r"$\bfTemperature$",
                          ylabel="TAL Solubility in Water",
                          xlabel="Temperature",
                          solubility_units= r"$\mathrm{g-TAL}\cdot\mathrm{L}^{-1}$",
                          temperature_units="\u00b0" + r"$\mathrm{C}$",
                          units_on_newline=(False, True), # (x, y)
                          legend_loc='center left',
                          fontname="Arial Unicode",
                          fontsize="14",
                          predicted_plot_color='black',
                          experimental_plot_color='#63C6CE', 
                          experimental_plot_marker='x',
                          dpi=600,
                          bbox_inches='tight',
                          xlabel_fontsize=16,
                          ylabel_fontsize=16,
                          ):


    plt.rcParams['font.sans-serif'] = fontname
    plt.rcParams['font.size'] = fontsize
    
    model_solubilities = [get_TAL_solubility_in_water_gpL(T) for T in model_Ts]
    
    fig, ax = plt.subplots()
    
    plt.plot(model_Ts-273.15, model_solubilities, color=predicted_plot_color, 
             label='predicted', zorder=1)
    
    plt.scatter(experimental_Ts-273.15, experimental_solubilities, marker=experimental_plot_marker, color=experimental_plot_color,
                label='experimental', zorder=2)
    
    plt.text(5, 180, f'R\u00b2 = {round(R_squared, 3)}')
    
    units_opening_brackets = ["[" , "[" ]
    for i in range(len(units_opening_brackets)):
        if units_on_newline[i]:
            units_opening_brackets[i] = "\n[" 
    
    plt.xlabel(xlabel +" " + units_opening_brackets[0] + temperature_units + "]", fontsize=xlabel_fontsize)
    plt.ylabel(ylabel +" " + units_opening_brackets[1] + solubility_units + "]",  fontsize=ylabel_fontsize)
    
    
    plt.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        direction='out',
        # right=True,
        width=0.65,
        labelsize=fontsize,
        # zorder=200,
        )
    
    plt.tick_params(
        axis='y',          
        which='major',      
        length=7,
        )

    plt.tick_params(
        axis='y',          
        which='minor',      
        length=3.5,
        )
    
    plt.tick_params(
        axis='x',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        direction='out',
        # right=True,
        width=0.65,
        labelsize=fontsize,
        # zorder=200,
        )

    plt.tick_params(
        axis='x',          
        which='major',      
        length=7,
        )

    plt.tick_params(
        axis='x',          
        which='minor',      
        length=3.5,
        )
    

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
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

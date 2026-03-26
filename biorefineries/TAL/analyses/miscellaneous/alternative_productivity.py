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

from biorefineries import TAL
from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import TAL_tea, TAL_lca, R302, spec, TAL_product, simulate_and_print, get_TAL_MPSP, theoretical_max_g_TAL_per_g_SA

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os

#%%
from matplotlib import pyplot as plt
def plot_multiple_metrics(x_axis_list, y_axis_list_of_lists, 
                          xlim=(0,100), 
                          ylims_list=[],
                          xticks=[0 + 10*i for i in range(0,11)],
                          yticks_list=[],
                          plot_type='line',
                          save_fig=True,
                          filename='alternative_productivities_plot.png',
                          ylabels_list=[r"$\bfMPSP$",
                                        r"$\bfFermentation bioreactors installed cost$".replace(' ', '\ '),
                                        r"$\bfTotal installed cost excl. fermentation bioreactors$".replace(' ', '\ '),
                                        r"$\bfFermentation bioreactors electricity demand$".replace(' ', '\ '),
                                        r"$\bfTotal positive electricity demand excl. fermentation bioreactors$".replace(' ', '\ '),
                                        ],
                          xlabel=r"$\bfFermentation$"+" " +r"$\bfProductivity$",
                          y_axis_units_list= [
                                              r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                                              "MM\$", 
                                              "MM\$", 
                                              "MW", 
                                              "MW", 
                                              r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                                              r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$",
                                              ],
                          x_axis_units= r"$\mathrm{g}\cdot\mathrm{L}^{-1}\cdot\mathrm{h}^{-1}$",
                          units_on_newline=False, # (x, y)
                          # legend_loc='center left',
                          fontname="Arial Unicode",
                          fontsize="12",
                          y_plotline_colors_list=[
                                                  "black",
                                                  "#A97802", 
                                                  '#607429',
                                                  '#A100A1',
                                                  '#f98f60',
                                                  ],
                          # experimental_plot_color='#63C6CE', 
                          # experimental_plot_marker='x',
                          dpi=600,
                          bbox_inches='tight',
                          ):

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=1)
    # twinaxes = [ax.twinx() for i in range(1, len(y_axis_list_of_lists))]
    plt.rcParams['font.sans-serif'] = fontname
    plt.rcParams['font.size'] = fontsize
    plt.xlim(xlim[0], xlim[1])
    
    p1 = None
    if plot_type=='line':
        p1, = ax.plot(x_axis_list, y_axis_list_of_lists[0],
                 color=y_plotline_colors_list[0],
             # label='predicted', 
             zorder=1)
    elif plot_type=='scatter':
        p1 = ax.scatter(x_axis_list, y_axis_list_of_lists[0],
                 color=y_plotline_colors_list[0],
             # label='predicted', 
             zorder=1)
    else:
        raise ValueError(f'Unknown plot_type {plot_type}.')
        
    if ylims_list:
        ax.set_ylim(ylims_list[0])
    ax.yaxis.label.set_color(y_plotline_colors_list[0])
    units_opening_brackets = ["[" , "[" ]
    for j in range(len(units_opening_brackets)):
        if units_on_newline:
            units_opening_brackets[j] = "\n[" 
            
    ax.set_ylabel(ylabels_list[0] +" " + units_opening_brackets[1] + y_axis_units_list[0] + "]")
    ax.tick_params(axis='y', colors=y_plotline_colors_list[0], direction='inout')
    ax.spines['right'].set_color(y_plotline_colors_list[0])
    
    ax.set_yticks(yticks_list[0])
    
    ax.tick_params(axis='x', direction='inout')
    
    for i in range(1, len(y_axis_list_of_lists)):
        twini = ax.twinx()
        
        pi = None
        if plot_type=='line':
            pi, = twini.plot(x_axis_list, y_axis_list_of_lists[i],
                         color=y_plotline_colors_list[i], 
                     # label='predicted', 
                     zorder=i+1)
        elif plot_type=='scatter':
            pi = twini.scatter(x_axis_list, y_axis_list_of_lists[i],
                         color=y_plotline_colors_list[i], 
                     # label='predicted', 
                     zorder=i+1)
        else:
            raise ValueError(f'Unknown plot_type {plot_type}.')
            
        twini.yaxis.label.set_color(y_plotline_colors_list[i])
        twini.spines['right'].set_color(y_plotline_colors_list[i])
        
        twini.spines.right.set_position(("axes", 1 + 0.2*(i-1)))
        twini.tick_params(axis='y', colors=y_plotline_colors_list[i], direction='inout')
        
        units_opening_brackets = ["[" , "[" ]
        for j in range(len(units_opening_brackets)):
            if units_on_newline:
                units_opening_brackets[j] = "\n[" 
                
        twini.set_ylabel(ylabels_list[i] +" " + units_opening_brackets[1] + y_axis_units_list[i] + "]")
        
        if ylims_list:
            twini.set_ylim(ylims_list[i][0], ylims_list[i][1])
        
        twini.set_yticks(yticks_list[i])
        
    ax.set_xlabel(xlabel +" " + units_opening_brackets[0] + x_axis_units + "]")
    ax.set_xticks(xticks)
    
    if save_fig:
        plt.savefig(filename, dpi=dpi, bbox_inches=bbox_inches,
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
    
    plt.show()
    
    


#%% 
simulate_and_plot = True

if simulate_and_plot:
    
    #%% Filepaths from biorefineries.TAL.models import models_TAL_solubility_exploit as models

    from biorefineries.TAL.models import models_TAL_solubility_exploit as models
    chdir = os.chdir

    dateTimeObj = datetime.now()

    ig = np.seterr(invalid='ignore')
    # bst.speed_up()

    product = TAL_product

    # search page for high end: https://www.alibaba.com/trade/search?spm=a2700.galleryofferlist.0.0.2a995827YzqZVg&fsb=y&IndexArea=product_en&assessmentCompany=true&keywords=590-00-1+sorbate&productTag=1200000228&ta=y&tab=all&
    SA_market_range=np.array([
                              6.74, # 2019 global high end from Sorbic Acid Market, Transparency Market Research
                              6.50 * 1.3397087, # $6.50/kg-potassium-sorbate from https://www.alibaba.com/product-detail/Lifecare-Supply-Potassium-Sorbate-High-Quality_1600897125355.html?spm=a2700.galleryofferlist.p_offer.d_title.1bc15827eAs1TL&s=p
                              ]) 

    TAL_maximum_viable_market_range = SA_market_range / theoretical_max_g_TAL_per_g_SA

    # TAL_maximum_viable_market_range = np.array([5.99, 7.74])

    TAL_filepath = TAL.__file__.replace('\\__init__.py', '')

    # ## Change working directory to biorefineries\\TAL\\analyses\\results
    # chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'

    model = models.TAL_model
    system = TAL_sys = models.TAL_sys
    
    system.flowsheet.M401.base_neutralizes_acids = False
    
    modes = [
                'A_FGI_sugarcane',
             ]
    
    
    parameter_distributions_filenames = [
                                        'parameter-distributions_TAL_' + mode +'.xlsx' 
                                        for mode in modes
                                        ]
    mode = modes[0]
    
    
    parameter_distributions_filename = TAL_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+parameter_distributions_filenames[0]
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)
    
    # load_additional_params()
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=2000, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    # ## Change working directory to biorefineries\\TAL\\analyses\\results
    # chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    #%% Baseline
    model.exception_hook = 'warn'
    print('\n\nSimulating baseline ...')
    baseline_initial = model.metrics_at_baseline()
    
    #%%
    simulate_and_print()
    
    #%%  Metrics
    broth = R302.outs[1]
    # SA_price_range = [6500, 7500]
    
    product_chemical_IDs = ['TAL',]
    get_product_MPSP = lambda: TAL_tea.solve_price(product) / get_product_purity() # USD / pure-kg
    get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
    get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])
    
    get_product_recovery = lambda: sum([product.imass[i] for i in product_chemical_IDs])/sum([broth.imass[i] for i in product_chemical_IDs])
    get_TAL_AOC = lambda: TAL_tea.AOC / 1e6 # million USD / y
    get_TAL_TCI = lambda: TAL_tea.TCI / 1e6 # million USD
    
    get_TAL_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol
    
    get_TAL_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol
    
    TAL_metrics = [get_product_MPSP, 
                    lambda: TAL_lca.GWP,
                    # lambda: TAL_lca.GWP - TAL_lca.net_electricity_GWP, 
                    lambda: TAL_lca.FEC, 
                    # lambda: TAL_lca.FEC - TAL_lca.net_electricity_FEC,
                   get_TAL_AOC, get_TAL_TCI, 
                   get_product_purity]
    
    #%%
    import biosteam as bst
    s,u = bst.main_flowsheet.stream, bst.main_flowsheet.unit
    R303 = u.R303
    
    #%% Run analysis
    steps = 50
    productivities = np.linspace(0.024, 0.600, steps)
    MPSPs, ferm_ICs, sys_ICs, ferm_pows, sys_pows = [], [], [], [], []
    recoveries = []
    
    for p in productivities:
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2,
                                 spec_3=p)
        MPSPs.append(get_TAL_MPSP())
        ferm_ICs.append(R302.installed_cost + R303.installed_cost)
        sys_ICs.append(TAL_tea.installed_equipment_cost)
        ferm_pows.append(R302.power_utility.rate + R303.power_utility.consumption)
        sys_pows.append(TAL_sys.power_utility.consumption)
    
    MPSPs, ferm_ICs, sys_ICs, ferm_pows, sys_pows =\
        np.array(MPSPs), np.array(ferm_ICs), np.array(sys_ICs), np.array(ferm_pows), np.array(sys_pows)
    
#%%

    plot_multiple_metrics(productivities,
                          [
                           MPSPs,
                           ferm_ICs/1e6, 
                           (sys_ICs - ferm_ICs)/1e6,
                           ferm_pows/1e3,
                           (sys_pows - ferm_pows)/1e3,
                           ],
                          
                          xlim=(0,0.6), 
                          ylims_list=[
                                      (0,8),
                                      (0,200),
                                      (0,200),
                                      (0,12),
                                      (0,12),
                                      ],
                          xticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                          yticks_list = [
                                     np.arange(0, 8.1, 1.0),
                                     np.arange(0, 200.1, 50),
                                     np.arange(0, 200.1, 50),
                                     np.arange(0, 12.1, 2),
                                     np.arange(0, 12.1, 2),
                                     ]
                          )
    



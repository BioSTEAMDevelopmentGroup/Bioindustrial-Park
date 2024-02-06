# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: sarangbhagwat
"""

# Run this cell first
from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors
import numpy as np

from biorefineries import TAL
from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import TAL_tea, TAL_lca, R302, spec, TAL_product, simulate_and_print, theoretical_max_g_TAL_per_g_SA

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os

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


#%% Filepaths
TAL_filepath = TAL.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\TAL\\analyses\\results
# chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'

#%% Load baseline
model = models.TAL_model
system = TAL_sys = models.TAL_sys

modes = [
            # 'A',
            'B',
         ]

parameter_distributions_filenames = [
                                    # 'parameter-distributions_A.xlsx',
                                    'parameter-distributions_B.xlsx',
                                    ]
mode = modes[0]


parameter_distributions_filename = TAL_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+parameter_distributions_filenames[0]
print(f'\n\nLoading parameter distributions ({mode}) ...')
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename)

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

#%%
splits = np.linspace(0., 0.95, 20)
MPSPs, GWPs, FECs = [], [], []
recoveries = []
for spl in splits:
    u.S403.split = spl
    simulate_and_print()
    MPSPs.append(get_product_MPSP())
    GWPs.append(TAL_lca.GWP)
    FECs.append(TAL_lca.FEC)
    recoveries.append(get_product_recovery())

#%%
from matplotlib import pyplot as plt
def plot_multiple_metrics(x_axis_list, y_axis_list_of_lists, 
                          xlim=(0,100), 
                          ylims_list=[],
                          xticks=[0 + 10*i for i in range(0,11)],
                          yticks_list=[],
                          save_fig=True,
                          filename='supernatant_recycling_fixed_decarboxylation_plot.png',
                          ylabels_list=[r"$\bfTAL$"+" " +r"$\bfRecovery$",
                                        r"$\bfMPSP$",
                                        r"$\mathrm{\bfGWP}_{\bf100}$",
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
                          # legend_loc='center left',
                          fontname="Arial Unicode",
                          fontsize="12",
                          y_plotline_colors_list=[
                                                  "black",
                                                  "#A97802", 
                                                  '#607429',
                                                  '#A100A1',
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
    p1, = ax.plot(x_axis_list, y_axis_list_of_lists[0],
             color=y_plotline_colors_list[0],
         # label='predicted', 
         zorder=1)
    if ylims_list:
        ax.set_ylim(ylims_list[0])
    ax.yaxis.label.set_color(p1.get_color())
    units_opening_brackets = ["[" , "[" ]
    for j in range(len(units_opening_brackets)):
        if units_on_newline:
            units_opening_brackets[j] = "\n[" 
            
    ax.set_ylabel(ylabels_list[0] +" " + units_opening_brackets[1] + y_axis_units_list[0] + "]")
    ax.tick_params(axis='y', colors=p1.get_color(), direction='inout')
    ax.spines['right'].set_color(p1.get_color())
    
    ax.set_yticks(yticks_list[0])
    
    ax.tick_params(axis='x', direction='inout')
    
    for i in range(1, len(y_axis_list_of_lists)):
        twini = ax.twinx()
        pi, = twini.plot(x_axis_list, y_axis_list_of_lists[i],
                     color=y_plotline_colors_list[i], 
                 # label='predicted', 
                 zorder=i+1)
        twini.yaxis.label.set_color(pi.get_color())
        twini.spines['right'].set_color(pi.get_color())
        
        twini.spines.right.set_position(("axes", 1 + 0.2*(i-1)))
        twini.tick_params(axis='y', colors=pi.get_color(), direction='inout')
        
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
plot_multiple_metrics(100*np.array(splits),
                      [
                       100*np.array(recoveries), 
                       MPSPs,
                       GWPs, 
                       FECs,
                       ],
                      
                      xlim=(0,100), 
                      ylims_list=[
                                  (0,100),
                                  (5,10),
                                  (6,16),
                                  (-30,60),
                                  ],
                      xticks=[0 + 10*i for i in range(0,11)],
                      yticks_list = [[0 + 10*i for i in range(0,11)],
                                 np.arange(5, 10.1, 1),
                                 np.arange(6, 16.1, 2),
                                 np.arange(-30, 60.1, 10),
                                 ]
                      )

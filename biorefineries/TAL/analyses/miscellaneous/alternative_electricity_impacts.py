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
from biorefineries.TAL.systems.system_TAL_solubility_exploit_ethanol_sugarcane import TAL_tea, TAL_lca, R302, spec, TAL_product, simulate_and_print, theoretical_max_g_TAL_per_g_SA, fermentation_group,\
    F301, F301_P, M304

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os

# from biorefineries.TAL.models import models_TAL_solubility_exploit_sugarcane as models
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
            # 'A_FGI_sugarcane',
            'A',
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

get_product_recovery_FGI = lambda: product.imass[product_chemical_IDs].sum()/broth.imass['TAL']

def get_F301_heat_utility_duty():
    if F301.heat_utilities:
        return abs(F301.heat_utilities[0].duty) + abs(F301.heat_utilities[1].duty)
    else:
        return 0

get_sugar_conc_utility_demand_proxy = lambda: get_F301_heat_utility_duty() +\
                                              3600*(F301_P.power_utility.consumption + M304.power_utility.consumption)

get_sugar_conc_TCI = lambda: F301.installed_cost + F301_P.installed_cost + M304.installed_cost

#%%
from biorefineries.TAL.process_settings import CFs
from matplotlib import pyplot as plt
    
#%%#%%

baseline_CFs = CFs['GWP_100']['Electricity'], CFs['FEC']['Electricity']

steps = 100

electricty_GWP_CFs = np.linspace(0., 1, steps)
electricty_FEC_CFs = np.linspace(0., 10, steps)

GWPs = []
for gcf in electricty_GWP_CFs:
    CFs['GWP_100']['Electricity'] = gcf
    GWPs.append(TAL_lca.GWP)
GWPs = np.array(GWPs)

FECs = []
for fcf in electricty_FEC_CFs:
    CFs['FEC']['Electricity'] = fcf
    FECs.append(TAL_lca.FEC)
FECs = np.array(FECs)

#%%
def plot_metric(x_axis_list, y_axis_list, 
                          xlim=(0,100), 
                          ylim=[],
                          xticks=[0 + 10*i for i in range(0,11)],
                          yticks=[],
                          plot_type='line',
                          save_fig=True,
                          filename='alternative_electricity_impacts.png',
                          # ylabel=r"$\bfMPSP$",
                                        # r"$\mathrm{\bfGWP}_{\bf100}$",
                          ylabel=r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$",
                          # ylabel=r"$\bfFEC$" ,
                          xlabel=r"$\bfSupernatant$"+" " +r"$\bfRecycling$",
                          # y_axis_units=r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                          y_axis_units=r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$",
                          # y_axis_units=r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$",
                          x_axis_units=r"$\mathrm{\%}$",
                          units_on_newline=False, # (x, y)
                          # legend_loc='center left',
                          fontname="Arial Unicode",
                          fontsize="12",
                          # linecolor="#A97802", 
                          linecolor='#607429',
                          # linecolor='#A100A1',
                          
                          # experimental_plot_color='#63C6CE', 
                          # experimental_plot_marker='x',
                          dpi=600,
                          bbox_inches='tight',
                          ):

    fig, ax = plt.subplots()
    
    # twinaxes = [ax.twinx() for i in range(1, len(y_axis_list))]
    plt.rcParams['font.sans-serif'] = fontname
    plt.rcParams['font.size'] = fontsize
    plt.xlim(xlim[0], xlim[1])
    
    p1 = None
    if plot_type=='line':
        p1, = ax.plot(x_axis_list, y_axis_list,
                 color=linecolor,
             # label='predicted', 
             zorder=1)
    elif plot_type=='scatter':
        p1 = ax.scatter(x_axis_list, y_axis_list,
                 color=linecolor,
             # label='predicted', 
             zorder=1)
    else:
        raise ValueError(f'Unknown plot_type {plot_type}.')
        
    if ylim:
        ax.set_ylim(ylim)
    # ax.yaxis.label.set_color(linecolor)
    units_opening_brackets = ["[" , "[" ]
    for j in range(len(units_opening_brackets)):
        if units_on_newline:
            units_opening_brackets[j] = "\n[" 
            
    ax.set_ylabel(ylabel +" " + units_opening_brackets[1] + y_axis_units + "]")
    ax.tick_params(axis='y', direction='inout')
    # ax.spines['right'].set_color(linecolor)
    
    ax.set_yticks(yticks)
    
    ax.tick_params(axis='x', direction='inout')
    
    
    ax.set_xlabel(xlabel +" " + units_opening_brackets[0] + x_axis_units + "]")
    ax.set_xticks(xticks)
    
    if save_fig:
        plt.savefig(filename, dpi=dpi, bbox_inches=bbox_inches,
                    # facecolor=plt.get_facecolor(),
                    # transparent=False,
                    )
    # breakpoint()
    plt.show()


#%%

plot_metric(electricty_GWP_CFs,
            GWPs,
                      
                      xlim=(0,1), 
                      ylim=(-10, 15),
                      xticks=[np.arange(0, 1.01, 0.1)],
                      yticks = np.arange(-10, 15.01, 5),
                      save_fig=False
                      )

#%%

plot_multiple_metrics(electricty_FEC_CFs,
                      [
                       FECs, 
                       0*electricty_GWP_CFs,
                       ],
                      
                      xlim=(0,1), 
                      ylim=[
                                  (-10, 15),
                                  (0, 0)
                                  ],
                      xticks=[np.arange(0, 1.01, 0.1)],
                      yticks = [np.arange(-10, 15.01, 5),
                                 ]
                      )

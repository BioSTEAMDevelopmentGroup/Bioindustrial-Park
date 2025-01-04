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
from biorefineries.TAL.systems.system_SA_solubility_exploit_ethanol_corn import TAL_tea, TAL_lca, R302, spec, KSA_product, simulate_and_print, theoretical_max_g_TAL_per_g_SA,\
    F301, F301_P, M304

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os

from biorefineries.TAL.models import models_SA_solubility_exploit_corn as models

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')
# bst.speed_up()

product = KSA_product

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
            'A_FGI_corn',
         ]


parameter_distributions_filenames = [
                                    'parameter-distributions_Sorbate_' + mode +'.xlsx' 
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
# simulate_and_print()

#%%  Metrics
broth = R302.outs[1]
# SA_price_range = [6500, 7500]

product_chemical_IDs = ['KSA',]
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


# %% Generate 3-specification meshgrid and set specification loading functions

steps = (50, 50, 1)

# Yield, titer, productivity (rate)
spec_1 = yields = np.linspace(0.05, 0.5724, steps[0]) # yield
spec_2 = titers = np.linspace(10., 
                              100., # although sugar concentration limit of 600 g/L would allow as high as 230 g-TAL/L, we set an upper limit of 100 g/L
                                   # based on achieved (50-68 g/L using E.coli, Candida) and targeted (100 g/L) titers for adipic acid, another organic solid with low water solubility
                                   # Skoog et al., 2018 ( https://doi.org/10.1016/j.biotechadv.2018.10.012 )
                                steps[1]) # titer

# spec_3 = productivities =\
#     np.linspace(0.05, 1.5, steps[2])

which_fig = 'insights'

#%%
additional_points = {}
spec_3 = productivities = None
TAL_metrics = None

if which_fig is None: # if you want to run at your own productivities and use your own metrics
    spec_3 = productivities =\
        np.array([
                    # 0.2*spec.baseline_productivity,
                    1.*spec.baseline_productivity,
                    # 5.*spec.baseline_productivity,
                  ])
    TAL_metrics = [get_product_MPSP, 
                    lambda: TAL_lca.GWP,
                    lambda: TAL_lca.FEC, 
                   get_TAL_AOC, get_TAL_TCI, 
                   get_product_purity]
    additional_points = {(40.5, 35.9):('D', 'w', 6)}
else:
    if which_fig=='5':
        spec_3 = productivities =\
            np.array([
                        1.*spec.baseline_productivity,
                      ])
        TAL_metrics = [get_product_MPSP, 
                        lambda: TAL_lca.GWP,
                        lambda: TAL_lca.FEC, 
                       get_TAL_AOC, get_TAL_TCI, 
                       get_product_purity]
        additional_points = {(40.5, 35.9):('D', 'w', 6)}
    elif which_fig=='S4':
        spec_3 = productivities =\
            np.array([
                        0.2*spec.baseline_productivity,
                      ])
        TAL_metrics = [get_product_MPSP, 
                        lambda: TAL_lca.GWP,
                        lambda: TAL_lca.FEC, 
                       get_TAL_AOC, get_TAL_TCI, 
                       get_product_purity]
    elif which_fig=='S5':
        spec_3 = productivities =\
            np.array([
                        5.*spec.baseline_productivity,
                      ])
        TAL_metrics = [get_product_MPSP, 
                        lambda: TAL_lca.GWP,
                        lambda: TAL_lca.FEC, 
                       get_TAL_AOC, get_TAL_TCI, 
                       get_product_purity]
    elif which_fig=='S11':
        spec_3 = productivities =\
            np.array([
                        1.*spec.baseline_productivity,
                      ])
        TAL_metrics = [get_product_MPSP, 
                        # lambda: TAL_lca.GWP,
                        lambda: TAL_lca.GWP - TAL_lca.net_electricity_GWP, 
                        # lambda: TAL_lca.FEC, 
                        lambda: TAL_lca.FEC - TAL_lca.net_electricity_FEC,
                       get_TAL_AOC, get_TAL_TCI, 
                       get_product_purity]
        additional_points = {(40.5, 35.9):('D', 'w', 6)}
    elif which_fig=='insights':
        spec_3 = productivities =\
            np.array([
                        5.*spec.baseline_productivity,
                      ])
        TAL_metrics = [get_product_MPSP, 
                        lambda: TAL_lca.GWP,
                        lambda: TAL_lca.FEC, 
                       # get_sugar_conc_utility_demand_proxy,
                       # get_sugar_conc_TCI,
                       get_product_recovery_FGI,
                       get_TAL_AOC, get_TAL_TCI,]
    
#%% Plot stuff

# Parameters analyzed across

x_label = r"$\bfYield$" # title of the x axis
x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
x_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

y_label = r"$\bfTiter$" # title of the y axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


z_label = r"$\bfProductivity$" # title of the z axis
z_units =  r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$"
z_ticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

# Metrics
MPSP_w_label = r"$\bfMPSP$" # title of the color axis
MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"

# GWP_w_label = r"$\mathrm{\bfGWP}_{\bf100}$"
GWP_w_label = r"$\mathrm{\bfCarbon}$" + " " + r"$\mathrm{\bfIntensity}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"

FEC_w_label = r"$\bfFEC$" # title of the color axis
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"

AOC_w_label = r"$\bfAOC$" # title of the color axis
AOC_units = r"$\mathrm{MM\$}\cdot\mathrm{y}^{-1}$"

TCI_w_label = r"$\bfTCI$" # title of the color axis
TCI_units = r"$\mathrm{MM\$}$"

Purity_w_label = r"$\bfPurity$" # title of the color axis
Purity_units = "wt " + r"$\mathrm{\%}$"

#%% Colors

marketrange_shadecolor = (*colors.neutral.shade(50).RGBn, 0.3)

oversaccharine_shadecolor_raw = colors.CABBI_teal_green.tint(40)
# oversaccharine_shadecolor_raw = colors.CABBI_green.shade(45)
inhibited_shadecolor_raw = colors.CABBI_grey.shade(60)
oversaccharine_shadecolor = (*oversaccharine_shadecolor_raw.RGBn, 1)
inhibited_shadecolor = (*inhibited_shadecolor_raw.RGBn, 1)
# overlap_color = (*(colors.CABBI_teal_green.tint(20).RGBn + colors.CABBI_black.tint(20).RGBn)/2, 1)
overlap_color = (*(oversaccharine_shadecolor_raw.RGBn + inhibited_shadecolor_raw.RGBn)/2, 1)
linecolor_dark = (*colors.CABBI_black.shade(40).RGBn, 0.95)
linecolor_light = (*colors.neutral_tint.RGBn, 0.85)
markercolor = (*colors.CABBI_orange.shade(5).RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)


def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    CABBI_colors = (colors.CABBI_orange.RGBn,
                    colors.CABBI_yellow.RGBn,

                    colors.CABBI_green.RGBn,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)


def CABBI_grey_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's grey colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_grey.RGBn,
                    colors.grey_shade.RGBn,
                    colors.grey_dark.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

def CABBI_blue_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's blue colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_blue_light.RGBn,
                    colors.CABBI_teal.RGBn,
                    colors.CABBI_teal_green.shade(25).RGBn,
                    colors.CABBI_green_dirty.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

#%% Tickmark utils (unused)

def tickmarks_from_data(data, accuracy=50, N_points=5):
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    dmin = floor(dmin/accuracy) * accuracy
    dmax = ceil(dmax/accuracy) * accuracy
    step = (dmax - dmin) / (N_points - 1)
    return [dmin + step * i for i in range(N_points)]

#%%
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = None
if which_fig is None:
    file_to_save = f'_{steps}_steps_'+'TAL_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
else:
    file_to_save = 'Fig_'+which_fig+'_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
  
#%% Create meshgrid
spec_1, spec_2 = np.meshgrid(spec_1, spec_2)

#%% Initial simulation
# simulate_and_print()

print('\n\nSimulating the initial point to avoid bugs ...')
# spec.byproduct_yields_decrease_policy = 'simultaneous, from 0 product yield'
spec.load_specifications(yields[0], titers[0], productivities[0])
# spec.set_production_capacity(desired_annual_production=spec.desired_annual_production)
# simulate_and_print()
for i in range(3): TAL_sys.simulate()
print(get_product_MPSP())

# %% Run TRY analysis 
system = TAL_sys
HXN = spec.HXN
results_metric_1, results_metric_2, results_metric_3 = [], [], []
results_metric_4, results_metric_5, results_metric_6 = [], [], []

def print_status(curr_no, total_no, s1, s2, s3, HXN_qbal_error, results=None, exception_str=None,):
    print('\n\n')
    print(f'{curr_no}/{total_no}')
    print('\n')
    print(s1, s2, s3)
    print('\n')
    print(f'HXN Qbal error = {round(HXN_qbal_error, 2)} %.')
    print('\n')
    print(results)
    print('\nError: ', exception_str)

max_HXN_qbal_percent_error = 0.

curr_no = 0
total_no = len(yields)*len(titers)*len(productivities)

print_status_every_n_simulations = 50

titers_mol_per_mol_total = [] # for fermentation generalizable insights work only

for p in productivities:
    # data_1 = TAL_data = spec.evaluate_across_specs(
    #         TAL_sys, spec_1, spec_2, TAL_metrics, [p])
    
    d1_Metric1, d1_Metric2, d1_Metric3 = [], [], []
    d1_Metric4, d1_Metric5, d1_Metric6 = [], [], []
    
    for y in yields:
        d1_Metric1.append([])
        d1_Metric2.append([])
        d1_Metric3.append([])
        d1_Metric4.append([])
        d1_Metric5.append([])
        d1_Metric6.append([])
        titers_mol_per_mol_total.append([])
        for t in titers:
            curr_no +=1
            error_message = None
            try:
                spec.load_specifications(spec_1=y, spec_2=t, spec_3=p)
                
                # spec.set_production_capacity(desired_annual_production=spec.desired_annual_production)
                
                for i in range(3): system.simulate()
                titers_mol_per_mol_total[-1].append(broth.imol['TAL']/broth.imol['TAL', 'Water'].sum())
                
                d1_Metric1[-1].append(TAL_metrics[0]())
                d1_Metric2[-1].append(TAL_metrics[1]())
                d1_Metric3[-1].append(TAL_metrics[2]())
                d1_Metric4[-1].append(TAL_metrics[3]())
                d1_Metric5[-1].append(TAL_metrics[4]())
                d1_Metric6[-1].append(TAL_metrics[5]())
                
                HXN_qbal_error = HXN.energy_balance_percent_error
                if abs(max_HXN_qbal_percent_error)<abs(HXN_qbal_error): max_HXN_qbal_percent_error = HXN_qbal_error
            
            except RuntimeError as e1:
                d1_Metric1[-1].append(np.nan)
                d1_Metric2[-1].append(np.nan)
                d1_Metric3[-1].append(np.nan)
                d1_Metric4[-1].append(np.nan)
                d1_Metric5[-1].append(np.nan)
                d1_Metric6[-1].append(np.nan)
                error_message = str(e1)
            
            except ValueError as e1:
                d1_Metric1[-1].append(np.nan)
                d1_Metric2[-1].append(np.nan)
                d1_Metric3[-1].append(np.nan)
                d1_Metric4[-1].append(np.nan)
                d1_Metric5[-1].append(np.nan)
                d1_Metric6[-1].append(np.nan)
                error_message = str(e1)
            
            if curr_no%print_status_every_n_simulations==0 or (error_message and not 'sugar' in error_message and not 'opposite sign' in error_message):
                print_status(curr_no, total_no,
                             y, t, p, 
                             results=[d1_Metric1[-1][-1], d1_Metric2[-1][-1], d1_Metric3[-1][-1],
                                      d1_Metric4[-1][-1], d1_Metric5[-1][-1], d1_Metric6[-1][-1],],
                             HXN_qbal_error=HXN.energy_balance_percent_error,
                             exception_str=error_message)
    
    d1_Metric1, d1_Metric2, d1_Metric3 = np.array(d1_Metric1), np.array(d1_Metric2), np.array(d1_Metric3)
    d1_Metric4, d1_Metric5, d1_Metric6 = np.array(d1_Metric4), np.array(d1_Metric5), np.array(d1_Metric6)
    
    d1_Metric1, d1_Metric2, d1_Metric3 = d1_Metric1.transpose(), d1_Metric2.transpose(), d1_Metric3.transpose()
    d1_Metric4, d1_Metric5, d1_Metric6 = d1_Metric4.transpose(), d1_Metric5.transpose(), d1_Metric6.transpose()
    
    results_metric_1.append(d1_Metric1)
    results_metric_2.append(d1_Metric2)
    results_metric_3.append(d1_Metric3)
    results_metric_4.append(d1_Metric4)
    results_metric_5.append(d1_Metric5)
    results_metric_6.append(d1_Metric6)


    # %% Save generated data
    
    csv_file_to_save = file_to_save + f'_prod_{p}'
    pd.DataFrame(d1_Metric1).to_csv(TAL_results_filepath+'MPSP-'+csv_file_to_save+'.csv')
    pd.DataFrame(d1_Metric2).to_csv(TAL_results_filepath+'GWP-'+csv_file_to_save+'.csv')
    pd.DataFrame(d1_Metric3).to_csv(TAL_results_filepath+'FEC-'+csv_file_to_save+'.csv')
    pd.DataFrame(d1_Metric4).to_csv(TAL_results_filepath+'AOC-'+csv_file_to_save+'.csv')
    pd.DataFrame(d1_Metric5).to_csv(TAL_results_filepath+'TCI-'+csv_file_to_save+'.csv')
    pd.DataFrame(d1_Metric6).to_csv(TAL_results_filepath+'Purity-'+csv_file_to_save+'.csv')
    

#%% Report maximum HXN energy balance error
print(f'Max HXN Q bal error was {round(max_HXN_qbal_percent_error, 3)} %.')

#%% Plot metrics vs titer, yield, and productivity

chdir(TAL_results_filepath)

results_metric_1 = np.array(results_metric_1)
results_metric_2 = np.array(results_metric_2)
results_metric_3 = np.array(results_metric_3)

#%% Save generated numpy file
np.save(TAL_results_filepath+'MPSP-'+file_to_save, results_metric_1)
np.save(TAL_results_filepath+'GWP-'+file_to_save, results_metric_2)
np.save(TAL_results_filepath+'FEC-'+file_to_save, results_metric_3)


#%% More plot utils

from math import floor, log
get_median = np.median
np_round = np.round

def get_OOM(number):
    return floor(log(number+1., 10))

def my_round(x, base=5):
    return base * round(x/base)

def get_contour_info_from_metric_data(
                                        metric_data, # numpy array
                                        round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM=1,
                                        n_levels_between_cbar_ticks=5,
                                        n_stdevs_for_bounds=1,
                                        lb=None,
                                        ub=None,
                                        multiply_step_size_by=0.5,
                                        w_ticks_round_to_base_divisor=2,
                                        remove_w_ticks_if_greater_than_this_fraction_of_successor=0.8,
                                      ):
    if not type(metric_data) == np.ndarray: metric_data = np.array(metric_data)
    metric_data = metric_data[~np.isnan(metric_data)]
    median, stdev = get_median(metric_data), metric_data.std()
    bound_diff = n_stdevs_for_bounds*stdev
    ub_temp = max(abs(median-bound_diff), abs(median+bound_diff))
    # breakpoint()
    OOM = get_OOM(ub_temp)
    log_ub_temp = log(ub_temp, 10)
    
    round_to_decimal_place = None
    
    if log_ub_temp < 1.: 
        round_to_decimal_place = int(round(ub_temp/(10.*round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM), 0))
    else:
        round_to_decimal_place = -OOM -1 + round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM
    
    if lb==None: lb = np_round(median - bound_diff, round_to_decimal_place)
    if ub==None: ub = np_round(median + bound_diff, round_to_decimal_place)
    
    cbar_ticks_step_size = (10**OOM)/(round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM) * multiply_step_size_by
    w_levels_step_size = cbar_ticks_step_size/n_levels_between_cbar_ticks
    cbar_ticks = np.arange(lb, ub+cbar_ticks_step_size, cbar_ticks_step_size)
    w_levels = np.arange(lb, ub+w_levels_step_size, w_levels_step_size)
    
    w_ticks_round_to_base = 10**-round_to_decimal_place / 2
    w_ticks = [*set([lb, 
               my_round(median-0.6*stdev, w_ticks_round_to_base),
               my_round(median-0.4*stdev, w_ticks_round_to_base),  
               my_round(median-0.2*stdev, w_ticks_round_to_base), 
               np_round(median, round_to_decimal_place),
               my_round(median+0.2*stdev, ), 
               my_round(median+0.4*stdev, w_ticks_round_to_base),  
               my_round(median+0.6*stdev, w_ticks_round_to_base), 
               ub])]
    
    # w_ticks.sort()
    # w_ticks_to_remove = []
    # for i in range(len(w_ticks)-1):
    #     if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
    #         w_ticks_to_remove.append(w_ticks[i])
    # for j in w_ticks_to_remove:
    #     w_ticks.remove(j)
    #     w_ticks.append(j-w_levels_step_size)
        
    w_ticks.sort()
    w_ticks_to_remove = []
    for i in range(len(w_ticks)-1):
        if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
            w_ticks_to_remove.append(w_ticks[i])
    for j in w_ticks_to_remove:
        w_ticks.remove(j)
    
    w_ticks.sort()
    return w_levels, w_ticks, cbar_ticks

#%% More plot stuff

fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5
keep_frames = True

print('\nCreating and saving contour plots ...\n')

#%% Plots
plot = False

if plot: 
    
    #%% MPSP
    
    # MPSP_w_levels, MPSP_w_ticks, MPSP_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
    MPSP_w_levels = np.arange(0., 12.25, 0.25)
    MPSP_cbar_ticks = np.arange(0., 12.1, 1.)
    MPSP_w_ticks = [1., 1.5, 2., 2.5, 2.75, 3., 3.5, 4, 5, 10, 12]
    # MPSP_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_1, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=MPSP_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=MPSP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=MPSP_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=MPSP_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=MPSP_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='MPSP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    comparison_range=TAL_maximum_viable_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    # manual_clabels_regular = {
                                    #     MPSP_w_ticks[5]: (45,28),
                                    #     },
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    # comparison_range=[MPSP_w_levels[-2], MPSP_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',
                                    
                                    # manual_clabels_regular = {
                                    #     2.50: (90,90),
                                    #     2.75: (85,60),
                                    #     3.00: (70,60),
                                    #     3.50: (85,30),
                                    #     4.00: (45,55),
                                    #     5.00: (50,20),
                                    #     10.0: (20,25),
                                    #     12.0: (20,15),
                                        
                                    # #     MPSP_w_ticks[0]: (70,70),
                                    # #     MPSP_w_ticks[1]: (60,58),
                                    # #     MPSP_w_ticks[2]: (60,45),
                                    # #     MPSP_w_ticks[3]: (55,30),
                                    # #     MPSP_w_ticks[4]: (50,25),
                                    # #     MPSP_w_ticks[5]: (12,12),
                                    #     },
                                    # manual_clabels_comparison_range =\
                                    #     {TAL_maximum_viable_market_range[0]:(33,20), 
                                    #       TAL_maximum_viable_market_range[1]:(32,10)},
                                    
                                    
                                    # contourplot_facecolor = colors.grey_dark.shade(8).RGBn,
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    bottom_fill_bounds = ((0,0), 
                                                          (1,11.),
                                                          (99,11.)),
                                    # zoom_data_scale=5,
                                    text_boxes = {'>12.0': [(5,5), 'white']},
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )
    
    #%% GWP
    
    # GWP_w_levels, GWP_w_ticks, GWP_cbar_ticks = get_contour_info_from_metric_data(results_metric_2,)
    GWP_w_levels = np.arange(-2, 14.1, 0.5)
    GWP_cbar_ticks = np.arange(-2, 14.1, 2.)
    GWP_w_ticks = [0, 1, 2, 3, 4, 6, 8, 10, 14]
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_2, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., GWP
                                    x_data=100*yields, # x axis values
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=GWP_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=GWP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=GWP_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=GWP_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=GWP_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='GWP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    additional_points = additional_points,
                                    
                                    # manual_clabels_regular = {
                                    #     # MPSP_w_ticks[0]: (80,300),
                                    #     # MPSP_w_ticks[1]: (50,320),
                                    #     # MPSP_w_ticks[2]: (25,310),
                                    #     # MPSP_w_ticks[3]: (18,300),
                                    #     # MPSP_w_ticks[4]: (10,50),
                                        
                                    #     GWP_w_ticks[0]: (60,70),
                                    #     GWP_w_ticks[1]: (60,60),
                                    #     GWP_w_ticks[2]: (60,50),
                                    #     GWP_w_ticks[3]: (60,40),
                                    #     GWP_w_ticks[4]: (60,30),
                                    #     GWP_w_ticks[5]: (60,20),
                                    #     GWP_w_ticks[6]: (60,10),
                                    #     }
                                    
                                    # cbar_n_minor_ticks=1,
                                    # comparison_range=[6.5, 7.5],
                                    # comparison_range=[GWP_w_levels[-2], GWP_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',
                                    
                                    
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    bottom_fill_bounds = ((0,0), 
                                                          (1,11.),
                                                          (99,11.)),
                                    # zoom_data_scale=5,
                                    text_boxes = {'>14.0': [(80,5), 'white']},
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    
                                    )
    
    
    #%% FEC
    
    # FEC_w_levels, FEC_w_ticks, FEC_cbar_ticks = get_contour_info_from_metric_data(results_metric_3,)
    FEC_w_levels = np.arange(-100, 101, 10)
    FEC_cbar_ticks = np.arange(-100, 101, 20)
    FEC_w_ticks = [-100, -60, -30, 0, 30, 60, 100]
    # FEC_w_ticks = [40, 50, 70, 80, 100]
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_3, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., FEC
                                    x_data=100*yields, # x axis values
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=FEC_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=FEC_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=FEC_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=FEC_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(200), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    cmap_under_color = colors.CABBI_orange.shade(1).RGBn,
                                    extend_cmap='both',
                                    cbar_ticks=FEC_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='FEC_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 1,
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    # comparison_range=[6.5, 7.5],
                                    # comparison_range=[FEC_w_levels[-2], FEC_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',
                                    # manual_clabels_regular = {
                                    #     FEC_w_ticks[0]: (15,18),
                                    #     FEC_w_ticks[1]: (32,50),
                                    #     FEC_w_ticks[2]: (55,45),
                                    #     FEC_w_ticks[3]: (55,30),
                                    #     FEC_w_ticks[4]: (70,28),
                                    #     FEC_w_ticks[5]: (80,20),
                                    #     FEC_w_ticks[6]: (60,20),
                                    #     },
                                    
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    bottom_fill_bounds = ((0,0), 
                                                          (1,5.),
                                                          (99,7.)),
                                    # bottom_fill_appears_above_other_areas_to_fill=True,
                                    
                                    # other_areas_to_fill_color = colors.CABBI_orange.shade(1).RGBn,
                                    # other_areas_to_fill_bounds= ((0,0), 
                                    #                       (5,10.),
                                    #                       (95,7.)),
                                    
                                    text_boxes = {'>100': [(80, 5), 'white']},
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )
    
    #%% AOC
    
    # AOC_w_levels, AOC_w_ticks, AOC_cbar_ticks = get_contour_info_from_metric_data(results_metric_4,)
    AOC_w_levels = np.arange(0, 151., 5.)
    AOC_cbar_ticks = AOC_w_levels = np.arange(0, 151., 10.)
    AOC_w_ticks = [0, 25, 30, 40, 50, 100, 150]
    # AOC_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_4, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., AOC
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=AOC_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=AOC_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=AOC_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=AOC_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=AOC_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='AOC_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    # comparison_range=TAL_maximum_viable_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 4,
                                    # comparison_range=[AOC_w_levels[-2], AOC_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',fill_bottom_with_cmap_over_color=False, # for TRY
                                    
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    bottom_fill_bounds = ((0,0), 
                                                          (1,11.),
                                                          (99,11.)),
                                    # zoom_data_scale=5,
                                    text_boxes = {'>150': [(80,3), 'white']},
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )
    
    #%% TCI
    
    TCI_w_levels, TCI_w_ticks, TCI_cbar_ticks = get_contour_info_from_metric_data(results_metric_5,)
    # TCI_w_levels = np.arange(2, 8.1, 0.2)
    # TCI_cbar_ticks = np.arange(2, 8.1, 1.)
    
    # TCI_w_ticks = [150, 200, 300, 400,]
    TCI_w_levels = np.arange(150, 401, 10)
    TCI_cbar_ticks = [150, 200, 250, 300, 350, 400]
    TCI_w_ticks = [150, 200, 220, 250, 300, 350, 400]
    # TCI_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_5, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., TCI
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=TCI_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=TCI_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=TCI_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=TCI_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=TCI_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='TCI_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    # comparison_range=TAL_maximum_viable_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 4,
                                    # comparison_range=[TCI_w_levels[-2], TCI_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',fill_bottom_with_cmap_over_color=False, # for TRY
                                    
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    fill_bottom_with_cmap_over_color=False, # for TRY
                                    bottom_fill_bounds = ((0,0), 
                                                          (1,11.),
                                                          (99,11.)),
                                    # zoom_data_scale=5,
                                    text_boxes = {'>400': [(80,5), 'white']},
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )
    
    #%% Purity
    
    # Purity_w_levels, Purity_w_ticks, Purity_cbar_ticks = get_contour_info_from_metric_data(results_metric_6,)
    Purity_w_levels = 100.*np.arange(0., 1.01, 0.025)
    Purity_cbar_ticks = 100.*np.arange(0., 1.01, 0.1)
    # Purity_cbar_ticks = np.arange(2, 8.1, 1.)
    Purity_w_ticks = 100.*np.array([0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95])
    # Purity_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=100.*np.array(results_metric_6), # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., Purity
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=productivities, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=Purity_w_label, # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=Purity_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=Purity_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=Purity_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='neither',
                                    cbar_ticks=Purity_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='Purity_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    # comparison_range=TAL_maximum_viable_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 4,
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    # comparison_range=[Purity_w_levels[-2], Purity_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )
    
    #%% Relative impact of yield and titer on MPSP
    rel_impact = []
    MPSPs = results_metric_1[0]
    t_step = 1
    y_step = 1
    tot_t_steps, tot_y_steps = len(MPSPs), len(MPSPs[0]) # only works if tot_t_steps == tot_y_steps
    max_increment = 1
    for i in range(tot_t_steps-1):
        rel_impact.append([])
        for j in range(tot_y_steps-1):
            t_increments = range(1, min(max_increment+1, tot_t_steps-i))
            y_increments = range(1, min(max_increment+1, tot_t_steps-j))
            ri_sum = 0.
            for inc in zip(t_increments, y_increments): # ensures stopping at the lower maximum increment
                t_inc, y_inc = inc
                curr_MPSP = MPSPs[i][j]
                t_MPSP = MPSPs[i+t_inc][j]
                y_MPSP = MPSPs[i][j+y_inc]
                ri_sum += (curr_MPSP-t_MPSP)/(curr_MPSP-y_MPSP)
            rel_impact[-1].append(ri_sum/len(t_increments)*len(y_increments))
    
    
    for i in range(len(MPSPs)-1): rel_impact[i].append(np.nan)
    rel_impact.append([])
    for i in range(len(MPSPs[i])): rel_impact[-1].append(np.nan)
    
    
    rel_impact = np.array(rel_impact)
    
    Rel_impact_MPSP_w_levels = np.arange(0., 10.01, 0.25)
    Rel_impact_MPSP_cbar_ticks = np.arange(0., 10.01, 0.5)
    # Rel_impact_MPSP_cbar_ticks = np.arange(2, 8.1, 1.)
    Rel_impact_MPSP_w_ticks = np.arange(0., 10.01, 0.5)
    # Rel_impact_MPSP_w_levels = np.arange(0., 15.5, 0.5)
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=np.array([rel_impact]), # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., Rel_impact_MPSP
                                    x_data=100*yields, # x axis values
                                    # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                    y_data=titers, # y axis values
                                    z_data=[productivities[0]], # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label='Rel impact by titer:yield on MPSP', # title of the color axis
                                    x_ticks=100*x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=Rel_impact_MPSP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=Rel_impact_MPSP_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units='',
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='neither',
                                    cbar_ticks=Rel_impact_MPSP_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='Rel_impact_MPSP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    # comparison_range=TAL_maximum_viable_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 4,
                                    additional_points ={(40.5, 35.9):('D', 'w', 6)},
                                    # comparison_range=[Rel_impact_MPSP_w_levels[-2], Rel_impact_MPSP_w_levels[-1]],
                                    # comparison_range_hatch_pattern='////',
                                    
                                    add_shapes = {
                                        # coords as tuple of tuples: (color, zorder),
                                        ((1,0), (47,100), (1,100)): ('white', 2), # infeasible region smoothing
                                        }
                                    )

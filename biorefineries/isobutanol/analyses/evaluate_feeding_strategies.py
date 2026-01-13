#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
from biorefineries import isobutanol

from matplotlib import pyplot as plt


from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os


import biosteam as bst

system = corn_EtOH_IBO_sys = isobutanol.system.corn_EtOH_IBO_sys
load_simulate_get_EtOH_MPSP = isobutanol.system.load_simulate_get_EtOH_MPSP
model_specification = isobutanol.system.model_specification
baseline_spec = isobutanol.system.baseline_spec
tea = corn_EtOH_IBO_sys_tea = isobutanol.system.corn_EtOH_IBO_sys_tea
fbs_spec = isobutanol.system.fbs_spec
optimize_tau_for_MPSP = isobutanol.system.optimize_tau_for_MPSP
optimize_max_n_glu_spikes_for_MPSP = isobutanol.system.optimize_max_n_glu_spikes_for_MPSP

f = system.flowsheet

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')

ferm_reactor = f.V406
HXN = f.HXN1001

product = f.ethanol
broth = ferm_reactor.outs[1]

EtOH_market_range=np.array([0.7, 1.0]) 
                
#%% Filepaths
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
isobutanol_results_filepath = isobutanol_filepath + '\\analyses\\results\\'


#%% Load baseline

#%% Baseline -- simulate and solve TEA
model_specification(**baseline_spec,
    n_sims=3,
    n_tea_solves=3,
    plot=True,
    )

#%%  Metrics
product_chemical_IDs = ['Ethanol',]
get_product_MPSP = lambda: tea.solve_price(product) / get_product_purity() # USD / pure-kg
get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])
get_product_recovery = lambda: sum([product.imol[i] for i in product_chemical_IDs])/sum([broth.imol[i] for i in product_chemical_IDs])
get_AOC = lambda: tea.AOC / 1e6 # million USD / y
get_TCI = lambda: tea.TCI / 1e6 # million USD

get_yield_nsk = lambda: ferm_reactor.results_specific_tau_dict['y_EtOH_glu_added']
get_titer_nsk = lambda: ferm_reactor.results_specific_tau_dict['[s_EtOH]']
get_prod_nsk = lambda: ferm_reactor.results_specific_tau_dict['prod_EtOH']

get_tau = lambda: ferm_reactor.tau

metrics = [get_product_MPSP, 
            get_AOC,
            get_TCI,
            get_yield_nsk,
            get_titer_nsk,
            get_prod_nsk,
            get_tau,]

#%%
results = {i: [] for i in range(len(metrics))}

# %% Generate 3-specification meshgrid and set specification loading functions

steps = (10, 10, 1)

spec_1 = threshold_conc_sugarses = np.linspace(0., 500., steps[0])

spec_2 = target_conc_sugarses = np.linspace(10., 500., steps[1])


spec_3 = conc_sugars_feed_spikes =\
    np.array([
              1.*baseline_spec['conc_sugars_feed_spike'],
              ])

#%% Plot stuff

# Parameters analyzed across

x_label = r"$\bfThreshold glucose concentration$" # title of the x axis
x_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
x_ticks = [0, 100, 200, 300, 400, 500]

y_label = r"$\bfTarget glucose concentration$" # title of the x axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 100, 200, 300, 400, 500]

z_label = r"$\bfSpike feed glucose concentration$" # title of the x axis
z_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
z_ticks = [0, 200, 400, 600, 800]

# Metrics
MPSP_w_label = r"$\bfMPSP$" # title of the color axis
MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
# MPSP_units = r"$\mathrm{\$/kg}$"

AOC_w_label = r"$\bfAOC$" # title of the color axis
AOC_units = r"$\mathrm{MM\$}\cdot\mathrm{y}^{-1}$"
# AOC_units = r"$\mathrm{MM\$/y}$"

TCI_w_label = r"$\bfTCI$" # title of the color axis
TCI_units = r"$\mathrm{MM\$}$"

Yield_w_label = r"$\bfYield$" # title of the color axis
Yield_units = r"$\mathrm{g}\cdot\mathrm{g}^{-1}$"

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

def JBEI_UCB_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    JBEI_orange = (233/255, 83/255, 39/255)
    UCB_blue = (0/255, 38/255, 118/255)
    UCB_yellow = (253/255, 181/255, 21/255)
    cmap_colors = (
                    UCB_yellow,
                    JBEI_orange,
                    UCB_blue,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)

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
file_to_save = f'_{steps}_steps_'+'etoh_fbs_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)

#%% Initial simulation

print('\n\nSimulating the initial point to avoid bugs ...')
curr_spec = fbs_spec.current_specifications
curr_spec.update({'threshold_conc_sugars':threshold_conc_sugarses[0],})
curr_spec.update({'target_conc_sugars':target_conc_sugarses[0],})

model_specification(**curr_spec,
    n_sims=3,
    n_tea_solves=3,
    plot=True,
    )

# %% Run TRY analysis 

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
total_no = len(spec_1)*len(spec_2)*len(spec_3)

print_status_every_n_simulations = 1

for s3 in spec_3:
    for v in list(results.values()): v.append([])
    
    for s1 in spec_1:
        for v in list(results.values()): v[-1].append([])
        for s2 in spec_2:
            curr_no +=1
            error_message = None
            try:
                curr_spec = {k:v for k,v in baseline_spec.items()}
                curr_spec.update({'threshold_conc_sugars':s1,})
                curr_spec.update({'target_conc_sugars':s2,})
                curr_spec.update({'conc_sugars_feed_spike':s3,})
                
                model_specification(**curr_spec,
                    n_sims=3,
                    n_tea_solves=3,
                    plot=False,
                    )
                optimize_max_n_glu_spikes_for_MPSP()
                
                assert s1<s2
                for k, v in list(results.items()): 
                    v[-1][-1].append(metrics[k]())
                
                HXN_qbal_error = HXN.energy_balance_percent_error
                if abs(max_HXN_qbal_percent_error)<abs(HXN_qbal_error): max_HXN_qbal_percent_error = HXN_qbal_error
            
            except Exception as e:
                str_e = str(e).lower()
                print('Error in model spec: %s'%str_e)
                for v in list(results.values()): v[-1][-1].append(np.nan)
                error_message = str_e
            
            if curr_no%print_status_every_n_simulations==0 or error_message:
                print_status(curr_no, total_no,
                             s1, s2, s3, 
                             results=[v[-1][-1][-1] for v in list(results.values())],
                             HXN_qbal_error=HXN.energy_balance_percent_error,
                             exception_str=error_message)

    # Convert last 2D list to array and transpose
    for k in results.keys(): 
        results[k][-1] = np.array(results[k][-1]).transpose()

    # Save generated data
    for k, v in results.items():
        csv_file_to_save = file_to_save + f'_metric_{k}'
        pd.DataFrame(v[-1]).to_csv(isobutanol_results_filepath+'MPSP-'+csv_file_to_save+'.csv')

#%% Report maximum HXN energy balance error
print(f'Max HXN Q bal error was {round(max_HXN_qbal_percent_error, 3)} %.')

#%% 

chdir(isobutanol_results_filepath)

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

#%% Smoothing
smoothing = False

if smoothing:
    for arr in list(results.values()):
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                for k in range(arr.shape[2]):
                    if j>0 and k>0 and j<arr.shape[1]-1 and k<arr.shape[2]-1 :
                        if np.isnan(arr[i,j,k]):
                            manhattan_neighbors = np.array([
                                         # arr[i][j-2][k],
                                         # arr[i][j+2][k],
                                         arr[i][j][k-1],
                                         arr[i][j][k+1]
                                         ])
                            if not np.any(np.isnan(manhattan_neighbors)):
                                arr[i,j,k] = np.mean(manhattan_neighbors)
                        # else:
                        #     manhattan_neighbors = np.array([
                        #                  arr[i][j-1][k],
                        #                  arr[i][j+1][k],
                        #                  arr[i][j][k-1],
                        #                  arr[i][j][k+1]
                        #                  ])
                        #     if not np.any(np.isnan(manhattan_neighbors)):
                        #         if not round(arr[i,j,k]/np.mean(manhattan_neighbors),0)==1:
                        #             print(i,j,k)
                    
#%% Plots
plot = True

if plot: 
    
    #%% MPSP
    
    # MPSP_w_levels, MPSP_w_ticks, MPSP_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
    MPSP_w_levels = np.arange(0.75, 1.51, 0.025)
    MPSP_cbar_ticks = np.arange(0.75, 1.51, 0.25)
    MPSP_w_ticks = [0.7, 0.8, 0.9]
    # MPSP_w_levels = np.arange(0., 15.5, 0.5)
    
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[0], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                    x_data=spec_1, # x axis values
                                    # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                    y_data=spec_2, # y axis values
                                    z_data=spec_3, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=MPSP_w_label, # title of the color axis
                                    x_ticks=x_ticks,
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
                                    cmap=JBEI_UCB_colormap(), # can use 'viridis' or other default matplotlib colormaps
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
                                    # comparison_range=EtOH_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    units_on_newline = (False, False, False, False), # x,y,z,w
                                    units_opening_brackets = [" (",] * 4,
                                    units_closing_brackets = [")",] * 4,
                                    )
    
    #%% Yield
    
    # Yield_w_levels, Yield_w_ticks, Yield_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
    Yield_w_levels = np.arange(0., 0.5, 0.01)
    Yield_cbar_ticks = np.arange(0., 0.5, 0.05)
    Yield_w_ticks = [0.1, 0.2, 0.3, 0.4, 0.5]
    # Yield_w_levels = np.arange(0., 15.5, 0.5)
    
    
    contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[3], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., Yield
                                    x_data=spec_1, # x axis values
                                    # x_data = yields/theoretical_max_g_HP_acid_per_g_glucose,
                                    y_data=spec_2, # y axis values
                                    z_data=spec_3, # z axis values
                                    x_label=x_label, # title of the x axis
                                    y_label=y_label, # title of the y axis
                                    z_label=z_label, # title of the z axis
                                    w_label=Yield_w_label, # title of the color axis
                                    x_ticks=x_ticks,
                                    y_ticks=y_ticks,
                                    z_ticks=z_ticks,
                                    w_levels=Yield_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                    w_ticks=Yield_w_ticks, # labeled, lined contours; a subset of w_levels
                                    x_units=x_units,
                                    y_units=y_units,
                                    z_units=z_units,
                                    w_units=Yield_units,
                                    # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                    cmap=JBEI_UCB_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                    cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                    extend_cmap='max',
                                    cbar_ticks=Yield_cbar_ticks,
                                    z_marker_color='g', # default matplotlib color names
                                    fps=fps, # animation frames (z values traversed) per second
                                    n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                    animated_contourplot_filename='yield_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                    keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                    axis_title_fonts=axis_title_fonts,
                                    clabel_fontsize = clabel_fontsize,
                                    default_fontsize = default_fontsize,
                                    axis_tick_fontsize = axis_tick_fontsize,
                                    # comparison_range=EtOH_market_range,
                                    n_minor_ticks = 1,
                                    cbar_n_minor_ticks = 3,
                                    units_on_newline = (False, False, False, False), # x,y,z,w
                                    units_opening_brackets = [" (",] * 4,
                                    units_closing_brackets = [")",] * 4,
                                    )
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np

from matplotlib import pyplot as plt

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

from biorefineries import isobutanol

import math

fbs_spec = isobutanol.models.fbs_spec

#%%
EtOH_market_range=np.array([
    0.52, # 1.5475 $/gal/(3.7854 L/gal * 0.789 kg/L)
    1.15, # 3.4500 $/gal/(3.7854 L/gal * 0.789 kg/L)
    ]) # Jan 2021 - Dec 2025 5-year low and high from https://tradingeconomics.com/commodity/ethanol
               
#%% Pub results main filepath
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'

#%%  Metrics
metrics_units = {"MPSP":  r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
                "AOC": "MM$/y",
                "TCI": "MM$",
                "Combined Yield": "g-EtOH-and-IBO/g-sugars",
                "EtOH Titer": "g-EtOH/L-broth",
                "EtOH Productivity": "g-EtOH/L-broth/h",
                "Number of glucose spikes": "",
                "Fermentation time": "h",
                "Total Q sugar evap": "kJ/h",
                "Target sugars concentration": "g-sugars/L-broth",
                "Cell loading": "g-cell/L-broth",
                "Active cell loading": "g-cell/L-broth",
                "EtOH Yield":"g-EtOH/g-sugars",
                "IBO Yield": "g-IBO/g-sugars",
                "IBO Titer": "g-IBO/L-broth",
                "IBO Productivity": "g-IBO/L-broth/h",
                'Actual aeration required': 'kmol-O2/h',
                }

#%% Input Details
# !!!

# Spec names
x_label = "Threshold glucose concentration" # title of the x axis
x_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
x_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

y_label = "Target glucose concentration" # title of the x axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

z_label = "Spike feed glucose concentration" # title of the x axis
z_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
z_ticks = [0, 200, 400, 600, 800]

# Misc
steps = (25, 25, 1)

#
spec_1 = threshold_conc_sugarses = np.linspace(1., 400., steps[0])

spec_2 = target_conc_sugarses = np.linspace(10., 400., steps[1])


spec_3 = conc_sugars_feed_spikes =\
    np.array([
              # 1.*baseline_spec['conc_sugars_feed_spike'],
              fbs_spec.conc_sugars_feed_spike,
              ])
    
#%% Intermediate details

# Misc
subfolder_name = f'Feed-strat\\'

#%% Metric names for plots

metrics_plot_names = {k: k for k in metrics_units.keys()}
metrics_plot_names["Total Q sugar evap"] = "Slurry evaporation duty"
metrics_plot_names["Actual aeration required"] = "Aeration required"

for k in metrics_plot_names.keys():
    metrics_plot_names[k] = metrics_plot_names[k].replace("IBO", "Isobutanol").replace("EtOH", "Ethanol")

metrics_plot_names["Combined Yield"] = "Ethanol Yield"
metrics_plot_names["Cell loading"] = "Cell density"

for k in metrics_plot_names.keys():
    if not k in ('TCI', 'AOC', 'MPSP',):
        metrics_plot_names[k] = metrics_plot_names[k].lower()
        
#%% Chdir

os.chdir(isobutanol_results_pub_filepath + subfolder_name)

#%% Filename
file_to_load = 'ibo_(25, 25, 1)_Thres_Targe_Spike_'

#%% Load results
results = {}
for k in metrics_units.keys():
    try:
        df = pd.read_csv(file_to_load + f'_{k}.csv')
        results[k] = np.array([df[df.columns[1:]].to_numpy(dtype='float', na_value=np.nan)])
    except Exception as e:
        if 'no such file' in str(e).lower():
            print(f'Metric {k} file not found.')
        else:
            raise e
            
#%% Colors

def JBEI_UCB_colormap(N_levels=90, reverse=False):
    JBEI_orange = (233/255, 83/255, 39/255)
    UCB_blue = (0/255, 38/255, 118/255)
    UCB_yellow = (253/255, 181/255, 21/255)
    cmap_colors = [
                    UCB_yellow,
                    JBEI_orange,
                    UCB_blue,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn]
    if reverse: cmap_colors.reverse()
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

#%% More plot stuff

fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5
keep_frames = True
keep_gifs = False

x_label_for_plot = r"$\bf{"+x_label.replace(' ', '\ ')+"}$"
y_label_for_plot = r"$\bf{"+y_label.replace(' ', '\ ')+"}$"

#%% Plots
plot_all_generic = True

if plot_all_generic: 
    print('\nCreating and saving contour plots ...\n')
    
    # Coords for markers

    # metrics_to_opt = [
    #                   'Combined Yield', 
    #                   'EtOH Titer', 
    #                   'EtOH Productivity', 
    #                   'Cell loading',
    #                   'TCI',
    #                   'AOC',
    #                   'MPSP', 
    #                   'Total Q sugar evap',
    #                   'Actual aeration required',
    #                   ]
    metrics_to_opt = [
                      'Cell loading',
                      'EtOH Titer', 
                      'EtOH Productivity', 
                      'Combined Yield', 
                      'Total Q sugar evap',
                      'Actual aeration required',
                      'TCI',
                      'AOC',
                      'MPSP', 
                      ]

    opt_coords = {}
    for m in  metrics_to_opt:
        m_non_nans = np.array(results[m])[np.where(~np.isnan(np.array(results[m])))]
        if m in ('MPSP', 'AOC', 'TCI', 'Total Q sugar evap', 'Fermentation time', 'Actual aeration required'):
            m_opt_coords = np.where(results[m][0]==m_non_nans.min())
        else:
            m_opt_coords = np.where(results[m][0]==m_non_nans.max())
        opt_s1 = spec_1[m_opt_coords[1][0]]
        opt_s2 = spec_2[m_opt_coords[0][0]]
        opt_coords[m] = (opt_s1, opt_s2)
        # print(results['MPSP'][0][opt_coords])
        # print(opt_s1, opt_s2)
    
    additional_points = {}
    
    if x_label in ('k_1e',):
        baseline_coords = (47.1, 0.04)
        additional_points[baseline_coords] = ('D', 'gray', 6)
    elif x_label in ('k_13',):
        # baseline_coords = (5.81, 0.04)
        baseline_coords = (0.0, 0.04)
        # additional_points[baseline_coords] = ('D', 'gray', 6)
    
    opt_marker_shapes = ['o', '^', 's', 'p', 'v', '<', '>', 'h', 
                         # 'P', 'X',
                         ]
    if x_label in ("k_1e",):
        for i in ['v', '<', '>']:
            opt_marker_shapes.remove(i)
            
    metric_optima_markers = {}
    
    shapes_i = 0
    for m in metrics_to_opt:
        if m=='MPSP':
            marker_shape='*'
            marker_color='#33ccff'
            marker_size = 8
        elif m=='TCI':
            marker_shape='s'
            marker_color='#33ccff'
            marker_size = 6
        elif m=='AOC':
            marker_shape='p'
            marker_color='#33ccff'
            marker_size = 6
        elif m=='Total Q sugar evap':
            marker_shape='P'
            marker_color = 'w'
            marker_size = 6
        elif m=='Actual aeration required':
            marker_shape='X'
            marker_color = 'w'
            marker_size = 6
        else:
            marker_shape = opt_marker_shapes[shapes_i]
            shapes_i += 1
            marker_color = 'w'
            marker_size = 6
            
        metric_opt_name = metrics_plot_names[m]
        if not metric_opt_name in ('TCI', 'AOC', 'MPSP',):
            metric_opt_name = metric_opt_name.lower()
            
        metric_optima_markers[metric_opt_name] = (marker_shape, marker_color, marker_size)
        additional_points[opt_coords[m]] = (marker_shape, marker_color, marker_size)
        print(m, marker_shape, marker_color, marker_size)
    
    fig, ax = plt.subplots()
    contourplots.utils.marker_legend(ax, metric_optima_markers, title="Optima", loc="upper right")
    fig.set_figwidth(5)
    fig.set_figheight(10)
    plt.savefig(
        fname=file_to_load+'_marker_legend.png',
        transparent=False,  
        facecolor='white',
        bbox_inches='tight',
        # figheight=20, figwidth=10,
        dpi=600,)
    plt.close()
        
        
    #%% All metrics
    for curr_metric, val in metrics_units.items():
        if not curr_metric in results.keys(): continue
        cbar_n_minor_ticks = 3
        lccm = curr_metric.lower()

        if 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'loading' in lccm:
            cmap = JBEI_UCB_colormap(reverse=True)
            cmap_over_color = colors.yellow_tint.RGBn
        
        else:
            cmap = JBEI_UCB_colormap(reverse=False)
            cmap_over_color = colors.grey_dark.shade(8).RGBn
            
        # curr_metric_w_levels, curr_metric_w_ticks, curr_metric_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
        curr_metric_non_nans = np.array(results[curr_metric])[np.where(~np.isnan(np.array(results[curr_metric])))]
        
        curr_metric_w_levels = np.arange(curr_metric_non_nans.min(), 
                                      curr_metric_non_nans.max()*1.001, 
                                      (curr_metric_non_nans.max()-curr_metric_non_nans.min())/80
                                      )
        curr_metric_cbar_ticks = np.arange(curr_metric_non_nans.min(), 
                                      curr_metric_non_nans.max()*1.001, 
                                      (curr_metric_non_nans.max()-curr_metric_non_nans.min())/5
                                      )
        
        curr_metric_w_ticks = list(set([np.percentile(curr_metric_non_nans, 25),
                            np.percentile(curr_metric_non_nans, 50),
                            np.percentile(curr_metric_non_nans, 75),
                            curr_metric_non_nans.max()]))
        curr_metric_w_ticks.sort(reverse=False)
        # curr_metric_w_levels = np.arange(0., 15.5, 0.5)

        # else:
        #     break
        
        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[curr_metric], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., curr_metric
                                        x_data=spec_1, # x axis values
                                        # x_data = curr_metrics/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=spec_2, # y axis values
                                        z_data=spec_3, # z axis values
                                        x_label=x_label_for_plot, # title of the x axis
                                        y_label=y_label_for_plot, # title of the y axis
                                        z_label=r"$\bf"+z_label+"$", # title of the z axis
                                        w_label=r"$\bf"+curr_metric+"$", # title of the color axis
                                        x_ticks=x_ticks,
                                        y_ticks=y_ticks,
                                        z_ticks=z_ticks,
                                        w_levels=curr_metric_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                        w_ticks=curr_metric_w_ticks, # labeled, lined contours; a subset of w_levels
                                        x_units=x_units,
                                        y_units=y_units,
                                        z_units=z_units,
                                        w_units=val,
                                        # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                        fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                        cmap=cmap, # can use 'viridis' or other default matplotlib colormaps
                                        # cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                        cmap_over_color=cmap_over_color,
                                        extend_cmap='max',
                                        cbar_ticks=curr_metric_cbar_ticks,
                                        z_marker_color='g', # default matplotlib color names
                                        fps=fps, # animation frames (z values traversed) per second
                                        n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                        animated_contourplot_filename=file_to_load+f'_{curr_metric}', # file name to save animated contourplot as (no extensions)
                                        keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                        keep_gifs=keep_gifs,
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        # comparison_range=EtOH_market_range,
                                        n_minor_ticks = 3,
                                        additional_points=additional_points,
                                        cbar_n_minor_ticks = cbar_n_minor_ticks,
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" [",] * 4,
                                        units_closing_brackets = ["]",] * 4,
                                        round_xticks_to=0,
                                        round_yticks_to=1,
                                        )

#%% Tailor-made plots

# !!!

#%% MPSP
curr_metric = 'MPSP'
curr_metric_w_levels = np.arange(0.7, 0.9, 0.01)
curr_metric_cbar_ticks = np.arange(0.7, 0.9, 0.05)
curr_metric_w_ticks = []
cbar_n_minor_ticks = 4
lccm = curr_metric.lower()

if 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'loading' in lccm:
    cmap = JBEI_UCB_colormap(reverse=True)
    cmap_over_color = colors.yellow_tint.RGBn

else:
    cmap = JBEI_UCB_colormap(reverse=False)
    cmap_over_color = colors.grey_dark.shade(8).RGBn


val = metrics_units[curr_metric]

where_batch_mode = np.where(results['Number of glucose spikes'][0]==0)
coords_batch_mode = [(spec_1[where_batch_mode[1][i]], spec_2[where_batch_mode[0][i]]) for i in range(len(where_batch_mode[0]))]

# # Reorder coords to form a closed polygon
# # 1. Compute centroid
# cx = sum(x for x, y in coords_batch_mode) / len(coords_batch_mode)
# cy = sum(y for x, y in coords_batch_mode) / len(coords_batch_mode)

# # 2. Sort by angle from centroid
# ordered_coords_batch_mode = sorted(coords_batch_mode, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))

# ordered_coords_batch_mode = tuple(ordered_coords_batch_mode)

def get_zeroth_spec_2_val_for_condition_for_all_spec_1_vals(metric_arr, condition):
    s2s = []
    for s1i in range(len(metric_arr[0])):
        success = False
        for s2i in range(len(metric_arr)):
            if condition(metric_arr[s2i][s1i]):
                s2s.append(spec_2[s2i])
                success = True
                break
        if not success:
            s2s.append(np.nan)
    return s2s

line_first_app_n_glu_spikes_0 = tuple(get_zeroth_spec_2_val_for_condition_for_all_spec_1_vals(results['Number of glucose spikes'][0],
                                                                  condition = lambda i: i==0))

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[curr_metric], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., curr_metric
                                x_data=spec_1, # x axis values
                                # x_data = curr_metrics/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=spec_2, # y axis values
                                z_data=spec_3, # z axis values
                                x_label=x_label_for_plot, # title of the x axis
                                y_label=y_label_for_plot, # title of the y axis
                                z_label=r"$\bf"+z_label+"$", # title of the z axis
                                w_label=r"$\bf"+curr_metric+"$", # title of the color axis
                                x_ticks=x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=curr_metric_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=curr_metric_w_ticks, # labeled, lined contours; a subset of w_levels
                                x_units=x_units,
                                y_units=y_units,
                                z_units=z_units,
                                w_units=val,
                                # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                cmap=cmap, # can use 'viridis' or other default matplotlib colormaps
                                # cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                cmap_over_color=cmap_over_color,
                                extend_cmap='max',
                                cbar_ticks=curr_metric_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename=file_to_load+f'_{curr_metric}', # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                keep_gifs=keep_gifs,
                                axis_title_fonts=axis_title_fonts,
                                clabel_fontsize = clabel_fontsize,
                                default_fontsize = default_fontsize,
                                axis_tick_fontsize = axis_tick_fontsize,
                                additional_points=additional_points,
                                # comparison_range=EtOH_market_range,
                                n_minor_ticks = 1,
                                cbar_n_minor_ticks = cbar_n_minor_ticks,
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                round_xticks_to=0,
                                round_yticks_to=1,
                                add_lines={line_first_app_n_glu_spikes_0: {'color': 'white', 'linewidth': 0.6, 'alpha': 0.6}},
                                )

#%% All metric values at opt_coords
opt_coords_metric_vals = {}

for m1 in metrics_to_opt:
    opt_s1, opt_s2 = opt_coords[m1]
    opt_s1_ind, opt_s2_ind = np.where(spec_1==opt_s1)[0][0], np.where(spec_2==opt_s2)[0][0]
    opt_coords_metric_vals[m1] = {'Coords': (opt_s1, opt_s2),}
    for m2 in metrics_to_opt:
        opt_coords_metric_vals[m1][m2] = results[m2][0, opt_s2_ind, opt_s1_ind]

#%%
round_off = contourplots.utils.round_off
def get_optima_comparisons(opt_coords_metric_values, rel_to_m='MPSP'):
    print(f"\n\nRelative to the optimum for '{rel_to_m}', at the optimum for:")
    print('\n')
    i = 0
    for m1, v in opt_coords_metric_values.items():
        i+=1
        print(f"{i}. '{m1}' ({v['Coords']}),")
        for m2 in v.keys():
            if not m2=='Coords':
                try:
                    rel_diff = v[m2]/opt_coords_metric_values[rel_to_m][m2] - 1
                    sign = '+' if rel_diff>0 else '-'
                    print(f"'{m2}' is {round_off(v[m2],3)}, which is {sign} {abs(int(100*rel_diff))}%.")
                except Exception as e:
                    if 'divide' in str(e).lower():
                        print(f"'{m2}' is zero.")
                    else:
                        breakpoint()
        print('\n')
    print('\n Note that optima overlap for:')
    print('\n')
    coords_metrics = {}
    for metric, v in opt_coords_metric_values.items():
        if not v['Coords'] in coords_metrics.keys():
            coords_metrics[v['Coords']] = [metric]
        else:
            coords_metrics[v['Coords']].append(metric)
    for coords, metrics in coords_metrics.items():
        if len(metrics)>1:
            print(f'{metrics}: optimum coordinates are {coords}.')
        
#%%
get_optima_comparisons(opt_coords_metric_vals)

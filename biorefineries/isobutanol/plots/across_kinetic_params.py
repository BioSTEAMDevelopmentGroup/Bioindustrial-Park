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
                "AOC": "MM\$/y",
                "TCI": "MM\$",
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
                "Actual aeration required": "kmol-O2/h",
                }

#%% Input Details
# !!!

# Spec names
x_label = "k_13"
y_label = "k_7ii"
z_label = "Spike feed glucose concentration"

# Misc
strategy='adaptive_fed-batch'
steps = (25, 25, 1)

#%% Intermediate details

# Misc
subfolder_name = f'Kinetics_{x_label}_{y_label}\\{strategy}\\'
perform_feeding_strategy_opt = 'adaptive' in subfolder_name
max_n = 0 if not 'fed-batch' in subfolder_name else 21

## Spec arrays, units, and ticks
spec_1, spec_2, spec_3 = None, None, None
x_units, y_units, z_units = None, None, None
x_ticks, y_ticks, z_ticks = None, None, None

# spec_1
if x_label in ('k_1e',):
    spec_1 = np.linspace(1., 300., steps[0])
    x_units = r"$\mathrm{g} \cdot \mathrm{g}^{-1} \cdot \mathrm{h}^{-1}$"
    x_ticks = np.array([0, 100, 200, 300])

elif x_label in ('k_13',):
    spec_1 = np.linspace(0.0, 40., steps[0])
    x_units = r"$\mathrm{g} \cdot \mathrm{g}^{-1} \cdot \mathrm{h}^{-1}$"
    x_ticks = np.array([0, 10, 20, 30, 40])

# spec_2
if y_label in ('k_1ie', 'k_1ii', 'k_7ie', 'k_7ii',):
    spec_2 = np.linspace(0.0001, 0.5, steps[1])
    y_units = r"$\mathrm{L} \cdot \mathrm{g}^{-1}$"
    y_ticks = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

# spec_3
if z_label in ('Spike feed glucose concentration',):
    spec_3 = conc_sugars_feed_spikes =\
        np.array([
                  fbs_spec.conc_sugars_feed_spike,
                  ])
    z_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
    z_ticks = [0, 200, 400, 600, 800]

#%% Metric names for plots

metrics_plot_names = {k: k for k in metrics_units.keys()}
metrics_plot_names["Total Q sugar evap"] = "Slurry evaporation duty"
metrics_plot_names["Actual aeration required"] = "Aeration required"

for k in metrics_plot_names.keys():
    metrics_plot_names[k] = metrics_plot_names[k].replace("IBO", "Isobutanol").replace("EtOH", "Ethanol")
    
if x_label in ("k_1e",):
    metrics_plot_names["Combined Yield"] = "Ethanol Yield"
   
metrics_plot_names["Cell loading"] = "Cell density"

for k in metrics_plot_names.keys():
    if not k in ('TCI', 'AOC', 'MPSP',):
        metrics_plot_names[k] = metrics_plot_names[k].lower()
    
#%% Chdir

os.chdir(isobutanol_results_pub_filepath + subfolder_name)

#%% Filename
file_to_load = f'ibo_{steps}_{x_label[:5]}_{y_label[:5]}_{z_label[:5]}_opt={perform_feeding_strategy_opt}_max_n={max_n}_'

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

x_label_for_plot = r"$\bf{" + x_label[:x_label.index('_')] + '}_{' + x_label[x_label.index('_')+1:] + '}' +"$"
y_label_for_plot = r"$\bf{" + y_label[:y_label.index('_')] + '}_{' + y_label[y_label.index('_')+1:] + '}' +"$"

#%% Plots
plot_all_generic = True

if plot_all_generic: 
    print('\nCreating and saving contour plots ...\n')
    
    
    # Coords for markers
    
    if x_label in ('k_1e',):
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
    
    elif x_label in ('k_13',):
        metrics_to_opt = [
                          'Cell loading',
                          'EtOH Titer', 
                          'EtOH Productivity', 
                          'EtOH Yield',
                          'IBO Titer', 
                          'IBO Productivity', 
                          'IBO Yield',
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
            marker_size = 12
        elif m=='TCI':
            marker_shape='s'
            marker_color='#33ccff'
            marker_size = 8
        elif m=='AOC':
            marker_shape='p'
            marker_color='#33ccff'
            marker_size = 8
        elif m=='Total Q sugar evap':
            marker_shape='P'
            marker_color = 'w'
            marker_size = 8
        elif m=='Actual aeration required':
            marker_shape='X'
            marker_color = 'w'
            marker_size = 8
        else:
            marker_shape = opt_marker_shapes[shapes_i]
            shapes_i += 1
            marker_color = 'w'
            marker_size = 8
            
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
        if 'spike' in lccm or 'total q' in lccm or 'target sugars' in lccm:
            if not perform_feeding_strategy_opt: 
                continue
            else: 
                if 'spike' in lccm:
                    if max_n == 0.:
                        continue
                else:
                    pass
        elif 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'loading' in lccm:
            cmap = JBEI_UCB_colormap(reverse=True)
            cmap_over_color = colors.yellow_tint.RGBn
        
        else:
            cmap = JBEI_UCB_colormap(reverse=False)
            cmap_over_color = colors.grey_dark.shade(8).RGBn
        
        # curr_metric_w_levels, curr_metric_w_ticks, curr_metric_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
        curr_metric_non_nans = np.array(results[curr_metric])[np.where(~np.isnan(np.array(results[curr_metric])))]
        
        w_level_top_val = curr_metric_non_nans.max()*0.8
        extend_cmap = 'max'
        if 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'cell loading' in lccm:
            w_level_top_val = curr_metric_non_nans.max()
            extend_cmap = 'neither'
        arange_w_level_top_val = w_level_top_val*1.00000000001
        curr_metric_w_levels = np.arange(curr_metric_non_nans.min(), 
                                      arange_w_level_top_val, 
                                      (w_level_top_val-curr_metric_non_nans.min())/80
                                      )
        curr_metric_cbar_ticks = np.arange(curr_metric_non_nans.min(), 
                                      arange_w_level_top_val, 
                                      (w_level_top_val-curr_metric_non_nans.min())/5
                                      )
        
        curr_metric_w_ticks = list(set([np.percentile(curr_metric_non_nans, 25),
                            np.percentile(curr_metric_non_nans, 50),
                            np.percentile(curr_metric_non_nans, 75),
                            w_level_top_val,
                            ]))
        curr_metric_w_ticks.sort(reverse=False)
        # curr_metric_w_levels = np.arange(0., 15.5, 0.5)
        
        if 'mpsp' in lccm:
            curr_metric_w_levels = np.arange(0.25, 5.001, 0.1)
            curr_metric_cbar_ticks = np.arange(0.25, 5.001, 0.25)
            curr_metric_w_ticks = [0.4, 0.9, 2.5, 5.0]
            cbar_n_minor_ticks = 4
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
                                        w_label=r"$\bf"+metrics_plot_names[curr_metric].replace(' ', '\ ')+"$", # title of the color axis
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
                                        extend_cmap=extend_cmap,
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
                                        # comparison_range=np.array(
                                        #     [max(EtOH_market_range[0], curr_metric_non_nans.min()),
                                        #      min(EtOH_market_range[1], curr_metric_non_nans.max())]),
                                        n_minor_ticks = 3,
                                        cbar_n_minor_ticks = cbar_n_minor_ticks,
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" [",] * 4,
                                        units_closing_brackets = ["]",] * 4,
                                        round_xticks_to=0,
                                        round_yticks_to=1,
                                        additional_points=additional_points,
                                        )

#%% Tailor-made plots

# !!!

#%% MPSP vs k_1e, k_7ie, spike feed conc
curr_metric = 'MPSP'
# if x_label in ('k_1e',):
#     curr_metric_w_levels = np.arange(0.6, 3.001, 0.05)
#     curr_metric_cbar_ticks = np.arange(0.6, 3.001, 0.4)
#     curr_metric_w_ticks = [0.75, 2.0, 3.0]
#     cbar_n_minor_ticks = 3
if x_label in ('k_13', 'k_1e'):
    # curr_metric_w_levels = np.arange(0.4, 4.001, 0.1)
    # curr_metric_cbar_ticks = np.arange(0.4, 4.001, 0.4)
    # curr_metric_w_ticks = [0.55, 2.5, 4.00]
    curr_metric_w_levels = np.arange(0.2, 2.401, 0.05)
    curr_metric_cbar_ticks = np.arange(0.2, 2.401, 0.2)
    curr_metric_w_ticks = [2.4]
    cbar_n_minor_ticks = 3
lccm = curr_metric.lower()

if 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'loading' in lccm:
    cmap = JBEI_UCB_colormap(reverse=True)
    cmap_over_color = colors.yellow_tint.RGBn

else:
    cmap = JBEI_UCB_colormap(reverse=False)
    cmap_over_color = colors.grey_dark.shade(8).RGBn


val = metrics_units[curr_metric]

curr_metric_non_nans = np.array(results[curr_metric])[np.where(~np.isnan(np.array(results[curr_metric])))]
 
contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[curr_metric], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., curr_metric
                                x_data=spec_1, # x axis values
                                # x_data = curr_metrics/theoretical_max_g_HP_acid_per_g_glucose,
                                y_data=spec_2, # y axis values
                                z_data=spec_3, # z axis values
                                x_label=x_label_for_plot, # title of the x axis
                                y_label=y_label_for_plot, # title of the y axis
                                z_label=r"$\bf"+z_label+"$", # title of the z axis
                                w_label=r"$\bf"+metrics_plot_names[curr_metric].replace(' ', '\ ')+"$", # title of the color axis
                                x_ticks=x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=curr_metric_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=curr_metric_w_ticks, # labeled, lined contours; a subset of w_levels
                                x_units=x_units,
                                y_units=y_units,
                                z_units=z_units,
                                w_units=val,
                                fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.2f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                # fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
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
                                comparison_range=np.array(
                                    [max(EtOH_market_range[0], curr_metric_non_nans.min()),
                                     min(EtOH_market_range[1], curr_metric_non_nans.max())]),
                                n_minor_ticks = 3,
                                cbar_n_minor_ticks = cbar_n_minor_ticks,
                                units_on_newline = (False, False, False, False), # x,y,z,w
                                units_opening_brackets = [" [",] * 4,
                                units_closing_brackets = ["]",] * 4,
                                round_xticks_to=0,
                                round_yticks_to=1,
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

def get_metric_given_specs(metric, spec_1, spec_2):
    return results[metric][0][spec_2][spec_1]

#%%
get_optima_comparisons(opt_coords_metric_vals)

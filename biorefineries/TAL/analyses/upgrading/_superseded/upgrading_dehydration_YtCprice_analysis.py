# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: sarangbhagwat
"""

# Run this cell first
from warnings import filterwarnings
import copy
filterwarnings('ignore')
from biosteam.process_tools import UnitGroup
from biosteam.digraph import digraph_from_units, save_digraph
from biosteam.utils import streams_from_units, filter_out_missing_streams, colors
from biosteam.report import unit_result_tables, tables_to_excel, stream_table
import numpy as np
import biosteam as bst
from biosteam.utils import colors
# from matplotlib.ticker import AutoMinorLocator as AML

# from TAL.system_solubility_exploit import TAL_sys, TAL_tea, R302, spec
from biorefineries import TAL
from biorefineries.TAL.system_SA_adsorption_sugarcane import TAL_sys, TAL_tea, TAL_lca, R302, KSA_product, u, TAL_chemicals
# from biorefineries.TAL.system_TAL_adsorption_glucose import TAL_sys, TAL_tea, R302, spec, SA
# from biorefineries.TAL.system_ethyl_esters import TAL_sys, TAL_tea, R302, spec, Mixed_esters
# get_GWP, get_non_bio_GWP, get_FEC, get_SPED
# from TAL.system_glucose_to_TAL_adsorb_evap_dry import SA as product

from matplotlib import pyplot as plt
from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from biosteam.plots import  MetricBar, plot_scatter_points, plot_contour_2d #, CABBI_green_colormap
from math import floor, ceil
from datetime import datetime
import flexsolve as flx
from math import log
from biosteam.utils import colors
from biosteam.utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
import matplotlib.colors as mcolors

from datetime import datetime
import os

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')

#%% Filepaths

TAL_filepath = TAL.__file__.replace('\\__init__.py', '')
TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'

#%% Initialize process specification

R401, R402, R403 = u.R401, u.R402, u.R403

fresh_R402_cat_stream = R402.ins[3]
spec_upgrading = TAL._general_process_specification.GeneralProcessSpecification(
    system=TAL_sys,
    baseline_spec_values=[R401.TAL_to_HMTHP_rxn.X, R402.HMTHP_to_PSA_rxn.X, R403.PSA_to_SA_rxn.X],
    HXN = u.HXN1001,
    )

def load_dehydration_PSA_yield(PSA_yield):
    R402.HMTHP_to_PSA_rxn.X = PSA_yield

def load_dehydration_time(dehydration_time):
    R402.tau = dehydration_time
    
def load_dehydration_cat_price(dehydration_cat_price):
    fresh_R402_cat_stream.price = dehydration_cat_price
    if type(dehydration_cat_price) == np.ndarray:
        R402.Amberlyst70_catalyst_price = dehydration_cat_price[0]
    else:
        R402.Amberlyst70_catalyst_price = dehydration_cat_price

spec_upgrading.load_spec_1, spec_upgrading.load_spec_2, spec_upgrading.load_spec_3 = \
    load_dehydration_PSA_yield, load_dehydration_time, load_dehydration_cat_price


spec = spec_upgrading

# %% Generate 3-specification meshgrid and set specification loading functions

steps = (5, 5, 20)

# Yield, titer, productivity (rate)
spec_1 = PSA_yields = np.linspace(0.3, 0.999, steps[0]) # TAL->HMTHP conversion in the hydrogenation reactor
spec_2 = dehydration_times = np.linspace(0.1, 20., steps[1]) # HMTHP->PSA conversion in the dehydration reactor
spec_3 = dehydration_cat_prices = np.linspace(10, 1000, steps[2])# PSA->SA conversion in the ring opening & hydrolysis reactor


spec_1, spec_2 = np.meshgrid(spec_1, spec_2)
results_metric_1, results_metric_2, results_metric_3 = [], [], []

#%% Metrics

product = KSA_product
product_chemical_IDs = ['PotassiumSorbate']
per_kg_KSA_to_per_kg_SA = TAL_chemicals.PotassiumSorbate.MW/TAL_chemicals.SorbicAcid.MW
get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass

get_product_MPSP = lambda: TAL_tea.solve_price(product) * per_kg_KSA_to_per_kg_SA / get_product_purity() # USD / pure kg SA-eq.
get_product_GWP = lambda: TAL_lca.GWP * per_kg_KSA_to_per_kg_SA
get_product_FEC = lambda: TAL_lca.FEC * per_kg_KSA_to_per_kg_SA

get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs]) * TAL_tea.operating_hours
get_TAL_VOC = lambda: TAL_tea.VOC / 1e6 # million USD / yr
get_TAL_FCI = lambda: TAL_tea.FCI / 1e6 # million USD

TAL_metrics = [get_product_MPSP, get_product_GWP, get_product_FEC]

#%% Plot stuff

# Parameters analyzed across

x_label = r"$\bfDehydration$" + " " + r"$\bfPSA$" + " " + r"$\bfYield$" # title of the x axis
x_units = r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$"
x_ticks=np.arange(0.3, 1.1, 0.1)

y_label = r"$\bfDehydration$" + " " + r"$\bfReaction$" + " " + r"$\bfTime$" # title of the y axis
y_units = r"$\mathrm{h}$"
y_ticks=np.arange(0., 22., 2)


z_label = r"$\bfDehydration$" + " " + r"$\bfCatalyst$" + " " + r"$\bfPrice$" # title of the z axis
z_units =  r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
z_ticks=np.arange(0., 1100, 100)

# Metrics
MPSP_w_label = r"$\bfMPSP$" # title of the color axis
MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"

GWP_w_label = r"$\mathrm{\bfGWP}_{\bf100}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"

FEC_w_label = r"$\bfFEC$" # title of the color axis
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"


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

# %% Run TRY analysis 
for p in dehydration_cat_prices:
    data_1 = TAL_data = spec.evaluate_across_specs(
            TAL_sys, spec_1, spec_2, TAL_metrics, [p])
    
    # %% Save generated data
    
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    file_to_save = f'_{steps}_steps_'+'TAL_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    np.save(TAL_results_filepath+file_to_save, data_1)
    
    pd.DataFrame(data_1[:, :, 0, :][:,:,0]).to_csv(TAL_results_filepath+'MPSP-'+file_to_save+'.csv')
    pd.DataFrame(data_1[:, :, 1, :][:,:,0]).to_csv(TAL_results_filepath+'GWP-'+file_to_save+'.csv')
    pd.DataFrame(data_1[:, :, 2, :][:,:,0]*TAL_tea.operating_hours).to_csv(TAL_results_filepath+'FEC-'+file_to_save+'.csv')
    
    
    # %% Load previously saved data
    file_to_load = TAL_results_filepath+file_to_save
    data_1 = np.load(file_to_load+'.npy')
    
    # data_1_copy = copy.deepcopy(data_1)

    results_metric_1.append(data_1[:, :, 0, :][:,:,0])
    results_metric_2.append(data_1[:, :, 1, :][:,:,0])
    results_metric_3.append(data_1[:, :, 2, :][:,:,0])

#%% Plot metrics vs titer, yield, and productivity

import contourplots
chdir(TAL_results_filepath)

results_metric_1 = np.array(results_metric_1)
results_metric_2 = np.array(results_metric_2)
results_metric_3 = np.array(results_metric_3)

#%% More plot stuff

fps = 6
axis_title_fonts={'size': {'x': 8, 'y':8, 'z':8, 'w':8},}
clabel_fontsize = 8.
default_fontsize = 8.
keep_frames = False

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
    median, stdev = get_median(metric_data), metric_data.std()
    bound_diff = n_stdevs_for_bounds*stdev
    ub_temp = max(abs(median-bound_diff), abs(median+bound_diff))
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

#%% MPSP

MPSP_w_levels, MPSP_w_ticks, MPSP_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=0)


contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_1, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*PSA_yields, # x axis values
                                y_data=dehydration_times, # y axis values
                                z_data=dehydration_cat_prices, # z axis values
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
                                fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
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
                                # comparison_range=[6.5, 7.5],
                                # comparison_range=[MPSP_w_levels[-2], MPSP_w_levels[-1]],
                                # comparison_range_hatch_pattern='////',
                                )

#%% GWP

GWP_w_levels, GWP_w_ticks, GWP_cbar_ticks = get_contour_info_from_metric_data(results_metric_2,)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_2, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*PSA_yields, # x axis values
                                y_data=dehydration_times, # y axis values
                                z_data=dehydration_cat_prices, # z axis values
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
                                fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                cbar_ticks=GWP_cbar_ticks,
                                z_marker_color='g', # default matplotlib color namesaxis_title_fonts={'size': {'x': 8, 'y':8, 'z':8, 'w':8},},
                                axis_title_fonts=axis_title_fonts,
                                clabel_fontsize = clabel_fontsize,
                                default_fontsize = default_fontsize,
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='GWP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                )


#%% FEC

FEC_w_levels, FEC_w_ticks, FEC_cbar_ticks = get_contour_info_from_metric_data(results_metric_3,)

contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_3, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*PSA_yields, # x axis values
                                y_data=dehydration_times, # y axis values
                                z_data=dehydration_cat_prices, # z axis values
                                x_label=x_label, # title of the x axis
                                y_label=y_label, # title of the y axis
                                z_label=z_label, # title of the z axis
                                w_label=MPSP_w_label, # title of the color axis
                                x_ticks=100*x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=FEC_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=FEC_w_ticks, # labeled, lined contours; a subset of w_levels
                                x_units=x_units,
                                y_units=y_units,
                                z_units=z_units,
                                w_units=FEC_units,
                                fmt_clabel=lambda cvalue: "{:.0f}".format(cvalue), # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cbar_ticks=FEC_cbar_ticks,
                                z_marker_color='g', # default matplotlib color names
                                axis_title_fonts=axis_title_fonts,
                                clabel_fontsize = clabel_fontsize,
                                default_fontsize = default_fontsize,
                                fps=fps, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='FEC_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                )

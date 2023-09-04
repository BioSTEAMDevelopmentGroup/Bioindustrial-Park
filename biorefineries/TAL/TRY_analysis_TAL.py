# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: sarangbhagwat
"""

# Run this cell first
from warnings import filterwarnings
import copy
filterwarnings('ignore')
# from biosteam.process_tools import UnitGroup
from biosteam.digraph import digraph_from_units, save_digraph
from biosteam.utils import streams_from_units, filter_out_missing_streams, colors
from biosteam.report import unit_result_tables, tables_to_excel, stream_table
import numpy as np
import biosteam as bst
# from biosteam.utils import colors
# from matplotlib.ticker import AutoMinorLocator as AML

# from TAL.system_solubility_exploit import TAL_sys, TAL_tea, R302, spec
from biorefineries import TAL
from biorefineries.TAL.system_TAL_solubility_exploit_ethanol_sugarcane import TAL_sys, TAL_tea, TAL_lca, R302, spec, TAL_product, simulate_and_print
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
# bst.speed_up()

product = TAL_product

#%% Filepath
TAL_filepath = TAL.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\TAL\\analyses\\results
# chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'

# %%Colors

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    # CABBI_colors = (colors.CABBI_yellow.tint(50).RGBn,
    #                 colors.CABBI_yellow.shade(10).RGBn,
    #                 colors.CABBI_green.RGBn,
    #                 colors.CABBI_green.shade(30).RGBn)

    CABBI_colors = (colors.CABBI_orange.RGBn,
                    colors.CABBI_yellow.RGBn,

                    colors.CABBI_green.RGBn,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)


def CABBI_grey_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_grey.RGBn,
                    colors.grey_shade.RGBn,
                    colors.grey_dark.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

def CABBI_blue_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_blue_light.RGBn,
                    colors.CABBI_teal.RGBn,
                    colors.CABBI_teal_green.shade(25).RGBn,
                    colors.CABBI_green_dirty.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

#%% 
SA_price_range = [6500, 7500]

product_chemical_IDs = ['TAL',]
get_product_MPSP = lambda: TAL_tea.solve_price(product) * 907.185 / get_product_purity() # USD / ton
get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])
# get_product_MPSP = lambda: TAL_tea.solve_price(product) / get_product_purity() # USD / kg
get_TAL_VOC = lambda: TAL_tea.VOC / 1e6 # million USD / yr
get_TAL_FCI = lambda: TAL_tea.FCI / 1e6 # million USD


get_TAL_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol

get_TAL_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol

# get_rel_impact_t_y = lambda: rel_impact_fn(steps)

TAL_metrics = [get_product_MPSP, lambda: TAL_lca.GWP, lambda: TAL_lca.FEC]
# TAL_metrics = [get_TAL_MPSP, get_GWP, get_FEC]

# %% Generate 3-specification meshgrid and set specification loading functions
steps = 20

# Yield, titer, productivity (rate)
spec_1 = yields = np.linspace(0.1, 0.5, steps) # yield
spec_2 = titers = np.linspace(10., 50., steps) # titer


# spec_3 = productivities =\
#     np.array([0.2*spec.baseline_productivity, spec.baseline_productivity, 5.*spec.baseline_productivity])
    
spec_3 = productivities =\
    np.array([spec.baseline_productivity,])


#%%
spec_1, spec_2 = np.meshgrid(spec_1, spec_2)

#%%
# simulate_and_print()
spec.load_specifications(yields[0], titers[0], productivities[0])
spec.set_production_capacity()
simulate_and_print()

# %% Run TRY analysis
system = TAL_sys
HXN = spec.HXN
results_metric_1, results_metric_2, results_metric_3 = [], [], []

def print_status(curr_no, total_no, s1, s2, s3, HXN_qbal_error, results=None, exception_str=None,):
    print('\n\n')
    print(f'{curr_no}/{total_no}')
    print('\n')
    print(s1, s2, s3)
    print('\n')
    print(f'HXN qbal error = {round(HXN_qbal_error, 2)} %.')
    print('\n')
    print(results)

max_HXN_qbal_percent_error = 0.

curr_no = 0
total_no = len(yields)*len(titers)*len(productivities)

for p in productivities:
    # data_1 = TAL_data = spec.evaluate_across_specs(
    #         TAL_sys, spec_1, spec_2, TAL_metrics, [p])
    
    d1_Metric1, d1_Metric2, d1_Metric3 = [], [], []
    for y in yields:
        d1_Metric1.append([])
        d1_Metric2.append([])
        d1_Metric3.append([])
        for t in titers:
            curr_no +=1
            error_message = None
            try:
                spec.load_specifications(spec_1=y, spec_2=t, spec_3=p)
                spec.set_production_capacity()
                # system.simulate()
                d1_Metric1[-1].append(TAL_metrics[0]())
                d1_Metric2[-1].append(TAL_metrics[1]())
                d1_Metric3[-1].append(TAL_metrics[2]())
                
                HXN_qbal_error = HXN.energy_balance_percent_error
                if max_HXN_qbal_percent_error<HXN_qbal_error: max_HXN_qbal_percent_error = HXN_qbal_error
            
            except RuntimeError as e1:
                d1_Metric1[-1].append(np.nan)
                d1_Metric2[-1].append(np.nan)
                d1_Metric3[-1].append(np.nan)
                error_message = str(e1)
            
            except ValueError as e1:
                d1_Metric1[-1].append(np.nan)
                d1_Metric2[-1].append(np.nan)
                d1_Metric3[-1].append(np.nan)
                error_message = str(e1)
            
            print_status(curr_no, total_no,
                         y, t, p, 
                         results=[d1_Metric1[-1][-1], d1_Metric2[-1][-1], d1_Metric3[-1][-1],],
                         HXN_qbal_error=HXN.energy_balance_percent_error,
                         exception_str=error_message)
    
    d1_Metric1, d1_Metric2, d1_Metric3 = np.array(d1_Metric1), np.array(d1_Metric2), np.array(d1_Metric3)
    d1_Metric1, d1_Metric2, d1_Metric3 = d1_Metric1.transpose(), d1_Metric2.transpose(), d1_Metric3.transpose()
    
    results_metric_1.append(d1_Metric1/907.185)
    results_metric_2.append(d1_Metric2)
    results_metric_3.append(d1_Metric3)

print(f'Max HXN Q bal error was {round(max_HXN_qbal_percent_error, 3)} %.')

    # %% Save generated data
    
    
    # minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    # file_to_save = f'_{steps}_steps_'+'TAL_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    # np.save(TAL_results_filepath+file_to_save, data_1)
    
    # pd.DataFrame(data_1[:, :, 0, :][:,:,0]/907.185).to_csv(TAL_results_filepath+'MPSP-'+file_to_save+'.csv')
    # pd.DataFrame(data_1[:, :, 1, :][:,:,0]).to_csv(TAL_results_filepath+'GWP-'+file_to_save+'.csv')
    # pd.DataFrame(data_1[:, :, 2, :][:,:,0]*TAL_tea.operating_hours).to_csv(TAL_results_filepath+'FEC-'+file_to_save+'.csv')
    
    
    # %% Load previously saved data
    # file_to_load = TAL_results_filepath+file_to_save
    # # file_to_load = 'C:/Users/saran/Documents/Academia/Spring 2020/BioSTEAM/Bioindustrial-Park/BioSTEAM 2.x.x/biorefineries/TAL/TAL_TRY_2020.9.22-16.44'
    # data_1 = np.load(file_to_load+'.npy')
    # # data_1 = np.load('TAL_TRY_2022.3.18-22.15.npy')
    # data_1_copy = copy.deepcopy(data_1)
    
    # data_2 = data_1
    
    # d1_Metric1 = data_1[:, :, 0, :]
    # d1_Metric2 = data_1[:, :, 1, :]
    # d1_Metric3 = data_1[:, :, 2, :]
    
    # d2_Metric1 = data_2[:, :, 0, :]
    # d2_Metric2 = data_2[:, :, 1, :]
    # d2_Metric3 = data_2[:, :, 2, :]
    
    
    # results_metric_1.append(d1_Metric1[:,:,0]/907.185)
    # results_metric_2.append(d1_Metric2[:,:,0])
    # results_metric_3.append(d2_Metric3[:,:,0])

#%% Plot metrics vs titer, yield, and productivity

import contourplots

chdir(TAL_results_filepath)

MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"

# yields = np.linspace(0.4, 0.9, steps) # x axis values
# titers = np.linspace(40., 120., steps) # y axis values
# productivities = np.arange(0.1, 2.1, 0.1) # z axis values
# MPSPs = results_metric_1 # too big to show here; shape = z * x * y # color (or w) axis values


x_ticks=[10, 20, 30, 40, 50]
y_ticks=[10, 20, 30, 40, 50]
z_ticks=[0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.30]
                                
keep_frames=True

minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = f'_{steps}_steps_'+'TAL_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)

#%% MPSP
# MPSP_w_levels = np.array([0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25.])
MPSP_w_levels = np.arange(0., 15.5, 0.5)
contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_1, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*yields, # x axis values
                                # x_data = yields/theoretical_max_g_TAL_acid_per_g_glucose,
                                y_data=titers, # y axis values
                                z_data=productivities, # z axis values
                                x_label=r"$\bfYield$", # title of the x axis
                                y_label=r"$\bfTiter$", # title of the y axis
                                z_label=r"$\bfProductivity$", # title of the z axis
                                w_label=r"$\bfMPSP$", # title of the color axis
                                x_ticks=x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=MPSP_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([1., 2., 3., 4., 5., 6., 8., 10., ]), # labeled, lined contours; a subset of w_levels
                                # w_ticks=[],
                                x_units=r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$",
                                # x_units=r"$\mathrm{g} \cdot \mathrm{g}$" + " " + r"$\mathrm{glucose}$" + " " + r"$\mathrm{eq.}^{-1}$",
                                y_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}$",
                                z_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$",
                                w_units=MPSP_units,
                                # fmt_clabel=lambda cvalue: "{:.1f}".format(cvalue), # format of contour labels
                                fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                extend_cmap='max',
                                # cbar_ticks=[0., 5., 10., 15., 20., 25.,],
                                cbar_ticks=np.array([0., 2.5, 5., 7.5, 10., 12.5, 15., 17.5]),
                                z_marker_color='g', # default matplotlib color names
                                axis_title_fonts={'size': {'x': 14, 'y':14, 'z':14, 'w':14},},
                                fps=3, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='part2_MPSP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                # comparison_range=[6.51 * 112.12652/126.11004, 7.43* 112.12652/126.11004], # maximum allowable TAL price = market price of SA ($/kg-SA) * molar mass of SA / molar mass of TAL
                                # comparison_range_hatch_pattern='////',
                                )

#%% GWP
contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_2, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*yields, # x axis values
                                y_data=titers, # y axis values
                                z_data=productivities, # z axis values
                                x_label=r"$\bfYield$", # title of the x axis
                                y_label=r"$\bfTiter$", # title of the y axis
                                z_label=r"$\bfProductivity$", # title of the z axis
                                w_label=r"$\mathrm{\bfGWP}_{\bf100}$", # title of the color axis
                                x_ticks=x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=np.arange(0., 5.25, 0.25), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([2, 3, 5,]), # labeled, lined contours; a subset of w_levels
                                x_units=r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$",
                                # x_units=r"$\mathrm{g} \cdot \mathrm{g}$" + " " + r"$\mathrm{glucose}$" + " " + r"$\mathrm{eq.}^{-1}$",
                                y_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}$",
                                z_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$",
                                w_units=GWP_units,
                                fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cbar_ticks=np.arange(0.0, 6., 1.),
                                z_marker_color='g', # default matplotlib color names
                                axis_title_fonts={'size': {'x': 14, 'y':14, 'z':14, 'w':14},},
                                fps=3, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='GWP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                )


#%% FEC
contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_3, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*yields, # x axis values
                                y_data=titers, # y axis values
                                z_data=productivities, # z axis values
                                x_label=r"$\bfYield$", # title of the x axis
                                y_label=r"$\bfTiter$", # title of the y axis
                                z_label=r"$\bfProductivity$", # title of the z axis
                                w_label=r"$\bfFEC$", # title of the color axis
                                x_ticks=x_ticks,
                                y_ticks=y_ticks,
                                z_ticks=z_ticks,
                                w_levels=np.arange(-30, 100, 10), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([-20, -10, 0, 10, 20, 50, 90]), # labeled, lined contours; a subset of w_levels
                                x_units=r"$\mathrm{\%}$" + " " + r"$\mathrm{theoretical}$",
                                # x_units=r"$\mathrm{g} \cdot \mathrm{g}$" + " " + r"$\mathrm{glucose}$" + " " + r"$\mathrm{eq.}^{-1}$",
                                y_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}$",
                                z_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$",
                                w_units=FEC_units,
                                fmt_clabel=lambda cvalue: "{:.0f}".format(cvalue), # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cbar_ticks=np.arange(-30, 100, 30),
                                z_marker_color='g', # default matplotlib color names
                                axis_title_fonts={'size': {'x': 14, 'y':14, 'z':14, 'w':14},},
                                fps=3, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='FEC_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                )
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
from matplotlib.ticker import AutoMinorLocator as AML

# from biorefineries.succinic.system_sc import succinic_sys, succinic_tea, R302, spec, product_stream

from biorefineries import succinic
from biorefineries.succinic.models import model, succinic_sys, succinic_tea, succinic_LCA, R302, spec, product_stream, theoretical_max_g_succinic_acid_per_g_glucose
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
ig = np.seterr(invalid='ignore')
# bst.speed_up()
import os

chdir = os.chdir

product = product_stream

modes = ['lab_batch', 'lab_fed-batch', 'pilot_batch']

parameter_distributions_filenames = ['parameter-distributions_lab-scale_batch.xlsx',
                                    'parameter-distributions_lab-scale_fed-batch.xlsx',
                                    'parameter-distributions_pilot-scale_batch.xlsx',
                                    ]


dateTimeObj = datetime.now()

#%% Set mode
# ## Change working directory to biorefineries\\succinic
# chdir(succinic.__file__.replace('\\__init__.py', ''))
# ##
mode = 'pilot_batch'

succinic_filepath = succinic.__file__.replace('\\__init__.py', '')
parameter_distributions_filename = succinic_filepath+'\\analyses\\parameter_distributions\\'+parameter_distributions_filenames[modes.index(mode)]
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename)
model.metrics_at_baseline()

# ## Change working directory to biorefineries\\succinic\\analyses\\results
# chdir(succinic.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
succinic_results_filepath = succinic_filepath + '\\analyses\\results\\'
# %%Colors
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

# def CABBI_green_colormap(N_levels=90):
#     """
#     Return a matplotlib.colors.LinearSegmentedColormap object
#     that serves as CABBI's green colormap theme for contour plots.
    
#     """
#     # CABBI_colors = (colors.CABBI_yellow.tint(50).RGBn,
#     #                 colors.CABBI_yellow.shade(10).RGBn,
#     #                 colors.CABBI_green.RGBn,
#     #                 colors.CABBI_green.shade(30).RGBn)

#     CABBI_colors = (colors.CABBI_yellow.RGBn,
#                     colors.CABBI_green.RGBn,
#                     colors.CABBI_teal_green.shade(90).RGBn)
#     return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

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
# def plot_contour_2d(X_grid, Y_grid, Z_1d, data, 
#                     xlabel, ylabel, xticks, yticks, 
#                     metric_bars, Z_label=None,
#                     Z_value_format=lambda Z: str(Z),
#                     fillblack=True):
#     """Create contour plots and return the figure and the axes."""
#     nrows = len(metric_bars)
#     ncols = len(Z_1d)
#     assert data.shape == (*X_grid.shape, nrows, ncols), (
#         "data shape must be (X, Y, M, Z), where (X, Y) is the shape of both X_grid and Y_grid, "
#         "M is the number of metrics, and Z is the number of elements in Z_1d"
#     )
#     widths = np.ones(ncols + 1)
#     widths[-1] /= 4
#     gs_kw = dict(width_ratios=widths)
#     fig, axes = plt.subplots(ncols=ncols + 1, nrows=nrows, gridspec_kw=gs_kw)
#     row_counter = 0
#     for row in range(nrows):
#         if row  == 2:
#             metric_bar = metric_bars[row]
#             for col in range(ncols):
#                 ax = axes[row, col]
#                 plt.sca(ax)
#                 style_plot_limits(xticks, yticks)
#                 yticklabels = col == 0
#                 xticklabels = row == nrows - 1
#                 if fillblack: fill_plot()
#             #     cp = plt.contourf(X_grid, Y_grid, data[:, :, row, col],
#             #                       levels=metric_bar.levels,
#             #                       cmap=metric_bar.cmap)
#             cbar_ax = axes[row, -1]
#             # style_axis(ax, xticks, yticks, xticklabels, yticklabels)
#             # # Only if you want a log-scale colorbar
            
#             # pcm = ax.pcolormesh(X_grid, Y_grid, data[:, :, row, col],
#             #            norm=mcolors.SymLogNorm(linthresh=0.003, linscale=1,
#             #                                   vmin=-4, vmax=4, base=10),
#             #            cmap=metric_bar.cmap)
#             pcm = plt.contourf(X_grid, Y_grid, data[:, :, row, col],
#                        norm=mcolors.SymLogNorm(linthresh=0.003, linscale=1,
#                             levels=metric_bar.levels, vmin=-4, vmax=4, base=10),
#                        cmap=metric_bar.cmap)
#             fig.colorbar(pcm, ax=cbar_ax, shrink = 0.8)
#             # else:
#             #     metric_bar.colorbar(fig, cbar_ax, colorplot = cp, shrink=0.8)
#         # plt.clim()
#     for col in range(ncols):
#         if not col and Z_label:
#             title = f"{Z_label}: {Z_value_format(Z_1d[col])}"
#         else:
#             title = Z_value_format(Z_1d[col])
#         ax = axes[0, col]
#         ax.set_title(title)
#     for ax in axes[:, -1]:
#         plt.sca(ax)
#         plt.axis('off')
#     set_axes_labels(axes[:, :-1], xlabel, ylabel)
#     plt.subplots_adjust(hspace=0.1, wspace=0.1)
#     return fig, axes

million_dollar = r"\mathrm{MM\$}"
MPSP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"

# VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
# FCI_units = f"${million_dollar}$"
VOC_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
# VOC_units = r"$\mathrm{kg CO2 eq.} \cdot \mathrm{kg MEK}^{-1}$"
# FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
# FCI_units = r"$\mathrm{MJ} \cdot \mathrm{kg MEK}^{-1}$"

def tickmarks_from_data(data, accuracy=50, N_points=5):
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    dmin = floor(dmin/accuracy) * accuracy
    dmax = ceil(dmax/accuracy) * accuracy
    step = (dmax - dmin) / (N_points - 1)
    return [dmin + step * i for i in range(N_points)]

# areas = {
#     'Fermentation': (ol.T102, ol.P102, ol.M102, ol.H101, ol.R101, ol.T105, ol.P104),
#     '3-Phase Decanter': (ol.C101, ol.P107),
#     'Dehydration': (ol.H102, ol.M103, ol.H103, ol.R102, ol.T106),
#     'Separation': (ol.H105, ol.P105, ol.C102, ol.P106, ol.H104, ol.D101,
#                    ol.D102, ol.M104, ol.M105, ol.P108, ol.P109, ol.H106),
#     'OSBL': (ol.T101, ol.T104, ol.T103, ol.T108,
#              ol.CT, ol.BT, ol.CWP, ol.T107),
# }
# area_groups = [UnitGroup(i, j) for i,j in areas.items()]
# area_groups.append(UnitGroup('Total', sum(areas.values(), ())))

# Contour stuff

target_spec_1 = .90 # yield
target_spec_2 = 180 # titer
target_spec_3 = 1.5 # productivity
lab_spec_1 = .80
lab_spec_2 = 109.9
lab_spec_3 = 1

# target_spec_2 = 180
# target_spec_1 = 1451.496
# target_spec_3 = 7.126
# lab_spec_2 = 109.9
# lab_spec_1 = 30
# lab_spec_3 = 71.26

# succinic_price_range = [2.3 * 907.185, 3.2 * 907.185] 
SA_price_range = [6500, 7500]
# temporary price range from https://www.alibaba.com/product-detail/hot-sale-C4H8O-butanon-mek_62345760689.html?spm=a2700.7724857.normalList.26.1d194486SbCyfR

product_chemical_IDs = ['SuccinicAcid',]
# get_product_MPSP = lambda: succinic_tea.solve_price(product) * 907.185 / get_product_purity() # USD / ton

def get_product_MPSP():
    model.specification()
    return succinic_tea.solve_price(product) * 907.185 / get_product_purity() # USD / ton

get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])
# get_product_MPSP = lambda: succinic_tea.solve_price(product) / get_product_purity() # USD / kg
get_succinic_VOC = lambda: succinic_tea.VOC / 1e6 # million USD / yr
get_succinic_FCI = lambda: succinic_tea.FCI / 1e6 # million USD



# def rel_impact_fn(steps):
#     rel_impact_yield_titer = None
#     lower_bound = 1/10
#     upper_bound = 10
#     middle = 1
#     rel_impact_titer_yield = None
#     max_theoretical_yield = 1
#     max_theoretical_titer = 300
#     d_MPSP_d_yield, d_MPSP_d_titer = None, None
    
#     fermentor = spec.reactor
#     curr_MPSP = succinic_tea.solve_price(product) * 907.185
#     curr_yield = fermentor.glucose_to_succinic_rxn.X
#     curr_titer = fermentor.outs[0].imass['TAL']/fermentor.outs[0].F_vol
    
#     multi_for_yield = 1
#     # d_yield = multi_for_yield * rel_step * max_theoretical_yield
    
#     multi_for_titer = 1
#     # d_titer = multi_for_titer * rel_step * max_theoretical_titer
    
#     d_yield, d_titer = (0.99-0.30)/steps, (210-70)/steps
#     next_yield = curr_yield + d_yield
#     next_titer = curr_titer + d_titer
    
#     # if (next_yield, next_titer) in TYM_dict.keys():
#     #     return 1
#     try:
#         if (curr_yield + d_yield > max_theoretical_yield - 0.01
#             and not curr_titer + d_titer > max_theoretical_titer):
#             if steps > 80:
#                 rel_impact_yield_titer = upper_bound
#             else:
#                 return rel_impact_fn(1.1*steps)
#         elif (curr_titer + d_titer > max_theoretical_titer
#             and not curr_yield + d_yield > max_theoretical_yield - 0.01):
#             if steps > 80:
#                 rel_impact_yield_titer = lower_bound
#             else:
#                 return rel_impact_fn(1.1*steps)
#         elif (curr_yield + d_yield > max_theoretical_yield - 0.01
#             and curr_titer + d_titer > max_theoretical_titer):
#             rel_impact_yield_titer = middle
#         else:
#             # d_yield = multi_for_yield * rel_step * max_theoretical_yield
#             spec.load_yield(curr_yield + d_yield)
#             succinic_sys.simulate()
#             MPSP_d_yield = succinic_tea.solve_price(product) * 907.185
#             d_MPSP_d_yield = (MPSP_d_yield - curr_MPSP) * multi_for_yield
#             spec.load_yield(curr_yield)
#             # d_titer = multi_for_titer * rel_step * max_theoretical_titer
#             spec.load_titer(curr_titer + d_titer)
#             succinic_sys.simulate()
#             MPSP_d_titer = succinic_tea.solve_price(product) * 907.185
#             d_MPSP_d_titer = (MPSP_d_titer - curr_MPSP) * multi_for_titer
            
#             rel_impact_yield_titer = d_MPSP_d_yield / d_MPSP_d_titer
#             rel_impact_titer_yield = 1/rel_impact_yield_titer
#             spec.load_titer(curr_titer)
#     except ValueError: # inadequate yield for given titer
#         rel_impact_titer_yield = lower_bound
#     if (rel_impact_titer_yield == None or rel_impact_titer_yield == np.nan
#         or rel_impact_titer_yield<0):
#         return rel_impact_fn(1.1*steps)
#     print(curr_yield, curr_titer, d_MPSP_d_yield, d_MPSP_d_titer, rel_impact_titer_yield)
#     # return 1 + log(rel_impact_yield_titer, 10)
#     return rel_impact_titer_yield

get_succinic_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol

get_succinic_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol

get_product_GWP = lambda: succinic_LCA.GWP
get_product_FEC = lambda: succinic_LCA.FEC
# get_rel_impact_t_y = lambda: rel_impact_fn(steps)

succinic_metrics = [get_product_MPSP, get_product_GWP, get_product_FEC]
# succinic_metrics = [get_succinic_MPSP, get_GWP, get_FEC]

# %% Generate 3-specification meshgrid and set specification loading functions
steps = 50

# Neutralization
spec.neutralization = False

# Yield, titer, productivity (rate)
spec_1 = yields = np.linspace(0.2, 0.8, steps) # yield
spec_2 = titers = np.linspace(20., 120., steps) # titer
# spec_3 = productivities = np.linspace(0.1, 1.5, 6) # productivity
spec_3 = productivities = [spec.baseline_productivity]

# spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer
# spec.load_spec_3 = spec.load_productivity
xlabel = "Yield"
ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
xticks = [0.1, 0.3, 0.5, 0.7, 0.9]
yticks = [0, 5, 10, 15, 20, 25, 30]
spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"

spec_1, spec_2 = np.meshgrid(spec_1, spec_2)

results_metric_1, results_metric_2, results_metric_3 = [], [], []

# %% Run TRY analysis 
for p in productivities:
    data_1 = succinic_data = spec.evaluate_across_specs(
            succinic_sys, spec_1, spec_2, succinic_metrics, [p])
    
    
    # %% Save generated data
    
    
    minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    file_to_save = mode+f'_{steps}_steps_'+f'neutralization_{spec.neutralization}_'+ 'succinic_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    np.save(succinic_results_filepath+file_to_save, data_1)
    
    pd.DataFrame(data_1[:, :, 0, :][:,:,0]/907.185).to_csv(succinic_results_filepath+'MPSP-'+file_to_save+'.csv')
    pd.DataFrame(data_1[:, :, 1, :][:,:,0]).to_csv(succinic_results_filepath+'GWP-'+file_to_save+'.csv')
    pd.DataFrame(data_1[:, :, 2, :][:,:,0]*succinic_tea.operating_hours).to_csv(succinic_results_filepath+'FEC-'+file_to_save+'.csv')
    
    
    # %% Load previously saved data
    file_to_load = succinic_results_filepath+file_to_save
    # file_to_load = 'C:/Users/saran/Documents/Academia/Spring 2020/BioSTEAM/Bioindustrial-Park/BioSTEAM 2.x.x/biorefineries/TAL/succinic_TRY_2020.9.22-16.44'
    data_1 = np.load(file_to_load+'.npy')
    # data_1 = np.load('succinic_TRY_2022.3.18-22.15.npy')
    data_1_copy = copy.deepcopy(data_1)
    
    data_2 = data_1
    
    d1_Metric1 = data_1[:, :, 0, :]
    d1_Metric2 = data_1[:, :, 1, :]
    d1_Metric3 = data_1[:, :, 2, :]
    
    d2_Metric1 = data_2[:, :, 0, :]
    d2_Metric2 = data_2[:, :, 1, :]
    d2_Metric3 = data_2[:, :, 2, :]
    
    # %% Functions to make regions with total sugars > 150 g/L and total inhibitors > 1000 mg/L
    # also infeasible
    
    
    def make_oversaccharine_region_infeasible():
        # mask = d1_Metric2>150.
        # data_1_copy[mask] = np.nan
        # # d1_Metric1[mask] = np.nan
        # # d1_Metric2[mask] = np.nan
        # # d1_Metric3[mask] = np.nan
        # # d2_Metric1[mask] = np.nan
        # # d2_Metric2[mask] = np.nan
        # # d2_Metric3[mask] = np.nan
        infeas = np.where(d1_Metric2>150.)
        data_1_copy[infeas] = np.nan
    # make_oversaccharine_region_infeasible()
    
    
    def make_inhibited_region_infeasible():
        infeas = np.where(d1_Metric3>1000.)
        data_1_copy[infeas] = np.nan
        # d1_Metric1[infeas] = np.nan
        # d1_Metric2[infeas] = np.nan
        # d1_Metric3[infeas] = np.nan
        # d2_Metric1[infeas] = np.nan
        # d2_Metric2[infeas] = np.nan
        # d2_Metric3[infeas] = np.nan
        
    # make_inhibited_region_infeasible()
    # %% Plot contours2
    # data_2 = data_1
    MPSP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"
    
    # VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
    # FCI_units = f"${million_dollar}$"
    VOC_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
    # VOC_units = r"$\mathrm{kg CO2 eq.} \cdot \mathrm{kg MEK}^{-1}$"
    # FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
    FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
    # FCI_units = r"$\mathrm{MJ eq.} \cdot \mathrm{kg MEK}^{-1}$"
    # data_1_copy = copy.deepcopy(data_1)
    
    Metric_1_tickmarks = tickmarks(
        dmin = min(d1_Metric1[~np.isnan(d1_Metric1)].min(), d2_Metric1[~np.isnan(d2_Metric1)].min()),
        dmax = max(d1_Metric1[~np.isnan(d1_Metric1)].max(), d2_Metric1[~np.isnan(d2_Metric1)].max())
    )
    Metric_2_tickmarks = tickmarks(
        dmin = min(d1_Metric2[~np.isnan(d1_Metric2)].min(), d2_Metric2[~np.isnan(d2_Metric2)].min()),
        dmax = max(d1_Metric2[~np.isnan(d1_Metric2)].max(), d2_Metric2[~np.isnan(d2_Metric2)].max())
    )
    Metric_3_tickmarks = tickmarks(
        dmin = min(d1_Metric3[~np.isnan(d1_Metric3)].min(), d2_Metric3[~np.isnan(d2_Metric3)].min()),
        dmax = max(d1_Metric3[~np.isnan(d1_Metric3)].max(), d2_Metric3[~np.isnan(d2_Metric3)].max())
    )
    
    # Metric_3_tickmarks = [0.0*1000, 0.24*1000, 0.48*1000, 0.72*1000, 0.96*1000, 1.2*1000]
    
    # Metric_1_tickmarks = [3000, 4500, 6000, 7500, 9000, 10500, 12000, 13500]
    # Metric_1_tickmarks = [2000, 2500, 3000, 3500, 4000, 4500]
    # Metric_2_tickmarks = [4, 5, 6, 7, 8, 9, 10]
    Metric_2_tickmarks = [0, 50, 100, 150, 200, 250, 300]
    # Metric_3_tickmarks = [60, 70, 80, 90, 100, 110, 120]
    Metric_3_tickmarks = [0, 300, 600, 900, 1200, 1500]
    
    results_metric_1.append(data_1[:, :, 0, :][:,:,0]/907.185)
    results_metric_2.append(data_1[:, :, 1, :][:,:,0])
    results_metric_3.append(data_1[:, :, 2, :][:,:,0])

#%% Plot metrics vs titer, yield, and productivity

import contourplots

chdir(succinic_results_filepath)

MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
GWP_units = r"$\mathrm{kg}$"+" "+ r"$\mathrm{CO}_{2}\mathrm{-eq.}\cdot\mathrm{kg}^{-1}$"
FEC_units = r"$\mathrm{MJ}\cdot\mathrm{kg}^{-1}$"

# yields = np.linspace(0.4, 0.9, steps) # x axis values
# titers = np.linspace(40., 120., steps) # y axis values
# productivities = np.arange(0.1, 2.1, 0.1) # z axis values
# MPSPs = results_metric_1 # too big to show here; shape = z * x * y # color (or w) axis values

#%% MPSP
contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results_metric_1, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., MPSP
                                x_data=100*yields, # x axis values
                                # x_data = yields/theoretical_max_g_succinic_acid_per_g_glucose,
                                y_data=titers, # y axis values
                                z_data=productivities, # z axis values
                                x_label=r"$\bfYield$", # title of the x axis
                                y_label=r"$\bfTiter$", # title of the y axis
                                z_label=r"$\bfProductivity$", # title of the z axis
                                w_label=r"$\bfMPSP$", # title of the color axis
                                x_ticks=[20, 30, 40, 50, 60, 70, 80],
                                y_ticks=[20, 40, 60, 80, 100, 120],
                                z_ticks=np.arange(0.0, 2.5, 0.5),
                                w_levels=np.arange(0.8, 2.3, 0.1), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([1.0, 1.1, 1.2, 1.3,  1.5, 1.8, 2.2,]), # labeled, lined contours; a subset of w_levels
                                x_units=r"$\mathrm{\% theoretical}$",
                                # x_units=r"$\mathrm{g} \cdot \mathrm{g}$" + " " + r"$\mathrm{glucose}$" + " " + r"$\mathrm{eq.}^{-1}$",
                                y_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}$",
                                z_units=r"$\mathrm{g} \cdot \mathrm{L}^{-1}  \cdot \mathrm{h}^{-1}$",
                                w_units=MPSP_units,
                                fmt_clabel=lambda cvalue: "{:.2f}".format(cvalue), # format of contour labels
                                cmap=CABBI_green_colormap(), # can use 'viridis' or other default matplotlib colormaps
                                cbar_ticks=np.arange(0.8, 2.3, 0.2),
                                z_marker_color='g', # default matplotlib color names
                                axis_title_fonts={'size': {'x': 14, 'y':14, 'z':14, 'w':14},},
                                fps=3, # animation frames (z values traversed) per second
                                n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                animated_contourplot_filename='MPSP_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                keep_frames=False, # leaves frame PNG files undeleted after running; False by default
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
                                x_ticks=[20, 30, 40, 50, 60, 70, 80],
                                y_ticks=[20, 40, 60, 80, 100, 120],
                                z_ticks=np.arange(0.0, 2.5, 0.5),
                                w_levels=np.arange(0., 5.25, 0.25), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([1.5, 2, 3, 5,]), # labeled, lined contours; a subset of w_levels
                                x_units=r"$\mathrm{\% theoretical}$",
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
                                keep_frames=False, # leaves frame PNG files undeleted after running; False by default
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
                                x_ticks=[20, 30, 40, 50, 60, 70, 80],
                                y_ticks=[20, 40, 60, 80, 100, 120],
                                z_ticks=np.arange(0.0, 2.5, 0.5),
                                w_levels=np.arange(-30, 100, 10), # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                w_ticks=np.array([-20, -10, 0, 10, 20, 50, 90]), # labeled, lined contours; a subset of w_levels
                                x_units=r"$\mathrm{\% theoretical}$",
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
                                keep_frames=False, # leaves frame PNG files undeleted after running; False by default
                                )
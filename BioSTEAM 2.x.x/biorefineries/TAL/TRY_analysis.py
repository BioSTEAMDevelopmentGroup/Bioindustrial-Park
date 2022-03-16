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

# from TAL.system_solubility_exploit import TAL_sys, TAL_tea, R302, spec
# from TAL.system_solubility_exploit import MEK as product
from TAL.system_TAL_adsorption_glucose import TAL_sys, TAL_tea, R302, spec
# get_GWP, get_non_bio_GWP, get_FEC, get_SPED
from TAL.system_TAL_adsorption_glucose import SA as product

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
bst.speed_up()
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

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    # CABBI_colors = (colors.CABBI_yellow.tint(50).RGBn,
    #                 colors.CABBI_yellow.shade(10).RGBn,
    #                 colors.CABBI_green.RGBn,
    #                 colors.CABBI_green.shade(30).RGBn)

    CABBI_colors = (colors.CABBI_yellow.RGBn,
                    colors.CABBI_green.RGBn,
                    colors.CABBI_teal_green.shade(30).RGBn)
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

# TAL_price_range = [2.3 * 907.185, 3.2 * 907.185] 
SA_price_range = [6500, 7500]
# temporary price range from https://www.alibaba.com/product-detail/hot-sale-C4H8O-butanon-mek_62345760689.html?spm=a2700.7724857.normalList.26.1d194486SbCyfR

product_chemical_ID = 'TAL'
get_SA_MPSP = lambda: TAL_tea.solve_price(product) * 907.185 / get_product_purity() # USD / ton
get_product_purity = lambda: product.imass[product_chemical_ID]/product.F_mass
# get_SA_MPSP = lambda: TAL_tea.solve_price(product) / get_product_purity() # USD / kg
get_TAL_VOC = lambda: TAL_tea.VOC / 1e6 # million USD / yr
get_TAL_FCI = lambda: TAL_tea.FCI / 1e6 # million USD



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
#     curr_MPSP = TAL_tea.solve_price(product) * 907.185
#     curr_yield = fermentor.glucose_to_TAL_rxn.X
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
#             TAL_sys.simulate()
#             MPSP_d_yield = TAL_tea.solve_price(product) * 907.185
#             d_MPSP_d_yield = (MPSP_d_yield - curr_MPSP) * multi_for_yield
#             spec.load_yield(curr_yield)
#             # d_titer = multi_for_titer * rel_step * max_theoretical_titer
#             spec.load_titer(curr_titer + d_titer)
#             TAL_sys.simulate()
#             MPSP_d_titer = TAL_tea.solve_price(product) * 907.185
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

get_TAL_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol

get_TAL_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol

# get_rel_impact_t_y = lambda: rel_impact_fn(steps)

TAL_metrics = [get_SA_MPSP, get_TAL_sugars_conc, get_TAL_inhibitors_conc]
# TAL_metrics = [get_TAL_MPSP, get_GWP, get_FEC]

# %% Generate 3-specification meshgrid and set specification loading functions
steps = 50

# Yield, titer, productivity (rate)
spec_1 = np.linspace(0.1, 0.9, steps) # yield
spec_2 = np.linspace(1., 30., steps) # titer
# spec_1 = np.linspace(0.2, 0.99, steps) # yield
# spec_2 = np.linspace(45, 225, steps) # titer
spec_3 = np.array([0.21,]) # productivity
spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity
xlabel = "Yield"
ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
# xticks = [0.33, 0.66, 0.99]
xticks = [0.1, 0.3, 0.5, 0.7, 0.9]
# yticks = [75, 150, 225]
yticks = [0, 5, 10, 15, 20, 25, 30]
# xticks = [0.2, 0.6, 0.99]
# yticks = [45, 135, 225]
spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"

# # Yield, titer, productivity (rate)
# spec_1 = np.linspace(0.70, 0.90, 10) # yield
# spec_2 = np.linspace(90, 130, 10) # titer
# spec_3 = np.array([1]) # productivity
# spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer
# spec.load_spec_3 = spec.load_productivity
# xlabel = "Yield"
# ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
# xticks = [.75, .80, .85, .90]
# yticks = [90, 100, 110, 120, 130]
# spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"

# # Dehydration conversion, titer, feedstock price
# spec_1 = np.linspace(0.50, 0.80, 3)
# spec_2 = np.linspace(40, 200, 3)
# spec_3 = np.array([7.126, 71.26, 712.26])
# spec.load_spec_1 = spec.load_dehydration_conversion
# spec.load_spec_2 = spec.load_titer
# spec.load_spec_3 = spec.load_feedstock_price
# xlabel = "Dehydration conversion [%]"
# ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
# xticks = [.50, .60, .70, .80]
# yticks = [40, 80, 120, 160, 200]
# spec_3_units = "$\mathrm{\$} \cdot \mathrm{dry-ton}^{-1}$"

# # Byproducts (Acetoin and IBA) selling price, titer, feedstock price
# spec_1 = np.linspace(1, 2000, 10)
# spec_2 = np.linspace(40, 200, 10)
# spec_3 = np.array([7.126, 71.26, 712.26])
# spec.load_spec_1 = spec.load_byproducts_price
# spec.load_spec_2 = spec.load_titer
# spec.load_spec_3 = spec.load_feedstock_price
# xlabel = "Acetoin and IBA selling price [$\mathrm{\$} \cdot \mathrm{ton}^{-1}$]"
# ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
# xticks = [0, 400, 800, 1200, 1600, 2000]
# yticks = [40, 80, 120, 160, 200]
# spec_3_units = "$\mathrm{\$} \cdot \mathrm{dry-ton}^{-1}$"

spec_1, spec_2 = np.meshgrid(spec_1, spec_2)

# titers_2 = np.linspace(70, 190, 5)
# yields_2 = np.linspace(0.70, 0.90, 5)
# productivities_2 = np.array([0.5, 1, 1.5])
# titers_2, yields_2 = np.meshgrid(titers_2, yields_2)

# utilizes_xylose = False


# %% Run TRY analysis 

data_1 = TAL_data = spec.evaluate_across_specs(
        TAL_sys, spec_1, spec_2, TAL_metrics, spec_3)

# spec.load_spec_1 = spec.load_dehydration_conversion
# spec.load_spec_2 = spec.load_titer
# spec.load_spec_3 = spec.load_feedstock_price
# data_1 = TAL_data = spec.evaluate_across_specs(
#         TAL_sys, spec_1, spec_2, TAL_metrics, spec_3)
# utilizes_xylose = True
# data_2 = TAL_data = spec.evaluate_across_TRY(
#         TAL_sys, titers_2, yields_2, TAL_metrics, productivities_2)




# %% Save generated data


dateTimeObj = datetime.now()
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
file_to_save = 'TAL_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
np.save(file_to_save, data_1)

pd.DataFrame(data_1[:, :, 0, :][:,:,0]).to_csv('MPSP-'+file_to_save+'.csv')
pd.DataFrame(data_1[:, :, 1, :][:,:,0]).to_csv('GWP-'+file_to_save+'.csv')
pd.DataFrame(data_1[:, :, 2, :][:,:,0]).to_csv('FEC-'+file_to_save+'.csv')


# %% Load previously saved data
file_to_load = file_to_save
# file_to_load = 'C:/Users/saran/Documents/Academia/Spring 2020/BioSTEAM/Bioindustrial-Park/BioSTEAM 2.x.x/biorefineries/TAL/TAL_TRY_2020.9.22-16.44'
data_1 = np.load(file_to_load+'.npy')
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


# # %% Plot contours1
# # data_2 = data_1


# # data_1_copy = copy.deepcopy(data_1)

# Metric_1_tickmarks = tickmarks(
#     dmin = min(d1_Metric1[~np.isnan(d1_Metric1)].min(), d2_Metric1[~np.isnan(d2_Metric1)].min()),
#     dmax = max(d1_Metric1[~np.isnan(d1_Metric1)].max(), d2_Metric1[~np.isnan(d2_Metric1)].max())
# )
# Metric_2_tickmarks = tickmarks(
#     dmin = min(d1_Metric2[~np.isnan(d1_Metric2)].min(), d2_Metric2[~np.isnan(d2_Metric2)].min()),
#     dmax = max(d1_Metric2[~np.isnan(d1_Metric2)].max(), d2_Metric2[~np.isnan(d2_Metric2)].max())
# )
# Metric_3_tickmarks = tickmarks(
#     dmin = min(d1_Metric3[~np.isnan(d1_Metric3)].min(), d2_Metric3[~np.isnan(d2_Metric3)].min()),
#     dmax = max(d1_Metric3[~np.isnan(d1_Metric3)].max(), d2_Metric3[~np.isnan(d2_Metric3)].max())
# )

# # Metric_3_tickmarks = [0.0*1000, 0.24*1000, 0.48*1000, 0.72*1000, 0.96*1000, 1.2*1000]

# # Metric_1_tickmarks = [2100, 2800, 3500, 4200, 4900, 5600, 6300]
# Metric_1_tickmarks = [1500, 2000, 2500, 3000, 3500, 4000, 4500]
# # Metric_2_tickmarks = [0,1,2,3,4,5,6]
# Metric_2_tickmarks = [0, 100, 200, 300, 400]
# # Metric_3_tickmarks = [48, 72, 96, 120, 144]
# Metric_3_tickmarks = [200, 400, 600, 800, 1000, 1200]


# def plot(data, titers, yields, productivities, 
#          Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks):
#     metric_bars = (MetricBar('MPSP\n', MPSP_units, CABBI_green_colormap(),
#                              Metric_1_tickmarks,
#                              1 + int((max(Metric_1_tickmarks) - min(Metric_1_tickmarks))/250)),
#                    MetricBar('Total sugars\n', VOC_units, CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
#                              Metric_2_tickmarks,
#                              1 + int((max(Metric_2_tickmarks) - min(Metric_2_tickmarks))/25)),
#                    MetricBar('Total inhibitors\n', FCI_units,
#                              plt.cm.get_cmap('bone_r'), # plt.cm.get_cmap('bone_r'),
#                              Metric_3_tickmarks,
#                              1 + int((max(Metric_3_tickmarks) - min(Metric_3_tickmarks))/100)),)
    
#     return plot_contour_2d(titers, yields, productivities, data, 
#                                 xlabel, ylabel, xticks, yticks, metric_bars, 
#                                 Z_value_format=lambda Z: f"{Z:.1f} [{spec_3_units}]",
#                                 fillblack=False)






# metric_bars = (MetricBar('MPSP', MPSP_units, CABBI_green_colormap(),
#                          Metric_1_tickmarks, 800),
#                MetricBar('Total sugars', VOC_units, CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
#                          Metric_2_tickmarks, 80),
#                MetricBar("Relative impact on MPSP\n[impact of titer : impact of yield]", FCI_units,
#                          plt.cm.get_cmap('BrBG'), # plt.cm.get_cmap('bone_r'),
#                          Metric_3_tickmarks, 100))
    
    
# # if TAL_metrics[1] is get_TAL_sugars_conc:
# make_oversaccharine_region_infeasible()
# # # if TAL_metrics[2] is get_TAL_inhibitors_conc:
# make_inhibited_region_infeasible()

# fig, axes = plot(data_1_copy, spec_1, spec_2, spec_3, Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks)
# # spec_2 = 100 * spec_2
# index_feas_1 = d1_Metric2[:, :, 0]<150.
# index_feas_2 = d1_Metric3[:, :, 0]<1000.
# CS1, CS2, CS3, CS4 = 0, 0, 0, 0
# for i, ax_col in enumerate(axes[:, :len(spec_3)].transpose()):
#     MSP = data_1_copy[:, :, 0, i]
#     sugars_orig = data_1[:, :, 1, i]
#     inhibitors_orig = data_1[:, :, 2, i]
#     MSP_orig = data_1[:, :, 0, i]
#     j=0
#     for ax in ax_col:
        
#         plt.sca(ax)
#         # ax.set_facecolor('black')
#         # ax.patch.set_facecolor(CABBI_brown)
#         ax.xaxis.set_minor_locator(AML())
#         ax.yaxis.set_minor_locator(AML())
#         ax.tick_params(which='both', top=True, bottom=True, left=True, right=True)
#         ax.tick_params(which='minor', direction = 'inout', length = 4)
#         ax.tick_params(which='major', length = 10)
#         # ax.tick_params(which='both', width=2)
#         # ax.tick_params(which='major', length=7)
#         # ax.tick_params(which='minor', length=4, color='r')

#         if j==0: # MPSP plot only
#             Z_infeas_1 = MSP_orig * 0. + 1.
#             Z_infeas_2 = MSP_orig * 0. + 1.
#             Z_infeas_3 = MSP_orig * 0. + 1.
#             # Z_infeas_1[index_infeas_2] = np.nan
#             # Z_infeas_1[index_infeas_1] = np.nan
            
#             Z_infeas_2[index_feas_2] = np.nan
#             Z_infeas_3[(index_feas_1 | index_feas_2)] = np.nan
            
#             CS2 = plt.contourf(spec_1, spec_2, Z_infeas_1, zorder=0,
#                               levels=1, colors=[oversaccharine_shadecolor])
#             CS3 = plt.contourf(spec_1, spec_2, Z_infeas_2, zorder=0,
#                               levels=1, colors=[inhibited_shadecolor])
#             CS4 = plt.contourf(spec_1, spec_2, Z_infeas_3, zorder=0,
#                               levels=1, colors=[overlap_color])
            

                
                
#             # CS2_lines = plt.contour(CS2, zorder=1e9, linewidths=1000.,
#             #             levels=[0, 2], colors=[oversaccharine_shadecolor])
#             # CS3_lines = plt.contour(CS3, zorder=1e9, linewidths=1000.,
#             #             levels=[0, 2], colors=[inhibited_shadecolor])
#             # make_oversaccharine_region_infeasible()
#             # make_inhibited_region_infeasible()
            
#             CS1 = plt.contourf(spec_1, spec_2, MSP, zorder=1e6,
#                               levels=SA_price_range, colors=[marketrange_shadecolor])
#             CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
#                         levels=SA_price_range, colors=[linecolor_dark])
#             plt.clabel(CS1_lines, levels=SA_price_range, inline_spacing = 0., \
#                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
#             CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
#                         levels=SA_price_range, colors=[linecolor_dark])
#             plt.clabel(CS1_lines, levels=SA_price_range, inline_spacing = 0., 
#                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
#             CS1a_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
#                         levels=[2500, 3000, 3500], colors=[linecolor_dark])
#             plt.clabel(CS1a_lines, levels=[2500, 3000, 3500], inline_spacing = 0., \
#                        fmt=lambda x: format(x,'.0f'), inline=True,
#                        manual = [(0.97, 205), (0.85,154), (0.75,120)], fontsize=12)
#             CS1b_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
#                         levels=[4000, 4500], colors=[linecolor_dark])
#             plt.clabel(CS1b_lines, levels=[4000, 4500], inline_spacing = 0., \
#                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
#         elif j==1: # sugars plot only
#             CS2_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='solid', linewidths=0.9,
#             levels=[50, 100], colors=[linecolor_dark])
#             plt.clabel(CS2_lines, levels=[50,100], inline_spacing = 0., \
#                        manual = [(0.9, 135), (0.75,135)], 
#                        fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
#             CS2a_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
#             levels=[150.], colors=[linecolor_dark])
#             plt.clabel(CS2a_lines, levels=[150.], inline_spacing = 0., \
#                         manual = [(0.65,135)],
#                         fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
#         elif j==2: # inhibitors plot only
#             CS3_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
#                         levels=[1000.], colors=[linecolor_light])
#             plt.clabel(CS3_lines, levels=[1000.], inline_spacing = 0., \
#                         manual = [(0.45,110)],
#                         fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
#             CS3a_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='solid', linewidths=.7,
#                         levels=[600, 800], colors=[linecolor_dark])
#             plt.clabel(CS3a_lines, levels=[600, 800], inline_spacing = 0., \
#                         manual = [(0.725,110), (0.575,110)],
#                         fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
                
            
#             # infeas_1 = spec_1.copy()
#             # infeas_1[d1_Metric3<1000.] = np.nan
#             # infeas_1 = spec_1[np.where(d1_Metric3>1000.)[0:2]], 100*spec_2[np.where(d1_Metric3>1000.)[0:2]]
#             # infeas_2 = spec_1[np.where(d1_Metric2>150.)[0:2]], 100*spec_2[np.where(d1_Metric2>150.)[0:2]]
            
            
            
#             # ax.fill(infeas_1[0], infeas_1[1], [oversaccharine_shadecolor], zorder = 1e9)
#             # ax.fill(infeas_2[0], infeas_2[1], [inhibited_shadecolor], zorder = 1e9)
#             # for i in range(len(infeas_1[0])):
#             #     ax.scatter(infeas_1[0][i], infeas_1[1][i], color=[oversaccharine_shadecolor], plotnonfinite =True)
            
            
#         j+=1
        
# add_markers = False

# if add_markers:
    
#     axes_lab_spec_3 = axes[:, 0]
#     for i, ax in enumerate(axes_lab_spec_3):
#         plt.sca(ax)
#         # plt.clabel(CS, fmt=lambda x: format(x,'.0f'), inline=1, fontsize=12)
#         plot_scatter_points([lab_spec_1], [lab_spec_2], marker='^', s=80, color=markercolor,
#                             edgecolor=edgecolor)
    
#     axes_target_spec_3 = axes[:, 0]
#     for ax in axes_target_spec_3:
#         plt.sca(ax)
#         plot_scatter_points([target_spec_1], [target_spec_2], marker='s', s=80, color=markercolor,
#                             edgecolor=edgecolor)
    
# plt.show()

# # fig_to_save = file_to_save + '.png'

# # plt.savefig(fig_to_save, format = 'png', dpi=500)


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


def plot(data, titers, yields, productivities, 
         Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks):
    metric_bars = (MetricBar('MPSP\n', MPSP_units, CABBI_green_colormap(),
                             Metric_1_tickmarks,
                             1 + int((max(Metric_1_tickmarks) - min(Metric_1_tickmarks))/250)),
                   MetricBar(' \n', ' ',
                            CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
                             Metric_2_tickmarks, 
                             1 + int((max(Metric_2_tickmarks) - min(Metric_2_tickmarks))/6.25)),
                   MetricBar(' \n', ' ',
                             plt.cm.get_cmap('bone_r'), # plt.cm.get_cmap('bone_r'),
                             Metric_3_tickmarks,
                             1 + int((max(Metric_3_tickmarks) - min(Metric_3_tickmarks))/62.5)))
    
    return plot_contour_2d(titers, yields, productivities, data, 
                                xlabel, ylabel, xticks, yticks, metric_bars, 
                                Z_value_format=lambda Z: f"{Z:.1f} [{spec_3_units}]",)
                                # fillblack=False)






metric_bars = (MetricBar('MPSP', MPSP_units, CABBI_green_colormap(),
                         Metric_1_tickmarks, 800),
               MetricBar('Total sugars', VOC_units, CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
                         Metric_2_tickmarks, 80),
               MetricBar("Relative impact on MPSP\n[impact of titer : impact of yield]", FCI_units,
                         plt.cm.get_cmap('BrBG'), # plt.cm.get_cmap('bone_r'),
                         Metric_3_tickmarks, 100))
    
    
# if TAL_metrics[1] is get_TAL_sugars_conc:
# make_oversaccharine_region_infeasible()
# # if TAL_metrics[2] is get_TAL_inhibitors_conc:
# make_inhibited_region_infeasible()

fig, axes = plot(data_1_copy, spec_1, spec_2, spec_3, Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks)
# spec_2 = 100 * spec_2
index_feas_1 = d1_Metric2[:, :, 0]<150.
index_feas_2 = d1_Metric3[:, :, 0]<1000.
CS1, CS2, CS3, CS4 = 0, 0, 0, 0
for i, ax_col in enumerate(axes[:, :len(spec_3)].transpose()):
    MSP = data_1_copy[:, :, 0, i]
    sugars_orig = data_1[:, :, 1, i]
    inhibitors_orig = data_1[:, :, 2, i]
    MSP_orig = data_1[:, :, 0, i]
    j=0
    for ax in ax_col:
        
        plt.sca(ax)
        # ax.set_facecolor('black')
        # ax.patch.set_facecolor(CABBI_brown)
        ax.xaxis.set_minor_locator(AML())
        ax.yaxis.set_minor_locator(AML())
        ax.tick_params(which='both', top=True, bottom=True, left=True, right=True)
        ax.tick_params(which='minor', direction = 'inout', length = 4)
        ax.tick_params(which='major', length = 10)
        # ax.tick_params(which='minor', left=True, bottom = True, direction = 'inout')
        if j==0: # MPSP plot only
            Z_infeas_1 = MSP_orig * 0. + 1.
            Z_infeas_2 = MSP_orig * 0. + 1.
            Z_infeas_3 = MSP_orig * 0. + 1.
            # Z_infeas_1[index_infeas_2] = np.nan
            # Z_infeas_1[index_infeas_1] = np.nan
            
            Z_infeas_2[index_feas_2] = np.nan
            Z_infeas_3[(index_feas_1 | index_feas_2)] = np.nan
            
            # CS2 = plt.contourf(spec_1, spec_2, Z_infeas_1, zorder=0,
            #                   levels=1, colors=[oversaccharine_shadecolor])
            # CS3 = plt.contourf(spec_1, spec_2, Z_infeas_2, zorder=0,
            #                   levels=1, colors=[inhibited_shadecolor])
            # CS4 = plt.contourf(spec_1, spec_2, Z_infeas_3, zorder=0,
            #                   levels=1, colors=[overlap_color])
            

                
                
            # CS2_lines = plt.contour(CS2, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[oversaccharine_shadecolor])
            # CS3_lines = plt.contour(CS3, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[inhibited_shadecolor])
            # make_oversaccharine_region_infeasible()
            # make_inhibited_region_infeasible()
            
            CS1 = plt.contourf(spec_1, spec_2, MSP, zorder=1e6,
                              levels=SA_price_range, colors=[marketrange_shadecolor])
            CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=SA_price_range, colors=[linecolor_dark])
            plt.clabel(CS1_lines, levels=SA_price_range, inline_spacing = 0., \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
                       manual = [(0.925, 30), (0.925, 40)])
            # CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
            #             levels=SA_price_range, colors=[linecolor_dark])
            # plt.clabel(CS1_lines, levels=SA_price_range, inline_spacing = 0., 
            #            fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
            # CS1a_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
            #             levels=[2500, 3000, 3500], colors=[linecolor_dark],)
            # plt.clabel(CS1a_lines, levels=[2500, 3000, 3500], inline_spacing = 0., \
            #            fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
            #            manual = [(0.925, 90), (0.925, 80), (0.925, 70)])
            CS1b_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
                        levels=[3500, 4000, 4500, 5000], colors=[linecolor_dark])
            plt.clabel(CS1b_lines, levels=[3500, 4000, 4500, 5000], inline_spacing = 0., \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
                       manual = [(0.925, 80), (0.925, 70), (0.925, 60), (0.925, 50)])
        elif j==1: # sugars plot only
            CS2_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='solid', linewidths=0.9,
            levels=[25, 50], colors=[linecolor_dark])
            plt.clabel(CS2_lines, levels=[25, 50], inline_spacing = 0.,
                       fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
            CS2a_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
            levels=[150.], colors=[linecolor_light])
            plt.clabel(CS2a_lines, levels=[150.], inline_spacing = 0., \
                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
        elif j==2: # inhibitors plot only
            CS3_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=[1000.], colors=[linecolor_light])
            plt.clabel(CS3_lines, levels=[1000.], inline_spacing = 0., \
                        fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
            CS3a_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='solid', linewidths=.7,
                        levels=[250, 500, 750], colors=[linecolor_dark])
            plt.clabel(CS3a_lines, levels=[250, 500, 750], inline_spacing = 0.,
                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            # infeas_1 = spec_1.copy()
            # infeas_1[d1_Metric3<1000.] = np.nan
            # infeas_1 = spec_1[np.where(d1_Metric3>1000.)[0:2]], 100*spec_2[np.where(d1_Metric3>1000.)[0:2]]
            # infeas_2 = spec_1[np.where(d1_Metric2>150.)[0:2]], 100*spec_2[np.where(d1_Metric2>150.)[0:2]]
            
            
            
            # ax.fill(infeas_1[0], infeas_1[1], [oversaccharine_shadecolor], zorder = 1e9)
            # ax.fill(infeas_2[0], infeas_2[1], [inhibited_shadecolor], zorder = 1e9)
            # for i in range(len(infeas_1[0])):
            #     ax.scatter(infeas_1[0][i], infeas_1[1][i], color=[oversaccharine_shadecolor], plotnonfinite =True)
            
            
        j+=1
        
add_markers = False

if add_markers:
    
    axes_lab_spec_3 = axes[:, 0]
    for i, ax in enumerate(axes_lab_spec_3):
        plt.sca(ax)
        # plt.clabel(CS, fmt=lambda x: format(x,'.0f'), inline=1, fontsize=12)
        plot_scatter_points([lab_spec_1], [lab_spec_2], marker='^', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
    axes_target_spec_3 = axes[:, 0]
    for ax in axes_target_spec_3:
        plt.sca(ax)
        plot_scatter_points([target_spec_1], [target_spec_2], marker='s', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
plt.show()

# fig_to_save = file_to_save + '.png'

# plt.savefig(fig_to_save, format = 'png', dpi=500)

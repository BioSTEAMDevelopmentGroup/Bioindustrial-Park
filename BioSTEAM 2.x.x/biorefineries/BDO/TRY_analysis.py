# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: saran
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
from BDO.system import BDO_sys, BDO_tea, R302
from BDO.system import spec
from BDO.system import BDO as product
from matplotlib import pyplot as plt
from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from biosteam.plots import plot_contour_2d, MetricBar, plot_scatter_points#, CABBI_green_colormap
from math import floor, ceil
from datetime import datetime

ig = np.seterr(invalid='ignore')

# Colors
marketrange_shadecolor = (*colors.neutral.RGBn, 0.15)

oversaccharine_shadecolor_raw = colors.CABBI_teal_green.tint(40)
inhibited_shadecolor_raw = colors.CABBI_grey.shade(60)
oversaccharine_shadecolor = (*oversaccharine_shadecolor_raw.RGBn, 1)
inhibited_shadecolor = (*inhibited_shadecolor_raw.RGBn, 1)
# overlap_color = (*(colors.CABBI_teal_green.tint(20).RGBn + colors.CABBI_black.tint(20).RGBn)/2, 1)
overlap_color = (*(oversaccharine_shadecolor_raw.RGBn + inhibited_shadecolor_raw.RGBn)/2, 1)
linecolor_dark = (*colors.neutral_shade.RGBn, 0.85)
linecolor_light = (*colors.neutral_tint.RGBn, 0.85)
markercolor = (*colors.CABBI_orange.shade(5).RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)

def CABBI_green_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_yellow.RGBn,
                    colors.CABBI_green.RGBn,
                    colors.CABBI_teal_green.shade(62).RGBn)
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




million_dollar = r"\mathrm{MM\$}"
MPSP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"
productivity_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"
# VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
# FCI_units = f"${million_dollar}$"
VOC_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"


def tickmarks_from_data(data, accuracy=50, N_points=5):
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    dmin = floor(dmin/50) * 50
    dmax = ceil(dmax/50) * 50
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

target_yield = 90
target_titer = 180
target_productivity = 1.5
lab_yield = 80
lab_titer = 109.9
lab_productivity = 1


BDO_price_range = [2.3 * 907.185, 3.2 * 907.185] 
# temporary price range from https://www.alibaba.com/product-detail/Supply-2-3-butanediol-With-Good_62437411997.html?spm=a2700.7724857.normalList.2.a0b56ac0EDVZ7h&s=p&fullFirstScreen=true
# previous TEA for 2,3-BDO from monosaccharides: http://dx.doi.org/10.1016/j.biortech.2015.12.005

get_BDO_MPSP = lambda: BDO_tea.solve_price(product) * 907.185 # To USD / ton
get_BDO_VOC = lambda: BDO_tea.VOC / 1e6 # million USD / yr
get_BDO_FCI = lambda: BDO_tea.FCI / 1e6 # million USD

get_BDO_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol

get_BDO_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol

BDO_metrics = [get_BDO_MPSP, get_BDO_sugars_conc, get_BDO_inhibitors_conc]


# %% Generate TRY meshgrid

titers_1 = np.linspace(40, 200, 180)
yields_1 = np.linspace(0.45, 0.99, 180)
productivities_1 = np.array([1])
titers_1, yields_1 = np.meshgrid(titers_1, yields_1)

# titers_2 = np.linspace(70, 190, 5)
# yields_2 = np.linspace(0.70, 0.90, 5)
# productivities_2 = np.array([0.5, 1, 1.5])
# titers_2, yields_2 = np.meshgrid(titers_2, yields_2)

# utilizes_xylose = False


# %% Run TRY analysis 
data_1 = BDO_data = spec.evaluate_across_TRY(
        BDO_sys, titers_1, yields_1, BDO_metrics, productivities_1)

# utilizes_xylose = True
# data_2 = BDO_data = spec.evaluate_across_TRY(
#         BDO_sys, titers_2, yields_2, BDO_metrics, productivities_2)




# %% Save generated data


dateTimeObj = datetime.now()
file_to_save = 'BDO_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)
np.save(file_to_save, data_1)

# %% Load previously saved data
# file_to_load = file_to_save
file_to_load = 'BDO_TRY_2020.8.19-2.44'
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



# %% Plot contours
# data_2 = data_1

xlabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
ylabel = "Yield [%]"

xticks = [40, 80, 120, 160, 200]
yticks = [45, 54, 63, 72, 81, 90, 99]


def plot(data, titers, yields, productivities, 
         Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks):
    metric_bars = (MetricBar('MPSP', MPSP_units, CABBI_green_colormap(),
                             Metric_1_tickmarks, 200),
                   MetricBar('Total sugars', VOC_units, CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
                             Metric_2_tickmarks, 45),
                   MetricBar("Total inhibitors", FCI_units,
                             CABBI_grey_colormap(), # plt.cm.get_cmap('bone_r'),
                             Metric_3_tickmarks, 60))
    
    return plot_contour_2d(titers, 100.*yields, productivities, data, 
                                xlabel, ylabel, xticks, yticks, metric_bars, 
                                Z_value_format=lambda Z: f"{Z:.1f} [{productivity_units}]",
                                fillblack=False)




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

Metric_1_tickmarks = [1800, 2700, 3600, 4500, 5400]
Metric_2_tickmarks = [0, 75, 150, 225, 300]
Metric_3_tickmarks = [0, 300, 600, 900, 1200]
make_oversaccharine_region_infeasible()
make_inhibited_region_infeasible()

fig, axes = plot(data_1_copy, titers_1, yields_1, productivities_1, Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks)
percent_yields_1 = 100 * yields_1
index_feas_1 = d1_Metric2[:, :, 0]<150.
index_feas_2 = d1_Metric3[:, :, 0]<1000.
CS1, CS2, CS3, CS4 = 0, 0, 0, 0
for i, ax_col in enumerate(axes[:, :len(productivities_1)].transpose()):
    MSP = data_1_copy[:, :, 0, i]
    sugars_orig = data_1[:, :, 1, i]
    inhibitors_orig = data_1[:, :, 2, i]
    MSP_orig = data_1[:, :, 0, i]
    j=0
    for ax in ax_col:
        plt.sca(ax)
        # ax.patch.set_facecolor(CABBI_brown)
        
        if j==0: # MPSP plot only
            Z_infeas_1 = MSP_orig * 0. + 1.
            Z_infeas_2 = MSP_orig * 0. + 1.
            Z_infeas_3 = MSP_orig * 0. + 1.
            # Z_infeas_1[index_infeas_2] = np.nan
            # Z_infeas_1[index_infeas_1] = np.nan
            Z_infeas_2[index_feas_2] = np.nan
            Z_infeas_3[(index_feas_1 | index_feas_2)] = np.nan
            
            CS2 = plt.contourf(titers_1, percent_yields_1, Z_infeas_1, zorder=0,
                              levels=1, colors=[oversaccharine_shadecolor])
            CS3 = plt.contourf(titers_1, percent_yields_1, Z_infeas_2, zorder=0,
                              levels=1, colors=[inhibited_shadecolor])
            CS4 = plt.contourf(titers_1, percent_yields_1, Z_infeas_3, zorder=0,
                              levels=1, colors=[overlap_color])
            

                
                
            # CS2_lines = plt.contour(CS2, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[oversaccharine_shadecolor])
            # CS3_lines = plt.contour(CS3, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[inhibited_shadecolor])
            # make_oversaccharine_region_infeasible()
            # make_inhibited_region_infeasible()
            
            CS1 = plt.contourf(titers_1, percent_yields_1, MSP, zorder=1e6,
                              levels=BDO_price_range, colors=[marketrange_shadecolor])
            CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=BDO_price_range, colors=[linecolor_dark])
            plt.clabel(CS1_lines, levels=BDO_price_range, inline_spacing = 0.3, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=BDO_price_range, colors=[linecolor_dark])
            plt.clabel(CS1_lines, levels=BDO_price_range, inline_spacing = 0.3, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            
            CS1a_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.5,
                        levels=[1800, 2400, 3600, 4200], colors=[linecolor_dark])
            plt.clabel(CS1a_lines, levels=[1800, 2400, 3600, 4200], inline_spacing = 0.1, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
        elif j==1:
            CS2_lines = plt.contour(titers_1, percent_yields_1, sugars_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
            levels=[150.], colors=[linecolor_dark])
            plt.clabel(CS2_lines, levels=[150.], inline_spacing = 0.3, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            CS2a_lines = plt.contour(titers_1, percent_yields_1, sugars_orig, zorder=1e6, linestyles='solid', linewidths=.5,
            levels=[75.], colors=[linecolor_dark])
            plt.clabel(CS2a_lines, levels=[75.], inline_spacing = 0.1, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            
        elif j==2:
            CS3_lines = plt.contour(titers_1, percent_yields_1, inhibitors_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=[1000.], colors=[linecolor_light])
            plt.clabel(CS3_lines, levels=[1000.], inline_spacing = 0.3, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            CS3a_lines = plt.contour(titers_1, percent_yields_1, inhibitors_orig, zorder=1e6, linestyles='solid', linewidths=.5,
                        levels=[500.], colors=[linecolor_light])
            plt.clabel(CS3a_lines, levels=[500.], inline_spacing = 0.3, \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=10)
            # infeas_1 = titers_1.copy()
            # infeas_1[d1_Metric3<1000.] = np.nan
            # infeas_1 = titers_1[np.where(d1_Metric3>1000.)[0:2]], 100*yields_1[np.where(d1_Metric3>1000.)[0:2]]
            # infeas_2 = titers_1[np.where(d1_Metric2>150.)[0:2]], 100*yields_1[np.where(d1_Metric2>150.)[0:2]]
            
            
            
            # ax.fill(infeas_1[0], infeas_1[1], [oversaccharine_shadecolor], zorder = 1e9)
            # ax.fill(infeas_2[0], infeas_2[1], [inhibited_shadecolor], zorder = 1e9)
            # for i in range(len(infeas_1[0])):
            #     ax.scatter(infeas_1[0][i], infeas_1[1][i], color=[oversaccharine_shadecolor], plotnonfinite =True)
            
            
        j+=1
        
add_markers = True

if add_markers:
    
    axes_lab_productivity = axes[:, 0]
    for i, ax in enumerate(axes_lab_productivity):
        plt.sca(ax)
        # plt.clabel(CS, fmt=lambda x: format(x,'.0f'), inline=1, fontsize=10)
        plot_scatter_points([lab_titer], [lab_yield], marker='^', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
    axes_target_productivity = axes[:, 0]
    for ax in axes_target_productivity:
        plt.sca(ax)
        plot_scatter_points([target_titer], [target_yield], marker='s', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
plt.show()

fig_to_save = file_to_save + '.png'

plt.savefig(fig_to_save, format = 'png', dpi=500)



# fig, axes = plot(data_2, titers_2, yields_2, productivities_2, Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks)
# percent_yields_2 = 100 * yields_2
# for i, ax_col in enumerate(axes[:, :3].transpose()):
#     MSP = data_2[:, :, 0, i]
#     for ax in ax_col:
#         plt.sca(ax)
#         CS = plt.contourf(titers_2, percent_yields_2, MSP, zorder=1e6,
#                           levels=BDO_price_range, colors=[marketrange_shadecolor])
#         plt.contour(CS,zorder=1e6, linestyles='dashed', linewidths=1.,
#                     levels=BDO_price_range, colors=[linecolor_dark])

# axes_lab_productivity = axes[:, 1]
# for i, ax in enumerate(axes_lab_productivity):
#     plt.sca(ax)
#     # plt.clabel(CS, fmt=lambda x: format(x,'.0f'), inline=1, fontsize=10)
#     plot_scatter_points([lab_titer], [lab_yield], marker='p', s=125, color=markercolor,
#                         edgecolor=edgecolor)

# axes_target_productivity = axes[:, 1]
# for ax in axes_target_productivity:
#     plt.sca(ax)
#     plot_scatter_points([target_titer], [target_yield], marker='h', s=125, color=markercolor,
#                         edgecolor=edgecolor)
    


# plt.show()




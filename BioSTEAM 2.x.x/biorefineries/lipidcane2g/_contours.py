# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 23:44:10 2021

@author: yrc2
"""
from biorefineries import lipidcane2g as lc
import biosteam as bst
import numpy as np
import pandas as pd
from biosteam.utils import colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from biosteam.utils import colors
from biosteam.plots import (
    plot_contour_2d,
    MetricBar,
    plot_scatter_points,
    plot_contour_single_metric,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks
)
from math import floor, ceil
from biosteam import plots
from biosteam.utils import CABBI_colors
from thermosteam.units_of_measure import format_units
from biosteam.plots.utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
from biosteam import Metric
from warnings import filterwarnings
import os

__all__ = ('plot_extraction_efficiency_and_lipid_content_contours',
           'plot_ethanol_and_biodiesel_price_contours')

filterwarnings('ignore', category=bst.utils.DesignWarning)
    
shadecolor = (*colors.neutral.RGBn, 0.20)
linecolor = (*colors.neutral_shade.RGBn, 0.85)
markercolor = (*colors.orange_tint.RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)
    
CABBI_colors = (colors.CABBI_yellow.tint(75).RGBn, 
                colors.CABBI_yellow.RGBn,
                colors.CABBI_green.RGBn,
                colors.CABBI_teal_green.shade(60).RGBn)

CABBI_colors_x = (colors.CABBI_blue_light.tint(90).RGBn,
                  colors.CABBI_blue_light.tint(40).RGBn, 
                  colors.CABBI_blue_light.RGBn, 
#                   colors.CABBI_teal.RGBn,
#                   colors.CABBI_teal_green.tint(10).RGBn,
                  colors.CABBI_teal_green.tint(40).shade(15).RGBn,
                  colors.CABBI_teal_green.shade(45).RGBn)

diverging_colormaps = [
    plt.cm.get_cmap('RdYlGn')
]

colormaps = [
    LinearSegmentedColormap.from_list('CABBI', CABBI_colors, 25),
    LinearSegmentedColormap.from_list('CABBI', CABBI_colors_x, 25),
    plt.cm.get_cmap('inferno_r'),
    plt.cm.get_cmap('copper_r'),
    plt.cm.get_cmap('bone_r'),
] * 2

def plot_ethanol_and_biodiesel_price_contours(N=30, benefit=False, cache={}):
    ethanol_price = np.linspace(1., 3., N)
    biodiesel_price = np.linspace(2, 6, N)
    lipid_content = [5, 10, 15]
    N_rows = len(lipid_content)
    configuration = ['L1', 'L1*', 'L2', 'L2*']
    N_cols = len(configuration)
    if (N, benefit) in cache:
        Z = cache[N, benefit]
    else:
        Z = np.zeros([N, N, N_rows, N_cols])
        for i in range(N_rows):
            for j in range(N_cols):
                lc.load(configuration[j])
                lc.set_feedstock_lipid_content.setter(lipid_content[i])
                lc.lipidcane_sys.simulate()
                X, Y = np.meshgrid(ethanol_price, biodiesel_price)
                if benefit:
                    Z[:, :, i, j] = lc.evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices(X, Y)
                else:
                    Z[:, :, i, j] = lc.evaluate_MFPP_across_ethanol_and_biodiesel_prices(X, Y)
    xlabel = f"Ethanol price [{format_units('$/gal')}]"
    ylabels = [f"Biodiesel price [{format_units('$/gal')}]"] * 4
    xticks = [1., 1.5, 2.0, 2.5, 3.0]
    yticks = [2, 3, 4, 5, 6]
    marks = tickmarks(Z, 5, 5, center=0.) if benefit else tickmarks(Z, 5, 5)
    colormap = (diverging_colormaps if benefit else colormaps)[0]
    metric_bar = MetricBar('MFPP', format_units('$/ton'), colormap, marks, 12,
                           center=0.)
    if benefit:
        baseline = ['S1', 'S1*', 'S2', 'S2*']
        titles = [f"{lc.format_name(i)} - {lc.format_name(j)}"
                  for i, j in  zip(configuration, baseline)]
    else:
        titles = [lc.format_name(i) for i in configuration]
    fig, axes, CSs, CB = plot_contour_single_metric(
        X, Y, Z, xlabel, ylabels, xticks, yticks, metric_bar, 
        fillblack=False, styleaxiskw=dict(xtick0=False), label=True,
        titles=titles,
    )
    for i in range(N_rows):
        for j in range(N_cols):
            ax = axes[i, j]
            CS = CSs[i, j]
            plt.sca(ax)
            data = Z[:, :, i, j]
            lb = data.min()
            ub = data.max()
            levels = [i for i in CS.levels if lb <= i <= ub]
            CS = plt.contour(X, Y, data=data, zorder=1e16, linestyles='dashed', linewidths=1.,
                             levels=levels, colors=[linecolor])
            ax.clabel(CS, levels=CS.levels, inline=True, fmt=lambda x: f'{round(x):,}',
                      fontsize=10, colors=[linecolor], zorder=1e16)

    plt.show()
    
def plot_extraction_efficiency_and_lipid_content_contours(load=False, metric_index=0):
    # Generate contour data
    x = np.linspace(0.4, 1., 5)
    y = np.linspace(0.05, 0.15, 5)
    X, Y = np.meshgrid(x, y)
    # dollar_per_mt = format_units(r'\$/MT')
    metric = bst.metric
    kg_per_ton = 907.18474
    
    # NG = None
    # @metric(units=format_units(r'$/ton'))
    # def MFPP():
    #     return kg_per_ton * lc.lipidcane_tea.solve_price(lc.lipidcane)

    # @metric(units=format_units(r'10^6*$'))
    # def TCI():
    #     return lc.lipidcane_tea.TCI / 1e6 # 10^6*$
    folder = os.path.dirname(__file__)
    file = 'lipid_extraction_analysis.npy'
    file = os.path.join(folder, file)
    configurations = [1, 2]
    agile = [False, True]
    if load:
        data = np.load(file)
    else:
        data = lc.evaluate_configurations_across_extraction_efficiency_and_lipid_content(
            X, Y, 0.70, agile, configurations, 
        )
    np.save(file, data)
    data = data[:, :, :, :, metric_index]
    
    # Plot contours
    xlabel = 'Lipid extraction [%]'
    ylabel = "Lipid content [dry wt. %]"
    ylabels = [f'Lipid-cane only\n{ylabel}',
               f'Lipid-cane & lipid-sorghum\n{ylabel}']
    xticks = [40, 60, 80, 100]
    yticks = [5, 7.5, 10, 12.5, 15]
    metric = lc.all_metric_mockups[metric_index]
    metric_bar = MetricBar(metric.name, format_units(metric.units), colormaps[metric_index], tickmarks(data, 5, 5), 18)
    fig, axes, CSs, CB = plot_contour_single_metric(
        100.*X, 100.*Y, data, xlabel, ylabels, xticks, yticks, metric_bar, 
        fillblack=False, styleaxiskw=dict(xtick0=False), label=True,
        titles=['Configuration I', 'Configuration II'],
    )
    M = len(configurations)
    N = len(agile)
    for i in range(M):
        for j in range(N):
            ax = axes[i, j]
            CS = CSs[i, j]
            plt.sca(ax)
            metric_data = data[:, :, i, j]
            lb = metric_data.min()
            ub = metric_data.max()
            levels = [i for i in CS.levels if lb <= i <= ub]
            CS = plt.contour(100.*X, 100.*Y, data=metric_data, zorder=1e16, linestyles='dashed', linewidths=1.,
                             levels=levels, colors=[linecolor])
            ax.clabel(CS, levels=CS.levels, inline=True, fmt=lambda x: f'{round(x):,}',
                      fontsize=10, colors=[linecolor], zorder=1e16)
            if j == 0:
                lb = 47.5
                ub = 52.5
            else:
                lb = 75
                ub = 80
            baseline = (lb + ub) / 2.
            plt.fill_between([lb, ub], [2], [20], 
                             color=shadecolor,
                             linewidth=1)
            plot_vertical_line(lb, ls='-.',
                               color=linecolor,
                               linewidth=1.0)
            plot_vertical_line(ub, ls='-.',
                               color=linecolor,
                               linewidth=1.0)
            plot_scatter_points([baseline], [10], marker='*', s=125, color=markercolor,
                                edgecolor=edgecolor)

    plt.show()
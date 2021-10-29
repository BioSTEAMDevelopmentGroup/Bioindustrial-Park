# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 03:21:03 2021

@author: yrc2
"""
from biorefineries import actag
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
    plot_contour_1d,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks,
)
from math import floor, ceil
from biosteam import plots
from biosteam.utils import CABBI_colors
from thermosteam.units_of_measure import format_units
from biosteam.plots.utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
from biosteam import Metric
from warnings import filterwarnings
import os
from thermosteam.units_of_measure import format_units

__all__ = ('plot_yield_titer_selectivity_productivity_contours',
           'plot_purity_across_selectivity',
           'plot_metrics_across_titer',
           'plot_tea_across_titer',
           'plot_MPSP_across_titer_and_yield')

baselinecolor = (*colors.red.RGBn, 1)
shadecolor = (*colors.neutral.RGBn, 0.20)
linecolor = (*colors.neutral_shade.RGBn, 0.85)
markercolor = (*colors.orange_tint.RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)
    
CABBI_colors = (colors.CABBI_yellow.tint(75).RGBn, 
                colors.CABBI_yellow.RGBn,
                colors.CABBI_green.RGBn,
                colors.CABBI_teal_green.shade(60).RGBn,
                colors.CABBI_teal_green.shade(90).RGBn)

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

### Defaults
rho = 0.900 # kg / L
f = 1 / rho * 907.1847 # L / ton
lubricating_oil_market_price = (0.92 * f, 1.35 * f) # USD / L to USD / ton

def plot_yield_titer_selectivity_productivity_contours(
        configuration=1, load=True, price_ranges=[lubricating_oil_market_price],
        metric_index=0,
    ):
    baseline = actag.baseline[configuration]
    target = actag.target[configuration]
    baseline = [[baseline['Yield']], [baseline['Titer']]]
    target  = [[target['Yield']], [target['Titer']]]
    
    # Generate contour data
    X, Y, z, w, data = actag.fermentation_data(configuration, load)
    
    data = data[:, :, :, :, metric_index]
    
    # Plot contours
    xlabel = "Yield\n[% theoretical]" 
    ylabel = f"Titer\n[{format_units('g/L')}]"
    xticks = [30, 45, 60, 75, 90]
    yticks = [5, 25, 50, 75, 100]
    if metric_index == 0:
        metric_bar = MetricBar(
            'MSP', format_units('$/ton'), 
            colormaps[metric_index],
            tickmarks(data, 5, 5), 15
        )
    elif metric_index == 1:
        metric_bar = MetricBar(
            'TCI', format_units('10^6*$'), 
            colormaps[metric_index], 
            tickmarks(data, 5, 5), 15
        )
    elif metric_index == 2:
        metric_bar = MetricBar(
            'Heating', format_units('GJ/hr'), 
            colormaps[metric_index], 
            tickmarks(data, 5, 5), 15
        )
    elif metric_index == 3:
        metric_bar = MetricBar(
            'Cooling', format_units('GJ/hr'), 
            colormaps[metric_index], 
            tickmarks(data, 5, 5), 15
        )
    X = X[:, :, 0]
    Y = Y[:, :, 0]
    fig, axes, CSs, CB = plot_contour_single_metric(
        X, Y, data, xlabel, ylabel, xticks, yticks, metric_bar, 
        fillblack=False, styleaxiskw=dict(xtick0=False), label=False,
    )
    
    *_, nrows, ncols = data.shape
    # colors = [linecolor]
    # hatches = ['//', r'\\']
    if price_ranges:
        for i, price_range in enumerate(price_ranges):
            for row in range(nrows):
                for col in range(ncols):
                    ax = axes[row, col]
                    plt.sca(ax)
                    metric_data = data[:, :, row, col]
                    # csf = ax.contourf(X, Y, metric_data, hatches=hatches[i],
                    #                  levels=np.linspace(*price_range, 5), colors='none')
                    # # For each level, set the color of its hatch 
                    # for collection in csf.collections: collection.set_edgecolor(colors[i])
                    # # Doing this also colors in the box around each level
                    # # We can remove the colored line around the levels by setting the linewidth to 0
                    # for collection in csf.collections: collection.set_linewidth(0.)
                    cs = plt.contour(X, Y, metric_data, zorder=1e6, linestyles='dashed', linewidths=1.,
                                     levels=price_range, colors=[linecolor])
                    # clabels = ax.clabel(cs, levels=cs.levels, inline=True, fmt=metric_bar.fmt,
                    #                     fontsize=12, colors=['k'], zorder=1e16)
                    # for i in clabels: i.set_rotation(0)
    plt.sca(axes[0, 0])
    plt.scatter(*baseline, color=baselinecolor, marker='o', s=75, edgecolor='black', zorder=1e6)
    plt.sca(axes[1, 1])
    plt.scatter(*target, color=baselinecolor, marker='*', s=125, edgecolor='black', zorder=1e6)
    plt.show()
    # breakpoint()

def plot_purity_across_selectivity(configuration=1):
    actag.load(configuration)
    x = np.linspace(50, 90, 20)
    def simulate_purity(selectivity):
        actag.set_selectivity(selectivity)
        actag.sys.simulate()
        return actag.acTAG.imass['AcetylDiOlein'] / actag.acTAG.F_mass
    y = [simulate_purity(i) for i in x]
    plt.plot(x, y)
    
def plot_metrics_across_titer(configuration=1):
    actag.load(configuration)
    x = np.linspace(5, 80, 20)
    def simulate_titer(titer):
        actag.fermentation.titer = titer
        actag.sys.simulate()
        return [i() for i in actag.model.metrics]
    data = np.array([simulate_titer(i) for i in x])
    M, N = data.shape
    fig, axes = plt.subplots(N)
    for i in range(N):
        ax = axes[i]
        y = data[:, i]
        plt.sca(ax)
        plt.plot(x, y)
    plt.show()
    
def plot_tea_across_titer(configuration=1):
    actag.load(configuration)
    x = np.linspace(5, 80, 30)
    actag.fermentation.product_yield = 0.70
    def simulate_titer(titer):
        actag.fermentation.titer = titer
        actag.sys.simulate()
        return [actag.MPSP(), actag.acTAG.F_mass / 1e3, actag.tea.TCI / 1e6, actag.tea.VOC / 1e6]
    data = np.array([simulate_titer(i) for i in x])
    M, N = data.shape
    fig, axes = plt.subplots(N)
    for i in range(N):
        ax = axes[i]
        y = data[:, i]
        plt.sca(ax)
        plt.plot(x, y)
    plt.show()
    
def plot_MPSP_across_titer_and_yield(configuration=1):
    actag.load(configuration)
    x = np.linspace(5, 80, 20)
    def simulate_titer(titer):
        actag.fermentation.titer = titer
        values = []
        for i in [0.50, 0.70, 0.90]:
            actag.fermentation.product_yield = i
            actag.sys.simulate()
            values.append(actag.MPSP())
        return values
    data = np.array([simulate_titer(i) for i in x])
    M, N = data.shape
    fig, axes = plt.subplots(N)
    for i in range(N):
        ax = axes[i]
        y = data[:, i]
        plt.sca(ax)
        plt.plot(x, y)
    plt.show()
    
def plot_recovery_across_selectivity(configuration=1):
    actag.load(configuration)
    x = np.linspace(50, 90, 20)
    def simulate_recovery(selectivity):
        actag.set_selectivity(selectivity)
        actag.sys.simulate()
        recovered = actag.acTAG.imass['AcetylDiOlein']
        total = actag.fermentation.outs[1].imass['AcetylDiOlein']
        return recovered / total
    y = [simulate_recovery(i) for i in x]
    plt.plot(x, y)
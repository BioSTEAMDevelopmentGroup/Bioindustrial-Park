# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 23:44:10 2021

@author: yrc2
"""
from biorefineries import oilcane as oc
import biosteam as bst
import numpy as np
from biosteam.utils import colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from biosteam.plots import (
    plot_contour_2d,
    MetricBar,
    plot_contour_single_metric,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks,
    plot_scatter_points,
)
from thermosteam.units_of_measure import format_units
from thermosteam.utils import set_figure_size, set_font
from biorefineries.oilcane._load_data import images_folder
from warnings import filterwarnings
import os
from ._distributions import (
    biodiesel_prices_quarter,
    ethanol_prices,
)

__all__ = (
    'plot_relative_sorghum_oil_content_and_cane_oil_content_contours_manuscript',
    'plot_extraction_efficiency_and_oil_content_contours_manuscript',
    'plot_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_enhanced_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_benefit_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_extraction_efficiency_and_oil_content_contours',
    'plot_relative_sorghum_oil_content_and_cane_oil_content_contours',
    'plot_ethanol_and_biodiesel_price_contours'
)

filterwarnings('ignore', category=bst.utils.DesignWarning)
    
shadecolor = (*colors.neutral.RGBn, 0.20)
linecolor = (*colors.neutral_shade.RGBn, 0.85)
targetcolor = (*colors.red_tint.RGBn, 1)
startcolor = (*colors.purple.tint(10).RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)
    
CABBI_colors = (colors.CABBI_yellow.tint(75).RGBn, 
                colors.CABBI_yellow.tint(30).RGBn,
                colors.CABBI_yellow.RGBn,
                colors.CABBI_green.tint(5).RGBn,
                colors.CABBI_teal_green.shade(40).RGBn,
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

light_letter_color = colors.neutral.tint(98).RGBn
letter_color = colors.neutral.tint(80).RGBn
dark_letter_color = colors.neutral.shade(80).RGBn

# %% Plot functions for publication

def _add_letter_labels(axes, xpos, ypos, colors):
    M, N = shape = colors.shape
    if shape == (2, 2):
        letters=np.array([['A', 'C'], ['B', 'D']])
    elif shape == (3, 2):
        letters=np.array([['A', 'D'], ['B', 'E'], ['C', 'F']])
    for i in range(M):
        for j in range(N):
            ax = axes[i, j]
            letter = letters[i, j]
            xlb, xub = ax.get_xlim()
            ylb, yub = ax.get_ylim()
            if hasattr(ax, '_cached_ytwin'):                
                ax = ax._cached_ytwin
            ax.text((xlb + xub) * xpos, (yub + ylb) * ypos, letter, color=colors[i, j],
                     horizontalalignment='center',verticalalignment='center',
                     fontsize=12, fontweight='bold', zorder=1e17)

def plot_extraction_efficiency_and_oil_content_contours_manuscript():
    set_font(size=8)
    set_figure_size()
    fig, axes = plot_extraction_efficiency_and_oil_content_contours(
        load=True,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[dark_letter_color, dark_letter_color],
                 [dark_letter_color, dark_letter_color]]
    _add_letter_labels(axes, 0.68, 0.70, colors)
    plt.subplots_adjust(right=0.92, wspace=0.1, top=0.9, bottom=0.10)
    file = os.path.join(images_folder, 'extraction_efficiency_and_oil_content_contours.svg')
    plt.savefig(file, transparent=True)

def plot_relative_sorghum_oil_content_and_cane_oil_content_contours_manuscript():
    set_font(size=8)
    set_figure_size()
    fig, axes = plot_relative_sorghum_oil_content_and_cane_oil_content_contours(
        load=True,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [dark_letter_color, light_letter_color]]
    _add_letter_labels(axes, 0.82, 0.70, colors)
    plt.subplots_adjust(right=0.92, wspace=0.1, top=0.9, bottom=0.12)
    file = os.path.join(images_folder, 'relative_sorghum_oil_content_and_cane_oil_content_contours.svg')
    plt.savefig(file, transparent=True)

def plot_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        titles=['CONVENTIONAL', 'CELLULOSIC'],
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[letter_color, letter_color],
                 [letter_color, letter_color],
                 [letter_color, letter_color]]
    _add_letter_labels(axes, 0.7, 0.65, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    file = os.path.join(images_folder, 'ethanol_and_biodiesel_price_contours.svg')
    plt.savefig(file, transparent=True)

def plot_enhanced_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        enhanced_cellulosic_performance=True,
        titles=['CONVENTIONAL', 'CELLULOSIC'],
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[letter_color, letter_color],
                 [letter_color, letter_color],
                 [letter_color, letter_color]]
    _add_letter_labels(axes, 0.7, 0.65, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    file = os.path.join(images_folder, 'enhanced_ethanol_and_biodiesel_price_contours.svg')
    plt.savefig(file, transparent=True)
    
def plot_benefit_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        benefit=True,
        titles=['CONVENTIONAL', 'CELLULOSIC'],
        dist=True,
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[dark_letter_color, dark_letter_color],
                 [dark_letter_color, dark_letter_color],
                 [dark_letter_color, dark_letter_color]]
    _add_letter_labels(axes, 0.7, 0.65, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    file = os.path.join(images_folder, 'benefit_ethanol_and_biodiesel_price_contours.svg')
    plt.savefig(file, transparent=True)

# %%

def plot_ethanol_and_biodiesel_price_contours(N=30, benefit=False, cache={}, 
                                              enhanced_cellulosic_performance=False,
                                              titles=None, load=True, save=True,
                                              dist=False):
    ethanol_price = np.linspace(1., 3., N)
    biodiesel_price = np.linspace(2, 6, N)
    oil_content = [5, 10, 15]
    N_rows = len(oil_content)
    configuration = ['O1', 'O2']
    N_cols = len(configuration)
    folder = os.path.dirname(__file__)
    file = 'price_contours'
    if benefit:
        file += "_benefit"
    if enhanced_cellulosic_performance:
        file += '_enhanced'
    file = os.path.join(folder, file + '.npy')
    X, Y = np.meshgrid(ethanol_price, biodiesel_price)
    if load:
        try: Z = np.load(file)
        except: load = True
        else: load = False
    if load:
        if (N, benefit, enhanced_cellulosic_performance) in cache:
            Z = cache[N, benefit, enhanced_cellulosic_performance]
        else:
            Z = np.zeros([N, N, N_rows, N_cols])
            for i in range(N_rows):
                for j in range(N_cols):
                    oc.load(configuration[j], enhanced_cellulosic_performance=enhanced_cellulosic_performance)
                    oc.set_cane_oil_content(oil_content[i])
                    oc.set_relative_sorghum_oil_content(0)
                    oc.sys.simulate()
                    if benefit:
                        Z[:, :, i, j] = oc.evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices(X, Y)
                    else:
                        Z[:, :, i, j] = oc.evaluate_MFPP_across_ethanol_and_biodiesel_prices(X, Y)
        if save:
            np.save(file, Z)
    xlabel = f"Ethanol price\n[{format_units('$/gal')}]"
    ylabels = [f"Biodiesel price\n[{format_units('$/gal')}]"] * 4
    xticks = [1., 1.5, 2.0, 2.5, 3.0]
    yticks = [2, 3, 4, 5, 6]
    marks = tickmarks(Z, 5, 5, center=0.) if benefit else tickmarks(Z, 5, 5)
    colormap = (diverging_colormaps if benefit else colormaps)[0]
    metric_bar = MetricBar('MFPP', format_units('$/ton'), colormap, marks, 15,
                           forced_size=0.3)
    if benefit:
        baseline = ['S1', 'S2']
        if titles is None:
            titles = [f"{oc.format_name(i)} - {oc.format_name(j)}"
                      for i, j in  zip(configuration, baseline)]
    elif titles is None:
        titles = [oc.format_name(i) for i in configuration]
    fig, axes, CSs, CB = plot_contour_single_metric(
        X, Y, Z, xlabel, ylabels, xticks, yticks, metric_bar, 
        styleaxiskw=dict(xtick0=False), label=True,
        titles=titles,
    )
    # for i in range(N_rows):
    #     for j in range(N_cols):
    #         ax = axes[i, j]
    #         CS = CSs[i, j]
    #         plt.sca(ax)
    #         data = Z[:, :, i, j]
    #         lb = data.min()
    #         ub = data.max()
    #         levels = [i for i in CS.levels if lb <= i <= ub]
    #         CS = plt.contour(X, Y, data=data, zorder=1, linestyles='dashed', linewidths=1.,
    #                          levels=levels, colors=[linecolor])
    #         ax.clabel(CS, levels=CS.levels, inline=True, fmt=lambda x: f'{round(x):,}',
    #                   colors=[linecolor], zorder=1)
    for ax in axes.flatten():
        try: fig.sca(ax)
        except: continue
        plot_scatter_points(ethanol_prices, biodiesel_prices_quarter, 
                            marker='o', s=2, color=dark_letter_color,
                            edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
    
def relative_sorghum_oil_content_and_cane_oil_content_data(load, relative):
    # Generate contour data
    y = np.linspace(0.05, 0.15, 20)
    x = np.linspace(-0.03, 0., 20) if relative else np.linspace(0.02, 0.15, 20)
    X, Y = np.meshgrid(x, y)
    folder = os.path.dirname(__file__)
    file = 'oil_content_analysis.npy'
    if relative: file = 'relative_' + file
    file = os.path.join(folder, file)
    configurations = [1, 2]
    if load:
        data = np.load(file)
    else:
        data = oc.evaluate_configurations_across_sorghum_and_cane_oil_content(
            X, Y, configurations, relative,
        )
    np.save(file, data)
    return X, Y, data
    
def plot_relative_sorghum_oil_content_and_cane_oil_content_contours(
        load=False, configuration_index=..., relative=False
    ):
    # Generate contour data
    X, Y, data = relative_sorghum_oil_content_and_cane_oil_content_data(load, relative)
    
    data = data[:, :, configuration_index, [0, 6]]
    
    # Plot contours
    xlabel = "Sorghum oil content\n[dry wt. %]" 
    if relative: xlabel = ('relative ' + xlabel).capitalize()
    ylabel = 'Cane oil content\n[dry wt. %]'
    yticks = [5, 7.5, 10, 12.5, 15]
    xticks = [-3, -2, -1, 0] if relative else [2, 5, 7.5, 10, 12.5, 15]
    MFPP = oc.all_metric_mockups[0]
    TCI = oc.all_metric_mockups[6]
    if configuration_index == 0:
        Z = np.array(["AGILE-CONVENTIONAL"])
        data = data[:, :, :, np.newaxis]
    elif configuration_index == 1:
        Z = np.array(["AGILE-CELLULOSIC"])
        data = data[:, :, :, np.newaxis]
    elif configuration_index == ...:
        Z = np.array(["AGILE-CONVENTIONAL", "AGILE-CELLULOSIC"])
        data = np.swapaxes(data, 2, 3)
    else:
        raise ValueError('configuration index must be either 0 or 1')
    metric_bars = [
        MetricBar(MFPP.name, format_units(MFPP.units), colormaps[0], tickmarks(data[:, :, 0], 5, 5), 15, 1),
        MetricBar(TCI.name, format_units(TCI.units), colormaps[1], tickmarks(data[:, :, 1], 5, 5), 29),
    ]
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X, 100.*Y, Z, data, xlabel, ylabel, xticks, yticks, metric_bars, 
        styleaxiskw=dict(xtick0=True), label=True,
    )
    return fig, axes
    
def plot_extraction_efficiency_and_oil_content_contours(
        load=False, metric_index=0, N_decimals=0,
    ):
    # Generate contour data
    x = np.linspace(0.4, 1., 8)
    y = np.linspace(0.05, 0.15, 8)
    X, Y = np.meshgrid(x, y)
    metric = bst.metric
    folder = os.path.dirname(__file__)
    file = 'oil_extraction_analysis.npy'
    file = os.path.join(folder, file)
    configurations = [1, 2]
    agile = [False, True]
    if load:
        data = np.load(file)
    else:
        data = oc.evaluate_configurations_across_extraction_efficiency_and_oil_content(
            X, Y, 0.70, agile, configurations, 
        )
    np.save(file, data)
    data = data[:, :, :, :, metric_index]
    
    # Plot contours
    xlabel = 'Oil extraction [%]'
    ylabel = "Oil content [dry wt. %]"
    ylabels = [f'Oilcane only\n{ylabel}',
               f'Oilcane & oilsorghum\n{ylabel}']
    xticks = [40, 60, 80, 100]
    yticks = [5, 7.5, 10, 12.5, 15]
    metric = oc.all_metric_mockups[metric_index]
    units = metric.units if metric.units == '%' else format_units(metric.units)
    metric_bar = MetricBar(
        metric.name, units, colormaps[metric_index], 
        tickmarks(data, 5, 5), 18, N_decimals=N_decimals,
        forced_size=0.3
    )
    fig, axes, CSs, CB = plot_contour_single_metric(
        100.*X, 100.*Y, data, xlabel, ylabels, xticks, yticks, metric_bar, 
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
        titles=['CONVENTIONAL', 'CELLULOSIC'],
    )
    M = len(configurations)
    N = len(agile)
    for i in range(N):
        for j in range(M):
            ax = axes[i, j]
            CS = CSs[i, j]
            plt.sca(ax)
            metric_data = data[:, :, i, j]
            lb = metric_data.min()
            ub = metric_data.max()
            levels = [i for i in CS.levels if lb <= i <= ub]
            CS = plt.contour(100.*X, 100.*Y, data=metric_data, zorder=1, linestyles='dashed', linewidths=1.,
                             levels=levels, colors=[linecolor])
            # ax.clabel(CS, levels=CS.levels, inline=True, fmt=lambda x: f'{round(x):,}',
            #           colors=[linecolor], zorder=1e16)
            if j == 0:
                lb = 50.0
                ub = 70.0
            else:
                lb = 70.0
                ub = 90.0
            plt.fill_between([lb, ub], [2], [20], 
                             color=shadecolor,
                             linewidth=1)
            plot_vertical_line(lb, ls='-.',
                               color=linecolor,
                               linewidth=1.0)
            plot_vertical_line(ub, ls='-.',
                               color=linecolor,
                               linewidth=1.0)
            if hasattr(ax, '_cached_ytwin'):                
                plt.sca(ax._cached_ytwin)
            plot_scatter_points([lb], [5], marker='o', s=50, color=startcolor,
                                edgecolor=edgecolor, clip_on=False, zorder=3)
            plot_scatter_points([ub], [15], marker='*', s=100, color=targetcolor,
                                edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
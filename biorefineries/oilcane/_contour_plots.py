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
from scipy.ndimage.filters import gaussian_filter
import os
from ._distributions import (
    biodiesel_prices,
    ethanol_no_RIN_prices,
    advanced_ethanol_prices,
    cellulosic_ethanol_prices,
)

__all__ = (
    'plot_relative_sorghum_oil_content_and_cane_oil_content_contours_manuscript',
    'plot_recovery_and_oil_content_contours_manuscript',
    'plot_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_enhanced_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_benefit_ethanol_and_biodiesel_price_contours_manuscript',
    'plot_recovery_and_oil_content_contours',
    'plot_relative_sorghum_oil_content_and_cane_oil_content_contours',
    'plot_ethanol_and_biodiesel_price_contours',
    'plot_recovery_and_oil_content_contours_biodiesel_only',
)

filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    
shadecolor = (*colors.neutral.RGBn, 0.20)
linecolor = (*colors.neutral_shade.RGBn, 0.85)
targetcolor = (*colors.red_tint.RGBn, 1)
startcolor = (*colors.red_tint.RGBn, 1)
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
    plt.cm.get_cmap('viridis'),
    plt.cm.get_cmap('copper_r'),
    # LinearSegmentedColormap.from_list('CABBI', CABBI_colors, 25),
    # LinearSegmentedColormap.from_list('CABBI', CABBI_colors_x, 25),
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

def plot_recovery_and_oil_content_contours_manuscript(load=True, fs=8, smooth=1):
    set_font(size=fs)
    set_figure_size()
    fig, axes = plot_recovery_and_oil_content_contours(
        load=load, 
        smooth=smooth,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 1 - 0.68, 0.7, colors)
    plt.subplots_adjust(right=0.92, wspace=0.1 * (fs/8) ** 2, top=0.9, bottom=0.10)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'recovery_and_oil_content_contours.{i}')
        plt.savefig(file, transparent=True)
        
def plot_recovery_and_oil_content_contours_biodiesel_only(load=True, fs=8, metric_indices=None):
    set_font(size=fs)
    set_figure_size()
    if metric_indices is None: metric_indices = (0, 2, 6, 10)
    for cmap, i in zip(colormaps, metric_indices):
        fig, axes = plot_recovery_and_oil_content_contours(
            load=load, configurations=np.array([[7, 8], [5, 6]]),
            N_points=20, yticks=[0, 2.5, 5, 7.5, 10, 12.5, 15],
            titles=['Batch', 'Fed-Batch'],
            metric_index=i, cmap=cmap,
        )
        load = True
        colors = np.zeros([2, 2], object)
        colors[:] = [[light_letter_color, light_letter_color],
                     [light_letter_color, light_letter_color]]
        _add_letter_labels(axes, 1 - 0.68, 0.85, colors)
        plt.subplots_adjust(right=0.92, wspace=0.1 * (fs/8) ** 2, top=0.9, bottom=0.10)
        for j in ('svg', 'png'):
            file = os.path.join(images_folder, f'recovery_and_oil_content_contours_biodiesel_only_{i}.{j}')
            plt.savefig(file, transparent=True)

def plot_relative_sorghum_oil_content_and_cane_oil_content_contours_manuscript(load=True, fs=8, smooth=0.9):
    set_font(size=fs)
    set_figure_size()
    fig, axes = plot_relative_sorghum_oil_content_and_cane_oil_content_contours(
        load=load, smooth=smooth,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 1 - 0.82, 0.70, colors)
    plt.subplots_adjust(right=0.92, wspace=0.1, top=0.9, bottom=0.12)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'relative_sorghum_oil_content_and_cane_oil_content_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        titles=['Direct Cogeneration', 'Integrated Co-Fermentation'],
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[letter_color, letter_color],
                 [letter_color, letter_color],
                 [letter_color, letter_color]]
    _add_letter_labels(axes, 0.7, 0.75, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'ethanol_and_biodiesel_price_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_enhanced_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        enhanced_cellulosic_performance=True,
        titles=['Direct Cogeneration', 'Integrated Co-Fermentation'],
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[letter_color, letter_color],
                 [letter_color, letter_color],
                 [letter_color, letter_color]]
    _add_letter_labels(axes, 0.7, 0.75, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'enhanced_ethanol_and_biodiesel_price_contours.{i}')
        plt.savefig(file, transparent=True)
    
def plot_benefit_ethanol_and_biodiesel_price_contours_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.7)
    fig, axes = plot_ethanol_and_biodiesel_price_contours(
        benefit=True,
        titles=['Direct Cogeneration', 'Integrated Co-Fermentation'],
        dist=True,
    )
    colors = np.zeros([3, 2], object)
    colors[:] = [[dark_letter_color, dark_letter_color],
                 [dark_letter_color, dark_letter_color],
                 [dark_letter_color, dark_letter_color]]
    _add_letter_labels(axes, 0.7, 0.75, colors)
    plt.subplots_adjust(left=0.2, right=0.92, wspace=0.1, top=0.94, bottom=0.12)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'benefit_ethanol_and_biodiesel_price_contours.{i}')
        plt.savefig(file, transparent=True)

# %% General    

def plot_ethanol_and_biodiesel_price_contours(N=30, benefit=False, cache={}, 
                                              enhanced_cellulosic_performance=False,
                                              titles=None, load=False, save=True,
                                              dist=False):
    ethanol_price = np.linspace(0.25, 0.5, N)
    biodiesel_price = np.linspace(0.65, 0.95, N)
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
    if load and (N, benefit, enhanced_cellulosic_performance) in cache:
        Z = cache[N, benefit, enhanced_cellulosic_performance]
    else:
        Z = np.zeros([N, N, N_rows, N_cols])
        for i in range(N_rows):
            for j in range(N_cols):
                oc.load(configuration[j],
                        reduce_chemicals=False,
                        enhanced_cellulosic_performance=enhanced_cellulosic_performance)
                oc.set_cane_oil_content(oil_content[i])
                oc.set_relative_sorghum_oil_content(0)
                oc.sys.simulate()
                if benefit:
                    Z[:, :, i, j] = oc.evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices(X, Y)
                else:
                    Z[:, :, i, j] = oc.evaluate_MFPP_across_ethanol_and_biodiesel_prices(X, Y)
    if save:
        np.save(file, Z)
    xlabel = f"Ethanol price\n[{format_units('$/L')}]"
    ylabels = [f"Biodiesel price\n[{format_units('$/L')}]"] * 4
    xticks = [0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
    yticks = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    marks = tickmarks(Z, 5, 5, center=0.) if benefit else tickmarks(Z, 5, 5)
    colormap = (diverging_colormaps if benefit else colormaps)[0]
    name = 'MFPP'
    if benefit: name = r'$\Delta$MFPP'
    metric_bar = MetricBar(name, format_units('$/MT'), colormap, marks, 15,
                           forced_size=0.3)
    if benefit:
        baseline = ['S1', 'S2']
        if titles is None:
            titles = [f"{oc.format_name(i)} - {oc.format_name(j)}"
                      for i, j in  zip(configuration, baseline)]
    elif titles is None:
        titles = [oc.format_name(i) for i in configuration]
    fig, axes, cps, cb, other_axes = plot_contour_single_metric(
        X, Y, Z, xlabel, ylabels, xticks, yticks, metric_bar, 
        styleaxiskw=dict(xtickf=False, ytickf=False), label=True,
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
        plot_scatter_points(ethanol_no_RIN_prices, biodiesel_prices, 
                            marker='o', s=2, color=(*colors.brown.RGBn, 1),
                            edgecolor=(*colors.brown.RGBn, 1), clip_on=True, zorder=3)
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
        load=False, configuration_index=..., relative=False, smooth=None,
    ):
    # Generate contour data
    X, Y, data = relative_sorghum_oil_content_and_cane_oil_content_data(load, relative)
    
    data = data[:, :, configuration_index, [0, 6]]
    
    # Plot contours
    xlabel = "Oil-sorghum oil content [dry wt. %]" 
    if relative: xlabel = ('relative ' + xlabel).capitalize()
    ylabel = 'Oilcane oil content\n[dry wt. %]'
    yticks = [5, 7.5, 10, 12.5, 15]
    xticks = [-3, -2, -1, 0] if relative else [2, 5, 7.5, 10, 12.5, 15]
    MFPP = oc.all_metric_mockups[0]
    TCI = oc.all_metric_mockups[6]
    if configuration_index == 0:
        Z = np.array(["Direct Cogeneration"])
        data = data[:, :, :, np.newaxis]
    elif configuration_index == 1:
        Z = np.array(["Integrated Co-Fermentation"])
        data = data[:, :, :, np.newaxis]
    elif configuration_index == ...:
        Z = np.array(["Direct Cogeneration", "Integrated Co-Fermentation"])
        data = np.swapaxes(data, 2, 3)
    else:
        raise ValueError('configuration index must be either 0 or 1')
    metric_bars = [
        [MetricBar(MFPP.name, format_units(MFPP.units), colormaps[0], tickmarks(data[:, :, 0, 0], 5, 1, expand=0, p=0.5), 10, 1),
         MetricBar(MFPP.name, format_units(MFPP.units), colormaps[0], tickmarks(data[:, :, 0, 1], 5, 1, expand=0, p=0.5), 10, 1)],
        [MetricBar(TCI.name, format_units(TCI.units), colormaps[1], tickmarks(data[:, :, 1, 0], 5, 5, expand=0, p=5), 10, 1),
         MetricBar(TCI.name, format_units(TCI.units), colormaps[1], tickmarks(data[:, :, 1, 1], 5, 5, expand=0, p=5), 10, 1)],
    ]
    
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, M, N = data.shape
        for m in range(M):
            for n in range(N):
                metric_data = data[:, :, m, n]
                data[:, :, m, n] = gaussian_filter(metric_data, smooth)
                # for a in range(A):
                #     values = metric_data[a, :]
                #     values.sort()
                #     p = np.arange(values.size)
                #     coeff = np.polyfit(p, values, 5)
                #     values[:] = np.polyval(coeff, p)
                # for b in range(B):
                #     values = metric_data[:, b]
                #     values.sort()
                #     p = np.arange(values.size)
                #     coeff = np.polyfit(p, values, 5)
                #     values[:] = np.polyval(coeff, p)
    
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X, 100.*Y, Z, data, xlabel, ylabel, xticks, yticks, metric_bars, 
        styleaxiskw=dict(xtick0=True), label=True,
    )
    for i in axes.flatten():
        plt.sca(i)
        plot_scatter_points([7], [10], marker='*', s=100, color=startcolor,
                            edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
    
def plot_recovery_and_oil_content_contours(
        load=False, metric_index=0, N_decimals=1, configurations=None,
        N_points=20, yticks=None, titles=None, cmap=None, smooth=None,
    ):
    if yticks is None: yticks = [5, 7.5, 10, 12.5, 15]
    if configurations is None:
        configurations = np.array([['O1', 'O1*'], ['O2', 'O2*']])
        if titles is None: titles = ['Oilcane Only', 'Oilcane & Oil-sorghum']
        
    # Generate contour data
    x = np.linspace(0.40, 1.0, N_points)
    y = np.linspace(yticks[0] / 100, yticks[-1] / 100, N_points)
    X, Y = np.meshgrid(x, y)
    metric = bst.metric
    folder = os.path.dirname(__file__)
    file = "oil_extraction_analysis.npy"
    file = os.path.join(folder, file)
    
    if load:
        data = np.load(file)
    else:
        data = oc.evaluate_configurations_across_recovery_and_oil_content(
            X, Y, configurations,
        )
    np.save(file, data)
    data = data[:, :, :, :, metric_index]
    # data = np.swapaxes(data, 2, 3)
    
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, M, N = data.shape
        for m in range(M):
            for n in range(N):
                metric_data = data[:, :, m, n]
                data[:, :, m, n] = gaussian_filter(metric_data, smooth)
                # for a in range(A):
                #     values = metric_data[a, :]
                #     # values.sort()
                #     p = np.arange(values.size)
                #     coeff = np.polyfit(p, values, 3)
                #     values[:] = np.polyval(coeff, p)
                # for b in range(B):
                #     values = metric_data[:, b]
                #     # values.sort()
                #     p = np.arange(values.size)
                #     coeff = np.polyfit(p, values, 3)
                #     values[:] = np.polyval(coeff, p)
    
    # Plot contours
    xlabel = 'Crushing mill oil recovery [%]'
    ylabel = "Oil content [dry wt. %]"
    ylabels = [f'Direct Cogeneration\n{ylabel}',
               f'Integrated Co-Fermentation\n{ylabel}']
    xticks = [40, 50, 60, 70, 80, 90, 100]
    
    metric = oc.all_metric_mockups[metric_index]
    units = metric.units if metric.units == '%' else format_units(metric.units)
    if cmap is None: cmap = colormaps[metric_index]
    # if metric_index == 10:
    #     breakpoint()
    mb = lambda x, name=None: MetricBar(
        metric.name if name is None else name, units if name is None else "", cmap, 
        tickmarks(data[:, :, x, :], 5, 0.1, expand=0, p=0.1, f=lambda x: round(x, N_decimals)),
        10, N_decimals=N_decimals
    )
    
    metric_bars = [mb(0), mb(1, "")]
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X, 100.*Y, titles, data, xlabel, ylabels, xticks, yticks, metric_bars, 
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    M, N = configurations.shape
    for i in range(N):
        for j in range(M):
            ax = axes[i, j]
            plt.sca(ax)
            plt.fill_between([60, 90], [yticks[0]], [yticks[-1]], 
                              color=shadecolor,
                              linewidth=1)
            plot_vertical_line(60, ls='-.',
                                color=linecolor,
                                linewidth=1.0)
            plot_vertical_line(90, ls='-.',
                               color=linecolor,
                               linewidth=1.0)
            if hasattr(ax, '_cached_ytwin'):                
                plt.sca(ax._cached_ytwin)
            plot_scatter_points([60], [10], marker='*', s=100, color=startcolor,
                                edgecolor=edgecolor, clip_on=False, zorder=3)
            # plot_scatter_points([ub], [15], marker='*', s=100, color=targetcolor,
            #                     edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
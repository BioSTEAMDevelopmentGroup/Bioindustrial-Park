# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 23:44:10 2021

@author: yrc2
"""
from biorefineries import cane
import biosteam as bst
import numpy as np
from biosteam.utils import colors, GG_colors
import matplotlib.pyplot as plt
from biosteam.plots import (
    rounded_linspace,
    generate_contour_data,
    plot_contour_2d,
    MetricBar,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks,
    plot_scatter_points,
    plot_contour_single_metric,
)
from thermosteam.units_of_measure import format_units
from thermosteam.utils import set_figure_size, set_font
from .results import images_folder
from . import feature_mockups as feature
from .data import microbial_oil_baseline as perf
from warnings import filterwarnings
from scipy.ndimage.filters import gaussian_filter
from colorpalette import ColorWheel
from scipy.interpolate import LinearNDInterpolator
import os

__all__ = (
    'plot_sorghum_oil_content_and_cane_oil_content_contours_manuscript',
    'plot_metrics_across_composition_manuscript',
    'plot_metrics_across_biomass_yield_manuscript',
    'plot_sorghum_oil_content_and_cane_oil_content_contours_seminar',
    'plot_recovery_and_oil_content_contours_manuscript',
    'plot_oil_recovery_integration_manuscript',
    'plot_recovery_and_oil_content_contours',
    'plot_sorghum_oil_content_and_cane_oil_content_contours',
    'plot_recovery_and_oil_content_contours_biodiesel_only',
    'plot_recovery_and_oil_content_contours_with_oilsorghum_only',
    'plot_metrics_across_biomass_yield',
    'plot_metrics_across_composition',
    'plot_oil_recovery_integration',
    'plot_lines_biomass_yield',
)

filterwarnings('ignore', category=bst.exceptions.DesignWarning)
line_color_wheel = ColorWheel([
    GG_colors.orange,
    GG_colors.blue,
    GG_colors.purple,
    GG_colors.green,
    GG_colors.yellow,
    colors.CABBI_brown,
    colors.CABBI_grey,
    colors.CABBI_teal,
])
for i, color in enumerate(line_color_wheel.colors): line_color_wheel.colors[i] = color.tint(20)
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
default_box_color = colors.neutral.RGBn

# %% Plot functions for publication

def _add_letter_labels(axes, xpos, ypos, colors, box_color=None):
    M, N = colors.shape
    if box_color is None:
        box_color = default_box_color
    txtbox = dict(boxstyle='round', facecolor=box_color, 
                  edgecolor='None', alpha=0.6)
    letters = 'ABCDEFG'
    n = -1
    for i in range(M):
        for j in range(N):
            n += 1
            letter = letters[n]
            ax = axes[i, j]
            xlb, xub = ax.get_xlim()
            ylb, yub = ax.get_ylim()
            if hasattr(ax, '_cached_ytwin'):                
                ax = ax._cached_ytwin
            ax.text((xlb + xub) * xpos, (yub + ylb) * ypos, letter, color=colors[i, j],
                     horizontalalignment='center',verticalalignment='center',
                     fontsize=10, fontweight='bold', zorder=1, bbox=txtbox)

def plot_metrics_across_composition_manuscript(load=True, fs=8, smooth=1):
    set_font(size=fs)
    width = None
    aspect_ratio = 0.5
    #x = 0.58
    set_figure_size(width=width, aspect_ratio=aspect_ratio)
    fig, axes = plot_metrics_across_composition(
        load=load, 
        smooth=smooth,
    )
    colors = np.zeros([1, 2], object)
    colors[:] = [[light_letter_color, light_letter_color]]
    # _add_letter_labels(axes, 1 - 0.9, 0.62, colors)
    plt.subplots_adjust(left=0.1, right=0.95, wspace=0.15, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'metrics_across_composition_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_metrics_across_biomass_yield_manuscript(load=True, fs=8, smooth=1):
    set_font(size=fs)
    width = None
    aspect_ratio = 0.8
    #x = 0.58
    set_figure_size(width=width, aspect_ratio=aspect_ratio)
    fig, axes = plot_metrics_across_biomass_yield(
        load=load, 
        smooth=smooth,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 0.9, 0.8, colors)
    plt.subplots_adjust(left=0.1, right=0.95, wspace=0.15, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'metrics_biomass_yield_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_oil_recovery_integration_manuscript(load=True, fs=8):
    set_font(size=fs)
    width = None
    aspect_ratio = 0.8
    #x = 0.58
    set_figure_size(width=width, aspect_ratio=aspect_ratio)
    fig, axes, other_axes = plot_oil_recovery_integration(
        load=load, 
    )
    colors = np.zeros([2, 3], object)
    colors[:] = [[light_letter_color, light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 0.62, 0.62, colors)
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.15, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oil_recovery_integration_contours.{i}')
        plt.savefig(file, transparent=True)

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

def plot_recovery_and_oil_content_contours_with_oilsorghum_only(fs=10, smooth=1):
    set_font(size=fs)
    set_figure_size(4, 1.1)
    fig, axes = plot_recovery_and_oil_content_contours(
        load=True, 
        smooth=smooth,
        with_oilsorghum_only=True
    )
    colors = np.zeros([2, 1], object)
    colors[:] = [[light_letter_color], [light_letter_color]]
    # _add_letter_labels(axes, 1 - 0.68, 0.7, colors)
    plt.subplots_adjust(left=.2, right=0.92, wspace=0.1 * (fs/8) ** 2, top=0.9, bottom=0.10)
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

def plot_sorghum_oil_content_and_cane_oil_content_contours_manuscript(load=True, fs=8, smooth=0.9, configuration_index=None):
    set_font(size=fs)
    set_figure_size()
    fig, axes = plot_sorghum_oil_content_and_cane_oil_content_contours(
        load=load, smooth=smooth, configuration_index=configuration_index,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 1 - 0.75, 0.75, colors)
    plt.subplots_adjust(left=0.12, right=0.92, wspace=0.12, hspace=0.12, top=0.9, bottom=0.12)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'relative_sorghum_oil_content_and_cane_oil_content_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_sorghum_oil_content_and_cane_oil_content_contours_seminar(load=True):
    plot_sorghum_oil_content_and_cane_oil_content_contours_manuscript(load, 10, 1.3, [2, 6])

# %% General    

def plot_metrics_across_composition(
        load=False, N_decimals=1, 
        yticks=None, titles=None, 
        cmap=None, smooth=None,
    ):
    xticks = [0,   2,   4,   6,   8,   10]
    xlim = np.array([xticks[0], xticks[-1]])
    yticks = [35,  45,  55,  65]
    ylim = np.array([yticks[0], yticks[-1]])
    metrics = CBY, = [cane.competitive_biomass_yield]
    # metrics = CBY, NEP = [cane.competitive_biomass_yield, cane.net_energy_production]
    metric_indices = [cane.all_metric_mockups.index(i) for i in metrics] 
    cane.YRCP2023()
    X, Y, data = generate_contour_data(
        cane.evaluate_metrics_at_composition, 
        file=contour_file('configurations_compositional_analysis'),
        xlim=xlim / 100, ylim=ylim / 100, load=load, n = 10,
    )
    data = data[:, :, metric_indices, :]
    # data[:, 0, 0, 0] = data[:, 1, 0, 0] + (data[:, 1, 0, 0] - data[:, 2, 0, 0])
    # data[:, 0, 1, 0] = data[:, 1, 1, 0] + (data[:, 1, 1, 0] - data[:, 2, 1, 0])
    
    # Plot contours
    xlabel = 'Oil content [dw %]'
    ylabel = "Fiber content [dw %]"
    if titles is None: titles = np.array(['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'])
    d0 = data[:, :, 0, :]
    # d1 = data[:, :, 1, :]
    CBY_units = format_units('DMT/ha/y')
    metric_bars = [
        # MetricBar(MFPP.name, format_units(MFPP.units), colormaps[1], tickmarks(data[:, :, 0, :], 5, 1, expand=0, p=0.5), 18, 1),
        MetricBar('Competitive biomass yield', CBY_units, 
                  plt.cm.get_cmap('viridis_r'), 
                  tickmarks(d0[~np.isnan(d0)], 5, 1, expand=0, p=0.5),
                  15, 1),
        # MetricBar('Biod. prod.', format_units(BP.units), plt.cm.get_cmap('copper'), tickmarks(data[:, :, 1, :], 8, 1, expand=0, p=0.5), 10, 1),
        # MetricBar('MBSP', format_units(MBSP.units), plt.cm.get_cmap('inferno_r'), tickmarks(d1[~np.isnan(d1)], 5, 0.1, expand=0, p=0.5), 20, 1),
        # MetricBar('Net energy prod.', format_units(NEP.units), plt.cm.get_cmap('inferno'), tickmarks(d1[~np.isnan(d1)], 5, 1, expand=0, p=2), 20, 1),
    ]
    fig, axes, CSs, CB, other_axes = plot_contour_2d(
        100.*X, 100.*Y, data, xlabel, ylabel, xticks, yticks, 
        metric_bars=metric_bars, titles=titles, 
        styleaxiskw=dict(xtick0=True), label=True, wbar=2.1
    )
    try:
        df = cane.get_composition_data()
    except:
        pass
    else:
        names = df.index
        lines = []
        for name, color in zip(names, line_color_wheel):
            data = df.loc[name]
            name = str(name)
            if name in ('19B', '316', '1565'): continue
            lines.append(
                (name, 
                 data['Stem oil (dw)']['Mean'] * 100, 
                 data['Fiber (dw)']['Mean'] * 100,
                 data['Water (wt)']['Mean'] * 100,
                 data['Biomass yield (dry MT/ha)']['Mean'],
                 color.RGBn)
            )
        for i in axes.flat: 
            if hasattr(i, '_cached_ytwin'):
                plt.sca(i); plt.xlim(xlim)
                plt.sca(i._cached_ytwin); plt.xlim(xlim)
        txtbox = dict(boxstyle='round', facecolor=colors.neutral.shade(20).RGBn, 
                      edgecolor='None', alpha=0.99, pad=0.1)
        
        for i, (*axes_columns, _) in enumerate(axes):
            for (name, lipid, fiber, moisture, biomass_yield, color) in lines:
                for j, ax in enumerate(axes_columns):
                    plt.sca(ax._cached_ytwin)
                    feedstock = 'Sugarcane' if name == 'WT' else 'Oilcane'
                    if name == '1580':
                        text = f"{feedstock} {name}:\n{biomass_yield:.0f} {CBY_units}"
                        txtbox_1580 = txtbox.copy()
                        txtbox_1580['pad'] = 0.3
                        plt.text(
                            lipid + 0.15, fiber, text, c=color,
                            verticalalignment='bottom',
                            bbox=txtbox_1580, 
                        )
                    else:
                    
                        text = f"{feedstock} {name}: {biomass_yield:.0f} {CBY_units}"
                        plt.text(
                            lipid + 0.15, fiber, text, c=color,
                            verticalalignment='bottom',
                            bbox=txtbox, 
                        )
                    plot_scatter_points(
                        [lipid], [fiber], marker='o', s=20,  linewidth=0.6,
                        color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
                    )
                    
    return fig, axes

def plot_lines_biomass_yield():
    fig, ax = plt.subplots(nrows=1, ncols=1)
    df = cane.get_composition_data(ignored=set())
    names = df.index
    key_lines = {'WT', 1566, 1580}
    names = sorted(names, key=lambda x: x not in key_lines)
    lines = []
    for name, color in zip(names, line_color_wheel):
        data = df.loc[name]
        name = str(name)
        lines.append(
            (name, 
             data['Stem oil (dw)']['Mean'] * 100, 
             data['Fiber (dw)']['Mean'] * 100,
             data['Water (wt)']['Mean'] * 100,
             data['Biomass yield (dry MT/ha)']['Mean'],
             color)
        )
    CBY_units = format_units('DMT/ha/y')
    oil_units = 'dw %'
    default_options = dict(
        dx=0, dy=0, 
        dx_text=0.05, dy_text=0, 
        horizontalalignment='left', 
        verticalalignment='bottom',
    )
    options = {
        '19B': {'dx': 0.3, 'dy': 3, 
                'dy_text': 0.05, 'dx_text': 0.02,
                'verticalalignment': 'center'},
        '1566': {'verticalalignment': 'top'},
        '1580': {'horizontalalignment':'right', 'dx_text':-0.05},
    }
    lipids = []
    biomass_yields = []
    for (name, lipid, fiber, moisture, biomass_yield, color) in lines:
        if name == 'WT':
            feedstock = 'Sugarcane'
            x0 = lipid
            y0 = biomass_yield
        else:
            feedstock = 'Oilcane'
        text = (
            f"{feedstock} {name}: {lipid:.0f} {oil_units}, "
            f"{biomass_yield:.0f} {CBY_units}"
        )
        bst.plots.annotate_point(
            text, lipid, biomass_yield,
            textcolor=color.shade(50).RGBn,
            linecolor=color.shade(50).RGBn,
            fontsize=10,
            **(default_options | options.get(name, {})), 
        )    
        plot_scatter_points(
            [lipid], [biomass_yield], marker='o', s=30,  linewidth=0.6,
            color=color.RGBn, edgecolor=color.shade(50).RGBn,
            clip_on=False, zorder=1e6,
        )
        lipids.append(lipid)
        biomass_yields.append(biomass_yield)
    
    # # Linear regression that passes through point WT
    # lipids = np.array(lipids)
    # biomass_yields = np.array(biomass_yields)
    # x = lipids - x0
    # y = biomass_yields - y0
    # m = (y * x).sum() / (x * x).sum()
    # b = y0 - m * x0
    # plt.plot(lipids, m * lipids + b, color=colors.CABBI_grey.shade(20).RGBn, 
    #          zorder=-1e9)
    plt.ylim([0, 30])
    plt.yticks([5, 10, 15, 20, 25, 30])
    plt.xlim([0, 6])
    plt.xlabel(f'Oil content [{oil_units}]')
    plt.ylabel(f'Biomass yield [{CBY_units}]')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'lines_oil_content_and_biomass_yield.{i}')
        plt.savefig(file, transparent=True)

def plot_metrics_across_biomass_yield(
        load=False, N_decimals=1, 
        yticks=None, titles=None, 
        cmap=None, smooth=None,
    ):
    xticks = [0,   2,   4,   6,   8,   10]
    xlim = (xticks[0], xticks[-1])
    yticks = [4, 10, 16, 22, 28]
    ylim = (yticks[0], yticks[-1])
    metrics = MBSP, GWP = [cane.MBSP, cane.GWP_biofuel_allocation]
    metric_indices = [cane.all_metric_mockups.index(i) for i in metrics] 
    cane.YRCP2023()
    X, Y, data = generate_contour_data(
        cane.evaluate_metrics_at_biomass_yield,
        xlim, ylim, file=contour_file('configurations_biomass_yield_analysis'),
        smooth=smooth,
        load=load,
        n=10,
    )
    data = data[:, :, metric_indices, :]
    # Plot contours
    xlabel = 'Oil content [dw %]'
    BY_units = format_units('DMT/ha/y')
    ylabel = f"Biomass yield [{BY_units}]"
    if titles is None: titles = np.array(['Direct Cogeneration', 'Integrated Co-Fermentation'])
    d0 = data[:, :, 0, :]
    metric_bars = [
        MetricBar('MBSP', format_units('USD/L'), plt.cm.get_cmap('viridis_r'), tickmarks(d0[~np.isnan(d0)], 5, 1, expand=0, p=0.5), 15, 1),
        MetricBar('GWP', format_units(GWP.units).replace('CO2', 'CO_2'), plt.cm.get_cmap('inferno_r'), tickmarks(d0[~np.isnan(d0)], 5, 0.1, expand=0, p=0.1), 15, 1),
    ]
    fig, axes, CSs, CB, other_axes = plot_contour_2d(
        X, Y, data, xlabel, ylabel, xticks, yticks, 
        metric_bars=metric_bars, titles=titles,
        styleaxiskw=dict(xtick0=True), 
        label=True, wbar=2.1
    )
    try:
        df = cane.get_composition_data()
    except:
        pass
    else:
        names = df.index
        lines = []
        for name, color in zip(names, line_color_wheel):
            data = df.loc[name]
            if str(name) in ('19B', '316', '1565'): continue
            oil = data['Stem oil (dw)']['Mean'] * 100
            lines.append(
                (name, 
                 oil, 
                 data['Biomass yield (dry MT/ha)']['Mean'],
                 color.RGBn)
            )
        for i in axes.flat: 
            if hasattr(i, '_cached_ytwin'):
                plt.sca(i); plt.xlim(xlim)
                plt.sca(i._cached_ytwin); plt.xlim(xlim)
        txtbox = dict(boxstyle='round', facecolor=colors.neutral.shade(20).RGBn, 
                      edgecolor='None', alpha=0.999, pad=0.3)
        
        for (name, lipid, biomass_yield, color) in lines:
            feedstock = 'Sugarcane' if name == 'WT' else 'Oilcane'
            line_results = cane.evaluate_metrics_at_biomass_yield(lipid, biomass_yield)[metric_indices]
            for i, (*axes_columns, _) in enumerate(axes):
                mb = metric_bars[i]
                units = mb.units.replace('CO2', 'CO_2')
                x = lipid + 0.2
                y = biomass_yield
                if biomass_yield > 25:
                    verticalalignment = 'top'
                elif biomass_yield < 10:
                    verticalalignment = 'bottom'
                else:
                    verticalalignment = 'center'
                    
                for j, ax in enumerate(axes_columns):
                    plt.sca(ax._cached_ytwin)
                    value = line_results[i, j]
                    text = f"{feedstock} {name}:\n{value:.1f} {units}"
                    plt.text(
                        x, y, text, c=color,
                        bbox=txtbox, verticalalignment=verticalalignment,
                    )
                    plot_scatter_points(
                        [lipid], [biomass_yield], marker='o', s=20,  linewidth=0.6,
                        color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
                    )
                    
    return fig, axes

def relative_sorghum_oil_content_and_cane_oil_content_data(load, configurations):
    # Generate contour data
    y = np.linspace(0.01, 0.10, 20)
    x = np.linspace(0.01, 0.10, 20)
    X, Y = np.meshgrid(x, y)
    folder = os.path.dirname(__file__)
    config = ''.join([str(i) for i in configurations])
    file = f'oil_content_analysis_{config}.npy'
    file = os.path.join(folder, file)
    if configurations is None: configurations = [1, 2]
    if load:
        data = np.load(file)
    else:
        from warnings import filterwarnings; filterwarnings('ignore')
        data = cane.evaluate_configurations_across_sorghum_and_cane_oil_content(
            X, Y, configurations, 
        )
        np.save(file, data)
    return X, Y, data
    
def contour_file(name):
    folder = os.path.dirname(__file__)
    file = name + '.npy'
    return os.path.join(folder, file)

# def plot_fermentation_performance(
#         load=False, 
#     ):
#     metrics = [feature.MBSP, feature.GWP_biodiesel]
#     metrics_index = [feature.all_metric_mockups.index(i) for i in metrics]
#     X, Y, Z = generate_contour_data(
#         cane.evaluate_fermentation_performance, # TODO
#         xlim=[50, 90],
#         ylim=[100 * perf.min_lipid_yield_glucose, 100 * perf.max_lipid_yield_glucose],
#         file=contour_file('oil_recovery_integration'),
#         n=10,
#         load=load,
#     )
#     Z = Z[:, :, metrics_index, :]
#     # Plot contours
#     xlabel = 'Microbial oil recovery [%]'
#     ylabel = "Microbial oil yield [wt %]"
#     xticks = [50, 60, 70, 80, 90]
#     yticks = rounded_linspace(
#         100 * perf.min_lipid_yield_glucose, 100 * perf.max_lipid_yield_glucose, 5, 
#     )
#     titles = ['Mechanical oil recovery', 'Integrated oil recovery']
    
#     d0 = Z[:, :, 0, :]
#     d1 = Z[:, :, 1, :]
#     d2 = Z[:, :, 2, :]
#     metric_bars = [
#         MetricBar('MBSP', format_units('USD/L'), plt.cm.get_cmap('viridis_r'), tickmarks(d0[~np.isnan(d0)], 5, 1, expand=0, p=0.5), 15, 1),
#         MetricBar('ROI', format_units('% USD'), plt.cm.get_cmap('viridis'), tickmarks(d1[~np.isnan(d1)], 5, 1, expand=0, p=1), 15, 1),
#         MetricBar("GWP", format_units('USD/L'), plt.cm.get_cmap('inferno_r'), tickmarks(d2[~np.isnan(d2)], 5, 0.1, expand=0, p=0.1), 15, 1),
#     ]
#     fig, axes, CSs, CB, other_axes = plot_contour_2d(
#         X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  titles,
#         fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
#     )
#     return fig, other_axes

def plot_oil_recovery_integration(
        load=False, metric=None,
    ):
    productivity = np.array([
        perf.hydrolysate_productivity, 
        perf.batch_productivity_mean, 
    ])
    titer = np.array([
        perf.hydrolysate_titer, 
        perf.batch_titer_mean, 
    ])
    cane.YRCP2023()
    X, Y, Z = generate_contour_data(
        cane.evaluate_metrics_oil_recovery_integration,
        xlim=[50, 90],
        ylim=[100 * perf.min_lipid_yield_glucose, 100 * perf.fed_batch_lipid_yield_mean],
        args=(
            productivity, titer,
        ),
        file=contour_file('oil_recovery_integration'),
        load=load,
        n=10,
    )
    if metric is None:
        metric = feature.MBSP
    elif isinstance(metric, str):
        metric = getattr(feature, metric)
    metric_index = cane.all_metric_mockups.index(metric)
    Z = Z[..., metric_index]
    # Plot contours
    xlabel = 'Microbial oil recovery [%]'
    ylabel = "Microbial oil yield [wt %]"
    titer_units = r'$g \cdot L^{-1}$'
    productivity_units = r'$g \cdot L^{-1} \cdot h^{-1}$'
    ylabels = [f"{ylabel}\n"
               f"productivity = {round(productivity[i], 2)} {productivity_units}\n"
               f"titer = {round(titer[i], 2)} {titer_units}"
               for i in range(len(productivity))]
    xticks = [50, 60, 70, 80, 90]
    yticks = [12, 14, 16, 18, 20, 22]
    titles = ['Direct Cogeneration',
              'Integrated\nCo-Fermentation', 
              'Integrated Co-Fermentation\n& Recovery']
    units = format_units(metric.units)
    metric_bar = MetricBar(
        metric.name, units, plt.cm.get_cmap('viridis_r'), tickmarks(Z[~np.isnan(Z)], 5, 0.5, expand=0, p=0.5), 25, 2
    )
    fig, axes, CSs, CB, other_axes = plot_contour_single_metric(
        X, Y, Z, xlabel, ylabels, xticks, yticks, metric_bar,  titles,
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    txtbox = dict(boxstyle='round', facecolor=colors.neutral.shade(20).RGBn, 
                  edgecolor='None', alpha=0.99, pad=0.1)
    color = GG_colors.orange.RGBn
    x = np.array([70, 70, 50])
    y = 100 * perf.hydrolysate_lipid_yield
    Z0 = Z[..., 0, :]
    new_shape = [Z0.shape[0] * Z0.shape[1], Z0.shape[-1]]
    interp = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), Z0.reshape(new_shape))
    baseline_values = interp(x, y)
    for i, ax in enumerate(axes[0, :3]):
        plt.sca(ax._cached_ytwin)
        value = baseline_values[i, i]
        text = f'Baseline:\n{value:.2f} {units}'
        plt.text(
            x[i] + 0.9, y + 0.7, text, c=color,
            verticalalignment='center',
            bbox=txtbox, 
        )
        plot_scatter_points(
            [x[i]], [y], 
            marker='o', s=20, linewidth=0.6,
            color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
        )
    color = GG_colors.yellow.RGBn
    y = 100 * perf.batch_lipid_yield_mean
    Z1 = Z[..., 1, :]
    interp = LinearNDInterpolator(list(zip(X.flatten(), Y.flatten())), Z1.reshape(new_shape))
    target_values = interp(x, y)
    for i, ax in enumerate(axes[1, :3]):
        plt.sca(ax._cached_ytwin)
        value = target_values[i, i]
        text = f'Target:\n{value:.2f} {units}'
        plt.text(
            x[i] + 0.9, y + 0.7, text, c=color,
            verticalalignment='center',
            bbox=txtbox, 
        )
        plot_scatter_points(
            [x[i]], [y], 
            marker='o', s=20, linewidth=0.6,
            color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
        )
    return fig, axes, other_axes

def plot_sorghum_oil_content_and_cane_oil_content_contours(
        load=False, configuration_index=None, smooth=None,
    ):
    if configuration_index is None: configuration_index = [1, 2]
    # Generate contour data
    X, Y, data = relative_sorghum_oil_content_and_cane_oil_content_data(load, configuration_index)
    MFPP = cane.MFPP
    TCI = cane.TCI
    metrics = [MFPP, TCI]
    metric_indices = [cane.all_metric_mockups.index(i) for i in metrics]
    data = data[:, :, :, metric_indices]
    
    # Plot contours
    xlabel = "Oil-sorghum oil content [dw %]" 
    ylabel = 'Oilcane oil content\n[dw %]'
    yticks = [2, 4, 6, 8, 10]
    xticks = [2, 4, 6, 8, 10]
    if configuration_index == [1, 2]:
        Z = np.array(["Direct Cogeneration", "Integrated Co-Fermentation"])
        data = np.swapaxes(data, 2, 3)
    elif configuration_index == [2, 6]:
        Z = np.array(["Bioethanol", "Microbial oil"])
        data = np.swapaxes(data, 2, 3)
    else:
        raise ValueError('configuration index must be either [1, 2] or [2, 6]')
    
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, M, N = data.shape
        for m in range(M):
            for n in range(N):
                metric_data = data[:, :, m, n]
                if isinstance(smooth, int):
                    for a in range(A):
                        values = metric_data[a, :]
                        values.sort()
                        p = np.arange(values.size)
                        coeff = np.polyfit(p, values, smooth)
                        values[:] = np.polyval(coeff, p)
                    for b in range(B):
                        values = metric_data[:, b]
                        values.sort()
                        p = np.arange(values.size)
                        coeff = np.polyfit(p, values, smooth)
                        values[:] = np.polyval(coeff, p)
                else:
                    data[:, :, m, n] = gaussian_filter(metric_data, smooth)
                    
    metric_bars = [
        # [MetricBar(MFPP.name, format_units(MFPP.units), colormaps[0], tickmarks(data[:, :, 0, 0], 5, 1, expand=0, p=0.5), 10, 1),
         MetricBar(MFPP.name, format_units(MFPP.units), colormaps[0], tickmarks(data[:, :, 0, [0, 1]], 5, 1, expand=0, p=0.5), 12, 1),
        # [MetricBar(TCI.name, format_units(TCI.units), colormaps[1], tickmarks(data[:, :, 1, 0], 5, 5, expand=0, p=5), 10, 1),
         MetricBar(TCI.name, format_units(TCI.units), colormaps[1], tickmarks(data[:, :, 1, [0, 1]], 5, 5, expand=0, p=5), 15, 1),
    ]
    
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X, 100.*Y, Z, data, xlabel, ylabel, xticks, yticks, metric_bars, 
        styleaxiskw=dict(xtick0=True), label=True,
    )
    for *cols, _ in axes:
        for i in cols:
            plt.sca(i)
            plt.xlim([1, 10])
    # for i in axes.flatten():
    #     plt.sca(i)
    #     plot_scatter_points([7], [10], marker='*', s=100, color=startcolor,
    #                         edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
    
def plot_recovery_and_oil_content_contours(
        load=False, metric_index=0, N_decimals=1, configurations=None,
        N_points=20, yticks=None, titles=None, cmap=None, smooth=None,
        with_oilsorghum_only=False,
    ):
    if yticks is None: yticks = [5, 7.5, 10, 12.5, 15]
    if configurations is None:
        if with_oilsorghum_only:
            configurations = np.array([['O1*'], ['O2*']])
            if titles is None: titles = ['Oilcane & Oil-sorghum']
        else:
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
        data = cane.evaluate_configurations_across_recovery_and_oil_content(
            X, Y, configurations,
        )
    np.save(file, data)
    data = data[:, :, :, :, metric_index]
    
    if with_oilsorghum_only:
        data = data[:, :, :, -1:]
    
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, M, N = data.shape
        for m in range(M):
            for n in range(N):
                metric_data = data[:, :, m, n]
                data[:, :, m, n] = gaussian_filter(metric_data, smooth)
    
    # Plot contours
    xlabel = 'Crushing mill oil recovery [%]'
    ylabel = "Oil content [dw %]"
    ylabels = [f'Direct Cogeneration\n{ylabel}',
               f'Integrated Co-Fermentation\n{ylabel}']
    xticks = [40, 50, 60, 70, 80, 90, 100]
    
    metric = cane.all_metric_mockups[metric_index]
    units = metric.units if metric.units == '%' else format_units(metric.units)
    if cmap is None: cmap = colormaps[metric_index]
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
    for i in range(M):
        for j in range(N):
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
            # plot_scatter_points([60], [10], marker='*', s=100, color=startcolor,
            #                     edgecolor=edgecolor, clip_on=False, zorder=3)
            # plot_scatter_points([ub], [15], marker='*', s=100, color=targetcolor,
            #                     edgecolor=edgecolor, clip_on=False, zorder=3)
    return fig, axes
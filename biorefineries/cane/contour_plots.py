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
    plot_contour_2d,
    MetricBar,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks,
    plot_scatter_points,
)
from thermosteam.units_of_measure import format_units
from thermosteam.utils import set_figure_size, set_font
from .results import images_folder
from warnings import filterwarnings
from scipy.ndimage.filters import gaussian_filter
from colorpalette import ColorWheel
import os

__all__ = (
    'plot_sorghum_oil_content_and_cane_oil_content_contours_manuscript',
    'plot_metrics_across_composition_configurations_manuscript',
    'plot_sorghum_oil_content_and_cane_oil_content_contours_seminar',
    'plot_metrics_across_composition_configurations',
    'plot_recovery_and_oil_content_contours_manuscript',
    'plot_recovery_and_oil_content_contours',
    'plot_sorghum_oil_content_and_cane_oil_content_contours',
    'plot_recovery_and_oil_content_contours_biodiesel_only',
    'plot_recovery_and_oil_content_contours_with_oilsorghum_only',
    'plot_metrics_across_composition',
    'plot_metrics_across_composition_manuscript',
)

filterwarnings('ignore', category=bst.exceptions.DesignWarning)
line_color_wheel = ColorWheel([
    GG_colors.orange,
    GG_colors.purple,
    GG_colors.green,
    GG_colors.blue,
    GG_colors.yellow,
    colors.CABBI_teal,
    colors.CABBI_grey,
    colors.CABBI_brown,
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
    M, N = shape = colors.shape
    if shape == (2, 2):
        letters=np.array([['A', 'C'], ['B', 'D']])
    elif shape == (3, 2):
        letters=np.array([['A', 'D'], ['B', 'E'], ['C', 'F']])
    elif shape == (2, 3):
        letters=np.array([['A', 'C', 'E'], ['B', 'D', 'F']])
    elif shape == (2, 1):
        letters=np.array([['A'], ['B']])
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
                     fontsize=12, fontweight='bold', zorder=1, bbox=txtbox)

def plot_metrics_across_composition_manuscript(load=True, configuration=None, fs=8, smooth=1):
    set_font(size=fs)
    set_figure_size()
    fig, axes = plot_metrics_across_composition(
        configuration=configuration,
        load=load, 
        smooth=smooth,
    )
    colors = np.zeros([2, 3], object)
    colors[:] = [[light_letter_color, light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 1 - 0.88, 0.63, colors)
    plt.subplots_adjust(right=0.92, hspace=0.1, wspace=0.1, top=0.9, bottom=0.10)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'metrics_across_composition_contours.{i}')
        plt.savefig(file, transparent=True)

def plot_metrics_across_composition_configurations_manuscript(load=True, fs=8, smooth=1):
    set_font(size=fs)
    set_figure_size()
    fig, axes = plot_metrics_across_composition_configurations(
        load=load, 
        smooth=smooth,
    )
    colors = np.zeros([2, 2], object)
    colors[:] = [[light_letter_color, light_letter_color],
                 [light_letter_color, light_letter_color]]
    _add_letter_labels(axes, 1 - 0.86, 0.61, colors)
    plt.subplots_adjust(right=0.92, hspace=0.1, wspace=0.1, top=0.9, bottom=0.10)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'metrics_across_composition_contours.{i}')
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

def metrics_across_oil_and_fiber_content(configuration, load, xlim):
    # Generate contour data
    x0, xf = xlim
    x = np.linspace(x0 / 100, xf / 100, 10)
    y = np.linspace(0.35, 0.75, 5)
    z = np.array([0.6, 0.65, 0.7])
    X, Y, Z = np.meshgrid(x, y, z)
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    file = f'{configuration}_composition_analysis.npy'
    file = os.path.join(folder, file)
    if load:
        data = np.load(file, allow_pickle=True)
    else:
        from warnings import filterwarnings
        bst.System.strict_convergence = False
        filterwarnings('ignore')
        # This returns data of all metrics for the given configuration,
        # but we are mainly interested in MFPP and productivity in L biodiesel per MT cane.
        data = cane.evaluate_metrics_across_composition(
            X, Y, Z, configuration,
        ) 
    np.save(file, data)
    return X, Y, Z, data

def metrics_across_oil_and_fiber_content_configurations(load, xlim):
    # Generate contour data
    x0, xf = xlim
    x = np.linspace(x0 / 100, xf / 100, 10)
    y = np.linspace(0.35, 0.75, 5)
    X, Y = np.meshgrid(x, y)
    folder = os.path.dirname(__file__)
    folder = os.path.join(folder, 'results')
    file = 'configurations_composition_analysis.npy'
    file = os.path.join(folder, file)
    if load:
        data = np.load(file, allow_pickle=True)
    else:
        from warnings import filterwarnings
        bst.System.strict_convergence = False
        filterwarnings('ignore')
        # This returns data of all metrics for the given configuration,
        # but we are mainly interested in MFPP and productivity in L biodiesel per MT cane.
        data = cane.evaluate_metrics_across_composition_configurations(
            X, Y,
        ) 
    np.save(file, data)
    return X, Y, data

def plot_metrics_across_composition_configurations(
        configuration=None, load=False, N_decimals=1, 
        yticks=None, titles=None, 
        cmap=None, smooth=None,
    ):
    xticks = [0,   2,   4,   6,   8,   10]
    xlim = (0, 10)
    # xticks_right = [0,   2,   4,   6,   8,   10]
    # xlim_right = (0, 10)
    # xticks_left = [2, 4, 6, 8, 10]
    # xlim_left = (0.5, 10)
    yticks = [35,  45,  55,  65,  75]
    metrics = COBY, NEP = [cane.competitive_biomass_yield, cane.net_energy_production]
    metric_indices = [cane.all_metric_mockups.index(i) for i in metrics] 
    X, Y, data = metrics_across_oil_and_fiber_content_configurations(load, xlim)
    data = data[:, :, :, metric_indices]
    data[:, 0, 0, 0] = data[:, 1, 0, 0] + (data[:, 1, 0, 0] - data[:, 2, 0, 0])
    data[:, 0, 0, 1] = data[:, 1, 0, 1] + (data[:, 1, 0, 1] - data[:, 2, 0, 1])
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, C, D = data.shape
        for i in range(C):
            for j in range(D):
                try:
                    data[:, :, i, j] = gaussian_filter(np.asarray(data[:, :, i, j], dtype=float), smooth)
                except:
                    breakpoint()
    data = np.swapaxes(data, 2, 3)
    # Plot contours
    xlabel = 'Oil content [dry wt. %]'
    ylabel = "Fiber content [dry wt. %]"
    if titles is None: titles = np.array(['Bioethanol', 'Microbial oil'])
    d0 = data[:, :, 0, :]
    d1 = data[:, :, 1, :]
    metric_bars = [
        # MetricBar(MFPP.name, format_units(MFPP.units), colormaps[1], tickmarks(data[:, :, 0, :], 5, 1, expand=0, p=0.5), 18, 1),
        MetricBar('Comp. biomass yield', 'Dry-'+format_units('MT/ha'), plt.cm.get_cmap('viridis_r'), tickmarks(d0[~np.isnan(d0)], 5, 1, expand=0, p=0.5), 10, 1),
        # MetricBar('Biod. prod.', format_units(BP.units), plt.cm.get_cmap('copper'), tickmarks(data[:, :, 1, :], 8, 1, expand=0, p=0.5), 10, 1),
        MetricBar('Net energy prod.', format_units(NEP.units), plt.cm.get_cmap('inferno'), tickmarks(d1[~np.isnan(d1)], 5, 1, expand=0, p=2), 20, 1),
    ]
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X, 100.*Y, titles, data, xlabel, ylabel, xticks, yticks, metric_bars, 
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
            if str(name) in ('19B', '316'): continue
            oil = data['Stem oil (dw)']['Mean'] * 100
            lines.append(
                (name, 
                 oil, 
                 data['Fiber (dw)']['Mean'] * 100,
                 data['Water (wt)']['Mean'] * 100,
                 color.RGBn)
            )
        for i in axes.flat: 
            if hasattr(i, '_cached_ytwin'):
                plt.sca(i); plt.xlim(xlim)
                plt.sca(i._cached_ytwin); plt.xlim(xlim)
        txtbox = dict(boxstyle='round', facecolor=colors.neutral.shade(20).RGBn, 
                      edgecolor='None', alpha=0.999, pad=0.2)
        
        for *axes_columns, _ in axes:
            for (name, lipid, fiber, moisture, color) in lines:
                left, right = axes_columns
                plt.sca(right._cached_ytwin)
                plt.text(
                    lipid + 0.1, fiber + 1, name, weight='bold', c=color,
                    bbox=txtbox,
                )
                plot_scatter_points(
                    [lipid], [fiber], marker='X', s=50, 
                    color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
                )
                
                plt.sca(left)
                plt.xlim([0.5, 10])
                if lipid > 0.5:
                    plt.text(
                        lipid + 0.1, fiber + 1, name, weight='bold', c=color,
                        bbox=txtbox,
                    )
                    plot_scatter_points(
                        [lipid], [fiber], marker='X', s=50, 
                        color=color, edgecolor=edgecolor, clip_on=False, zorder=1e6,
                    )
                    
    return fig, axes

def plot_metrics_across_composition(
        configuration=None, load=False, N_decimals=1, 
        yticks=None, titles=None, 
        cmap=None, smooth=None,
    ):
    if configuration is None: configuration = 'O2'
    config = cane.parse_configuration(configuration)
    if config.number > 4 or config.number < 0:
        xticks = [0,   2,   4,   6,   8,   10]
        xlim = (0, 10)
    else:
        xticks = [2, 4, 6, 8, 10]
        xlim = (0.5, 10)
    yticks = [35,  45,  55,  65,  75]
    metrics = COBY, NEP = [cane.competitive_biomass_yield, cane.net_energy_production]
    metric_indices = [cane.all_metric_mockups.index(i) for i in metrics] 
    X, Y, Z, data = metrics_across_oil_and_fiber_content(configuration, load, xlim)
    data = data[:, :, :, metric_indices]
    if smooth: # Smooth curves due to heat exchanger network and discontinuities in design decisionss
        A, B, C, D = data.shape
        for i in range(C):
            for j in range(D):
                try:
                    data[:, :, i, j] = gaussian_filter(np.asarray(data[:, :, i, j], dtype=float), smooth)
                except:
                    breakpoint()
    data = np.swapaxes(data, 2, 3)
    # Plot contours
    xlabel = 'Oil content [dry wt. %]'
    ylabel = "Fiber content [dry wt. %]"
    if titles is None: titles = np.array(['60% moisture', '65% moisture', '70% moisture'])
    metric_bars = [
        # MetricBar(MFPP.name, format_units(MFPP.units), colormaps[1], tickmarks(data[:, :, 0, :], 5, 1, expand=0, p=0.5), 18, 1),
        MetricBar('Comp. biomass yield', format_units(COBY.units), colormaps[0], tickmarks(data[:, :, 0, :], 5, 1, expand=0, p=0.5), 10, 1),
        # MetricBar('Biod. prod.', format_units(BP.units), plt.cm.get_cmap('copper'), tickmarks(data[:, :, 1, :], 8, 1, expand=0, p=0.5), 10, 1),
        MetricBar('Net energy prod.', format_units(NEP.units), colormaps[2], tickmarks(data[:, :, 1, :], 5, 1, expand=0, p=2), 20, 1),
    ]
    fig, axes, CSs, CB = plot_contour_2d(
        100.*X[:, :, 0], 100.*Y[:, :, 0], titles, data, xlabel, ylabel, xticks, yticks, metric_bars, 
        styleaxiskw=dict(xtick0=True), label=True, wbar=2.1
    )
    def determine_axis_column(moisture_content):
        for column, axis_moisture in enumerate([60, 65, 70]):
            if abs(moisture_content - axis_moisture) < 2.5:
                return column
        raise RuntimeError('could not determine axis with similar moisture content')
    try:
        df = cane.get_composition_data()
    except:
        pass
    else:
        names = df.index
        lines = []
        for name, color in zip(names, line_color_wheel):
            data = df.loc[name]
            oil = data['Stem oil (dw)']['Mean'] * 100
            if oil < 1.5: continue
            lines.append(
                (name, 
                 oil, 
                 data['Fiber (dw)']['Mean'] * 100,
                 data['Water (wt)']['Mean'] * 100,
                 color.RGBn)
            )
        for i in axes.flat: 
            if hasattr(i, '_cached_ytwin'):
                plt.sca(i); plt.xlim(xlim)
                plt.sca(i._cached_ytwin); plt.xlim(xlim)
        txtbox = dict(boxstyle='round', facecolor=colors.neutral.shade(20).RGBn, 
                      edgecolor='None', alpha=0.99, pad=0.2)
        for *axes_columns, _ in axes:
            for (name, lipid, fiber, moisture, color) in lines:
                index = determine_axis_column(moisture)
                plt.sca(axes_columns[index]._cached_ytwin)
                plt.text(
                    lipid + 0.1, fiber + 1, name, weight='bold', c=color,
                    bbox=txtbox,
                )
                plot_scatter_points(
                    [lipid], [fiber], marker='X', s=50, 
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
    xlabel = "Oil-sorghum oil content [dry wt. %]" 
    ylabel = 'Oilcane oil content\n[dry wt. %]'
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
    ylabel = "Oil content [dry wt. %]"
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
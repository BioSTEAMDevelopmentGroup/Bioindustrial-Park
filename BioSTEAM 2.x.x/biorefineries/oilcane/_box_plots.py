# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:34:00 2021

@author: yrc2
"""
import biosteam as bst
import biorefineries.oilcane as oc
from biosteam.utils import CABBI_colors, colors
from thermosteam.utils import set_figure_size, set_font
from thermosteam.units_of_measure import format_units
from colorpalette import Palette
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from warnings import warn
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from . import _variable_mockups as variables
from ._variable_mockups import (
    tea_monte_carlo_metric_mockups,
    tea_monte_carlo_derivative_metric_mockups,
    lca_monte_carlo_metric_mockups, 
    lca_monte_carlo_derivative_metric_mockups,
    MFPP, ethanol_production, biodiesel_production,
    GWP_economic
)
from ._load_data import (
    images_folder,
    get_monte_carlo,
    spearman_file,
)
import os
from._parse_configuration import format_name
from math import log10, floor

__all__ = (
    'plot_all',
    'plot_breakdown_main_manuscript',
    'plot_montecarlo_main_manuscript',
    'plot_montecarlo_agile_main_manuscript',
    'plot_montecarlo_derivative_main_manuscript',
    'plot_montecarlo_SI',
    'plot_spearman_tea_main_manuscript',
    'plot_spearman_lca_main_manuscript',
    'plot_monte_carlo_across_coordinate',
    'monte_carlo_box_plot',
    'monte_carlo_results',
    'plot_monte_carlo',
    'plot_spearman',
    'plot_configuration_breakdown',
    'plot_TCI_areas_across_oil_content',
)

area_colors = {
    'Feedstock handling': CABBI_colors.teal, 
    'Juicing': CABBI_colors.green_dirty,
    'EtOH prod.': CABBI_colors.blue,
    'Ethanol production': CABBI_colors.blue,
    'Oil ext.': CABBI_colors.brown,
    'Oil extraction': CABBI_colors.brown,
    'Biod. prod.': CABBI_colors.orange,
    'Biodiesel production': CABBI_colors.orange,
    'Pretreatment': CABBI_colors.green,
    'Wastewater treatment': colors.purple,
    'CH&P': CABBI_colors.yellow,
    'Co-Heat and Power': CABBI_colors.yellow,
    'Utilities': colors.red,
    'Storage': CABBI_colors.grey,
    'HXN': colors.orange,
    'Heat exchanger network': colors.orange,
}

area_hatches = {
    'Feedstock handling': 'x', 
    'Juicing': '-',
    'EtOH prod.': '/',
    'Ethanol production': '/',
    'Oil ext.': '\\',
    'Oil extraction': '\\',
    'Biod. prod.': '/|',
    'Biodiesel production': '/|',
    'Pretreatment': '//',
    'Wastewater treatment': r'\\',
    'CH&P': '',
    'Co-Heat and Power': '',
    'Utilities': '\\|',
    'Storage': '',
    'HXN': '+',
    'Heat exchanger network': '+',
}

for i in area_colors: area_colors[i] = area_colors[i].tint(20)
palette = Palette(**area_colors)
letter_color = colors.neutral.shade(25).RGBn

# %% Publication

def plot_all():
    plot_montecarlo_main_manuscript()
    plot_montecarlo_derivative_main_manuscript()
    plot_montecarlo_SI()
    plot_spearman_tea_main_manuscript()
    plot_spearman_lca_main_manuscript()
    plot_breakdown_main_manuscript()

def plot_montecarlo_main_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.75)
    fig, axes = plot_monte_carlo(
        derivative=False, absolute=False, comparison=True,
        tickmarks=None, agile=False, split=True,
        labels=[
            'Conventional',
            'Cellulosic',
            '     Cellu. vs. conven.',
            # 'Conventional',
            # 'Cellulosic',
        ],
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 
            'orange', 'orange', 'orange'
        ])
    )
    for ax, letter in zip(axes, 'ABCDEF'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(2.6, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=12, fontweight='bold')
        if letter in 'CF':
            tub = ylb
            x = 0.5
            plt.text(x, ylb - (yub - ylb) * 0.21, 
                     'Impact of processing\noilcane over sugarcane', 
                     horizontalalignment='center',verticalalignment='center',
                     fontsize=8)
            y = [ylb - (yub - ylb) * 0.4, tub]
            
            line = Line2D([1.5], y, lw=0.8, color='k', alpha=1, linestyle='--')
            line.set_clip_on(False)
            ax.add_line(line)
            
            plt.text(2.17, ylb - (yub - ylb) * 0.26, 
                     'Impact of integrating\ncellulosic ethanol\nproduction', 
                     horizontalalignment='center',verticalalignment='center',
                     fontsize=8)
            
            # line = Line2D([2.5], y, lw=0.8, color='k', alpha=1, linestyle='--')
            # line.set_clip_on(False)
            # ax.add_line(line)
            
            # plt.text(3.5, ylb - (yub - ylb) * 0.35, 
            #          'Impact of integrating oilsorghum\nat an agile oilcane biorefinery', 
            #          horizontalalignment='center',verticalalignment='center',
            #          fontsize=8)
            
            
    plt.subplots_adjust(right=0.96, left=0.105, wspace=0.38, top=0.98, bottom=0.12)
    file = os.path.join(images_folder, 'montecarlo_main_manuscript.svg')
    plt.savefig(file, transparent=True)
    
def plot_montecarlo_agile_main_manuscript():
    set_font(size=8)
    set_figure_size(width=3.3071, aspect_ratio=1.0)
    fig, axes = plot_monte_carlo(
        derivative=False, absolute=False, comparison=True,
        tickmarks=None, agile_only=True, split=False,
        labels=[
            'Conventional',
            'Cellulosic'
        ]
    )
    for ax, letter in zip(axes, 'AB'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(1.6, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=12, fontweight='bold')
        if letter == 'B':
            plt.text(0.5, ylb - (yub - ylb) * 0.25, 
                      'Impact of integrating oilsorghum\nat an agile oilcane biorefinery', 
                      horizontalalignment='center',verticalalignment='center',
                      fontsize=8)
    plt.subplots_adjust(right=0.9, left=0.2, wspace=0.5, top=0.98, bottom=0.15)
    file = os.path.join(images_folder, 'montecarlo_agile_main_manuscript.svg')
    plt.savefig(file, transparent=True)
    
def plot_montecarlo_derivative_main_manuscript():
    set_font(size=8)
    set_figure_size(width=3.3071, aspect_ratio=1.85)
    fig, axes = plot_monte_carlo(
        derivative=True, absolute=True, 
        comparison=False, agile=False,
        split=False,
        # tickmarks=np.array([
        #     [-3, -2, -1, 0, 1, 2, 3, 4, 5],
        #     [-9, -6, -3,  0, 3, 6, 9, 12, 15],
        #     [-2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2],
        #     [-16, -8, 0, 8, 16, 24, 32, 40, 48],
        #     [-400, -300, -200, -100, 0, 100, 200, 300, 400],
        #     [-300, -225, -150, -75, 0, 75, 150, 225, 300]
        # ], dtype=object),
        labels=['Conventional', 'Cellulosic'],
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 'teal', 'brown', 
            'orange', 
        ])
    )
    for ax, letter in zip(axes, 'ABCDEFGH'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(1.6, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=12, fontweight='bold')
    plt.subplots_adjust(left=0.25, right=0.9, wspace=0.4, top=0.98, bottom=0.05)
    file = os.path.join(images_folder, 'montecarlo_derivative_main_manuscript.svg')
    plt.savefig(file, transparent=True)

def plot_montecarlo_SI():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.9)
    fig, axes = plot_monte_carlo(
        absolute=True, comparison=False, split=True,
        expand=0.1, allocated=True,
        # labels=['Sugarcane\nconventional', 'Oilcane\nconventional',
        #         'Sugarcane\ncellulosic', 'Oilcane\ncellulosic',
        #         'Sugarcane\nconventional\nagile', 'Oilcane\nconventional\nagile',
        #         'Sugarcane\ncellulosic\nagile', 'Oilcane\ncellulosic\nagile'],
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJ'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(7.8, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=12, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.40, top=0.98, bottom=0.04)
    file = os.path.join(images_folder, 'montecarlo_SI.svg')
    plt.savefig(file, transparent=True)
    
def plot_spearman_tea_main_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.80)
    plot_spearman(
        configurations=[
            'O1', 'O1*',
            'O2', 'O2*',
        ],
        labels=[
            'Conventional', 'Agile-Conventional',
            'Cellulosic', 'Agile-Cellulosic',
        ],
        cutoff=0.03,
        kind='TEA',
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.08)
    file = os.path.join(images_folder, 'spearman_tea_main_manuscript.svg')
    plt.savefig(file, transparent=True)

def plot_spearman_lca_main_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.65)
    plot_spearman(
        configurations=[
            'O1', 'O1*',
            'O2', 'O2*',
        ],
        labels=[
            'Conventional', 'Agile-Conventional',
            'Cellulosic', 'Agile-Cellulosic',
        ],
        cutoff=0.03,
        kind='LCA',
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.10)
    file = os.path.join(images_folder, 'spearman_lca_main_manuscript.svg')
    plt.savefig(file, transparent=True)

def plot_breakdown_main_manuscript():
    set_font(size=8)
    set_figure_size(aspect_ratio=0.68)
    fig, axes = plt.subplots(nrows=1, ncols=2)
    plt.sca(axes[0])
    plot_configuration_breakdown('O1', ax=axes[0], legend=False)
    plt.sca(axes[1])
    plot_configuration_breakdown('O2', ax=axes[1], legend=True)
    yticks = axes[1].get_yticks()
    plt.yticks(yticks, ['']*len(yticks))
    plt.ylabel('')
    plt.subplots_adjust(left=0.09, right=0.96, wspace=0., top=0.84, bottom=0.31)
    for ax, letter in zip(axes, ['(A) Conventional', '(B) Cellulosic']):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        xlb, xub = plt.xlim()
        plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.2, letter, color=letter_color,
                  horizontalalignment='center',verticalalignment='center',
                  fontsize=12, fontweight='bold')
    
    file = os.path.join(images_folder, 'breakdown_main_manuscript.svg')
    plt.savefig(file, transparent=True)

# %% General

def plot_monte_carlo_across_coordinate(coordinate, data, color_wheel):
    if isinstance(data, list):
        return [plot_monte_carlo_across_coordinate(coordinate, i, color_wheel) for i in data]
    else:
        color = color_wheel.next()
        return bst.plots.plot_montecarlo_across_coordinate(
            coordinate, data,
            light_color=color.tint(50).RGBn,
            dark_color=color.shade(50).RGBn,
        )

def monte_carlo_box_plot(data, positions, light_color, dark_color, width=None):
    if width is None: width = 0.8
    return plt.boxplot(x=data, positions=positions, patch_artist=True,
                     widths=width, whis=[5, 95],
                     boxprops={'facecolor':light_color,
                               'edgecolor':dark_color},
                     medianprops={'color':dark_color,
                                  'linewidth':1.5},
                     flierprops = {'marker':'D',
                                   'markerfacecolor': light_color,
                                   'markeredgecolor': dark_color,
                                   'markersize':3})

def monte_carlo_results(with_units=False):
    results = {}
    ethanol_over_biodiesel = bst.MockVariable('Ethanol over biodiesel', 'Gal/ton', 'Biorefinery')
    for name in oc.configuration_names + oc.comparison_names + oc.other_comparison_names:
        try: 
            df = get_monte_carlo(name)
        except:
            warn(f'could not load {name}', RuntimeWarning)
            continue
        results[name] = dct = {}
        if name in ('O1', 'O2'):
            index = ethanol_over_biodiesel.index
            key = index[1] if with_units else index[1].split(' [')[0]
            data = df[ethanol_production.index].values / df[biodiesel_production.index].values
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            dct[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
        for metric in (*tea_monte_carlo_metric_mockups, *tea_monte_carlo_derivative_metric_mockups,
                       *lca_monte_carlo_metric_mockups, *lca_monte_carlo_derivative_metric_mockups,
                       variables.GWP_ethanol_displacement, variables.GWP_ethanol_allocation):
            index = metric.index
            try:
                data = df[index].values
            except:
                breakpoint()
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            key = index[1] if with_units else index[1].split(' [')[0]
            if key in dct: key = f'{key}, {index[0]}'
            dct[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
    try:
        df_O2O1 = get_monte_carlo('O2 - O1')
        df_O1 = get_monte_carlo('O1')
    except:
        warn('could not load O2 - O1', RuntimeWarning)
    else:
        results['(O2 - O1) / O1'] = relative_results = {}
        for metric in (biodiesel_production, ethanol_production):
            index = metric.index
            key = index[1] if with_units else index[1].split(' [')[0]
            data = (df_O2O1[index].values / df_O1[index].values)
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            relative_results[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
    return results

def plot_monte_carlo(derivative=False, absolute=True, comparison=True,
                     configuration_names=None, comparison_names=None,
                     labels=None, tickmarks=None, agile=True, split=True,
                     expand=None, step_min=None,
                     agile_only=False, xrot=None,
                     color_wheel=None):
    if configuration_names is None: configuration_names = oc.configuration_names
    if comparison_names is None: comparison_names = oc.comparison_names
    factors = []
    GWP_units_gal = '$\\mathrm{kg} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{gal}^{-1}$'
    # GWP_units_economic = '$\\mathrm{g} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{USD}^{-1}$'
    # GWP_units_energy = GWP_units_economic.replace('USD', 'GGE')
    if derivative:
        configuration_names = ['O1', 'O2']
        comparison_names = ['O2 - O1']
        GWP_economic, GWP_ethanol, GWP_biodiesel, GWP_electricity, GWP_crude_glycerol, = lca_monte_carlo_derivative_metric_mockups
        MFPP, TCI, *production, electricity_production, natural_gas_consumption = tea_monte_carlo_derivative_metric_mockups
        rows = [
            MFPP, 
            TCI, 
            production,
            electricity_production,
            natural_gas_consumption,
            GWP_ethanol,
        ]
        factors = [(-1, 1000)] # For changing units of measure; not currently in use
    else:
        GWP_economic, GWP_ethanol, GWP_biodiesel, GWP_electricity, GWP_crude_glycerol, = lca_monte_carlo_metric_mockups
        MFPP, TCI, *production, electricity_production, natural_gas_consumption = tea_monte_carlo_metric_mockups
        GWP_ethanol_displacement = variables.GWP_ethanol_displacement
        if absolute:
            GWP_economic = [variables.GWP_ethanol, variables.GWP_biodiesel]
            GWP_energy = [variables.GWP_ethanol_allocation, variables.GWP_biodiesel_allocation]
            rows = [
                MFPP, 
                TCI, 
                production,
                electricity_production,
                natural_gas_consumption,
                GWP_ethanol_displacement,
                GWP_economic, 
                GWP_energy
            ]
            ylabels = [
                f"MFPP\n[{format_units('USD/ton')}]",
                f"TCI\n[{format_units('10^6*USD')}]",
                f"Production\n[{format_units('Gal/ton')}]",
                f"Elec. prod.\n[{format_units('kWhr/ton')}]",
                f"NG cons.\n[{format_units('cf/ton')}]",
                "GWP$_{\\mathrm{displacement}}$" f"\n[{GWP_units_gal}]",
                "GWP$_{\\mathrm{energy}}$" f"\n[{GWP_units_gal}]",
                "GWP$_{\\mathrm{economic}}$" f"\n[{GWP_units_gal}]",
            ]
        else:
            GWP_ethanol_allocation = variables.GWP_ethanol_allocation
            if agile_only:
                rows = [
                    MFPP, 
                    TCI, 
                ]
                ylabels = [
                    f"MFPP\n[{format_units('USD/ton')}]",
                    f"TCI\n[{format_units('10^6*USD')}]",
                ]
            else:
                rows = [
                    MFPP, 
                    TCI, 
                    production,
                    GWP_ethanol_displacement,
                    GWP_ethanol_allocation,
                    GWP_ethanol, # economic
                ]
                ylabels = [
                    f"MFPP\n[{format_units('USD/ton')}]",
                    f"TCI\n[{format_units('10^6*USD')}]",
                    f"Production\n[{format_units('Gal/ton')}]",
                    "GWP$_{\\mathrm{displacement}}$" f"\n[{GWP_units_gal}]",
                    "GWP$_{\\mathrm{energy}}$" f"\n[{GWP_units_gal}]",
                    "GWP$_{\\mathrm{economic}}$" f"\n[{GWP_units_gal}]",
                ]
    if derivative:
        ylabels = [
            r"$\Delta$" + format_units(r"MFPP/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('USD/ton')}]",
            r"$\Delta$" + format_units(r"TCI/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('10^6*USD')}]",
            r"$\Delta$" + format_units(r"Prod./OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('Gal/ton')}]",
            r"$\Delta$" + format_units(r"EP/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('kWhr/ton')}]",
            r"$\Delta$" + format_units(r"NGC/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('cf/ton')}]",
            r"$\Delta$" + r"GWP$_{\mathrm{econ.}} \cdot \Delta \mathrm{OC}^{-1}$" f"\n[{GWP_units_gal.replace('kg','g')}]",
        ]
    elif comparison and not absolute:
        ylabels = [r"$\Delta$" + i for i in ylabels]
    if color_wheel is None:
        color_wheel = CABBI_colors.wheel()
    N_rows = len(rows)
    if split:
        ncols = 2
        nrows = int(round(N_rows / 2))
    else:
        ncols = 1
        nrows = N_rows
    fig, axes_box = plt.subplots(ncols=ncols, nrows=nrows)
    plt.subplots_adjust(wspace=0.45)
    axes = axes_box.transpose()
    axes = axes.flatten()
    
    combined = absolute and comparison
    if agile_only:
        configuration_names = [i for i in configuration_names if '*' in i]
        comparison_names = [i for i in comparison_names if '*' in i]
    elif not agile:
        configuration_names = [i for i in configuration_names if '*' not in i]
        comparison_names = [i for i in comparison_names if '*' not in i]
    if combined:
        columns = configurations = configuration_names + comparison_names
    elif absolute:
        columns = configurations = configuration_names
    elif comparison:
        columns = configurations = comparison_names
    else:
        columns = configurations = []
    N_cols = len(columns)    
    xtext = labels or [format_name(i).replace(' ', '') for i in configurations]
    N_marks = len(xtext)
    xticks = tuple(range(N_marks))
    
    def get_data(metric, name):
        try:
            df = get_monte_carlo(name, metric)
        except:
            return np.zeros([1, 1])
        else:
            values = df.values
            return values
    
    def plot(arr, position):
        if arr.ndim == 2:
            N = arr.shape[1]
            width = 0.618 / N
            boxwidth = 0.618 / (N + 1/N)
            plots = []
            for i in range(N):
                color = color_wheel.next()
                boxplot = monte_carlo_box_plot(
                    data=arr[:, i], positions=[position + (i-1/N)*width], 
                    light_color=color.RGBn, 
                    dark_color=color.shade(60).RGBn,
                    width=boxwidth,
                )
                plots.append(boxplot)
            return plots
        else:
            color = color_wheel.next()
            return monte_carlo_box_plot(
                data=arr, positions=[position], 
                light_color=color.RGBn, 
                dark_color=color.shade(60).RGBn,
                width=0.618,
            )
    
    data = np.zeros([N_rows, N_cols], dtype=object)
    data[:] = [[get_data(i, j) for j in columns] for i in rows]
    for i, j in factors: data[i, :] *= j
    def round2sigfigs(x):
        try:
            value = round(x, 1-int(floor(log10(abs(x)))))
        except:
            return int(x)
        if int(value) == value:
            return int(value)
        else:
            return value
        
    if tickmarks is None: 
        tickmarks = [
            bst.plots.rounded_tickmarks_from_data(
                i, step_min=step_min, N_ticks=8, lb_max=0, center=0,
                f=round2sigfigs, expand=expand,
            ) 
            for i in data
        ]

    x0 = len(configuration_names) - 0.5
    xf = len(columns) - 0.5
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        if combined:
            bst.plots.plot_vertical_line(x0)
            ax.axvspan(x0, xf, color=colors.purple_tint.tint(60).RGBn)
        plt.xlim(-0.5, xf)

    for j in range(N_cols):
        color_wheel.restart()
        for i in range(N_rows):
            ax = axes[i]
            plt.sca(ax)
            plot(data[i, j], j)
            plt.ylabel(ylabels[i])
    
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        yticks = tickmarks[i]
        plt.ylim([yticks[0], yticks[1]])
        if yticks[0] < 0.:
            bst.plots.plot_horizontal_line(0, color=CABBI_colors.black.RGBn, lw=0.8, linestyle='--')
        try:
            xticklabels = xtext if ax in axes_box[-1] else []
        except:
            xticklabels = xtext if i == N_rows - 1 else []
        bst.plots.style_axis(ax,  
            xticks = xticks,
            yticks = yticks,
            xticklabels= xticklabels, 
            ytick0=False,
            ytickf=False,
            offset_xticks=True,
            xrot=xrot,
            ha='center',
        )
    fig.align_ylabels(axes)
    plt.subplots_adjust(hspace=0)
    plt.sca(axes[1])
    # legend = plt.legend(
    #     handles=[
    #         mpatches.Patch(facecolor=color_wheel[0].RGBn, 
    #                        edgecolor=CABBI_colors.black.RGBn,
    #                        label='Oilcane only'),
    #         mpatches.Patch(facecolor=color_wheel[1].RGBn, 
    #                        edgecolor=CABBI_colors.black.RGBn,
    #                        label='Oilcane & oilsorghum'),
    #     ], 
    #     bbox_to_anchor=(0, 1, 1, 0), 
    #     loc="lower right", 
    #     # mode="expand", 
    #     # ncol=2
    # )
    # legend.get_frame().set_linewidth(0.0)
    return fig, axes

def plot_spearman(configurations, top=None, labels=None, metric=None, cutoff=None,
                  kind='TEA'):
    if metric is None:
        if kind == 'TEA':
            metric = MFPP
            metric_name = metric.name
        elif kind == 'LCA':
            metric = GWP_economic
            metric_name = r'GWP$_{\mathrm{economic}}$'
        else:
            raise ValueError(f"invalid kind '{kind}'")
    else:
        if metric == 'MFPP':
            metric = MFPP
        elif metric == 'GWP':
            metric = GWP_economic
        metric_name = metric.name
    stream_price = format_units('USD/gal')
    USD_ton = format_units('USD/ton')
    ng_price = format_units('USD/cf')
    electricity_price = format_units('USD/kWhr')
    operating_days = format_units('day/yr')
    capacity = format_units('ton/hr')
    titer = format_units('g/L')
    productivity = format_units('g/L/hr')
    material_GWP = '$\\mathrm{kg} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{kg}^{-1}$'
    feedstock_GWP = '$\\mathrm{g} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{kg}^{-1}$'
    index, ignored_list = zip(*[
         ('Bagasse oil retention [40 $-$ 70 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Oil extraction efficiency [baseline + 0 $-$ 20 %]', ['S2', 'S1', 'S2*', 'S1*']),
        (f'Plant capacity [330 $-$ 404 {capacity}]', []),
        (f'Ethanol price [1.02, 1.80, 2.87 {stream_price}]', []),
        (f'Relative biodiesel price [0.31, 2.98, 4.11 {stream_price}]', []),
        (f'Natural gas price [3.71, 4.73, 6.18 {ng_price}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Electricity price [0.0583, 0.065, 0.069 {electricity_price}]', ['S2', 'O2', 'S2*', 'O2*']),
        (f'Operating days [180 $-$ 210 {operating_days}]', []),
         ('IRR [10 $-$ 15 %]', []),
        (f'Crude glycerol price [91 $-$ 200 {USD_ton}]', ['S2', 'S1', 'S2*', 'S1*']),
        (f'Pure glycerol price [501 $-$ 678 {USD_ton}]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Saccharification reaction time [54 $-$ 90 hr]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Cellulase price [144 $-$ 240 {USD_ton}]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Cellulase loading [1.5 $-$ 2.5 wt. % cellulose]', ['S1', 'O1', 'S1*', 'O1*']),
         ('PTRS base cost [14.9 $-$ 24.7 MMUSD]', ['S1', 'O1', 'S1*', 'O1*']),
         # ('Pretreatment reactor system base cost [14.9 $-$ 24.7 MMUSD]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Cane glucose yield [85 $-$ 97.5 %]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Sorghum glucose yield [85 $-$ 97.5 %]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Cane xylose yield [65 $-$ 97.5 %]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Sorghum xylose yield [65 $-$ 97.5 %]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Glucose to ethanol yield [90 $-$ 95 %]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Xylose to ethanol yield [50 $-$ 95 %]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Titer [65 $-$ 130 {titer}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Productivity [1.0 $-$ 2.0 {productivity}]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Cane PL content [7.5 $-$ 12.5 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Sorghum PL content [7.5 $-$ 12.5 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Cane FFA content [7.5 $-$ 12.5 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Sorghum FFA content [7.5 $-$ 12.5 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Cane oil content [5 $-$ 15 dry wt. %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Relative sorghum oil content [-3 $-$ 0 dry wt. %]', ['S2', 'S1', 'S2*', 'S1*', 'O2', 'O1']),
         ('TAG to FFA conversion [17.25 $-$ 28.75 % theoretical]', ['S1', 'O1', 'S1*', 'O1*']),
        # TODO: change lower upper values to baseline +- 10%
        (f'Feedstock GWPCF [26.3 $-$ 44.0 {feedstock_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Methanol GWPCF [0.338 $-$ 0.563 {material_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Pure glycerine GWPCF [1.25 $-$ 2.08 {material_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Cellulase GWPCF [6.05 $-$ 10.1 {material_GWP}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Natural gas GWPCF [0.297 $-$ 0.363 {material_GWP}]', ['S1', 'O1', 'S1*', 'O1*']),
    ])
    ignored_dct = {
        'S1': [],
        'O1': [],
        'S2': [],
        'O2': [],
        'S1*': [],
        'O1*': [],
        'S2*': [],
        'O2*': [],
    }
    for i, ignored in enumerate(ignored_list):
        for name in ignored: ignored_dct[name].append(i)
        index_name = index[i]
        if kind == 'LCA':
            for term in ('cost', 'price', 'IRR', 'time', 'capacity'):
                if term in index_name:
                    for name in ignored_dct: ignored_dct[name].append(i)
                    break
        elif kind == 'TEA':
            if 'GWP' in index_name:
                for name in ignored_dct: ignored_dct[name].append(i)
        else:
            raise ValueError(f"invalid kind '{kind}'")
    
    rhos = []
    for name in configurations:
        file = spearman_file(name)
        try: 
            df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
        except: 
            warning = RuntimeWarning(f"file '{file}' not found")
            warn(warning)
            continue
        s = df[metric.index]
        s.iloc[ignored_dct[name]] = 0.
        rhos.append(s)
    color_wheel = [CABBI_colors.orange, CABBI_colors.green_soft, CABBI_colors.blue, CABBI_colors.brown]
    fig, ax = bst.plots.plot_spearman_2d(rhos, top=top, index=index, cutoff=cutoff, 
                                         color_wheel=color_wheel,
                                         name=metric_name)
    plt.legend(
        handles=[
            mpatches.Patch(
                color=color_wheel[i].RGBn, 
                label=labels[i] if labels else format_name(configurations[i])
            )
            for i in range(len(configurations))
        ], 
        loc='lower left'
    )
    return fig, ax

def plot_configuration_breakdown(name, across_coordinate=False, **kwargs):
    oc.load(name)
    if across_coordinate:
        return bst.plots.plot_unit_groups_across_coordinate(
            oc.set_cane_oil_content,
            [5, 7.5, 10, 12.5],
            'Feedstock oil content [dry wt. %]',
            oc.unit_groups,
            colors=[area_colors[i.name].RGBn for i in oc.unit_groups],
            hatches=[area_hatches[i.name] for i in oc.unit_groups],
            **kwargs,
        )
    else:
        def format_total(x):
            if x < 1e3:
                return format(x, '.3g')
            else:
                x = int(x)
                n = 10 ** (len(str(x)) - 3)
                value = int(round(x / n) * n)
                return format(value, ',')
        for i in oc.unit_groups: 
            if i.name == 'EtOH prod.':
                i.name = 'Ethanol production'
            elif i.name == 'Oil ext.':
                i.name = 'Oil extraction'
            elif i.name == 'Biod. prod.':
                i.name = 'Biodiesel production'
            i.metrics[0].name = 'Inst. eq.\ncost'
            i.metrics[3].name = 'Elec.\ncons.'
            i.metrics[4].name = 'Mat.\ncost'
        
        return bst.plots.plot_unit_groups(
            oc.unit_groups,
            colors=[area_colors[i.name].RGBn for i in oc.unit_groups],
            hatches=[area_hatches[i.name] for i in oc.unit_groups],
            format_total=format_total,
            fraction=True,
            legend_kwargs=dict(
                loc='lower center',
                ncol=4,
                bbox_to_anchor=(0, -0.52),
                labelspacing=1.5, handlelength=2.8,
                handleheight=1, scale=0.8,
            ),
            **kwargs,
        )

def plot_TCI_areas_across_oil_content(configuration='O2'):
    oc.load(configuration)
    data = {i.name: [] for i in oc.unit_groups}
    increasing_areas = []
    decreasing_areas = []
    oil_contents = np.linspace(5, 15, 10)
    for i in oil_contents:
        oc.set_cane_oil_content(i)
        oc.sys.simulate()
        for i in oc.unit_groups: data[i.name].append(i.get_installed_cost())
    for name, group_data in data.items():
        lb, *_, ub = group_data
        if ub > lb: 
            increasing_areas.append(group_data)
        else:
            decreasing_areas.append(group_data)
    increasing_values = np.sum(increasing_areas, axis=0)
    increasing_values -= increasing_values[0]
    decreasing_values = np.sum(decreasing_areas, axis=0)
    decreasing_values -= decreasing_values[-1]
    plt.plot(oil_contents, increasing_values, label='Oil & fiber areas')
    plt.plot(oil_contents, decreasing_values, label='Sugar areas')
    
# def plot_monte_carlo_across_oil_content(kind=0, derivative=False):
#     MFPP, TCI, *production, electricity_production, natural_gas_consumption = tea_monte_carlo_metric_mockups
#     rows = [MFPP, TCI, production]
#     if kind == 0:
#         columns = across_oil_content_names
#     elif kind == 1:
#         columns = across_oil_content_agile_names
#     elif kind == 2:
#         columns = across_oil_content_comparison_names
#     elif kind == 3:
#         columns = across_oil_content_agile_comparison_names
#     elif kind == 4:
#         columns = across_oil_content_agile_direct_comparison_names
#     else:
#         raise NotImplementedError(str(kind))
#     if derivative:
#         x = 100 * (oil_content[:-1] + np.diff(oil_content) / 2.)
#         ylabels = [
#             f"MFPP der. [{format_units('USD/ton')}]",
#             f"TCI der. [{format_units('10^6*USD')}]",
#             f"Production der. [{format_units('gal/ton')}]"
#         ]
#     else:
#         x = 100 * oil_content
#         ylabels = [
#             f"MFPP$\backprime$ [{format_units('USD/ton')}]",
#             f"TCI [{format_units('10^6*USD')}]",
#             f"Production [{format_units('gal/ton')}]"
#         ]
#     N_cols = len(columns)
#     N_rows = len(rows)
#     fig, axes = plt.subplots(ncols=N_cols, nrows=N_rows)
#     data = np.zeros([N_rows, N_cols], dtype=object)
    
#     def get_data(metric, name):
#         if isinstance(metric, bst.Variable):
#             return get_monte_carlo_across_oil_content(name, metric, derivative)
#         else:
#             return [get_data(i, name) for i in metric]
    
#     data = np.array([[get_data(i, j) for j in columns] for i in rows])
#     tickmarks = [None] * N_rows
#     get_max = lambda x: max([i.max() for i in x]) if isinstance(x, list) else x.max()
#     get_min = lambda x: min([i.min() for i in x]) if isinstance(x, list) else x.min()
#     N_ticks = 5
#     for r in range(N_rows):
#         lb = min(min([get_min(i) for i in data[r, :]]), 0)
#         ub = max([get_max(i) for i in data[r, :]])
#         diff = 0.1 * (ub - lb)
#         ub += diff
#         if derivative:
#             lb = floor(lb)
#             ub = ceil(ub)
#             step = (ub - lb) / (N_ticks - 1)
#             tickmarks[r] = [0, 1] if step == 0 else [int(lb + step * i) for i in range(N_ticks)]
#         else:
#             if rows[r] is MFPP:
#                 if kind == 0 or kind == 1:
#                     tickmarks[r] = [-20, 0, 20, 40, 60]
#                 elif kind == 2:
#                     tickmarks[r] = [-20, -10, 0, 10, 20]
#                 elif kind == 3:
#                     tickmarks[r] = [-10, 0, 10, 20, 30]
#                 elif kind == 4:
#                     tickmarks[r] = [-5, 0, 5, 10, 15]
#                 continue
#             lb = floor(lb / 15) * 15
#             ub = ceil(ub / 15) * 15
#             step = (ub - lb) / (N_ticks - 1)
#             tickmarks[r] = [0, 1] if step == 0 else [int(lb + step * i) for i in range(N_ticks)]
#     color_wheel = CABBI_colors.wheel()
#     for j in range(N_cols):
#         color_wheel.restart()
#         for i in range(N_rows):
#             arr = data[i, j]
#             ax = axes[i, j]
#             plt.sca(ax)
#             percentiles = plot_monte_carlo_across_coordinate(x, arr, color_wheel)
#             if i == 0: ax.set_title(format_name(columns[j]))
#             xticklabels = i == N_rows - 1
#             yticklabels = j == 0
#             if xticklabels: plt.xlabel('Oil content [dry wt. %]')
#             if yticklabels: plt.ylabel(ylabels[i])
#             bst.plots.style_axis(ax,  
#                                  xticks = [5, 10, 15],
#                                  yticks = tickmarks[i],
#                                  xticklabels= xticklabels, 
#                                  yticklabels= yticklabels,
#                                  ytick0=False)
#     for i in range(N_cols): fig.align_ylabels(axes[:, i])
#     plt.subplots_adjust(hspace=0.1, wspace=0.1)
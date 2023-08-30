# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 01:34:00 2021

@author: yrc2
"""
from typing import NamedTuple
import biosteam as bst
from biosteam.utils import CABBI_colors, GG_colors, colors
from thermosteam.utils import set_figure_size, set_font, roundsigfigs
from thermosteam.units_of_measure import format_units
from colorpalette import Palette
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from warnings import warn
import numpy as np
import pandas as pd
from biorefineries import cane
import matplotlib.patches as patch
import flexsolve as flx
from . import feature_mockups as features
from .feature_mockups import (
    get_YRCP2023_spearman_names,
    tea_monte_carlo_metric_mockups,
    tea_monte_carlo_derivative_metric_mockups,
    lca_monte_carlo_metric_mockups, 
    lca_monte_carlo_derivative_metric_mockups,
    ROI, MFPP, MBSP, TCI, electricity_production, GWP_biofuel_allocation,
    # natural_gas_consumption,
    ethanol_production, biodiesel_production, biodiesel_yield, IRR,
    GWP_ethanol, GWP_biodiesel, GWP_electricity,
    GWP_ethanol_allocation, GWP_biodiesel_allocation, GWP_biodiesel_displacement,
    GWP_economic, MFPP_derivative, 
    TCI_derivative, 
    ethanol_production_derivative,
    biodiesel_production_derivative,
    electricity_production_derivative,
    # natural_gas_consumption_derivative,
    GWP_ethanol_derivative,
)
from .results import (
    images_folder,
    monte_carlo_file,
    get_monte_carlo,
    get_line_monte_carlo,
    spearman_file,
)
import os
from colorpalette import ColorWheel, Color
from.parse_configuration import format_name
from scipy import interpolate

__all__ = (
    'plot_all',
    # 'plot_montecarlo_main_manuscript',
    'plot_breakdowns',
    'plot_montecarlo_feedstock_comparison',
    'plot_montecarlo_configuration_comparison',
    'plot_montecarlo_agile_comparison',
    'plot_montecarlo_derivative',
    'plot_montecarlo_absolute',
    'plot_lines_monte_carlo_manuscript',
    'plot_spearman_tea',
    'plot_spearman_lca',
    'plot_spearman_tea_short',
    'plot_spearman_lca_short',
    'plot_monte_carlo_across_coordinate',
    'monte_carlo_box_plot',
    'plot_monte_carlo',
    'plot_spearman',
    'plot_configuration_breakdown',
    'plot_feedstock_conventional_comparison_kde',
    'plot_feedstock_cellulosic_comparison_kde',
    'plot_configuration_comparison_kde',
    'plot_open_comparison_kde',
    'plot_feedstock_comparison_kde',
    'plot_feedstock_comparison_kde_2023',
    'plot_crude_configuration_comparison_kde',
    'plot_agile_comparison_kde',
    'plot_separated_configuration_comparison_kde',
    'plot_unlabeled_feedstock_conventional_comparison_kde',
    'plot_competitive_biomass_yield_across_oil_content',
    'plot_competitive_microbial_oil_yield_across_oil_content',
    'plot_microbial_oil_feedstock_comparison_kde',
    'plot_microbial_oil_bioethanol_comparison_kde',
    'plot_microbial_oil_kde_2023',
    'plot_spearman_tea_YRCP2023',
    'area_colors',
    'area_colors_biodiesel',
    'area_hatches',
    'area_hatches_biodiesel',
    'plot_spearman_YRCP2023',
    'plot_spearman_lca_YRCP2023',
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
    'Oil prod. & ext.': CABBI_colors.blue,
    'Oil production and extraction': CABBI_colors.blue,
}

area_colors_biodiesel = {
    'Juicing': CABBI_colors.green_dirty,
    'Oil production': CABBI_colors.blue,
    'Pretreatment': CABBI_colors.teal,
    'Wastewater treatment': colors.purple,
    'CH&P': CABBI_colors.yellow,
    'Co-Heat and Power': CABBI_colors.yellow,
    'Utilities': colors.red,
    'Storage': CABBI_colors.grey,
    'HXN': colors.orange,
    'Heat exchanger network': colors.orange,
    'Biod. prod.': CABBI_colors.orange,
    'Biodiesel production': CABBI_colors.orange,
}
area_hatches_biodiesel = {
    'Juicing': '',
    'Oil production': '',
    'Biod. prod.': '',
    'Biodiesel production': '',
    'Pretreatment': '',
    'Wastewater treatment': r'',
    'CH&P': '',
    'Co-Heat and Power': '',
    'Utilities': '',
    'Storage': '',
    'HXN': '',
    'Heat exchanger network': '',
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
    'TAG prod.': '/',
    'Oil prod. & ext.': '/',
    'Oil production and extraction': '/',
}
line_color_wheel = ColorWheel([
    GG_colors.orange,
    GG_colors.blue,
    GG_colors.purple,
    GG_colors.green,
    GG_colors.yellow,
    colors.CABBI_teal,
    colors.CABBI_grey,
    colors.CABBI_brown,
])
for i, color in enumerate(line_color_wheel.colors): line_color_wheel.colors[i] = color.tint(20)
edgecolor = (*colors.CABBI_black.RGBn, 1)
for i in area_colors: area_colors[i] = area_colors[i].tint(20)
for i in area_colors_biodiesel: area_colors_biodiesel[i] = area_colors_biodiesel[i].tint(20)
palette = Palette(**area_colors)
letter_color = colors.neutral.shade(25).RGBn
GWP_units_L = '$\\mathrm{kg} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{L}^{-1}$'
GWP_units_GGE = GWP_units_L.replace('L', 'GGE')
GWP_units_L_small = GWP_units_L.replace('kg', 'g')
CABBI_colors.orange_hatch = CABBI_colors.orange.copy(hatch='////')
    
ethanol_over_biodiesel = bst.MockVariable('Ethanol over biodiesel', 'L/MT', 'Biorefinery')
GWP_ethanol_displacement = features.GWP_ethanol_displacement
production = (ethanol_production, biodiesel_production)

mc_metric_settings = {
    'MFPP': (MFPP, f"MFPP\n[{format_units('USD/MT')}]", None),
    'TCI': (TCI, f"TCI\n[{format_units('10^6*USD')}]", None),
    'production': (production, f"Production\n[{format_units('L/MT')}]", None),
    'electricity_production': (electricity_production, f"Elec. prod.\n[{format_units('kWhr/MT')}]", None),
    # 'natural_gas_consumption': (natural_gas_consumption, f"NG cons.\n[{format_units('m^3/MT')}]", None),
    'GWP_ethanol_displacement': (GWP_ethanol_displacement, "GWP$_{\\mathrm{displacement}}$" f"\n[{GWP_units_L}]", None),
    'GWP_economic': ((GWP_ethanol, GWP_biodiesel), "GWP$_{\\mathrm{economic}}$" f"\n[{GWP_units_L}]", None),
    'GWP_energy': ((GWP_ethanol_allocation, GWP_biodiesel_allocation), "GWP$_{\\mathrm{energy}}$" f"\n[{GWP_units_L}]", None),
}

mc_line_metric_settings = {
    # 'Biodiesel production': (biodiesel_production, f"Biodiesel production\n[{format_units('L/MT')}]"),
    'MBSP': (MBSP, f"MBSP\n[{format_units('USD/gal')}]"),
    'TCI': (TCI, f"TCI\n[{format_units('10^6*USD')}]"),
    'IRR': (IRR, "IRR\n[%]"),
    'ROI': (ROI, f"ROI\n[{format_units('%/y')}]"),
    # 'Electricity production': (electricity_production, f"Elec. prod.\n[{format_units('kWhr/MT')}]"),
    'Biodiesel yield': (biodiesel_yield, f"Biodiesel yield\n[{format_units('L/ha/y')}]"),
    'GWP energy': (GWP_biodiesel_allocation,"GWP$_{\\mathrm{energy}}$\n" f"[{GWP_units_L}]"),
    'GWP economic': (GWP_biodiesel, "GWP$_{\\mathrm{economic}}$\n" f"[{GWP_units_L}]"),
    'GWP displacement': (GWP_biodiesel_displacement, "GWP$_{\\mathrm{displacement}}$\n" f"[{GWP_units_L}]"),
}
line_metrics = dict(
    economic=['MBSP', 'TCI', 'Biodiesel yield'],
    environmental=['GWP energy', 'GWP economic', 'GWP displacement']
)
mc_comparison_settings = {
    'MFPP': (MFPP, r"$\Delta$" + f"MFPP\n[{format_units('USD/MT')}]", None),
    'MBSP': (MBSP, r"$\Delta$" + f"MBSP\n[{format_units('USD/L')}]", None),
    'TCI': (TCI, r"$\Delta$" + f"TCI\n[{format_units('10^6*USD')}]", None),
    'production': (production, r"$\Delta$" + f"Production\n[{format_units('L/MT')}]", None),
    'electricity_production': (electricity_production, r"$\Delta$" + f"Elec. prod.\n[{format_units('kWhr/MT')}]", None),
    'GWP biofuel': (GWP_biofuel_allocation, r"$\Delta$" +  f"GWP [{GWP_units_GGE}]", None),
    # 'natural_gas_consumption': (natural_gas_consumption, r"$\Delta$" + f"NG cons.\n[{format_units('m^3/MT')}]", None),
    'GWP_ethanol_displacement': (GWP_ethanol_displacement, r"$\Delta$" + "GWP$_{\\mathrm{displacement}}$" f"\n[{GWP_units_L}]", None),
    'GWP_economic': (GWP_ethanol, r"$\Delta$" + "GWP$_{\\mathrm{economic}}$" f"\n[{GWP_units_L}]", None),
    'GWP_energy': (GWP_ethanol_allocation, r"$\Delta$" + "GWP$_{\\mathrm{energy}}$" f"\n[{GWP_units_L}]", None),
    'GWP_property_allocation': ((GWP_ethanol, GWP_ethanol_allocation), r"$\Delta$" + f"GWP\n[{GWP_units_L}]", None),
}

mc_derivative_metric_settings = {
    'MFPP': (MFPP_derivative, r"$\Delta$" + format_units(r"MFPP/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('USD/MT')}]", None),
    'TCI': (TCI_derivative,  r"$\Delta$" + format_units(r"TCI/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('10^6*USD')}]", None),
    'production': ((ethanol_production_derivative, biodiesel_production_derivative), r"$\Delta$" + format_units(r"Prod./OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('L/MT')}]", None),
    'electricity_production': (electricity_production_derivative, r"$\Delta$" + format_units(r"EP/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('kWhr/MT')}]", None),
    # 'natural_gas_consumption': (natural_gas_consumption_derivative, r"$\Delta$" + format_units(r"NGC/OC").replace('cdot', r'cdot \Delta') + f"\n[{format_units('m^3/MT')}]", None),
    'GWP_economic': (GWP_ethanol_derivative, r"$\Delta$" + r"GWP $\cdot \Delta \mathrm{OC}^{-1}$" f"\n[{GWP_units_L_small}]", 1000),
}

kde_metric_settings = {j[0]: j for j in mc_metric_settings.values()}
kde_comparison_settings = {j[0]: j for j in mc_comparison_settings.values()}
kde_derivative_settings = {j[0]: j for j in mc_derivative_metric_settings.values()}

# %% Helpful classes

class Point(NamedTuple):
    x: float
    y: float
    c: np.array

# %% Plots for publication

def plot_all():
    # plot_montecarlo_main_manuscript()
    plot_montecarlo_absolute()
    plot_spearman_tea()
    plot_spearman_lca()
    plot_breakdowns()

def plot_montecarlo_feedstock_comparison(axes_box=None, letters=None, 
                                         single_column=True):
    if single_column:
        width = 'half'
        aspect_ratio = 2.25
        ncols = 1
        left = 0.255
        bottom = 0.05
    else:
        width = None
        aspect_ratio = 0.75
        left = 0.105
        bottom = 0.12
        ncols = 3
    if axes_box is None:
        set_font(size=8)
        set_figure_size(width=width, aspect_ratio=aspect_ratio)       
    fig, axes = plot_monte_carlo(
        derivative=False, absolute=False, comparison=True,
        tickmarks=None, agile=False, ncols=ncols, axes_box=axes_box,
        labels=[
            'Direct Cogeneration',
            'Integrated Co-Fermentation',
            # 'Direct Cogeneration',
            # 'Integrated Co-Fermentation',
        ],
        comparison_names=['O1 - S1', 'O2 - S2'],
        metrics = ['MFPP', 'TCI', 'production', 'GWP_property_allocation', 
                   'natural_gas_consumption', 'electricity_production'],
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 
            'orange', 'orange_hatch', 'grey', 'brown',
        ])
    )
    for ax, letter in zip(axes, 'ABCDEFGH' if letters is None else letters):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(1.65, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
        # if axes_box is None and letter in 'DH':
        #     x = 0.5
        #     plt.text(x, ylb - (yub - ylb) * 0.3, 
        #              'Impact of processing\noilcane over sugarcane', 
        #              horizontalalignment='center',verticalalignment='center',
        #              fontsize=8)
    if axes_box is None:
        plt.subplots_adjust(right=0.96, left=left, wspace=0.38, top=0.98, bottom=bottom)
        for i in ('svg', 'png'):
            file = os.path.join(images_folder, f'montecarlo_feedstock_comparison.{i}')
            plt.savefig(file, transparent=True)

    
def plot_montecarlo_configuration_comparison(axes_box=None, letters=None,
                                             single_column=True):
    if single_column:
        width = 'half'
        aspect_ratio = 2.25
        ncols = 1
        left = 0.255
        bottom = 0.05
        x = 1.65
        metrics= ['MFPP', 'TCI', 'production', 'GWP_property_allocation',
                  'natural_gas_consumption', 'electricity_production']
    else:
        width = None
        aspect_ratio = 0.75
        left = 0.105
        bottom = 0.12
        ncols = 2
        x = 0.58
        metrics= ['MFPP', 'TCI', 'production', 'GWP_property_allocation']
    if axes_box is None:
        set_font(size=8)
        set_figure_size(width=width, aspect_ratio=aspect_ratio)
    fig, axes = plot_monte_carlo(
        derivative=False, absolute=False, comparison=True,
        tickmarks=None, agile=False, ncols=ncols, axes_box=axes_box,
        labels=[
            'Oilcane',
            # 'Sugarcane',
        ],
        comparison_names=[
            'O2 - O1', 
            # 'S2 - S1'
        ],
        metrics=metrics,
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 
            'orange', 'orange_hatch', 
        ])
    )
    for ax, letter in zip(axes, 'ABCDEF' if letters is None else letters):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(x, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    if axes_box is None:
        plt.subplots_adjust(right=0.96, left=left, wspace=0.38, top=0.98, bottom=bottom)
        for i in ('svg', 'png'):
            file = os.path.join(images_folder, f'montecarlo_configuration_comparison.{i}')
            plt.savefig(file, transparent=True)
    
def plot_montecarlo_agile_comparison(axes_box=None, letters=None):
    if axes_box is None:
        set_font(size=8)
        set_figure_size(width=3.3071, aspect_ratio=1.0)
    fig, axes = plot_monte_carlo(
        derivative=False, absolute=False, comparison=True,
        tickmarks=None, agile_only=True, ncols=1,
        labels=[
            'Direct Cogeneration',
            'Integrated Co-Fermentation'
        ],
        metrics=['MFPP', 'TCI'],
        axes_box=axes_box,
    )
    for ax, letter in zip(axes, 'AB'  if letters is None else letters):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(1.65, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
        if axes_box is None and letter == 'B':
            plt.text(0.5, ylb - (yub - ylb) * 0.25, 
                      'Impact of integrating oilsorghum\nat an agile oilcane biorefinery', 
                      horizontalalignment='center',verticalalignment='center',
                      fontsize=8)
    if axes_box is None:
        plt.subplots_adjust(right=0.9, left=0.2, wspace=0.5, top=0.98, bottom=0.15)
        for i in ('svg', 'png'):
            file = os.path.join(images_folder, f'montecarlo_agile_comparison.{i}')
            plt.savefig(file, transparent=True)
    
def plot_montecarlo_derivative():
    set_font(size=8)
    set_figure_size(
        aspect_ratio=0.5,
        # width=3.3071, aspect_ratio=1.85
    )
    fig, axes = plot_monte_carlo(
        derivative=True, absolute=True, 
        comparison=False, agile=False,
        ncols=3,
        # tickmarks=np.array([
        #     [-3, -2, -1, 0, 1, 2, 3, 4, 5],
        #     [-9, -6, -3,  0, 3, 6, 9, 12, 15],
        #     [-2.0, -1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5, 2],
        #     [-16, -8, 0, 8, 16, 24, 32, 40, 48],
        #     [-400, -300, -200, -100, 0, 100, 200, 300, 400],
        #     [-300, -225, -150, -75, 0, 75, 150, 225, 300]
        # ], dtype=object),
        labels=['DC', 'ICF'],
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 'grey', 
            'orange', 
        ])
    )
    for ax, letter in zip(axes, 'ABCDEFGH'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(1.65, ylb + (yub - ylb) * 0.90, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(
        hspace=0, wspace=0.7,
        top=0.95, bottom=0.1,
        left=0.12, right=0.96
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'montecarlo_derivative.{i}')
        plt.savefig(file, transparent=True)

def plot_montecarlo_absolute():
    set_font(size=8)
    set_figure_size(aspect_ratio=1.05)
    fig, axes = plot_monte_carlo(
        absolute=True, comparison=False, ncols=2,
        expand=0.1, 
        labels=['Sugarcane\nDC', 'Oilcane\nDC',
                'Sugarcane\nICF', 'Oilcane\nICF',
                'Sugarcane &\nSorghum DC', 'Oilcane &\nOil-sorghum DC',
                'Sugarcane &\nSorghum ICF', 'Oilcane &\nOil-sorghum ICF'],
        xrot=90,
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 'grey',
            'orange', 'orange', 'green', 'orange', 'green',
        ])
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJ'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(7.8, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.40, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'montecarlo_absolute.{i}')
        plt.savefig(file, transparent=True)
    
def plot_montecarlo_absolute_YRCP2023():
    set_font(size=8)
    set_figure_size(aspect_ratio=1.05)
    fig, axes = plot_monte_carlo(
        absolute=True, comparison=False, ncols=2,
        expand=0.1, 
        labels=['Sugarcane\nDC', 'Oilcane\nDC',
                'Sugarcane\nICF', 'Oilcane\nICF'],
        xrot=90,
        color_wheel = CABBI_colors.wheel([
            'blue_light', 'green_dirty', 'orange', 'green', 'grey',
            'orange', 'orange', 'green', 'orange', 'green',
        ])
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJ'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(7.8, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.40, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'montecarlo_absolute.{i}')
        plt.savefig(file, transparent=True)
    
def plot_lines_monte_carlo_manuscript(fs=8):
    set_font(size=fs)
    set_figure_size(aspect_ratio=1.2)
    for name, metrics in line_metrics.items():
        fig, axes = plot_lines_monte_carlo(
            metrics=metrics,
            labels=['Sugarcane WT\n0.03 wt % Oil', 
                    'Oilcane 1566\n1.8 wt % Oil', 
                    'Oilcane 1580\n5.4 wt % Oil',
                    'Oilcane Target\n10 wt % Oil'],
            xrot=90,
            # configurations=['O7'],
            # tickmarks=[[0, 3, 6, 9, 12, 15],
            #             [0, 2, 4, 6, 8, 10]], 
            color_wheel = line_color_wheel,
        )
        for ax, letter in zip(axes.flat, 'ABCDEFGHIJKL'):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            plt.text(-0.25, ylb + (yub - ylb) * 0.9, letter, color=letter_color,
                     horizontalalignment='center',verticalalignment='center',
                     fontsize=10, fontweight='bold')
        plt.subplots_adjust(left=0.12, right=0.95, wspace=0, top=0.95, bottom=0.15)
        for i in ('svg', 'png'):
            file = os.path.join(images_folder, f'montecarlo_lines_{name}.{i}')
            plt.savefig(file, transparent=True)
    
def plot_spearman_tea(with_units=None, aspect_ratio=0.8, **kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=aspect_ratio)
    plot_spearman(
        configurations=[
            'O1', 'O1*',
            'O2', 'O2*',
        ],
        labels=[
            'DC', 'Oil-sorghum int., DC',
            'ICF', 'Oil-sorghum int., ICF',
        ],
        kind='TEA',
        with_units=with_units,
        cutoff=0.05,
        **kwargs
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.08)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_tea.{i}')
        plt.savefig(file, transparent=True)

def plot_spearman_tea_YRCP2023(aspect_ratio=0.8, **kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=aspect_ratio)
    plot_spearman_YRCP2023(
        configurations=[
            'O7', 'O9',
        ],
        labels=[
            'DC', 'ICF'
        ],
        kind='TEA',
        top=10,
        **kwargs
    )
    plt.subplots_adjust(left=0.3, right=0.975, top=0.98, bottom=0.08)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_tea_YRCP2023.{i}')
        plt.savefig(file, transparent=True)

def plot_spearman_lca_YRCP2023(aspect_ratio=0.8, **kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=aspect_ratio)
    plot_spearman_YRCP2023(
        configurations=[
            'O7', 'O9',
        ],
        labels=[
            'DC', 'ICF'
        ],
        kind='LCA',
        top=10,
        **kwargs
    )
    plt.subplots_adjust(left=0.35, right=0.975, top=0.98, bottom=0.08)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_lca_YRCP2023.{i}')
        plt.savefig(file, transparent=True)

def plot_spearman_tea_short(**kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=0.65, width=6.6142 * 2/3)
    plot_spearman(
        configurations=[
            'O1', 
            'O2', 
        ],
        labels=[
            'DC', 
            'ICF',
        ],
        kind='TEA',
        with_units=False,
        cutoff=0.05,
        top=5,
        legend=True,
        legend_kwargs={'loc': 'upper left'},
        **kwargs
    )
    plt.subplots_adjust(left=0.35, right=0.975, top=0.98, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_tea.{i}')
        plt.savefig(file, transparent=True)

def plot_spearman_lca_short(with_units=False, aspect_ratio=0.65, **kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=aspect_ratio, width=6.6142 * 2/3)
    plot_spearman(
        configurations=[
            'O1', 
            'O2', 
        ],
        labels=[
            'DC', 
            'ICF',
        ],
        kind='LCA',
        with_units=with_units,
        cutoff=0.025,
        top=5,
        legend=False,
        **kwargs
    )
    plt.subplots_adjust(left=0.35, right=0.975, top=0.98, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_lca.{i}')
        plt.savefig(file, transparent=True)

def plot_spearman_lca(with_units=None, aspect_ratio=0.65, **kwargs):
    set_font(size=8)
    set_figure_size(aspect_ratio=aspect_ratio)
    plot_spearman(
        configurations=[
            'O1', 'O1*',
            'O2', 'O2*',
        ],
        labels=[
            'DC', 'Oil-sorghum int., DC',
            'ICF', 'Oil-sorghum int., ICF',
        ],
        kind='LCA',
        with_units=with_units,
        cutoff=0.025,
        **kwargs
    )
    plt.subplots_adjust(left=0.45, right=0.975, top=0.98, bottom=0.10)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_lca.{i}')
        plt.savefig(file, transparent=True)

def plot_breakdowns(biodiesel_only=False):
    set_font(size=8)
    set_figure_size(aspect_ratio=0.68)
    if biodiesel_only:
        cane.YRCP2023()
        fig, axes = plt.subplots(nrows=1, ncols=3)
        configurations = ('O7', 'O9', 'O8')
        yticks = [-20, 0, 20, 40, 60, 80, 100]
        ytexts = ['']*len(yticks)
        for i, config in enumerate(configurations):
            ax = axes[i]
            plt.sca(ax)
            plot_configuration_breakdown(config, ax=ax, legend=i == 2, biodiesel_only=biodiesel_only)
            if i: 
                plt.yticks(yticks, ytexts)
                plt.ylabel('')
            else:
                plt.yticks(yticks)
        plt.subplots_adjust(left=0.07, right=0.98, wspace=0., top=0.84, bottom=0.31)
        titles = ['(A) Direct\nCogeneration', '(B) Integrated\nCo-Fermentation', '(C) Integrated\nCo-Fermentation & Recovery']
        for ax, letter in zip(axes, titles):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            xlb, xub = plt.xlim()
            plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.2, letter, color=letter_color,
                      horizontalalignment='center',verticalalignment='center',
                      fontsize=10, fontweight='bold')
        for i in ('svg', 'png'):
            name = f'breakdowns_biodiesel_only.{i}'
            file = os.path.join(images_folder, name)
            plt.savefig(file, transparent=True)
    else:
        fig, axes = plt.subplots(nrows=1, ncols=2)
        plt.sca(axes[0])
        c1, c2 = ('O1', 'O2')
        plot_configuration_breakdown(c1, ax=axes[0], legend=False, biodiesel_only=biodiesel_only)
        plt.sca(axes[1])
        plot_configuration_breakdown(c2, ax=axes[1], legend=True, biodiesel_only=biodiesel_only)
        yticks = axes[1].get_yticks()
        plt.yticks(yticks, ['']*len(yticks))
        plt.ylabel('')
        plt.subplots_adjust(left=0.09, right=0.96, wspace=0., top=0.84, bottom=0.31)
        for ax, letter in zip(axes, ['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation']):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            xlb, xub = plt.xlim()
            plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.2, letter, color=letter_color,
                      horizontalalignment='center',verticalalignment='center',
                      fontsize=10, fontweight='bold')
        for i in ('svg', 'png'):
            name = f'breakdowns.{i}'
            file = os.path.join(images_folder, name)
            plt.savefig(file, transparent=True)

# %% KDE


def plot_kde(name, metrics=(GWP_ethanol, MFPP), xticks=None, yticks=None,
             xbox_kwargs=None, ybox_kwargs=None, top_left='',
             top_right='Tradeoff', bottom_left='Tradeoff',
             bottom_right='', fs=None, ticklabels=True, aspect_ratio=1.1):
    set_font(size=fs or 8)
    set_figure_size(width='half', aspect_ratio=aspect_ratio)
    Xi, Yi = [i.index for i in metrics]
    df = get_monte_carlo(name, metrics)
    y = df[Yi].values
    x = df[Xi].values
    sX, sY = [kde_comparison_settings[i] for i in metrics]
    _, xlabel, fx = sX
    _, ylabel, fy = sY
    if fx: x *= fx
    if fy: y *= fy
    ax = bst.plots.plot_kde(
        y=y, x=x, xticks=xticks, yticks=yticks,
        xticklabels=ticklabels, yticklabels=ticklabels,
        xbox_kwargs=xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        ybox_kwargs=ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        aspect_ratio=1.2,
    )
    plt.sca(ax)
    plt.xlabel(xlabel.replace('\n', ' '))
    plt.ylabel(ylabel.replace('\n', ' '))
    bst.plots.plot_quadrants(data=[x, y], text=[top_left, top_right, bottom_left, bottom_right])
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.2, right=0.98,
    )

def plot_kde_fake_scenarios_ethanol_price(name, xticks=None, yticks=None,
             xbox_kwargs=None, ybox_kwargs=None, top_left='',
             top_right='Tradeoff', bottom_left='Tradeoff',
             bottom_right='', fs=None, ticklabels=True, aspect_ratio=1.1):
    set_font(size=fs or 8)
    set_figure_size(width='half', aspect_ratio=aspect_ratio)
    ethanol_price = features.set_ethanol_price
    relative_biodiesel_price = features.set_biodiesel_price
    metrics = (GWP_ethanol, MFPP)
    all_features = (*metrics, ethanol_price, relative_biodiesel_price, 
                    ethanol_production, biodiesel_production)
    Xi, Yi = [i.index for i in metrics]
    df = get_monte_carlo(name, all_features)
    x = df[Xi].values
    y = df[Yi].values
    ethanol_price_range = [0.269, 0.758 * 2]
    biodiesel_price_range = [0.0819, 1.09 * 2]
    settings = kde_metric_settings if '-' in name else  kde_comparison_settings
    sX, sY = [settings[i] for i in metrics]
    _, xlabel, fx = sX
    _, ylabel, fy = sY
    if fx: x *= fx
    if fy: y *= fy
    xs = (x, x)
    ys = tuple([y + (i - df[ethanol_price.index]) * df[ethanol_production.index]
                  + (i + j - df[ethanol_price.index] + df[relative_biodiesel_price.index]) * df[biodiesel_production.index]
                for (i, j) in zip(ethanol_price_range, biodiesel_price_range)])
    ax = bst.plots.plot_kde(
        y=ys, x=xs, xticks=xticks, yticks=yticks,
        xticklabels=ticklabels, yticklabels=ticklabels,
        xbox_kwargs=xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        ybox_kwargs=ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        cmaps=['inferno', 'viridis'],
        aspect_ratio=1.2,
    )
    plt.sca(ax)
    plt.xlabel(xlabel.replace('\n', ' '))
    plt.ylabel(ylabel.replace('\n', ' '))
    bst.plots.plot_quadrants(
        data=None, 
        text=['Market\nCompetitive', 'Market\nCompetitive',
              'Financial\nRisk', 'Financial\nRisk'],
        y=50,
    )
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.15, right=0.98,
    )

def plot_kde_2d(name, metrics=(GWP_biofuel_allocation, MFPP), xticks=None, yticks=None,
                top_left='', top_right='', bottom_left='',
                bottom_right='', xbox_kwargs=None, ybox_kwargs=None, titles=None,
                fs=None, ticklabels=True, rotate_quadrants=0,
                x_center=None, y_center=None, 
                xlabel=None, ylabel=None):
    set_font(size=fs or 8)
    set_figure_size(aspect_ratio=0.6)
    if isinstance(name, str): name = (name,)
    Xi, Yi = [i.index for i in metrics]
    dfs = [get_monte_carlo(i, metrics) for i in name]
    xs = np.array([[df[Xi].values for df in dfs]])
    ys = np.array([[df[Yi].values for df in dfs]])
    if not (xlabel or ylabel):
        sX, sY = [kde_comparison_settings[i] for i in metrics]
        _, xlabel, fx = sX
        _, ylabel, fy = sY
        if fx: xs *= fx
        if fy: ys *= fy
    ticklabels = [True, True] if ticklabels else [False, False]
    axes = bst.plots.plot_kde_2d(
        xs=xs, ys=ys,
        xticks=xticks, yticks=yticks,
        xticklabels=ticklabels, yticklabels=ticklabels,
        xbox_kwargs=2*[xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn)],
        ybox_kwargs=[ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn)],
        aspect_ratio=1.,
    )
    M, N = axes.shape
    text = [top_left, top_right, bottom_left, bottom_right]
    for i in range(M):
        for j in range(N):
            ax = axes[i, j]
            plt.sca(ax)
            if i == M - 1: plt.xlabel(xlabel.replace('\n', ' '))
            if j == 0: plt.ylabel(ylabel.replace('\n', ' '))
            df = dfs[j]
            x = df[Xi]
            y = df[Yi]
            bst.plots.plot_quadrants(data=[x, y], x=x_center, y=y_center, 
                                     text=text, rotate=rotate_quadrants)
    plt.subplots_adjust(
        hspace=0, wspace=0,
        top=0.98, bottom=0.15,
        left=0.1, right=0.98,
    )
    if titles:
        plt.subplots_adjust(
            top=0.90,
        )
        for ax, letter in zip(axes[0, :], titles):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            xlb, xub = plt.xlim()
            plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.17, letter, color=letter_color,
                      horizontalalignment='center', verticalalignment='center',
                      fontsize=10, fontweight='bold')

def plot_unlabeled_feedstock_conventional_comparison_kde(fake_scenarios=True):
    mfpp = kde_comparison_settings[MFPP]
    gwp = kde_comparison_settings[GWP_ethanol]
    gwp_dummy = gwp[1].replace('$_{\\mathrm{economic}}$', '')
    if fake_scenarios: 
        gwp_dummy = gwp_dummy.replace(r"$\Delta$", '')
    kde_comparison_settings[GWP_ethanol] = [
        gwp[0], 
        gwp_dummy, 
        gwp[2]
    ]
    kde_comparison_settings[MFPP] = [
        mfpp[0], 
        'IRR [%]', 
        mfpp[2]
    ]
    (plot_kde_fake_scenarios_ethanol_price if fake_scenarios else plot_kde)(
        'O1' if fake_scenarios else 'O1 - S1',
        yticks=[-75, 0, 75, 150, 225, 300] if fake_scenarios else [-15, 0, 15, 30, 45, 60],
        xticks=[0.25, 0.3, 0.35, 0.4] if fake_scenarios else [-0.12, -0.09, -0.06, -0.03, 0, 0.03, 0.06],
        top_left='Conf. A\nFavored()',
        bottom_right='Conf. B\nFavored()',
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        ticklabels=False,
        aspect_ratio=0.8,
    )
    kde_comparison_settings[GWP_ethanol] = gwp
    kde_comparison_settings[MFPP] = mfpp
    for i in ('svg', 'png'):
        file = os.path.join(
            images_folder, 
            f'prelim_scenario_comparison_kde.{i}' if fake_scenarios else f'prelim_comparison_kde.{i}'
        )
        plt.savefig(file, transparent=True)

def plot_feedstock_conventional_comparison_kde():
    plot_kde(
        'O1 - S1',
        yticks=[-15, 0, 15, 30, 45, 60],
        xticks=[-0.12, -0.09, -0.06, -0.03, 0, 0.03, 0.06],
        top_left='Oilcane Favored',
        bottom_right='Sugarcane\nFavored',
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'feedstock_conventional_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_feedstock_cellulosic_comparison_kde(fs=None):
    plot_kde(
        'O2 - S2',
        yticks=[-15, 0, 15, 30, 45, 60],
        xticks=[-0.2, -0.1, 0, 0.1, 0.2],
        top_left='Oilcane\nFavored()',
        bottom_right='Sugarcane\nFavored()',
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        fs=fs,
        # fx=1000.,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'feedstock_cellulosic_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_feedstock_comparison_kde(fs=None):
    plot_kde_2d(
        ('O1 - S1', 'O2 - S2'),
        yticks=[[-18, 0, 18, 36, 54]],
        xticks=[[-0.12, -0.09, -0.06, -0.03, 0, 0.03, 0.06],
                [-0.06, -0.03, 0., 0.03, 0.06, 0.09]],
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Oilcane\nFavored()',
        bottom_right='\nSugarcane\nFavored()',
        titles=['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'],
        fs=fs,
    )
    plt.subplots_adjust(
        wspace=0,
        
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'feedstock_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_feedstock_comparison_kde_2023(fs=None):
    plot_kde_2d(
        ('O1|constant biomass yield - S1.WT', 'O2|constant biomass yield - S2.WT'),
        yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='TCI\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Oilcane\nFavored()',
        bottom_right='\nSugarcane\nFavored()',
        titles=['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'],
        metrics=[TCI, MFPP],
        fs=fs,
    )
    plt.subplots_adjust(
        wspace=0,
        
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'feedstock_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_microbial_oil_kde_2023(fs=None):
    plot_kde_2d(
        ('O7.WT', 'O9.WT'),
        yticks=[[0, 800, 1600, 2400, 3200, 4000, 4800]],
        xticks=[[0, 0.5, 1.0, 1.5, 2.0, 2.5],
                [0, 0.5, 1.0, 1.5, 2.0, 2.5]],
        top_right='GWP\nTradeoff()',
        bottom_left='Yield\nTradeoff()',
        top_left='Microbial oil\nFavored()',
        bottom_right='Soybean oil\nFavored()',
        titles=['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'],
        metrics=[GWP_biodiesel_allocation, biodiesel_yield],
        x_center=1.29,
        y_center=617,
        xlabel="GWP$_{\\mathrm{energy}}$\n" f"[{GWP_units_L}]",
        ylabel=f"Biodiesel yield\n[{format_units('L/ha/y')}]",
        fs=fs,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'microbial_oil_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_configuration_comparison_kde(fs=None):
    plot_kde(
        'O2 - O1',
        yticks=sorted([-1 * i for i in [-39, -26, -13, 0, 13, 26]]),
        xticks=sorted([-1 * i for i in [-0.06, -0.04, -0.02, 0, 0.02, 0.04]]),
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='ICF\nFavored()',
        bottom_right='DC Favored()',
        fs=fs,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'configuration_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_microbial_oil_feedstock_comparison_kde(fs=None):
    plot_kde_2d(
        ('O7 - O7.WT', 'O9 - O9.WT'),
        yticks=[[-90, -60, -30, 0, 30, 60, 90]],
        xticks=[[-1.5, -1, -0.5, 0, 0.5, 1],
                [-3, -2, -1, 0, 1, 2]],
        # yticks=sorted([-1 * i for i in [-39, -26, -13, 0, 13, 26]]),
        # xticks=sorted([-1 * i for i in [-0.06, -0.04, -0.02, 0, 0.02, 0.04]]),
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Plant oil\nFavored()',
        bottom_right='Microb. oil\nFavored()',
        titles=['(A) Direct Cogeneration',
                '(B) Integrated Co-Fermentation'],
        fs=fs,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'microbial_oil_feedstock_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_microbial_oil_bioethanol_comparison_kde(fs=None):
    plot_kde_2d(
        ('O7.WT - S1.WT', 'O9.WT - S2.WT'),
        yticks=[[-150, -120, -90, -60, -30, 0, 30]],
        xticks=[[-0.24, -0.12, 0, 0.12, 0.24, .36],
                [-0.5, 0, 0.5, 1, 1.5, 2.]],
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Microb. Oil\nFavored()',
        bottom_right='Bioethanol\nFavored()',
        titles=['(A) Direct Cogeneration',
                '(B) Integrated Co-Fermentation'],
        fs=fs,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'microbial_oil_bioethanol_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_separated_configuration_comparison_kde():
    plot_kde_2d(
        ('O1', 'O2'),
        yticks=[[-20, 0, 20, 40, 60]],
        xticks=[
            [0, 0.5, 1, 1.5],
            [0, 2, 4, 6, 8, 10]
        ],
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='DC Favored()',
        bottom_right='ICF\nFavored()',
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'separated_configuration_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_crude_configuration_comparison_kde():
    plot_kde_2d(
        ('O1 - O3', 'O2 - O4'),
        yticks=[[0, 10, 20, 30, 40]],
        xticks=[
            [-0.3, -0.24, -0.18, -0.12, -0.06, 0],
            [-0.15, -0.12, -0.09, -0.06, -0.03, 0]
        ],
        top_right='GWP\nTradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Biodiesel\nProduction Favored()',
        bottom_right='Crude Oil\nProduction Favored()',
        titles=['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'],
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'crude_configuration_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_agile_comparison_kde(fs=None):
    plot_kde_2d(
        ('O1* - O1', 'O2* - O2'),
        metrics=[TCI, MFPP],
        yticks=[[0, 3, 6, 9, 12, 15]], 
        xticks=[
            [-100, -80, -60, -40, -20, 0],
            [-150, -125, -100, -75, -50, -25, 0],
        ],
        top_right='TCI-Tradeoff()',
        bottom_left='MFPP\nTradeoff()',
        top_left='Oil-sorghum\nIntegration Favored()',
        bottom_right='Cane-only\nFavored()',
        xbox_kwargs=dict(light=CABBI_colors.green_dirty.RGBn, 
                         dark=CABBI_colors.green_dirty.shade(60).RGBn),
        titles=['(A) Direct Cogeneration', '(B) Integrated Co-Fermentation'],
        fs=fs,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'agile_conventional_comparison_kde.{i}')
        plt.savefig(file, transparent=True)

def plot_open_comparison_kde(overlap=False):
    metrics = [MFPP, TCI, GWP_ethanol, biodiesel_production]
    df_conventional_oc = get_monte_carlo('O1', metrics)
    df_cellulosic_oc = get_monte_carlo('O2', metrics)
    df_conventional_sc = get_monte_carlo('S1', metrics)
    df_cellulosic_sc = get_monte_carlo('S2', metrics)
    MFPPi = MFPP.index
    TCIi = TCI.index
    if overlap:
        ys = np.zeros([1, 2], dtype=object)
        xs = np.zeros([1, 2], dtype=object)    
        ys[0, 0] = (df_conventional_oc[MFPPi], df_cellulosic_oc[MFPPi])
        ys[0, 1] = (df_conventional_sc[MFPPi], df_cellulosic_sc[MFPPi])
        xs[0, 0] = (df_conventional_oc[TCIi], df_cellulosic_oc[TCIi])
        xs[0, 1] = (df_conventional_sc[TCIi], df_cellulosic_sc[TCIi])
        yticks = [[-30, -15, 0, 15, 30, 45, 60, 75]]
        xticks = 2*[[200, 300, 400, 500, 600]]
    else:
        ys = np.array([
            [df_conventional_oc[MFPPi], df_conventional_sc[MFPPi]],
            [df_cellulosic_oc[MFPPi], df_cellulosic_sc[MFPPi]]
        ])
        xs = np.array([
            [df_conventional_oc[TCIi], df_conventional_sc[TCIi]],
            [df_cellulosic_oc[TCIi], df_cellulosic_sc[TCIi]]
        ])
        yticks = 2*[[-30, -15, 0, 15, 30, 45, 60, 75]]
        xticks = 2*[[200, 300, 400, 500, 600]]
    bst.plots.plot_kde_2d(
        ys=ys, xs=xs, xticks=xticks, yticks=yticks,
        xbox_kwargs=[dict(position=1), dict(position=1)],
        ybox_kwargs=[dict(position=0), dict(position=0)],
    )

#%%  General Monte Carlo box plots

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

def monte_carlo_box_plot(data, positions, light_color, dark_color, width=None, 
                         hatch=None, outliers=False, **kwargs):
    if width is None: width = 0.8
    if outliers:
        flierprops = {'marker':'D',
                  'markerfacecolor': light_color,
                  'markeredgecolor': dark_color,
                  'markersize':3}
    else:
        flierprops = {'marker':''}
    bp = plt.boxplot(
        x=data, positions=positions, patch_artist=True,
        widths=width, whis=[5, 95],
        boxprops={'facecolor':light_color,
                  'edgecolor':dark_color},
        medianprops={'color':dark_color,
                     'linewidth':1.5},
        flierprops=flierprops,
        **kwargs
    )
    if hatch:
        for box in bp['boxes']:
            box.set(hatch = hatch)

def plot_monte_carlo(derivative=False, absolute=True, comparison=True,
                     configuration_names=None, comparison_names=None,
                     metrics=None, labels=None, tickmarks=None, agile=True, 
                     ncols=1, expand=None, step_min=None,
                     agile_only=False, xrot=None,
                     color_wheel=None, axes_box=None):
    if derivative:
        default_configuration_names = ['O1', 'O2']
        default_comparison_names = ['O2 - O1']
        metric_info = mc_derivative_metric_settings
        default_metrics = list(metric_info)
    else:
        default_configuration_names = configuration_names[:-2]
        default_comparison_names = comparison_names
        if comparison:
            metric_info = mc_comparison_settings
        else:
            metric_info = mc_metric_settings
        if agile_only:
            default_configuration_names = [i for i in default_configuration_names if '*' in i]
            default_comparison_names = [i for i in default_comparison_names if '*' in i]
            default_metrics = ['MFPP', 'TCI', 'production']
        else:
            default_metrics = list(metric_info)
    if configuration_names is None: configuration_names = default_configuration_names
    if comparison_names is None: comparison_names = default_comparison_names
    if metrics is None: metrics = default_metrics
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
    rows, ylabels, factors = zip(*[metric_info[i] for i in metrics])
    factors = [(i, j) for i, j in enumerate(factors) if j is not None]
    if color_wheel is None: color_wheel = CABBI_colors.wheel()
    N_rows = len(rows)
    if axes_box is None:
        nrows = int(round(N_rows / ncols))
        fig, axes_box = plt.subplots(ncols=ncols, nrows=nrows)
        plt.subplots_adjust(wspace=0.45)
    else:
        nrows = N_rows
        fig = None
    axes = axes_box.transpose()
    axes = axes.flatten()
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
                    data=arr[:, i], positions=[position + (i-(N-1)/2)*width], 
                    light_color=color.RGBn, 
                    dark_color=color.shade(60).RGBn,
                    width=boxwidth,
                    hatch=getattr(color, 'hatch', None),
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
        
    if tickmarks is None: 
        tickmarks = [
            bst.plots.rounded_tickmarks_from_data(
                i, step_min=step_min, N_ticks=8, lb_max=0, center=0,
                f=roundsigfigs, expand=expand,
                f_min=lambda x: np.percentile(x, 5),
                f_max=lambda x: np.percentile(x, 95),
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
            xticklabels = xtext if (ax in axes_box[-1] or i == N_rows - 1) else []
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
        )
    for i in range(N_rows, nrows * ncols):
        ax = axes[i]
        plt.sca(ax)
        plt.axis('off')
    if fig is None:
        fig = plt.gcf()
    else:
        plt.subplots_adjust(hspace=0)
    fig.align_ylabels(axes)
    return fig, axes

def plot_competitive_biomass_yield_across_oil_content(
        fs=None,
    ):
    if fs is None: fs = 8
    set_figure_size(aspect_ratio=0.8)
    fig, axes = plt.subplots(ncols=1, nrows=2)
    set_font(size=fs)
    ax1, ax2 = axes.flatten()
    plt.sca(ax1)
    xs, fys = _plot_competitive_biomass_yield_across_oil_content('O7')
    plt.xlim([0, 15])
    plt.ylim([0, 40])
    bst.plots.style_axis(
        xticks=[0, 3, 6, 9, 12, 15],
        xticklabels=False,
        ytick0=False,
        ytickf=True,
    )
    plt.text(
        15 - 0.25, 36, '(A) Direct Cogeneration', weight='bold',
        ha='right', va='center', fontsize=10,
        c=bst.plots.title_color,
    )
    _add_lines_biomass_yield_vs_oil_content(*fys)
    plt.sca(ax2)
    xs, fys = _plot_competitive_biomass_yield_across_oil_content('O9')
    plt.xlim([0, 15])
    plt.ylim([0, 40])
    plt.xlabel('Oil content [dw %]')
    bst.plots.style_axis(
        xticks=[0, 3, 6, 9, 12, 15],
        ytick0=False,
        ytickf=True,
    )
    plt.text(
        15 - 0.25, 36, '(B) Integrated Co-Fermentation', weight='bold', 
        ha='right', va='center', fontsize=10,
        c=bst.plots.title_color,
    )
    _add_lines_biomass_yield_vs_oil_content(*fys)
    plt.subplots_adjust(hspace=0.05, left=0.125, right=0.96, bottom=0.10, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'competitive_biomass_yield_MCAC.{i}')
        plt.savefig(file, transparent=True)

def _add_lines_biomass_yield_vs_oil_content(f5, f50, f95):
    df = cane.get_composition_data()
    points = {i: None for i in ('WT', '1566', '1580', 'Target')}
    for name, color in zip(df.index, line_color_wheel):
        line = df.loc[name]
        name = str(name)
        if name not in points: continue
        oil = line['Stem oil (dw)']['Mean'] * 100
        biomass = line['Biomass yield (dry MT/ha)']['Mean']
        points[name] = Point(oil, biomass, color.shade(50).RGBn)
        bst.plots.plot_scatter_points(
            [oil], [biomass], marker='o', s=20, color=color.RGBn,
            edgecolor=edgecolor, clip_on=False, zorder=3,
            linewidth=0.6,
        )
    arrow = r'$\rightarrow$'
    name = 'WT'
    pnt = points[name]
    y50_target = f50(pnt.x)
    x_units = f"{format_units('dw')} % oil"
    y_units = f"{format_units('DMT/ha/y')}"
    text = (
        f"Sugarcane {name}, {pnt.x:.2g} {x_units}\n"
        f"Yield: {pnt.y:.1f} {y_units}"
    )
    bst.plots.annotate_point(
        text, pnt.x, pnt.y, horizontalalignment='left', 
        textcolor=pnt.c, linecolor=pnt.c,
        dx=0.3, dy=7, dx_text=-0.1, dy_text=0, fontsize=8,
    )
    
    name = '1566'
    pnt = points[name]
    y5_target, y50_target, y95_target = f5(pnt.x), f50(pnt.x), f95(pnt.x)
    ax = plt.gca()
    ax.annotate('',
        xy=(pnt.x, y50_target + 0.25), 
        xytext=(pnt.x, pnt.y),
        arrowprops=dict(arrowstyle="->", color=pnt.c),
    )
    text = (
        f"Oilcane {name}, {pnt.x:.2g} {x_units}\n"
        f"Yield: {pnt.y:.1f} {arrow} {y50_target:.1f} [{y5_target:.1f}, {y95_target:.1f}] {y_units}\n"
    )
    bst.plots.annotate_point(
        text, pnt.x, pnt.y, horizontalalignment='left', verticalalignment='center',
        textcolor=pnt.c, linecolor=pnt.c,
        dx=2.5, dy=13.5, dx_text=0.1, dy_text=0, fontsize=8,
    )
    
    name = '1580'
    pnt = points[name]
    y5_target, y50_target, y95_target = f5(pnt.x), f50(pnt.x), f95(pnt.x)
    ax = plt.gca()
    ax.annotate('',
        xy=(pnt.x, y50_target + 0.25), 
        xytext=(pnt.x, pnt.y),
        arrowprops=dict(arrowstyle="->", color=pnt.c),
    )
    text = (
        f"Oilcane {name}, {pnt.x:.2g} {x_units}\n"
        f"Yield: {pnt.y:.1f} {arrow} {y50_target:.1f} [{y5_target:.1f}, {y95_target:.1f}] {y_units}\n"
    )
    bst.plots.annotate_point(
        text, pnt.x, pnt.y, horizontalalignment='left',
        verticalalignment='center',
        textcolor=pnt.c, linecolor=pnt.c,
        dx=0.25, dy=0, dx_text=0.1, dy_text=0, fontsize=8,
    )
    
    name = 'Target'
    pnt = points[name]
    y5_target, y50_target, y95_target = f5(pnt.x), f50(pnt.x), f95(pnt.x)
    ax = plt.gca()
    text = (
        f"Oilcane {name}, {pnt.x:.2g} {x_units}\n"
        f"Yield: {pnt.y:.1f} {y_units}"
    )
    bst.plots.annotate_point(
        text, pnt.x, pnt.y, horizontalalignment='center',
        verticalalignment='bottom',
        textcolor=pnt.c, linecolor=pnt.c,
        dx=0, dy=4.6, dx_text=0, dy_text=0.1, fontsize=8,
    )

def _plot_theoretical_competitive_biomass_yield_across_oil_content():
    f = 0.93 # Maximum fraction of fiber/sugar
    Yb = 25.62 # Baseline biomass yield
    Ym = 0.32 # Microbial oil yield
    x = np.linspace(0, 0.15) # Fraction of oil
    Yo = f * Yb * Ym / ((f - x) * Ym + x) # Theoretical oilcane yield
    # print(
    #     ((1 - Yo / Yb)[1:] / x[1:]).mean()
    # )
    plt.plot(100 * x, Yo, color=bst.plots.title_color, lw=1, zorder=0,
             linestyle='--')
    x = 0.14
    y = f * Yb * Ym / ((f - x) * Ym + x) # Theoretical oilcane yield
    bst.plots.annotate_point(
        'Theoretical\nmaximum\ntargets', 100 * x, y,
        dy=3, dy_text=0.5, fontsize=8,
        linecolor=bst.plots.title_color,
        textcolor=bst.plots.title_color,
        verticalalignment='bottom',
    )
    
def _plot_competitive_biomass_yield_across_oil_content(
        configuration,
    ):
    file = monte_carlo_file(configuration, across_lines=False, across_oil_content='oilcane vs sugarcane')
    df = pd.read_excel(file, sheet_name=features.competitive_biomass_yield.short_description, index_col=0)
    df = df.dropna()
    oil_content = np.array(df.columns) * 100
    plt.ylabel(f"Competitive biomass yield\n[{format_units('DMT/ha/y')}]")
    _plot_theoretical_competitive_biomass_yield_across_oil_content()
    biomass_yields = bst.plots.plot_montecarlo_across_coordinate(
        oil_content, df, 
        fill_color=CABBI_colors.green_dirty.tint(50).RGBn,
        median_color=CABBI_colors.green_dirty.shade(10).RGBn,
        p5_color=CABBI_colors.green_dirty.RGBn,
        smooth=1,
    )
    index = [0, 2, 4]
    f5, f50, f95 = [interpolate.interp1d(oil_content, biomass_yields[i]) for i in index]
    bst.plots.annotate_point(
        'Financially\ncompetitive\ntargets', 13.5, f50(13.5),
        dy=-4, dy_text=-0.5, fontsize=8,
        linecolor=CABBI_colors.green_dirty.shade(40).RGBn,
        textcolor=CABBI_colors.green_dirty.shade(40).RGBn,
        verticalalignment='top',
    )
    return (oil_content, (f5, f50, f95))

def plot_competitive_microbial_oil_yield_across_oil_content(
        configuration=None,
    ):
    if configuration is None: configuration = 'O7'
    file = monte_carlo_file(configuration, across_lines=False, across_oil_content='microbial oil vs bioethanol')
    df = pd.read_excel(file, sheet_name=features.competitive_microbial_oil_yield.short_description, index_col=0)
    df = df.dropna(axis=0)
    fig = plt.figure()
    oil_fraction = np.array(df.columns) * 100
    plt.ylabel('Competitive microbial oil yield [wt %]')
    bst.plots.plot_montecarlo_across_coordinate(oil_fraction, df)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'competitive_microbial_oil_yield_MCAC_{configuration}.{i}')
        plt.savefig(file, transparent=True)

def plot_lines_monte_carlo(
        configurations=None, metrics=None, labels=None, tickmarks=None, 
        expand=None, step_min=None, xrot=None, color_wheel=None,
    ):
    if configurations is None: configurations = ['O7', 'O9']
    if metrics is None: metrics = list(mc_line_metric_settings)
    df = cane.get_composition_data()
    columns = df.index
    rows, ylabels = zip(*[mc_line_metric_settings[i] for i in metrics])
    if color_wheel is None: color_wheel = line_color_wheel
    nrows = len(rows)
    ncols = len(configurations)
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows)
    axes = axes.reshape([nrows, ncols])
    plt.subplots_adjust(wspace=0.45)
    nmc = len(columns)    
    xtext = labels or columns
    N_marks = len(xtext)
    xticks = tuple(range(N_marks))
    
    def get_data(metric, line, configuration):
        df = get_line_monte_carlo(line, configuration, metric)
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
                boxplot = bst.plots.plot_montecarlo(
                    data=arr[:, i], positions=[position + (i-(N-1)/2)*width], 
                    light_color=color.RGBn, 
                    dark_color=color.shade(60).RGBn,
                    width=boxwidth,
                    hatch=getattr(color, 'hatch', None),
                )
                plots.append(boxplot)
            return plots
        else:
            color = color_wheel.next()
            return bst.plots.plot_montecarlo(
                data=arr, positions=[position], 
                light_color=color.RGBn, 
                dark_color=color.shade(60).RGBn,
                hatch=getattr(color, 'hatch', None),
                width=0.618,
            )
    
    data = np.zeros([nrows, nmc], dtype=object)
    data[:] = [[np.array([get_data(i, j, k) for k in configurations]) for j in columns] for i in rows]
    
    if tickmarks is None: 
        tickmarks = [
            bst.plots.rounded_tickmarks_from_data(
                i, step_min=step_min, N_ticks=6, lb_max=0, center=0,
                f=roundsigfigs, expand=expand,
                f_min=lambda x: np.percentile(x, 0.1),
                f_max=lambda x: np.percentile(x, 99.9),
            ) 
            for i in data
        ]

    xf = len(columns) - 0.5
    for i in range(nrows):
        for k in range(ncols):
            color_wheel.restart()
            for j in range(nmc):
                ax = axes[i, k]
                plt.sca(ax)
                plt.xlim(-0.5, xf)
                plot(data[i, j][k], j)
                if k == 0: plt.ylabel(ylabels[i])
    
    titles = ['Direct Cogeneration', 'Integrated Co-Fermentation']
    soybean_color = CABBI_colors.grey.shade(30).RGBn
    MBSP_soybean = 617
    GWP_soybean = 1.29
    for i in range(nrows):
        for j in range(ncols):
            ax = axes[i, j]
            plt.sca(ax)
            if i == 0: ax.set_title(titles[j], fontsize=10, weight='bold', color=letter_color)
            yticks = tickmarks[i]
            plt.ylim([yticks[0], yticks[1]])
            if yticks[0] < 0.:
                bst.plots.plot_horizontal_line(0, color=CABBI_colors.black.RGBn,
                                               lw=0.8, linestyle='--')
            if metrics[i] == 'MBSP':
                if j == 0:
                    lb = 0.80
                    ub = 1.63
                    plt.fill_between([-1, ncols + 2], lb, ub,
                                     color=CABBI_colors.grey.tint(30).RGBn,
                                     linewidth=1.0,
                                     zorder=0)
                    plt.text(1.5, (lb + ub) * 0.5, 'Market price (D4 RIN)',
                             color=letter_color, 
                             verticalalignment='center',
                             horizontalalignment='center',
                             fontweight='bold',
                             fontsize=8,
                             zorder=1)
                else:
                    lb = 1.16
                    ub = 2.16
                    plt.fill_between([-1, ncols + 2], lb, ub,
                                     color=CABBI_colors.grey.tint(30).RGBn,
                                     linewidth=1.0,
                                     zorder=0)
                    plt.text(1.5, (lb + ub) * 0.5, 'Market price (D3 RIN)',
                             color=letter_color, 
                             fontweight='bold',
                             verticalalignment='center',
                             horizontalalignment='center',
                             fontsize=8,
                             zorder=1)
            elif metrics[i] == 'Biodiesel yield':
                if j == 0:
                    ax.annotate('Soybean',
                        xy=(1.5, MBSP_soybean), 
                        xytext=(1.6, MBSP_soybean * 2.5),
                        arrowprops=dict(arrowstyle="->", color=soybean_color),
                        color=soybean_color,
                        fontsize=8,
                        fontweight='bold',
                    )
                bst.plots.plot_horizontal_line(MBSP_soybean, color=soybean_color,
                                               linestyle='-', zorder=0)
            elif metrics[i] == 'GWP economic':
                bst.plots.plot_horizontal_line(GWP_soybean, color=soybean_color,
                                               linestyle='-', zorder=0)
                if j == 0:
                    ax.annotate('Soybean',
                        xy=(0.5, GWP_soybean), 
                        xytext=(0.6, GWP_soybean * 1.5),
                        arrowprops=dict(arrowstyle="->", color=soybean_color),
                        color=soybean_color,
                        fontsize=8,
                        fontweight='bold',
                    )
                    
            xticklabels = xtext if i == nrows - 1 else []
            yticklabels = j == 0
            bst.plots.style_axis(ax,  
                xticks = xticks,
                yticks = yticks,
                xticklabels= xticklabels, 
                yticklabels=yticklabels,
                ytick0=False,
                ytickf=False,
                offset_xticks=False,
                xrot=xrot,
            )
    if fig is None:
        fig = plt.gcf()
    else:
        plt.subplots_adjust(hspace=0)
        plt.subplots_adjust(wspace=0)
    fig.align_ylabels(axes)
    return fig, axes

#%% Spearman

def plot_spearman(configurations, labels=None, metric=None, 
                  kind=None, with_units=None, legend=None, 
                  legend_kwargs=None, **kwargs):
    if kind is None: kind = 'TEA'
    if with_units is None: with_units = True
    if legend is None: legend = True
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
    stream_price = format_units('USD/L')
    USD_MT = format_units('USD/MT')
    ng_price = format_units('USD/m^3')
    electricity_price = format_units('USD/kWhr')
    operating_days = format_units('day/yr')
    capacity = format_units('10^6 MT/yr')
    titer = format_units('g/L')
    productivity = format_units('g/L/hr')
    material_GWP = '$\\mathrm{kg} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{kg}^{-1}$'
    feedstock_GWP = '$\\mathrm{g} \\cdot \\mathrm{CO}_{2}\\mathrm{eq} \\cdot \\mathrm{kg}^{-1}$'
    index, ignored_list = zip(*[
         ('Crushing mill oil recovery [60 $-$ 90 %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Saccharification oil recovery [70 $-$ 90 %]', ['S2', 'S1', 'S2*', 'S1*', 'O1', 'O1*']),
        (f'Cane operating days [120 $-$ 180 {operating_days}]', []),
        (f'Sorghum operating days [30 $-$ 60 {operating_days}]', ['S2', 'S1', 'O1', 'O2']),
        (f'Crushing capacity [1.2 $-$ 2.0 {capacity}]', []),
        (f'Cellulosic ethanol price [0.655, 1.02, 1.03 {stream_price}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Advanced ethanol price [0.414, 0.599, 0.705 {stream_price}]', []),
        (f'Biodiesel price [0.800, 0.801, 1.57 {stream_price}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Natural gas price [0.105, 0.122, 0.175 {ng_price}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Electricity price [0.0583, 0.065, 0.069 {electricity_price}]', ['S2', 'O2', 'S2*', 'O2*']),
         ('IRR [10 $-$ 15 %]', []),
        (f'Crude glycerol price [100 $-$ 220 {USD_MT}]', ['S2', 'S1', 'S2*', 'S1*']),
        (f'Pure glycerol price [488 $-$ 812 {USD_MT}]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Saccharification reaction time [54 $-$ 90 hr]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Cellulase price [159 $-$ 265 {USD_MT}]', ['S1', 'O1', 'S1*', 'O1*']),
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
         ('Cane oil content [5 $-$ 15 dw %]', ['S2', 'S1', 'S2*', 'S1*']),
         ('Relative sorghum oil content [-3 $-$ 0 dw %]', ['S2', 'S1', 'S2*', 'S1*', 'O2', 'O1']),
         ('TAG to FFA conversion [17.25 $-$ 28.75 % theoretical]', ['S1', 'O1', 'S1*', 'O1*']),
        # TODO: change lower upper values to baseline +- 10%
        (f'Feedstock GWPCF [26.3 $-$ 44.0 {feedstock_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Methanol GWPCF [0.338 $-$ 0.563 {material_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Pure glycerine GWPCF [1.25 $-$ 2.08 {material_GWP}]', ['S1', 'S2', 'S1*', 'S2*']),
        (f'Cellulase GWPCF [6.05 $-$ 10.1 {material_GWP}]', ['S1', 'O1', 'S1*', 'O1*']),
        (f'Natural gas GWPCF [0.297 $-$ 0.363 {material_GWP}]', ['S1', 'O1', 'S1*', 'O1*']),
         ('Corporate income tax [21 $-$ 28 %]', []),
    ])
    if not with_units: index = [i.split(' [')[0] for i in index]
    ignored_dct = {
        'S1': [],
        'O1': [],
        'S2': [],
        'O2': [],
        'S1*': [],
        'O1*': [],
        'S2*': [],
        'O2*': [],
        'O7': [],
        'O8': [],
    }
    for i, ignored in enumerate(ignored_list):
        for name in ignored: ignored_dct[name].append(i)
        index_name = index[i]
        if kind == 'LCA':
            for term in ('cost', 'price', 'IRR', 'time', 'capacity', 'operating'):
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
    fig, ax = bst.plots.plot_spearman_2d(rhos, index=index,
                                         color_wheel=color_wheel,
                                         name=metric_name,
                                         **kwargs)
    if legend:
        if legend_kwargs is None:
            legend_kwargs = {'loc': 'lower left'}
        plt.legend(
            handles=[
                mpatches.Patch(
                    color=color_wheel[i].RGBn, 
                    label=labels[i] if labels else format_name(configurations[i])
                )
                for i in range(len(configurations))
            ], 
            **legend_kwargs,
        )
    return fig, ax

def plot_spearman_YRCP2023(configurations=None, labels=None, 
                  kind=None, with_units=None, legend=None, 
                  legend_kwargs=None, **kwargs):
    cane.YRCP2023()
    if configurations is None: configurations = ('O7', 'O9')
    if kind is None: kind = 'TEA'
    if with_units is None: with_units = True
    if legend is None: legend = True
    if kind == 'TEA':
        metric = MBSP
        metric_name = metric.name
    elif kind == 'LCA':
        metric = GWP_biofuel_allocation
        metric_name = r'GWP$_{\mathrm{energy}}$'
    else:
        raise ValueError(f"invalid kind '{kind}'")
    rhos = []
    for configuration in configurations:
        names = get_YRCP2023_spearman_names(configuration, kind)
        file = spearman_file(configuration)
        try: 
            df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
        except: 
            warning = RuntimeWarning(f"file '{file}' not found")
            warn(warning)
            continue
        s = df[metric.index]
        rhos.append(
            {j: s[i] for i, j in names.items()}
        )
    rhos = pd.DataFrame(rhos).T
    rhos.fillna(0, inplace=True)
    index = rhos.index
    rhos = [rhos[i] for i in rhos]
    color_wheel = [GG_colors.orange, GG_colors.blue]
    fig, ax = bst.plots.plot_spearman_2d(rhos, index=index,
                                         color_wheel=color_wheel,
                                         name=metric_name,
                                         **kwargs)
    if legend:
        if legend_kwargs is None:
            legend_kwargs = {'loc': 'lower right'}
        plt.legend(
            handles=[
                mpatches.Patch(
                    color=color_wheel[i].RGBn, 
                    label=labels[i] if labels else format_name(configurations[i])
                )
                for i in range(len(configurations))
            ], 
            **legend_kwargs,
        )
    return fig, ax

# %% Other

def plot_configuration_breakdown(name, across_coordinate=False, 
                                 biodiesel_only=None,
                                 **kwargs):
    oc = cane.Biorefinery(name)
    if biodiesel_only:
        colors = [area_colors_biodiesel[i.name].RGBn for i in oc.unit_groups]
        hatches = [area_hatches_biodiesel[i.name] for i in oc.unit_groups]
        for i in oc.unit_groups:
            if len(i.metrics) == 4: continue
            i.metrics = i.metrics[:4]
        units = sum([i.units for i in oc.unit_groups], [])
        joint_group = bst.UnitGroup(None, units)
        joint_group.autofill_metrics(shorthand=True, material_cost=False)
    else:
        colors = [area_colors[i.name].RGBn for i in oc.unit_groups]
        hatches = [area_hatches[i.name] for i in oc.unit_groups]
        joint_group = None
    if across_coordinate:
        return bst.plots.plot_unit_groups_across_coordinate(
            cane.set_cane_oil_content,
            [5, 7.5, 10, 12.5],
            'Feedstock oil content [dw %]',
            oc.unit_groups,
            colors=colors,
            hatches=hatches,
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
            if len(i.metrics) == 5: i.metrics[4].name = 'Mat.\ncost'
        
        return bst.plots.plot_unit_groups(
            oc.unit_groups,
            colors=colors,
            hatches=hatches,
            format_total=format_total,
            fraction=True,
            joint_group=joint_group,
            legend_kwargs=dict(
                loc='lower right' if biodiesel_only else 'lower center',
                ncol=4,
                bbox_to_anchor=(0.9 if biodiesel_only else 0, -0.52),
                labelspacing=1.5, handlelength=2.8,
                handleheight=1, scale=0.8,
            ),
            **kwargs,
        )
    
    
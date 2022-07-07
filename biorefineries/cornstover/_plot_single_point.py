# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:44:23 2019

@author: yoelr
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from biosteam.plots import (plot_scatter_points, plot_horizontal_line,
                           plot_montecarlo, plot_vertical_line)
import os
from biorefineries.cornstover.model import cornstover_model
import thermosteam as tmo
from thermosteam.utils import colors
from colorpalette import Color
                    
red = Color('red', '#ed5a6a')
red_shade = Color('red_shade', '#d45063')
orange = Color('orange', '#f98f60')
orange_shade = Color('orange', '#de7e55')
green = Color('green', '#7ac083')
green_shade = Color('green_shade', '#6ba873')
blue = Color('blue', '#60c1cf')
blue_shade = Color('blue_shade', '#55a8b5')
purple = Color('purple', '#a381ba')
purple_shade = Color('purple_shade', '#8b6f9d')
black = Color('black', '#3a4551')

# %% Setup of subplots

# light_color = colors.blue_tint.RGBn
# dark_color = colors.blue_shade.RGBn
# dot_color = colors.purple_shade.RGBn
nums = tuple(range(1, 9))
tmo.utils.plots.set_font(size=8)
tmo.utils.plots.set_figure_size(width=6.6142, aspect_ratio=0.45)
fig, axes = plt.subplots(ncols=1, nrows=2, constrained_layout=True, gridspec_kw=dict(height_ratios=[1, 4]))
magnitude_ax, metric_over_benchmark_ax = axes
plt.sca(metric_over_benchmark_ax)

# %% Plot MESP and ethanol sales

positions_other = (16, 17, 18)
other_index = [('Biorefinery', 'Steam demand [kg/hr]'),
               ('Biorefinery', 'Ethanol production [kg/hr]'),
               ('Biorefinery', 'Minimum ethanol selling price [USD/gal]')]
data = cornstover_model.metrics_at_baseline()
other_data = np.array(data[other_index])
other_data[0] *= 100./234784.
other_data[1] *= 100./22273.
other_data[2] *= 100./2.15
bx_other = plt.bar(positions_other, other_data, 0.5,
        align='center',
        color=blue.RGBn,
        edgecolor=blue_shade.RGBn)

# %% Plot electricity

plot_vertical_line(15.5, color=colors.grey_tint.RGBn)
positions_electricity = tuple(range(9))
areas = [f'Area {i}00' for i in nums]
units = 'MW'
positions = np.arange(0, 9)
electricity_cols = [(i, 'Electricity [MW]') for i in areas]
humbird_electricity = 41 * np.array((0.02, 0.14, 0.06, 0.05,
                                     0.18, 0.003, 0.03, 0.08,
                                     0.44))
electricity_data = data[electricity_cols]  #/ humbird_electricity
electricity_data[('Biorefinery', 'Excess electricity [MW]')] = data[('Biorefinery', 'Excess electricity [MW]')]
electricity_data_humbird_normalized = electricity_data * (100/humbird_electricity)
bx_electricity = plt.bar(positions_electricity, electricity_data_humbird_normalized, 0.5,
                         align='center',
                         color=orange.RGBn,
                         edgecolor=orange_shade.RGBn)

# plot_vertical_line(7.5, color=colors.orange_tint.shade(15).RGBn, ls='-.')
# plot_vertical_line(9.5, color=colors.orange_tint.shade(15).RGBn, ls='-.')

# %% Plot installed cost

units = '10^6 USD'
plot_vertical_line(8.5, color=colors.grey_tint.RGBn)
positions_installed = tuple(range(9, 16))
installed_cols = [(i, 'Installed equipment cost [10^6 USD]') for i in areas[1:]]
humbird_installed = np.array([24.2, 32.9, 31.2, 22.3, 49.4, 5, 66, 6.9])
installed_data = data[installed_cols]
# installed_data[('Biorefinery', 'Installed cost')] = installed_data.sum(1)
installed_data_humbird_normalized = installed_data * (100/humbird_installed[1:])
bx_installed = plt.bar(
    positions_installed, installed_data_humbird_normalized, 0.5,
    align='center',
    color=red.RGBn,
    edgecolor=red_shade.RGBn,
)
plot_horizontal_line(100, ls='--')
u_lb = 0; y_ub = 200
plt.ylim(0, 200)
yticks = np.arange(0, 201, 50)
plt.yticks(yticks)
y_text = 0.885*y_ub
y_letter = 0.875 * y_ub
plt.text(4, y_text, "Electricity demand", color=orange_shade.RGBn,
          horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(-0.25, y_letter, "C", color=black.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(12, y_text, "Installed cost", color=red_shade.RGBn,
          horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(8.75, y_letter, "D", color=black.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(15.75, y_letter, "E", color=black.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
# plt.text(9.25, y_letter, "D", color=black.RGBn,
#          horizontalalignment='center', fontsize=14, fontweight='bold')

# plt.text(16, y_text, "Ethanol\nsales", color=colors.red_shade.RGBn,
#           horizontalalignment='center', fontsize=12, fontweight='bold')
# plt.text(16.0, y_text, "MESP", color=colors.blue_shade.RGBn,
#           horizontalalignment='center', fontsize=12, fontweight='bold')

plt.xlim(-0.5, 18.5)
plt.ylabel("Metric over benchmark [%]")
area_marks = [i.replace(' ', '\n') for i in areas]
xmarks = area_marks + ['Excess'] + area_marks[1:] + ['Steam\ndemand', '   EtOH\n    prod.', '    MESP']
xticks = positions_electricity + positions_installed + positions_other
plt.xticks(xticks, xmarks)
metric_over_benchmark_ax.set_zorder(1e6)

plt.sca(magnitude_ax)
plt.fill_between([-0.5, 15.5], 0, 1, color=colors.neutral_tint.tint(85).RGBn)

electricity_areas = humbird_electricity # electricity_data.median()
electricity_areas /= max(electricity_areas)
plt.bar(positions_electricity, electricity_areas, 0.5,
        align='center', label="Electricity demand",
        color=orange_shade.RGBn,
        edgecolor=orange_shade.tint(10).shade(15).RGBn)
plot_vertical_line(8.5, color=colors.grey_tint.RGBn)

installed_areas = humbird_installed[1:]# installed_data.median()
installed_areas /= max(installed_areas)
plt.bar(positions_installed, installed_areas, 0.5,
        align='center', label="Installed equipment cost",
        color=red_shade.RGBn,
        edgecolor=red_shade.tint(10).shade(15).RGBn)

plot_vertical_line(15.5, color=colors.grey_tint.RGBn)
# plot_vertical_line(15.5, color='k', lw=0.8)
# plt.bar(positions_MESP, [1], 0.5,
#         align='center', label="MESP",
#         color=colors.blue.tint(30).shade(15).RGBn,
#         edgecolor=colors.blue_shade.RGBn)

plot_vertical_line(-0.5, color='k')
plt.hlines([1], [-0.5], [15.5], color='k')
magnitude_ax.spines['top'].set_visible(False)
magnitude_ax.spines['right'].set_visible(False)
magnitude_ax.tick_params(axis="x", direction="in", length=0)
magnitude_ax.set_zorder(2)
metric_over_benchmark_ax.set_zorder(1)
plt.yticks([], [])
plt.xticks(xticks[:-3], [])
plt.ylim(0, 1)
plt.xlim(-0.5, 18.5)
plt.text(4, 0.65, "Benchmark magnitude", color=black.RGBn,
         horizontalalignment='center', fontsize=12, fontweight='bold')
plt.text(-0.25, 0.60, "A", color=black.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.text(8.75, 0.60, "B", color=black.RGBn,
         horizontalalignment='left', fontsize=16, fontweight='bold')
plt.subplots_adjust(hspace=.0)

metric_over_benchmark_ax.tick_params(axis='x', direction="inout", length=4)
for ax in axes:
    ax.tick_params(axis='y', right=False, direction="inout", length=4)
ax2 = metric_over_benchmark_ax.twinx()
plt.sca(ax2)
plt.ylim(0, 200)
plt.yticks(yticks, [])
ax2.zorder = 1000
ax2.tick_params(direction="in")
    
xlabels = metric_over_benchmark_ax.get_xticklabels()    

plt.subplots_adjust(
    hspace=0, wspace=0.05,
    top=0.95, bottom=0.15,
    left=0.08, right=0.98,
)

images_folder = os.path.join(os.path.dirname(__file__), 'images')
for i in ('svg', 'png'):
    file = os.path.join(images_folder, f'single_point_evaluation.{i}')
    plt.savefig(file, transparent=True)
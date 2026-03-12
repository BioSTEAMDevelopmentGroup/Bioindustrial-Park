#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multi-panel contour plots across kinetic parameter pairs and strategies.

This script is adapted from:
- across_kinetic_params.py for file discovery, metric metadata, and spec grids
- HP_all_TRY_plot.py for the multi-panel subplot layout pattern

Rows correspond to (x_label, y_label) pairs.
Columns correspond to strategy names.
Each panel loads the matching CSV and plots a contour map for one metric.

Example layout:
    row_parameter_pairs = [
        ('k_1e', 'k_1ie'),
        ('k_1e', 'k_7ie'),
        ('k_13', 'k_7ii'),
    ]
    strategies = [
        'fixed_batch',
        'adaptive_batch',
        'fixed_fed-batch',
        'adaptive_fed-batch',
    ]
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

import contourplots
from biosteam.utils import colors
from biorefineries import isobutanol

# -----------------------------------------------------------------------------
# User inputs
# -----------------------------------------------------------------------------
metric = 'MPSP'
z_label = 'Spike feed glucose concentration'
steps = (25, 25, 1)

row_parameter_pairs = [
    ('k_1e', 'k_1ie'),
    ('k_1e', 'k_7ie'),
    ('k_13', 'k_7ii'),
]

strategies = [
    'fixed_batch',
    'adaptive_batch',
    'fixed_fed-batch',
    # 'adaptive_fed-batch',
]

output_filename = 'MPSP_multi_panel_kinetics.png'

# Optional: set to None to auto-compute from all loaded panels.
manual_w_levels = np.arange(0.2, 2.401, 0.05) if metric == 'MPSP' else None
manual_cbar_ticks = np.arange(0.2, 2.401, 0.2) if metric == 'MPSP' else None
manual_w_ticks = [2.4] if metric == 'MPSP' else None

# -----------------------------------------------------------------------------
# Shared metadata from across_kinetic_params.py
# -----------------------------------------------------------------------------
get_rounded_str = contourplots.utils.get_rounded_str
fbs_spec = isobutanol.models.fbs_spec

metrics_units = {
    'MPSP': r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$",
    'AOC': r'MM\$/y',
    'TCI': r'MM\$',
    'Combined Yield': 'g-EtOH-and-IBO/g-sugars',
    'EtOH Titer': 'g-EtOH/L-broth',
    'EtOH Productivity': 'g-EtOH/L-broth/h',
    'Number of glucose spikes': '',
    'Fermentation time': 'h',
    'Total Q sugar evap': 'kJ/h',
    'Target sugars concentration': 'g-sugars/L-broth',
    'Cell loading': 'g-cell/L-broth',
    'Active cell loading': 'g-cell/L-broth',
    'EtOH Yield': 'g-EtOH/g-sugars',
    'IBO Yield': 'g-IBO/g-sugars',
    'IBO Titer': 'g-IBO/L-broth',
    'IBO Productivity': 'g-IBO/L-broth/h',
    'Actual aeration required': 'kmol-O2/h',
}

metrics_plot_names = {k: k for k in metrics_units.keys()}
metrics_plot_names['Total Q sugar evap'] = 'Slurry evaporation duty'
metrics_plot_names['Actual aeration required'] = 'Aeration required'
for k in list(metrics_plot_names.keys()):
    metrics_plot_names[k] = (
        metrics_plot_names[k]
        .replace('IBO', 'Isobutanol')
        .replace('EtOH', 'Ethanol')
    )
metrics_plot_names['Cell loading'] = 'Cell density'

# Preserve the original special casing from across_kinetic_params.py:
# if x_label == k_1e, Combined Yield becomes Ethanol Yield.
# We handle that panel-by-panel during plotting.

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------
def JBEI_UCB_colormap(N_levels=90, reverse=False):
    jbei_orange = (233 / 255, 83 / 255, 39 / 255)
    ucb_blue = (0 / 255, 38 / 255, 118 / 255)
    ucb_yellow = (253 / 255, 181 / 255, 21 / 255)
    cmap_colors = [ucb_yellow, jbei_orange, ucb_blue, colors.grey_dark.RGBn]
    if reverse:
        cmap_colors.reverse()
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)


def get_metric_plot_name(metric_name, x_label):
    name = metrics_plot_names[metric_name]
    if x_label == 'k_1e' and metric_name == 'Combined Yield':
        name = 'Ethanol Yield'
    if metric_name not in ('TCI', 'AOC', 'MPSP'):
        name = name.lower()
    return name


def format_param_label(label):
    return r"$\bf{" + label[:label.index('_')] + '}_{' + label[label.index('_') + 1:] + '}"$'


def get_specs_and_ticks(x_label, y_label, z_label, steps):
    spec_1 = spec_2 = spec_3 = None
    x_units = y_units = z_units = None
    x_ticks = y_ticks = z_ticks = None

    if x_label in ('k_1e',):
        spec_1 = np.linspace(1.0, 300.0, steps[0])
        x_units = r"$\mathrm{g} \cdot \mathrm{g}^{-1} \cdot \mathrm{h}^{-1}$"
        x_ticks = np.array([0, 100, 200, 300])
    elif x_label in ('k_13',):
        spec_1 = np.linspace(0.0, 40.0, steps[0])
        x_units = r"$\mathrm{g} \cdot \mathrm{g}^{-1} \cdot \mathrm{h}^{-1}$"
        x_ticks = np.array([0, 10, 20, 30, 40])
    else:
        raise ValueError(f'Unsupported x_label: {x_label}')

    if y_label in ('k_1ie', 'k_1ii', 'k_7ie', 'k_7ii'):
        spec_2 = np.linspace(0.0001, 0.5, steps[1])
        y_units = r"$\mathrm{L} \cdot \mathrm{g}^{-1}$"
        y_ticks = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    else:
        raise ValueError(f'Unsupported y_label: {y_label}')

    if z_label == 'Spike feed glucose concentration':
        spec_3 = np.array([fbs_spec.conc_sugars_feed_spike])
        z_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
        z_ticks = [0, 200, 400, 600, 800]
    else:
        raise ValueError(f'Unsupported z_label: {z_label}')

    return spec_1, spec_2, spec_3, x_units, y_units, z_units, x_ticks, y_ticks, z_ticks


def build_subfolder_name(x_label, y_label, strategy):
    return f'Kinetics_{x_label}_{y_label}\\{strategy}\\'


def build_file_prefix(x_label, y_label, z_label, strategy, steps):
    subfolder_name = build_subfolder_name(x_label, y_label, strategy)
    perform_feeding_strategy_opt = 'adaptive' in subfolder_name
    max_n = 0 if 'fed-batch' not in subfolder_name else 200
    file_prefix = (
        f"ibo_{steps}_{x_label[:5]}_{y_label[:5]}_{z_label[:5]}"
        f"_opt={perform_feeding_strategy_opt}_max_n={max_n}_"
    )
    return subfolder_name, file_prefix


def load_metric_array(x_label, y_label, strategy, z_label, metric, steps):
    subfolder_name, file_prefix = build_file_prefix(x_label, y_label, z_label, strategy, steps)
    folder = isobutanol_results_pub_filepath + subfolder_name
    filepath = os.path.join(folder, file_prefix + f'_{metric}.csv')

    if not os.path.exists(filepath):
        return None, filepath

    df = pd.read_csv(filepath)
    arr_2d = df[df.columns[1:]].to_numpy(dtype=float, na_value=np.nan)
    arr_3d = np.array([arr_2d])
    return arr_3d, filepath


def choose_metric_colormap(metric_name):
    lccm = metric_name.lower()
    if any(key in lccm for key in ('yield', 'titer', 'productivity', 'loading')):
        cmap = JBEI_UCB_colormap(reverse=True)
        cmap_over_color = colors.yellow_tint.RGBn
        extend_cmap = 'neither'
    else:
        cmap = JBEI_UCB_colormap(reverse=False)
        cmap_over_color = colors.grey_dark.shade(8).RGBn
        extend_cmap = 'max'
    return cmap, cmap_over_color, extend_cmap


def compute_levels_from_arrays(arrays, metric_name):
    valid = []
    for arr in arrays:
        if arr is None:
            continue
        vals = arr[np.where(~np.isnan(arr))]
        if vals.size:
            valid.append(vals)

    if not valid:
        raise ValueError(f'No valid data found for metric {metric_name!r}.')

    values = np.concatenate(valid)

    if manual_w_levels is not None and manual_cbar_ticks is not None:
        return manual_w_levels, manual_cbar_ticks, manual_w_ticks

    lccm = metric_name.lower()
    min_val = float(np.nanmin(values))
    max_val = float(np.nanmax(values))

    if any(key in lccm for key in ('yield', 'titer', 'productivity', 'cell loading')):
        top_val = max_val
    else:
        top_val = max_val * 0.8

    if np.isclose(top_val, min_val):
        top_val = min_val + 1e-6

    levels = np.arange(min_val, top_val * 1.00000000001, (top_val - min_val) / 80)
    cbar_ticks = np.arange(min_val, top_val * 1.00000000001, (top_val - min_val) / 5)
    w_ticks = sorted(set([
        float(np.percentile(values, 25)),
        float(np.percentile(values, 50)),
        float(np.percentile(values, 75)),
        top_val,
    ]))
    return levels, cbar_ticks, w_ticks


# -----------------------------------------------------------------------------
# Load all requested data first so the whole figure shares one color scale.
# -----------------------------------------------------------------------------
panel_data = {}
all_arrays = []
missing_files = []

for x_label, y_label in row_parameter_pairs:
    panel_data[(x_label, y_label)] = {}
    for strategy in strategies:
        arr, filepath = load_metric_array(x_label, y_label, strategy, z_label, metric, steps)
        panel_data[(x_label, y_label)][strategy] = arr
        if arr is None:
            missing_files.append(filepath)
        else:
            all_arrays.append(arr)

if not all_arrays:
    raise FileNotFoundError(
        'No matching metric CSV files were found for the requested row/column layout.'
    )

w_levels, cbar_ticks, w_ticks = compute_levels_from_arrays(all_arrays, metric)
cmap, cmap_over_color, extend_cmap = choose_metric_colormap(metric)

# -----------------------------------------------------------------------------
# Figure layout, based on HP_all_TRY_plot.py style.
# -----------------------------------------------------------------------------
nrows = len(row_parameter_pairs)
ncols = len(strategies)
fig, axs = plt.subplots(nrows, ncols, constrained_layout=True)

# Normalize axs indexing for 1-row or 1-col cases
if nrows == 1 and ncols == 1:
    axs = np.array([[axs]])
elif nrows == 1:
    axs = np.array([axs])
elif ncols == 1:
    axs = np.array([[ax] for ax in axs])

axis_title_fonts = {'size': {'x': 11, 'y': 11, 'z': 11, 'w': 11}}
default_fontsize = 11.0
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5

for i, (x_label, y_label) in enumerate(row_parameter_pairs):
    (
        spec_1,
        spec_2,
        spec_3,
        x_units,
        y_units,
        z_units,
        x_ticks,
        y_ticks,
        z_ticks,
    ) = get_specs_and_ticks(x_label, y_label, z_label, steps)

    x_label_for_plot = format_param_label(x_label)
    y_label_for_plot = format_param_label(y_label)
    w_label = r"$\bf" + get_metric_plot_name(metric, x_label).replace(' ', '\\ ') + "$"

    for j, strategy in enumerate(strategies):
        ax = axs[i, j]
        arr = panel_data[(x_label, y_label)][strategy]

        if arr is None:
            ax.set_axis_off()
            ax.text(
                0.5,
                0.5,
                'File not found',
                transform=ax.transAxes,
                ha='center',
                va='center',
                fontsize=10,
            )
            continue

        show_x_ticklabels = i == nrows - 1
        show_y_ticklabels = j == 0

        contourplots.animated_contourplot(
            w_data_vs_x_y_at_multiple_z=arr,
            x_data=spec_1,
            y_data=spec_2,
            z_data=spec_3,
            x_label=x_label_for_plot,
            y_label=y_label_for_plot,
            z_label=r"$\bf" + z_label + "$",
            w_label=w_label,
            x_ticks=x_ticks,
            y_ticks=y_ticks,
            z_ticks=z_ticks,
            w_levels=w_levels,
            w_ticks=w_ticks,
            x_units=x_units,
            y_units=y_units,
            z_units=z_units,
            w_units=metrics_units[metric],
            fmt_clabel=lambda cvalue: get_rounded_str(cvalue, 3),
            cmap=cmap,
            cmap_over_color=cmap_over_color,
            extend_cmap=extend_cmap,
            cbar_ticks=cbar_ticks,
            z_marker_color='g',
            fps=1,
            n_loops='inf',
            animated_contourplot_filename='ignore_this_filename',
            keep_frames=False,
            keep_gifs=False,
            axis_title_fonts=axis_title_fonts,
            clabel_fontsize=clabel_fontsize,
            default_fontsize=default_fontsize,
            axis_tick_fontsize=axis_tick_fontsize,
            n_minor_ticks=3,
            cbar_n_minor_ticks=3,
            units_on_newline=(False, False, False, False),
            units_opening_brackets=[' ['] * 4,
            units_closing_brackets=[']'] * 4,
            round_xticks_to=0,
            round_yticks_to=1,
            include_top_bar=False,
            include_cbar=False,
            include_axis_labels=False,
            include_x_axis_ticklabels=show_x_ticklabels,
            include_y_axis_ticklabels=show_y_ticklabels,
            fig_ax_to_use=(fig, ax),
        )

        if i == 0:
            ax.set_title(strategy.replace('_', '\n'), fontsize=11, fontweight='bold')
        if j == 0:
            ax.set_ylabel(f'{x_label} vs\n{y_label}', fontsize=11, fontweight='bold')

# -----------------------------------------------------------------------------
# Shared labels and shared colorbar
# -----------------------------------------------------------------------------
fig.supxlabel('Kinetic parameter on x-axis', fontsize=12, fontweight='bold')
fig.supylabel('Kinetic parameter on y-axis', fontsize=12, fontweight='bold')

norm = Normalize(vmin=float(w_levels[0]), vmax=float(w_levels[-1]))
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs.ravel().tolist(), shrink=0.95, pad=0.02)
cbar.set_label(f"{get_metric_plot_name(metric, row_parameter_pairs[0][0]).title()} [{metrics_units[metric]}]", fontsize=11)
cbar.set_ticks(cbar_ticks)

fig.set_figwidth(3.2 * ncols + 1.2)
fig.set_figheight(2.7 * nrows + 0.8)

plt.savefig(
    output_filename,
    transparent=False,
    facecolor='white',
    bbox_inches='tight',
    dpi=300,
)
plt.close(fig)

print(f'Saved multi-panel figure to: {output_filename}')
if missing_files:
    print('\nMissing files:')
    for fp in missing_files:
        print(f'  - {fp}')

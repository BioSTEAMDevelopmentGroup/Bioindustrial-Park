#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2025-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.



import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

import contourplots
from biosteam.utils import colors
from biorefineries import isobutanol

#%%
# -----------------------------------------------------------------------------
# User inputs
# -----------------------------------------------------------------------------
metric = 'MPSP'
x_label = "Threshold glucose concentration" # title of the x axis
x_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
x_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

y_label = "Target glucose concentration" # title of the y axis
y_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
y_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

z_label = "Max. no. of glucose spikes" # title of the z axis
z_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
z_ticks = [0, 5, 10, 15, 20]

steps = (25, 25, 5)

spec_1 = threshold_conc_sugarses = np.linspace(1., 400., steps[0])

spec_2 = target_conc_sugarses = np.linspace(10., 400., steps[1])

spec_3 = row_max_ns = max_n_glu_spikes = np.linspace(0, 20, steps[2])

stage_1_times = [10., 15., 20., 25.]

output_filename = 'MPSP_multi_panel_feed_strat.png'

# Optional: set to None to auto-compute from all loaded panels.
manual_w_levels = np.arange(0.65, 1.2, 0.01) if metric == 'MPSP' else None
manual_cbar_ticks = np.arange(0.65, 1.2, 0.05) if metric == 'MPSP' else None
# manual_w_ticks = [2.4] if metric == 'MPSP' else None
manual_w_ticks = [0.8, 1.2]

comparison_range = []
# if metric == 'MPSP':
#     comparison_range = EtOH_market_range = np.array([
#         0.52, # 1.5475 $/gal/(3.7854 L/gal * 0.789 kg/L)
#         1.15, # 3.4500 $/gal/(3.7854 L/gal * 0.789 kg/L)
#         ]) # Jan 2021 - Dec 2025 5-year low and high from https://tradingeconomics.com/commodity/ethanol

get_rounded_str = contourplots.utils.get_rounded_str
fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3)
# if metric == 'MPSP':
#     fmt_clabel = lambda cvalue: rf"$\${cvalue:.2f}\cdot\mathrm{{kg}}^{{-1}}$"
    
# -----------------------------------------------------------------------------
# Shared metadata from across_kinetic_params.py
# -----------------------------------------------------------------------------
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


def get_metrics_to_opt_for_x_label(x_label):
    return [
        'Cell loading',
        'EtOH Titer',
        'EtOH Productivity',
        'Combined Yield',
        'Total Q sugar evap',
        'Actual aeration required',
        'TCI',
        'AOC',
        'MPSP',
        ]


def metric_is_minimized(metric_name):
    return metric_name in (
        'MPSP', 'AOC', 'TCI', 'Total Q sugar evap',
        'Fermentation time', 'Actual aeration required'
    )


def format_param_label(label):
    return r"$\bf{label.replace(' ', '\ ')}$"


def format_axis_label_with_units(label, units):
    base = format_param_label(label)
    if units:
        return f"{base} [{units}]"
    return base



def load_metric_array(x_label, y_label, s1t, max_n, metric_name, steps):
    subfolder_name = 'Feed-strat\\'
    file_prefix = f'ibo_{steps}_{x_label[:5]}_{y_label[:5]}_Spike_s1t={s1t}_max_n={max_n}_'
    folder = isobutanol_results_pub_filepath + subfolder_name
    filepath = os.path.join(folder, file_prefix + f'_{metric_name}.csv')
    if not os.path.exists(filepath):
        return None, filepath
    df = pd.read_csv(filepath)
    arr_2d = df[df.columns[1:]].to_numpy(dtype=float, na_value=np.nan)
    return np.array([arr_2d]), filepath


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


def build_metric_marker_styles(metrics_to_opt, row_max_ns):
    opt_marker_shapes = ['o', '^', 's', 'p', 'v', '<', '>', 'h']
    if x_label == 'Threshold glucose concentration':
        for shape in ['v', '<', '>']:
            if shape in opt_marker_shapes:
                opt_marker_shapes.remove(shape)

    styles = {}
    legend_labels = {}
    generic_i = 0
    for metric_name in metrics_to_opt:
        if metric_name == 'MPSP':
            style = ('*', '#33ccff', 16)
        elif metric_name == 'TCI':
            style = ('s', '#33ccff', 12)
        elif metric_name == 'AOC':
            style = ('p', '#33ccff', 12)
        elif metric_name == 'Total Q sugar evap':
            style = ('P', 'w', 12)
        elif metric_name == 'Actual aeration required':
            style = ('X', 'w', 12)
        else:
            style = (opt_marker_shapes[generic_i], 'w', 12)
            generic_i += 1
        styles[metric_name] = style

        label = metrics_plot_names[metric_name]
        if metric_name == 'Combined Yield' and x_label == 'Threshold glucose concentration':
            label = 'Ethanol Yield'
        if label not in ('TCI', 'AOC', 'MPSP'):
            label = label.lower()
        legend_labels[metric_name] = label
    return styles, legend_labels


def compute_global_optima(metrics_to_opt, row_max_ns, stage_1_times, z_label, steps):
    global_optima = {}
    missing_opt_metric_files = []
    for opt_metric in metrics_to_opt:
        best_value = None
        best_record = None
        for row_i, max_n in enumerate(row_max_ns):
            for s1t in stage_1_times:
                arr, filepath = load_metric_array(x_label, y_label, s1t, max_n, opt_metric, steps)
                if arr is None:
                    missing_opt_metric_files.append(filepath)
                    continue
                arr2d = arr[0]
                finite_mask = np.isfinite(arr2d)
                if not finite_mask.any():
                    continue
                candidate_value = float(np.nanmin(arr2d)) if metric_is_minimized(opt_metric) else float(np.nanmax(arr2d))
                is_better = best_value is None
                if not is_better:
                    is_better = candidate_value < best_value if metric_is_minimized(opt_metric) else candidate_value > best_value
                if is_better:
                    idx = np.where(arr2d == candidate_value)
                    row_idx = int(idx[0][0])
                    col_idx = int(idx[1][0])
                    best_value = candidate_value
                    best_record = {
                        'metric': opt_metric,
                        'row_index': row_i,
                        'x_label': x_label,
                        'y_label': y_label,
                        'max_n': max_n,
                        's1t': s1t,
                        'x': float(spec_1[col_idx]),
                        'y': float(spec_2[row_idx]),
                        'value': candidate_value,
                    }
        if best_record is not None:
            global_optima[opt_metric] = best_record
    return global_optima, sorted(set(missing_opt_metric_files))

        
def get_zeroth_spec_2_val_for_condition_for_all_spec_1_vals(metric_arr, condition):
    s2s = []
    for s1i in range(len(metric_arr[0])):
        success = False
        for s2i in range(len(metric_arr)):
            if condition(metric_arr[s2i][s1i]):
                s2s.append(spec_2[s2i])
                success = True
                break
        if not success:
            s2s.append(np.nan)
    return s2s


round_off = contourplots.utils.round_off

def get_optima_metric_values(global_optima, metrics_to_opt, steps):
    opt_coords_metric_values = {}
    missing_metric_value_files = []
    for metric_name in metrics_to_opt:
        if metric_name not in global_optima:
            continue

        opt_record = global_optima[metric_name]
        opt_x = float(opt_record['x'])
        opt_y = float(opt_record['y'])
        opt_s1t = opt_record['s1t']
        opt_max_n = opt_record['max_n']

        x_matches = np.where(np.isclose(spec_1, opt_x))[0]
        y_matches = np.where(np.isclose(spec_2, opt_y))[0]
        if not len(x_matches) or not len(y_matches):
            continue
        opt_x_ind = int(x_matches[0])
        opt_y_ind = int(y_matches[0])

        coords = (opt_x, opt_y, opt_s1t, opt_max_n)
        opt_coords_metric_values[metric_name] = {'Coords': coords}
        for comparison_metric in metrics_to_opt:
            arr, filepath = load_metric_array(x_label, y_label, opt_s1t, opt_max_n, comparison_metric, steps)
            if arr is None:
                missing_metric_value_files.append(filepath)
                opt_coords_metric_values[metric_name][comparison_metric] = np.nan
                continue
            opt_coords_metric_values[metric_name][comparison_metric] = arr[0][opt_y_ind, opt_x_ind]
    return opt_coords_metric_values, sorted(set(missing_metric_value_files))


def get_optima_comparisons(opt_coords_metric_values, rel_to_m='MPSP'):
    if rel_to_m not in opt_coords_metric_values:
        print(f"\nNo optima data available for '{rel_to_m}'.")
        return

    print(f"\n\nRelative to the optimum for '{rel_to_m}', at the optimum for:")
    print('\n')
    i = 0
    for m1, v in opt_coords_metric_values.items():
        i += 1
        print(f"{i}. '{m1}' ({v['Coords']}),")
        for m2 in v.keys():
            if m2 == 'Coords':
                continue

            curr_val = v[m2]
            ref_val = opt_coords_metric_values[rel_to_m].get(m2, np.nan)
            if np.isnan(curr_val):
                print(f"'{m2}' is unavailable.")
                continue
            if np.isnan(ref_val):
                print(f"'{m2}' is {round_off(curr_val, 3)}.")
                continue
            if np.isclose(ref_val, 0.0):
                if np.isclose(curr_val, 0.0):
                    print(f"'{m2}' is zero.")
                else:
                    print(f"'{m2}' is {round_off(curr_val, 3)}.")
                continue

            rel_diff = curr_val / ref_val - 1
            sign = '+' if rel_diff > 0 else '-'
            print(f"'{m2}' is {round_off(curr_val, 3)}, which is {sign} {abs(int(100 * rel_diff))}%.")
        print('\n')

    print('\n Note that optima overlap for:')
    print('\n')
    coords_metrics = {}
    for metric, v in opt_coords_metric_values.items():
        coords_metrics.setdefault(v['Coords'], []).append(metric)
    any_overlaps = False
    for coords, metrics in coords_metrics.items():
        if len(metrics) > 1:
            any_overlaps = True
            print(f'{metrics}: optimum coordinates are {coords}.')
    if not any_overlaps:
        print('No overlapping optimum coordinates were found.')

# -----------------------------------------------------------------------------
# Determine which optimization metrics to show.
# -----------------------------------------------------------------------------
# metrics_to_opt = []
# for x_label, _ in row_max_ns:
#     for m in get_metrics_to_opt_for_x_label(x_label):
#         if m not in metrics_to_opt:
#             metrics_to_opt.append(m)
metrics_to_opt = [
    'Cell loading',
    'EtOH Titer',
    'EtOH Productivity',
    'Combined Yield',
    'Total Q sugar evap',
    'Actual aeration required',
    'TCI',
    'AOC',
    'MPSP',
    ]

metric_marker_styles, _metric_legend_labels = build_metric_marker_styles(metrics_to_opt, row_max_ns)

global_optima, missing_opt_metric_files = compute_global_optima(
    metrics_to_opt, row_max_ns, stage_1_times, z_label, steps
)

opt_coords_metric_values, missing_opt_value_files = get_optima_metric_values(
    global_optima, metrics_to_opt, steps
)

# -----------------------------------------------------------------------------
# Load panel data for the plotted metric.
# -----------------------------------------------------------------------------
panel_data = {}
all_arrays = []
missing_files = []
panel_data[(x_label, y_label)] = {}
for max_n in row_max_ns:
    for s1t in stage_1_times:
        arr, filepath = load_metric_array(x_label, y_label, s1t, max_n, metric, steps)
        panel_data[(x_label, y_label)][(s1t, max_n)] = arr
        if arr is None:
            missing_files.append(filepath)
        else:
            all_arrays.append(arr)

if not all_arrays:
    raise FileNotFoundError('No matching metric CSV files were found for the requested row/column layout.')

w_levels, cbar_ticks, w_ticks = compute_levels_from_arrays(all_arrays, metric)
cmap, cmap_over_color, extend_cmap = choose_metric_colormap(metric)

# -----------------------------------------------------------------------------
# Figure layout
# -----------------------------------------------------------------------------
nrows = len(row_max_ns)
ncols = len(stage_1_times)
fig, axs = plt.subplots(nrows, ncols, constrained_layout=True)
if nrows == 1 and ncols == 1:
    axs = np.array([[axs]])
elif nrows == 1:
    axs = np.array([axs])
elif ncols == 1:
    axs = np.array([[ax] for ax in axs])

axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 15.
clabel_fontsize = 12
axis_tick_fontsize = 15

for i, max_n in enumerate(row_max_ns):
    # spec_1, spec_2, spec_3, x_units, y_units, z_units, x_ticks, y_ticks, z_ticks = get_specs_and_ticks(
    #     x_label, y_label, z_label, steps
    # )
    x_label_for_plot = format_param_label(x_label)
    y_label_for_plot = format_param_label(y_label)
    w_label = r"$\bf" + get_metric_plot_name(metric, x_label).replace(' ', '\\ ') + "$"

    for j, s1t in enumerate(stage_1_times):
        ax = axs[i, j]
        arr = panel_data[(x_label, y_label)][(s1t, max_n)]

        if arr is None:
            ax.set_axis_off()
            ax.text(0.5, 0.5, 'File not found', transform=ax.transAxes, ha='center', va='center', fontsize=10)
            continue

        additional_points = {}
        for opt_metric, opt_record in global_optima.items():
            if opt_record['row_index'] == i and opt_record['s1t'] == s1t:
                additional_points[(opt_record['x'], opt_record['y'])] = metric_marker_styles[opt_metric]
        
        add_lines = {}
        if not arr is None and not max_n==0.0:
            arr_n_spikes, filepath = load_metric_array(x_label, y_label, s1t, max_n, 'Number of glucose spikes', steps)
        
            line_first_app_n_glu_spikes_0 = tuple(get_zeroth_spec_2_val_for_condition_for_all_spec_1_vals(arr_n_spikes[0],
                                                                              condition = lambda i: i==0))
            add_lines = {line_first_app_n_glu_spikes_0: {'color': 'white', 'linewidth': 1.0, 'alpha': 1.0}}
            
        contourplots.animated_contourplot(
            w_data_vs_x_y_at_multiple_z=arr,
            x_data=spec_1,
            y_data=spec_2,
            # z_data=spec_3,
            z_data=[max_n,],
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
            fmt_clabel=fmt_clabel,
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
            include_x_axis_ticklabels=(i == nrows - 1),
            # include_x_axis_ticklabels=True,
            include_last_x_axis_ticklabel=(j==ncols-1),
            # include_last_x_axis_ticklabel=True,
            include_y_axis_ticklabels=(j == 0),
            include_last_y_axis_ticklabel=(i==0),
            additional_points=additional_points,
            fig_ax_to_use=(fig, ax),
            comparison_range=comparison_range,
            inline_spacing=0.,
            label_over_color='black',
            add_lines=add_lines,
        )

        # if i == 0:
        #     ax.set_title(str(s1t).replace('_', '\n'), fontsize=11, fontweight='bold')

# -----------------------------------------------------------------------------
# Shared colorbar and layout
# -----------------------------------------------------------------------------
plt.subplots_adjust(wspace=0., hspace=0.)

norm = Normalize(vmin=float(w_levels[0]), vmax=float(w_levels[-1]))
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs.ravel().tolist(), shrink=0.95, pad=0.02)
cbar.set_label(f"{get_metric_plot_name(metric, x_label).title()} [{metrics_units[metric]}]", fontsize=11)
cbar.set_ticks(cbar_ticks)

fig.set_figwidth(3.2 * ncols + 1.2)
fig.set_figheight(2.7 * nrows + 1.6)

plt.savefig(output_filename, transparent=False, facecolor='white', bbox_inches='tight', dpi=600)
plt.close(fig)

print(f'Saved multi-panel figure to: {output_filename}')
if global_optima:
    print('\nGlobal optima used for markers:')
    for metric_name in metrics_to_opt:
        if metric_name not in global_optima:
            continue
        opt_record = global_optima[metric_name]
        print(
            '  {metric}: value={value}, x={x}, y={y}, panel=({xl}, {yl}, max_n={mn}, s1t={st})'.format(
                metric=metric_name,
                value=get_rounded_str(opt_record['value'], 6),
                x=get_rounded_str(opt_record['x'], 6),
                y=get_rounded_str(opt_record['y'], 6),
                xl=opt_record['x_label'],
                yl=opt_record['y_label'],
                mn=get_rounded_str(opt_record['max_n'], 6),
                st=get_rounded_str(opt_record['s1t'], 6),
            )
        )
if missing_files:
    print('\nMissing plotted-metric files:')
    for fp in sorted(set(missing_files)):
        print(f'  - {fp}')
if missing_opt_metric_files:
    print('\nMissing optimization-metric files:')
    for fp in missing_opt_metric_files:
        print(f'  - {fp}')
if missing_opt_value_files:
    print('\nMissing comparison-value files:')
    for fp in missing_opt_value_files:
        print(f'  - {fp}')

get_optima_comparisons(opt_coords_metric_values, rel_to_m='MPSP')

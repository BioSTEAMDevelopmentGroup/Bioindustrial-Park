#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2025-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

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

#%%
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
    'adaptive_fed-batch',
]


#%% Global variables
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
    
#%%
def plot_save_kinetic_multipanel(metric, row_parameter_pairs=None, strategies=None, steps=None,
                                 z_label=None,):
    if z_label is None:
        z_label = 'Spike feed glucose concentration'
    if steps is None:
        steps = (25, 25, 1)
    if row_parameter_pairs is None:
        row_parameter_pairs = [
            ('k_1e', 'k_1ie'),
            ('k_1e', 'k_7ie'),
            ('k_13', 'k_7ii'),
        ]
    if strategies is None:
        strategies = [
            'fixed_batch',
            'adaptive_batch',
            'fixed_fed-batch',
            'adaptive_fed-batch',
        ]
        
    output_filename = f'{metric}_{row_parameter_pairs}.png'
    
    # Optional: set to None to auto-compute from all loaded panels.
    manual_w_levels = np.arange(0.2, 2.401, 0.05) if metric == 'MPSP' else None
    manual_cbar_ticks = np.arange(0.2, 2.401, 0.2) if metric == 'MPSP' else None
    manual_w_ticks = [2.4] if metric == 'MPSP' else None
    
    comparison_range = []
    if metric == 'MPSP':
        comparison_range = EtOH_market_range = np.array([
            0.52, # 1.5475 $/gal/(3.7854 L/gal * 0.789 kg/L)
            1.15, # 3.4500 $/gal/(3.7854 L/gal * 0.789 kg/L)
            ]) # Jan 2021 - Dec 2025 5-year low and high from https://tradingeconomics.com/commodity/ethanol
    
    get_rounded_str = contourplots.utils.get_rounded_str
    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3)
    # if metric == 'MPSP':
    #     fmt_clabel = lambda cvalue: rf"$\${cvalue:.2f}\cdot\mathrm{{kg}}^{{-1}}$"
        
    # -----------------------------------------------------------------------------
    # Shared metadata from across_kinetic_params.py
    # -----------------------------------------------------------------------------
    fbs_spec = isobutanol.models.fbs_spec
    
    
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
    
    os.chdir(isobutanol_results_pub_filepath)
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
        if x_label == 'k_1e':
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
        if x_label == 'k_13':
            return [
                'Cell loading',
                'EtOH Titer',
                'EtOH Productivity',
                'EtOH Yield',
                'IBO Titer',
                'IBO Productivity',
                'IBO Yield',
                'Combined Yield',
                'Total Q sugar evap',
                'Actual aeration required',
                'TCI',
                'AOC',
                'MPSP',
            ]
        return []
    
    
    def metric_is_minimized(metric_name):
        return metric_name in (
            'MPSP', 'AOC', 'TCI', 'Total Q sugar evap',
            'Fermentation time', 'Actual aeration required'
        )
    
    
    def format_param_label(label):
        return r"$\bf{" + label[:label.index('_')] + '}_{' + label[label.index('_') + 1:] + '}"$'
    
    
    def format_axis_label_with_units(label, units):
        base = format_param_label(label)
        if units:
            return f"{base} [{units}]"
        return base
    
    
    def get_specs_and_ticks(x_label, y_label, z_label, steps):
        spec_1 = spec_2 = spec_3 = None
        x_units = y_units = z_units = None
        x_ticks = y_ticks = z_ticks = None
    
        if x_label == 'k_1e':
            spec_1 = np.linspace(1.0, 300.0, steps[0])
            x_units = r"$\mathrm{g} \cdot \mathrm{g}^{-1} \cdot \mathrm{h}^{-1}$"
            x_ticks = np.array([0, 100, 200, 300])
        elif x_label == 'k_13':
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
        fed_batch_max_n = 21
        if 'fixed' in strategy:
            fed_batch_max_n = 3
        max_n = 0 if not 'fed-batch' in subfolder_name else fed_batch_max_n
        file_prefix = (
            f"ibo_{steps}_{x_label[:5]}_{y_label[:5]}_{z_label[:5]}"
            f"_opt={perform_feeding_strategy_opt}_max_n={max_n}_"
        )
        return subfolder_name, file_prefix
    
    
    def load_metric_array(x_label, y_label, strategy, z_label, metric_name, steps):
        subfolder_name, file_prefix = build_file_prefix(x_label, y_label, z_label, strategy, steps)
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
    
    
    def build_metric_marker_styles(metrics_to_opt, row_parameter_pairs):
        opt_marker_shapes = ['^', 's', 'p', 'v', '<', '>', 'h']
        if all(x_label == 'k_1e' for x_label, _ in row_parameter_pairs):
            for shape in ['v', '<', '>']:
                if shape in opt_marker_shapes:
                    opt_marker_shapes.remove(shape)
    
        styles = {}
        legend_labels = {}
        generic_i = 0
        for metric_name in metrics_to_opt:
            if metric_name == 'Cell loading':
                style = ('o', 'w', 8)
            elif metric_name == 'EtOH Titer':
                style = ('^', 'w', 8)
            elif metric_name == 'EtOH Productivity':
                style = ('s', 'w', 8)
            elif metric_name == 'EtOH Yield':
                style = ('p', 'w', 8)
            elif metric_name == 'IBO Titer':
                style = ('v', 'w', 8)
            elif metric_name == 'IBO Productivity':
                style = ('<', 'w', 8)
            elif metric_name == 'IBO Yield':
                style = ('>', 'w', 8)
            elif metric_name == 'Combined Yield':
                style = ('h', 'w', 8)
            elif metric_name == 'Total Q sugar evap':
                style = ('P', 'w', 8)
            elif metric_name == 'Actual aeration required':
                style = ('X', 'w', 8)
            elif metric_name == 'TCI':
                style = ('s', '#33ccff', 8)
            elif metric_name == 'AOC':
                style = ('p', '#33ccff', 8)
            elif metric_name == 'MPSP':
                style = ('*', '#33ccff', 12)
            else:
                style = (opt_marker_shapes[generic_i], 'w', 8)
                generic_i += 1
            styles[metric_name] = style
    
            label = metrics_plot_names[metric_name]
            if metric_name == 'Combined Yield' and any(xl == 'k_1e' for xl, _ in row_parameter_pairs):
                label = 'Ethanol Yield'
            if label not in ('TCI', 'AOC', 'MPSP'):
                label = label.lower()
            legend_labels[metric_name] = label
        return styles, legend_labels
    
    
    def compute_rowwise_optima(metrics_to_opt, row_parameter_pairs, strategies, z_label, steps):
        rowwise_optima = {}
        missing_opt_metric_files = []
        for row_i, (x_label, y_label) in enumerate(row_parameter_pairs):
            rowwise_optima[row_i] = {}
            spec_1, spec_2, *_ = get_specs_and_ticks(x_label, y_label, z_label, steps)
            for opt_metric in metrics_to_opt:
                best_value = None
                best_record = None
                found_any = False
                for strategy in strategies:
                    
                    # For data with x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie'), only ethanol is produced; however, etoh yield is saved as 'combined yield' rather than 'etoh yield'
                    # if opt_metric=='EtOH Yield' and x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie'): 
                    #     arr, filepath = load_metric_array(x_label, y_label, strategy, z_label, 'Combined Yield', steps)
                    # elif opt_metric=='Combined Yield' and x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie'):
                    #     pass
                    # else: # business as usual
                    arr, filepath = load_metric_array(x_label, y_label, strategy, z_label, opt_metric, steps)
                    if arr is None:
                        missing_opt_metric_files.append(filepath)
                        continue
                    arr2d = arr[0]
                    finite_mask = np.isfinite(arr2d)
                    if not finite_mask.any():
                        continue
                    found_any = True
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
                            'strategy': strategy,
                            'x': float(spec_1[col_idx]),
                            'y': float(spec_2[row_idx]),
                            'value': candidate_value,
                        }
                if found_any and best_record is not None:
                    rowwise_optima[row_i][opt_metric] = best_record
        return rowwise_optima, sorted(set(missing_opt_metric_files))
    
    
    round_off = contourplots.utils.round_off
    
    
    def get_metrics_at_coordinate_for_panel(x_label, y_label, strategy, z_label, steps, x_coord, y_coord, metrics_to_opt):
        spec_1, spec_2, *_ = get_specs_and_ticks(x_label, y_label, z_label, steps)
        x_idx = int(np.abs(spec_1 - x_coord).argmin())
        y_idx = int(np.abs(spec_2 - y_coord).argmin())
    
        values = {}
        missing_summary_files = []
        for metric_name in metrics_to_opt:
            arr, filepath = load_metric_array(x_label, y_label, strategy, z_label, metric_name, steps)
            if arr is None:
                missing_summary_files.append(filepath)
                values[metric_name] = np.nan
                continue
            try:
                values[metric_name] = float(arr[0][y_idx, x_idx])
            except Exception:
                values[metric_name] = np.nan
        return values, missing_summary_files
    
    
    def build_rowwise_opt_coords_metric_values(rowwise_optima, row_parameter_pairs, metrics_to_opt, z_label, steps):
        rowwise_opt_coords_metric_values = {}
        missing_summary_files = []
        for row_i, (x_label, y_label) in enumerate(row_parameter_pairs):
            row_summary = {}
            for metric_name, opt_record in rowwise_optima.get(row_i, {}).items():
                metric_values, missing_for_metric = get_metrics_at_coordinate_for_panel(
                    x_label=x_label,
                    y_label=y_label,
                    strategy=opt_record['strategy'],
                    z_label=z_label,
                    steps=steps,
                    x_coord=opt_record['x'],
                    y_coord=opt_record['y'],
                    metrics_to_opt=metrics_to_opt,
                )
                missing_summary_files.extend(missing_for_metric)
                row_summary[metric_name] = {
                    'Coords': (opt_record['x'], opt_record['y'], opt_record['strategy']),
                    **metric_values,
                }
            rowwise_opt_coords_metric_values[row_i] = row_summary
        return rowwise_opt_coords_metric_values, sorted(set(missing_summary_files))
    
    
    def get_optima_comparisons(opt_coords_metric_values, rel_to_m='MPSP'):
        if rel_to_m not in opt_coords_metric_values:
            print(f"Reference optimum '{rel_to_m}' is not available for this row.")
            return
    
        print(f"\n\nRelative to the optimum for '{rel_to_m}', at the optimum for:")
        print('\n')
        i = 0
        for m1, v in opt_coords_metric_values.items():
            i += 1
            print(f"{i}. '{m1}' ({v['Coords']}),")
            for m2 in v.keys():
                if m2 != 'Coords':
                    try:
                        ref_val = opt_coords_metric_values[rel_to_m][m2]
                        rel_diff = v[m2] / ref_val - 1
                        sign = '+' if rel_diff > 0 else '-'
                        sign_other = '+' if sign=='-' else '-'
                        print(f"'{m2}' is {round_off(v[m2],3)}, which is {sign} {abs(int(100*rel_diff))}% (so {rel_to_m} opt is {sign_other} {abs(int(100*(rel_diff/(1+rel_diff))))}%).")
                    except Exception as e:
                        if 'divide' in str(e).lower() or (isinstance(v[m2], float) and np.isnan(v[m2])):
                            print(f"'{m2}' is zero or unavailable.")
                        else:
                            print(f"'{m2}' is {v[m2]}.")
            print('\n')
    
        print('\n Note that optima overlap for:')
        print('\n')
        coords_metrics = {}
        for metric, v in opt_coords_metric_values.items():
            if v['Coords'] not in coords_metrics.keys():
                coords_metrics[v['Coords']] = [metric]
            else:
                coords_metrics[v['Coords']].append(metric)
        found_overlap = False
        for coords, metrics in coords_metrics.items():
            if len(metrics) > 1:
                found_overlap = True
                print(f"{metrics}: optimum coordinates are {coords}.")
        if not found_overlap:
            print('None.')
    
    
    # -----------------------------------------------------------------------------
    # Determine which optimization metrics to show.
    # -----------------------------------------------------------------------------
    metrics_to_opt = []
    for x_label, _ in row_parameter_pairs:
        for m in get_metrics_to_opt_for_x_label(x_label):
            if m not in metrics_to_opt:
                metrics_to_opt.append(m)
    
    metric_marker_styles, _metric_legend_labels = build_metric_marker_styles(metrics_to_opt, row_parameter_pairs)
    
    rowwise_optima, missing_opt_metric_files = compute_rowwise_optima(
        metrics_to_opt, row_parameter_pairs, strategies, z_label, steps
    )
    
    # -----------------------------------------------------------------------------
    # Load panel data for the plotted metric.
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
        raise FileNotFoundError('No matching metric CSV files were found for the requested row/column layout.')
    
    w_levels, cbar_ticks, w_ticks = compute_levels_from_arrays(all_arrays, metric)
    if not metric=='MPSP':
        w_ticks = []
    cmap, cmap_over_color, extend_cmap = choose_metric_colormap(metric)
    
    # -----------------------------------------------------------------------------
    # Figure layout
    # -----------------------------------------------------------------------------
    nrows = len(row_parameter_pairs)
    ncols = len(strategies)
    fig, axs = plt.subplots(nrows, ncols, constrained_layout=True)
    if nrows == 1 and ncols == 1:
        axs = np.array([[axs]])
    elif nrows == 1:
        axs = np.array([axs])
    elif ncols == 1:
        axs = np.array([[ax] for ax in axs])
    
    axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
    default_fontsize = 11.
    clabel_fontsize = 11
    axis_tick_fontsize = 11
    
    for i, (x_label, y_label) in enumerate(row_parameter_pairs):
        spec_1, spec_2, spec_3, x_units, y_units, z_units, x_ticks, y_ticks, z_ticks = get_specs_and_ticks(
            x_label, y_label, z_label, steps
        )
        x_label_for_plot = format_param_label(x_label)
        y_label_for_plot = format_param_label(y_label)
        w_label = r"$\bf" + get_metric_plot_name(metric, x_label).replace(' ', '\\ ') + "$"
    
        for j, strategy in enumerate(strategies):
            ax = axs[i, j]
            arr = panel_data[(x_label, y_label)][strategy]
    
            if arr is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, 'File not found', transform=ax.transAxes, ha='center', va='center', fontsize=10)
                continue
    
            additional_points = {}
            for opt_metric, opt_record in rowwise_optima.get(i, {}).items():
                if opt_record['strategy'] == strategy:
                    # For data with x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie'), only ethanol is produced; however, etoh yield is saved as 'combined yield' rather than 'etoh yield'
                    if opt_metric=='Combined Yield' and x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie'):
                        try:
                            additional_points[(opt_record['x'], opt_record['y'])] = metric_marker_styles['EtOH Yield']
                        except:
                            additional_points[(opt_record['x'], opt_record['y'])] = ('p', 'w', 8)
                    else:
                        additional_points[(opt_record['x'], opt_record['y'])] = metric_marker_styles[opt_metric]
            
            # add baseline marker
            baseline_coords = None
            if x_label in ('k_1e',) and y_label in ('k_1ie', 'k_7ie',):
                baseline_coords = (47.1, 0.04)
            elif x_label in ('k_13',) and y_label in ('k_7ii',):
                # baseline_coords = (5.81, 0.04)
                # baseline_coords = (0.0, 0.04)
                pass
            if baseline_coords is not None:
                additional_points[baseline_coords] = ('D', 'gray', 6)
                
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
                # include_x_axis_ticklabels=(i == nrows - 1),
                include_x_axis_ticklabels=True,
                # include_last_x_axis_ticklabel=(j==ncols-1),
                include_last_x_axis_ticklabel=True,
                include_y_axis_ticklabels=(j == 0),
                additional_points=additional_points,
                fig_ax_to_use=(fig, ax),
                comparison_range=comparison_range,
                show_comparison_range_clabels=False,
                inline_spacing=0.1,
                label_over_color='black',
            )
    
            if i == 0:
                ax.set_title(strategy.replace('_', '\n'), fontsize=11, fontweight='bold')
    
    # -----------------------------------------------------------------------------
    # Shared colorbar and layout
    # -----------------------------------------------------------------------------
    plt.subplots_adjust(wspace=0.15, hspace=0.6)
    
    norm = Normalize(vmin=float(w_levels[0]), vmax=float(w_levels[-1]))
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axs.ravel().tolist(), shrink=0.95, pad=0.02)
    cbar.set_label(f"{get_metric_plot_name(metric, row_parameter_pairs[0][0]).title()} [{metrics_units[metric]}]", fontsize=11)
    cbar.set_ticks(cbar_ticks)
    
    fig.set_figwidth(3.2 * ncols + 1.2)
    fig.set_figheight(2.7 * nrows + 1.6)
    
    plt.savefig(output_filename, transparent=False, facecolor='white', bbox_inches='tight', dpi=600)
    plt.close(fig)
    
    print(f'Saved multi-panel figure to: {output_filename}')
    if any(rowwise_optima.values()):
        print('\nRow-wise optima used for markers:')
        for row_i, (x_label, y_label) in enumerate(row_parameter_pairs):
            if not rowwise_optima.get(row_i):
                continue
            print(f'  Row {row_i + 1} ({x_label}, {y_label}):')
            for metric_name in metrics_to_opt:
                if metric_name not in rowwise_optima[row_i]:
                    continue
                opt_record = rowwise_optima[row_i][metric_name]
                print(
                    '    {metric}: value={value}, x={x}, y={y}, panel=({xl}, {yl}, {st})'.format(
                        metric=metric_name,
                        value=get_rounded_str(opt_record['value'], 6),
                        x=get_rounded_str(opt_record['x'], 6),
                        y=get_rounded_str(opt_record['y'], 6),
                        xl=opt_record['x_label'],
                        yl=opt_record['y_label'],
                        st=opt_record['strategy'],
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
    
    rowwise_opt_coords_metric_values, missing_summary_files = build_rowwise_opt_coords_metric_values(
        rowwise_optima=rowwise_optima,
        row_parameter_pairs=row_parameter_pairs,
        metrics_to_opt=metrics_to_opt,
        z_label=z_label,
        steps=steps,
    )
    
    for row_i, (x_label, y_label) in enumerate(row_parameter_pairs):
        row_summary = rowwise_opt_coords_metric_values.get(row_i, {})
        if not row_summary:
            continue
        print(f"\n===== Optima comparison summary for row {row_i + 1} ({x_label}, {y_label}) =====")
        get_optima_comparisons(row_summary, rel_to_m='MPSP')
    
    if missing_summary_files:
        print('\nMissing summary files:')
        for fp in missing_summary_files:
            print(f'  - {fp}')
    
    return panel_data

#%% Plot for all metrics
results = {}
if __name__ == '__main__':
    for metric in list(metrics_units.keys()):
        results[metric] = plot_save_kinetic_multipanel(metric)
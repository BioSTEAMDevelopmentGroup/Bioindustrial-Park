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
from matplotlib.lines import Line2D, TICKDOWN, TICKLEFT
from matplotlib.transforms import ScaledTranslation
import matplotlib.text
import matplotlib.patheffects
import contourpy

import contourplots
from biosteam.utils import colors
from biorefineries import isobutanol
# No isobutanol.load() here: this script only reads the sweep CSVs and needs
# the package path, and the import is build-free (2026-08-30 lazy-load
# refactor), so it never touches the simulation / numba cache.

#%%
# -----------------------------------------------------------------------------
# User inputs
# -----------------------------------------------------------------------------
# Metric to plot; override with IBO_FEED_STRAT_METRIC=IRR (or 'IBO MPSP', ...).
metric = os.environ.get('IBO_FEED_STRAT_METRIC', 'MPSP')
x_label = "Threshold glucose concentration" # title of the x axis
# units stay in mathtext (rendered in the figure typeface, see
# apply_font_rcparams) with a middle dot: Arial has no dot operator (\cdot,
# U+22C5) or superscript minus, which would fall back to another face
x_units =r"$\mathrm{g·L}^{-1}$"
x_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

y_label = "Target glucose concentration" # title of the y axis
y_units =r"$\mathrm{g·L}^{-1}$"
y_ticks = [0, 100, 200, 300, 400,
           # 300, 400, 500,
           ]

z_label = "Max. no. of glucose spikes" # optimized per grid point in these sweeps; the z axis is not drawn
z_units = ""
z_ticks = [0, 1]

steps = (25, 25, 1)

spec_1 = threshold_conces = np.linspace(1., 400., steps[0])

spec_2 = target_conces = np.linspace(10., 400., steps[1])

# Panels: one column per scenario. Each reads the per-metric CSVs that a run of
# analyses/evaluate_feeding_strategies.py (spike cap optimized for ethanol MPSP
# at every threshold/target point) saved, copied into
# analyses/results/publication/Feed-strat/ under their
# 'ibo_(25, 25, 1)_Thres_Targe_Max n_<scenario>_' prefix.
scenarios = ['B'] # override with IBO_FEED_STRAT_SCENARIOS=A or =A,B
scenarios = os.environ.get('IBO_FEED_STRAT_SCENARIOS', ','.join(scenarios)).split(',')

output_filename = f'{metric}_multi_panel_feed_strat_{"".join(scenarios)}.png'

# Optional: set to None to auto-compute from all loaded panels.
if metric == 'MPSP' and scenarios == ['A']:
    # scenario A spans ~0.87-1.10 $/kg on the 2026-09-02 model (optimum = the
    # 0.866 baseline strategy)
    manual_w_levels = np.arange(0.85, 1.1001, 0.005)
    manual_cbar_ticks = np.arange(0.85, 1.1001, 0.05)
    manual_w_ticks = [0.90, 0.95, 1.00]
elif metric == 'MPSP':
    # scenario B spans ~0.66-2.06 $/kg (baseline ~1.08); keep the level count
    # within the colormap's 90 colors
    manual_w_levels = np.arange(0.6, 2.0001, 0.02)
    manual_cbar_ticks = np.arange(0.6, 2.0001, 0.2)
    manual_w_ticks = [0.8, 1.2, 1.6]
elif metric == 'IBO MPSP':
    # scenario B spans ~1.29-3.01 $/kg (baseline ~1.80); A makes no isobutanol
    manual_w_levels = np.arange(1.2, 3.0001, 0.025)
    manual_cbar_ticks = np.arange(1.2, 3.0001, 0.3)
    manual_w_ticks = [1.5, 2.0, 2.5]
else:
    # IRR (and any other metric): fitted to the loaded panels in
    # compute_levels_from_arrays
    manual_w_levels = manual_cbar_ticks = manual_w_ticks = None

comparison_range = []
# if metric == 'MPSP':
#     comparison_range = EtOH_market_range = np.array([
#         0.52, # 1.5475 $/gal/(3.7854 L/gal * 0.789 kg/L)
#         1.15, # 3.4500 $/gal/(3.7854 L/gal * 0.789 kg/L)
#         ]) # Jan 2021 - Dec 2025 5-year low and high from https://tradingeconomics.com/commodity/ethanol

get_rounded_str = contourplots.utils.get_rounded_str
if metric in ('MPSP', 'IBO MPSP', 'IRR'):
    # fixed two decimals: the colorbar ticks use the same formatter, so the
    # bar and the inline labels agree (0.80 on both, not 0.8 vs 0.80)
    fmt_clabel = lambda cvalue: f"{cvalue:.2f}"
else:
    fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3)

# -----------------------------------------------------------------------------
# Typeface, font sizes, tick style
# -----------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'panel_title': 12, 'cbar_title': 12,
         'clabel': 10, 'legend': 9}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out


def apply_font_rcparams():
    # contourplots resets font.sans-serif to the uninstalled 'Arial Unicode'
    # (silent DejaVu Sans fallback) and font.size on every call, so this is
    # applied before the figure is made and again after each panel call.
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def restyle_panel_text(ax):
    # text created inside animated_contourplot (tick labels, contour labels)
    # picked up the library's rcParams; force the figure typeface on it
    for t in ax.findobj(matplotlib.text.Text):
        t.set_fontfamily(FONT_FAMILY)


def style_ticks(ax):
    # top/right: inward only; left/bottom: in and out, the same length each
    # way. Call after fig.canvas.draw() so every tick object exists.
    for which, L in TICK_LEN.items():
        ax.tick_params(axis='both', which=which, direction='inout', length=2*L,
                       top=True, right=True, labelsize=FONTS['tick'])
        get = 'get_major_ticks' if which == 'major' else 'get_minor_ticks'
        for tick in getattr(ax.xaxis, get)():
            tick.tick2line.set_marker(TICKDOWN)
            tick.tick2line.set_markersize(L)
        for tick in getattr(ax.yaxis, get)():
            tick.tick2line.set_marker(TICKLEFT)
            tick.tick2line.set_markersize(L)

# -----------------------------------------------------------------------------
# Shared metadata from across_kinetic_params.py
# -----------------------------------------------------------------------------
metrics_units = {
    'MPSP': r"$\mathrm{\$·kg}^{-1}$",
    'IBO MPSP': r"$\mathrm{\$·kg}^{-1}$",
    'IRR': '',
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
    if metric_name not in ('TCI', 'AOC', 'MPSP', 'IBO MPSP', 'IRR'):
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
        'MPSP', 'IBO MPSP', 'AOC', 'TCI', 'Total Q sugar evap',
        'Fermentation time', 'Actual aeration required'
    )


def format_param_label(label):
    return r"$\bf{label.replace(' ', '\ ')}$"


def format_axis_label_with_units(label, units):
    base = format_param_label(label)
    if units:
        return f"{base} [{units}]"
    return base



def load_metric_array(x_label, y_label, scenario, metric_name, steps):
    subfolder_name = 'Feed-strat\\'
    file_prefix = f'ibo_{steps}_{x_label[:5]}_{y_label[:5]}_Max n_{scenario}_'
    folder = isobutanol_results_pub_filepath + subfolder_name
    filepath = os.path.join(folder, file_prefix + f'_{metric_name}.csv')
    if not os.path.exists(filepath):
        return None, filepath
    df = pd.read_csv(filepath)
    arr_2d = df[df.columns[1:]].to_numpy(dtype=float, na_value=np.nan)
    return np.array([arr_2d]), filepath


def choose_metric_colormap(metric_name):
    lccm = metric_name.lower()
    cmap_under_color = None
    if 'irr' in lccm:
        # higher is better (reversed map: yellow = high, dark grey = 0). The
        # colorbar starts at 0 (compute_levels_from_arrays), so every
        # money-losing point (IRR < 0, incl. -inf = unsolvable) is painted
        # in one flat "under zero" colour: a grey darker than the map's
        # darkest colour (grey_dark, 0.28 -> shade(20), 0.23).
        cmap = JBEI_UCB_colormap(reverse=True)
        cmap_over_color = colors.yellow_tint.RGBn
        cmap_under_color = colors.grey_dark.shade(20).RGBn
        extend_cmap = 'both'
    elif any(key in lccm for key in ('yield', 'titer', 'productivity', 'loading')):
        cmap = JBEI_UCB_colormap(reverse=True)
        cmap_over_color = colors.yellow_tint.RGBn
        extend_cmap = 'neither'
    else:
        cmap = JBEI_UCB_colormap(reverse=False)
        cmap_over_color = colors.grey_dark.shade(8).RGBn
        extend_cmap = 'max'
    return cmap, cmap_over_color, cmap_under_color, extend_cmap


def compute_levels_from_arrays(arrays, metric_name):
    valid = []
    for arr in arrays:
        if arr is None:
            continue
        # finite values only: NaN = not simulated / not solved, and IRR is
        # -inf where no real IRR exists (solve_TEA, since 2026-09-03)
        vals = arr[np.isfinite(arr)]
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
    if 'irr' in lccm:
        # Fit the IRR colorbar to the loaded panels (same recipe as the
        # sweep script's IRR block): finest level step that keeps the filled
        # contours within the colormap's 90 colors, bounds rounded outward
        # to that step, colorbar ticks every 10 steps. The lower bound is
        # floored at 0 so the money-losing corners (IRR reaches -0.7 at high
        # target concs in scenario B) don't stretch the scale: every IRR < 0
        # cell takes the flat "under zero" colour (choose_metric_colormap),
        # and the 0.00 level is the boundary of that region. Scenario A
        # spans ~-0.17-0.12, B ~-0.7-0.21.
        floor = 0.0
        for step in (0.005, 0.01, 0.02, 0.025, 0.05, 0.1):
            lb = max(np.floor(min_val/step)*step, floor)
            ub = np.ceil(max_val/step)*step
            if (ub - lb)/step <= 80: break
        levels = np.arange(lb, ub + step/10, step)
        cbar_ticks = np.arange(lb, ub + step/10, 10*step)
        # labeled contours at the colorbar ticks inside the data range
        # (break-even 0.00 included); quartile-based labels sit too close
        # together for scenario A's narrow 0.04-0.12 spread
        w_ticks = sorted(set(
            [0.0] + [float(round(t, 3)) for t in cbar_ticks if lb < t < max_val]))
        return levels, cbar_ticks, w_ticks
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


def build_metric_marker_styles(metrics_to_opt):
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
        elif metric_name == 'IBO MPSP':
            style = ('*', 'w', 16)
        elif metric_name == 'IRR':
            style = ('D', '#33ccff', 12)
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


def compute_global_optima(metrics_to_opt, scenarios, steps):
    global_optima = {}
    missing_opt_metric_files = []
    for opt_metric in metrics_to_opt:
        best_value = None
        best_record = None
        for scenario in scenarios:
            arr, filepath = load_metric_array(x_label, y_label, scenario, opt_metric, steps)
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
                    'x_label': x_label,
                    'y_label': y_label,
                    'scenario': scenario,
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

# minimum clearance a label wants from each obstacle class, in axes-span units
# (a 10 pt "0.20" label is ~0.1 of a 3.2 in panel wide, a 9 pt marker ~0.04)
CLEARANCE = {'marker': 0.13, 'label': 0.12, 'edge': 0.06, 'line': 0.05}

def place_contour_labels(arr2d, levels, marker_xy=(), line_y=None):
    """One inline label per level: the contour vertex with the largest
    clearance (distance / CLEARANCE, in axes-span units) from the optimum
    markers, the labels already placed, the panel edges and the white
    line. Vertices within 6 % of arc length of a segment end are skipped so
    the label is not cut by the feasible-region edge. Returns
    {level: (x, y)} for the levels that have a contour in this panel (the
    others are dropped, so the library draws no automatic labels)."""
    x0, x1 = x_ticks[0], x_ticks[-1]
    y0, y1 = y_ticks[0], y_ticks[-1]
    to_n = lambda x, y: np.array([(x - x0)/(x1 - x0), (y - y0)/(y1 - y0)])
    markers = [to_n(*p) for p in marker_xy]
    line_pts = []
    if line_y is not None:
        line_pts = [to_n(x, y) for x, y in zip(spec_1, line_y) if np.isfinite(y)]
    gen = contourpy.contour_generator(spec_1, spec_2, np.ma.masked_invalid(arr2d))
    placed, positions = [], {}

    def clearance(p):
        d = min(p[0], 1 - p[0], p[1], 1 - p[1]) / CLEARANCE['edge']
        for cls, pts in (('marker', markers), ('label', placed), ('line', line_pts)):
            for q in pts:
                d = min(d, np.hypot(*(p - q)) / CLEARANCE[cls])
        return d

    for level in levels:
        segs = [seg for seg in gen.lines(level) if len(seg) >= 2]
        seg_npts = [np.array([to_n(x, y) for x, y in seg]) for seg in segs]
        seg_arcs = [np.concatenate([[0.], np.cumsum(np.hypot(*np.diff(npts, axis=0).T))])
                    for npts in seg_npts]
        # candidates: only segments at least half as long as the level's
        # longest one (and never shorter than a ~0.12-wide label when a
        # proper segment exists), so a clear little loop never outranks the
        # main contour line
        longest = max((arc[-1] for arc in seg_arcs), default=0.)
        best_seg, best_k, best_score = None, None, -np.inf
        for seg, npts, arc in zip(segs, seg_npts, seg_arcs):
            if arc[-1] == 0 or arc[-1] < 0.5*longest or (arc[-1] < 0.12 and longest >= 0.12):
                continue
            frac = arc / arc[-1]
            inner = (frac >= 0.06) & (frac <= 0.94)
            if not inner.any():
                inner[len(inner)//2] = True
            for k in np.flatnonzero(inner):
                # clearance first, then a preference for long segments and
                # for the middle of a segment
                score = (clearance(npts[k]) + 0.5*min(arc[-1], 0.5)
                         + 0.05*(1 - abs(frac[k] - 0.5)*2))
                if score > best_score:
                    best_seg, best_k, best_score = seg, k, score
        if best_seg is not None:
            positions[level] = (float(best_seg[best_k][0]), float(best_seg[best_k][1]))
            placed.append(to_n(*positions[level]))
    return positions


def get_optima_metric_values(global_optima, metrics_to_opt, steps):
    opt_coords_metric_values = {}
    missing_metric_value_files = []
    for metric_name in metrics_to_opt:
        if metric_name not in global_optima:
            continue

        opt_record = global_optima[metric_name]
        opt_x = float(opt_record['x'])
        opt_y = float(opt_record['y'])
        opt_scenario = opt_record['scenario']

        x_matches = np.where(np.isclose(spec_1, opt_x))[0]
        y_matches = np.where(np.isclose(spec_2, opt_y))[0]
        if not len(x_matches) or not len(y_matches):
            continue
        opt_x_ind = int(x_matches[0])
        opt_y_ind = int(y_matches[0])

        coords = (opt_x, opt_y, opt_scenario)
        opt_coords_metric_values[metric_name] = {'Coords': coords}
        for comparison_metric in metrics_to_opt:
            arr, filepath = load_metric_array(x_label, y_label, opt_scenario, comparison_metric, steps)
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
            if np.isinf(curr_val):
                print(f"'{m2}' is {curr_val} (no real IRR on the valid domain).")
                continue
            if not np.isfinite(ref_val):
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
    'IBO MPSP',
    'IRR',
    ]

metric_marker_styles, _metric_legend_labels = build_metric_marker_styles(metrics_to_opt)

global_optima, missing_opt_metric_files = compute_global_optima(
    metrics_to_opt, scenarios, steps
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
for scenario in scenarios:
    arr, filepath = load_metric_array(x_label, y_label, scenario, metric, steps)
    panel_data[scenario] = arr
    if arr is None:
        missing_files.append(filepath)
    else:
        all_arrays.append(arr)

if not all_arrays:
    raise FileNotFoundError('No matching metric CSV files were found for the requested row/column layout.')

w_levels, cbar_ticks, w_ticks = compute_levels_from_arrays(all_arrays, metric)
cmap, cmap_over_color, cmap_under_color, extend_cmap = choose_metric_colormap(metric)

# -----------------------------------------------------------------------------
# Figure layout
# -----------------------------------------------------------------------------
nrows = 1
ncols = len(scenarios)
apply_font_rcparams()
fig, axs = plt.subplots(nrows, ncols, constrained_layout=True)
if nrows == 1 and ncols == 1:
    axs = np.array([[axs]])
elif nrows == 1:
    axs = np.array([axs])
elif ncols == 1:
    axs = np.array([[ax] for ax in axs])

axis_title_fonts={'size': {k: FONTS['axis_title'] for k in 'xyzw'},}
default_fontsize = FONTS['tick']
clabel_fontsize = FONTS['clabel']
axis_tick_fontsize = FONTS['tick']

for i in range(nrows):
    x_label_for_plot = format_param_label(x_label)
    y_label_for_plot = format_param_label(y_label)
    w_label = r"$\bf" + get_metric_plot_name(metric, x_label).replace(' ', '\\ ') + "$"

    for j, scenario in enumerate(scenarios):
        ax = axs[i, j]
        arr = panel_data[scenario]

        if arr is None:
            ax.set_axis_off()
            ax.text(0.5, 0.5, 'File not found', transform=ax.transAxes, ha='center', va='center', fontsize=10)
            continue
        if 'irr' in metric.lower():
            # unsolvable IRRs (-inf: no real root on the valid domain) belong
            # in the under-zero region, but contourf masks non-finite cells;
            # draw them just below the lower bound instead
            arr = np.where(np.isneginf(arr),
                           w_levels[0] - (w_levels[1] - w_levels[0]), arr)

        # Optima sharing a grid point keep the first metric's marker (in
        # metrics_to_opt order), so the MPSP star is never hidden under a
        # later co-located marker; overlaps are listed by get_optima_comparisons.
        additional_points = {}
        for opt_metric in metrics_to_opt:
            opt_record = global_optima.get(opt_metric)
            if opt_record is None or opt_record['scenario'] != scenario:
                continue
            additional_points.setdefault((opt_record['x'], opt_record['y']),
                                         metric_marker_styles[opt_metric])

        # White line: for each threshold, the lowest target at which the
        # MPSP-optimized spike cap lands on batch mode (zero spikes).
        add_lines = {}
        arr_n_spikes, filepath = load_metric_array(x_label, y_label, scenario, 'Number of glucose spikes', steps)
        if arr_n_spikes is not None:
            line_first_app_n_glu_spikes_0 = tuple(get_zeroth_spec_2_val_for_condition_for_all_spec_1_vals(arr_n_spikes[0],
                                                                              condition = lambda i: i==0))
            add_lines = {line_first_app_n_glu_spikes_0: {'color': 'white', 'linewidth': 1.0, 'alpha': 1.0}}

        # One label per labelled level, placed clear of markers, other labels,
        # the panel edges and the white line (the library's automatic
        # placement repeated levels and stacked labels on steep gradients).
        marker_xy = [(rec['x'], rec['y']) for rec in global_optima.values()
                     if rec['scenario'] == scenario]
        manual_clabels = place_contour_labels(
            arr[0], w_ticks, marker_xy,
            line_first_app_n_glu_spikes_0 if arr_n_spikes is not None else None)

        contourplots.animated_contourplot(
            w_data_vs_x_y_at_multiple_z=arr,
            x_data=spec_1,
            y_data=spec_2,
            z_data=[0.,],
            x_label=x_label_for_plot,
            y_label=y_label_for_plot,
            z_label=r"$\bf" + z_label + "$",
            w_label=w_label,
            x_ticks=x_ticks,
            y_ticks=y_ticks,
            z_ticks=z_ticks,
            w_levels=w_levels,
            w_ticks=list(manual_clabels.keys()), # only levels present in this panel
            manual_clabels_regular=manual_clabels,
            x_units=x_units,
            y_units=y_units,
            z_units=z_units,
            w_units=metrics_units[metric],
            fmt_clabel=fmt_clabel,
            cmap=cmap,
            cmap_over_color=cmap_over_color,
            cmap_under_color=cmap_under_color,
            extend_cmap=extend_cmap,
            cbar_ticks=cbar_ticks,
            fontname={'fontname': FONT_FAMILY},
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
            include_last_x_axis_ticklabel=True, # panels are visibly separated (layout wspace below)
            include_y_axis_ticklabels=(j == 0),
            include_last_y_axis_ticklabel=(i==0),
            additional_points=additional_points,
            fig_ax_to_use=(fig, ax),
            comparison_range=comparison_range,
            inline_spacing=0.,
            label_over_color='black',
            add_lines=add_lines,
        )
        apply_font_rcparams()
        restyle_panel_text(ax)

        if ncols > 1:
            # pad keeps the title clear of optimum markers sitting on the top edge
            ax.set_title(f'Scenario {scenario}', fontsize=FONTS['panel_title'],
                         fontweight='bold', pad=8)
        else:
            ax.set_title('') # drop the library's blank 20 pt spacer title

# -----------------------------------------------------------------------------
# Shared colorbar and layout
# -----------------------------------------------------------------------------
# (plt.subplots_adjust is a no-op under constrained layout; set the gap here)
fig.get_layout_engine().set(wspace=0.08, hspace=0.)

norm = Normalize(vmin=float(w_levels[0]), vmax=float(w_levels[-1]))
# carry the panels' under/over colours onto the shared colorbar as extension
# triangles (e.g. IRR's flat "under zero" grey, MPSP's dark over-range grey)
cbar_cmap = cmap.copy()
if cmap_under_color is not None:
    cbar_cmap.set_under(cmap_under_color)
if cmap_over_color is not None:
    cbar_cmap.set_over(cmap_over_color)
sm = ScalarMappable(norm=norm, cmap=cbar_cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs.ravel().tolist(), shrink=0.95, pad=0.02,
                    extend=extend_cmap)
_cbar_name = get_metric_plot_name(metric, x_label)
if not _cbar_name.isupper(): # capitalize plain names; keep acronyms (MPSP, TCI, AOC) intact
    _cbar_name = _cbar_name[0].upper() + _cbar_name[1:]
_cbar_units = metrics_units[metric]
cbar.set_label(f"{_cbar_name} [{_cbar_units}]" if _cbar_units else _cbar_name,
               fontsize=FONTS['cbar_title'])
cbar.set_ticks(cbar_ticks)
cbar.set_ticklabels([fmt_clabel(t) for t in cbar_ticks]) # same precision as the inline labels
cbar.ax.tick_params(labelsize=FONTS['tick'])
restyle_panel_text(cbar.ax)

fig.supxlabel(f"{x_label} [{x_units}]", fontsize=FONTS['axis_title'])
fig.supylabel(f"{y_label} [{y_units}]", fontsize=FONTS['axis_title'])

fig.set_figwidth(3.2 * ncols + 1.2)
fig.set_figheight(2.7 * nrows + 1.6)

fig.canvas.draw() # materialize every tick before restyling them
for ax in axs.ravel():
    style_ticks(ax)

output_filepath = os.path.join(isobutanol_results_pub_filepath, 'Feed-strat', output_filename)
plt.savefig(output_filepath, transparent=False, facecolor='white', bbox_inches='tight', dpi=600)
plt.close(fig)

print(f'Saved multi-panel figure to: {output_filepath}')
if global_optima:
    print('\nGlobal optima used for markers:')
    for metric_name in metrics_to_opt:
        if metric_name not in global_optima:
            continue
        opt_record = global_optima[metric_name]
        print(
            '  {metric}: value={value}, x={x}, y={y}, panel=({xl}, {yl}, scenario={sc})'.format(
                metric=metric_name,
                value=get_rounded_str(opt_record['value'], 6),
                x=get_rounded_str(opt_record['x'], 6),
                y=get_rounded_str(opt_record['y'], 6),
                xl=opt_record['x_label'],
                yl=opt_record['y_label'],
                sc=opt_record['scenario'],
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

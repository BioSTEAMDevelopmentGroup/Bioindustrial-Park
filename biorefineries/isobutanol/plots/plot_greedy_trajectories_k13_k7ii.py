#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
IRR contour of the scenario-B k_13 x k_7ii kinetic sweep, overlaid with greedy
(8-neighbor hill-climb) trajectories from the scenario-B baseline kinetics to
the local maxima of IRR and of the isobutanol and ethanol titer, yield and
productivity (seven trajectories; the metric set is configurable below).

Consumes the per-metric CSVs written by analyses/evaluate_EtOH_k13_k7ii.py
(20 x 20 grid, scenario B, no feeding-strategy optimization) and runs no
biorefinery simulation: the plotting module is loaded standalone, never via
`import biorefineries.isobutanol`.

Run:  & "$py" biorefineries/isobutanol/plots/plot_greedy_trajectories_k13_k7ii.py
Writes the PNG and a CSV of the trajectory points next to the sweep CSVs in
analyses/results/.
"""
import os
import importlib.util

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from biosteam.utils import colors

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, '..', 'analyses', 'results')

#%% Sweep being plotted (must match analyses/evaluate_EtOH_k13_k7ii.py)

# `file_to_save` prefix of the sweep run: ibo_{steps}_{x}_{y}_{z}_opt=..._max_n=..._
SWEEP_PREFIX = 'ibo_(20, 20, 1)_k_13_k_7ii_Spike_opt=False_max_n=13_'
STEPS = (20, 20)
SPEC_1 = np.linspace(0.0, 20.0, STEPS[0])       # k_13  (x axis, sweep spec_1)
SPEC_2 = np.linspace(0.0001, 0.2, STEPS[1])     # k_7ii (y axis, sweep spec_2)

# Scenario-B baseline kinetics (Baseline column of
# analyses/full/parameter_distributions/parameter-distributions_corn_IBO_EtOH_B.xlsx);
# snapped to the nearest grid cell by the plotting function.
BASELINE_K13 = 5.81
BASELINE_K7II = 0.15

COLOR_METRIC = 'IRR'
TRAJECTORY_METRICS = ['IRR',
                      'IBO Titer', 'IBO Yield', 'IBO Productivity',
                      'EtOH Titer', 'EtOH Yield', 'EtOH Productivity']
SENSES = {m: 'max' for m in TRAJECTORY_METRICS}
# Styling on the reversed (yellow = high IRR) colormap. Greedy climbs from one
# baseline share their first segments, so every trajectory needs a distinct
# look even where it lies on top of another:
#   color      = metric type   (black titer, white yield, magenta productivity)
#   line style = product       (solid isobutanol, dashed ethanol)
#   width      = draw order    (IRR widest and first, then titer > yield >
#                               productivity, so each narrower line stays
#                               visible inside a wider coincident one)
_TITER, _YIELD, _PROD = '#1a1a1a', '#ffffff', '#ff1aff'
TRAJECTORY_COLORS = {'IRR': '#00bfff',
                     'IBO Titer': _TITER, 'IBO Yield': _YIELD, 'IBO Productivity': _PROD,
                     'EtOH Titer': _TITER, 'EtOH Yield': _YIELD, 'EtOH Productivity': _PROD}
TRAJECTORY_LINESTYLES = {'IRR': '-',
                         'IBO Titer': '-', 'IBO Yield': '-', 'IBO Productivity': '-',
                         'EtOH Titer': '--', 'EtOH Yield': '--', 'EtOH Productivity': '--'}
TRAJECTORY_LINEWIDTHS = {'IRR': 5.0,
                         'IBO Titer': 3.4, 'IBO Yield': 2.0, 'IBO Productivity': 1.2,
                         'EtOH Titer': 3.4, 'EtOH Yield': 2.0, 'EtOH Productivity': 1.2}

#%% Plot styling (shared with the sweep script's IRR contour)

x_label = r"$\mathbf{k}_{13}$"
y_label = r"$\mathbf{k}_{7ii}$"
xy_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
x_ticks = [0, 5, 10, 15, 20]
y_ticks = [0.0, 0.05, 0.1, 0.15, 0.2]

# same bounds as evaluate_EtOH_k13_k7ii.py (grid IRR spans -0.12 to 0.19)
IRR_w_levels = np.arange(-0.1, 0.2001, 0.005)
IRR_cbar_ticks = np.arange(-0.1, 0.2001, 0.05)
IRR_w_ticks = [0.0, 0.05, 0.10, 0.15, 0.18]

axis_title_fonts = {'size': {'x': 11, 'y': 11, 'z': 11, 'w': 11}}
default_fontsize = 11.
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5


#%% Helpers

def _load_module():
    path = os.path.join(HERE, 'trajectory_contourplot.py')
    spec = importlib.util.spec_from_file_location('trajectory_contourplot', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_metric_csv(metric):
    """(1, ny, nx) array from '<SWEEP_PREFIX>_<metric>.csv'; rows are spec_2
    (k_7ii) values and columns spec_1 (k_13) values, as the sweep saves them."""
    path = os.path.join(RESULTS_DIR, f'{SWEEP_PREFIX}_{metric}.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f'{path}\nRun analyses/evaluate_EtOH_k13_k7ii.py first (scenario B, '
            f'steps {STEPS}).')
    arr = pd.read_csv(path, index_col=0).to_numpy(dtype=float)
    if arr.shape != (STEPS[1], STEPS[0]):
        raise ValueError(f'{metric}: expected shape {(STEPS[1], STEPS[0])}, '
                         f'got {arr.shape}')
    return arr[None, :, :]


#%% Main

def main():
    tc = _load_module()
    results = {m: _load_metric_csv(m)
               for m in dict.fromkeys([COLOR_METRIC, *TRAJECTORY_METRICS])}

    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, SPEC_1, SPEC_2,
        color_metric=COLOR_METRIC,
        color_metric_label=r"$\mathbf{IRR}$",
        color_metric_units='',
        baseline_point=(BASELINE_K13, BASELINE_K7II),
        trajectory_metrics=TRAJECTORY_METRICS,
        senses=SENSES,
        trajectory_colors=TRAJECTORY_COLORS,
        trajectory_linestyles=TRAJECTORY_LINESTYLES,
        trajectory_linewidths=TRAJECTORY_LINEWIDTHS,
        # eight entries: park the legend below the axes, three columns; the
        # mid-grey face keeps both the black and the white lines visible
        legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, -0.16),
                       'ncol': 3, 'fontsize': 8.5, 'framealpha': 1.0,
                       'facecolor': '#a9a9a9', 'edgecolor': 'k'},
        x_label=x_label, y_label=y_label,
        x_units=xy_units, y_units=xy_units,
        x_ticks=x_ticks, y_ticks=y_ticks,
        cmap=tc.JBEI_UCB_colormap(reverse=True),
        w_levels=IRR_w_levels, cbar_ticks=IRR_cbar_ticks, w_ticks=IRR_w_ticks,
        extend_cmap='both',
        cmap_over_color=colors.yellow_tint.RGBn,
        cmap_under_color=colors.grey_dark.shade(40).RGBn,
        # passed through to contourplots.animated_contourplot
        fmt_clabel=lambda cvalue: f'{cvalue:.2f}',
        axis_title_fonts=axis_title_fonts,
        clabel_fontsize=clabel_fontsize,
        default_fontsize=default_fontsize,
        axis_tick_fontsize=axis_tick_fontsize,
        n_minor_ticks=1,
        cbar_n_minor_ticks=4,
        round_yticks_to=2,
        units_on_newline=(False, False, False, False),
        # x, y, z, w: no brackets around the (empty) IRR units
        units_opening_brackets=[" (", " (", " (", ""],
        units_closing_brackets=[")", ")", ")", ""],
    )

    stem = f'{COLOR_METRIC}_greedy_trajectories_{SWEEP_PREFIX}'
    png = os.path.join(RESULTS_DIR, stem + '.png')
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    print('Saved:', png)

    # trajectory points, one row per step, for the record
    rows = []
    for m, d in traj.items():
        grid = results[m][0]
        for step, ((x, y), (iy, ix)) in enumerate(zip(d['path_xy'], d['path_ij'])):
            rows.append({'metric': m, 'sense': SENSES[m], 'step': step,
                         'k_13': x, 'k_7ii': y, 'value': grid[iy, ix]})
    csv = os.path.join(RESULTS_DIR, stem + '.csv')
    pd.DataFrame(rows).to_csv(csv, index=False)
    print('Saved:', csv)

    ix0 = tc._snap_index(SPEC_1, BASELINE_K13)
    iy0 = tc._snap_index(SPEC_2, BASELINE_K7II)
    print(f'\nBaseline (k_13, k_7ii) = ({BASELINE_K13}, {BASELINE_K7II}) '
          f'snapped to grid cell ({SPEC_1[ix0]:.4g}, {SPEC_2[iy0]:.4g})')
    for m in TRAJECTORY_METRICS:
        print(f'  {m:18s} {SENSES[m]}: {results[m][0][iy0, ix0]:.4g} at baseline')
    print()
    for m, d in traj.items():
        print(f'{m:18s} ({SENSES[m]}): {d["n_steps"]} steps -> optimum '
              f'{d["optimum_value"]:.4g} at (k_13, k_7ii) = '
              f'({d["optimum_xy"][0]:.4g}, {d["optimum_xy"][1]:.4g})')
    return fig, ax, traj


if __name__ == '__main__':
    main()

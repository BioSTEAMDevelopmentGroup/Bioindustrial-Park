#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
IRR contour of the scenario-A k_1e x inhib_ethanol-multiplier kinetic sweep
(enzyme burden ON), overlaid with greedy (8-neighbor hill-climb) trajectories
from the scenario-A baseline cell (k_1e = 47.1, inhib_ethanol multiplier = 1.0)
to the local maxima of IRR and of the ethanol titer, combined yield and
ethanol productivity (four trajectories; the metric set is configurable below).

Consumes the per-metric CSVs written by
analyses/evaluate_EtOH_k1e_inhib_ethanol.py (20 x 20 grid, scenario A, enzyme
burden ON, no feeding-strategy optimization) and runs no biorefinery
simulation: the plotting engine (plots/trajectory_contourplot.py) is loaded
standalone, never via `import biorefineries.isobutanol`.

The high-k_1e region is blank (NaN): those points are enzyme-burden-infeasible
(the ethanol-production capacity pushes the modeled proteome pool over the
flexible-sector cap). Unsolvable, money-losing cells (IRR = -inf) are drawn in
the colormap's under-color; the greedy max-climb never steps onto them.

Run:  & "$py" biorefineries/isobutanol/plots/plot_greedy_trajectories_k1e_inhib_ethanol.py
Writes the PNG and a CSV of the trajectory points next to the sweep CSVs in
analyses/results/.
"""
import os
import importlib.util

import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
from biosteam.utils import colors

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, '..', 'analyses', 'results')

#%% Sweep being plotted (must match analyses/evaluate_EtOH_k1e_inhib_ethanol.py)

# `file_to_save` prefix of the sweep run: ibo_{steps}_{x}_{y}_{z}_opt=..._max_n=..._
# (x_label[:5]='k_1e', y_label[:5]='inhib', z_label[:5]='Spike'; max_n = 16 for
# scenario A's fed-batch feeding strategy). The loader appends '_<metric>.csv'.
SWEEP_PREFIX = 'ibo_(20, 20, 1)_k_1e_inhib_Spike_opt=False_max_n=16_'
STEPS = (20, 20)
SPEC_1 = np.linspace(1.0, 300.0, STEPS[0])     # k_1e   (x axis, sweep spec_1)
SPEC_2 = np.linspace(0.2, 2.0, STEPS[1])       # inhib_ethanol multiplier (y, spec_2)

# Trajectory start point, snapped to the nearest grid cell by the plotting
# function: the scenario-A baseline -- k_1e = 47.1 (Baseline column of
# analyses/full/parameter_distributions/parameter-distributions_corn_IBO_EtOH_A.xlsx,
# read live off the model during the sweep) with the inhib_ethanol family at
# its baseline (multiplier = 1.0).
BASELINE_K1E = 47.1
BASELINE_MULT = 1.0

COLOR_METRIC = 'IRR'
TRAJECTORY_METRICS = ['IRR', 'EtOH Titer', 'Combined Yield', 'EtOH Productivity']
SENSES = {m: 'max' for m in TRAJECTORY_METRICS}
# Styling on the reversed (yellow = high IRR) colormap: IRR white, the three
# fermentation metrics one grey and told apart by their optimum markers. All
# lines are dashed at the same width; greedy climbs from one start share their
# first segment, so the IRR dashes are phase-shifted by one dash length against
# the grey ones and the two alternate where the paths coincide.
_GREY = '#a0a0a0'
TRAJECTORY_COLORS = {'IRR': '#ffffff', 'EtOH Titer': _GREY,
                     'Combined Yield': _GREY, 'EtOH Productivity': _GREY}
_DASH_IRR, _DASH_OTHER = (0, (4, 4)), (4, (4, 4))   # (offset, (on, off)) in pt
TRAJECTORY_LINESTYLES = {m: (_DASH_IRR if m == 'IRR' else _DASH_OTHER)
                         for m in TRAJECTORY_METRICS}
TRAJECTORY_LINEWIDTHS = {m: 1.5 for m in TRAJECTORY_METRICS}
TRAJECTORY_MARKERS = {'IRR': '*', 'EtOH Titer': 'o',
                      'Combined Yield': 's', 'EtOH Productivity': '^'}
TRAJECTORY_MARKER_SIZES = {'IRR': 12, 'EtOH Titer': 7,
                           'Combined Yield': 6.5, 'EtOH Productivity': 7.5}   # pt
BASELINE_MARKER = ('D', 'white', 6)   # (shape, fill color, size)

#%% Plot styling (shared with the sweep script's IRR contour)

x_label = r"$\mathbf{k}_{1e}$"
y_label = r"$\mathbf{inhib_{ethanol}\ multiplier}$"
x_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
y_units = r""    # dimensionless (x scenario-A baseline of each member)
x_ticks = [0, 100, 200, 300]
y_ticks = [0.2, 0.6, 1.0, 1.4, 1.8]

# IRR is plotted in percent. The grid's finite range is ~ -68 % to 13 %; the
# colorbar covers -10 % to 14 % and the extend arrows catch both ends (the
# under-color also absorbs the -inf money-losing cells cleaned in main()).
IRR_w_levels = np.arange(-10.0, 14.001, 0.5)
IRR_cbar_ticks = np.arange(-10.0, 14.001, 5.0)     # -10, -5, 0, 5, 10
IRR_w_ticks = [-5.0, 0.0, 5.0, 10.0, 12.0]         # black labeled contour lines
fmt_percent = lambda v, pos=None: f'{v:g}%'

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
    (multiplier) values and columns spec_1 (k_1e) values, as the sweep saves
    them."""
    path = os.path.join(RESULTS_DIR, f'{SWEEP_PREFIX}_{metric}.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f'{path}\nRun analyses/evaluate_EtOH_k1e_inhib_ethanol.py first '
            f'(scenario A, steps {STEPS}).')
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
    results['IRR'] = 100. * results['IRR']   # fraction -> percent

    # Money-losing cells solve_TEA reports as IRR = -inf: push them just below
    # the lowest level so contourf fills them with the under-color instead of
    # tripping on a non-finite value. For a max-climb they never improve, so
    # cleaning them to a finite low value leaves every trajectory unchanged.
    irr = results['IRR']
    if np.isneginf(irr).any():
        irr[np.isneginf(irr)] = IRR_w_levels[0] - (IRR_w_levels[1] - IRR_w_levels[0])

    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, SPEC_1, SPEC_2,
        color_metric=COLOR_METRIC,
        color_metric_label=r"$\mathbf{Internal\ Rate\ of\ Return}$",
        color_metric_units='',
        baseline_point=(BASELINE_K1E, BASELINE_MULT),
        trajectory_metrics=TRAJECTORY_METRICS,
        senses=SENSES,
        trajectory_colors=TRAJECTORY_COLORS,
        trajectory_linestyles=TRAJECTORY_LINESTYLES,
        trajectory_linewidths=TRAJECTORY_LINEWIDTHS,
        trajectory_markers=TRAJECTORY_MARKERS,
        trajectory_marker_sizes=TRAJECTORY_MARKER_SIZES,
        baseline_marker=BASELINE_MARKER,
        # a step must beat the sweep's own convergence tolerance so a plateau
        # of round-off does not read as a climb (same 1e-4 as the k13_k7ii plot)
        min_rel_improvement=1e-4,
        # legend below the axes in one row; the light-grey face keeps both the
        # white and the grey lines visible
        legend_kwargs={'loc': 'upper center', 'bbox_to_anchor': (0.5, -0.16),
                       'ncol': 5, 'fontsize': 8.5, 'framealpha': 1.0,
                       'facecolor': '#c8c8c8', 'edgecolor': 'k'},
        x_label=x_label, y_label=y_label,
        x_units=x_units, y_units=y_units,
        x_ticks=x_ticks, y_ticks=y_ticks,
        cmap=tc.JBEI_UCB_colormap(reverse=True),
        w_levels=IRR_w_levels, cbar_ticks=IRR_cbar_ticks, w_ticks=IRR_w_ticks,
        extend_cmap='both',
        cmap_over_color=colors.yellow_tint.RGBn,
        cmap_under_color=colors.grey_dark.shade(40).RGBn,
        # passed through to contourplots.animated_contourplot
        fmt_clabel=fmt_percent,
        axis_title_fonts=axis_title_fonts,
        clabel_fontsize=clabel_fontsize,
        default_fontsize=default_fontsize,
        axis_tick_fontsize=axis_tick_fontsize,
        n_minor_ticks=1,
        cbar_n_minor_ticks=4,
        round_yticks_to=1,
        units_on_newline=(False, False, False, False),
        # x, y, z, w: brackets only around the (present) k_1e units
        units_opening_brackets=[" (", "", "", ""],
        units_closing_brackets=[")", "", "", ""],
    )
    # percent ticks on the colorbar (the last axes contourplots added)
    cbar_ax = [a for a in fig.axes if a is not ax][-1]
    cbar_ax.yaxis.set_major_formatter(FuncFormatter(fmt_percent))

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
                         'k_1e': x, 'inhib_ethanol_multiplier': y,
                         'value': grid[iy, ix]})
    csv = os.path.join(RESULTS_DIR, stem + '.csv')
    pd.DataFrame(rows).to_csv(csv, index=False)
    print('Saved:', csv)

    ix0 = tc._snap_index(SPEC_1, BASELINE_K1E)
    iy0 = tc._snap_index(SPEC_2, BASELINE_MULT)
    print(f'\nBaseline (k_1e, multiplier) = ({BASELINE_K1E}, {BASELINE_MULT}) '
          f'snapped to grid cell ({SPEC_1[ix0]:.4g}, {SPEC_2[iy0]:.4g}); '
          f'IRR values in %')
    for m in TRAJECTORY_METRICS:
        print(f'  {m:18s} {SENSES[m]}: {results[m][0][iy0, ix0]:.4g} at baseline')
    print()
    for m, d in traj.items():
        print(f'{m:18s} ({SENSES[m]}): {d["n_steps"]} steps -> optimum '
              f'{d["optimum_value"]:.4g} at (k_1e, multiplier) = '
              f'({d["optimum_xy"][0]:.4g}, {d["optimum_xy"][1]:.4g})')
    return fig, ax, traj


if __name__ == '__main__':
    main()

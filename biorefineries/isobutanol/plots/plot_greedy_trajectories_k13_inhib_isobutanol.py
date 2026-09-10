#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
IRR contour of the opt_IRR k_13 x inhib_isobutanol-multiplier kinetic sweep
(enzyme burden ON), overlaid with greedy (8-neighbor hill-climb) trajectories
from an ethanol-only start point (k_13 = 0, inhib_isobutanol multiplier = 1.0)
to the local maxima of six metrics: isobutanol yield/titer/productivity and
ethanol yield/titer/productivity (all sense 'max'). The financial (IRR) optimum
is NOT hill-climbed -- a greedy IRR climb is trapped at the money-losing
ethanol-only start (every 8-neighbor is also money-losing) -- so instead the
opt_IRR point (k_13 = opt_IRR baseline, multiplier = 1.0) is marked with a star
in the same blue the kinetic-optimization parameter-sets figure gives the IRR
study (HUE_COLORS[0] in plot_kin_opt_parameter_sets.py).

Consumes the per-metric CSVs written by
analyses/evaluate_EtOH_k13_inhib_isobutanol.py (20 x 20 grid, opt_IRR baseline,
no feeding-strategy optimization) and runs no biorefinery simulation: the
plotting module is loaded standalone, never via `import biorefineries.isobutanol`.

IRR styling matches the sweep script's IRR contour: shown in percent on a HARD
0-25 % scale, a grey under-color for money-losing cells (< 0 %, including the
unsolvable -inf corners, which are pushed just below the lowest level so they
fill with the under-color instead of rendering blank), a white labeled 0 %
break-even contour line, black labeled 5/10/15/20 % lines, and % contour labels.

Run:  & "$py" biorefineries/isobutanol/plots/plot_greedy_trajectories_k13_inhib_isobutanol.py
Writes the PNG and a CSV of the trajectory points next to the sweep CSVs in
analyses/results/.
"""
import os
import importlib.util

import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from biosteam.utils import colors

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(HERE, '..', 'analyses', 'results')

#%% Sweep being plotted (must match analyses/evaluate_EtOH_k13_inhib_isobutanol.py)

# `file_to_save` prefix of the sweep run: ibo_{steps}_{x}_{y}_{z}_opt=..._max_n=..._
# (x_label[:5]="k_13", y_label[:5]="inhib", z_label[:5]="Spike"; opt_IRR's
# default_max_n_glu_spikes is 18). The loader adds the second '_' before the
# metric name, as the sweep's `csv_file_to_save = file_to_save + f'_{k}'` does.
SWEEP_PREFIX = 'ibo_(20, 20, 1)_k_13_inhib_Spike_opt=False_max_n=18_'
STEPS = (20, 20)
SPEC_1 = np.linspace(0.0, 6.5, STEPS[0])    # k_13                    (x, spec_1)
SPEC_2 = np.linspace(0.2, 2.0, STEPS[1])    # inhib_isobutanol mult.  (y, spec_2)

# Trajectory start point, snapped to the nearest grid cell by the plotting
# function: the ethanol-only strain -- the isobutanol Ehrlich entry switched
# off (k_13 = 0) at the unscaled inhib_isobutanol family (multiplier = 1.0).
BASELINE_K13 = 0.0
BASELINE_MULT = 1.0

COLOR_METRIC = 'IRR'
# IRR is the contour color metric but is NOT hill-climbed (see the module
# docstring); only the six fermentation metrics get greedy trajectories.
TRAJECTORY_METRICS = ['IBO Yield', 'IBO Titer', 'IBO Productivity',
                      'EtOH Yield', 'EtOH Titer', 'EtOH Productivity']
SENSES = {m: 'max' for m in TRAJECTORY_METRICS}

# The six fermentation metrics take distinct high-contrast colors that read
# against the grey->blue->orange->yellow map, each with its own optimum marker.
# All lines are dashed at one width with staggered dash offsets, so where greedy
# climbs share a segment (they all start at the same cell) the colors interleave
# rather than one hiding the rest.
# EtOH Productivity is black (not teal) to stay clear of the opt_IRR marker's
# blue (#18C4DC), which the requested styling fixes.
TRAJECTORY_COLORS = {'IBO Yield': '#66dd55', 'IBO Titer': '#33ccff',
                     'IBO Productivity': '#b266ff',
                     'EtOH Yield': '#ff66b3', 'EtOH Titer': '#ff4d4d',
                     'EtOH Productivity': '#111111'}
TRAJECTORY_MARKERS = {'IBO Yield': 's', 'IBO Titer': 'o', 'IBO Productivity': '^',
                      'EtOH Yield': 'D', 'EtOH Titer': 'v', 'EtOH Productivity': 'P'}
TRAJECTORY_MARKER_SIZES = {'IBO Yield': 6.5, 'IBO Titer': 7, 'IBO Productivity': 7.5,
                           'EtOH Yield': 6, 'EtOH Titer': 7, 'EtOH Productivity': 7.5}
_PERIOD = 8.0   # dash period (4 on + 4 off), pt
TRAJECTORY_LINESTYLES = {m: (i * _PERIOD / len(TRAJECTORY_METRICS), (4, 4))
                         for i, m in enumerate(TRAJECTORY_METRICS)}
TRAJECTORY_LINEWIDTHS = {m: 1.6 for m in TRAJECTORY_METRICS}
BASELINE_MARKER = ('D', 'white', 7)   # (shape, fill color, size)

# The financial (IRR) optimum: the opt_IRR point at (opt_IRR baseline k_13,
# multiplier 1.0), drawn as a star in the IRR study's blue from the
# kinetic-optimization parameter-sets figure (HUE_COLORS[0]). Not hill-climbed;
# k_13 is opt_IRR's live baseline printed by the sweep (scenarios' opt_IRR
# workbook), placed at its true data coordinates (not snapped to a grid cell).
OPT_IRR_K13 = 0.6421409697576548
OPT_IRR_MULT = 1.0
OPT_IRR_COLOR = '#18C4DC'
OPT_IRR_MARKER = ('*', OPT_IRR_COLOR, 15)   # (shape, fill color, size)

#%% Plot styling (shared with the sweep script's IRR contour)

x_label = r"$\mathbf{k}_{13}$"
y_label = r"$\mathbf{inhib\_isobutanol\ multiplier}$"
x_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
y_units = r""   # dimensionless (x opt_IRR baseline of each family member)
x_ticks = [0, 1, 2, 3, 4, 5, 6]
y_ticks = [0.2, 0.6, 1.0, 1.4, 1.8, 2.0]

# IRR in percent on a HARD 0-25 % scale (matches the sweep's IRR branch):
#  - money-losing cells (< 0 %, incl. -inf) render with the grey under-color;
#  - no over-color (grid max ~23 % < 25 %), so extend only the min side;
#  - 0 % break-even is a WHITE labeled contour line (comparison_lines);
#  - 5/10/15/20 % are black labeled contour lines (w_ticks).
IRR_w_levels = np.arange(0.0, 25.0001, 25.0 / 80)
IRR_cbar_ticks = np.arange(0.0, 25.0001, 5.0)          # 0,5,...,25
IRR_w_ticks = [5.0, 10.0, 15.0, 20.0]                  # black labeled lines
IRR_comparison_lines = [0.0]                           # white 0 % break-even
fmt_percent = lambda v, pos=None: f'{v:.0f}%'

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
    (inhib_isobutanol multiplier) values and columns spec_1 (k_13) values, as
    the sweep saves them."""
    path = os.path.join(RESULTS_DIR, f'{SWEEP_PREFIX}_{metric}.csv')
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f'{path}\nRun analyses/evaluate_EtOH_k13_inhib_isobutanol.py first '
            f'(opt_IRR, steps {STEPS}).')
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

    # IRR fraction -> percent, then push -inf (unsolvable, money-losing cells)
    # to just below the lowest filled-contour level so contourf fills them with
    # the grey under-color instead of masking them blank. Finite negative IRRs
    # are left as-is: extend_cmap='min' already colors anything below level 0.
    # This also keeps the IRR greedy climb well-behaved (no inf arithmetic).
    irr = 100. * results['IRR']
    step = IRR_w_levels[1] - IRR_w_levels[0]
    irr[np.isneginf(irr)] = IRR_w_levels[0] - step
    results['IRR'] = irr

    fig, ax, traj = tc.plot_metric_with_trajectories(
        results, SPEC_1, SPEC_2,
        color_metric=COLOR_METRIC,
        color_metric_label=r"$\mathbf{Internal\ Rate\ of\ Return}$",
        color_metric_units='%',
        baseline_point=(BASELINE_K13, BASELINE_MULT),
        trajectory_metrics=TRAJECTORY_METRICS,
        senses=SENSES,
        trajectory_colors=TRAJECTORY_COLORS,
        trajectory_linestyles=TRAJECTORY_LINESTYLES,
        trajectory_linewidths=TRAJECTORY_LINEWIDTHS,
        trajectory_markers=TRAJECTORY_MARKERS,
        trajectory_marker_sizes=TRAJECTORY_MARKER_SIZES,
        baseline_marker=BASELINE_MARKER,
        # a step must beat the sweep's own convergence tolerance (sim_rtol =
        # 1e-4) so a climb does not wander across a round-off plateau
        min_rel_improvement=1e-4,
        # the legend is built below (it must include the opt_IRR marker, which
        # the helper does not know about), so suppress the helper's own legend
        show_legend=False,
        x_label=x_label, y_label=y_label,
        x_units=x_units, y_units=y_units,
        x_ticks=x_ticks, y_ticks=y_ticks,
        cmap=tc.JBEI_UCB_colormap(reverse=True),
        w_levels=IRR_w_levels, cbar_ticks=IRR_cbar_ticks, w_ticks=IRR_w_ticks,
        extend_cmap='min',
        cmap_under_color=colors.grey_dark.shade(40).RGBn,
        # passed through to contourplots.animated_contourplot
        comparison_lines=IRR_comparison_lines,       # white 0 % break-even line
        comparison_lines_colors='white',
        fmt_clabel=fmt_percent,
        axis_title_fonts=axis_title_fonts,
        clabel_fontsize=clabel_fontsize,
        default_fontsize=default_fontsize,
        axis_tick_fontsize=axis_tick_fontsize,
        n_minor_ticks=1,
        cbar_n_minor_ticks=4,
        round_yticks_to=1,
        units_on_newline=(False, False, False, False),
        # x, y, z, w: bracket only the x units (y/w are dimensionless / percent)
        units_opening_brackets=[" (", "", "", ""],
        units_closing_brackets=[")", "", "", ""],
    )
    # percent ticks on the colorbar (the last axes contourplots added)
    cbar_ax = [a for a in fig.axes if a is not ax][-1]
    cbar_ax.yaxis.set_major_formatter(FuncFormatter(fmt_percent))

    # opt_IRR point (the financial optimum) at its true data coordinates, above
    # the trajectory lines/markers
    om_shape, om_color, om_size = OPT_IRR_MARKER
    ax.plot(OPT_IRR_K13, OPT_IRR_MULT, linestyle='None', marker=om_shape,
            markerfacecolor=om_color, markeredgecolor='k', markeredgewidth=0.8,
            markersize=om_size, zorder=700, clip_on=False)

    # Legend below the axes (8 handles: baseline + opt_IRR + 6 metrics) in two
    # rows; the light-grey face keeps the white baseline marker visible.
    b_shape, b_color, b_size = BASELINE_MARKER
    legend_handles = [
        Line2D([0], [0], color='none', linestyle='None', marker=b_shape,
               markerfacecolor=b_color, markeredgecolor='k',
               markersize=b_size, label='baseline (ethanol-only)'),
        Line2D([0], [0], color='none', linestyle='None', marker=om_shape,
               markerfacecolor=om_color, markeredgecolor='k',
               markersize=om_size, label='opt_IRR optimum'),
    ]
    for m in TRAJECTORY_METRICS:
        c = TRAJECTORY_COLORS[m]
        legend_handles.append(Line2D(
            [0], [0], color=c, marker=TRAJECTORY_MARKERS[m], markerfacecolor=c,
            markeredgecolor='k', markersize=TRAJECTORY_MARKER_SIZES[m],
            linestyle=TRAJECTORY_LINESTYLES[m],
            linewidth=TRAJECTORY_LINEWIDTHS[m], label=m))
    ax.legend(handles=legend_handles, loc='upper center',
              bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=8, framealpha=1.0,
              facecolor='#c8c8c8', edgecolor='k')

    stem = f'{COLOR_METRIC}_greedy_trajectories_{SWEEP_PREFIX}'
    png = os.path.join(RESULTS_DIR, stem + '.png')
    fig.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
    print('Saved:', png)

    # trajectory points, one row per step, for the record
    rows = []
    for m, d in traj.items():
        grid = results[m][0]
        for st, ((x, y), (iy, ix)) in enumerate(zip(d['path_xy'], d['path_ij'])):
            rows.append({'metric': m, 'sense': SENSES[m], 'step': st,
                         'k_13': x, 'inhib_isobutanol_multiplier': y,
                         'value': grid[iy, ix]})
    csv = os.path.join(RESULTS_DIR, stem + '.csv')
    pd.DataFrame(rows).to_csv(csv, index=False)
    print('Saved:', csv)

    ix0 = tc._snap_index(SPEC_1, BASELINE_K13)
    iy0 = tc._snap_index(SPEC_2, BASELINE_MULT)
    print(f'\nBaseline (k_13, inhib_isobutanol mult) = ({BASELINE_K13}, '
          f'{BASELINE_MULT}) snapped to grid cell '
          f'({SPEC_1[ix0]:.4g}, {SPEC_2[iy0]:.4g}); IRR values in %')
    for m in TRAJECTORY_METRICS:
        print(f'  {m:18s} {SENSES[m]}: {results[m][0][iy0, ix0]:.4g} at baseline')
    print()
    for m, d in traj.items():
        print(f'{m:18s} ({SENSES[m]}): {d["n_steps"]} steps -> optimum '
              f'{d["optimum_value"]:.4g} at (k_13, mult) = '
              f'({d["optimum_xy"][0]:.4g}, {d["optimum_xy"][1]:.4g})')
    return fig, ax, traj


if __name__ == '__main__':
    main()

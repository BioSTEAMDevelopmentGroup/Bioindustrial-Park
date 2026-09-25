#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Two-panel feeding-strategy figure for scenario A.

A  Minimum ethanol selling price (MESP = the sweep's purity-adjusted ethanol
   MPSP at a 15 % IRR, converted from $/kg to $/GGE exactly as in
   plots/plot_uncertainty_MPSP_vs_TCI.py) over the threshold x target
   glucose-concentration sweep of analyses/evaluate_feeding_strategies.py
   (spike cap optimized for MPSP at every grid point), with the grid optimum
   of cell density, ethanol titer, productivity, yield, TCI and MESP marked
   and annotated, and the white line = lowest target at which the optimized
   spike cap is zero (batch above, fed-batch below).
B  Concentration vs time of the scenario-A baseline fermentation
   (scenarios.SCENARIOS['A'] feeding strategy), with the ends of aeration and
   of fermentation marked.

Two stages:
  1. (simulation) one isobutanol.load() + scenarios.load_scenario('A'); the
     V406 time course is cached to
     analyses/results/feeding_strategy_baseline_A_trajectory.{csv,json}.
     Runs only when that cache is missing or with --resimulate.
  2. (sim-safe) reads the sweep CSVs + the cached trajectory and draws.
     With a cache present the package is never imported, so it is safe
     alongside any running simulation.

Output: analyses/results/publication/Feed-strat/feeding_strategy_baseline_A.{png,pdf}
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.lines import TICKDOWN, TICKLEFT
from matplotlib.patches import Polygon
from matplotlib.ticker import AutoMinorLocator, FixedLocator, MultipleLocator
import matplotlib.patheffects as pe

#%% Paths and settings

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PACKAGE_DIR, 'analyses', 'results')
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Feed-strat')
OUTPUT_STEM = 'feeding_strategy_baseline_A'
TRAJECTORY_STEM = os.path.join(RESULTS_DIR, 'feeding_strategy_baseline_A_trajectory')

SCENARIO = 'A'
# Grid of the sweep to read (evaluate_feeding_strategies.py, IBO_SWEEP_STEPS);
# its spec_1 / spec_2 below must match that script's linspaces.
SWEEP_STEPS = (40, 40, 1)
SPEC_1 = np.linspace(1., 400., SWEEP_STEPS[0])   # threshold glucose conc. [g/L]
SPEC_2 = np.linspace(10., 400., SWEEP_STEPS[1])  # target glucose conc. [g/L]
SWEEP_CSV = os.path.join(
    RESULTS_DIR,
    f'ibo_{SWEEP_STEPS}_Thres_Targe_Max n_{SCENARIO}__{{metric}}.csv')

# MESP in $ per gasoline gallon equivalent, with the conversion of
# plots/plot_uncertainty_MPSP_vs_TCI.py: $/GGE = $/kg x KG_PER_GAL /
# GGE_PER_GAL; ethanol at 0.789 kg/L (20 C) x 3.785411784 L/gal = 2.987 kg/gal,
# 1 gal ethanol = 0.67 GGE (AFDC), so ~4.458 ($/GGE)/($/kg)
ETHANOL_DENSITY_KG_PER_L = 0.789
L_PER_GAL = 3.785411784
KG_PER_GAL = ETHANOL_DENSITY_KG_PER_L * L_PER_GAL
GGE_PER_GAL = 0.67
USD_PER_KG_TO_USD_PER_GGE = KG_PER_GAL / GGE_PER_GAL

# Panel A colour scale (MESP, $/GGE): the sweep spans 3.82-5.03
MESP_LEVELS = np.arange(3.75, 5.25001, 0.025)
MESP_CBAR_TICKS = np.arange(3.75, 5.25001, 0.25)
MESP_CBAR_MINOR_STEP = 0.05

# Optimum markers: (sweep metric, 'min'/'max', label, marker, face colour,
# size [pt], label offset from the marker [pt], arrow curvature). Offsets are
# set by eye for the current sweep; retune them if the optima move. A label is
# left-/right-aligned by the sign of its x offset unless LABEL_HA overrides it.
# Aeration, AOC and slurry evaporation duty are deliberately not marked.
OPTIMA = [
    ('EtOH Titer',        'max', 'titer',        '^', 'white',   10, (11, -15),  0.3),
    ('Cell loading',      'max', 'cell density', 'o', 'white',   10, (-6, 16),  -0.3),
    ('EtOH Productivity', 'max', 'productivity', 's', 'white',    9, (8, 17),   -0.2),
    ('EtOH Yield',        'max', 'yield',        'p', 'white',   10, (12, -14),  0.3),
    ('TCI',               'min', 'TCI',          'p', '#33ccff', 10, (0, -22),   0.3),
    ('MPSP',              'min', 'MESP',         '*', '#33ccff', 14, (12, -14),  0.3),
]
# productivity: centred above its optimum, clear of the batch label to its
# right; TCI: right-aligned just right of its optimum, so the label sits
# wholly in the white (threshold > target) region
LABEL_HA = {'productivity': 'center', 'TCI': 'right'}
# label text colour (default black); titer sits on the dark high-MESP band
LABEL_COLOR = {'titer': 'white'}

# Both batch / fed-batch labels: anchored on the boundary line at threshold
# 150 g/L, then shifted by this many points (x, y)
BATCH_LABEL_SHIFT = (-8., -5.)

# Panel B series: (nskinetics results key, legend label); default colour cycle
TRAJECTORY_SERIES = [
    ('[x]', 'cell density'),
    ('curr_a', 'active cell density'),
    ('[s_glu]', 'glucose'),
    ('[s_EtOH]', 'ethanol'),
    ('[s_acetate]', 'acetate'),
    ('[s_IBO]', 'isobutanol'),
]
TIME_LIM = (0., 80.)
CONC_LIM = (0., 250.)

#%% Formatting (see plots/across_feed_strat_multipanel_global_with_optima_summary.py)

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 13, 'annotation': 12, 'legend': 12,
         'panel_letter': 20}
TICK_LEN = {'major': 4.0, 'minor': 2.0}  # pt, each way for left/bottom
G_PER_L = r'$\mathrm{g·L}^{-1}$'
USD_PER_GGE = r'$\mathrm{\$·GGE}^{-1}$'


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'
    plt.rcParams['pdf.fonttype'] = 42


def bold_title(name, units):
    return r'$\bf{' + name.replace(' ', r'\ ') + '}$' + f' [{units}]'


def style_ticks(ax):
    # all four sides; top/right inward, left/bottom in and out equally.
    # Call after fig.canvas.draw() so every tick object exists.
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


def JBEI_UCB_colormap(N_levels=90):
    # low = UCB yellow -> JBEI orange -> UCB blue -> dark grey = high
    try:
        from biosteam.utils import colors
        grey_dark = colors.grey_dark.RGBn
    except Exception:  # keep stage 2 importable without biosteam
        grey_dark = (0.28, 0.28, 0.28)
    cmap_colors = [(253/255, 181/255, 21/255), (233/255, 83/255, 39/255),
                   (0/255, 38/255, 118/255), grey_dark]
    return LinearSegmentedColormap.from_list('JBEI_UCB', cmap_colors, N_levels)

#%% Stage 1: baseline fermentation time course (simulation)

def simulate_baseline_trajectory():
    """One load + the scenario-A baseline; caches V406's time course."""
    from biorefineries import isobutanol
    isobutanol.load()
    from biorefineries.isobutanol import scenarios
    bundle = scenarios.load_scenario(SCENARIO)  # runs one baseline simulation
    results = bundle['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'))
    V406 = bundle['V406']
    d = V406.nsk_results_dict
    df = pd.DataFrame({'time': np.asarray(d['time'], dtype=float)})
    for key, _ in TRAJECTORY_SERIES:
        df[key] = np.asarray(d[key], dtype=float)
    df.to_csv(TRAJECTORY_STEM + '.csv', index=False)
    specific = V406.nsk_results_specific_tau_dict
    meta = {
        'scenario': SCENARIO,
        'feeding_kwargs': bundle['feeding_kwargs'],
        'max_n_spikes': int(bundle['fbs_spec'].max_n_spikes),
        'n_glu_spikes': float(specific['curr_n_glu_spikes']),
        'tau': float(V406.tau),
        'tau_stop_aeration': float(V406.tau_stop_aeration),
        'EtOH_titer': float(specific['[s_EtOH]']),
        'cell_density': float(specific['[x]']),
        'ethanol_MPSP': float(results['MPSPs']['ethanol']),
        'IRR': float(results['IRR']),
    }
    with open(TRAJECTORY_STEM + '.json', 'w') as file:
        json.dump(meta, file, indent=2)
    print(f'Cached the baseline trajectory: {TRAJECTORY_STEM}.csv/.json')
    return df, meta


def load_baseline_trajectory(resimulate=False):
    if resimulate or not (os.path.exists(TRAJECTORY_STEM + '.csv')
                          and os.path.exists(TRAJECTORY_STEM + '.json')):
        return simulate_baseline_trajectory()
    with open(TRAJECTORY_STEM + '.json') as file:
        meta = json.load(file)
    print(f'Using the cached baseline trajectory {TRAJECTORY_STEM}.csv '
          '(--resimulate to refresh).')
    return pd.read_csv(TRAJECTORY_STEM + '.csv'), meta

#%% Stage 2: sweep data

def load_sweep_metric(metric):
    """(n_target, n_threshold) array: row = target (SPEC_2), column =
    threshold (SPEC_1); NaN where threshold > target (not simulated)."""
    path = SWEEP_CSV.format(metric=metric)
    arr = pd.read_csv(path).iloc[:, 1:].to_numpy(dtype=float)
    if arr.shape != (len(SPEC_2), len(SPEC_1)):
        raise ValueError(f'{path}: shape {arr.shape} does not match the '
                         f'{SWEEP_STEPS} grid.')
    return arr


def grid_optimum(arr, sense):
    finite = np.isfinite(arr)
    value = np.nanmin(arr[finite]) if sense == 'min' else np.nanmax(arr[finite])
    i, j = np.argwhere(arr == value)[0]
    return SPEC_1[j], SPEC_2[i], value


def fill_past_diagonal(arr):
    """Carry each row's last finite value across the unsimulated cells right
    of it (threshold > target), so the filled contours reach the exact
    threshold = target edge that clips them. Interior NaNs are reported."""
    out = arr.copy()
    for row in out:
        idx = np.flatnonzero(np.isfinite(row))
        if not idx.size:
            continue
        interior = ~np.isfinite(row[:idx[-1]+1])
        if interior.any():
            print(f'  WARNING: {interior.sum()} interior NaN cell(s) in a row '
                  '(failed sweep points) -- filled from the left.')
            for k in np.flatnonzero(interior):
                row[k] = row[k-1] if k else row[idx[0]]
        row[idx[-1]+1:] = row[idx[-1]]
    return out


def batch_boundary(n_spikes):
    """For each threshold: the lowest target whose MPSP-optimized spike cap
    is zero (batch above). Truncated where it meets the threshold = target
    edge, since from there on it would run along that edge."""
    xs, ys = [], []
    for j, x in enumerate(SPEC_1):
        col = n_spikes[:, j]
        rows = np.flatnonzero(np.isfinite(col) & (col == 0))
        if not rows.size:
            continue
        y = SPEC_2[rows[0]]
        first_row = np.flatnonzero(np.isfinite(col))[0]
        xs.append(x); ys.append(y)
        if rows[0] == first_row:
            break
    return np.array(xs), np.array(ys)

#%% Figure

def draw_panel_A(fig, ax, cax):
    mesp = load_sweep_metric('MPSP') * USD_PER_KG_TO_USD_PER_GGE
    filled = fill_past_diagonal(mesp)
    # pad a threshold = 0 column so the fill reaches the y axis
    x = np.concatenate([[0.], SPEC_1])
    z = np.column_stack([filled[:, 0], filled])
    cmap = JBEI_UCB_colormap()
    norm = BoundaryNorm(MESP_LEVELS, cmap.N)
    cs = ax.contourf(x, SPEC_2, z, levels=MESP_LEVELS, cmap=cmap, norm=norm,
                     zorder=1)
    cs.set_edgecolor('face')  # no hairline seams between bands in the PDF
    domain = Polygon([(0., SPEC_2[0]), (0., 400.), (400., 400.),
                      (SPEC_2[0], SPEC_2[0])], closed=True,
                     transform=ax.transData, facecolor='none', edgecolor='none')
    ax.add_patch(domain)
    cs.set_clip_path(domain)

    ax.set_xlim(0., 400.)
    ax.set_ylim(0., 400.)

    # batch / fed-batch boundary, labelled on either side along its local
    # slope (rotation in screen space; axes position and limits are final)
    bx, by = batch_boundary(load_sweep_metric('Number of glucose spikes'))
    ax.plot(bx, by, color='white', linewidth=0.8, alpha=0.85, zorder=3,
            clip_path=domain)
    k = int(np.argmin(np.abs(bx - 150.)))
    k0, k1 = max(k - 5, 0), min(k + 5, len(bx) - 1)
    slope = np.polyfit(bx[k0:k1+1], by[k0:k1+1], 1)[0]
    p0 = ax.transData.transform((bx[k], by[k]))
    p1 = ax.transData.transform((bx[k] + 1., by[k] + slope))
    angle = np.arctan2(p1[1] - p0[1], p1[0] - p0[0])
    normal = np.array([-np.sin(angle), np.cos(angle)])
    # (gap from the line in points: the batch label sits over the steps)
    for text, gap in (('↑ batch', 12.), ('↓ fed-batch', -9.)):
        ax.annotate(text, xy=(bx[k], by[k]),
                    xytext=gap * normal + np.array(BATCH_LABEL_SHIFT),
                    textcoords='offset points', ha='center', va='center',
                    rotation=np.degrees(angle), rotation_mode='anchor',
                    color='white', fontsize=FONTS['annotation'], zorder=4)

    # optimum markers + labels
    optima = {}
    for metric, sense, label, marker, color, size, offset, rad in OPTIMA:
        ox, oy, value = grid_optimum(load_sweep_metric(metric), sense)
        if metric == 'MPSP':
            value *= USD_PER_KG_TO_USD_PER_GGE  # MESP, $/GGE
        optima[label] = (ox, oy, value)
        ax.plot(ox, oy, linestyle='none', marker=marker, markersize=size,
                markerfacecolor=color, markeredgecolor='black',
                markeredgewidth=0.8, zorder=10, clip_on=False)
        ax.annotate(label, xy=(ox, oy), xytext=offset,
                    textcoords='offset points', fontsize=FONTS['annotation'],
                    color=LABEL_COLOR.get(label, 'black'),
                    ha=LABEL_HA.get(label, 'left' if offset[0] >= 0 else 'right'),
                    va='bottom' if offset[1] >= 0 else 'top', zorder=11,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=9, color='black', lw=0.9,
                                    shrinkA=1, shrinkB=size/2 + 1,
                                    connectionstyle=f'arc3,rad={rad}'))

    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(MultipleLocator(100.))
        axis.set_minor_locator(AutoMinorLocator(2))
    ax.set_xlabel(bold_title('Threshold glucose concentration', G_PER_L),
                  fontsize=FONTS['axis_title'])
    ax.set_ylabel(bold_title('Target glucose concentration', G_PER_L),
                  fontsize=FONTS['axis_title'])

    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cax, spacing='proportional')
    cbar.set_ticks(MESP_CBAR_TICKS)
    cbar.set_ticklabels([f'{t:.2f}' for t in MESP_CBAR_TICKS])
    cbar.ax.yaxis.set_minor_locator(FixedLocator(
        [v for v in np.arange(MESP_LEVELS[0], MESP_LEVELS[-1] + 1e-9,
                              MESP_CBAR_MINOR_STEP)
         if not np.any(np.isclose(v, MESP_CBAR_TICKS))]))
    cbar.ax.tick_params(which='major', labelsize=FONTS['tick'],
                        length=TICK_LEN['major'])
    cbar.ax.tick_params(which='minor', length=TICK_LEN['minor'])
    cbar.set_label(bold_title('MESP', USD_PER_GGE), fontsize=FONTS['axis_title'])
    return optima


def draw_panel_B(fig, ax, df, meta):
    for key, label in TRAJECTORY_SERIES:
        ax.plot(df['time'], df[key], label=label, linewidth=1.6)
    ax.set_xlim(*TIME_LIM)
    ax.set_ylim(*CONC_LIM)
    ax.xaxis.set_major_locator(MultipleLocator(20.))
    ax.yaxis.set_major_locator(MultipleLocator(50.))
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_minor_locator(AutoMinorLocator(4))
    ax.set_xlabel(bold_title('Time', 'h'), fontsize=FONTS['axis_title'])
    ax.set_ylabel(bold_title('Concentration', G_PER_L),
                  fontsize=FONTS['axis_title'])
    ax.legend(loc='upper right', fontsize=FONTS['legend'], frameon=True,
              fancybox=False, edgecolor='black', framealpha=1.,
              handlelength=2.2, borderaxespad=0.6)

    for t, text in ((meta['tau_stop_aeration'], 'end of aeration'),
                    (meta['tau'], 'end of fermentation')):
        ax.axvline(t, color='grey', linestyle='--', linewidth=1.0, zorder=0)
        ax.annotate(text, xy=(t, CONC_LIM[1]), xytext=(4, 12),
                    textcoords='offset points', ha='left', va='bottom',
                    fontsize=FONTS['annotation'], annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=9, color='black', lw=0.9,
                                    shrinkA=1, shrinkB=1,
                                    connectionstyle='arc3,rad=0.4'))


def make_figure(df, meta):
    apply_font_rcparams()
    fig = plt.figure(figsize=(11.0, 4.2))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[1.0, 0.045, 0.42, 0.86],
                  wspace=0.05, left=0.08, right=0.98, bottom=0.14, top=0.86)
    axA = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    axB = fig.add_subplot(gs[0, 3])
    optima = draw_panel_A(fig, axA, cax)
    draw_panel_B(fig, axB, df, meta)
    for ax, letter in ((axA, 'A'), (axB, 'B')):
        ax.text(-0.2, 1.07, letter, transform=ax.transAxes, ha='left',
                va='bottom', fontsize=FONTS['panel_letter'], fontweight='bold')
    fig.canvas.draw()
    style_ticks(axA)
    style_ticks(axB)
    return fig, optima


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--resimulate', action='store_true',
                        help='re-run the scenario-A baseline simulation even '
                             'if the trajectory cache exists')
    args = parser.parse_args()

    df, meta = load_baseline_trajectory(args.resimulate)
    fig, optima = make_figure(df, meta)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUTPUT_DIR, f'{OUTPUT_STEM}.{ext}')
        fig.savefig(path, dpi=600, bbox_inches='tight', facecolor='white')
        print(f'Saved {path}')
    plt.close(fig)

    print('\nGrid optima (threshold, target [g/L]: value):')
    for label, (ox, oy, value) in optima.items():
        print(f'  {label:14s} {ox:6.1f}, {oy:6.1f}: {value:.5g}')
    print(f"\nBaseline: feeding {meta['feeding_kwargs']}, "
          f"{meta['n_glu_spikes']:g} spikes, tau {meta['tau']:.2f} h, end of "
          f"aeration {meta['tau_stop_aeration']:.2f} h, ethanol titer "
          f"{meta['EtOH_titer']:.1f} g/L, MESP "
          f"{meta['ethanol_MPSP'] * USD_PER_KG_TO_USD_PER_GGE:.4f} $/GGE "
          f"({meta['ethanol_MPSP']:.5f} $/kg), IRR {meta['IRR']:.4f}")


if __name__ == '__main__':
    main()

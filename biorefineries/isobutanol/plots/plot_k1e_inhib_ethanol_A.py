#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""k_1e x ethanol-inhibition sweep figure for scenario A.

Minimum ethanol selling price (MESP = the sweep's purity-adjusted ethanol MPSP
at a 15 % IRR, converted from $/kg to $/GGE exactly as in
plots/plot_uncertainty_MPSP_vs_TCI.py) over the k_1e x inhib_ethanol-multiplier
grid of analyses/evaluate_EtOH_k1e_inhib_ethanol.py (scenario A, enzyme burden
ON, A's fed-batch feeding strategy, no feeding-strategy optimization), with the
grid optimum of cell density, ethanol titer, productivity, yield, TCI and MESP
marked and annotated. Format and style follow
plots/plot_feeding_strategy_baseline_A.py (panel A), whose helpers are reused
by file path.

Sim-safe: reads the sweep CSVs only and never imports the package.

Output: analyses/results/publication/Kinetic-sweeps/k1e_inhib_ethanol_A.{png,pdf}
"""

import importlib.util
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator, FixedLocator, MultipleLocator
from matplotlib.transforms import ScaledTranslation

#%% Shared style (plots/plot_feeding_strategy_baseline_A.py, by file path)

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_feeding_strategy_style', os.path.join(HERE, 'plot_feeding_strategy_baseline_A.py'))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

FONTS = fs.FONTS
TICK_LEN = fs.TICK_LEN
USD_PER_KG_TO_USD_PER_GGE = fs.USD_PER_KG_TO_USD_PER_GGE

#%% Paths and settings

PACKAGE_DIR = os.path.dirname(HERE)
RESULTS_DIR = os.path.join(PACKAGE_DIR, 'analyses', 'results')
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Kinetic-sweeps')
OUTPUT_STEM = 'k1e_inhib_ethanol_A'

# Grid of the sweep to read (analyses/evaluate_EtOH_k1e_inhib_ethanol.py);
# SPEC_1 / SPEC_2 must match that script's linspaces: k_1e over the
# metabolic_14d glycolysis band 0.2x-4x of scenario A's k_1e (47.1 g/L/h, the
# nskinetics antimony value), the inhib_ethanol multiplier over the default
# group band 0.75x-1.5x.
SWEEP_STEPS = (40, 40, 1)
BASELINE_K1E = 47.1
SPEC_1 = np.linspace(0.2 * BASELINE_K1E, 4.0 * BASELINE_K1E, SWEEP_STEPS[0])  # k_1e [g/L/h]
SPEC_2 = np.linspace(0.75, 1.5, SWEEP_STEPS[1])  # inhib_ethanol multiplier [-]
SWEEP_CSV = os.path.join(
    RESULTS_DIR,
    f'ibo_{SWEEP_STEPS}_k_1e_inhib_Spike_opt=False_max_n=16__{{metric}}.csv')

# Colour scale (MESP, $/GGE): the same levels as the feeding-strategy figure;
# the sweep spans 3.79-8.48, so the high-k_1e wall past 5.25 takes the
# over-colour (extend='max')
MESP_LEVELS = fs.MESP_LEVELS
MESP_CBAR_TICKS = fs.MESP_CBAR_TICKS
MESP_CBAR_MINOR_STEP = fs.MESP_CBAR_MINOR_STEP

# Optimum markers: (sweep metric, 'min'/'max', label, marker, face colour,
# size [pt], label offset from the marker [pt], arrow curvature), as in the
# feeding-strategy figure. The sweep has no separate ethanol-yield column:
# 'Combined Yield' (ethanol + isobutanol per g sugars added) is the ethanol
# yield here, since scenario A makes no isobutanol. Offsets are set by eye for
# the current sweep; retune them if the optima move. A label is left-/right-
# aligned by the sign of its x offset unless LABEL_HA overrides it.
OPTIMA = [
    ('EtOH Titer',        'max', 'titer',        '^', 'white',   10, (-10, 22),  0.3),
    ('Cell loading',      'max', 'cell density', 'o', 'white',   10, (6, 32),   -0.3),
    ('EtOH Productivity', 'max', 'productivity', 's', 'white',    9, (0, 22),   -0.2),
    ('Combined Yield',    'max', 'yield',        'p', 'white',   10, (-14, 16),  0.3),
    ('TCI',               'min', 'TCI',          'p', '#33ccff', 10, (24, 6),    0.3),
    ('MPSP',              'min', 'MESP',         '*', '#33ccff', 14, (10, 34),  -0.3),
]
LABEL_HA = {'productivity': 'center'}
LABEL_COLOR = {}
# Co-located optima are fanned apart by this many points (x, y): cell density
# and TCI both sit at the lowest-k_1e, lowest-multiplier corner
MARKER_NUDGE = {'cell density': (0., 5.), 'TCI': (5., 0.)}

X_MAJOR, X_MINOR_DIV = 50., 2
Y_MAJOR, Y_MINOR_DIV = 0.25, 5
G_PER_L_PER_H = r'$\mathrm{g·L}^{-1}\mathrm{·h}^{-1}$'

#%% Sweep data

def load_sweep_metric(metric):
    """(n_multiplier, n_k_1e) array: row = inhib_ethanol multiplier (SPEC_2),
    column = k_1e (SPEC_1)."""
    path = SWEEP_CSV.format(metric=metric)
    arr = pd.read_csv(path).iloc[:, 1:].to_numpy(dtype=float)
    if arr.shape != (len(SPEC_2), len(SPEC_1)):
        raise ValueError(f'{path}: shape {arr.shape} does not match the '
                         f'{SWEEP_STEPS} grid.')
    return arr


def grid_optimum(arr, sense):
    finite = np.isfinite(arr)
    value = np.min(arr[finite]) if sense == 'min' else np.max(arr[finite])
    i, j = np.argwhere(arr == value)[0]
    return SPEC_1[j], SPEC_2[i], value


def fill_failed_cells(arr):
    """Replace each non-finite cell (failed sweep point) by the mean of its
    finite 4-neighbours, so the filled contours have no hole."""
    out = arr.copy()
    bad = np.argwhere(~np.isfinite(arr))
    for i, j in bad:
        neighbours = [arr[a, b] for a, b in ((i-1, j), (i+1, j), (i, j-1), (i, j+1))
                      if 0 <= a < arr.shape[0] and 0 <= b < arr.shape[1]
                      and np.isfinite(arr[a, b])]
        out[i, j] = np.mean(neighbours)
        print(f'  WARNING: failed sweep point at k_1e = {SPEC_1[j]:.1f}, '
              f'multiplier = {SPEC_2[i]:.3f} -- filled from its neighbours.')
    return out

#%% Figure

def draw_panel(fig, ax, cax):
    mesp = fill_failed_cells(load_sweep_metric('MPSP') * USD_PER_KG_TO_USD_PER_GGE)
    cmap = fs.JBEI_UCB_colormap()
    norm = BoundaryNorm(MESP_LEVELS, cmap.N, extend='max')
    cs = ax.contourf(SPEC_1, SPEC_2, mesp, levels=MESP_LEVELS, cmap=cmap,
                     norm=norm, extend='max', zorder=1)
    cs.set_edgecolor('face')  # no hairline seams between bands in the PDF
    ax.set_xlim(SPEC_1[0], SPEC_1[-1])
    ax.set_ylim(SPEC_2[0], SPEC_2[-1])

    # optimum markers + labels
    optima = {}
    for metric, sense, label, marker, color, size, offset, rad in OPTIMA:
        ox, oy, value = grid_optimum(load_sweep_metric(metric), sense)
        if metric == 'MPSP':
            value *= USD_PER_KG_TO_USD_PER_GGE  # MESP, $/GGE
        optima[label] = (ox, oy, value)
        dx, dy = MARKER_NUDGE.get(label, (0., 0.))
        where = ax.transData + ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
        ax.plot(ox, oy, linestyle='none', marker=marker, markersize=size,
                markerfacecolor=color, markeredgecolor='black',
                markeredgewidth=0.8, zorder=10, clip_on=False, transform=where)
        ax.annotate(label, xy=(ox, oy), xycoords=where, xytext=offset,
                    textcoords='offset points', fontsize=FONTS['annotation'],
                    color=LABEL_COLOR.get(label, 'black'),
                    ha=LABEL_HA.get(label, 'left' if offset[0] >= 0 else 'right'),
                    va='bottom' if offset[1] >= 0 else 'top', zorder=11,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=9, color='black', lw=0.9,
                                    shrinkA=1, shrinkB=size/2 + 1,
                                    connectionstyle=f'arc3,rad={rad}'))

    ax.xaxis.set_major_locator(MultipleLocator(X_MAJOR))
    ax.xaxis.set_minor_locator(AutoMinorLocator(X_MINOR_DIV))
    ax.yaxis.set_major_locator(MultipleLocator(Y_MAJOR))
    ax.yaxis.set_minor_locator(AutoMinorLocator(Y_MINOR_DIV))
    ax.set_xlabel(r'$\bf{Glycolytic\ capacity}$, $k_{\mathrm{1e}}$' + f' [{G_PER_L_PER_H}]',
                  fontsize=FONTS['axis_title'])
    ax.set_ylabel(r'$\bf{Ethanol\ inhibition\ multiplier}$',
                  fontsize=FONTS['axis_title'])

    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cax, spacing='proportional', extend='max',
                        extendfrac=0.04)
    cbar.set_ticks(MESP_CBAR_TICKS)
    cbar.set_ticklabels([f'{t:.2f}' for t in MESP_CBAR_TICKS])
    cbar.ax.yaxis.set_minor_locator(FixedLocator(
        [v for v in np.arange(MESP_LEVELS[0], MESP_LEVELS[-1] + 1e-9,
                              MESP_CBAR_MINOR_STEP)
         if not np.any(np.isclose(v, MESP_CBAR_TICKS))]))
    cbar.ax.tick_params(which='major', labelsize=FONTS['tick'],
                        length=TICK_LEN['major'])
    cbar.ax.tick_params(which='minor', length=TICK_LEN['minor'])
    cbar.set_label(fs.bold_title('MESP', fs.USD_PER_GGE), fontsize=FONTS['axis_title'])
    return optima


def make_figure():
    fs.apply_font_rcparams()
    fig = plt.figure(figsize=(5.6, 4.2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 0.045], wspace=0.05,
                  left=0.15, right=0.86, bottom=0.14, top=0.86)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    optima = draw_panel(fig, ax, cax)
    fig.canvas.draw()
    fs.style_ticks(ax)
    return fig, optima


def main():
    fig, optima = make_figure()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUTPUT_DIR, f'{OUTPUT_STEM}.{ext}')
        fig.savefig(path, dpi=600, bbox_inches='tight', facecolor='white')
        print(f'Saved {path}')
    plt.close(fig)

    print('\nGrid optima (k_1e [g/L/h], inhib_ethanol multiplier: value):')
    for label, (ox, oy, value) in optima.items():
        print(f'  {label:14s} {ox:6.1f}, {oy:6.3f}: {value:.5g}')


if __name__ == '__main__':
    main()

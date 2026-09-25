#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Bivariate uncertainty plot of the purity-adjusted ethanol MPSP (y) against
the total capital investment (x), from an uncertainty-analysis results
workbook (*_1_full_evaluation.xlsx as written by
analyses/full/uncertainties_IBO_EtOH.py).

Joint panel: every Monte Carlo sample as a teal dot at 25 % opacity and the
baseline (the 'initial' row of the companion *_0_baseline.xlsx) as a white
diamond (unlabelled: name it in the caption), over a light grey band showing
the ethanol market price range (ETHANOL_MARKET_RANGE), with contour lines
of a Gaussian KDE enclosing 5 / 25 / 50 / 75 / 95 % of the samples (the
highest-density regions; each line is the density quantile AT the samples,
so it holds that share of them). Marginal box plots
outside the panel: box = 25th-75th percentile, line = median, whiskers = 5th-95th
percentile (the whis=[5, 95] of contourplots.box_and_whiskers_plot), dots =
1st and 99th percentiles.

Sim-safe: pure pandas/matplotlib/scipy, never imports biorefineries.

Usage:
    python plot_uncertainty_MPSP_vs_TCI.py [<results.xlsx>]
        [--baseline-file PATH] [--out-dir DIR] [--stem STEM] [--dpi N]

With no workbook, the newest scenario-A *_1_full_evaluation.xlsx under
analyses/results is used.
"""

import os
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.lines import TICKDOWN, TICKLEFT
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from scipy import stats

__all__ = ('plot_uncertainty_MPSP_vs_TCI',)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(os.path.dirname(HERE), 'analyses', 'results')
DEFAULT_OUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Uncertainty')

MPSP_COL = ('Biorefinery', 'Purity-adjusted ethanol MPSP [$/kg]')
TCI_COL = ('Biorefinery', 'Total capital investment [10^6 $]')

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'clabel': 10}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# the TRY-informed profitability campaign's teal (RELAY_COLOR in
# plots/plot_kin_opt_parameter_sets.py)
TEAL = '#0B6E7A'
SAMPLE_ALPHA = 0.25   # joint-panel samples: 75 % transparent
SAMPLE_SIZE = 6 # pt^2


def _mix(color, other, t):
    """`color` moved a fraction `t` of the way towards `other` (RGB)."""
    a, b = np.array(to_rgb(color)), np.array(to_rgb(other))
    return to_hex((1 - t)*a + t*b)


BOX_FACE = TEAL
BOX_EDGE = _mix(TEAL, 'black', 0.5)
BOX_MEDIAN = 'black'
# KDE contour lines: share of samples each encloses, drawn in dark teal
CONTOUR_SHARES = (0.05, 0.25, 0.50, 0.75, 0.95)
CONTOUR_COLOR = BOX_EDGE
CONTOUR_LW = 1.0
# on-screen bearing [deg] from the density peak at which each contour is
# labelled: the inner contours are too close together to label along one ray,
# so the labels alternate sides (and keep clear of the baseline marker)
CONTOUR_LABEL_BEARINGS = {0.05: 0, 0.25: 180, 0.50: 0, 0.75: 180, 0.95: 0}
INK = '#0b0b0b'

# ethanol market price range, the same one analyses/full/uncertainties_IBO_EtOH.py
# draws on its MPSP box plot: Jan 2021 - Dec 2025 five-year low and high,
# 1.5475 and 3.4500 $/gal / (3.7854 L/gal * 0.789 kg/L), from
# https://tradingeconomics.com/commodity/ethanol
ETHANOL_MARKET_RANGE = (0.52, 1.15) # $/kg
MPSP_AXIS_LIMITS = (0.0, 1.2) # $/kg
# a light shade of the baseline grey of plots/plot_kin_opt_parameter_sets.py
# (BASELINE_COLOR); the band is unlabelled, named in the caption
BASELINE_GRAY = '#90918e'
MARKET_BAND_COLOR = _mix(BASELINE_GRAY, 'white', 0.7)
BOX_PERCENTILES = {'whis': (5, 95), 'dots': (1, 99)}


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


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


#%% Data

def newest_results_file(scenario='A'):
    files = [f for f in glob.glob(os.path.join(glob.escape(RESULTS_DIR), '*_1_full_evaluation.xlsx'))
             if f"['{scenario}']" in os.path.basename(f)]
    if not files:
        raise FileNotFoundError(f'no scenario-{scenario} *_1_full_evaluation.xlsx in {RESULTS_DIR}')
    return max(files, key=os.path.getmtime)


def load_samples(results_file):
    df = pd.read_excel(results_file, sheet_name='TEA results',
                       header=[0, 1], index_col=0)
    xy = df[[TCI_COL, MPSP_COL]].astype(float)
    finite = np.isfinite(xy.values).all(axis=1)
    return xy.values[finite, 0], xy.values[finite, 1], int((~finite).sum())


def load_baseline(baseline_file):
    df = pd.read_excel(baseline_file, header=[0, 1], index_col=0)
    row = df.loc['initial']
    return float(row[TCI_COL]), float(row[MPSP_COL])


#%% Drawing

def nice_ticks(values, nbins=7):
    """'Nice' major ticks enclosing the data; the axis limits sit on the
    first and last one."""
    ticks = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(
        values.min(), values.max())
    return ticks[(ticks <= values.min()).nonzero()[0][-1]:
                 (ticks >= values.max()).nonzero()[0][0] + 1]


def draw_samples(ax, x, y):
    # every sample; rasterized so the PDF stays small (axes, text and the
    # baseline marker stay vector)
    ax.scatter(x, y, s=SAMPLE_SIZE, color=TEAL, alpha=SAMPLE_ALPHA, lw=0,
               rasterized=True, zorder=2)


def draw_hdr_contours(ax, x, y, xlim, ylim, n_grid=300, pad=0.15):
    """KDE contour lines enclosing CONTOUR_SHARES of the samples, labelled
    inline with the share, each where its contour crosses the ray from the
    density peak at that share's CONTOUR_LABEL_BEARINGS angle (on screen)."""
    kde = stats.gaussian_kde(np.vstack([x, y]))
    dx, dy = np.ptp(x), np.ptp(y)
    GX, GY = np.meshgrid(np.linspace(x.min() - pad*dx, x.max() + pad*dx, n_grid),
                         np.linspace(y.min() - pad*dy, y.max() + pad*dy, n_grid))
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    at_samples = kde(np.vstack([x, y]))
    # the region holding share p is {density >= the (1 - p) quantile of the
    # density at the samples}; larger share -> lower level
    levels = {np.quantile(at_samples, 1 - p): p for p in CONTOUR_SHARES}
    cs = ax.contour(GX, GY, Z, levels=sorted(levels), colors=CONTOUR_COLOR,
                    linewidths=CONTOUR_LW, zorder=3)

    # label positions: per level, the contour vertex whose on-screen bearing
    # from the peak is closest to the ray (axes size in inches converts the
    # data offsets to screen proportions)
    i_peak = np.unravel_index(np.argmax(Z), Z.shape)
    px, py = GX[i_peak], GY[i_peak]
    fig_w, fig_h = ax.figure.get_size_inches()
    box = ax.get_position()
    sx = box.width*fig_w/(xlim[1] - xlim[0])
    sy = box.height*fig_h/(ylim[1] - ylim[0])
    positions = []
    for lev, segs in zip(cs.levels, cs.allsegs):
        ray = np.deg2rad(CONTOUR_LABEL_BEARINGS[levels[lev]])
        verts = np.concatenate([seg for seg in segs if len(seg)])
        bearing = np.arctan2((verts[:, 1] - py)*sy, (verts[:, 0] - px)*sx)
        k = np.argmin(np.abs(np.angle(np.exp(1j*(bearing - ray)))))
        positions.append(tuple(verts[k]))
    labels = ax.clabel(cs, fmt={lev: f'{round(100*p)}%' for lev, p in levels.items()},
                       manual=positions, fontsize=FONTS['clabel'], inline=True,
                       inline_spacing=3)
    for t in labels:
        t.set_color(INK)


def draw_box(ax, values, orientation):
    lo_w, hi_w = BOX_PERCENTILES['whis']
    ax.boxplot(values, whis=[lo_w, hi_w], orientation=orientation,
               widths=0.6, showfliers=False, patch_artist=True,
               boxprops={'facecolor': BOX_FACE, 'edgecolor': BOX_EDGE, 'linewidth': 1.0},
               medianprops={'color': BOX_MEDIAN, 'linewidth': 1.4},
               whiskerprops={'color': BOX_EDGE, 'linewidth': 0.8},
               capprops={'color': BOX_EDGE, 'linewidth': 0.8})
    dots = np.percentile(values, BOX_PERCENTILES['dots'])
    ones = np.ones_like(dots)
    xy = (dots, ones) if orientation == 'horizontal' else (ones, dots)
    ax.plot(*xy, 'o', ms=5, mfc=BOX_FACE, mec='none', clip_on=False)
    ax.set_axis_off()


def plot_uncertainty_MPSP_vs_TCI(results_file=None, baseline_file=None,
                                 out_dir=DEFAULT_OUT_DIR, stem=None, dpi=600):
    results_file = results_file or newest_results_file('A')
    if baseline_file is None:
        candidate = results_file.replace('_1_full_evaluation', '_0_baseline')
        baseline_file = candidate if os.path.exists(candidate) else None
    print(f'Results:  {results_file}')
    print(f'Baseline: {baseline_file}')

    tci, mpsp, n_dropped = load_samples(results_file)
    rho, p = stats.spearmanr(tci, mpsp)
    base = load_baseline(baseline_file) if baseline_file else None
    print(f'{tci.size} samples ({n_dropped} non-finite dropped); '
          f'Spearman rho = {rho:.3f} (p = {p:.1e})')
    for name, v in (('TCI [MM$]', tci), ('MPSP [$/kg]', mpsp)):
        q = np.percentile(v, [1, 5, 25, 50, 75, 95, 99])
        print(f'  {name}: percentiles 1/5/25/50/75/95/99 = '
              + ' / '.join(f'{i:.4g}' for i in q))
    if base:
        print(f'  baseline: TCI {base[0]:.4g} MM$, MPSP {base[1]:.5g} $/kg')

    apply_font_rcparams()
    fig = plt.figure(figsize=(5.6, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=(4.2, 0.55), height_ratios=(0.55, 4.2),
                          wspace=0.03, hspace=0.03,
                          left=0.19, right=0.96, bottom=0.15, top=0.97)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    # fixed MPSP axis from zero (spans the market band and every sample)
    xticks = nice_ticks(tci)
    yticks = nice_ticks(np.array(MPSP_AXIS_LIMITS))
    if not (yticks[0] <= mpsp.min() and mpsp.max() <= yticks[-1]):
        raise ValueError(f'MPSP_AXIS_LIMITS {MPSP_AXIS_LIMITS} clip samples '
                         f'({mpsp.min():.3g}-{mpsp.max():.3g} $/kg)')
    xlim, ylim = (xticks[0], xticks[-1]), (yticks[0], yticks[-1])
    ax.axhspan(*ETHANOL_MARKET_RANGE, color=MARKET_BAND_COLOR, lw=0, zorder=0)
    draw_samples(ax, tci, mpsp)
    draw_hdr_contours(ax, tci, mpsp, xlim, ylim)
    if base:
        ax.plot(*base, 'D', ms=8, mfc='w', mec=INK, mew=1.2, zorder=5)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel('Total capital investment [MM\\$]',
                  fontsize=FONTS['axis_title'])
    ax.set_ylabel('Minimum ethanol selling price '
                  r'[$\mathrm{\$·kg}^{-1}$]', fontsize=FONTS['axis_title'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    draw_box(ax_top, tci, 'horizontal')
    draw_box(ax_right, mpsp, 'vertical')
    ax_top.set_ylim(0.4, 1.6)
    ax_right.set_xlim(0.4, 1.6)

    fig.canvas.draw()
    style_ticks(ax)

    os.makedirs(out_dir, exist_ok=True)
    if stem is None:
        tag = os.path.basename(results_file).split('_A_1_full_evaluation')[0]
        tag = tag.lstrip('_').replace("['A']_", '').replace('IBO_', '')
        stem = f'MPSP_vs_TCI_A_{tag}'
    paths = []
    for ext in ('png', 'pdf'):
        path = os.path.join(out_dir, f'{stem}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        paths.append(path)
        print(f'Saved {path}')
    plt.close(fig)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('results_file', nargs='?', default=None)
    parser.add_argument('--baseline-file', default=None)
    parser.add_argument('--out-dir', default=DEFAULT_OUT_DIR)
    parser.add_argument('--stem', default=None)
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    plot_uncertainty_MPSP_vs_TCI(args.results_file, args.baseline_file,
                                 args.out_dir, args.stem, args.dpi)


if __name__ == '__main__':
    main()

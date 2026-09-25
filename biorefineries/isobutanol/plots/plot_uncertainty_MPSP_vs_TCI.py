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

Joint panel: every Monte Carlo sample as a translucent dot, the 50 % and
90 % highest-density regions of a Gaussian KDE (the contours that enclose
that share of the samples), the baseline (the 'initial' row of the
companion *_0_baseline.xlsx) and the Spearman rank correlation. Marginal
panels: the KDE of each metric with its median.

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
from matplotlib.lines import Line2D, TICKDOWN, TICKLEFT
from matplotlib.ticker import AutoMinorLocator
from scipy import stats

__all__ = ('plot_uncertainty_MPSP_vs_TCI',)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(os.path.dirname(HERE), 'analyses', 'results')
DEFAULT_OUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Uncertainty')

MPSP_COL = ('Biorefinery', 'Purity-adjusted ethanol MPSP [$/kg]')
TCI_COL = ('Biorefinery', 'Total capital investment [10^6 $]')

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'clabel': 10, 'legend': 9}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# dataviz reference palette: one series -> sequential blue (dots at step 450,
# density contours at step 650); text stays in ink tokens
DOT_COLOR = '#2a78d6'
CONTOUR_COLOR = '#104281'
INK = '#0b0b0b'
INK_SECONDARY = '#52514e'
HDR_LEVELS = (0.90, 0.50) # share of samples enclosed, outer first


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
    pattern = os.path.join(RESULTS_DIR, f"*['{scenario}']*_{scenario}_1_full_evaluation.xlsx")
    files = [f for f in glob.glob(glob.escape(RESULTS_DIR) + os.sep + '*_1_full_evaluation.xlsx')
             if f"['{scenario}']" in os.path.basename(f)]
    if not files:
        raise FileNotFoundError(f'no scenario-{scenario} results workbook matching {pattern}')
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


#%% Density

def kde_hdr(x, y, levels=HDR_LEVELS, n_grid=200, pad=0.12):
    """Gaussian KDE on a grid plus the density thresholds whose superlevel
    sets hold `levels` of the samples (density quantiles AT the samples)."""
    kde = stats.gaussian_kde(np.vstack([x, y]))
    dx, dy = np.ptp(x), np.ptp(y)
    gx = np.linspace(x.min() - pad*dx, x.max() + pad*dx, n_grid)
    gy = np.linspace(y.min() - pad*dy, y.max() + pad*dy, n_grid)
    GX, GY = np.meshgrid(gx, gy)
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    at_samples = kde(np.vstack([x, y]))
    thresholds = [np.quantile(at_samples, 1 - p) for p in levels]
    return GX, GY, Z, thresholds


def draw_marginal(ax, values, median, baseline, orientation):
    kde = stats.gaussian_kde(values)
    lo, hi = values.min(), values.max()
    pad = 0.12 * (hi - lo)
    grid = np.linspace(lo - pad, hi + pad, 400)
    dens = kde(grid)
    if orientation == 'horizontal': # top panel: x = metric
        ax.fill_between(grid, 0, dens, color=DOT_COLOR, alpha=0.15, lw=0)
        ax.plot(grid, dens, color=DOT_COLOR, lw=2)
        ax.axvline(median, color=INK_SECONDARY, lw=1)
        ax.plot([baseline], [kde(baseline)[0]], 'D', ms=7, mfc='w', mec=INK, mew=1.2,
                zorder=5)
        ax.set_ylim(0, dens.max() * 1.15)
    else: # right panel: y = metric
        ax.fill_betweenx(grid, 0, dens, color=DOT_COLOR, alpha=0.15, lw=0)
        ax.plot(dens, grid, color=DOT_COLOR, lw=2)
        ax.axhline(median, color=INK_SECONDARY, lw=1)
        ax.plot([kde(baseline)[0]], [baseline], 'D', ms=7, mfc='w', mec=INK, mew=1.2,
                zorder=5)
        ax.set_xlim(0, dens.max() * 1.15)
    for side in ('top', 'right', 'left', 'bottom'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom' if orientation == 'horizontal' else 'left'].set_visible(True)
    ax.spines['bottom' if orientation == 'horizontal' else 'left'].set_color(INK_SECONDARY)
    ax.tick_params(axis='both', which='both', length=0,
                   labelbottom=False, labelleft=False)


#%% Figure

def plot_uncertainty_MPSP_vs_TCI(results_file=None, baseline_file=None,
                                 out_dir=DEFAULT_OUT_DIR, stem=None, dpi=600):
    results_file = results_file or newest_results_file('A')
    if baseline_file is None:
        candidate = results_file.replace('_1_full_evaluation', '_0_baseline')
        baseline_file = candidate if os.path.exists(candidate) else None
    print(f'Results:  {results_file}')
    print(f'Baseline: {baseline_file}')

    tci, mpsp, n_dropped = load_samples(results_file)
    n = tci.size
    rho, p = stats.spearmanr(tci, mpsp)
    base = load_baseline(baseline_file) if baseline_file else None
    print(f'{n} samples ({n_dropped} non-finite dropped); Spearman rho = {rho:.3f} (p = {p:.1e})')
    for name, v in (('TCI [MM$]', tci), ('MPSP [$/kg]', mpsp)):
        q = np.percentile(v, [5, 50, 95])
        print(f'  {name}: median {q[1]:.4g}, 5-95 % [{q[0]:.4g}, {q[2]:.4g}]')
    if base:
        print(f'  baseline: TCI {base[0]:.4g} MM$, MPSP {base[1]:.5g} $/kg')

    apply_font_rcparams()
    fig = plt.figure(figsize=(6.0, 6.0))
    gs = fig.add_gridspec(2, 2, width_ratios=(4.2, 1.0), height_ratios=(1.0, 4.2),
                          wspace=0.04, hspace=0.04,
                          left=0.14, right=0.97, bottom=0.18, top=0.97)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    # joint panel
    ax.scatter(tci, mpsp, s=5, color=DOT_COLOR, alpha=0.22, lw=0, rasterized=True,
               zorder=2)
    GX, GY, Z, thresholds = kde_hdr(tci, mpsp)
    cs = ax.contour(GX, GY, Z, levels=sorted(thresholds), colors=CONTOUR_COLOR,
                    linewidths=(1.2, 1.2), zorder=3)
    fmt = {lev: f'{int(round(100*share))}%'
           for lev, share in zip(sorted(thresholds), sorted(HDR_LEVELS, reverse=True))}
    labels = ax.clabel(cs, fmt=fmt, fontsize=FONTS['clabel'], inline=True,
                       inline_spacing=4)
    for t in labels:
        t.set_color(INK)
    if base:
        ax.plot(*base, 'D', ms=9, mfc='w', mec=INK, mew=1.4, zorder=5)
    ax.text(0.97, 0.97, rf'Spearman $\rho$ = {rho:.2f}'.replace('-', '−'),
            transform=ax.transAxes, ha='right', va='top', fontsize=FONTS['legend'],
            color=INK)

    ax.set_xlabel(r'Total capital investment [MM\$]', fontsize=FONTS['axis_title'])
    ax.set_ylabel(r'Ethanol MPSP [$\mathrm{\$·kg}^{-1}$]', fontsize=FONTS['axis_title'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # marginals
    draw_marginal(ax_top, tci, np.median(tci), base[0] if base else np.nan, 'horizontal')
    draw_marginal(ax_right, mpsp, np.median(mpsp), base[1] if base else np.nan, 'vertical')

    # legend under the panels
    handles = [
        Line2D([], [], ls='none', marker='o', ms=5, mfc=DOT_COLOR, mec='none',
               alpha=0.6, label=f'Monte Carlo sample (n = {n:,})'),
        Line2D([], [], color=CONTOUR_COLOR, lw=1.2,
               label='Region holding 50% / 90% of samples'),
        Line2D([], [], color=INK_SECONDARY, lw=1, label='Median (marginals)'),
    ]
    if base:
        handles.insert(1, Line2D([], [], ls='none', marker='D', ms=7, mfc='w',
                                 mec=INK, mew=1.2, label='Baseline'))
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
               fontsize=FONTS['legend'], bbox_to_anchor=(0.5, 0.0),
               handletextpad=0.5, columnspacing=1.5)

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
        fig.savefig(path, dpi=dpi, transparent=False, facecolor='white')
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

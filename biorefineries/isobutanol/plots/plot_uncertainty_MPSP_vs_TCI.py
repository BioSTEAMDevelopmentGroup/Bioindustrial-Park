#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Bivariate uncertainty plot of the purity-adjusted ethanol MPSP (y, converted
from the workbook's $/kg to $/GGE, see USD_PER_KG_TO_USD_PER_GGE) against
the total capital investment (x), from an uncertainty-analysis results
workbook (*_1_full_evaluation.xlsx as written by
analyses/full/uncertainties_IBO_EtOH.py).

Joint panel: every Monte Carlo sample as a teal dot at 25 % opacity and the
baseline (the 'initial' row of the companion *_0_baseline.xlsx) as a white
diamond (unlabelled: name it in the caption), over a light grey band for the
ethanol market price range (ETHANOL_MARKET_RANGE) spanning the typical corn
ethanol biorefinery TCI (TYPICAL_CORN_ETHANOL_TCI), the MPSP-TCI Pareto
frontier of the samples (lower-left, both minimized) as a blue staircase, and dark grey dashed lines
at the ends of the gasoline price range (GASOLINE_PRICE_RANGE), with contour
lines
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
FONTS = {'tick': 12, 'axis_title': 12}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# the uninformed profitability campaign's colour (HUE_COLORS[0] of
# plots/plot_kin_opt_parameter_sets.py), for the MPSP-TCI Pareto frontier
PARETO_COLOR = '#18C4DC'
PARETO_LW = 1.5
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
INK = '#0b0b0b'

# MPSP is plotted per gasoline gallon equivalent: $/GGE = $/kg * KG_PER_GAL
# / GGE_PER_GAL. kg per gal from ethanol's density at 20 C, 0.789 kg/L (the
# value analyses/full/uncertainties_IBO_EtOH.py uses to put the market range
# in $/kg) x 3.785411784 L/gal = 2.987 kg/gal; 1 gal ethanol = 0.67 GGE from
# the AFDC fuel properties comparison, https://afdc.energy.gov/fuels/properties
ETHANOL_DENSITY_KG_PER_L = 0.789
L_PER_GAL = 3.785411784
KG_PER_GAL = ETHANOL_DENSITY_KG_PER_L * L_PER_GAL
GGE_PER_GAL = 0.67
USD_PER_KG_TO_USD_PER_GGE = KG_PER_GAL / GGE_PER_GAL # ~4.458

# ethanol market price range, the same one uncertainties_IBO_EtOH.py draws on
# its MPSP box plot: Jan 2021 - Dec 2025 five-year low and high, 1.5475 and
# 3.4500 $/gal, from https://tradingeconomics.com/commodity/ethanol; converted
# from $/gal directly (not via its rounded $/kg values)
ETHANOL_MARKET_RANGE_PER_GAL = (1.5475, 3.4500) # $/gal
ETHANOL_MARKET_RANGE = tuple(v / GGE_PER_GAL for v in ETHANOL_MARKET_RANGE_PER_GAL) # $/GGE
# US retail gasoline price range, Jan 2021 - Dec 2025 low and high, per EIA
# Gasoline and Diesel Fuel Update, https://www.eia.gov/petroleum/gasdiesel/
# (a gallon of gasoline is 1 GGE by definition)
GASOLINE_PRICE_RANGE = (2.16, 4.84) # $/GGE
MPSP_AXIS_LIMITS = (0.0, 6.0) # $/GGE
TCI_AXIS_LIMITS = (75.0, 200.0) # MM$
TCI_TICK_STEP = 25.0 # MM$
N_MINOR_PER_MAJOR = 4 # minor ticks between adjacent major ticks, both axes
# the ethanol market band spans only the typical total capital investment of a
# corn ethanol biorefinery
# TODO: cite the source of the 100-150 MM$ typical corn-ethanol TCI range
TYPICAL_CORN_ETHANOL_TCI = (100.0, 150.0) # MM$
# the ethanol range is a band in a light shade of the baseline grey of
# plots/plot_kin_opt_parameter_sets.py (BASELINE_COLOR); the gasoline range,
# which almost coincides with it, is two dashed lines in a dark shade of the
# same grey. Both unlabelled, named in the caption.
BASELINE_GRAY = '#90918e'
MARKET_BAND_COLOR = _mix(BASELINE_GRAY, 'white', 0.7)
GASOLINE_LINE_COLOR = _mix(BASELINE_GRAY, 'black', 0.45)
GASOLINE_LINE_STYLE = dict(lw=1.0, ls=(0, (5, 3)))
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


def pareto_frontier(x, y):
    """The non-dominated samples when minimizing both x and y, sorted by
    ascending x (so descending y)."""
    order = np.lexsort((y, x)) # by x, ties by y
    front, best = [], np.inf
    for i in order:
        if y[i] < best:
            front.append(i)
            best = y[i]
    front = np.array(front)
    return x[front], y[front]


def draw_pareto_frontier(ax, x, y):
    # staircase: between two frontier samples the lowest attainable MPSP is
    # the left one's, so step horizontally, then down
    fx, fy = pareto_frontier(x, y)
    ax.step(fx, fy, where='post', color=PARETO_COLOR, lw=PARETO_LW, zorder=4)
    return fx, fy


def draw_hdr_contours(ax, x, y, n_grid=300, pad=0.15):
    """KDE contour lines enclosing CONTOUR_SHARES of the samples (unlabelled:
    name the shares in the caption)."""
    kde = stats.gaussian_kde(np.vstack([x, y]))
    dx, dy = np.ptp(x), np.ptp(y)
    GX, GY = np.meshgrid(np.linspace(x.min() - pad*dx, x.max() + pad*dx, n_grid),
                         np.linspace(y.min() - pad*dy, y.max() + pad*dy, n_grid))
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    at_samples = kde(np.vstack([x, y]))
    # the region holding share p is {density >= the (1 - p) quantile of the
    # density at the samples}; larger share -> lower level
    levels = {np.quantile(at_samples, 1 - p): p for p in CONTOUR_SHARES}
    ax.contour(GX, GY, Z, levels=sorted(levels), colors=CONTOUR_COLOR,
                    linewidths=CONTOUR_LW, zorder=3)


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

    tci, mpsp_kg, n_dropped = load_samples(results_file)
    mpsp = mpsp_kg * USD_PER_KG_TO_USD_PER_GGE
    rho, p = stats.spearmanr(tci, mpsp)
    base = load_baseline(baseline_file) if baseline_file else None
    if base:
        base = (base[0], base[1] * USD_PER_KG_TO_USD_PER_GGE)
    print(f'{tci.size} samples ({n_dropped} non-finite dropped); '
          f'Spearman rho = {rho:.3f} (p = {p:.1e})')
    print(f'MPSP conversion: {KG_PER_GAL:.4f} kg/gal / {GGE_PER_GAL} GGE/gal = '
          f'{USD_PER_KG_TO_USD_PER_GGE:.4f} ($/GGE)/($/kg); market range '
          + ' - '.join(f'{v:.3f}' for v in ETHANOL_MARKET_RANGE) + ' $/GGE')
    for name, v in (('TCI [MM$]', tci), ('MPSP [$/GGE]', mpsp)):
        q = np.percentile(v, [1, 5, 25, 50, 75, 95, 99])
        print(f'  {name}: percentiles 1/5/25/50/75/95/99 = '
              + ' / '.join(f'{i:.4g}' for i in q))
    if base:
        print(f'  baseline: TCI {base[0]:.4g} MM$, MPSP {base[1]:.5g} $/GGE')

    apply_font_rcparams()
    fig = plt.figure(figsize=(5.6, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=(4.2, 0.55), height_ratios=(0.55, 4.2),
                          wspace=0.03, hspace=0.03,
                          left=0.19, right=0.96, bottom=0.15, top=0.97)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)

    # fixed MPSP axis from zero (spans the market band and every sample)
    xticks = np.arange(TCI_AXIS_LIMITS[0], TCI_AXIS_LIMITS[1] + TCI_TICK_STEP/2,
                       TCI_TICK_STEP)
    yticks = nice_ticks(np.array(MPSP_AXIS_LIMITS))
    if not (xticks[0] <= tci.min() and tci.max() <= xticks[-1]):
        raise ValueError(f'TCI_AXIS_LIMITS {TCI_AXIS_LIMITS} clip samples '
                         f'({tci.min():.4g}-{tci.max():.4g} MM$)')
    if not (yticks[0] <= mpsp.min() and mpsp.max() <= yticks[-1]):
        raise ValueError(f'MPSP_AXIS_LIMITS {MPSP_AXIS_LIMITS} clip samples '
                         f'({mpsp.min():.3g}-{mpsp.max():.3g} $/GGE)')
    xlim, ylim = (xticks[0], xticks[-1]), (yticks[0], yticks[-1])
    ax.fill_between(TYPICAL_CORN_ETHANOL_TCI, *ETHANOL_MARKET_RANGE,
                    color=MARKET_BAND_COLOR, lw=0, zorder=0)
    for price in GASOLINE_PRICE_RANGE:
        ax.axhline(price, color=GASOLINE_LINE_COLOR, zorder=1, **GASOLINE_LINE_STYLE)
    draw_samples(ax, tci, mpsp)
    draw_hdr_contours(ax, tci, mpsp)
    fx, fy = draw_pareto_frontier(ax, tci, mpsp)
    print(f'Pareto frontier: {fx.size} non-dominated samples, TCI '
          f'{fx[0]:.4g}-{fx[-1]:.4g} MM$, MPSP {fy[0]:.4g}-{fy[-1]:.4g} $/GGE')
    if base:
        ax.plot(*base, 'D', ms=8, mfc='w', mec=INK, mew=1.2, zorder=5)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel('Total capital investment [MM\\$]',
                  fontsize=FONTS['axis_title'])
    ax.set_ylabel('Minimum ethanol selling price '
                  r'[$\mathrm{\$·GGE}^{-1}$]', fontsize=FONTS['axis_title'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))

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

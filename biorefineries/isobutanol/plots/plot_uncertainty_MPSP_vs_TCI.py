#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Multi-panel bivariate uncertainty figure (2 x 2 grid, fourth cell reserved)
from an uncertainty-analysis results workbook (*_1_full_evaluation.xlsx as
written by analyses/full/uncertainties_IBO_EtOH.py):

  A  purity-adjusted ethanol MPSP (y, converted from the workbook's $/kg to
     $/GGE, see USD_PER_KG_TO_USD_PER_GGE) vs total capital investment (x)
  B  ethanol titer (y, g/L-water) vs ethanol yield (x, g/g sugars added)
  C  ethanol sale revenue (y) vs DDGS sale revenue (x), MM$/yr at the default
     product prices; ethanol revenue is not a workbook column and is derived
     as the annual product sale minus the DDGS, crude-oil and isobutanol sale
     revenues (exact: the model's tea.sales is the sum of those four)
  D  empty (a fourth distribution will be added)

Every joint panel, in its own hue (A teal, B purple, C yellow): every Monte
Carlo sample as a dot at 35 % opacity,
the baseline (the 'initial' row of the companion *_0_baseline.xlsx) as a
white diamond (unlabelled: name it in the caption), contour lines of a
Gaussian KDE enclosing 5 / 25 / 50 / 75 / 95 % of the samples (the
highest-density regions; each line is the density quantile AT the samples,
so it holds that share of them), and marginal box plots outside the panel:
box = 25th-75th percentile, line = median, whiskers = 5th-95th percentile
(the whis=[5, 95] of contourplots.box_and_whiskers_plot), dots = 1st and 99th
percentiles.

Panel A also carries a light grey band for the ethanol market price range
(ETHANOL_MARKET_RANGE) spanning the typical corn ethanol biorefinery TCI
(TYPICAL_CORN_ETHANOL_TCI), the MPSP-TCI Pareto frontier of the samples
(lower-left, both minimized) as a red dashed staircase, and dark grey dashed
lines at the ends of the gasoline price range (GASOLINE_PRICE_RANGE).

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
ETOH_TITER_COL = ('Fermentation', 'Et OH titer [g-EtOH/L-water]')
ETOH_YIELD_COL = ('Fermentation', 'Et OH yield [g-EtOH/g-sugars-added]')
SALES_COL = ('Biorefinery', 'Annual product sale (excl. electricity) [10^6 $/yr]')
DDGS_REVENUE_COL = ('Coproducts', 'DDGS sale revenue [$/y]')
# every product sale revenue other than ethanol's; ethanol's = SALES_COL
# minus these
COPRODUCT_REVENUE_COLS = (DDGS_REVENUE_COL,
                          ('Coproducts', 'Crude oil sale revenue [$/y]'),
                          ('Coproducts', 'Isobutanol sale revenue [$/y]'))

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'panel_letter': 14}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# the MPSP-TCI Pareto frontier: dashed, in the red of the hue palette of
# plots/plot_kin_opt_parameter_sets.py (HUE_COLORS[3])
PARETO_COLOR = '#ED586F'
PARETO_LW = 1.5
PARETO_LS = (0, (4, 2))
# the TRY-informed profitability campaign's teal (RELAY_COLOR in
# plots/plot_kin_opt_parameter_sets.py)
TEAL = '#0B6E7A'
# panels B and C: the ethanol-yield purple and ethanol-titer yellow of that
# figure's hue palette (HUE_COLORS[4], HUE_COLORS[5])
PURPLE = '#a280b9'
YELLOW = '#f3c354'
SAMPLE_ALPHA = 0.35   # joint-panel samples: 65 % transparent
SAMPLE_SIZE = 6 # pt^2


def _mix(color, other, t):
    """`color` moved a fraction `t` of the way towards `other` (RGB)."""
    a, b = np.array(to_rgb(color)), np.array(to_rgb(other))
    return to_hex((1 - t)*a + t*b)


# each panel has one hue: samples and box faces in it, box edges / whiskers
# and the KDE contour lines in a dark shade of it (dark_shade)
dark_shade = lambda color: _mix(color, 'black', 0.5)
BOX_MEDIAN = 'black'
# KDE contour lines: share of samples each encloses
CONTOUR_SHARES = (0.05, 0.25, 0.50, 0.75, 0.95)
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
MPSP_AXIS_LIMITS = (2.0, 6.0) # $/GGE
TCI_AXIS_LIMITS = (75.0, 200.0) # MM$
TCI_TICK_STEP = 25.0 # MM$
N_MINOR_PER_MAJOR = 4 # minor ticks between adjacent major ticks, both axes
# the ethanol market band spans only the typical total capital investment of a
# conventional dry-grind corn ethanol biorefinery:
# - low end, ~40 MGY: Kurambhatti et al. (2019) report $83.95 million for the
#   conventional dry-grind process at 40.2 million gal/yr, a TCI (direct fixed
#   capital + working capital at 5 % of DFC)
# - high end: $135 million, from U.S. DOE EERE (2010), Current State of the
#   U.S. Ethanol Industry, https://www1.eere.energy.gov/bioenergy/pdfs/
#   current_state_of_the_us_ethanol_industry.pdf (raised from Tanzil et al.
#   (2021)'s $115 million for a typical 230 ML/yr (~60.8 MGY) mill)
# Each bound is escalated from its source's year to the model's 2023 dollars by
# CEPCI ratio (annual averages, Chemical Engineering magazine; 2023 = 797.9,
# process_settings.CEPCI). Years are the sources' publication years.
CEPCI_ANNUAL = {2010: 550.8, 2019: 607.5, 2023: 797.9}
TCI_COST_YEAR = 2023
TYPICAL_CORN_ETHANOL_TCI_SOURCE = ((83.95, 2019), # MM$, Kurambhatti et al.
                                   (135.0, 2010)) # MM$, DOE EERE
TYPICAL_CORN_ETHANOL_TCI = tuple(
    v * CEPCI_ANNUAL[TCI_COST_YEAR] / CEPCI_ANNUAL[year]
    for v, year in TYPICAL_CORN_ETHANOL_TCI_SOURCE) # MM$ (2023$): ~110.3, ~195.6
# the ethanol range is a band in a light shade of the baseline grey of
# plots/plot_kin_opt_parameter_sets.py (BASELINE_COLOR); the gasoline range,
# which almost coincides with it, is two dashed lines in a dark shade of the
# same grey. Both unlabelled, named in the caption.
BASELINE_GRAY = '#90918e'
MARKET_BAND_COLOR = _mix(BASELINE_GRAY, 'white', 0.7)
GASOLINE_LINE_COLOR = _mix(BASELINE_GRAY, 'black', 0.45)
GASOLINE_LINE_STYLE = dict(lw=1.0, ls=(0, (5, 3)))
BOX_PERCENTILES = {'whis': (5, 95), 'dots': (1, 99)}

# the joint panels of the 2 x 2 grid, row-major (the fourth cell is left empty
# for now): (x outcome, y outcome) keyed by _outcomes' names, axis titles, and
# panel hue, fixed axis ticks / limits (neither = 'nice' ticks enclosing the samples and
# the baseline)
PANELS = (
    dict(x='TCI', y='MPSP',
         xlabel='Total capital investment [MM\\$]',
         ylabel=r'Minimum ethanol selling price [$\mathrm{\$·GGE}^{-1}$]',
         xticks=np.arange(TCI_AXIS_LIMITS[0], TCI_AXIS_LIMITS[1] + TCI_TICK_STEP/2,
                          TCI_TICK_STEP),
         ylim=MPSP_AXIS_LIMITS, market=True, pareto=True, color=TEAL),
    dict(x='EtOH yield', y='EtOH titer', color=PURPLE,
         xlabel=r'Ethanol yield [$\mathrm{g·g}^{-1}$ sugars]',
         ylabel=r'Ethanol titer [$\mathrm{g·L}^{-1}$]'),
    dict(x='DDGS revenue', y='EtOH revenue', color=YELLOW,
         xlabel=r'DDGS sale revenue [$\mathrm{MM\$·yr}^{-1}$]',
         ylabel=r'Ethanol sale revenue [$\mathrm{MM\$·yr}^{-1}$]'),
)
GRID_SHAPE = (2, 2)


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


def _outcomes(df):
    """The plotted outcomes of a TEA-results frame, as float Series, keyed by
    the names used in PANELS."""
    col = lambda c: df[c].astype(float)
    return {
        'TCI': col(TCI_COL),
        'MPSP': col(MPSP_COL) * USD_PER_KG_TO_USD_PER_GGE,
        'EtOH titer': col(ETOH_TITER_COL),
        'EtOH yield': col(ETOH_YIELD_COL),
        'DDGS revenue': col(DDGS_REVENUE_COL) / 1e6,
        'EtOH revenue': col(SALES_COL) - sum(col(c) for c in COPRODUCT_REVENUE_COLS) / 1e6,
    }


def load_samples(results_file):
    df = pd.read_excel(results_file, sheet_name='TEA results',
                       header=[0, 1], index_col=0)
    return pd.DataFrame(_outcomes(df))


def load_baseline(baseline_file):
    df = pd.read_excel(baseline_file, header=[0, 1], index_col=0)
    return {k: float(v.iloc[0]) for k, v in _outcomes(df.loc[['initial']]).items()}


#%% Drawing

def nice_ticks(values, nbins=7):
    """'Nice' major ticks enclosing the data; the axis limits sit on the
    first and last one."""
    ticks = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(
        values.min(), values.max())
    return ticks[(ticks <= values.min()).nonzero()[0][-1]:
                 (ticks >= values.max()).nonzero()[0][0] + 1]


def draw_samples(ax, x, y, color):
    # every sample; rasterized so the PDF stays small (axes, text and the
    # baseline marker stay vector)
    ax.scatter(x, y, s=SAMPLE_SIZE, color=color, alpha=SAMPLE_ALPHA, lw=0,
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
    ax.step(fx, fy, where='post', color=PARETO_COLOR, lw=PARETO_LW, ls=PARETO_LS,
            zorder=4)
    return fx, fy


def draw_hdr_contours(ax, x, y, color, n_grid=300, pad=0.15):
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
    ax.contour(GX, GY, Z, levels=sorted(levels), colors=dark_shade(color),
                    linewidths=CONTOUR_LW, zorder=3)


def draw_box(ax, values, orientation, color):
    edge = dark_shade(color)
    lo_w, hi_w = BOX_PERCENTILES['whis']
    ax.boxplot(values, whis=[lo_w, hi_w], orientation=orientation,
               widths=0.6, showfliers=False, patch_artist=True,
               boxprops={'facecolor': color, 'edgecolor': edge, 'linewidth': 1.0},
               medianprops={'color': BOX_MEDIAN, 'linewidth': 1.4},
               whiskerprops={'color': edge, 'linewidth': 0.8},
               capprops={'color': edge, 'linewidth': 0.8})
    dots = np.percentile(values, BOX_PERCENTILES['dots'])
    ones = np.ones_like(dots)
    xy = (dots, ones) if orientation == 'horizontal' else (ones, dots)
    ax.plot(*xy, 'o', ms=5, mfc=color, mec='none', clip_on=False)
    ax.set_axis_off()


def _panel_ticks(panel, axis, values):
    """Major ticks of one axis of a panel: the fixed `<axis>ticks`, 'nice'
    ticks spanning the fixed `<axis>lim`, or 'nice' ticks enclosing `values`.
    Fixed ticks that would clip a sample raise."""
    ticks = panel.get(f'{axis}ticks')
    if ticks is None and panel.get(f'{axis}lim') is not None:
        ticks = nice_ticks(np.array(panel[f'{axis}lim']))
    if ticks is None:
        return nice_ticks(values, nbins=6)
    if not (ticks[0] <= values.min() and values.max() <= ticks[-1]):
        raise ValueError(f'{panel[axis]} axis {ticks[0]:.4g}-{ticks[-1]:.4g} clips '
                         f'samples ({values.min():.4g}-{values.max():.4g})')
    return np.asarray(ticks)


def draw_joint_panel(fig, cell, panel, samples, base, letter):
    """One joint panel (samples, KDE contours, baseline, marginal boxes) in
    the grid cell `cell`; returns the joint axes."""
    gs = cell.subgridspec(2, 2, width_ratios=(4.2, 0.55), height_ratios=(0.55, 4.2),
                          wspace=0.03, hspace=0.03)
    ax = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax)
    x, y = samples[panel['x']].values, samples[panel['y']].values

    with_base = lambda k, v: np.append(v, base[panel[k]]) if base else v
    xticks = _panel_ticks(panel, 'x', with_base('x', x))
    yticks = _panel_ticks(panel, 'y', with_base('y', y))
    if panel.get('market'):
        ax.fill_between(TYPICAL_CORN_ETHANOL_TCI, *ETHANOL_MARKET_RANGE,
                        color=MARKET_BAND_COLOR, lw=0, zorder=0)
        for price in GASOLINE_PRICE_RANGE:
            ax.axhline(price, color=GASOLINE_LINE_COLOR, zorder=1, **GASOLINE_LINE_STYLE)
    draw_samples(ax, x, y, panel['color'])
    draw_hdr_contours(ax, x, y, panel['color'])
    if panel.get('pareto'):
        fx, fy = draw_pareto_frontier(ax, x, y)
        print(f'  Pareto frontier: {fx.size} non-dominated samples, '
              f"{panel['x']} {fx[0]:.4g}-{fx[-1]:.4g}, {panel['y']} {fy[0]:.4g}-{fy[-1]:.4g}")
    if base:
        ax.plot(base[panel['x']], base[panel['y']], 'D', ms=8, mfc='w', mec=INK,
                mew=1.2, zorder=5)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    ax.set_xlabel(panel['xlabel'], fontsize=FONTS['axis_title'])
    ax.set_ylabel(panel['ylabel'], fontsize=FONTS['axis_title'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))

    draw_box(ax_top, x, 'horizontal', panel['color'])
    draw_box(ax_right, y, 'vertical', panel['color'])
    ax_top.set_ylim(0.4, 1.6)
    ax_right.set_xlim(0.4, 1.6)
    # panel letter in the cell's top-left corner, above the y-axis title
    ax_top.text(-0.22, 1.0, letter, transform=ax_top.transAxes, ha='left', va='top',
                fontsize=FONTS['panel_letter'], fontweight='bold')
    return ax


def plot_uncertainty_MPSP_vs_TCI(results_file=None, baseline_file=None,
                                 out_dir=DEFAULT_OUT_DIR, stem=None, dpi=600):
    results_file = results_file or newest_results_file('A')
    if baseline_file is None:
        candidate = results_file.replace('_1_full_evaluation', '_0_baseline')
        baseline_file = candidate if os.path.exists(candidate) else None
    print(f'Results:  {results_file}')
    print(f'Baseline: {baseline_file}')

    samples = load_samples(results_file)
    used = list(dict.fromkeys(panel[k] for panel in PANELS for k in ('x', 'y')))
    finite = np.isfinite(samples[used].values).all(axis=1)
    n_dropped = int((~finite).sum())
    samples = samples[finite]
    base = load_baseline(baseline_file) if baseline_file else None
    print(f'{len(samples)} samples ({n_dropped} with a non-finite plotted outcome dropped)')
    print(f'MPSP conversion: {KG_PER_GAL:.4f} kg/gal / {GGE_PER_GAL} GGE/gal = '
          f'{USD_PER_KG_TO_USD_PER_GGE:.4f} ($/GGE)/($/kg); market range '
          + ' - '.join(f'{v:.3f}' for v in ETHANOL_MARKET_RANGE) + ' $/GGE')
    for name in used:
        q = np.percentile(samples[name], [1, 5, 25, 50, 75, 95, 99])
        print(f'  {name}: percentiles 1/5/25/50/75/95/99 = '
              + ' / '.join(f'{i:.4g}' for i in q)
              + (f'; baseline {base[name]:.5g}' if base else ''))

    apply_font_rcparams()
    n_rows, n_cols = GRID_SHAPE
    fig = plt.figure(figsize=(5.6*n_cols, 5.4*n_rows))
    grid = fig.add_gridspec(n_rows, n_cols, wspace=0.38, hspace=0.26,
                            left=0.1, right=0.98, bottom=0.075, top=0.985)
    axes = []
    for i, panel in enumerate(PANELS):
        letter = chr(ord('A') + i)
        rho, p = stats.spearmanr(samples[panel['x']], samples[panel['y']])
        print(f"Panel {letter}: {panel['y']} vs {panel['x']}, "
              f'Spearman rho = {rho:.3f} (p = {p:.1e})')
        axes.append(draw_joint_panel(fig, grid[divmod(i, n_cols)], panel,
                                     samples, base, letter))

    fig.canvas.draw()
    for ax in axes:
        style_ticks(ax)

    os.makedirs(out_dir, exist_ok=True)
    if stem is None:
        tag = os.path.basename(results_file).split('_A_1_full_evaluation')[0]
        tag = tag.lstrip('_').replace("['A']_", '').replace('IBO_', '')
        stem = f'uncertainty_bivariate_A_{tag}'
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

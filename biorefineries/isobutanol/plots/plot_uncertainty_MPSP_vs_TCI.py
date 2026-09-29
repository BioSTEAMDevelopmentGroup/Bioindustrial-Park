#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Multi-panel bivariate uncertainty figure (2 x 3 grid)
from an uncertainty-analysis results workbook (*_1_full_evaluation.xlsx as
written by analyses/full/uncertainties_IBO_EtOH.py). MESP = the
purity-adjusted ethanol MPSP, converted from the workbook's $/kg to $/GGE (see
USD_PER_KG_TO_USD_PER_GGE):

  A  MESP (y) vs total capital investment (x)
  B  total capital investment (y) vs ethanol production (x, million gal of
     pure ethanol per year: the purity-adjusted production rate /
     ETHANOL_DENSITY_KG_PER_L / L_PER_GAL)
  C  MESP (y) vs feed corn starch content (x, kg/kg dry corn; a sampled
     parameter, read from the workbook's 'Parameters' sheet; its baseline from
     the scenario's parameter-distribution workbook, PARAMETER_DISTRIBUTIONS)
  D  MESP (y) vs ethanol yield (x, g/g sugars added)
  E  MESP (y) vs ethanol titer (x, g/L-water)
  F  MESP (y) vs ethanol productivity (x, g/L-water/h)

Every joint panel, in its own hue (A teal, B green, C orange, D purple,
E yellow, F blue): the
Gaussian-KDE density of the Monte Carlo samples as filled contours (one
sequential ramp of the hue, light = sparse to dark = dense; the lowest band is
left unfilled, so the panel background stays white),
the baseline (the 'initial' row of the companion *_0_baseline.xlsx) as a
white diamond (unlabelled: name it in the caption), and marginal box plots
outside the panel:
box = 25th-75th percentile, line = median, whiskers = 5th-95th percentile
(the whis=[5, 95] of contourplots.box_and_whiskers_plot), dots = the most
outlying samples (minimum and maximum).

Panel A also carries a light grey band for the ethanol market price range
(ETHANOL_MARKET_RANGE) spanning the typical corn ethanol biorefinery TCI
(TYPICAL_CORN_ETHANOL_TCI) and dark grey dashed lines at the ends of the
gasoline price range (GASOLINE_PRICE_RANGE). Panel B carries a box in the
same grey spanning a typical corn ethanol biorefinery's ethanol production
(TYPICAL_CORN_ETHANOL_PRODUCTION) and TCI (TYPICAL_CORN_ETHANOL_TCI), also
unlabelled.

Panels A and B also carry the Pareto frontier of their samples as a solid red
staircase, in the sense set by the panel's `pareto` entry: A lower-left (TCI
and MESP minimized), B lower-right (production maximized, TCI minimized).
Panels C-F carry a binned-median trend line instead (`trend`): the samples
split into N_TREND_BINS equal-count bins of x, a solid black line through each
bin's median x and median MESP. Their x is not a design choice that trades
against MESP (starch content is a feedstock property; yield, titer and
productivity are outcomes of the uncertain kinetics, and all four go WITH
lower MESP), so a frontier would only trace the lucky corn-price / starch /
capacity draws that happen to sit at high x; the binned medians show how much
MESP actually moves with x.

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
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
from matplotlib.lines import TICKDOWN, TICKLEFT
from matplotlib.transforms import ScaledTranslation
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
ETOH_PRODUCTIVITY_COL = ('Fermentation', 'Et OH productivity [g-EtOH/L-water/h]')
ETOH_PRODUCTION_COL = ('Biorefinery', 'Adjusted production rate [10^6 kg/yr]') # pure ethanol
# the feed corn starch content is a sampled parameter, not a metric: the
# samples are in the results workbook's 'Parameters' sheet, and the baseline
# workbook (metrics only) lacks it, so its baseline comes from the Baseline
# column of the scenario's parameter-distribution workbook
STARCH_COL = ('Feedstock', 'Feedstock starch content [kg/dry-kg]') # 'Parameters' sheet
STARCH_PARAMETER_NAME = 'Feedstock starch content'
PARAMETER_DISTRIBUTIONS = os.path.join(
    os.path.dirname(HERE), 'analyses', 'full', 'parameter_distributions',
    'parameter-distributions_corn_IBO_EtOH_A.xlsx')

FONT_FAMILY = 'Arial'
FONTS = {'tick': 13.8, 'axis_title': 13.8, 'panel_letter': 16.1} # 1.15 x (12, 12, 14) since 2026-09-29
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# the Pareto frontiers: solid, in the red of the hue palette of
# plots/plot_kin_opt_parameter_sets.py (HUE_COLORS[3])
PARETO_COLOR = '#ED586F'
PARETO_LW = 1.5
PARETO_LS = '-'
# panels C-F: binned-median trend line, solid black (no markers, no halo)
N_TREND_BINS = 10 # equal-count bins of x (deciles)
TREND_COLOR = '#0b0b0b'
TREND_LW = 1.5
TREND_LS = '-'
# the TRY-informed profitability campaign's teal (RELAY_COLOR in
# plots/plot_kin_opt_parameter_sets.py)
TEAL = '#0B6E7A'
# panels B-F, from that figure's hue palette: B the isobutanol-titer green
# (HUE_COLORS[2]), C the isobutanol-yield orange ([1]), and D-F the
# ethanol-yield purple, ethanol-titer yellow and ethanol-productivity blue
# ([4], [5], [6]), so each ethanol metric keeps its hue across figures
GREEN = '#79bf82'
ORANGE = '#f98f60'
PURPLE = '#a280b9'
YELLOW = '#f3c354'
BLUE = '#5a6bcc'


def _mix(color, other, t):
    """`color` moved a fraction `t` of the way towards `other` (RGB)."""
    a, b = np.array(to_rgb(color)), np.array(to_rgb(other))
    return to_hex((1 - t)*a + t*b)


# each panel has one hue: box faces in it, box edges / whiskers in a dark
# shade of it (dark_shade), the filled KDE on a ramp of it (density_ramp)
dark_shade = lambda color: _mix(color, 'black', 0.5)
BOX_MEDIAN = 'black'
# filled KDE: N_DENSITY_LEVELS equal steps from 0 to the peak density, coloured
# on ONE sequential ramp of the panel's hue (white-mixed tints up to the hue
# itself, then darker shades of it)
density_ramp = lambda color: tuple(
    [_mix(color, 'white', t) for t in (0.94, 0.78, 0.58, 0.36, 0.16)]
    + [color] + [_mix(color, 'black', t) for t in (0.35, 0.65)])
N_DENSITY_LEVELS = 10
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
TCI_AXIS_LIMITS = (100.0, 200.0) # MM$
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
# typical dry-grind corn ethanol biorefinery capacity, for panel B's box
TYPICAL_CORN_ETHANOL_PRODUCTION = (40.0, 60.0) # MM gal/y
# ethanol yield, titer and productivity reported for high-gravity fermentation
# in U.S. corn dry-grind facilities, for the boxes of panels D-F (each x the
# ethanol market price range, as panel A's box):
# - Gomes, D. et al. Very High Gravity Bioethanol Revisited: Main Challenges
#   and Advances. Fermentation 7, 38 (2021).
# - Deparis, Q., Claes, A., Foulquie-Moreno, M. R. & Thevelein, J. M.
#   Engineering tolerance to industrially relevant stress factors in yeast
#   cell factories. FEMS Yeast Res. 17 (2017).
# - Tsegaye, K. N., Alemnew, M. & Berhane, N. Saccharomyces cerevisiae for
#   lignocellulosic ethanol production: a look at key attributes and genome
#   shuffling. Front. Bioeng. Biotechnol. 12, 1466644 (2024).
# - Devantier, R., Pedersen, S. & Olsson, L. Characterization of very high
#   gravity ethanol fermentation of corn mash. Effect of glucoamylase dosage,
#   pre-saccharification and yeast strain. Appl. Microbiol. Biotechnol. 68,
#   622-629 (2005).
HIGH_GRAVITY_ETHANOL_YIELD = (0.42, 0.48) # g/g
HIGH_GRAVITY_ETHANOL_TITER = (100.0, 142.0) # g/L
HIGH_GRAVITY_ETHANOL_PRODUCTIVITY = (2.0, 4.4) # g/L/h
# the ethanol range is a band in a light shade of the baseline grey of
# plots/plot_kin_opt_parameter_sets.py (BASELINE_COLOR); the gasoline range,
# which almost coincides with it, is two dashed lines in a dark shade of the
# same grey. Both unlabelled, named in the caption.
BASELINE_GRAY = '#90918e'
MARKET_BAND_COLOR = _mix(BASELINE_GRAY, 'white', 0.7)
GASOLINE_LINE_COLOR = _mix(BASELINE_GRAY, 'black', 0.45)
GASOLINE_LINE_STYLE = dict(lw=1.0, ls=(0, (5, 3)))
BOX_PERCENTILES = {'whis': (5, 95), 'dots': (0, 100)} # dots: min and max

# the joint panels of the 2 x 3 grid, row-major: (x outcome, y outcome) keyed
# by _outcomes' names, axis titles, and panel hue, Pareto sense per axis
# ('min' / 'max', x then y; None = no frontier), binned-median trend line
# (`trend`), fixed axis ticks / limits (neither = 'nice' ticks enclosing the
# samples and the baseline)
TCI_TICKS = np.arange(TCI_AXIS_LIMITS[0], TCI_AXIS_LIMITS[1] + TCI_TICK_STEP/2,
                      TCI_TICK_STEP)
MESP_LABEL = r'Minimum ethanol selling price [$\mathrm{\$·GGE}^{-1}$]'
# panels C-F: axes fitted closely to the samples, so the marginal boxes'
# min / max dots sit near the limits (6000-sim scenario-A run, min - max:
# MESP 3.350-4.851 $/GGE, yield 0.4167-0.4707 g/g, titer 93.89-129.8 g/L,
# productivity 0.832-3.812 g/L/h); since 2026-09-29 their MESP axes are
# panel A's MPSP_AXIS_LIMITS (3.3-4.9 every 0.2 before), so the ethanol
# market-price rows of the D-F boxes show
ETOH_YIELD_TICKS = np.round(np.linspace(0.41, 0.48, 8), 10) # g/g, every 0.01
ETOH_TITER_TICKS = np.linspace(90.0, 130.0, 5) # g/L, every 10
ETOH_PRODUCTIVITY_TICKS = np.round(np.linspace(0.8, 4.0, 9), 10) # g/L/h, every 0.4
# panels C-F: MESP (y, MPSP_AXIS_LIMITS as panel A) vs one driver (x), with a binned-median
# trend line and no Pareto frontier (see the module docstring)
_mesp_panel = lambda x, xlabel, color, **kw: {
    **dict(x=x, y='MPSP', xlabel=xlabel, ylabel=MESP_LABEL, ylim=MPSP_AXIS_LIMITS,
           pareto=None, trend=True, color=color),
    **kw}
PANELS = (
    dict(x='TCI', y='MPSP',
         xlabel='Total capital investment [MM\\$]',
         ylabel=MESP_LABEL,
         xticks=TCI_TICKS,
         ylim=MPSP_AXIS_LIMITS, market=True, pareto=('min', 'min'), color=TEAL,
         callouts=True),
    dict(x='EtOH production', y='TCI', color=GREEN, pareto=('max', 'min'),
         xlabel=r'Ethanol production [$\mathrm{MM\ gal·y}^{-1}$]',
         ylabel='Total capital investment [MM\\$]',
         yticks=TCI_TICKS,
         box=(TYPICAL_CORN_ETHANOL_PRODUCTION, TYPICAL_CORN_ETHANOL_TCI)),
    _mesp_panel('Starch content',
                r'Corn starch content [$\mathrm{kg·kg}^{-1}$ dry corn]', ORANGE),
    _mesp_panel('EtOH yield', r'Ethanol yield [$\mathrm{g·g}^{-1}$ sugars]',
                PURPLE, xticks=ETOH_YIELD_TICKS,
                box=(HIGH_GRAVITY_ETHANOL_YIELD, ETHANOL_MARKET_RANGE)),
    _mesp_panel('EtOH titer', r'Ethanol titer [$\mathrm{g·L}^{-1}$]',
                YELLOW, xticks=ETOH_TITER_TICKS,
                box=(HIGH_GRAVITY_ETHANOL_TITER, ETHANOL_MARKET_RANGE)),
    _mesp_panel('EtOH productivity',
                r'Ethanol productivity [$\mathrm{g·L}^{-1}\mathrm{·h}^{-1}$]', BLUE,
                xticks=ETOH_PRODUCTIVITY_TICKS,
                box=(HIGH_GRAVITY_ETHANOL_PRODUCTIVITY, ETHANOL_MARKET_RANGE)),
)
GRID_SHAPE = (2, 3)
# figure size per grid cell and the figure margins (inches, so the cells keep
# their size and spacing whatever the grid shape), gaps between grid cells
# (fractions of the mean cell width / height), and the panel letter's
# position in its joint axes' coordinates
CELL_SIZE = (5.6, 5.4) # in
MARGINS_IN = {'left': 1.12, 'right': 0.224, 'bottom': 0.81, 'top': 0.162}
GRID_SPACE = {'w': 0.18, 'h': 0.145}
PANEL_LETTER_XY = (-0.15, 1.12)

# callouts (panel A, `callouts` in PANELS): bold labels with curved arrows,
# after the reference layout of the former 2x2
# figure; offsets in points from the annotated point
CALLOUT_FONTSIZE = 12.0
CALLOUT_ARROW = dict(arrowstyle='-|>', lw=1.0, mutation_scale=12, shrinkA=2)
# `tip` (points) moves the arrow head past the annotated point so it overlaps
# the item slightly; the baseline arrow instead stops `shrinkB` points short
# of the diamond's centre, just inside its edge
BASELINE_CALLOUT = dict(text='baseline', offset=(24, 22), rad=0.35, shrinkB=2.5)
PARETO_CALLOUT = dict(text='Pareto frontier', at=0.55, offset=(12, -22), rad=-0.35,
                      tip=(0, 1.5))
GASOLINE_CALLOUT = dict(text='gasoline market price range', x=105.0, tip=1.5)
CONVENTIONAL_CALLOUT = dict(text='conventional', at=(185.0, None), # None = box top
                            offset=(-20, 14), rad=-0.35, tip=(0, -2.0))
# marginal box axes: size relative to the joint axes' 4.2, box width in its axes
MARGINAL_RATIO = 0.45
BOX_WIDTH = 0.55


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def bold_title(label):
    """An axis title with its name in bold and its ' [units]' part regular:
    the name is set as mathtext bold (mathtext.bf = the bold of FONT_FAMILY)."""
    name, sep, units = label.partition(' [')
    name = r'$\mathbf{' + name.replace(' ', r'\ ') + '}$'
    return name + sep + units


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
    """The plotted metrics of a TEA-results frame, as float Series, keyed by
    the names used in PANELS (the sampled starch content is added by
    load_samples / load_baseline)."""
    col = lambda c: df[c].astype(float)
    return {
        'TCI': col(TCI_COL),
        'MPSP': col(MPSP_COL) * USD_PER_KG_TO_USD_PER_GGE,
        'EtOH titer': col(ETOH_TITER_COL),
        'EtOH yield': col(ETOH_YIELD_COL),
        'EtOH productivity': col(ETOH_PRODUCTIVITY_COL),
        'EtOH production': col(ETOH_PRODUCTION_COL) / ETHANOL_DENSITY_KG_PER_L / L_PER_GAL, # MM gal/y
    }


def load_samples(results_file):
    sheets = pd.read_excel(results_file, sheet_name=['TEA results', 'Parameters'],
                           header=[0, 1], index_col=0)
    samples = pd.DataFrame(_outcomes(sheets['TEA results']))
    samples['Starch content'] = sheets['Parameters'][STARCH_COL].astype(float)
    return samples


def load_baseline(baseline_file, parameter_distributions=PARAMETER_DISTRIBUTIONS):
    df = pd.read_excel(baseline_file, header=[0, 1], index_col=0)
    base = {k: float(v.iloc[0]) for k, v in _outcomes(df.loc[['initial']]).items()}
    dist = pd.read_excel(parameter_distributions)
    row = dist[dist['Parameter name'] == STARCH_PARAMETER_NAME]
    if len(row) != 1:
        raise ValueError(f'{len(row)} {STARCH_PARAMETER_NAME!r} rows in '
                         f'{parameter_distributions}')
    base['Starch content'] = float(row['Baseline'].iloc[0])
    return base


#%% Drawing

def nice_ticks(values, nbins=7):
    """'Nice' major ticks enclosing the data; the axis limits sit on the
    first and last one."""
    ticks = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(
        values.min(), values.max())
    return ticks[(ticks <= values.min()).nonzero()[0][-1]:
                 (ticks >= values.max()).nonzero()[0][0] + 1]


def draw_density(ax, x, y, xlim, ylim, color, n_grid=300):
    kde = stats.gaussian_kde(np.vstack([x, y]))
    GX, GY = np.meshgrid(np.linspace(*xlim, n_grid), np.linspace(*ylim, n_grid))
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    # the lowest band (0 to the first level) would fill the whole panel; it
    # is left unfilled so the background (and panel A's market band) stays
    # visible, the other bands keep their ramp colours
    levels = np.linspace(0, Z.max(), N_DENSITY_LEVELS + 1)
    cmap = LinearSegmentedColormap.from_list('density', density_ramp(color))
    colors = ['none'] + [cmap(i/(N_DENSITY_LEVELS - 1))
                         for i in range(1, N_DENSITY_LEVELS)]
    ax.contourf(GX, GY, Z, levels=levels, colors=colors, antialiased=True, zorder=2)


def pareto_frontier(x, y, sense=('min', 'min')):
    """The non-dominated samples under `sense` ('min' / 'max' for x, then y),
    ordered from the best x to the best y."""
    sx, sy = (1 if s == 'min' else -1 for s in sense)
    x_, y_ = sx*x, sy*y # both minimized
    order = np.lexsort((y_, x_)) # by x, ties by y
    front, best = [], np.inf
    for i in order:
        if y_[i] < best:
            front.append(i)
            best = y_[i]
    front = np.array(front)
    return x[front], y[front]


def draw_pareto_frontier(ax, x, y, sense):
    # staircase: between two consecutive frontier samples the best attainable
    # y is the first one's, so step along x first, then along y
    fx, fy = pareto_frontier(x, y, sense)
    ax.step(fx, fy, where='post', color=PARETO_COLOR, lw=PARETO_LW, ls=PARETO_LS,
            zorder=4)
    return fx, fy


def binned_medians(x, y, n_bins=N_TREND_BINS):
    """Median x and median y of each of `n_bins` equal-count bins of x (split
    at its quantiles), in increasing x."""
    edges = np.quantile(x, np.linspace(0, 1, n_bins + 1))
    bins = np.digitize(x, edges[1:-1]) # 0 .. n_bins - 1
    bx = np.array([np.median(x[bins == i]) for i in range(n_bins)])
    by = np.array([np.median(y[bins == i]) for i in range(n_bins)])
    return bx, by


def draw_binned_medians(ax, x, y):
    bx, by = binned_medians(x, y)
    ax.plot(bx, by, color=TREND_COLOR, lw=TREND_LW, ls=TREND_LS, zorder=4)
    return bx, by


def draw_box(ax, values, orientation, color):
    edge = dark_shade(color)
    lo_w, hi_w = BOX_PERCENTILES['whis']
    ax.boxplot(values, whis=[lo_w, hi_w], orientation=orientation,
               widths=BOX_WIDTH, showfliers=False, patch_artist=True,
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
    gs = cell.subgridspec(2, 2, width_ratios=(4.2, MARGINAL_RATIO), height_ratios=(MARGINAL_RATIO, 4.2),
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
    if panel.get('box'):
        (x0, x1), (y0, y1) = panel['box']
        ax.fill_between((x0, x1), y0, y1, color=MARKET_BAND_COLOR, lw=0, zorder=0)
    draw_density(ax, x, y, (xticks[0], xticks[-1]), (yticks[0], yticks[-1]),
                 panel['color'])
    fx = fy = None
    if panel.get('pareto'):
        fx, fy = draw_pareto_frontier(ax, x, y, panel['pareto'])
        print(f'  Pareto frontier: {fx.size} non-dominated samples, '
              f"{panel['x']} {fx[0]:.4g}-{fx[-1]:.4g}, {panel['y']} {fy[0]:.4g}-{fy[-1]:.4g}")
    if panel.get('trend'):
        bx, by = draw_binned_medians(ax, x, y)
        print(f'  binned medians ({bx.size} equal-count bins): '
              f"{panel['x']} {bx[0]:.4g} -> {bx[-1]:.4g}, "
              f"{panel['y']} {by[0]:.4g} -> {by[-1]:.4g} ({by[-1] - by[0]:+.3f})")
    if base:
        ax.plot(base[panel['x']], base[panel['y']], 'D', ms=8, mfc='w', mec=INK,
                mew=1.2, zorder=5)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_ylim(yticks[0], yticks[-1])
    ax.set_xlabel(bold_title(panel['xlabel']), fontsize=FONTS['axis_title'])
    ax.set_ylabel(bold_title(panel['ylabel']), fontsize=FONTS['axis_title'])
    ax.xaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(N_MINOR_PER_MAJOR + 1))

    if panel.get('callouts'):
        draw_panel_callouts(ax, panel, base, fx, fy)

    draw_box(ax_top, x, 'horizontal', panel['color'])
    draw_box(ax_right, y, 'vertical', panel['color'])
    ax_top.set_ylim(0.4, 1.6)
    ax_right.set_xlim(0.4, 1.6)
    # panel letter just above the y-axis title, level with the top box
    ax.text(PANEL_LETTER_XY[0], PANEL_LETTER_XY[1], letter, transform=ax.transAxes,
            ha='left', va='top',
                fontsize=FONTS['panel_letter'], fontweight='bold')
    return ax


def _nudged(ax, dx, dy):
    """Data coordinates shifted by (dx, dy) points."""
    return ax.transData + ScaledTranslation(dx/72, dy/72, ax.figure.dpi_scale_trans)


def _callout(ax, text, xy, offset, rad, color, shrinkB=0, tip=(0, 0), **kw):
    return ax.annotate(text, xy=xy, xycoords=_nudged(ax, *tip), xytext=offset,
                       textcoords='offset points',
                       fontsize=CALLOUT_FONTSIZE, fontweight='bold', color=color,
                       arrowprops=dict(**CALLOUT_ARROW, color=color, shrinkB=shrinkB,
                                       connectionstyle=f'arc3,rad={rad}'),
                       zorder=6, **kw)


def draw_panel_callouts(ax, panel, base, fx, fy):
    """Panel-A callouts: the baseline diamond, the Pareto frontier, the grey
    conventional-facility box and the gasoline price range between its two
    dashed lines."""
    if base:
        c = BASELINE_CALLOUT
        _callout(ax, c['text'], (base[panel['x']], base[panel['y']]), c['offset'],
                 c['rad'], INK, shrinkB=c['shrinkB'], ha='left', va='bottom')
    if fx is not None:
        c = PARETO_CALLOUT
        # a point on the staircase, a fraction `at` along its x span
        xs = fx[0] + c['at']*(fx[-1] - fx[0])
        ys = fy[np.searchsorted(fx, xs, side='right') - 1]
        _callout(ax, c['text'], (xs, ys), c['offset'], c['rad'], PARETO_COLOR, tip=c['tip'],
                 ha='left', va='top')
    if panel.get('market'):
        c = CONVENTIONAL_CALLOUT
        cx, cy = c['at']
        _callout(ax, c['text'], (cx, ETHANOL_MARKET_RANGE[1] if cy is None else cy),
                 c['offset'], c['rad'], GASOLINE_LINE_COLOR, tip=c['tip'],
                 ha='right', va='bottom')
        c = GASOLINE_CALLOUT
        lo, hi = GASOLINE_PRICE_RANGE
        ax.annotate('', xy=(c['x'], hi), xycoords=_nudged(ax, 0, c['tip']),
                    xytext=(c['x'], lo), textcoords=_nudged(ax, 0, -c['tip']),
                    arrowprops=dict(**{**CALLOUT_ARROW, 'arrowstyle': '<|-|>', 'shrinkA': 0},
                                    shrinkB=0, color=GASOLINE_LINE_COLOR), zorder=2)
        ax.text(c['x'], 0.5*(lo + hi), c['text'], rotation=90, ha='center', va='center',
                fontsize=CALLOUT_FONTSIZE, fontweight='bold', color=GASOLINE_LINE_COLOR,
                bbox=dict(facecolor='white', edgecolor='none', pad=2), zorder=3)


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
    width, height = CELL_SIZE[0]*n_cols, CELL_SIZE[1]*n_rows
    fig = plt.figure(figsize=(width, height))
    m = MARGINS_IN
    grid = fig.add_gridspec(n_rows, n_cols, wspace=GRID_SPACE['w'], hspace=GRID_SPACE['h'],
                            left=m['left']/width, right=1 - m['right']/width,
                            bottom=m['bottom']/height, top=1 - m['top']/height)
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

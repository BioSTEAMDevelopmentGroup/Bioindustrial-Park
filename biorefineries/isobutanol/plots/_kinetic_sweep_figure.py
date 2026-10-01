#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Shared engine of the single-panel 2-D kinetic-sweep MESP figures.

Draws one `evaluate_*` sweep as the minimum ethanol selling price (MESP = the
sweep's purity-adjusted ethanol MPSP at a 15 % IRR in $/GGE) with the grid
optimum of several sweep metrics marked and annotated, in the format and
style of plots/plot_feeding_strategy_baseline_A.py (panel A), whose helpers
are loaded by file path. Each figure script fills a `SweepFigure` and calls
`main(figure)`.

Sim-safe: reads the sweep CSVs only and never imports the package. Load this
module by file path (plots/ is not a package).
"""

import importlib.util
import os
from dataclasses import dataclass, field

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

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_feeding_strategy_style', os.path.join(HERE, 'plot_feeding_strategy_baseline_A.py'))
fs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fs)

FONTS = fs.FONTS
TICK_LEN = fs.TICK_LEN
USD_PER_KG_TO_USD_PER_GGE = fs.USD_PER_KG_TO_USD_PER_GGE

RESULTS_DIR = os.path.join(os.path.dirname(HERE), 'analyses', 'results')
OUTPUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Kinetic-sweeps')

# Default colour scale (MESP, $/GGE): the feeding-strategy figure's levels, so
# the figures read on one scale; values past the top take the over-colour. A
# figure may override all three (SweepFigure.levels / cbar_ticks /
# cbar_minor_step); the colormap is sized to the levels (at least 90 colours).
MESP_LEVELS = fs.MESP_LEVELS
MESP_CBAR_TICKS = fs.MESP_CBAR_TICKS
MESP_CBAR_MINOR_STEP = fs.MESP_CBAR_MINOR_STEP

# Default IRR colour scale [%]: -15 to 15 % (the sweeps top out below the 15 %
# hurdle; the finite negative IRRs are a thin rim down to ~-75 %, and a scale
# reaching that far would squeeze the profitable range into a sixth of the
# colours); the colormap is the MESP one reversed (high IRR = yellow, as low
# MESP). Cells below -15 % and the -inf of an outright money-loser (solve_TEA
# finds no IRR at all) take the dark-grey under-colour; a white contour marks
# break-even (IRR = 0).
IRR_LEVELS = np.arange(-15., 15.00001, 0.5)
IRR_CBAR_TICKS = np.arange(-15., 15.00001, 5.)
IRR_CBAR_MINOR_STEP = 1.
# = the MESP figures' over-colour: the colormap's dark-grey end
IRR_UNDER_COLOR = fs.JBEI_UCB_colormap()(1.)
IRR_REFERENCE_LINES = (0.,)

# Colour axes: the sweep metric read, its unit conversion, default scale,
# which end extends, the colorbar title, its tick format and the values
# marked by a white contour line (and a matching line on the colourbar)
COLOR_AXES = {
    'MESP': dict(sweep_metric='MPSP', scale=USD_PER_KG_TO_USD_PER_GGE,
                 levels=MESP_LEVELS, cbar_ticks=MESP_CBAR_TICKS,
                 cbar_minor_step=MESP_CBAR_MINOR_STEP, extend='max',
                 reverse=False, under_color=None,
                 title=fs.bold_title('MESP', fs.USD_PER_GGE), tick_fmt='{:.2f}',
                 reference_lines=()),
    'IRR': dict(sweep_metric='IRR', scale=100.,
                levels=IRR_LEVELS, cbar_ticks=IRR_CBAR_TICKS,
                cbar_minor_step=IRR_CBAR_MINOR_STEP, extend='min',
                reverse=True, under_color=IRR_UNDER_COLOR,
                title=fs.bold_title('IRR', '%'), tick_fmt='{:.0f}',
                reference_lines=IRR_REFERENCE_LINES),
}

# Ethanol market price range, the one plots/plot_uncertainty_MPSP_vs_TCI.py
# (ETHANOL_MARKET_RANGE) and uncertainties_IBO_EtOH.py use: Jan 2021 - Dec 2025
# five-year low and high, 1.5475 and 3.4500 $/gal, from
# https://tradingeconomics.com/commodity/ethanol, / 0.67 GGE/gal -> $/GGE
ETHANOL_MARKET_RANGE_PER_GAL = (1.5475, 3.4500) # $/gal
ETHANOL_MARKET_RANGE = tuple(v / fs.GGE_PER_GAL for v in ETHANOL_MARKET_RANGE_PER_GAL) # $/GGE

# Comparison-range style (contourplots' comparison_range convention: white
# '///' hatch over the in-range region + a white bounding contour, the same
# hatch on the colourbar)
COMPARISON_HATCH = '///'
COMPARISON_HATCH_RGBA = (1., 1., 1., 0.45)
COMPARISON_LINE_RGBA = (1., 1., 1., 0.9)
COMPARISON_LINE_WIDTH = 0.9
COMPARISON_HATCH_LINE_WIDTH = 0.6

BASELINE_MARKER_SIZE = 7 # pt

G_PER_L_PER_H =r'$\mathrm{g·L}^{-1}\mathrm{·h}^{-1}$'


@dataclass
class SweepFigure:
    """One sweep figure. `csv_prefix` is the sweep's `file_to_save` prefix
    (the loader appends '_<metric>.csv'); `spec_1` / `spec_2` are its x / y
    linspaces. `optima` rows: (sweep metric, 'min'/'max', label, marker, face
    colour, size [pt], label offset from the marker [pt], arrow curvature); a
    label is left-/right-aligned by the sign of its x offset unless
    `label_ha` overrides it; `marker_nudge` fans co-located optima apart
    [pt]. `color_axis` picks the coloured metric ('MESP' or 'IRR', see
    COLOR_AXES); `levels` / `cbar_ticks` / `cbar_minor_step` override its
    default colour scale; `comparison_range` = (low, high), in the colour
    axis's units, hatches the region whose value lies in that range (e.g.
    ETHANOL_MARKET_RANGE on MESP);
    `baseline` = (x, y) marks the baseline point with a white diamond, and
    `baseline_callout` = (label, label offset [pt], arrow curvature) labels it
    in the optima's callout style (None = no label)."""
    output_stem: str
    csv_prefix: str
    spec_1: np.ndarray
    spec_2: np.ndarray
    xlabel: str
    ylabel: str
    optima: list
    xlim: tuple = None
    ylim: tuple = None
    x_major: float = 1.
    x_minor_div: int = 2
    y_major: float = 1.
    y_minor_div: int = 2
    label_ha: dict = field(default_factory=dict)
    label_color: dict = field(default_factory=dict)
    marker_nudge: dict = field(default_factory=dict)
    figsize: tuple = (5.6, 4.2)
    color_axis: str = 'MESP'
    levels: np.ndarray = None
    cbar_ticks: np.ndarray = None
    cbar_minor_step: float = None
    comparison_range: tuple = None
    baseline: tuple = None
    baseline_callout: tuple = None

#%% Sweep data

def load_sweep_metric(figure, metric):
    """(n_y, n_x) array: row = spec_2, column = spec_1."""
    path = os.path.join(RESULTS_DIR, f'{figure.csv_prefix}_{metric}.csv')
    arr = pd.read_csv(path).iloc[:, 1:].to_numpy(dtype=float)
    if arr.shape != (len(figure.spec_2), len(figure.spec_1)):
        raise ValueError(f'{path}: shape {arr.shape} does not match the '
                         f'{len(figure.spec_1)} x {len(figure.spec_2)} grid.')
    return arr


def grid_optimum(figure, arr, sense):
    finite = np.isfinite(arr)
    value = np.min(arr[finite]) if sense == 'min' else np.max(arr[finite])
    i, j = np.argwhere(arr == value)[0]
    return figure.spec_1[j], figure.spec_2[i], value


def fill_failed_cells(figure, arr):
    """Replace each ISOLATED non-finite cell (a failed sweep point: every
    in-grid 4-neighbour finite) by the mean of its neighbours, so the filled
    contours have no pin-hole. Contiguous non-finite regions (e.g. a
    no-ethanol edge where the MESP is undefined) are left blank."""
    out = arr.copy()
    n_blank = 0
    for i, j in np.argwhere(~np.isfinite(arr)):
        neighbours = [arr[a, b] for a, b in ((i-1, j), (i+1, j), (i, j-1), (i, j+1))
                      if 0 <= a < arr.shape[0] and 0 <= b < arr.shape[1]]
        if all(np.isfinite(neighbours)):
            out[i, j] = np.mean(neighbours)
            print(f'  WARNING: failed sweep point at x = {figure.spec_1[j]:.4g}, '
                  f'y = {figure.spec_2[i]:.4g} -- filled from its neighbours.')
        else:
            n_blank += 1
    if n_blank:
        print(f'  NOTE: {n_blank} non-finite cell(s) in contiguous regions left blank.')
    return out

#%% Figure

def draw_panel(figure, fig, ax, cax):
    axis = COLOR_AXES[figure.color_axis]
    levels = axis['levels'] if figure.levels is None else figure.levels
    cbar_ticks = axis['cbar_ticks'] if figure.cbar_ticks is None else figure.cbar_ticks
    minor_step = figure.cbar_minor_step or axis['cbar_minor_step']
    extend = axis['extend']
    values = load_sweep_metric(figure, axis['sweep_metric']) * axis['scale']
    if axis['under_color'] is not None:
        # contourf masks non-finite cells: push -inf (no IRR, money-losing)
        # below the lowest level so it takes the under-colour, not a blank
        values[np.isneginf(values)] = levels[0] - (levels[1] - levels[0])
    values = fill_failed_cells(figure, values)
    # at least one colour per band plus the extend colour (BoundaryNorm with
    # extend); the default 90 unless a figure's levels need more
    cmap = fs.JBEI_UCB_colormap(max(90, len(levels)))
    if axis['reverse']:
        cmap = cmap.reversed()
    if axis['under_color'] is not None:
        cmap.set_under(axis['under_color'])
    norm = BoundaryNorm(levels, cmap.N, extend=extend)
    cs = ax.contourf(figure.spec_1, figure.spec_2, values, levels=levels,
                     cmap=cmap, norm=norm, extend=extend, zorder=1)
    cs.set_edgecolor('face')  # no hairline seams between bands in the PDF
    ax.set_xlim(*(figure.xlim or (figure.spec_1[0], figure.spec_1[-1])))
    ax.set_ylim(*(figure.ylim or (figure.spec_2[0], figure.spec_2[-1])))

    # comparison range: hatch the in-range region, bound it with a contour
    if figure.comparison_range is not None:
        with matplotlib.rc_context({'hatch.color': COMPARISON_HATCH_RGBA,
                                    'hatch.linewidth': COMPARISON_HATCH_LINE_WIDTH}):
            ax.contourf(figure.spec_1, figure.spec_2, values,
                        levels=list(figure.comparison_range), colors='none',
                        hatches=[COMPARISON_HATCH], zorder=2)
        ax.contour(figure.spec_1, figure.spec_2, values,
                   levels=list(figure.comparison_range),
                   colors=[COMPARISON_LINE_RGBA], linewidths=COMPARISON_LINE_WIDTH,
                   zorder=3)

    def callout(label, xy, where, offset, rad, size):
        ax.annotate(label, xy=xy, xycoords=where, xytext=offset,
                    textcoords='offset points', fontsize=FONTS['annotation'],
                    color=figure.label_color.get(label, 'black'),
                    ha=figure.label_ha.get(label, 'left' if offset[0] >= 0 else 'right'),
                    va='bottom' if offset[1] >= 0 else 'top', zorder=11,
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle='-|>', mutation_scale=9, color='black', lw=0.9,
                                    shrinkA=1, shrinkB=size/2 + 1,
                                    connectionstyle=f'arc3,rad={rad}'))

    # baseline point: white diamond, black outline, + its callout
    if figure.baseline is not None:
        ax.plot(*figure.baseline, linestyle='none', marker='D', markersize=BASELINE_MARKER_SIZE,
                markerfacecolor='white', markeredgecolor='black',
                markeredgewidth=0.8, zorder=9, clip_on=False)
        if figure.baseline_callout is not None:
            label, offset, rad = figure.baseline_callout
            callout(label, figure.baseline, 'data', offset, rad, BASELINE_MARKER_SIZE)

    # reference values (e.g. IRR break-even), in the comparison line's style
    if axis['reference_lines']:
        ax.contour(figure.spec_1, figure.spec_2, values,
                   levels=list(axis['reference_lines']),
                   colors=[COMPARISON_LINE_RGBA], linewidths=COMPARISON_LINE_WIDTH,
                   zorder=3)

    # optimum markers + labels
    optima = {}
    for metric, sense, label, marker, color, size, offset, rad in figure.optima:
        ox, oy, value = grid_optimum(figure, load_sweep_metric(figure, metric), sense)
        if metric == 'MPSP':
            value *= USD_PER_KG_TO_USD_PER_GGE  # MESP, $/GGE
        optima[label] = (ox, oy, value)
        dx, dy = figure.marker_nudge.get(label, (0., 0.))
        where = ax.transData + ScaledTranslation(dx/72, dy/72, fig.dpi_scale_trans)
        ax.plot(ox, oy, linestyle='none', marker=marker, markersize=size,
                markerfacecolor=color, markeredgecolor='black',
                markeredgewidth=0.8, zorder=10, clip_on=False, transform=where)
        callout(label, (ox, oy), where, offset, rad, size)

    ax.xaxis.set_major_locator(MultipleLocator(figure.x_major))
    ax.xaxis.set_minor_locator(AutoMinorLocator(figure.x_minor_div))
    ax.yaxis.set_major_locator(MultipleLocator(figure.y_major))
    ax.yaxis.set_minor_locator(AutoMinorLocator(figure.y_minor_div))
    ax.set_xlabel(figure.xlabel, fontsize=FONTS['axis_title'])
    ax.set_ylabel(figure.ylabel, fontsize=FONTS['axis_title'])

    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cax, spacing='proportional', extend=extend,
                        extendfrac=0.04)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([axis['tick_fmt'].format(t).replace('-', '−')
                         for t in cbar_ticks])
    cbar.ax.yaxis.set_minor_locator(FixedLocator(
        [v for v in np.arange(levels[0], levels[-1] + 1e-9, minor_step)
         if not np.any(np.isclose(v, cbar_ticks))]))
    if figure.comparison_range is not None:
        with matplotlib.rc_context({'hatch.color': COMPARISON_HATCH_RGBA,
                                    'hatch.linewidth': COMPARISON_HATCH_LINE_WIDTH}):
            cbar.ax.fill_betweenx(figure.comparison_range, 0., 1., facecolor='none',
                                  edgecolor=COMPARISON_HATCH_RGBA,
                                  hatch=COMPARISON_HATCH, linewidth=0., zorder=3)
    for v in axis['reference_lines']:
        cbar.ax.axhline(v, color=COMPARISON_LINE_RGBA, linewidth=COMPARISON_LINE_WIDTH,
                        zorder=3)
    cbar.ax.tick_params(which='major', labelsize=FONTS['tick'],
                        length=TICK_LEN['major'])
    cbar.ax.tick_params(which='minor', length=TICK_LEN['minor'])
    cbar.set_label(axis['title'], fontsize=FONTS['axis_title'])
    return optima


def make_figure(figure):
    fs.apply_font_rcparams()
    fig = plt.figure(figsize=figure.figsize)
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.0, 0.045], wspace=0.05,
                  left=0.15, right=0.86, bottom=0.14, top=0.86)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    optima = draw_panel(figure, fig, ax, cax)
    fig.canvas.draw()
    fs.style_ticks(ax)
    return fig, optima


def main(figure):
    fig, optima = make_figure(figure)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(OUTPUT_DIR, f'{figure.output_stem}.{ext}')
        fig.savefig(path, dpi=600, bbox_inches='tight', facecolor='white')
        print(f'Saved {path}')
    plt.close(fig)

    print('\nGrid optima (x, y: value):')
    for label, (ox, oy, value) in optima.items():
        print(f'  {label:14s} {ox:8.4g}, {oy:8.4g}: {value:.5g}')

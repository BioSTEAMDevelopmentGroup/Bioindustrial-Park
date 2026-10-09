#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Objective-landscape figure, "needle on a plateau" layout.

Two side-by-side panels over ONE x-axis quantity: the Euclidean distance, in
the normalized (unit-cube) 12-D `metabolic_split_12d` decision space, of every
recorded trial from the most profitable trial of the pooled campaigns
(ibo_yield campaign #842, PI 0.504).

(a) raw profitability index PI = NPV(15 % hurdle rate)/TCI of every trial vs
    that distance: a needle at d ~ 0 whose binned median collapses by
    d ~ 0.2;
(b) isobutanol yield of the SAME trials vs the same distance: a plateau out to
    d ~ 0.8, then a decay.

Same anchor, same trials in both panels, so sampling bias cannot explain the
different profile shapes. Each panel overlays the binned median and
interquartile range (0.1-wide distance bins, drawn only where a bin holds at
least MIN_BIN_N trials; the sparser bins nearest the anchor are left to the
individual points). Panel (a) marks PI = 0 (break-even at the 15 % hurdle
rate) and the scenario-A strain's PI, and pins every trial below
PI_FLOOR to an overflow strip under the axis floor (never dropped). The x-axis
stops at D_MAX; the number of trials beyond it is printed in the figure. The
neighbour ibo_yield #1882 (d = 0.135: higher IBO yield than #842, negative PI)
is ringed in both panels.

Data: the seven 2026-09-23 seed-350 `ethanol_isobutanol` x
`metabolic_split_12d` GP campaigns (`_rs350`), pooled and de-duplicated by
`plots/_objective_landscape_data.py` (13,622 unique COMPLETE trials). PI is
the raw tracked `'PI'` column, not the campaigns' log-tail objective.

Sim-safe: reads trajectory CSVs through the data layer, which is loaded BY
FILE PATH; never imports the biorefineries package, never load()s, so it can
run alongside any simulation on any numba-cache state.

Usage:
    python plot_objective_landscape_needle_plateau.py [--out-dir DIR]
        [--stem STEM] [--dpi 300]
"""
import os
import argparse
import importlib.util

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D, TICKDOWN, TICKLEFT
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_data_layer():
    path = os.path.join(_HERE, '_objective_landscape_data.py')
    spec = importlib.util.spec_from_file_location('_objective_landscape_data',
                                                  path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


old = _load_data_layer()

# -----------------------------------------------------------------------------
# Typeface, font sizes, tick style (from
# across_feed_strat_multipanel_global_with_optima_summary.py)
# -----------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'panel_label': 12, 'legend': 9,
         'annotation': 9}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out


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

# -----------------------------------------------------------------------------
# Shared encoding (identical across the objective-landscape figures)
# -----------------------------------------------------------------------------
GREY = '#B8B8B8'
ACCENT = '#C0392B'
ALL_KW = dict(s=4, c=GREY, alpha=0.5, linewidths=0, rasterized=True,
              zorder=1)
POS_KW = dict(s=9, c=ACCENT, alpha=0.9, linewidths=0, zorder=3)
STAR_KW = dict(marker='*', markersize=13, markerfacecolor=ACCENT,
               markeredgecolor='black', markeredgewidth=0.8,
               linestyle='none', zorder=6)

# This layout's own lines
MEDIAN_COLOR = '#1F3A5F'
MEDIAN_KW = dict(color=MEDIAN_COLOR, lw=1.4, marker='o', ms=3.0, zorder=4)
IQR_KW = dict(color=MEDIAN_COLOR, alpha=0.16, lw=0, zorder=2)
RING_KW = dict(marker='o', markersize=8, markerfacecolor='none',
               markeredgecolor='black', markeredgewidth=1.0,
               linestyle='none', zorder=5)
CLIP_KW = dict(marker='v', s=9, c='#8C8C8C', alpha=0.6, linewidths=0,
               rasterized=True, zorder=1)
TEXT = '#1A1A1A'
BREAKEVEN_KW = dict(color='black', lw=0.8, zorder=2.5)
PI_A_KW = dict(color='#4D4D4D', lw=0.8, ls=(0, (4, 2)), zorder=2.5)

# -----------------------------------------------------------------------------
# Figure settings
# -----------------------------------------------------------------------------
D_MAX = 1.5          # x-axis end (distance); trials beyond are counted, not drawn
BIN_WIDTH = 0.1      # distance-bin width of the median / IQR overlay
MIN_BIN_N = 10       # a bin is summarized only if it holds this many trials
PI_FLOOR = -3.5      # PI below this is pinned to the overflow strip
PI_STRIP = -3.68     # y of the pinned (overflow) markers
PI_YLIM = (-3.86, 0.66)
YIELD_YLIM = (-0.012, 0.42)
NEIGHBOUR = ('ibo_yield', 1882)   # the annotated neighbour of #842


def minus(s):
    """Typographic minus signs in a formatted number string."""
    return s.replace('-', '−')


def binned_profile(d, y, edges, min_n):
    """Per-bin (centre, n, q25, median, q75); NaN statistics where n < min_n."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (d >= lo) & (d < hi)
        n = int(s.sum())
        if n >= min_n:
            q25, med, q75 = np.quantile(y[s], [0.25, 0.5, 0.75])
        else:
            q25 = med = q75 = np.nan
        rows.append((0.5 * (lo + hi), n, q25, med, q75))
    return np.array(rows, dtype=float)


def draw_profile(ax, prof):
    ok = np.isfinite(prof[:, 3])
    x, q25, med, q75 = prof[ok, 0], prof[ok, 2], prof[ok, 3], prof[ok, 4]
    ax.fill_between(x, q25, q75, **IQR_KW)
    ax.plot(x, med, **MEDIAN_KW)


def inch_axes(fig, left, bottom, width, height):
    W, H = fig.get_size_inches()
    return fig.add_axes([left / W, bottom / H, width / W, height / H])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out-dir', default=os.path.join(
        old.RESULTS_DIR, 'publication', 'Objective-landscape'))
    parser.add_argument('--stem', default='objective_landscape_needle_plateau')
    parser.add_argument('--dpi', type=int, default=300)
    args = parser.parse_args()

    # ---- data ---------------------------------------------------------------
    data = old.load_landscape()
    f = data.frame
    a = data.anchor
    dist = data.distance_from(a)
    PI = f['PI'].to_numpy(float)
    yield_col, yield_label, yield_unit = old.OBJECTIVES[old.PROCESS_EXAMPLE]
    Y = f[yield_col].to_numpy(float)
    n_all = len(f)
    pos = PI > 0
    n_pos = int(pos.sum())
    anchor_name = f"#{int(f['trial_number'][a])}"
    PI_anchor = PI[a]

    nb = np.flatnonzero((f['campaign'] == NEIGHBOUR[0])
                        & (f['trial_number'] == NEIGHBOUR[1]))
    if len(nb) != 1:
        raise ValueError(f'neighbour {NEIGHBOUR} not found uniquely')
    nb = int(nb[0])

    shown = dist <= D_MAX
    n_shown = int(shown.sum())
    n_beyond = n_all - n_shown
    n_pos_beyond = int((pos & ~shown).sum())
    PI_max_beyond = float(PI[~shown].max())
    clipped = shown & (PI < PI_FLOOR)
    n_clipped = int(clipped.sum())

    edges = np.round(np.arange(0.0, D_MAX + BIN_WIDTH / 2, BIN_WIDTH), 10)
    prof_PI = binned_profile(dist, PI, edges, MIN_BIN_N)
    prof_Y = binned_profile(dist, Y, edges, MIN_BIN_N)

    # ---- figure -------------------------------------------------------------
    apply_font_rcparams()
    fig = plt.figure(figsize=(7.0, 4.06))
    ax_a = inch_axes(fig, 0.66, 1.48, 2.72, 2.40)
    ax_b = inch_axes(fig, 4.18, 1.48, 2.72, 2.40)

    for ax, y, prof in ((ax_a, np.where(clipped, np.nan, PI), prof_PI),
                        (ax_b, Y, prof_Y)):
        s = shown & np.isfinite(y)
        ax.scatter(dist[s], y[s], **ALL_KW)
        draw_profile(ax, prof)
        sp = s & pos
        ax.scatter(dist[sp], y[sp], **POS_KW)
        ax.plot(dist[a], y[a], **STAR_KW)
        ax.plot(dist[nb], y[nb], **RING_KW)
        ax.set_xlim(-0.06, D_MAX + 0.02)   # left margin: room for the star
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    # (a) overflow strip, reference lines
    ax_a.scatter(dist[clipped], np.full(n_clipped, PI_STRIP), **CLIP_KW)
    ax_a.axhline(PI_FLOOR, color='#9A9A9A', lw=0.5, ls=(0, (1, 1.5)), zorder=0)
    # Labelled in the legend: every in-plot spot next to these two lines, at
    # any distance, is occupied by PI > 0 trials (above) or the IQR band and
    # trials (below).
    ax_a.axhline(0.0, **BREAKEVEN_KW)
    ax_a.axhline(old.PI_A, **PI_A_KW)
    ax_a.set_ylim(*PI_YLIM)
    ax_a.yaxis.set_major_locator(MultipleLocator(1.0))
    ax_a.set_ylabel('Profitability index (PI)', fontsize=FONTS['axis_title'])

    # (b)
    ax_b.set_ylim(*YIELD_YLIM)
    ax_b.yaxis.set_major_locator(MultipleLocator(0.1))
    ax_b.set_ylabel(f'{yield_label} ({yield_unit})',
                    fontsize=FONTS['axis_title'])

    # neighbour annotations
    arrow = dict(arrowstyle='-', color=TEXT, lw=0.6, shrinkA=1, shrinkB=4.5)
    ax_a.annotate(f'#{NEIGHBOUR[1]}: PI {minus(f"{PI[nb]:.2f}")}\n'
                  f'{anchor_name}: PI {PI_anchor:.2f}',
                  xy=(dist[nb], PI[nb]), xytext=(0.15, -1.75),
                  ha='left', va='center', fontsize=FONTS['annotation'],
                  color=TEXT, arrowprops=arrow, zorder=7)
    ax_b.annotate(f'#{NEIGHBOUR[1]}: yield {Y[nb]:.3f}\n'
                  f'{anchor_name}: yield {Y[a]:.3f}',
                  xy=(dist[nb], Y[nb]), xytext=(0.15, 0.075),
                  ha='left', va='center', fontsize=FONTS['annotation'],
                  color=TEXT, arrowprops=arrow, zorder=7)

    # panel labels, shared x label, beyond-range note
    for ax, lab in ((ax_a, 'a'), (ax_b, 'b')):
        ax.text(-0.02, 1.02, lab, transform=ax.transAxes, ha='right',
                va='bottom', fontsize=FONTS['panel_label'], fontweight='bold')
    W, H = fig.get_size_inches()
    fig.text(0.5, 1.0 / H,
             f'Distance from the most profitable trial ({anchor_name}), '
             'normalized 12-D decision space',
             ha='center', va='bottom', fontsize=FONTS['axis_title'])

    # legend
    handles = [
        Line2D([], [], marker='o', ls='none', ms=3.2, mfc=GREY, mec='none',
               alpha=0.8, label=f'All trials (n = {n_all:,})'),
        Line2D([], [], marker='o', ls='none', ms=4.2, mfc=ACCENT, mec='none',
               label=f'PI > 0 (n = {n_pos})'),
        Line2D([], [], label=f'Most profitable trial ({anchor_name}, '
               f'PI {PI_anchor:.2f})', **STAR_KW),
        Line2D([], [], label=f'Neighbour #{NEIGHBOUR[1]} '
               f'(d = {dist[nb]:.2f})', **RING_KW),
        Line2D([], [], color=MEDIAN_COLOR, lw=1.4, marker='o', ms=3.0,
               label=f'Binned median ({BIN_WIDTH:g}-wide bins, n ≥ {MIN_BIN_N})'),
        Patch(fc=MEDIAN_COLOR, alpha=0.16, ec='none',
              label='Interquartile range'),
        Line2D([], [], marker='v', ls='none', ms=3.6, mfc='#8C8C8C',
               mec='none', alpha=0.8,
               label=f'PI < {minus(f"{PI_FLOOR:g}")}, drawn below the '
                     f'dotted floor (n = {n_clipped})'),
        Line2D([], [], label='Break-even: PI = 0 (NPV = 0 at a 15 % '
               'hurdle rate)', **BREAKEVEN_KW),
        Line2D([], [], label=f'Scenario-A strain (PI = '
               f'{minus(f"{old.PI_A:.2f}")})', **PI_A_KW),
        # the x-range note, as a handle-less last entry of the right column
        Line2D([], [], ls='none', label=f'Not shown: {n_beyond:,} trials at '
               f'd > {D_MAX:g} ({n_pos_beyond} with PI > 0)'),
    ]
    # column-major fill: trial markers left, summary / reference lines right
    order = [0, 1, 2, 3, 6, 4, 5, 7, 8, 9]
    fig.legend(handles=[handles[i] for i in order], loc='lower center',
               ncol=2, bbox_to_anchor=(0.5, 0.03 / H), fontsize=FONTS['legend'],
               frameon=False, handlelength=2.2, columnspacing=1.6,
               handletextpad=0.5, labelspacing=0.35)

    fig.canvas.draw()
    for ax in (ax_a, ax_b):
        style_ticks(ax)

    os.makedirs(args.out_dir, exist_ok=True)
    out_png = os.path.join(args.out_dir, args.stem + '.png')
    out_pdf = os.path.join(args.out_dir, args.stem + '.pdf')
    fig.savefig(out_png, dpi=args.dpi)
    fig.savefig(out_pdf, dpi=args.dpi)
    plt.close(fig)

    # ---- console summary ----------------------------------------------------
    print(f'{data.n_raw} COMPLETE rows -> {n_all} unique trials; '
          f'anchor {data.describe(a)}: PI {PI_anchor:.4f}, '
          f'{yield_col} {Y[a]:.4f}')
    print(f'PI > 0: {n_pos} ({n_pos - n_pos_beyond} shown); PI > 0.3: '
          f'{int((PI > 0.3).sum())}; PI_A = {old.PI_A}')
    print(f'shown d <= {D_MAX:g}: {n_shown}; beyond: {n_beyond} '
          f'({n_pos_beyond} with PI > 0, max PI {PI_max_beyond:.3f})')
    print(f'PI < {PI_FLOOR:g} pinned to the overflow strip: {n_clipped} '
          f'(min PI at d <= {D_MAX:g}: {PI[shown].min():.3f})')
    print(f'neighbour {NEIGHBOUR[0]} #{NEIGHBOUR[1]}: d {dist[nb]:.3f}, '
          f'PI {PI[nb]:.3f}, {yield_col} {Y[nb]:.3f}, IBO titer '
          f"{f['IBO titer'][nb]:.1f} vs {f['IBO titer'][a]:.1f} g/L, tau "
          f"{f['tau'][nb]:.1f} vs {f['tau'][a]:.1f} h")
    print(f'\n{"bin":>9} {"n":>5} | {"PI q25":>7} {"med":>7} {"q75":>7} | '
          f'{"Y q25":>6} {"med":>6} {"q75":>6}')
    for (c, n, p25, pm, p75), (_, _, y25, ym, y75), lo in zip(
            prof_PI, prof_Y, edges[:-1]):
        tag = '' if np.isfinite(pm) else '  (n < %d: points only)' % MIN_BIN_N
        print(f'{lo:4.1f}-{lo + BIN_WIDTH:3.1f} {int(n):5d} | {p25:7.3f} '
              f'{pm:7.3f} {p75:7.3f} | {y25:6.3f} {ym:6.3f} {y75:6.3f}{tag}')
    print(f'\nsaved {out_png}\nsaved {out_pdf}')


if __name__ == '__main__':
    main()

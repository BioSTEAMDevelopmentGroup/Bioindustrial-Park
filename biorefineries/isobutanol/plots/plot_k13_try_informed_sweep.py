#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Profitability index and isobutanol yield along a 1-D k_13 sweep through
the TRY-informed profitability campaign's optimum.

One panel over a log k_13 (ALS capacity) axis: the profitability index PI
(left y-axis) and the isobutanol yield (right y-axis) of every sweep point,
all other decision variables held at the TRY-informed (relay) campaign's best
trial #1912 (k_13 = 4.0 g/L/h, the top of the campaign band). Marked: the
trial itself, PI = 0 (break-even at the 15 % hurdle rate), the region beyond
the campaign's k_13 band, and the k_13 above which the enzyme burden derates
growth. Points that did not simulate (INFEASIBLE / ERROR / LOST) leave gaps
and are counted in the legend. The console prints a roughness summary of both
curves along the sweep (normalized total variation, number of local maxima).

Data: results/evaluate_k13_try_informed_optimum.csv + _anchor.json, written by
analyses/evaluate_k13_try_informed_optimum.py (the simulation stage, ask-first).
Sim-safe: reads the CSV / JSON only; never imports the biorefineries package.

Usage:
    python plot_k13_try_informed_sweep.py [--csv PATH] [--out-dir DIR]
        [--stem STEM] [--dpi 300]
"""
import os
import json
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D, TICKDOWN
from matplotlib.patches import Patch

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
SWEEP_STEM = 'evaluate_k13_try_informed_optimum'

# -----------------------------------------------------------------------------
# Typeface, font sizes, tick style (the package's publication conventions)
# -----------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'legend': 9, 'annotation': 9}
TICK_LEN = {'major': 4.0, 'minor': 2.0}   # pt; labelled sides extend in AND out


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def style_twin_ticks(ax, ax_right):
    """Bottom / left / right (the three labelled sides) in and out, top
    inward only; major and minor. Call after fig.canvas.draw()."""
    for which, L in TICK_LEN.items():
        ax.tick_params(axis='x', which=which, direction='inout', length=2*L,
                       top=True, labelsize=FONTS['tick'])
        ax.tick_params(axis='y', which=which, direction='inout', length=2*L,
                       left=True, right=False, labelsize=FONTS['tick'],
                       color=PI_COLOR)
        ax_right.tick_params(axis='y', which=which, direction='inout',
                             length=2*L, left=False, right=True,
                             labelsize=FONTS['tick'], color=YIELD_COLOR)
        get = 'get_major_ticks' if which == 'major' else 'get_minor_ticks'
        for tick in getattr(ax.xaxis, get)():
            tick.tick2line.set_marker(TICKDOWN)
            tick.tick2line.set_markersize(L)

# -----------------------------------------------------------------------------
# Encoding
# -----------------------------------------------------------------------------
PI_COLOR = '#C0392B'       # the objective-landscape figures' accent
YIELD_COLOR = '#1F5FA8'
TEXT = '#1A1A1A'
PI_KW = dict(color=PI_COLOR, lw=1.4, marker='o', ms=2.0, zorder=4)
YIELD_KW = dict(color=YIELD_COLOR, lw=1.4, marker='s', ms=1.7, zorder=3)
STAR_KW = dict(marker='*', markersize=13, markerfacecolor=PI_COLOR,
               markeredgecolor='black', markeredgewidth=0.8, linestyle='none',
               zorder=6)
RING_KW = dict(marker='o', markersize=8, markerfacecolor='none',
               markeredgecolor=YIELD_COLOR, markeredgewidth=1.3,
               linestyle='none', zorder=6)
BAND_KW = dict(facecolor='#EDEDED', edgecolor='none', zorder=0)
DERATE_KW = dict(color='#4D4D4D', lw=0.9, ls=(0, (2, 2)), zorder=1)
BREAKEVEN_KW = dict(color='black', lw=0.8, zorder=2)


def minus(s):
    return s.replace('-', '−')


def roughness(x, y):
    """(normalized total variation, n strict interior local maxima) of y
    along the sweep, over its finite points. TV/(max-min) = 1 for a
    monotone curve, 2 for a single clean peak; spikes add to it."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    span = y.max() - y.min()
    tv = np.abs(np.diff(y)).sum()/span if span > 0 else np.nan
    n_max = int(np.sum((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:])))
    return tv, n_max


def load(csv_path):
    frame = pd.read_csv(csv_path).sort_values('k_13').reset_index(drop=True)
    anchor_path = csv_path[:-len('.csv')] + '_anchor.json'
    with open(anchor_path) as fh:
        anchor = json.load(fh)
    return frame, anchor


def plot(frame, anchor, out_stem, dpi):
    apply_font_rcparams()
    ok = frame['state'] == 'OK'
    k = frame['k_13'].to_numpy(float)
    PI = np.where(ok, frame['PI'], np.nan)
    Y = np.where(ok, frame['IBO yield'], np.nan)
    d = frame['burden_factor'].to_numpy(float)
    k_opt = float(frame.loc[frame['is_anchor'] == 1, 'k_13'].iloc[0])
    i_opt = int(np.flatnonzero(k == k_opt)[0])
    band_high = float(anchor['grid']['band_high'])
    derated = np.flatnonzero(np.isfinite(d) & (d < 1 - 1e-9))
    k_derate = k[derated[0]] if derated.size else None

    fig, ax = plt.subplots(figsize=(4.9, 4.8))
    fig.subplots_adjust(left=0.17, right=0.83, top=0.96, bottom=0.42)
    ax_r = ax.twinx()
    ax.set_zorder(ax_r.get_zorder() + 1)   # PI (and its star) on top
    ax.patch.set_visible(False)

    x_lo, x_hi = k.min()/1.25, k.max()*1.25
    # on the BOTTOM (right-axis) layer, so it never covers either curve
    ax_r.axvspan(band_high, x_hi, **BAND_KW)
    if k_derate is not None:
        ax.axvline(k_derate, **DERATE_KW)
    ax.axhline(0.0, **BREAKEVEN_KW)

    ax_r.plot(k, Y, **YIELD_KW)
    ax_r.plot([k_opt], [Y[i_opt]], **RING_KW)
    ax.plot(k, PI, **PI_KW)
    ax.plot([k_opt], [PI[i_opt]], **STAR_KW)

    ax.set_xscale('log')
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel(r'ALS capacity, $k_{13}$ ($\mathrm{g·L}^{-1}·\mathrm{h}^{-1}$)',
                  fontsize=FONTS['axis_title'])
    ax.set_ylabel('Profitability index (PI)', fontsize=FONTS['axis_title'])
    ax_r.set_ylabel(r'Isobutanol yield ($\mathrm{g·g}^{-1}$)',
                    fontsize=FONTS['axis_title'], rotation=270, labelpad=16)
    y_top = np.nanmax(Y)
    ax_r.set_ylim(-0.03*y_top, 1.12*y_top)
    pi_lo, pi_hi = np.nanmin(PI), np.nanmax(PI)
    pad = 0.08*(pi_hi - pi_lo)
    ax.set_ylim(min(pi_lo - pad, -pad), pi_hi + 1.5*pad)
    ax.spines['left'].set_color(PI_COLOR)
    ax_r.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax_r.spines['right'].set_color(YIELD_COLOR)
    for side in ('left', 'right'):
        ax.spines[side].set_linewidth(1.2)
        ax_r.spines[side].set_linewidth(1.2)
    fig.canvas.draw()
    style_twin_ticks(ax, ax_r)
    for label in ax.get_yticklabels() + ax_r.get_yticklabels():
        label.set_text(minus(label.get_text()))

    n_bad = int((~ok).sum())
    handles = [
        Line2D([], [], **{k_: v for k_, v in PI_KW.items() if k_ != 'zorder'},
               label='Profitability index (left axis)'),
        Line2D([], [], **{k_: v for k_, v in YIELD_KW.items() if k_ != 'zorder'},
               label='Isobutanol yield (right axis)'),
        Line2D([], [], **{k_: v for k_, v in STAR_KW.items() if k_ != 'zorder'},
               label=f"TRY-informed optimum (#{anchor['trial_number']}, "
                     f"PI {PI[i_opt]:.2f})"),
        Line2D([], [], **BREAKEVEN_KW, label='PI = 0 (break-even at a '
                                             '15 % hurdle rate)'),
        Patch(**BAND_KW, label=f'Beyond the campaign band '
                               f'($k_{{13}}$ > {band_high:g})'),
    ]
    if k_derate is not None:
        handles.append(Line2D([], [], **DERATE_KW,
                              label=f'Growth derated by enzyme burden '
                                    f'($k_{{13}}$ ≥ {k_derate:.2g})'))
    if n_bad:
        states = frame.loc[~ok, 'state'].value_counts().to_dict()
        handles.append(Line2D([], [], color='none', label='Not simulated: ' +
                              ', '.join(f'{v} {s.lower()}'
                                        for s, v in states.items())))
    fig.legend(handles=handles, loc='lower center', ncol=1, frameon=False,
               fontsize=FONTS['legend'], bbox_to_anchor=(0.5, 0.0),
               handlelength=2.2)

    os.makedirs(os.path.dirname(out_stem), exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out_stem}.{ext}', dpi=dpi)
    plt.close(fig)

    # --- console summary ---
    print(f"{len(frame)} points ({int(ok.sum())} simulated): k_13 "
          f"{k.min():g}-{k.max():g}")
    print(f"anchor #{anchor['trial_number']}: reproduced PI "
          f"{anchor['reproduced_PI']:.5f} vs recorded {anchor['recorded_PI']:.5f}")
    i_max = int(np.nanargmax(PI))
    print(f"PI max {PI[i_max]:.4f} at k_13 {k[i_max]:.4g}; at the optimum "
          f"{PI[i_opt]:.4f}; IBO yield max {np.nanmax(Y):.4f} at k_13 "
          f"{k[int(np.nanargmax(Y))]:.4g}, at the optimum {Y[i_opt]:.4f}")
    if k_derate is not None:
        print(f'growth derated from k_13 {k_derate:.4g}')
    for name, y in (('PI', PI), ('IBO yield', Y)):
        tv, n_max = roughness(k, y)
        print(f'{name:10s} normalized total variation {tv:.2f}, '
              f'{n_max} local maxima')
    print(f'wrote {out_stem}.png / .pdf')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--csv', default=os.path.join(RESULTS_DIR,
                                                  SWEEP_STEM + '.csv'))
    ap.add_argument('--out-dir', default=os.path.join(
        RESULTS_DIR, 'publication', 'Objective-landscape'))
    ap.add_argument('--stem', default='k13_try_informed_sweep')
    ap.add_argument('--dpi', type=int, default=300)
    args = ap.parse_args(argv)
    frame, anchor = load(args.csv)
    plot(frame, anchor, os.path.join(args.out_dir, args.stem), args.dpi)


if __name__ == '__main__':
    main()

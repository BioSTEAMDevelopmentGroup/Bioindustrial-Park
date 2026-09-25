#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Profitability index and isobutanol yield along 1-D sweeps of any
metabolic_split_12d decision variable through the flagship profitability
campaign's optimum (trial #1912).

Two modes:

* --screen -- ranks the coarse slices of results/
  evaluate_axis_flagship_optimum_screen.csv (plus the in-band part of the
  full k_13 sweep, subsampled to the screen's density, as a reference) by how
  rough PI is against how smooth the isobutanol yield is, prints the table
  and draws a small-multiples overview (one twin-axis panel per variable).
  Roughness of a curve y along a slice, over its simulated points in axis
  order: normalized total variation TV/(max-min) (1 = monotone, 2 = one clean
  peak; every extra up-and-down adds 2), strict interior local maxima, and
  the largest single step as a fraction of the range. Score = excess TV of
  PI minus excess TV of the yield (excess = TV - 1).
* --param NAME -- the single-panel publication figure of a full-resolution
  sweep (results/evaluate_axis_flagship_optimum_<NAME>.csv), in the layout of
  plot_k13_flagship_sweep.py: PI on the left axis, isobutanol yield on the
  right, the trial marked, PI = 0, any part beyond the campaign band shaded,
  the onset of burden derating marked. --yield ethanol puts the ethanol yield
  on the right axis instead (stem suffix _etoh_yield). A third, outboard
  y axis carries the IRR (%) (--no-irr drops it), scaled so that the 15 %
  hurdle rate sits on the PI = 0 line (NPV at 15 % = 0 <=> IRR = 15 %).
  Each y axis is coloured like its curve; the legend lists only markers and
  reference lines.

Data: written by analyses/evaluate_axis_flagship_optimum.py (the simulation
stage, ask-first). Sim-safe: reads CSV / JSON only; reuses the typeface and
tick helpers of plot_k13_flagship_sweep.py by file path; never imports the
biorefineries package.

Usage:
    python plot_axis_flagship_sweep.py --screen [--dpi 300]
    python plot_axis_flagship_sweep.py --param threshold_conc [--dpi 300]
    python plot_axis_flagship_sweep.py --param k_3 --yield ethanol
"""
import os
import json
import argparse
import importlib.util

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import (FuncFormatter, LogLocator, MultipleLocator,
                               NullFormatter, NullLocator)

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    '_k13_sweep_plot', os.path.join(_here, 'plot_k13_flagship_sweep.py'))
k13p = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(k13p)

RESULTS_DIR = k13p.RESULTS_DIR
SWEEP_STEM = 'evaluate_axis_flagship_optimum'
K13_SWEEP_CSV = os.path.join(RESULTS_DIR, 'evaluate_k13_flagship_optimum.csv')
OUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Objective-landscape')
FONTS = k13p.FONTS
PI_COLOR, YIELD_COLOR = k13p.PI_COLOR, k13p.YIELD_COLOR
RATE = r'[$\mathrm{g·L}^{-1}·\mathrm{h}^{-1}$]'
CONC = r'[$\mathrm{g·L}^{-1}$]'

# decision variable -> (axis title, short title for the overview panels)
AXIS_LABELS = {
    'k_3': (rf'Pdc capacity, $k_{{3}}$ {RATE}', r'Pdc, $k_{3}$'),
    'k_6': (rf'Adh1 capacity, $k_{{6}}$ {RATE}', r'Adh1, $k_{6}$'),
    'k_13': (rf'ALS capacity, $k_{{13}}$ {RATE}', r'ALS, $k_{13}$'),
    'k_17': (rf'Adh6 capacity, $k_{{17}}$ {RATE}', r'Adh6, $k_{17}$'),
    'glycolysis': ('Glycolysis capacity multiplier', 'Glycolysis'),
    # the group's value IS the anchor k_14 (reference 1.0); k_15 / k_16
    # follow in the stoichiometric ratio 1.015 / 0.879
    'ehrlich_downstream': (rf'KARI–DHAD–Aro10 capacity, $k_{{14}}$ {RATE}',
                           'Ehrlich downstream'),
    'inhib_ethanol': ('Ethanol inhibition multiplier', 'Ethanol inhib.'),
    'inhib_isobutanol': ('Isobutanol inhibition multiplier',
                         'Isobutanol inhib.'),
    'inhib_acetate': ('Acetate inhibition multiplier', 'Acetate inhib.'),
    'threshold_conc': (f'Glucose spike threshold {CONC}', 'Spike threshold'),
    'target_delta': (f'Spike target above threshold {CONC}',
                     'Target − threshold'),
    'max_n_spikes': ('Maximum number of glucose spikes', 'Spike cap'),
}


#: steps below this fraction of a curve's range count as flat when local
#: maxima are counted (cross-process load-path drift is ~1e-4 relative)
FLAT_STEP = 1e-3


def roughness(y):
    """(normalized total variation, interior local maxima, largest single
    step / range) of y over its finite points, in order. A maximum is a
    rise followed by a fall once steps below FLAT_STEP of the range are
    treated as flat, so solver jitter on a plateau is not counted."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if y.size < 3:
        return np.nan, 0, np.nan
    span = y.max() - y.min()
    if span <= 0:
        return np.nan, 0, np.nan
    diff = np.diff(y)
    signs = np.sign(diff)[np.abs(diff) > FLAT_STEP*span]
    n_max = int(np.sum((signs[:-1] > 0) & (signs[1:] < 0)))
    return np.abs(diff).sum()/span, n_max, np.abs(diff).max()/span


def decoupled_step(PI, Y):
    """Largest single step of PI (as a fraction of its range) NOT matched by
    the isobutanol yield: max over neighbouring simulated points of
    |dPI|/range(PI) - |dY|/range(Y). Near 1 = a cliff in PI across which
    the yield does not move."""
    ok = np.isfinite(PI) & np.isfinite(Y)
    p, y = PI[ok], Y[ok]
    sp, sy = np.ptp(p), np.ptp(y)
    if p.size < 2 or sp <= 0 or sy <= 0:
        return np.nan
    return float(np.max(np.abs(np.diff(p))/sp - np.abs(np.diff(y))/sy))


def slice_summary(param, x, PI, Y, n_total):
    tv_p, nmax_p, jump_p = roughness(PI)
    tv_y, nmax_y, jump_y = roughness(Y)
    return dict(param=param, n=n_total, n_ok=int(np.isfinite(PI).sum()),
                PI_min=np.nanmin(PI), PI_max=np.nanmax(PI),
                PI_TV=tv_p, PI_maxima=nmax_p, PI_jump=jump_p,
                Y_min=np.nanmin(Y), Y_max=np.nanmax(Y),
                Y_TV=tv_y, Y_maxima=nmax_y, Y_jump=jump_y,
                decoupled=decoupled_step(PI, Y),
                score=(tv_p - 1) - (tv_y - 1))


#: the third (IRR) axis: colour, curve style, spine offset beyond the yield
#: spine, and the hurdle rate the PI is computed at (NPV at 15 % = 0 <=> IRR =
#: 15 %, so the IRR axis is scaled to put 15 % on the PI = 0 line)
IRR_COLOR = '#6A3D9A'
IRR_KW = dict(color=IRR_COLOR, lw=1.2, zorder=3)
#: the --param figure's curves are plain lines (only the optimum is marked)
PI_LINE_KW = {**k13p.PI_KW, 'marker': None}
YIELD_LINE_KW = {**k13p.YIELD_KW, 'marker': None}
IRR_AXIS_GAP_IN = 0.82
HURDLE_PCT = 15.0


def align_irr_axis(ax, ax_irr, IRR):
    """Scale the IRR axis (%) so HURDLE_PCT sits at PI = 0 on `ax` and every
    finite IRR fits; returns False (independent limits) when PI = 0 is not
    inside the PI axis or no IRR is finite."""
    lo, hi = ax.get_ylim()
    f0 = -lo/(hi - lo)   # axes fraction of PI = 0
    finite = IRR[np.isfinite(IRR)]
    if not finite.size:
        return False
    if not 0 < f0 < 1:
        pad = 0.08*max(np.ptp(finite), 1.0)
        ax_irr.set_ylim(finite.min() - pad, finite.max() + pad)
        return False
    up = max(finite.max() - HURDLE_PCT, 0.0)/(1 - f0)
    down = max(HURDLE_PCT - finite.min(), 0.0)/f0
    span = 1.06*max(up, down, 1e-9)   # % per unit axes fraction
    ax_irr.set_ylim(HURDLE_PCT - span*f0, HURDLE_PCT + span*(1 - f0))
    return True


#: --yield choice -> (sweep CSV column, axis / legend label, output-stem suffix)
YIELD_METRICS = {
    'isobutanol': ('IBO yield', 'Isobutanol yield', ''),
    'ethanol': ('EtOH yield', 'Ethanol yield', '_etoh_yield'),
}


def ok_series(frame, yield_column='IBO yield'):
    ok = frame['state'] == 'OK'
    return (frame['value'].to_numpy(float),
            np.where(ok, frame['PI'], np.nan),
            np.where(ok, frame[yield_column], np.nan))


def k13_reference(n_points):
    """The in-band part of the full k_13 sweep, subsampled to the points
    nearest a log grid of n_points (the screen's density)."""
    if not os.path.isfile(K13_SWEEP_CSV):
        return None
    d = pd.read_csv(K13_SWEEP_CSV).sort_values('k_13')
    d = d[d['k_13'] <= 4.0 + 1e-9].reset_index(drop=True)
    logk = np.log(d['k_13'].to_numpy(float))
    target = np.log(np.geomspace(d['k_13'].min(), 4.0, n_points))
    idx = sorted({int(np.argmin(np.abs(logk - t))) for t in target}
                 | set(np.flatnonzero(d['is_anchor'] == 1)))
    d = d.iloc[idx].rename(columns={'k_13': 'value'})
    d.insert(0, 'param', 'k_13')
    return d


# -----------------------------------------------------------------------------
# --screen: ranking table + small-multiples overview
# -----------------------------------------------------------------------------
def screen(dpi):
    path = os.path.join(RESULTS_DIR, f'{SWEEP_STEM}_screen.csv')
    with open(path[:-len('.csv')] + '_anchor.json') as fh:
        anchor = json.load(fh)
    frame = pd.read_csv(path)
    params = list(dict.fromkeys(frame['param']))
    n_screen = int(anchor['n_points'])
    ref = k13_reference(n_screen)
    slices = {p: frame[frame['param'] == p].sort_values('value')
              for p in params}
    if ref is not None:
        slices['k_13'] = ref
    rows = []
    for p, d in slices.items():
        x, PI, Y = ok_series(d)
        rows.append(slice_summary(p, x, PI, Y, len(d)))
    table = pd.DataFrame(rows).sort_values('score', ascending=False)
    with pd.option_context('display.width', 200, 'display.max_columns', 30,
                           'display.float_format', '{:.3f}'.format):
        print(f'anchor #{anchor["trial_number"]}: reproduced PI '
              f'{anchor["reproduced_PI"]:.5f} vs recorded '
              f'{anchor["recorded_PI"]:.5f}\n')
        print(table.to_string(index=False))
    states = frame['state'].value_counts().to_dict()
    print(f'\nscreen states: {states}')
    bad = frame[frame['state'] != 'OK']
    if len(bad):
        print(bad.groupby(['param', 'state']).size().to_string())
    table.to_csv(os.path.join(RESULTS_DIR, f'{SWEEP_STEM}_screen_ranking.csv'),
                 index=False)

    # --- overview figure, panels in ranking order ---
    k13p.apply_font_rcparams()
    order = list(table['param'])
    ncol = 4
    nrow = int(np.ceil(len(order)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(12.5, 2.75*nrow))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.07,
                        wspace=0.55, hspace=0.55)
    for ax, p in zip(axes.flat, order):
        d = slices[p]
        x, PI, Y = ok_series(d)
        ax_r = ax.twinx()
        ax.set_zorder(ax_r.get_zorder() + 1)
        ax.patch.set_visible(False)
        ax.axhline(0.0, color='black', lw=0.6)
        ax_r.plot(x, Y, color=YIELD_COLOR, lw=1.1, marker='s', ms=1.8)
        ax.plot(x, PI, color=PI_COLOR, lw=1.1, marker='o', ms=2.0)
        a = d['is_anchor'].to_numpy() == 1
        ax.plot(x[a], PI[a], **{**k13p.STAR_KW, 'markersize': 9})
        entry = anchor['bands'].get(p)
        if (entry and entry['log']) or p == 'k_13':
            ax.set_xscale('log')
            ratio = np.nanmax(x)/np.nanmin(x)
            if ratio < 3:        # the inhibition multipliers (0.75-1.5)
                ax.set_xticks([float(f'{v:.2g}') for v in
                               np.linspace(np.nanmin(x), np.nanmax(x), 4)])
                ax.xaxis.set_minor_locator(NullLocator())
            elif ratio < 30:     # under ~1.5 decades (glycolysis)
                ax.xaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
            if ratio < 30:
                ax.xaxis.set_minor_formatter(NullFormatter())
                ax.xaxis.set_major_formatter(
                    FuncFormatter(lambda v, _: f'{v:g}'))
        s = table.set_index('param').loc[p]
        ax.set_title(f"{AXIS_LABELS[p][1]}\nPI TV {s['PI_TV']:.1f}, "
                     f"yield TV {s['Y_TV']:.1f}", fontsize=9)
        ax.tick_params(labelsize=8)
        ax_r.tick_params(labelsize=8, colors=YIELD_COLOR)
        ax.tick_params(axis='y', colors=PI_COLOR)
    for ax in list(axes.flat)[len(order):]:
        ax.set_visible(False)
    fig.text(0.01, 0.5, 'Profitability index', color=PI_COLOR, rotation=90,
             va='center', fontsize=11)
    fig.text(0.99, 0.5, r'Isobutanol yield [$\mathrm{g·g}^{-1}$]',
             color=YIELD_COLOR, rotation=270, va='center', ha='right',
             fontsize=11)
    out = os.path.join(OUT_DIR, 'axis_flagship_screen')
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=dpi)
    plt.close(fig)
    print(f'wrote {out}.png / .pdf')


# -----------------------------------------------------------------------------
# --param: the publication figure of one full-resolution slice
# -----------------------------------------------------------------------------
def plot_param(param, dpi, stem=None, product='isobutanol', irr=True):
    yield_column, yield_label, stem_suffix = YIELD_METRICS[product]
    path = os.path.join(RESULTS_DIR, f'{SWEEP_STEM}_{param}.csv')
    with open(path[:-len('.csv')] + '_anchor.json') as fh:
        anchor = json.load(fh)
    frame = pd.read_csv(path).sort_values('value').reset_index(drop=True)
    band = anchor['bands'][param]
    x, PI, Y = ok_series(frame, yield_column)
    ok = frame['state'] == 'OK'
    IRR = frame['IRR'].to_numpy(float)*100
    IRR = np.where(ok & np.isfinite(IRR), IRR, np.nan)   # -inf: no IRR root
    dfac = frame['burden_factor'].to_numpy(float)
    i_opt = int(np.flatnonzero(frame['is_anchor'].to_numpy() == 1)[0])
    x_opt = x[i_opt]
    log = band['log']
    derated = np.flatnonzero(np.isfinite(dfac) & (dfac < 1 - 1e-9))

    # the legend lists only the markers / reference lines; the two curves are
    # identified by their colour-matched axes. Legend entries are known before
    # drawing, so the figure is sized to them (axes box fixed in inches).
    beyond = []
    if x.max() > band['high']*(1 + 1e-9):
        beyond.append('high')
    if x.min() < band['low']*(1 - 1e-9) - 1e-12:
        beyond.append('low')
    x_derate = (x[derated[0]] if derated.size and derated[0] > 0 else None)
    n_legend = (2 + bool(beyond) + (x_derate is not None)
                + bool((frame['state'] != 'OK').any()))
    axes_h, top_in, xlabel_in, entry_in = 2.59, 0.19, 0.72, 0.22
    axes_w, left_in, right_in = 3.234, 0.833, 0.833
    irr_gap_in = IRR_AXIS_GAP_IN if irr else 0.0   # yield spine -> IRR spine
    fig_w = left_in + axes_w + irr_gap_in + right_in
    fig_h = top_in + axes_h + xlabel_in + entry_in*n_legend + 0.1
    k13p.apply_font_rcparams()
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(left=left_in/fig_w, right=(left_in + axes_w)/fig_w,
                        top=1 - top_in/fig_h,
                        bottom=(fig_h - top_in - axes_h)/fig_h)
    ax_r = ax.twinx()
    ax.set_zorder(ax_r.get_zorder() + 2)   # PI (and its star) on top
    ax.patch.set_visible(False)
    if irr:   # a third y axis, outboard of the yield axis
        ax_irr = ax.twinx()
        ax_irr.set_zorder(ax_r.get_zorder() + 1)
        ax_irr.spines['right'].set_position(('axes', 1 + irr_gap_in/axes_w))
        for side in ('left', 'top', 'bottom'):
            ax_irr.spines[side].set_visible(False)
        ax_irr.spines['right'].set_color(IRR_COLOR)
        ax_irr.spines['right'].set_linewidth(1.2)

    if log:
        x_lo, x_hi = x.min()/1.25, x.max()*1.25
    else:
        pad = 0.03*(x.max() - x.min())
        x_lo, x_hi = x.min() - pad, x.max() + pad
    spans = {'high': (band['high'], x_hi), 'low': (x_lo, band['low'])}
    for side in beyond:   # bottom (right-axis) layer: never covers a curve
        ax_r.axvspan(*spans[side], **k13p.BAND_KW)
    if x_derate is not None:
        ax.axvline(x_derate, **k13p.DERATE_KW)
    ax.axhline(0.0, **k13p.BREAKEVEN_KW)

    ax_r.plot(x, Y, **YIELD_LINE_KW)
    ax_r.plot([x_opt], [Y[i_opt]], **k13p.RING_KW)
    if irr:   # no ring at the optimum: it would sit on the yield ring
        ax_irr.plot(x, IRR, **IRR_KW)
    ax.plot(x, PI, **PI_LINE_KW)
    ax.plot([x_opt], [PI[i_opt]], **k13p.STAR_KW)

    if log:
        ax.set_xscale('log')
    ax.set_xlim(x_lo, x_hi)
    ax.set_xlabel(AXIS_LABELS[param][0], fontsize=FONTS['axis_title'])
    ax.set_ylabel('Profitability index at 15% IRR', fontsize=FONTS['axis_title'],
                  color=PI_COLOR)
    ax_r.set_ylabel(yield_label + r' [$\mathrm{g·g}^{-1}$]',
                    fontsize=FONTS['axis_title'], rotation=270, labelpad=16,
                    color=YIELD_COLOR)
    y_top = np.nanmax(Y)
    ax_r.set_ylim(-0.03*y_top, 1.12*y_top)
    pi_lo, pi_hi = np.nanmin(PI), np.nanmax(PI)
    pad = 0.08*(pi_hi - pi_lo)
    ax.set_ylim(min(pi_lo - pad, -pad), pi_hi + 1.5*pad)
    aligned = False
    if irr:
        ax_irr.set_ylabel('IRR [%]',
                          fontsize=FONTS['axis_title'], rotation=270,
                          labelpad=16, color=IRR_COLOR)
        aligned = align_irr_axis(ax, ax_irr, IRR)
        if aligned:   # ticks through the hurdle rate: 0, 15, 30 %
            ax_irr.yaxis.set_major_locator(MultipleLocator(HURDLE_PCT))
            ax_irr.yaxis.set_minor_locator(MultipleLocator(HURDLE_PCT/3))
    ax.spines['left'].set_color(PI_COLOR)
    ax_r.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax_r.spines['right'].set_color(YIELD_COLOR)
    for side in ('left', 'right'):
        ax.spines[side].set_linewidth(1.2)
        ax_r.spines[side].set_linewidth(1.2)
    fig.canvas.draw()
    k13p.style_twin_ticks(ax, ax_r)
    ax.tick_params(axis='y', which='both', labelcolor=PI_COLOR)
    ax_r.tick_params(axis='y', which='both', labelcolor=YIELD_COLOR)
    if irr:
        for which, L in k13p.TICK_LEN.items():
            ax_irr.tick_params(axis='y', which=which, direction='inout',
                               length=2*L, left=False, right=True,
                               labelsize=FONTS['tick'], color=IRR_COLOR,
                               labelcolor=IRR_COLOR)
    for label in ax.get_yticklabels() + ax_r.get_yticklabels():
        label.set_text(k13p.minus(label.get_text()))

    def line(kw, label):
        return Line2D([], [], **{k: v for k, v in kw.items() if k != 'zorder'},
                      label=label)
    handles = [
        line(k13p.STAR_KW, f"Flagship optimum (#{anchor['trial_number']}, "
                           f"PI {PI[i_opt]:.2f})"),
        Line2D([], [], **k13p.BREAKEVEN_KW,
               label=('PI = 0 and IRR = 15 % (break-even at the hurdle rate)'
                      if aligned else
                      'PI = 0 (break-even at a 15 % hurdle rate)')),
    ]
    if beyond:
        handles.append(Patch(**k13p.BAND_KW, label='Beyond the campaign band'))
    if x_derate is not None:
        handles.append(Line2D([], [], **k13p.DERATE_KW,
                              label=f'Growth derated by enzyme burden '
                                    f'(from {x_derate:.3g})'))
    n_bad = int((~ok).sum())
    if n_bad:
        states = frame.loc[~ok, 'state'].value_counts().to_dict()
        handles.append(Line2D([], [], color='none', label='Not simulated: ' +
                              ', '.join(f'{v} {s.lower()}'
                                        for s, v in states.items())))
    fig.legend(handles=handles, loc='lower center', ncol=1, frameon=False,
               fontsize=FONTS['legend'], handlelength=2.2,
               bbox_to_anchor=((left_in + axes_w/2)/fig_w, 0.0))  # under the axes

    out = os.path.join(OUT_DIR, stem or f'{param}_flagship_sweep{stem_suffix}')
    os.makedirs(OUT_DIR, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}', dpi=dpi)
    plt.close(fig)

    s = slice_summary(param, x, PI, Y, len(frame))
    print(f"{len(frame)} points ({int(ok.sum())} simulated): {param} "
          f"{x.min():g}-{x.max():g}; anchor #{anchor['trial_number']} at "
          f"{x_opt:g} (reproduced PI {anchor['reproduced_PI']:.5f} vs "
          f"recorded {anchor['recorded_PI']:.5f})")
    print(f"PI {s['PI_min']:.3f}..{s['PI_max']:.3f}: TV {s['PI_TV']:.2f}, "
          f"{s['PI_maxima']} maxima, largest step {s['PI_jump']:.2f} of range")
    print(f"{yield_label} {s['Y_min']:.4f}..{s['Y_max']:.4f}: TV {s['Y_TV']:.2f}, "
          f"{s['Y_maxima']} maxima, largest step {s['Y_jump']:.2f} of range")
    if irr:
        print(f'IRR {np.nanmin(IRR):.2f}..{np.nanmax(IRR):.2f} %; '
              f'{int(np.isnan(IRR).sum())} points without an IRR; axis '
              f"{'aligned (15 % at PI = 0)' if aligned else 'NOT aligned'}")
    print(f'wrote {out}.png / .pdf')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--screen', action='store_true')
    mode.add_argument('--param')
    ap.add_argument('--yield', dest='product', choices=tuple(YIELD_METRICS),
                    default='isobutanol',
                    help='yield on the right axis of a --param figure')
    ap.add_argument('--no-irr', dest='irr', action='store_false',
                    help='drop the third (IRR) axis of a --param figure')
    ap.add_argument('--stem', default=None)
    ap.add_argument('--dpi', type=int, default=300)
    args = ap.parse_args(argv)
    if args.screen:
        screen(args.dpi)
    else:
        plot_param(args.param, args.dpi, args.stem, args.product, args.irr)


if __name__ == '__main__':
    main()

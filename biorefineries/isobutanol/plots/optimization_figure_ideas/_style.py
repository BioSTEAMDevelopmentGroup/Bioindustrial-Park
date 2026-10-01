#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Shared style layer of the kinetic-BO campaign figures: palette, fonts /
rcParams, tick styling, inch-based placement, IRR axes with the loss band,
direct-label keys, and the render checks (text overlap, font-size floor,
Arial glyph coverage, colour-vision palette check) plus the save helper.

Matplotlib only (Agg backend, set at import); sim-safe. Import it from a
sibling figure script with::

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import _common as C
    import _style as S
    fig = S.new_figure(*S.MAIN_SIZE)       # rcParams applied
    ...
    S.check_figure(fig)                    # draw + overlap / font / glyph
    S.save(fig, 'kinBO_main_scouts')       # PNG 300 dpi + PDF (+ _latest)

Public API
----------
Constants
    PALETTE {key: hex} (spec section 2), TEXT / NOTE colours, FS (font
    sizes, pt on the 10-in canvas), MIN_FONT_PT = 9, CANVAS_WIDTH_IN = 10,
    PRINT_WIDTH_MM = 180, PRINT_SCALE (0.709), MAIN_SIZE (10 x 7.9),
    S1_SIZE (10 x 3.9), S2_SIZE (10 x 5.6), DPI = 300
    IRR_LIM (-6, 31), LOSS_BAND (-6, 0), LOSS_CENTER -3, LOSS_JITTER 2.4
    PLATEAU_LINE / START_LINE / ZERO_LINE (line kwargs)
    UNIT_TITER 'g·L$^{-1}$', UNIT_PROD, UNIT_PROTEOME 'g·(g DCW)$^{-1}$'
    CVD_PAIRS, CVD_MIN_DE = 12, GRAY_PAIRS {pair: min dL*}
Set-up and placement
    apply_style()                  house rcParams (Arial, mathtext Arial,
                                   stixsans fallback, pdf/ps fonttype 42)
    new_figure(w, h)               apply_style() + plt.figure(figsize=(w, h))
    inch_axes(fig, x, y, w, h, **kw)  axes placed in inches (origin bottom-
                                   left)
    fig_text(fig, x, y, s, **kw)   figure text at (x, y) inches
    panel_letter(fig, x, y, letter) / panel_title(fig, x, y, text, dx=0.28)
    bold_axis_title(label)         bold name, regular [units] / (units)
    style_ticks(ax, x=True, y=True, minor_x=True, minor_y=True)
                                   ticks on all four sides; left/bottom
                                   in+out, top/right in; major 4 / minor 2;
                                   log axes keep log minors, logit axes get
                                   none; x=False / y=False removes that
                                   axis's ticks (categorical rows)
IRR axes
    irr_plot(v, rng=None)          IRR fraction -> plot coordinate (% if
                                   >= 0; a loss -> -3, jittered +-2.4 if rng)
    irr_axis(ax, which='y', step=5, title='IRR [%]', band=True,
             loss_fs=None, zero_label=True)
                                   limits -6..31, ticks 'loss', 0, step..30,
                                   minors outside the band, loss band + 0
                                   line (narrow horizontal axis: loss_fs=9)
    log_trial_axis(ax, which='x', lim=(1, 2200))  log trial axis, majors
                                   '1', '10', '100', '1,000', log minors
    logit_pct_axis(ax, which='x', ticks_pct=(1, 2, 5, ..., 90),
                   lim_pct=(0.8, 93))  logit share axis (plot FRACTIONS)
    loss_band(ax, which='y'), plateau_line(ax, U_pct, which='y', **kw),
    halo(text, lw=2) white stroke behind a label over dots (no bbox)
    start_line(ax, start_pct, which='y', **kw), tint_above(ax, U_pct,
    which='y'), panel_rng() (np.random.default_rng(0), one per panel)
Campaign encodings
    campaign_style(c)              {'color', 'dark', 'light', 'text',
                                   'marker'} for a _common.Campaign (or
                                   'base' for the starting strain)
Keys
    inline_key(fig, x, y, items, fontsize=9, ...)  one line of marker /
                                   swatch + text items in figure inches
                                   (glyphs drawn as markers, not Unicode)
Render checks (call after building the figure)
    text_overlaps(fig, exempt=(), tol_px=0.5) -> [str]  (empty = pass)
    tick_label_collisions(fig, min_gap_in=0.03) -> [str]  same-axis tick
                                   labels that overlap, or x-axis labels
                                   closer than 0.03 in (e.g. 'loss0')
    min_font_check(fig, min_pt=9) -> [str]
    glyph_check(fig) -> [str]      characters Arial cannot render
    check_figure(fig, size=None, raise_on_fail=True, exempt=()) -> dict
    FigureCheckError (AssertionError subclass)
Palette checks
    cvd_check(pairs=CVD_PAIRS, min_de=12, raise_on_fail=True,
              verbose=False, dataviz=True) -> (ok, report)
    delta_e76(h1, h2, kind=None), lstar(hex), darken(hex, dL=5)
Saving
    save(fig, stem, out_dir=OUT_DIR, stamp=None, latest=True,
         max_pdf_mb=10) -> {'png', 'pdf', 'png_latest', 'pdf_latest',
         'png_px'}: <stem>_<YYYY.MM.DD-HH.MM>.{png,pdf} at 300 dpi plus a
         stable copy <stem>_latest.{png,pdf}; never bbox_inches='tight'
"""
import glob
import importlib.util
import os
import re
import shutil
import struct
import sys
import tempfile
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt                       # noqa: E402
from matplotlib.lines import Line2D                         # noqa: E402
from matplotlib.patches import Rectangle                    # noqa: E402
from matplotlib.text import Text                            # noqa: E402
from matplotlib.ticker import (AutoMinorLocator, FixedFormatter,  # noqa: E402
                               FixedLocator, NullLocator)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import MAX_PATH_CHARS, OUT_DIR                 # noqa: E402

__all__ = [
    'PALETTE', 'TEXT', 'NOTE', 'FS', 'MIN_FONT_PT', 'CANVAS_WIDTH_IN',
    'PRINT_WIDTH_MM', 'PRINT_SCALE', 'MAIN_SIZE', 'S1_SIZE', 'S2_SIZE', 'DPI',
    'IRR_LIM', 'LOSS_BAND', 'LOSS_CENTER', 'LOSS_JITTER', 'PLATEAU_LINE',
    'START_LINE', 'ZERO_LINE', 'UNIT_TITER', 'UNIT_PROD', 'UNIT_PROTEOME',
    'apply_style', 'new_figure', 'inch_axes', 'fig_text', 'panel_letter',
    'panel_title', 'bold_axis_title', 'style_ticks', 'irr_plot', 'irr_axis',
    'loss_band', 'plateau_line', 'start_line', 'tint_above', 'panel_rng',
    'log_trial_axis', 'logit_pct_axis', 'tick_label_collisions', 'TINT_ALPHA',
    'campaign_style', 'inline_key', 'text_overlaps', 'min_font_check',
    'glyph_check', 'check_figure', 'FigureCheckError', 'CVD_PAIRS',
    'CVD_MIN_DE', 'GRAY_PAIRS', 'cvd_check', 'delta_e76', 'lstar', 'darken',
    'save',
]

# %% Palette (spec section 2) -----------------------------------------------------
# House hexes for the profitability campaigns and the baseline
# (plot_kin_opt_parameter_sets.py HUE_COLORS[0] / RELAY_COLOR /
# BASELINE_COLOR). The six per-campaign scout hues are replaced by two
# family hues (ethanol amber, isobutanol violet); the role (yield / titer /
# productivity) is the marker shape. All seven CVD_PAIRS pass cvd_check()
# unchanged (min CIE76 dE 17.7, adh1/adh6 under tritanopia), so no hex was
# darkened.
PALETTE = {
    'unin': '#18C4DC',          # uninformed lines / markers; dots alpha 0.65
    'unin_text': '#0E8FA1',     # cyan text (legible at 9-10 pt)
    'relay': '#0B6E7A',         # TRY-informed lines / markers / text
    'etoh': '#E0A030',          # ethanol product (d), Pdc (e)
    'etoh_dark': '#A86F0C',     # ethanol-scout markers (edge), headers
    'etoh_light': '#F0CF8E',    # ethanol-scout dots / strips (alpha 0.55)
    'ibo': '#7B5BA6',           # isobutanol product (d), ALS->Aro10 (e)
    'ibo_dark': '#5E4385',      # isobutanol-scout markers, headers
    'ibo_light': '#BFAEDA',     # isobutanol-scout dots / strips
    'adh1': '#F3D9A6',          # light shade of Pdc (e)
    'adh6': '#D3C6E8',          # light shade of ALS->Aro10 (e)
    'gly': '#8C8C8C',           # glycolysis (e)
    'tca': '#CFCFCF',           # TCA + acetate (e)
    'base': '#90918e',          # starting strain line / open diamond
    'seed': '#8F8F8F',          # seeds diamond (f)
    'pw': '#7A7A7A',            # price-weighted-yield replicate scout (S1)
    'pw_light': '#C4C4C4',      # its dots
    'loss_band': '#EDEDED',     # loss band fill
    'startup_band': '#F4F4F4',  # uninformed space-filling start-up (1-51)
    'tint': '#18C4DC',          # "above the plateau" region, alpha TINT_ALPHA
    'text': '#222222',
    'note': '#666666',
}
TINT_ALPHA = 0.07
TEXT = PALETTE['text']
NOTE = PALETTE['note']

# %% Canvas, fonts, units ---------------------------------------------------------
CANVAS_WIDTH_IN = 10.0
PRINT_WIDTH_MM = 180.0
PRINT_SCALE = PRINT_WIDTH_MM / 25.4 / CANVAS_WIDTH_IN       # 0.709
MAIN_SIZE = (10.0, 7.9)
S1_SIZE = (10.0, 3.9)
S2_SIZE = (10.0, 5.6)
DPI = 300
# canvas point sizes (print = x 0.709): tick 12 (8.5), axis title 12,
# panel letter 14 bold, panel title 12 bold, row label 11, annotation 10
# (7.1), note / key / table header 9 (6.4, the floor)
FS = {'tick': 12, 'axis': 12, 'letter': 14, 'title': 12, 'row': 11,
      'table': 11, 'annot': 10, 'note': 9, 'key': 9}
MIN_FONT_PT = 9
TICK_LABEL_GAP_IN = 0.03          # min gap between same-axis tick labels
HALO_LW = 2.0                     # white stroke behind text over dots (pt)
UNIT_TITER = 'g·L$^{-1}$'
UNIT_PROD = 'g·L$^{-1}$·h$^{-1}$'
UNIT_PROTEOME = 'g·(g DCW)$^{-1}$'

# %% IRR axes -----------------------------------------------------------------------
IRR_LIM = (-6.0, 31.0)
LOSS_BAND = (-6.0, 0.0)
LOSS_CENTER = -3.0
LOSS_JITTER = 2.4
PLATEAU_LINE = dict(color=PALETTE['unin'], lw=1.1, ls='--', zorder=1.2)
START_LINE = dict(color=PALETTE['base'], lw=1.1, ls=(0, (1.5, 1.2)),
                  zorder=1.1)
ZERO_LINE = dict(color='0.6', lw=0.6, zorder=0.6)


def apply_style():
    """House rcParams (pk.apply_fonts + the formatting preferences)."""
    rc = plt.rcParams
    rc['font.family'] = 'sans-serif'
    rc['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    rc['font.size'] = FS['tick']
    rc['xtick.labelsize'] = FS['tick']
    rc['ytick.labelsize'] = FS['tick']
    rc['axes.labelsize'] = FS['axis']
    rc['axes.titlesize'] = FS['title']
    rc['legend.fontsize'] = FS['note']
    rc['axes.linewidth'] = 0.8
    rc['hatch.linewidth'] = 0.6
    rc['mathtext.fontset'] = 'custom'
    rc['mathtext.rm'] = 'Arial'
    rc['mathtext.it'] = 'Arial:italic'
    rc['mathtext.bf'] = 'Arial:bold'
    rc['mathtext.fallback'] = 'stixsans'
    rc['pdf.fonttype'] = 42
    rc['ps.fonttype'] = 42
    rc['savefig.dpi'] = DPI
    rc['savefig.bbox'] = 'standard'       # never tight: the layout is in inches
    rc['text.color'] = TEXT
    rc['axes.labelcolor'] = TEXT
    rc['xtick.color'] = TEXT
    rc['ytick.color'] = TEXT
    rc['axes.edgecolor'] = TEXT
    rc['axes.unicode_minus'] = True       # Arial has U+2212


def new_figure(w, h):
    apply_style()
    return plt.figure(figsize=(w, h))


def inch_axes(fig, x, y, w, h, **kw):
    """Axes at (x, y) inches from the canvas's bottom-left, w x h inches."""
    W, H = fig.get_size_inches()
    return fig.add_axes((x / W, y / H, w / W, h / H), **kw)


def fig_text(fig, x, y, s, **kw):
    """Figure text at (x, y) inches."""
    kw.setdefault('color', TEXT)
    return fig.text(x, y, s, transform=fig.dpi_scale_trans, **kw)


def panel_letter(fig, x, y, letter, **kw):
    kw.setdefault('fontsize', FS['letter'])
    kw.setdefault('fontweight', 'bold')
    kw.setdefault('va', 'baseline')
    kw.setdefault('ha', 'left')
    return fig_text(fig, x, y, letter, **kw)


def panel_title(fig, x, y, text, dx=0.28, **kw):
    """Panel title starting `dx` inches right of the letter at (x, y)."""
    kw.setdefault('fontsize', FS['title'])
    kw.setdefault('fontweight', 'bold')
    kw.setdefault('va', 'baseline')
    kw.setdefault('ha', 'left')
    return fig_text(fig, x + dx, y, text, **kw)


# --- bold axis titles, unit annotations regular (pk._bold_axis_title rule) ----
_MATHTEXT_ESCAPE = {'%': r'\%', '&': r'\&', '#': r'\#', '_': r'\_',
                    '{': r'\{', '}': r'\}'}


def _bold_run(run):
    if not run.strip():
        return run
    parts = []
    for seg in re.split(r'(\$[^$]*\$)', run):
        if not seg:
            continue
        if seg.startswith('$') and seg.endswith('$'):
            inner = seg[1:-1]
            if inner:
                parts.append(r'\mathbf{' + inner + '}')
        else:
            esc = ''.join(_MATHTEXT_ESCAPE.get(c, c) for c in seg)
            parts.append(r'\mathbf{' + esc.replace(' ', r'\ ') + '}')
    return '$' + ''.join(parts) + '$'


def _split_units(line):
    segs, buf, close = [], '', None
    for c in line:
        if close is None and c in '[(':
            if buf:
                segs.append((False, buf))
            buf, close = c, (']' if c == '[' else ')')
        elif close is not None:
            buf += c
            if c == close:
                segs.append((True, buf))
                buf, close = '', None
        else:
            buf += c
    if buf:
        segs.append((close is not None, buf))
    return segs


def bold_axis_title(label):
    """'IRR [%]' -> bold 'IRR' + regular ' [%]' (copy of the house
    plot_kin_opt_parameter_sets._bold_axis_title rule; multi-line safe)."""
    out = []
    for line in label.split('\n'):
        out.append(''.join(seg if is_unit else _bold_run(seg)
                           for is_unit, seg in _split_units(line)))
    return '\n'.join(out)


# %% Ticks ---------------------------------------------------------------------------
def style_ticks(ax, x=True, y=True, minor_x=True, minor_y=True):
    """House ticks: all four sides, major + minor; left/bottom in-and-out,
    top/right inward; major 4 pt, minor 2 pt (plot_relay_campaign_progress
    .style_ticks). Call AFTER limits and locators: a linear axis whose minor
    locator is still the default NullLocator gets an AutoMinorLocator (a
    FixedLocator set by irr_axis is kept), a log axis keeps its log minors,
    a logit axis gets none. x=False / y=False removes that axis's ticks and
    tick labels (categorical rows)."""
    for axis, on, minor in ((ax.xaxis, x, minor_x), (ax.yaxis, y, minor_y)):
        scale = axis.get_scale()
        if not (on and minor) or scale == 'logit':
            axis.set_minor_locator(NullLocator())
        elif scale == 'linear' and isinstance(axis.get_minor_locator(),
                                              NullLocator):
            axis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', top=True, right=True, direction='in',
                   length=4)
    ax.tick_params(which='minor', top=True, right=True, direction='in',
                   length=2)
    ax.tick_params(axis='x', which='major', bottom=True, direction='inout',
                   length=4)
    ax.tick_params(axis='x', which='minor', bottom=True, direction='inout',
                   length=2)
    ax.tick_params(axis='y', which='major', left=True, direction='inout',
                   length=4)
    ax.tick_params(axis='y', which='minor', left=True, direction='inout',
                   length=2)
    # 'direction' is per axis, not per side: re-point top/right inward
    for tick in ax.xaxis.get_major_ticks() + ax.xaxis.get_minor_ticks():
        tick.tick2line.set_marker(3)        # TICKDOWN on the top spine
    for tick in ax.yaxis.get_major_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick2line.set_marker(0)        # TICKLEFT on the right spine
    if not x:
        ax.tick_params(axis='x', which='both', bottom=False, top=False,
                       labelbottom=False, labeltop=False)
    if not y:
        ax.tick_params(axis='y', which='both', left=False, right=False,
                       labelleft=False, labelright=False)


# %% IRR helpers ---------------------------------------------------------------------
def irr_plot(v, rng=None):
    """IRR fraction -> plot coordinate: 100*v if v >= 0; a loss (v < 0 or
    -inf) -> LOSS_CENTER, + U(-2.4, 2.4) jitter when `rng` is given (one
    draw per loss, in input order). NaN (not solved) stays NaN.
    Scalar in -> float out; array in -> array out."""
    a = np.atleast_1d(np.asarray(v, dtype=float))
    loss = (a < 0) | np.isneginf(a)
    out = np.where(loss, LOSS_CENTER, 100.0 * a)
    if rng is not None and loss.any():
        out[loss] = LOSS_CENTER + rng.uniform(-LOSS_JITTER, LOSS_JITTER,
                                              int(loss.sum()))
    return float(out[0]) if np.ndim(v) == 0 else out


def loss_band(ax, which='y'):
    """Grey loss band (-6..0) and the 0 line."""
    if which == 'y':
        ax.axhspan(*LOSS_BAND, color=PALETTE['loss_band'], lw=0, zorder=0.1)
        ax.axhline(0.0, **ZERO_LINE)
    else:
        ax.axvspan(*LOSS_BAND, color=PALETTE['loss_band'], lw=0, zorder=0.1)
        ax.axvline(0.0, **ZERO_LINE)


def irr_axis(ax, which='y', step=5, title='IRR [%]', band=True,
             loss_fs=None, zero_label=True):
    """An IRR axis: limits -6..31 %, major ticks 'loss' (at -3), 0, step, ..
    30, minors (step/5) outside the loss band, the loss band and 0 line, and
    the bold title (None = no title). Call style_ticks(ax) afterwards.

    On a HORIZONTAL IRR axis narrower than ~2.8 in the 'loss' and '0' tick
    labels collide at 12 pt (check_figure reports it as a tick collision):
    pass loss_fs=9 (fits down to ~2.0 in) or zero_label=False."""
    axis = ax.yaxis if which == 'y' else ax.xaxis
    majors = list(range(0, 31, step))
    labels = ['loss'] + [str(m) for m in majors]
    if not zero_label:
        labels[1] = ''
    axis.set_major_locator(FixedLocator([LOSS_CENTER] + majors))
    axis.set_major_formatter(FixedFormatter(labels))
    mstep = step / 5.0
    minors = [m for m in np.arange(0, IRR_LIM[1] + 1e-9, mstep)
              if min(abs(m - M) for M in majors) > 1e-9]
    axis.set_minor_locator(FixedLocator(minors))
    (ax.set_ylim if which == 'y' else ax.set_xlim)(*IRR_LIM)
    if loss_fs is not None:
        # every tick exists now (fixed locator), so this sticks
        t0 = axis.get_major_ticks()[0]
        t0.label1.set_fontsize(loss_fs)
        t0.label2.set_fontsize(loss_fs)
    if band:
        loss_band(ax, which)
    if title:
        (ax.set_ylabel if which == 'y' else ax.set_xlabel)(
            bold_axis_title(title))


def log_trial_axis(ax, which='x', lim=(1, 2200), title='Simulated trials'):
    """Log trial axis (panel a / S1a): limits 1..2,200, majors 1, 10, 100,
    1,000 labelled '1', '10', '100', '1,000', log minors (2..9). Call
    style_ticks(ax) afterwards (it keeps the log minors)."""
    from matplotlib.ticker import LogLocator
    axis = ax.xaxis if which == 'x' else ax.yaxis
    (ax.set_xscale if which == 'x' else ax.set_yscale)('log')
    (ax.set_xlim if which == 'x' else ax.set_ylim)(*lim)
    majors = [m for m in (1, 10, 100, 1000, 10000)
              if lim[0] <= m <= lim[1]]
    axis.set_major_locator(FixedLocator(majors))
    axis.set_major_formatter(FixedFormatter([f'{m:,}' for m in majors]))
    axis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    axis.set_minor_formatter(FixedFormatter([]))
    if title:
        (ax.set_xlabel if which == 'x' else ax.set_ylabel)(
            bold_axis_title(title))


def logit_pct_axis(ax, which='x', ticks_pct=(1, 2, 5, 10, 20, 50, 80, 90),
                   lim_pct=(0.8, 93), title='Trials making\nisobutanol [%]'):
    """Logit axis for a share in % (panel f's exploration): the data must be
    plotted as FRACTIONS (pct / 100); majors at `ticks_pct` with plain %
    labels, no minors (style_ticks enforces none on a logit axis)."""
    axis = ax.xaxis if which == 'x' else ax.yaxis
    (ax.set_xscale if which == 'x' else ax.set_yscale)('logit')
    (ax.set_xlim if which == 'x' else ax.set_ylim)(lim_pct[0] / 100,
                                                   lim_pct[1] / 100)
    axis.set_major_locator(FixedLocator([t / 100 for t in ticks_pct]))
    axis.set_major_formatter(FixedFormatter([f'{t:g}' for t in ticks_pct]))
    axis.set_minor_locator(NullLocator())
    if title:
        (ax.set_xlabel if which == 'x' else ax.set_ylabel)(
            bold_axis_title(title))


def plateau_line(ax, U_pct, which='y', **kw):
    """Dashed cyan line at the uninformed plateau U (in %)."""
    k = {**PLATEAU_LINE, **kw}
    return (ax.axhline if which == 'y' else ax.axvline)(U_pct, **k)


def start_line(ax, start_pct, which='y', **kw):
    """Dotted grey starting-strain line (in %)."""
    k = {**START_LINE, **kw}
    return (ax.axhline if which == 'y' else ax.axvline)(start_pct, **k)


def tint_above(ax, U_pct, which='y'):
    """Faint cyan tint over IRR > U (the "above the plateau" region)."""
    span = ax.axhspan if which == 'y' else ax.axvspan
    return span(U_pct, IRR_LIM[1], color=PALETTE['tint'], alpha=TINT_ALPHA,
                lw=0, zorder=0.05)


def halo(t, lw=HALO_LW, color='white'):
    """A thin white stroke behind a Text (instead of an opaque bbox): keeps a
    label legible over dots while hiding only the dots under its glyphs.
    Returns the Text."""
    from matplotlib import patheffects
    t.set_path_effects([patheffects.withStroke(linewidth=lw,
                                               foreground=color)])
    return t


def panel_rng():
    """The per-panel jitter generator (spec: default_rng(0), one per panel,
    fixed draw order)."""
    return np.random.default_rng(0)


# %% Campaign encodings ----------------------------------------------------------------
def campaign_style(c):
    """Encodings of a _common.Campaign (or the string 'base'): 'color' (lines
    / filled markers), 'dark' (marker edges, headers), 'light' (trial dots),
    'text' (direct-label colour), 'marker'."""
    if c == 'base':
        return {'color': PALETTE['base'], 'dark': PALETTE['base'],
                'light': PALETTE['base'], 'text': NOTE, 'marker': 'D'}
    if c.family == 'profit':
        col = PALETTE['relay'] if c.is_relay else PALETTE['unin']
        txt = PALETTE['relay'] if c.is_relay else PALETTE['unin_text']
        return {'color': col, 'dark': col, 'light': col, 'text': txt,
                'marker': c.marker}
    fam = {'etoh': 'etoh', 'ibo': 'ibo', 'pw': 'pw'}[c.family]
    dark = PALETTE[fam + '_dark'] if fam != 'pw' else PALETTE['pw']
    return {'color': PALETTE[fam], 'dark': dark,
            'light': PALETTE[fam + '_light'], 'text': dark,
            'marker': c.marker}


# %% Keys ----------------------------------------------------------------------------------
def inline_key(fig, x, y, items, fontsize=None, gap_in=0.05, item_gap_in=0.16,
               ha='left'):
    """Draw one line of key items at (x, y) inches (y = the text's vertical
    centre). Each item is a dict with 'text' (and optional 'color', 'style',
    'weight') plus either a Line2D marker spec ('marker', 'ms' [pt], 'mfc',
    'mec', 'mew') or a square swatch ('swatch': colour, 'size_in': 0.11).
    Glyphs are drawn as markers so they render in every font. ha='left'
    starts at x; ha='right' ends at x. Returns (x_end_in, artists)."""
    fontsize = FS['key'] if fontsize is None else fontsize
    renderer = fig.canvas.get_renderer()
    tr = fig.dpi_scale_trans
    widths, specs = [], []
    for it in items:                       # measure first (for ha='right')
        t = Text(0, 0, it.get('text', ''), fontsize=fontsize,
                 fontweight=it.get('weight', 'normal'),
                 fontstyle=it.get('style', 'normal'))
        t.set_figure(fig)
        tw = t.get_window_extent(renderer).width / fig.dpi if it.get(
            'text') else 0.0
        if 'swatch' in it:
            gw = it.get('size_in', 0.11)
        elif 'marker' in it:
            gw = it.get('ms', 6.0) / 72.0
        else:
            gw = 0.0
        widths.append((gw, tw))
        specs.append(it)
    total = sum(gw + (gap_in if gw and tw else 0) + tw for gw, tw in widths) \
        + item_gap_in * (len(items) - 1)
    cx = x if ha == 'left' else x - total
    artists = []
    for it, (gw, tw) in zip(specs, widths):
        if 'swatch' in it:
            s = it.get('size_in', 0.11)
            r = Rectangle((cx, y - s / 2), s, s, transform=tr,
                          facecolor=it['swatch'], edgecolor=it.get(
                              'edgecolor', 'none'), lw=it.get('lw', 0.0))
            fig.add_artist(r)
            artists.append(r)
        elif 'marker' in it:
            ln = Line2D([cx + gw / 2], [y], transform=tr, ls='none',
                        marker=it['marker'], ms=it.get('ms', 6.0),
                        mfc=it.get('mfc', TEXT), mec=it.get('mec', TEXT),
                        mew=it.get('mew', 1.0))
            fig.add_artist(ln)
            artists.append(ln)
        tx = cx + gw + (gap_in if gw and tw else 0)
        if it.get('text'):
            artists.append(fig.text(
                tx, y, it['text'], transform=tr, fontsize=fontsize,
                va='center', ha='left', color=it.get('color', TEXT),
                fontweight=it.get('weight', 'normal'),
                fontstyle=it.get('style', 'normal')))
        cx = tx + tw + item_gap_in
    return cx - item_gap_in, artists


# %% Render checks -------------------------------------------------------------------
class FigureCheckError(AssertionError):
    pass


def _safe_print(msg):
    """print() that never fails on a console codepage (cp1252) lacking a
    character of the message (e.g. U+2605): unencodable characters are
    backslash-escaped."""
    enc = getattr(sys.stdout, 'encoding', None) or 'utf-8'
    print(str(msg).encode(enc, errors='backslashreplace').decode(enc))


def _drawn_texts(fig):
    """[(Text, bbox_px, tick_axis_id or None, where)] of every visible,
    non-empty Text the figure draws (tick labels only for ticks inside the
    view interval; axis labels only when that axis is drawn)."""
    renderer = fig.canvas.get_renderer()
    owner, skip, where = {}, set(), {}
    for i, ax in enumerate(fig.axes):
        for name, axis in (('x', ax.xaxis), ('y', ax.yaxis)):
            shown = ax.axison and axis.get_visible() and ax.get_visible()
            try:
                drawn = {id(t) for t in axis._update_ticks()}
            except Exception:                    # private API fallback
                drawn = {id(t) for t in axis.get_major_ticks()
                         + axis.get_minor_ticks()}
            for t in axis.get_major_ticks() + axis.get_minor_ticks():
                for lab in (t.label1, t.label2):
                    owner[id(lab)] = id(axis)
                    where[id(lab)] = f'ax{i}.{name}tick'
                    if not shown or id(t) not in drawn:
                        skip.add(id(lab))
            where[id(axis.label)] = f'ax{i}.{name}label'
            if not shown:
                skip.add(id(axis.label))
                skip.add(id(axis.offsetText))
        for t in ax.texts:
            where[id(t)] = f'ax{i}.text'
    out = []
    for t in fig.findobj(Text):
        if id(t) in skip or not t.get_visible():
            continue
        if not t.get_text().strip() or t.get_alpha() == 0:
            continue
        ax = t.axes
        if ax is not None and not ax.get_visible():
            continue
        bb = t.get_window_extent(renderer)
        if bb.width <= 0 or bb.height <= 0:
            continue
        out.append((t, bb, owner.get(id(t)), where.get(id(t), 'fig.text')))
    return out


def text_overlaps(fig, exempt=(), tol_px=0.5):
    """Pairwise overlaps of the window extents of all drawn non-empty Text
    artists (tick labels of the SAME axis are exempt from each other), plus
    any text outside the canvas. `exempt`: Text artists to ignore. Draws the
    canvas first. Returns a list of messages; [] = pass."""
    fig.canvas.draw()
    ex = {id(t) for t in exempt}
    items = [x for x in _drawn_texts(fig) if id(x[0]) not in ex]
    msgs = []
    W, H = fig.bbox.width, fig.bbox.height
    for t, bb, _, w in items:
        if (bb.x0 < -tol_px or bb.y0 < -tol_px or bb.x1 > W + tol_px
                or bb.y1 > H + tol_px):
            msgs.append(f'outside canvas: {t.get_text()!r} [{w}] '
                        f'({bb.x0:.0f},{bb.y0:.0f})-({bb.x1:.0f},{bb.y1:.0f})'
                        f' px of {W:.0f}x{H:.0f}')
    for i in range(len(items)):
        ti, bi, oi, wi = items[i]
        for j in range(i + 1, len(items)):
            tj, bj, oj, wj = items[j]
            if oi is not None and oi == oj:
                continue
            dx = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            dy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if dx > tol_px and dy > tol_px:
                msgs.append(f'overlap: {ti.get_text()!r} [{wi}] x '
                            f'{tj.get_text()!r} [{wj}] ({dx:.1f} x {dy:.1f}'
                            f' px)')
    return msgs


def min_font_check(fig, min_pt=MIN_FONT_PT):
    """Every drawn Text artist must be >= min_pt (9 pt on the 10-in canvas
    = 6.4 pt in print). Returns a list of messages; [] = pass."""
    fig.canvas.draw()
    return [f'{t.get_fontsize():.1f} pt < {min_pt} pt: {t.get_text()!r} '
            f'[{w}]' for t, _, _, w in _drawn_texts(fig)
            if t.get_fontsize() < min_pt - 1e-9]


_ARIAL_CMAP = None


def _arial_cmap():
    global _ARIAL_CMAP
    if _ARIAL_CMAP is None:
        try:
            from fontTools.ttLib import TTFont
            from matplotlib import font_manager
            path = font_manager.findfont('Arial', fallback_to_default=False)
            _ARIAL_CMAP = set(TTFont(path, lazy=True).getBestCmap())
        except Exception:                         # fontTools / Arial missing
            _ARIAL_CMAP = False
    return _ARIAL_CMAP


def glyph_check(fig):
    """Characters of drawn text (math segments included) that Arial cannot
    render -- they would fall back to DejaVu / STIX (e.g. the Unicode
    triangle / star / diamond: draw those as markers, see inline_key).
    Returns a list of messages; [] = pass (or fontTools unavailable)."""
    cmap = _arial_cmap()
    if not cmap:
        return []
    fig.canvas.draw()
    msgs = []
    for t, _, _, w in _drawn_texts(fig):
        s = t.get_text()
        plain = re.sub(r'\\[A-Za-z]+', '', s)      # mathtext commands
        bad = sorted({c for c in plain if ord(c) > 127 and ord(c)
                      not in cmap})
        if bad:
            msgs.append(f'not in Arial: {[f"U+{ord(c):04X} {c}" for c in bad]}'
                        f' in {s!r} [{w}]')
    return msgs


def tick_label_collisions(fig, tol_px=0.5, min_gap_in=TICK_LABEL_GAP_IN):
    """Tick labels of the SAME axis that overlap or nearly touch each other.
    text_overlaps exempts same-axis tick labels from each other (spec 2);
    this stricter check catches the real collisions among them:
      * x-axis labels (side by side) collide when they overlap vertically
        and their horizontal gap is below `min_gap_in` inches (default 0.03
        in) -- e.g. 'loss' and '0' abutting as 'loss0' on a narrow
        horizontal IRR axis (round 2);
      * y-axis labels (stacked) collide only when they overlap (> tol_px in
        both directions), as before.
    Returns a list of messages; [] = pass."""
    fig.canvas.draw()
    gap = float(min_gap_in) * fig.dpi
    by_axis = {}
    for t, bb, owner, w in _drawn_texts(fig):
        if owner is not None:
            by_axis.setdefault(owner, []).append((t, bb, w))
    msgs = []
    for items in by_axis.values():
        for i in range(len(items)):
            ti, bi, wi = items[i]
            for tj, bj, wj in items[i + 1:]:
                dx = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
                dy = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                if 'xtick' in wi:
                    hit = dy > tol_px and dx > -gap
                else:
                    hit = dx > tol_px and dy > tol_px
                if hit:
                    msgs.append(f'tick labels collide: {ti.get_text()!r} x '
                                f'{tj.get_text()!r} [{wi}] ({dx:.1f} x '
                                f'{dy:.1f} px)')
    return msgs


def check_figure(fig, size=None, raise_on_fail=True, exempt=()):
    """Draw and run the render checks: text_overlaps, tick_label_collisions,
    min_font_check, glyph_check (and the canvas size if `size` = (w, h)
    inches is given). Prints every problem; raises FigureCheckError if any
    (unless raise_on_fail=False). Returns {'overlaps', 'ticks', 'fonts',
    'glyphs', 'size'} (lists of messages)."""
    res = {'overlaps': text_overlaps(fig, exempt=exempt),
           'ticks': tick_label_collisions(fig),
           'fonts': min_font_check(fig), 'glyphs': glyph_check(fig),
           'size': []}
    if size is not None:
        got = tuple(round(float(v), 4) for v in fig.get_size_inches())
        if got != tuple(round(float(v), 4) for v in size):
            res['size'].append(f'canvas {got} in != {tuple(size)} in')
    n = sum(len(v) for v in res.values())
    for k, v in res.items():
        for m in v:
            _safe_print(f'[{k}] {m}')
    print(f'figure checks: {n} problem(s) ('
          + ', '.join(f'{k} {len(v)}' for k, v in res.items()) + ')')
    if n and raise_on_fail:
        raise FigureCheckError(f'{n} figure-check problem(s)')
    return res


# %% Palette / colour-vision check -----------------------------------------------------
# Machado, Oliveira & Fernandes (2009) severity-1.0 matrices (linear RGB),
# the same as the dataviz skill's validator
_MACHADO = {
    'deutan': np.array([[0.367322, 0.860646, -0.227968],
                        [0.280085, 0.672501, 0.047413],
                        [-0.011820, 0.042940, 0.968881]]),
    'protan': np.array([[0.152286, 1.052583, -0.204868],
                        [0.114503, 0.786281, 0.099216],
                        [-0.003882, -0.048116, 1.051998]]),
    'tritan': np.array([[1.255528, -0.076749, -0.178779],
                        [-0.078411, 0.930809, 0.147602],
                        [0.004733, 0.691367, 0.303900]]),
}
_RGB2XYZ = np.array([[0.4124564, 0.3575761, 0.1804375],
                     [0.2126729, 0.7151522, 0.0721750],
                     [0.0193339, 0.1191920, 0.9503041]])
_D65 = np.array([0.95047, 1.0, 1.08883])
# co-occurring pairs (spec section 2) and the grayscale pairs (6.B.9)
CVD_PAIRS = (('unin', 'relay'), ('unin', 'etoh_light'), ('relay', 'ibo_light'),
             ('etoh', 'ibo'), ('adh1', 'adh6'), ('gly', 'tca'),
             ('etoh_light', 'ibo_light'))
CVD_MIN_DE = 12.0
GRAY_PAIRS = {('unin', 'relay'): 15.0, ('etoh', 'ibo'): 20.0}


def _hex2lin(h):
    h = h.lstrip('#')
    c = np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _lin2hex(lin):
    lin = np.clip(lin, 0, 1)
    c = np.where(lin <= 0.0031308, 12.92 * lin,
                 1.055 * lin ** (1 / 2.4) - 0.055)
    return '#' + ''.join(f'{int(round(v * 255)):02X}' for v in c)


def _lin2lab(lin):
    t = (_RGB2XYZ @ lin) / _D65
    d = 6 / 29
    f = np.where(t > d ** 3, np.cbrt(t), t / (3 * d * d) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]),
                     200 * (f[1] - f[2])])


def _lab2lin(lab):
    L, a, b = lab
    fy = (L + 16) / 116
    fx, fz = fy + a / 500, fy - b / 200
    d = 6 / 29
    f = np.array([fx, fy, fz])
    t = np.where(f > d, f ** 3, 3 * d * d * (f - 4 / 29))
    return np.linalg.solve(_RGB2XYZ, t * _D65)


def _sim(h, kind):
    lin = _hex2lin(h)
    return lin if kind is None else np.clip(_MACHADO[kind] @ lin, 0, 1)


def delta_e76(h1, h2, kind=None):
    """CIE76 dE between two hexes, under simulated deutan / protan / tritan
    vision (Machado 2009, severity 1) or normal vision (kind=None)."""
    return float(np.linalg.norm(_lin2lab(_sim(h1, kind))
                                - _lin2lab(_sim(h2, kind))))


def lstar(h):
    """CIELAB L* (the grayscale lightness of a colour)."""
    return float(_lin2lab(_hex2lin(h))[0])


def darken(h, dL=5.0):
    """The hex with CIELAB L* lowered by dL (the spec's fix for a failing
    pair: darken the lighter colour in 5-L* steps)."""
    lab = _lin2lab(_hex2lin(h))
    lab[0] -= dL
    return _lin2hex(_lab2lin(lab))


def _dataviz_validator():
    """The dataviz skill's validate_palette.py, loaded by file path if it is
    available ($DATAVIZ_VALIDATOR or the bundled-skills temp folder)."""
    cands = [os.environ.get('DATAVIZ_VALIDATOR', '')]
    cands += sorted(glob.glob(os.path.join(
        tempfile.gettempdir(), 'claude', 'bundled-skills', '*', '*',
        'dataviz', 'scripts', 'validate_palette.py')), reverse=True)
    for p in cands:
        if p and os.path.exists(p):
            try:
                spec = importlib.util.spec_from_file_location('dv_validate', p)
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                return mod, p
            except Exception:
                continue
    return None, None


def cvd_check(pairs=CVD_PAIRS, min_de=CVD_MIN_DE, gray_pairs=GRAY_PAIRS,
              raise_on_fail=True, verbose=False, dataviz=True, palette=None):
    """Colour-vision check of the co-occurring palette pairs: CIE76 dE >=
    `min_de` under normal vision and simulated deuteranopia / protanopia /
    tritanopia (Machado 2009, severity 1, linear RGB), and the grayscale
    lightness gaps GRAY_PAIRS (|dL*| >= 15 cyan/teal, >= 20 amber/violet).
    If the dataviz skill's validator is found, its OKLab pair deltas are
    reported too (ADVISORY: its lightness band / chroma floor target
    categorical series hues, not this palette's tints and greys).
    Returns (ok, report); raises AssertionError on failure unless
    raise_on_fail=False."""
    P = PALETTE if palette is None else palette
    report, ok = [], True
    for a, b in pairs:
        d = {k or 'normal': delta_e76(P[a], P[b], k)
             for k in (None, 'deutan', 'protan', 'tritan')}
        good = min(d.values()) >= min_de
        ok &= good
        report.append({'pair': (a, b), 'kind': 'cvd', 'ok': good, **d,
                       'min': min(d.values())})
    for (a, b), need in gray_pairs.items():
        dl = abs(lstar(P[a]) - lstar(P[b]))
        good = dl >= need
        ok &= good
        report.append({'pair': (a, b), 'kind': 'gray', 'ok': good,
                       'dL': dl, 'need': need})
    dv = None
    if dataviz:
        mod, path = _dataviz_validator()
        if mod is not None:
            dv = []
            for a, b in pairs:
                cvd = min(mod.deltaE(P[a], P[b], k)
                          for k in ('protan', 'deutan'))
                nor = mod.deltaE(P[a], P[b])
                dv.append({'pair': (a, b), 'cvd_oklab': cvd,
                           'normal_oklab': nor,
                           'ok': cvd >= mod.CVD_TARGET
                           and nor >= mod.NORMAL_FLOOR})
            report.append({'kind': 'dataviz', 'path': path, 'pairs': dv,
                           'ok': all(x['ok'] for x in dv)})
    if verbose:
        print('\npalette check (CIE76 dE, Machado 2009 severity 1; need '
              f'>= {min_de:g}):')
        for r in report:
            if r['kind'] == 'cvd':
                print(f"  {'ok  ' if r['ok'] else 'FAIL'} "
                      f"{r['pair'][0]:>10s} / {r['pair'][1]:<10s} "
                      f"normal {r['normal']:5.1f}  deutan {r['deutan']:5.1f}"
                      f"  protan {r['protan']:5.1f}  tritan "
                      f"{r['tritan']:5.1f}")
            elif r['kind'] == 'gray':
                print(f"  {'ok  ' if r['ok'] else 'FAIL'} grayscale "
                      f"{r['pair'][0]} / {r['pair'][1]}: |dL*| "
                      f"{r['dL']:.1f} (need >= {r['need']:g})")
        if dv is not None:
            print('  dataviz validator (advisory; OKLab x100, CVD target '
                  '8, normal floor 15):')
            for x in dv:
                print(f"    {'ok  ' if x['ok'] else 'warn'} "
                      f"{x['pair'][0]:>10s} / {x['pair'][1]:<10s} "
                      f"cvd {x['cvd_oklab']:5.1f}  normal "
                      f"{x['normal_oklab']:5.1f}")
        elif dataviz:
            print('  dataviz validator: not found (skipped)')
    if not ok and raise_on_fail:
        bad = [r['pair'] for r in report if r['kind'] in ('cvd', 'gray')
               and not r['ok']]
        raise AssertionError(f'palette check failed for {bad}')
    return ok, report


# %% Saving --------------------------------------------------------------------------
def _png_size(path):
    with open(path, 'rb') as fh:
        head = fh.read(24)
    return struct.unpack('>II', head[16:24])


def save(fig, stem, out_dir=OUT_DIR, stamp=None, latest=True, max_pdf_mb=10.0):
    """Write <stem>_<YYYY.MM.DD-HH.MM>.png (300 dpi) and .pdf to `out_dir`,
    plus stable copies <stem>_latest.png / .pdf (overwritten every run) so
    reviewers always find the newest render. Asserts every path <= 259
    characters, the PNG pixel size = canvas x 300 dpi (no cropping), the
    PDF <= max_pdf_mb and free of Type-3 fonts (fonttype 42). Never uses
    bbox_inches='tight'. Prints and returns the paths + the PNG size."""
    stamp = time.strftime('%Y.%m.%d-%H.%M') if stamp is None else stamp
    os.makedirs(out_dir, exist_ok=True)
    out = {}
    for ext in ('png', 'pdf'):
        p = os.path.abspath(os.path.join(out_dir, f'{stem}_{stamp}.{ext}'))
        q = os.path.abspath(os.path.join(out_dir, f'{stem}_latest.{ext}'))
        for path in (p, q):
            assert len(path) <= MAX_PATH_CHARS, \
                f'path {len(path)} > {MAX_PATH_CHARS} chars: {path}'
        fig.savefig(p, dpi=DPI, format=ext)
        out[ext] = p
        if latest:
            shutil.copyfile(p, q)
            out[ext + '_latest'] = q
    w, h = _png_size(out['png'])
    W, H = fig.get_size_inches()
    want = (int(round(W * DPI)), int(round(H * DPI)))
    assert (w, h) == want, f'PNG is {w} x {h} px, expected {want}'
    out['png_px'] = (w, h)
    size_mb = os.path.getsize(out['pdf']) / 1e6
    assert size_mb <= max_pdf_mb, f'PDF is {size_mb:.1f} MB > {max_pdf_mb}'
    with open(out['pdf'], 'rb') as fh:
        pdf = fh.read()
    assert b'/Type3' not in pdf, 'PDF holds Type-3 fonts (fonttype 42 off?)'
    out['pdf_mb'] = size_mb
    print(f'saved {out["png"]}  ({w} x {h} px)')
    print(f'saved {out["pdf"]}  ({size_mb:.2f} MB)')
    if latest:
        print(f'  latest copies: {os.path.basename(out["png_latest"])}, '
              f'{os.path.basename(out["pdf_latest"])}')
    return out

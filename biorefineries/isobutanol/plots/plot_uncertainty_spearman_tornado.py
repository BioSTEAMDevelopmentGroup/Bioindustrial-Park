#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Spearman tornado figure from an uncertainty-analysis results workbook
(*_1_full_evaluation.xlsx as written by analyses/full/uncertainties_IBO_EtOH.py):
Spearman's rank correlation coefficient (rho) of five outcomes with the
uncertain parameters that matter to each, one tornado per outcome:

  left column (economics)      A  MESP (purity-adjusted ethanol MPSP)
                               B  total capital investment (TCI)
  right column (fermentation)  C  ethanol titer
                               D  ethanol productivity
                               E  ethanol yield

A parameter is shown in a panel only if its rho with that panel's outcome is
significant: |rho| >= RHO_MIN (0.10) and p < P_MAX (0.05). With n = 6,000
samples the p cut is not the binding one (|rho| >= 0.10 has p ~ 1e-14). Bars
are sorted by |rho|, largest on top, and coloured by parameter category
(CATEGORIES: prices and tax, plant scale and uptime, corn composition,
fermentation kinetics); the light grey band marks |rho| < RHO_MIN. rho is read
from the workbook's 'Spearman' / 'Spearman p-values' sheets (biosteam
Model.spearman_r over all samples); MESP's rho is that of the $/kg MPSP column,
which a $/GGE conversion (a positive constant) leaves unchanged. A cross-check
recomputes each shown rho from the 'Parameters' and 'TEA results' sheets.

Sim-safe: pure pandas/matplotlib/scipy, never imports biorefineries.

Usage:
    python plot_uncertainty_spearman_tornado.py [<results.xlsx>]
        [--out-dir DIR] [--stem STEM] [--dpi N]

With no workbook, the newest scenario-A *_1_full_evaluation.xlsx under
analyses/results is used. Also writes <stem>.csv: every shown bar.
"""

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import TICKDOWN, TICKLEFT
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator
from scipy import stats

__all__ = ('plot_uncertainty_spearman_tornado',)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(os.path.dirname(HERE), 'analyses', 'results')
DEFAULT_OUT_DIR = os.path.join(RESULTS_DIR, 'publication', 'Uncertainty')

RHO_MIN = 0.10 # |rho| at or above it ...
P_MAX = 0.05   # ... and p below it -> shown

# (letter, title, workbook metric column), column by column
PANELS = (
    ('A', 'MESP', 'Purity-adjusted ethanol MPSP [$/kg]'),
    ('B', 'Total capital investment', 'Total capital investment [10^6 $]'),
    ('C', 'Ethanol titer', 'Et OH titer [g-EtOH/L-water]'),
    ('D', 'Ethanol productivity', 'Et OH productivity [g-EtOH/L-water/h]'),
    ('E', 'Ethanol yield', 'Et OH yield [g-EtOH/g-sugars-added]'),
)
COLUMNS = (('A', 'B'), ('C', 'D', 'E'))

# categorical hues, in this fixed order (validated with the dataviz skill's
# validate_palette.js on the white surface: lightness, chroma, CVD and
# normal-vision checks pass; the gold is below 3:1 contrast, relieved by the
# per-bar parameter labels and rho values)
CATEGORIES = (
    ('prices', 'Prices and tax', '#e0703e'),
    ('scale', 'Plant scale and uptime', '#5a6bcc'),
    ('feedstock', 'Corn composition', '#c99a1e'),
    ('kinetics', 'Fermentation kinetics', '#00949a'),
    ('other', 'Process design', '#8a8a8a'),
)
CATEGORY_COLOR = {key: color for key, _, color in CATEGORIES}
INK = '#0b0b0b'
INSIGNIFICANT_BAND = '#ececec'
GRID_COLOR = '#d9d9d9'

# display names of the non-kinetic parameters (workbook name without units)
PARAMETER_LABELS = {
    'Feedstock unit price': 'Corn price',
    'DDGS unit price (sale)': 'DDGS price',
    'Plant annual operating days': 'Annual operating days',
    'Feedstock capacity': 'Plant capacity (corn feed rate)',
    'Feedstock starch content': 'Corn starch content',
    'Federal corporate tax rate': 'Corporate income tax rate',
}
SCALE_PARAMETERS = ('Plant annual operating days', 'Feedstock capacity')
# short descriptions of the kinetic parameters (roles as in
# docs/reports/kinetic-parameters-scenario-A.md); unlisted ones show the
# symbol alone
KINETIC_DESCRIPTIONS = {
    'k_1l': 'Low-affinity glucose uptake',
    'k_1e': 'Acetaldehyde-enhanced glycolysis',
    'K_1i': 'Glycolysis signal saturation',
    'k_1ie': 'Ethanol inhibition of glycolysis',
    'k_4ie': 'Ethanol inhibition of acetate formation',
    'k_7': 'Growth capacity',
    'k_7ie': 'Ethanol inhibition of growth',
    'k_10': 'Active-biomass decay',
}

FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'label': 12, 'axis_title': 13, 'panel_title': 13, 'value': 10, 'legend': 11}
TICK_LEN = {'major': 4.0, 'minor': 2.0} # pt; left/bottom ticks extend this far in AND out

# layout, inches
BAR_PITCH_IN = 0.36
BAR_HEIGHT = 0.64 # of the pitch
PANEL_PAD_IN = 0.16 # above the first and below the last bar
AXES_W_IN = 3.3
LABEL_W_IN = 3.35 # room for the parameter labels left of each column
COL_GAP_IN = 0.45
PANEL_GAP_IN = 0.62 # between panels of a column (room for the title)
MARGINS_IN = {'top': 0.45, 'bottom': 0.75, 'right': 0.25}
XLIM = (-1.0, 1.0)
VALUE_INSIDE_ABOVE = 0.75 # |rho| above it: value printed inside the bar end


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
    # way. Call after fig.canvas.draw() so every tick object exists. y has no
    # minor ticks (categorical).
    for which, L in TICK_LEN.items():
        ax.tick_params(axis='x', which=which, direction='inout', length=2*L,
                       top=True, labelsize=FONTS['tick'])
        get = 'get_major_ticks' if which == 'major' else 'get_minor_ticks'
        for tick in getattr(ax.xaxis, get)():
            tick.tick2line.set_marker(TICKDOWN)
            tick.tick2line.set_markersize(L)
    ax.tick_params(axis='y', which='major', direction='inout', length=2*TICK_LEN['major'],
                   right=True, labelsize=FONTS['label'])
    for tick in ax.yaxis.get_major_ticks():
        tick.tick2line.set_marker(TICKLEFT)
        tick.tick2line.set_markersize(TICK_LEN['major'])


#%% Data

def newest_results_file(scenario='A'):
    files = [f for f in glob.glob(os.path.join(glob.escape(RESULTS_DIR), '*_1_full_evaluation.xlsx'))
             if f"['{scenario}']" in os.path.basename(f)]
    if not files:
        raise FileNotFoundError(f'no scenario-{scenario} *_1_full_evaluation.xlsx in {RESULTS_DIR}')
    return max(files, key=os.path.getmtime)


def _strip_units(name):
    return name.split(' [')[0].strip()


def kinetic_symbol(name):
    """'. k 1ie [g_per_l_per_h]' -> 'k_1ie'; 'K 1i [g_per_l]' -> 'K_1i'."""
    m = re.match(r'^\.?\s*([kK])\s+(\S+)', _strip_units(name))
    return f'{m.group(1)}_{m.group(2)}' if m else None


def parameter_category(element, name):
    if element.startswith('Kinetics'):
        return 'kinetics'
    base = _strip_units(name)
    if base in SCALE_PARAMETERS:
        return 'scale'
    if element == 'Feedstock':
        return 'feedstock'
    if element == 'TEA':
        return 'prices'
    return 'other'


def parameter_label(element, name):
    if element.startswith('Kinetics'):
        symbol = kinetic_symbol(name)
        if symbol:
            letter, sub = symbol.split('_', 1)
            math = rf'$\it{{{letter}}}_{{\mathrm{{{sub}}}}}$'
            desc = KINETIC_DESCRIPTIONS.get(symbol)
            return f'{desc}, {math}' if desc else math
    base = _strip_units(name)
    return PARAMETER_LABELS.get(base, base)


def load_spearman(results_file):
    """rho and p frames indexed by (element, parameter), the blank parameter
    dropped."""
    frames = []
    for sheet in ('Spearman', 'Spearman p-values'):
        df = pd.read_excel(results_file, sheet_name=sheet, header=0)
        df['Element'] = df['Element'].ffill()
        df = df[~df['Parameter'].str.startswith('Blank parameter')]
        frames.append(df.set_index(['Element', 'Parameter']))
    return frames


def significant(rho, p, metric):
    r, pv = rho[metric].astype(float), p[metric].astype(float)
    keep = (r.abs() >= RHO_MIN) & (pv < P_MAX)
    out = pd.DataFrame({'rho': r[keep], 'p': pv[keep]})
    return out.reindex(out['rho'].abs().sort_values(ascending=False).index)


def recompute_rho(results_file, shown):
    """Max |rho_sheet - rho_recomputed| over the shown bars, from the raw
    sample sheets (None if the sheets cannot be aligned)."""
    try:
        params = pd.read_excel(results_file, sheet_name='Parameters', header=[0, 1], index_col=0)
        tea = pd.read_excel(results_file, sheet_name='TEA results', header=[0, 1], index_col=0)
    except Exception as e:
        print(f'Cross-check skipped: {e}')
        return None
    pcols = {c[1]: c for c in params.columns}
    tcols = {c[1]: c for c in tea.columns}
    worst = 0.0
    for metric, table in shown.items():
        if metric not in tcols:
            return None
        y = tea[tcols[metric]].astype(float)
        for (element, name), row in table.iterrows():
            if name not in pcols:
                return None
            x = params[pcols[name]].astype(float)
            ok = np.isfinite(x) & np.isfinite(y)
            worst = max(worst, abs(stats.spearmanr(x[ok], y[ok])[0] - row['rho']))
    return worst, len(tea)


#%% Drawing

def draw_tornado(ax, table, letter, title):
    n = len(table)
    y = np.arange(n)[::-1] # largest |rho| on top
    ax.axvspan(-RHO_MIN, RHO_MIN, color=INSIGNIFICANT_BAND, lw=0, zorder=0)
    for xv in (-0.5, 0.5):
        ax.axvline(xv, color=GRID_COLOR, lw=0.6, zorder=0.5)
    colors = [CATEGORY_COLOR[parameter_category(e, nm)] for e, nm in table.index]
    ax.barh(y, table['rho'], height=BAR_HEIGHT, color=colors, lw=0, zorder=2)
    ax.axvline(0, color=INK, lw=0.9, zorder=3)
    for yi, r in zip(y, table['rho']):
        text = f'{r:.2f}'.replace('-', '−')
        inside = abs(r) > VALUE_INSIDE_ABOVE
        sign = np.sign(r)
        x = r - sign*0.02 if inside else r + sign*0.02
        ha = ('right' if r > 0 else 'left') if inside else ('left' if r > 0 else 'right')
        ax.text(x, yi, text, ha=ha, va='center', fontsize=FONTS['value'],
                color='white' if inside else INK,
                fontweight='bold' if inside else 'normal', zorder=4)
    ax.set_yticks(y, [parameter_label(e, nm) for e, nm in table.index])
    ax.set_ylim(-0.5 - PANEL_PAD_IN/BAR_PITCH_IN, n - 0.5 + PANEL_PAD_IN/BAR_PITCH_IN)
    ax.set_xlim(*XLIM)
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_title(rf'$\mathbf{{{letter}}}$   ' + r'$\mathbf{' + title.replace(' ', r'\ ') + '}$',
                 loc='left', fontsize=FONTS['panel_title'], pad=7,
                 x=-LABEL_W_IN/AXES_W_IN)
    for side in ('left', 'right', 'top', 'bottom'):
        ax.spines[side].set_linewidth(0.8)


def plot_uncertainty_spearman_tornado(results_file=None, out_dir=DEFAULT_OUT_DIR,
                                      stem=None, dpi=600):
    results_file = results_file or newest_results_file('A')
    print(f'Results: {results_file}')
    rho, p = load_spearman(results_file)
    tables = {}
    for letter, title, metric in PANELS:
        if metric not in rho.columns:
            raise KeyError(f'metric {metric!r} not in the Spearman sheet')
        tables[letter] = significant(rho, p, metric)
        t = tables[letter]
        print(f'{letter} {title}: {len(t)} parameters with |rho| >= {RHO_MIN} and p < {P_MAX}')
        for (e, nm), row in t.iterrows():
            print(f'    {row["rho"]:+.3f}  (p = {row["p"]:.1e})  {e}: {nm}')
    check = recompute_rho(results_file, {m: tables[l] for l, _, m in PANELS})
    if check:
        print(f'Cross-check over {check[1]} samples: max |rho_sheet - rho_recomputed| = {check[0]:.2e}')

    apply_font_rcparams()
    m = MARGINS_IN
    heights = {l: len(tables[l])*BAR_PITCH_IN + 2*PANEL_PAD_IN for l, _, _ in PANELS}
    col_heights = [sum(heights[l] for l in col) + (len(col) - 1)*PANEL_GAP_IN for col in COLUMNS]
    width = len(COLUMNS)*(LABEL_W_IN + AXES_W_IN) + (len(COLUMNS) - 1)*COL_GAP_IN + m['right']
    height = m['top'] + max(col_heights) + m['bottom']
    fig = plt.figure(figsize=(width, height))
    titles = {l: t for l, t, _ in PANELS}
    axes, bottom_axes = [], []
    for j, col in enumerate(COLUMNS):
        x0 = j*(LABEL_W_IN + AXES_W_IN + COL_GAP_IN) + LABEL_W_IN
        top = height - m['top']
        for k, letter in enumerate(col):
            h = heights[letter]
            ax = fig.add_axes([x0/width, (top - h)/height, AXES_W_IN/width, h/height])
            draw_tornado(ax, tables[letter], letter, titles[letter])
            if k < len(col) - 1:
                ax.tick_params(axis='x', labelbottom=False)
            else:
                # plain text: mathtext sets the apostrophe as a prime
                ax.set_xlabel('Spearman’s ρ', fontsize=FONTS['axis_title'],
                              fontweight='bold')
                bottom_axes.append((ax, top - h))
            axes.append(ax)
            top -= h + PANEL_GAP_IN

    # category legend under the shorter column, plus the grey band
    present = {parameter_category(e, nm) for t in tables.values() for e, nm in t.index}
    handles = [Patch(facecolor=c, lw=0, label=name) for key, name, c in CATEGORIES if key in present]
    handles.append(Patch(facecolor=INSIGNIFICANT_BAND, edgecolor='#bdbdbd', lw=0.6,
                         label=rf'$|\rho|$ < {RHO_MIN:.2f} (not shown)'))
    j_short = int(np.argmin(col_heights))
    ax_short, y_short = bottom_axes[j_short]
    x_left = j_short*(LABEL_W_IN + AXES_W_IN + COL_GAP_IN)
    fig.legend(handles=handles, loc='upper left', frameon=False, fontsize=FONTS['legend'],
               handlelength=1.4, handleheight=1.0, borderaxespad=0,
               bbox_to_anchor=((x_left + 0.25)/width, (y_short - 0.75)/height))

    fig.canvas.draw()
    for ax in axes:
        style_ticks(ax)

    os.makedirs(out_dir, exist_ok=True)
    if stem is None:
        tag = os.path.basename(results_file).split('_A_1_full_evaluation')[0]
        tag = tag.lstrip('_').replace("['A']_", '').replace('IBO_', '')
        stem = f'spearman_tornado_A_{tag}'
    paths = []
    for ext in ('png', 'pdf'):
        path = os.path.join(out_dir, f'{stem}.{ext}')
        fig.savefig(path, dpi=dpi, facecolor='white')
        paths.append(path)
        print(f'Saved {path}')
    plt.close(fig)
    rows = [dict(panel=l, metric=mt, element=e, parameter=nm,
                 label=parameter_label(e, nm), category=parameter_category(e, nm),
                 rho=r['rho'], p=r['p'])
            for l, _, mt in PANELS for (e, nm), r in tables[l].iterrows()]
    path = os.path.join(out_dir, f'{stem}.csv')
    pd.DataFrame(rows).to_csv(path, index=False)
    print(f'Saved {path}')
    return paths + [path]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('results_file', nargs='?', default=None)
    parser.add_argument('--out-dir', default=DEFAULT_OUT_DIR)
    parser.add_argument('--stem', default=None)
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    plot_uncertainty_spearman_tornado(args.results_file, args.out_dir, args.stem, args.dpi)


if __name__ == '__main__':
    main()

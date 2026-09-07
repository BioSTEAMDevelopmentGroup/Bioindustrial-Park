#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Grouped heatmap of a kinetic-BO study's best-objective trial against the
scenario-A baseline, in ABSOLUTE values: four strips (ethanol-production
rate constants, Ehrlich-branch rate constants, inhibition coefficients,
feeding parameters), each a 2 x N grid whose top sub-row is the scenario-A
baseline and bottom sub-row the best trial, colored by value on a log
scale with one sequential colorbar per strip. A zero (the Ehrlich block
under scenario A) has no log value and is drawn as a light-gray "0" cell.
The best-trial cells also carry the fold change vs. the baseline in small
italics.

Inhibition coefficients are shown per member (the applied_<member> CSV
columns of a grouped study), not per group multiplier. A member absent
from the scenario-A workbook (the isobutanol coefficients and the r16
terms) took its scenario-B value as the study baseline (the ethanol_
isobutanol preset samples the B workbook's rows starting at A); those
baseline cells are tagged "(B)".

Sim-safe: kinetic_optimization.py is loaded by file path (no biosteam
import, no load()), the trajectory CSV and the parameter-distribution
workbooks are plain file reads. Run:

    python plots/plot_kin_opt_best_vs_baseline_heatmap.py [STUDY_NAME]

STUDY_NAME defaults to the metabolic_minimal_subset IRR study. Writes
<study>_best_vs_baseline_heatmap_<stamp>.png/.pdf next to the CSV.
"""
import os
import sys
import importlib.util
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
from matplotlib.ticker import (LogLocator, FixedLocator, FuncFormatter,
                               NullFormatter)

# --- paths -------------------------------------------------------------------
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
DEFAULT_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset_irr'
                 '_rb0.001-10_ib0.2-2_burden')

_spec = importlib.util.spec_from_file_location(
    'ko', os.path.join(PKG_DIR, 'kinetic_optimization.py'))
ko = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ko)

# --- scenario-A baseline operating point (CLAUDE.md "Baseline feeding
# strategies"; hardcoded in every analysis script, no central constant) --------
BASELINE_A = {'threshold_conc': 217.125, 'target_conc': 221.25,
              'spike_conc': 600.0, 'max_n_spikes': 16, 'n_glu_spikes': 10,
              'IRR': 0.1225}

# --- typography (plot-formatting preferences: Arial, 12/12/10/9 hierarchy) ---
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'panel_title': 12, 'cbar_title': 12, 'cell': 10,
         'cell_sub': 8, 'suptitle': 12}


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


# --- color: one single-hue sequential map (light -> dark blue), shared by
# every strip; each strip gets its own log norm over its own values ----------
CMAP = LinearSegmentedColormap.from_list(
    'value', ['#e4edf7', '#9dbfe0', '#4a86c5', '#1f4e9c', '#0b2a5c'])
ZERO_COLOR = '#f2f2f2'
CMAP.set_bad(ZERO_COLOR)


def text_color_on(rgba):
    r, g, b = to_rgb(rgba)
    lum = 0.2126*r + 0.7152*g + 0.0722*b
    return 'white' if lum < 0.45 else '#1a1a1a'


def sub(name):
    """k_1e -> $k_{1e}$ style label; plain names pass through."""
    if '_' in name and name[0] in 'kK':
        head, tail = name.split('_', 1)
        return rf'$\mathit{{{head}}}_{{\mathrm{{{tail}}}}}$'
    return name


def fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return 'n/a'
    if float(v).is_integer() and abs(v) < 1e4:
        return f'{int(v)}'
    return f'{v:.3g}'


# --- data ---------------------------------------------------------------------
def load_best(study_name):
    csv_path = os.path.join(RESULTS_DIR, f'{study_name}_trajectory.csv')
    df = pd.read_csv(csv_path)
    ok = df[df['state'] == 'COMPLETE']
    best = ok.loc[ok['objective'].idxmax()]
    return best, csv_path


def build_rows(best, A, B):
    """Each row: dict(title, cbar_label, cells=[(label, baseline, best,
    tag)]) with `tag` a short marker printed under the baseline value
    ('(B)' = the baseline came from the scenario-B workbook)."""
    def rate(name):
        return (sub(name), A[name], float(best[name]), None)

    row_etoh = dict(
        title='Ethanol-production rate constants (glycolysis, PDC, ADH)',
        cbar_label='rate constant',
        cells=[rate('k_1l'), rate('k_1h'), rate('k_1e'), rate('k_3'),
               rate('k_6')])

    row_ehr = dict(
        title='Ehrlich-branch rate constants (branch off in scenario A)',
        cbar_label='rate constant',
        cells=[(sub(n), 0.0, float(best[n]), None)
               for n in ('k_13', 'k_14', 'k_15', 'k_16')])

    inh = []
    for gname, members in ko.METABOLIC_MINIMAL_SUBSET_GROUPS.items():
        for m in members:
            in_A = A.get(m) is not None
            base = A[m] if in_A else B[m]
            inh.append((sub(m), base, float(best[f'applied_{m}']),
                        None if in_A else '(B)'))
    row_inh = dict(
        title='Inhibition coefficients (one multiplier per effector family '
              'in the study; "(B)": absent in A, baseline = scenario B value)',
        cbar_label='inhibition coefficient', cells=inh, extra_bottom=0.3,
        group_spans=[(g, len(ms)) for g, ms in
                     ko.METABOLIC_MINIMAL_SUBSET_GROUPS.items()])

    threshold = float(best['threshold_conc'])
    target = min(ko.TARGET_CONC_MAX, threshold + float(best['target_delta']))
    row_feed = dict(
        title='Feeding parameters (spike feed pinned at 600 '
              r'$\mathrm{g·L}^{-1}$)',
        cbar_label=r'$\mathrm{g·L}^{-1}$ (conc.) / count (spikes)',
        cells=[
            ('threshold sugar\n' + r'conc. [$\mathrm{g·L}^{-1}$]',
             BASELINE_A['threshold_conc'], threshold, None),
            ('target sugar\n' + r'conc. [$\mathrm{g·L}^{-1}$]',
             BASELINE_A['target_conc'], target, None),
            ('max. glucose\nspikes', BASELINE_A['max_n_spikes'],
             float(best['max_n_spikes']), None),
            ('actual glucose\nspikes', BASELINE_A['n_glu_spikes'],
             float(best['n_glu_spikes']), None),
        ])
    return [row_etoh, row_ehr, row_inh, row_feed]


# --- figure -------------------------------------------------------------------
def plot(rows, best, study_name, out_stem):
    apply_font_rcparams()
    strip_w, cell_h = 10.0, 0.62               # inches; every strip spans strip_w
    left, cbar_gap, cbar_w, right = 1.35, 0.3, 0.18, 1.3
    title_h, gap = 0.4, 0.6                    # per-strip title band, inter-strip gap
    top_pad, bottom_pad = 0.6, 0.75
    fig_w = left + strip_w + cbar_gap + cbar_w + right
    n_sub = 2
    fig_h = (top_pad + len(rows)*(title_h + n_sub*cell_h + gap) - gap
             + bottom_pad + sum(r.get('extra_bottom', 0.0) for r in rows))
    fig = plt.figure(figsize=(fig_w, fig_h))
    sub_labels = ['scenario A\nbaseline',
                  f"best-IRR\ntrial {int(best['trial_number'])}"]

    y = fig_h - top_pad
    for row in rows:
        cells = row['cells']
        n = len(cells)
        data = np.array([[c[1] for c in cells], [c[2] for c in cells]],
                        dtype=float)
        positive = data[data > 0]
        norm = LogNorm(positive.min(), positive.max())
        masked = np.ma.masked_less_equal(data, 0.0)
        y -= title_h + n_sub*cell_h
        ax = fig.add_axes([left/fig_w, y/fig_h, strip_w/fig_w,
                           n_sub*cell_h/fig_h])
        ax.imshow(masked, cmap=CMAP, norm=norm, aspect='auto',
                  extent=(0, n, 0, n_sub), interpolation='nearest',
                  origin='upper')
        for i in range(1, n):                       # 2px surface gaps
            ax.axvline(i, color='white', lw=2)
        ax.axhline(1, color='white', lw=2)
        # Cells are pure color (values/multipliers read off the colorbars);
        # the only in-cell mark kept is the "(B)" baseline-source tag.
        for i, (label, base, bv, tag) in enumerate(cells):
            if tag:
                tc = (text_color_on(CMAP(norm(base))) if base > 0
                      else '#666666')
                ax.text(i + 0.5, n_sub - 0.5, tag, ha='center', va='center',
                        color=tc, fontsize=FONTS['cell_sub'])
        ax.set_xlim(0, n)
        ax.set_ylim(0, n_sub)
        ax.set_yticks([1.5, 0.5])
        ax.set_yticklabels(sub_labels, fontsize=FONTS['tick'])
        ax.set_xticks(np.arange(n) + 0.5)
        ax.set_xticklabels([c[0] for c in cells], fontsize=FONTS['tick'])
        ax.tick_params(axis='both', length=0, pad=4)
        # effector-family brackets under the inhibition members
        spans = row.get('group_spans')
        if spans:
            x0 = 0
            for gname, width in spans:
                x1 = x0 + width
                ax.plot([x0 + 0.08, x1 - 0.08], [-0.55, -0.55], color='#444444',
                        lw=0.9, clip_on=False)
                ax.text((x0 + x1)/2, -0.66,
                        gname.replace('inhib_', '') + ' inhibition',
                        ha='center', va='top', fontsize=FONTS['tick'] - 2,
                        color='#333333', clip_on=False)
                x0 = x1
        for s in ax.spines.values():
            s.set_edgecolor('#444444')
            s.set_linewidth(0.8)
        ax.set_title(row['title'], fontsize=FONTS['panel_title'],
                     fontweight='bold', loc='left', pad=6)
        # colorbar: this strip's own log norm
        cax = fig.add_axes([(left + strip_w + cbar_gap)/fig_w, y/fig_h,
                            cbar_w/fig_w, n_sub*cell_h/fig_h])
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=CMAP)
        cb = fig.colorbar(sm, cax=cax)
        lo, hi = norm.vmin, norm.vmax
        decades = [10.0**k for k in range(int(np.floor(np.log10(lo))),
                                          int(np.ceil(np.log10(hi))) + 1)
                   if lo <= 10.0**k <= hi]
        if len(decades) >= 2:
            cb.ax.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
        else:   # under a decade of range: label the 1-2-5 steps instead
            steps = [m*10.0**k for k in range(int(np.floor(np.log10(lo))),
                                              int(np.ceil(np.log10(hi))) + 1)
                     for m in (1, 2, 5) if lo <= m*10.0**k <= hi]
            cb.ax.yaxis.set_major_locator(FixedLocator(steps))
        cb.ax.yaxis.set_minor_locator(
            LogLocator(base=10, subs=(2, 3, 4, 5, 6, 7, 8, 9), numticks=20))
        cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: fmt(v)))
        cb.ax.yaxis.set_minor_formatter(NullFormatter())
        cb.ax.tick_params(which='major', labelsize=FONTS['cell_sub'] + 1,
                          length=3)
        cb.ax.tick_params(which='minor', length=1.5)
        cb.set_label(row['cbar_label'], fontsize=FONTS['cell_sub'] + 1)
        cb.outline.set_linewidth(0.8)
        y -= gap + row.get('extra_bottom', 0.0)

    fig.suptitle(
        f"Best-IRR trial {int(best['trial_number'])} "
        f"(IRR {best['IRR']:.3f}; isobutanol {best['IBO titer']:.0f}, "
        f"ethanol {best['EtOH titer']:.0f} "
        r"$\mathrm{g·L}^{-1}$" + ") vs. scenario A baseline "
        f"(IRR {BASELINE_A['IRR']:.3f})",
        fontsize=FONTS['suptitle'], fontweight='bold', y=1 - 0.2/fig_h)
    fig.text(0.5, 0.01, study_name, ha='center', va='bottom', fontsize=7,
             color='#777777')
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out_stem}.{ext}', dpi=300, bbox_inches='tight')
    return fig


def main(study_name=DEFAULT_STUDY):
    best, csv_path = load_best(study_name)
    A = ko.workbook_kinetic_baselines('A')
    B = ko.workbook_kinetic_baselines('B')
    rows = build_rows(best, A, B)
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    out_stem = os.path.join(RESULTS_DIR,
                            f'{study_name}_best_vs_baseline_heatmap_{stamp}')
    plot(rows, best, study_name, out_stem)
    print(f'best trial {int(best["trial_number"])}: IRR {best["IRR"]:.4f}')
    for row in rows:
        print(f'  [{row["title"]}]')
        for label, base, bv, tag in row['cells']:
            print(f'    {label!r:45s} A {base:10.4g}  best {bv:10.4g}  '
                  f'{tag or ""}')
    print(f'wrote {out_stem}.png / .pdf')
    return out_stem


if __name__ == '__main__':
    main(*sys.argv[1:2])

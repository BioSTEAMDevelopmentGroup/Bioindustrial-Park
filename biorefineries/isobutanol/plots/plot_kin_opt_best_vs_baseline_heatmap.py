#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Grouped heatmap of a kinetic-BO study's best-objective trial against the
scenario-A baseline: four single-row strips (ethanol-production rate
constants, Ehrlich-branch rate constants, inhibition-effector multipliers,
feeding parameters), one diverging colorbar per strip, cell color =
log10(best / baseline), cell text = the best value over the baseline.

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
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.ticker import FixedLocator, FuncFormatter

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
              'spike_conc': 600.0, 'max_n_spikes': 16, 'IRR': 0.1225}

# --- typography (plot-formatting preferences: Arial, 12/12/10/9 hierarchy) ---
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'panel_title': 12, 'cbar_title': 12, 'cell': 10,
         'cell_sub': 9, 'suptitle': 12}


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


# --- color: one diverging map (blue = below baseline, gray midpoint, red =
# above), shared by every strip; each strip gets its own symmetric norm -------
CMAP = LinearSegmentedColormap.from_list(
    'fold_change', ['#1f4e9c', '#7fa7d8', '#d9d9d9', '#e8967a', '#a8281c'])


def text_color_on(rgba):
    r, g, b = to_rgb(rgba)
    lum = 0.2126*r + 0.7152*g + 0.0722*b
    return 'white' if lum < 0.5 else '#1a1a1a'


def sub(name):
    """k_1e -> $k_{1e}$ style label; plain names pass through."""
    if '_' in name and name[0] in 'kK':
        head, tail = name.split('_', 1)
        return rf'$\mathit{{{head}}}_{{\mathrm{{{tail}}}}}$'
    return name


def fold_ticks(L):
    """Round fold-change tick values inside the symmetric log10 range +-L
    (at most 7; the finer 0.5/0.7/1.5/2 steps only for narrow ranges)."""
    coarse = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100]
    fine = [0.5, 0.7, 1.5, 2]
    within = lambda vals: [v for v in vals if abs(np.log10(v)) <= L + 1e-9]
    ticks = sorted(within(coarse) + within(fine))
    if len(ticks) > 7:
        ticks = within(coarse)
    return ticks


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
    """Each row: dict(title, cells=[(label, best, baseline, ref_note)],
    cbar_label). Color = log10(best/baseline)."""
    def rate(name):
        return (sub(name), float(best[name]), A[name], None)

    row_etoh = dict(
        title='Ethanol-production rate constants (glycolysis, PDC, ADH)',
        cells=[rate('k_1l'), rate('k_1h'), rate('k_1e'), rate('k_3'),
               rate('k_6')],
        cbar_label='fold change vs. scenario A')

    # Ehrlich block: scenario A = 0 (branch off), so the fold change is
    # taken against the scenario-B values and the cell says so.
    ehr = []
    for name in ('k_13', 'k_14', 'k_15', 'k_16'):
        ehr.append((sub(name), float(best[name]), B[name], 'A: 0'))
    row_ehr = dict(
        title='Ehrlich-branch rate constants (scenario A baseline = 0; '
              'colored vs. scenario B values)',
        cells=ehr, cbar_label='fold change vs. scenario B')

    groups = ko.METABOLIC_MINIMAL_SUBSET_GROUPS
    inh = []
    for gname, members in groups.items():
        effector = gname.replace('inhib_', '')
        mem = ', '.join(sub(m) for m in members)
        inh.append((f'{effector}\ninhibition', float(best[gname]), 1.0, mem))
    row_inh = dict(
        title='Inhibition-effector multipliers (one per family; members '
              'scale together, baseline 1×)',
        cells=inh, cbar_label='multiplier vs. scenario A',
        extra_bottom=0.3)

    threshold = float(best['threshold_conc'])
    target = min(ko.TARGET_CONC_MAX, threshold + float(best['target_delta']))
    feed = [
        ('threshold sugar\n' + r'conc. [$\mathrm{g·L}^{-1}$]', threshold,
         BASELINE_A['threshold_conc'], None),
        ('target sugar\n' + r'conc. [$\mathrm{g·L}^{-1}$]', target,
         BASELINE_A['target_conc'], None),
        ('max. glucose\nspikes', float(best['max_n_spikes']),
         BASELINE_A['max_n_spikes'], None),
    ]
    row_feed = dict(
        title='Feeding parameters (spike feed pinned at 600 '
              r'$\mathrm{g·L}^{-1}$; actual spikes: '
              f"{int(best['n_glu_spikes'])} vs. 10 at baseline)",
        cells=feed, cbar_label='fold change vs. scenario A')
    return [row_etoh, row_ehr, row_inh, row_feed]


# --- figure -------------------------------------------------------------------
def plot(rows, best, study_name, out_stem):
    apply_font_rcparams()
    n_max = max(len(r['cells']) for r in rows)
    cell_w, cell_h = 1.55, 0.95          # inches
    left, right_cbar = 0.25, 0.55        # inches
    title_h, gap = 0.42, 0.55            # per-strip title band, inter-strip gap
    top_pad, bottom_pad = 0.55, 0.95     # suptitle band, bottom tick-label band
    fig_w = left + n_max*cell_w + 0.35 + right_cbar + 1.35
    fig_h = (top_pad + len(rows)*(title_h + cell_h + gap) - gap + bottom_pad
             + sum(r.get('extra_bottom', 0.0) for r in rows))
    fig = plt.figure(figsize=(fig_w, fig_h))

    y = fig_h - top_pad
    for row in rows:
        cells = row['cells']
        n = len(cells)
        vals = np.array([c[1]/c[2] for c in cells])
        logs = np.log10(vals)
        L = max(0.2, float(np.ceil(np.max(np.abs(logs))*2)/2))
        norm = Normalize(-L, L)
        # heatmap strip
        y -= title_h + cell_h
        ax = fig.add_axes([left/fig_w, y/fig_h, n*cell_w/fig_w, cell_h/fig_h])
        ax.imshow(logs[None, :], cmap=CMAP, norm=norm, aspect='auto',
                  extent=(0, n, 0, 1), interpolation='nearest')
        # 2px surface gaps between cells
        for i in range(1, n):
            ax.axvline(i, color='white', lw=2)
        for i, (label, bv, base, note) in enumerate(cells):
            tc = text_color_on(CMAP(norm(logs[i])))
            mult = vals[i]
            ax.text(i + 0.5, 0.66, f'{fmt(bv)}', ha='center', va='center',
                    color=tc, fontsize=FONTS['cell'], fontweight='bold')
            base_txt = (f'A: {fmt(base)}' if note is None
                        else (f'{note}, B: {fmt(base)}' if note.startswith('A:')
                              else f'A: {fmt(base)}'))
            ax.text(i + 0.5, 0.36, base_txt, ha='center', va='center',
                    color=tc, fontsize=FONTS['cell_sub'])
            ax.text(i + 0.5, 0.12, f'{mult:.2g}×', ha='center', va='center',
                    color=tc, fontsize=FONTS['cell_sub'], style='italic')
        ax.set_xlim(0, n)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks(np.arange(n) + 0.5)
        ax.set_xticklabels([c[0] for c in cells], fontsize=FONTS['tick'])
        ax.tick_params(axis='x', length=0, pad=4)
        # member lists under the inhibition-group labels
        if any(c[3] and not c[3].startswith('A:') for c in cells):
            for i, c in enumerate(cells):
                ax.text(i + 0.5, -0.66, c[3], ha='center', va='top',
                        fontsize=8, color='#555555', transform=ax.transData)
        for s in ax.spines.values():
            s.set_edgecolor('#444444')
            s.set_linewidth(0.8)
        ax.set_title(row['title'], fontsize=FONTS['panel_title'],
                     fontweight='bold', loc='left', pad=6)
        # colorbar: its own symmetric norm, ticks labelled as fold changes
        cax = fig.add_axes([(left + n_max*cell_w + 0.35)/fig_w, y/fig_h,
                            0.18/fig_w, cell_h/fig_h])
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=CMAP)
        cb = fig.colorbar(sm, cax=cax)
        cb.locator = FixedLocator(np.log10(fold_ticks(L)))
        cb.formatter = FuncFormatter(lambda v, _: f'{10**v:.2g}×')
        cb.update_ticks()
        cb.ax.tick_params(labelsize=FONTS['cell_sub'], length=3)
        cb.set_label(row['cbar_label'], fontsize=FONTS['cell_sub'])
        cb.outline.set_linewidth(0.8)
        y -= gap + row.get('extra_bottom', 0.0)

    fig.suptitle(
        f"Best-IRR trial {int(best['trial_number'])} "
        f"(IRR {best['IRR']:.3f}; isobutanol {best['IBO titer']:.0f}, "
        f"ethanol {best['EtOH titer']:.0f} "
        r"$\mathrm{g·L}^{-1}$" + ") vs. scenario A baseline "
        f"(IRR {BASELINE_A['IRR']:.3f})",
        fontsize=FONTS['suptitle'], fontweight='bold', y=1 - 0.18/fig_h)
    fig.text(0.5, 0.012, study_name, ha='center', va='bottom', fontsize=7,
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
        for label, bv, base, note in row['cells']:
            print(f'    {label!r:45s} best {bv:10.4g}  base {base:10.4g}  '
                  f'x{bv/base:.3g}')
    print(f'wrote {out_stem}.png / .pdf')
    return out_stem


if __name__ == '__main__':
    main(*sys.argv[1:2])

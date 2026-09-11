#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""IRR vs total-capital-investment (TCI) Pareto frontier of a kinetic-
optimization study. One panel:

  * every COMPLETE trial with finite IRR and TCI as a light scatter cloud,
    coloured by isobutanol titer (the co-production story -- how much IBO a
    given (TCI, IRR) point makes);
  * the Pareto-optimal frontier -- the trials that are not dominated when
    IRR is MAXIMIZED and TCI is MINIMIZED -- highlighted and joined by a
    line, each labelled with its trial number;
  * the IRR = 0 break-even reference and, optionally, the scenario-A
    baseline point (hard-coded outcomes, simulated 2026-09-07) as a star,
    to show where the incumbent design sits relative to the frontier.

Sim-safe: a plain pandas read of one trajectory CSV under analyses/results.
No biosteam import, no load(); runnable while a study is in flight. Run:

    python plots/plot_kin_opt_pareto_irr_tci.py [--study <study name or CSV>]

With no --study it uses the 2026-09-07 metabolic_minimal_subset IRR study
(the run that found the high-isobutanol optimum, trial 774). Writes
<stem>_pareto_irr_tci_<stamp>.png and .pdf to --out-dir (analyses/results).
"""
import os
import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')

DEFAULT_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset'
                 '_irr_rb0.001-10_ib0.2-2_burden')

# Scenario-A baseline outcomes -- HARD-CODED (simulated 2026-09-07,
# smoke_test_1 protocol, IBO_2026), matching plot_kin_opt_parameter_sets.py.
BASELINE_A = {'IRR': 0.1230, 'TCI': 139.6, 'IBO titer': 0.0}

FONTS = {'tick': 9, 'legend': 10, 'axis': 12, 'title': 12, 'callout': 9}


def apply_fonts():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['xtick.labelsize'] = FONTS['tick']
    plt.rcParams['ytick.labelsize'] = FONTS['tick']
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def style_ticks(ax):
    """Ticks on all four sides, major + minor; top/right inward, left/bottom
    in-and-out (house preference)."""
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', top=True, right=True, length=4)
    ax.tick_params(which='minor', top=True, right=True, length=2)
    ax.tick_params(axis='x', which='major', bottom=True, direction='inout',
                   length=4)
    ax.tick_params(axis='x', which='minor', bottom=True, direction='inout',
                   length=2)
    ax.tick_params(axis='y', which='major', left=True, direction='inout',
                   length=4)
    ax.tick_params(axis='y', which='minor', left=True, direction='inout',
                   length=2)
    # top/right stay inward: redraw them as such (they inherit 'inout' above,
    # so re-set just those sides)
    ax.tick_params(which='major', top=True, right=True)
    ax.tick_params(which='minor', top=True, right=True)


def resolve_csv(study):
    """Accept a study name, a bare filename, or a path; return the CSV path."""
    if os.path.isfile(study):
        return study
    cand = study if study.endswith('.csv') else study + '_trajectory.csv'
    path = cand if os.path.isabs(cand) else os.path.join(RESULTS_DIR, cand)
    if not os.path.isfile(path):
        raise FileNotFoundError(f'no trajectory CSV for study {study!r} '
                                f'(looked for {path})')
    return path


def pareto_mask(irr, tci):
    """Boolean mask of non-dominated points: MAXIMIZE irr, MINIMIZE tci.
    Point i is dominated if some j is >= in IRR and <= in TCI and strictly
    better in at least one."""
    irr = np.asarray(irr, float)
    tci = np.asarray(tci, float)
    n = irr.size
    keep = np.ones(n, bool)
    for i in range(n):
        better = ((irr >= irr[i]) & (tci <= tci[i])
                  & ((irr > irr[i]) | (tci < tci[i])))
        if better.any():
            keep[i] = False
    return keep


def load_complete(csv_path):
    df = pd.read_csv(csv_path)
    d = df[(df['state'] == 'COMPLETE')
           & np.isfinite(df['IRR']) & np.isfinite(df['TCI'])].copy()
    if d.empty:
        raise ValueError(f'no COMPLETE finite-(IRR, TCI) trials in {csv_path}')
    return d


def make_figure(csv_path, show_baseline=True):
    d = load_complete(csv_path)
    irr = d['IRR'].to_numpy()
    tci = d['TCI'].to_numpy()
    ibo = d['IBO titer'].to_numpy()
    trial = d['trial_number'].to_numpy()

    keep = pareto_mask(irr, tci)
    front = d[keep].sort_values('TCI')

    apply_fonts()
    fig, ax = plt.subplots(figsize=(6.4, 5.0))

    # cloud, coloured by isobutanol titer
    sc = ax.scatter(tci, irr, c=ibo, cmap='viridis', s=22, alpha=0.55,
                    linewidths=0, zorder=2, vmin=0,
                    vmax=float(np.nanmax(ibo)))

    # break-even reference (label parked at the right, clear of the frontier)
    ax.axhline(0.0, color='0.55', lw=0.9, ls=(0, (5, 4)), zorder=1)
    ax.annotate('IRR = 0 (break-even)', xy=(1.0, 0.0),
                xycoords=('axes fraction', 'data'), xytext=(-4, 3),
                textcoords='offset points', va='bottom', ha='right',
                fontsize=FONTS['callout'], color='0.4')

    # Pareto frontier: joined line + outlined markers
    fx = front['TCI'].to_numpy()
    fy = front['IRR'].to_numpy()
    ax.plot(fx, fy, color='crimson', lw=1.6, zorder=3, solid_capstyle='round')
    ax.scatter(fx, fy, s=68, facecolor='white', edgecolor='crimson',
               linewidths=1.8, zorder=4)
    # stagger the trial labels above/below so close frontier points
    # (e.g. #400 and #1242) do not collide
    label_offsets = {925: (7, -3, 'left', 'top'),
                     1591: (7, -3, 'left', 'top'),
                     400: (-6, -4, 'right', 'top'),
                     1242: (7, 2, 'left', 'bottom'),
                     774: (8, 0, 'left', 'center')}
    for xi, yi, ti in zip(fx, fy, front['trial_number'].to_numpy()):
        dx, dy, ha, va = label_offsets.get(int(ti), (7, -3, 'left', 'top'))
        ax.annotate(f'#{int(ti)}', xy=(xi, yi), xytext=(dx, dy),
                    textcoords='offset points', ha=ha, va=va,
                    fontsize=FONTS['callout'], color='crimson')

    # scenario-A baseline
    handles = [
        Line2D([], [], marker='o', ls='none', mfc='0.7', mec='0.7', ms=6,
               label='Optimization trial (colour = IBO titer)'),
        Line2D([], [], marker='o', ls='-', color='crimson', mfc='white',
               mec='crimson', mew=1.8, ms=8, label='Pareto frontier'),
    ]
    if show_baseline:
        ax.scatter([BASELINE_A['TCI']], [BASELINE_A['IRR']], marker='*',
                   s=230, facecolor='#1f77b4', edgecolor='k', linewidths=0.8,
                   zorder=5)
        ax.annotate('Scenario-A\nbaseline', xy=(BASELINE_A['TCI'],
                    BASELINE_A['IRR']), xytext=(8, 8),
                    textcoords='offset points', ha='left', va='bottom',
                    fontsize=FONTS['callout'], color='#1f77b4')
        handles.append(Line2D([], [], marker='*', ls='none', mfc='#1f77b4',
                       mec='k', ms=13, label='Scenario-A baseline'))

    ax.set_xlabel(r'Total capital investment (MM\$)', fontsize=FONTS['axis'])
    ax.set_ylabel('Internal rate of return (IRR)', fontsize=FONTS['axis'])
    ax.set_title('IRR – capital Pareto frontier', fontsize=FONTS['title'])
    style_ticks(ax)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r'Isobutanol titer ($\mathrm{g·L}^{-1}$)',
                   fontsize=FONTS['axis'])
    cbar.ax.tick_params(labelsize=FONTS['tick'])

    ax.legend(handles=handles, loc='lower right', fontsize=FONTS['legend'],
              frameon=True, framealpha=0.9, edgecolor='0.8',
              handletextpad=0.5)
    fig.tight_layout()
    return fig, front


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study', default=DEFAULT_STUDY,
                   help='study name, bare CSV filename, or path')
    p.add_argument('--out-dir', default=RESULTS_DIR)
    p.add_argument('--no-baseline', action='store_true',
                   help='omit the scenario-A baseline star')
    args = p.parse_args()

    csv_path = resolve_csv(args.study)
    fig, front = make_figure(csv_path, show_baseline=not args.no_baseline)

    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    stem = stem.replace('_trajectory', '')
    base = os.path.join(args.out_dir, f'{stem}_pareto_irr_tci_{stamp}')
    os.makedirs(args.out_dir, exist_ok=True)
    fig.savefig(base + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    plt.close(fig)

    print('Pareto frontier (IRR maximized, TCI minimized):')
    cols = ['trial_number', 'IRR', 'TCI', 'IBO titer', 'EtOH titer', 'tau']
    print(front[cols].to_string(index=False))
    print('\nWrote:')
    print('  ' + base + '.png')
    print('  ' + base + '.pdf')


if __name__ == '__main__':
    main()

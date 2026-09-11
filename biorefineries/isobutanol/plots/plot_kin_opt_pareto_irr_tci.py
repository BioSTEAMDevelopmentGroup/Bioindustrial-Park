#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Two-objective Pareto-frontier scatter of a kinetic-optimization study.
One panel: every COMPLETE trial (finite IRR and TCI) as a scatter cloud
coloured by a third metric, with the non-dominated frontier drawn as a red
dashed line and each frontier trial labelled by number. The scenario-A
baseline (hard-coded outcomes, simulated 2026-09-07) is marked with a star,
and an IRR = 0 break-even reference is drawn whenever IRR is an axis.

Selectable via --view:

  irr_tci  (default)  x = total capital investment (minimized)
                      y = IRR (maximized), colour = isobutanol titer.
  ibo_irr             x = isobutanol titer (maximized)
                      y = IRR (maximized), colour = total capital investment.

Sim-safe: a plain pandas read of one trajectory CSV under analyses/results.
No biosteam import, no load(); runnable while a study is in flight. Run:

    python plots/plot_kin_opt_pareto_irr_tci.py [--view ibo_irr] \
        [--study <study name or CSV>]

With no --study it uses the 2026-09-07 metabolic_minimal_subset IRR study
(the run that found the high-isobutanol optimum, trial 774). Writes
<stem>_pareto_<view>_<stamp>.png and .pdf to --out-dir (analyses/results).
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
BASELINE_A = {'IRR': 0.1230, 'TCI': 139.6, 'IBO titer': 0.0, 'IBO yield': 0.0}

_TCI_LABEL = r'Total capital investment (MM\$)'
_IRR_LABEL = 'Internal rate of return (IRR)'
_IBO_LABEL = r'Isobutanol titer ($\mathrm{g·L}^{-1}$)'
_IBOY_LABEL = r'Isobutanol yield ($\mathrm{g·g}^{-1}$)'

# One entry per --view. `dir` is the optimization direction of each axis
# ('max'/'min'); the frontier is the set non-dominated under those. `cvmin`
# pins the colorbar floor (None = autoscale from data). `legend_loc` is the
# in-axes legend anchor.
VIEWS = {
    'irr_tci': dict(
        x='TCI', xdir='min', xlabel=_TCI_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
        color='IBO titer', clabel=_IBO_LABEL, cshort='IBO titer', cvmin=0.0,
        title='IRR – capital Pareto frontier', legend_loc='lower right'),
    'ibo_irr': dict(
        x='IBO titer', xdir='max', xlabel=_IBO_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
        color='TCI', clabel=_TCI_LABEL, cshort='TCI', cvmin=None,
        title='IRR – isobutanol titer Pareto frontier',
        legend_loc='lower right'),
    'iboyield_irr': dict(
        x='IBO yield', xdir='max', xlabel=_IBOY_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
        color='TCI', clabel=_TCI_LABEL, cshort='TCI', cvmin=None,
        title='IRR – isobutanol yield Pareto frontier',
        legend_loc='lower right'),
}

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
    # top/right stay inward (they inherit 'inout' above, so re-set them)
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


def _to_maximize(values, direction):
    """Sign-flip so larger is always better."""
    v = np.asarray(values, float)
    return v if direction == 'max' else -v


def pareto_mask(x, xdir, y, ydir):
    """Non-dominated mask for two objectives with per-axis directions.
    Point i is dominated if some j is at least as good on both objectives
    and strictly better on one."""
    a = _to_maximize(x, xdir)
    b = _to_maximize(y, ydir)
    n = a.size
    keep = np.ones(n, bool)
    for i in range(n):
        better = ((a >= a[i]) & (b >= b[i]) & ((a > a[i]) | (b > b[i])))
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


def make_figure(csv_path, view, show_baseline=True, show_frontier=True):
    v = VIEWS[view]
    d = load_complete(csv_path)
    x = d[v['x']].to_numpy()
    y = d[v['y']].to_numpy()
    c = d[v['color']].to_numpy()

    keep = pareto_mask(x, v['xdir'], y, v['ydir'])
    # order the frontier along x so the connecting line is monotone
    front = d[keep].sort_values(v['x'])

    apply_fonts()
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    # a little breathing room so edge frontier labels/stars are not clipped
    ax.margins(x=0.07, y=0.08)

    sc = ax.scatter(x, y, c=c, cmap='viridis', s=22, alpha=0.55,
                    linewidths=0, zorder=2,
                    vmin=v['cvmin'] if v['cvmin'] is not None
                    else float(np.nanmin(c)),
                    vmax=float(np.nanmax(c)))

    # break-even reference whenever IRR is an axis
    if v['y'] == 'IRR':
        ax.axhline(0.0, color='0.55', lw=0.9, ls=(0, (5, 4)), zorder=1)
        ax.annotate('IRR = 0 (break-even)', xy=(1.0, 0.0),
                    xycoords=('axes fraction', 'data'), xytext=(-4, 3),
                    textcoords='offset points', va='bottom', ha='right',
                    fontsize=FONTS['callout'], color='0.4')
    elif v['x'] == 'IRR':
        ax.axvline(0.0, color='0.55', lw=0.9, ls=(0, (5, 4)), zorder=1)

    handles = [
        Line2D([], [], marker='o', ls='none', mfc='0.7', mec='0.7', ms=6,
               label=f"Optimization trial (colour = {v['cshort']})"),
    ]

    # Pareto frontier: red dashed line, no markers; trial labels staggered
    # above/below along x, and pulled inward near an axis edge, so nothing
    # collides or clips
    if show_frontier:
        fx = front[v['x']].to_numpy()
        fy = front[v['y']].to_numpy()
        ax.plot(fx, fy, color='crimson', lw=1.6, ls=(0, (6, 4)), zorder=3,
                solid_capstyle='round')
        x0, x1 = ax.get_xlim()
        span = x1 - x0
        for k, (xi, yi, ti) in enumerate(
                zip(fx, fy, front['trial_number'].to_numpy())):
            frac = (xi - x0) / span if span else 0.5
            if frac > 0.88:            # near right edge: label to the left
                dx, ha = -7, 'right'
            elif frac < 0.12:          # near left edge: label to the right
                dx, ha = 7, 'left'
            else:
                dx, ha = 0, 'center'
            up = (k % 2 == 0)
            ax.annotate(f'#{int(ti)}', xy=(xi, yi),
                        xytext=(dx, 9 if up else -9),
                        textcoords='offset points', ha=ha,
                        va='bottom' if up else 'top',
                        fontsize=FONTS['callout'], color='crimson')
        handles.append(Line2D([], [], ls=(0, (6, 4)), color='crimson',
                       lw=1.6, label='Pareto frontier'))
    if show_baseline:
        bx, by = BASELINE_A[v['x']], BASELINE_A[v['y']]
        ax.scatter([bx], [by], marker='*', s=230, facecolor='#1f77b4',
                   edgecolor='k', linewidths=0.8, zorder=5)
        ax.annotate('Scenario-A\nbaseline', xy=(bx, by), xytext=(8, 8),
                    textcoords='offset points', ha='left', va='bottom',
                    fontsize=FONTS['callout'], color='#1f77b4')
        handles.append(Line2D([], [], marker='*', ls='none', mfc='#1f77b4',
                       mec='k', ms=13, label='Scenario-A baseline'))

    ax.set_xlabel(v['xlabel'], fontsize=FONTS['axis'])
    ax.set_ylabel(v['ylabel'], fontsize=FONTS['axis'])
    title = v['title'] if show_frontier \
        else v['title'].replace('Pareto frontier', 'trade-off')
    ax.set_title(title, fontsize=FONTS['title'])
    style_ticks(ax)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(v['clabel'], fontsize=FONTS['axis'])
    cbar.ax.tick_params(labelsize=FONTS['tick'])

    ax.legend(handles=handles, loc=v['legend_loc'], fontsize=FONTS['legend'],
              frameon=True, framealpha=0.9, edgecolor='0.8',
              handletextpad=0.5)
    fig.tight_layout()
    return fig, front


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--view', default='irr_tci', choices=sorted(VIEWS),
                   help='which two objectives / colour metric to plot')
    p.add_argument('--study', default=DEFAULT_STUDY,
                   help='study name, bare CSV filename, or path')
    p.add_argument('--out-dir', default=RESULTS_DIR)
    p.add_argument('--no-baseline', action='store_true',
                   help='omit the scenario-A baseline star')
    p.add_argument('--no-frontier', action='store_true',
                   help='draw only the scatter cloud, no Pareto frontier')
    args = p.parse_args()

    v = VIEWS[args.view]
    csv_path = resolve_csv(args.study)
    fig, front = make_figure(csv_path, args.view,
                             show_baseline=not args.no_baseline,
                             show_frontier=not args.no_frontier)

    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    stem = stem.replace('_trajectory', '')
    base = os.path.join(args.out_dir, f'{stem}_pareto_{args.view}_{stamp}')
    os.makedirs(args.out_dir, exist_ok=True)
    fig.savefig(base + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    plt.close(fig)

    _dir = {'max': 'maximized', 'min': 'minimized'}
    print(f"Pareto frontier ({v['x']} {_dir[v['xdir']]}, "
          f"{v['y']} {_dir[v['ydir']]}):")
    cols = ['trial_number', 'IRR', 'TCI', 'IBO titer', 'EtOH titer', 'tau']
    print(front[cols].to_string(index=False))
    print('\nWrote:')
    print('  ' + base + '.png')
    print('  ' + base + '.pdf')


if __name__ == '__main__':
    main()

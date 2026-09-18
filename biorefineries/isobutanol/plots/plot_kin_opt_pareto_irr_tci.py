#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Two-objective Pareto-frontier scatter across kinetic-optimization
campaigns. The COMPLETE trials (finite IRR and TCI) of one or more studies
are POOLED into one scatter cloud, each point coloured by which campaign it
came from (the same per-campaign hues as plot_kin_opt_parameter_sets.py);
the non-dominated frontier of the pool is drawn as a red dashed line with
each frontier point labelled <campaign>#<trial>. The scenario-A baseline
(hard-coded outcomes, simulated 2026-09-07) is a grey star, and an IRR = 0
break-even reference is drawn whenever IRR is an axis.

Selectable via --view (both axes; colour is always the campaign):

  irr_tci  (default)  x = total capital investment (minimized)
                      y = IRR (maximized).
  ibo_irr             x = isobutanol titer (maximized), y = IRR (maximized).
  iboyield_irr        x = isobutanol yield (maximized), y = IRR (maximized).

Sim-safe: plain pandas reads of trajectory CSVs under analyses/results. No
biosteam import, no load(); runnable while a study is in flight. Run:

    python plots/plot_kin_opt_pareto_irr_tci.py [--view ibo_irr] \
        [--studies <study name or CSV> ...]

With no --studies it pools the five metabolic_minimal_subset campaigns (one
per objective: IRR / ethanol titer / isobutanol titer / ethanol yield /
isobutanol yield) plotted by default in plot_kin_opt_parameter_sets.py, so
the frontier is the best achievable ACROSS all five objectives. Writes
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

# The five metabolic_minimal_subset campaigns plotted by default -- one per
# optimized objective -- matching plot_kin_opt_parameter_sets.py's default
# set. Pooling their COMPLETE trials gives the achievable frontier ACROSS all
# five objectives (a richer, cross-campaign frontier than any single study).
_MINIMAL_SUBSET_STUDY = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset'
                         '_%s_rb0.001-10_ib0.2-2_burden')
# objective slug -> short tag used to label pooled frontier points
STUDY_TAGS = {'irr': 'IRR', 'etoh_titer': 'EtT', 'ibo_titer': 'IBT',
              'ibo_yield': 'IBY', 'etoh_yield': 'EtY'}
DEFAULT_STUDIES = [_MINIMAL_SUBSET_STUDY % o for o in STUDY_TAGS]

# Per-campaign categorical colours, tag -> colour, matching the hue order
# plot_kin_opt_parameter_sets.py assigns its default sets (baseline grey,
# then irr / ibo_titer / etoh_titer / ibo_yield / etoh_yield taking
# HUE_COLORS[0..4]); trials are coloured by which campaign they came from.
CAMPAIGN_ORDER = ['IRR', 'IBT', 'EtT', 'IBY', 'EtY']
CAMPAIGN_COLORS = {'IRR': '#18C4DC', 'IBT': '#f98f60', 'EtT': '#79bf82',
                   'IBY': '#a280b9', 'EtY': '#f3c354'}
CAMPAIGN_LABELS = {'IRR': 'IRR', 'IBT': 'IBO titer', 'EtT': 'EtOH titer',
                   'IBY': 'IBO yield', 'EtY': 'EtOH yield'}
_UNKNOWN_CAMPAIGN_COLOR = '0.6'
# scenario-A baseline marker colour (plot_kin_opt_parameter_sets.py's
# BASELINE_COLOR) -- a neutral grey, distinct from the campaign hues
BASELINE_COLOR = '#90918e'

# Scenario-A baseline outcomes -- HARD-CODED (simulated 2026-09-07,
# smoke_test_1 protocol, IBO_2026), matching plot_kin_opt_parameter_sets.py.
BASELINE_A = {'IRR': 0.1230, 'TCI': 139.6, 'IBO titer': 0.0, 'IBO yield': 0.0}

_TCI_LABEL = r'Total capital investment (MM\$)'
_IRR_LABEL = 'Internal rate of return (IRR)'
_IBO_LABEL = r'Isobutanol titer ($\mathrm{g·L}^{-1}$)'
_IBOY_LABEL = r'Isobutanol yield ($\mathrm{g·g}^{-1}$)'

# One entry per --view. `xdir`/`ydir` are the optimization directions
# ('max'/'min'); the frontier is the set non-dominated under those.
# `legend_loc` is the in-axes legend anchor. Points are always coloured by
# source campaign (see CAMPAIGN_COLORS), so there is no colour metric here.
VIEWS = {
    'irr_tci': dict(
        x='TCI', xdir='min', xlabel=_TCI_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
        title='IRR – capital Pareto frontier', legend_loc='lower right'),
    'ibo_irr': dict(
        x='IBO titer', xdir='max', xlabel=_IBO_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
        title='IRR – isobutanol titer Pareto frontier',
        legend_loc='lower right'),
    'iboyield_irr': dict(
        x='IBO yield', xdir='max', xlabel=_IBOY_LABEL,
        y='IRR', ydir='max', ylabel=_IRR_LABEL,
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


def _tukey_upper(values, k=3.0):
    """Upper Tukey fence q3 + k*IQR -- a robust display cap. Pooling the
    ethanol-titer campaign brings in a detached high-TCI artifact cluster
    (~400+ MM\\$, non-physical high-ethanol-titer designs) that would
    otherwise blow out a TCI axis or colorbar; capping at this fence keeps
    the economically meaningful region legible while those points stay on
    the plot (clipped to the edge / saturated colour)."""
    q1, q3 = np.nanpercentile(values, [25, 75])
    return q3 + k * (q3 - q1)


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


def study_tag(study):
    """Short label for a study's source objective, matched from its name;
    falls back to a truncated basename for an unrecognized study."""
    base = os.path.basename(str(study))
    for slug, tag in STUDY_TAGS.items():
        if f'_{slug}_' in base or base.endswith(f'_{slug}'):
            return tag
    return base[:6]


def load_pool(studies):
    """Pool the COMPLETE finite-(IRR, TCI) trials of one or more studies into
    a single frame, tagging each row with its source campaign in 'study'."""
    frames = []
    for s in studies:
        csv_path = resolve_csv(s)
        df = pd.read_csv(csv_path)
        d = df[(df['state'] == 'COMPLETE')
               & np.isfinite(df['IRR']) & np.isfinite(df['TCI'])].copy()
        if d.empty:
            raise ValueError(f'no COMPLETE finite-(IRR, TCI) trials in '
                             f'{csv_path}')
        d['study'] = study_tag(s)
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def _place_frontier_labels(ax, front, v, n_studies):
    """Label each frontier point <campaign>#<trial> (plain #<trial> for a
    single campaign) in a de-collided vertical column on the side away from
    the cloud -- left for a min-x frontier, right for a max-x frontier --
    with a thin leader to the point. Frontier points often crowd a narrow
    band, so labels are spread to a minimum vertical gap and joined to their
    point by a line."""
    fx = front[v['x']].to_numpy()
    fy = front[v['y']].to_numpy()
    trials = front['trial_number'].to_numpy()
    tags = front['study'].to_numpy()
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xspan, yspan = x1 - x0, y1 - y0
    # open a gutter on the side the frontier hugs (left for a min-x frontier,
    # right for a max-x frontier) and drop the label column into it, so the
    # labels never land on the cloud or on the colorbar
    left = (v['xdir'] == 'min')
    ha = 'right' if left else 'left'
    gutter = 0.22 * xspan
    if left:
        ax.set_xlim(left=x0 - gutter)
        text_x = x0 - 0.03 * xspan
    else:
        ax.set_xlim(right=x1 + gutter)
        text_x = x1 + 0.03 * xspan

    order = np.argsort(fy)
    fracs = [(fy[i] - y0) / yspan for i in order]
    min_gap = 0.045
    for i in range(1, len(fracs)):
        if fracs[i] - fracs[i - 1] < min_gap:
            fracs[i] = fracs[i - 1] + min_gap
    overflow = fracs[-1] - 0.985
    if overflow > 0:
        fracs = [f - overflow for f in fracs]

    for rank, i in enumerate(order):
        lbl = (f'{tags[i]}#{int(trials[i])}' if n_studies > 1
               else f'#{int(trials[i])}')
        ax.annotate(lbl, xy=(fx[i], fy[i]),
                    xytext=(text_x, y0 + fracs[rank] * yspan),
                    textcoords='data', ha=ha, va='center',
                    fontsize=FONTS['callout'], color='crimson',
                    arrowprops=dict(arrowstyle='-', color='crimson',
                                    lw=0.5, alpha=0.6,
                                    shrinkA=1, shrinkB=1))


def make_figure(pool, view, show_baseline=True, show_frontier=True):
    v = VIEWS[view]
    d = pool
    n_studies = d['study'].nunique()
    x = d[v['x']].to_numpy()
    y = d[v['y']].to_numpy()

    keep = pareto_mask(x, v['xdir'], y, v['ydir'])
    # order the frontier along x so the connecting line is monotone
    front = d[keep].sort_values(v['x'])

    apply_fonts()
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    # a little breathing room so edge frontier labels/stars are not clipped
    ax.margins(x=0.07, y=0.08)

    # colour each point by its source campaign (categorical). Plot the
    # larger campaigns first so the sparse ones stay visible on top; build
    # the legend in the canonical CAMPAIGN_ORDER regardless of plot order.
    tags_present = list(d['study'].unique())
    ordered = ([t for t in CAMPAIGN_ORDER if t in tags_present]
               + [t for t in tags_present if t not in CAMPAIGN_ORDER])
    by_size = sorted(ordered, key=lambda t: int((d['study'] == t).sum()),
                     reverse=True)
    for t in by_size:
        sub = d[d['study'] == t]
        ax.scatter(sub[v['x']], sub[v['y']],
                   color=CAMPAIGN_COLORS.get(t, _UNKNOWN_CAMPAIGN_COLOR),
                   s=18, alpha=0.55, linewidths=0, zorder=2)
    campaign_handles = [
        Line2D([], [], marker='o', ls='none', ms=7, mec='none',
               mfc=CAMPAIGN_COLORS.get(t, _UNKNOWN_CAMPAIGN_COLOR),
               label=CAMPAIGN_LABELS.get(t, t)) for t in ordered]

    # break-even reference whenever IRR is an axis
    if v['y'] == 'IRR':
        ax.axhline(0.0, color='0.55', lw=0.9, ls=(0, (5, 4)), zorder=1)
        ax.annotate('IRR = 0 (break-even)', xy=(1.0, 0.0),
                    xycoords=('axes fraction', 'data'), xytext=(-4, 3),
                    textcoords='offset points', va='bottom', ha='right',
                    fontsize=FONTS['callout'], color='0.4')
    elif v['x'] == 'IRR':
        ax.axvline(0.0, color='0.55', lw=0.9, ls=(0, (5, 4)), zorder=1)

    handles = list(campaign_handles)

    # Pareto frontier: red dashed line, no markers. The trial labels are
    # placed further below (once the axis limits are final) so they can be
    # de-collided against the finished scale.
    if show_frontier:
        fx = front[v['x']].to_numpy()
        fy = front[v['y']].to_numpy()
        ax.plot(fx, fy, color='crimson', lw=1.6, ls=(0, (6, 4)), zorder=3,
                solid_capstyle='round')
        handles.append(Line2D([], [], ls=(0, (6, 4)), color='crimson',
                       lw=1.6, label='Pareto frontier'))
    if show_baseline:
        bx, by = BASELINE_A[v['x']], BASELINE_A[v['y']]
        ax.scatter([bx], [by], marker='*', s=230, facecolor=BASELINE_COLOR,
                   edgecolor='k', linewidths=0.8, zorder=6)
        ax.annotate('Scenario-A\nbaseline', xy=(bx, by), xytext=(8, 8),
                    textcoords='offset points', ha='left', va='bottom',
                    fontsize=FONTS['callout'], color='0.35')
        handles.append(Line2D([], [], marker='*', ls='none', mfc=BASELINE_COLOR,
                       mec='k', ms=13, label='Scenario-A baseline'))

    # robust cap for a TCI x-axis: crop the detached high-TCI artifact
    # cluster so the meaningful region fills the panel (points stay drawn,
    # clipped at the right edge; a note reports how many are off-scale)
    if v['x'] == 'TCI':
        cap = _tukey_upper(x)
        if cap < float(np.nanmax(x)):
            need = [cap]
            if show_baseline:
                need.append(BASELINE_A['TCI'])
            if show_frontier and len(front):
                need.append(float(front['TCI'].max()))
            xhi = max(need) * 1.03
            n_off = int((x > xhi).sum())
            ax.set_xlim(right=xhi)
            if n_off:
                ax.annotate(f'+{n_off} trials off-scale\n(TCI up to '
                            f'{float(np.nanmax(x)):.0f})', xy=(0.985, 0.5),
                            xycoords='axes fraction', ha='right', va='center',
                            fontsize=FONTS['callout'], color='0.45')

    # frontier labels last, so they de-collide against the finished limits
    if show_frontier and len(front):
        _place_frontier_labels(ax, front, v, n_studies)

    ax.set_xlabel(v['xlabel'], fontsize=FONTS['axis'])
    ax.set_ylabel(v['ylabel'], fontsize=FONTS['axis'])
    title = v['title'] if show_frontier \
        else v['title'].replace('Pareto frontier', 'trade-off')
    ax.set_title(title, fontsize=FONTS['title'])
    style_ticks(ax)

    ax.legend(handles=handles, loc=v['legend_loc'], fontsize=FONTS['legend'],
              frameon=True, framealpha=0.9, edgecolor='0.8',
              handletextpad=0.5, title='Campaign', title_fontsize=FONTS['legend'])
    fig.tight_layout()
    return fig, front


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--view', default='irr_tci', choices=sorted(VIEWS),
                   help='which two objectives to plot (colour is the campaign)')
    p.add_argument('--studies', nargs='+', default=DEFAULT_STUDIES,
                   metavar='STUDY',
                   help='one or more study names, bare CSV filenames, or '
                        'paths; their COMPLETE trials are pooled. Default: '
                        'the five metabolic_minimal_subset campaigns.')
    p.add_argument('--out-dir', default=RESULTS_DIR)
    p.add_argument('--no-baseline', action='store_true',
                   help='omit the scenario-A baseline star')
    p.add_argument('--no-frontier', action='store_true',
                   help='draw only the scatter cloud, no Pareto frontier')
    args = p.parse_args()

    v = VIEWS[args.view]
    pool = load_pool(args.studies)
    fig, front = make_figure(pool, args.view,
                             show_baseline=not args.no_baseline,
                             show_frontier=not args.no_frontier)

    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    if len(args.studies) == 1:
        stem = os.path.splitext(os.path.basename(
            resolve_csv(args.studies[0])))[0].replace('_trajectory', '')
    else:
        stem = ('kin_opt_ethanol_isobutanol_metabolic_minimal_subset'
                f'_{len(args.studies)}campaigns')
    base = os.path.join(args.out_dir, f'{stem}_pareto_{args.view}_{stamp}')
    os.makedirs(args.out_dir, exist_ok=True)
    fig.savefig(base + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(base + '.pdf', bbox_inches='tight')
    plt.close(fig)

    _dir = {'max': 'maximized', 'min': 'minimized'}
    print(f"Pooled {len(args.studies)} campaign(s), {len(pool)} COMPLETE "
          f"trials. Pareto frontier ({v['x']} {_dir[v['xdir']]}, "
          f"{v['y']} {_dir[v['ydir']]}):")
    cols = ['study', 'trial_number', 'IRR', 'TCI', 'IBO titer', 'EtOH titer',
            'tau']
    print(front[cols].to_string(index=False))
    print('\nWrote:')
    print('  ' + base + '.png')
    print('  ' + base + '.pdf')


if __name__ == '__main__':
    main()

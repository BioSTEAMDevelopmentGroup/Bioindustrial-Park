#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Objective-landscape figure, fitness-distance layout: the profitability
index (PI) has a spiky ("needle") landscape over the metabolic_split_12d
decision space, whereas isobutanol yield -- the best process-level objective
-- has a smooth landscape whose high region passes through the most profitable
designs.

Data: the 13,622 unique COMPLETE trials pooled over the seven 2026-09-23
seed-350 `ethanol_isobutanol` x `metabolic_split_12d` GP campaigns (one per
objective: PI (log-tail), isobutanol / ethanol yield, titer, productivity),
read through the shared data layer `plots/_objective_landscape_data.py`
(loaded by file path). Distances are Euclidean in the campaigns' own 12-D
unit cube (log-scale axes log-normalized, linear axes min-max) -- the distance
the GP itself sees.

Panels (one row; (a) and (b) the same size, (c) narrower):

  (a) PI of every trial vs its unit-cube distance from the most profitable
      trial (ibo_yield campaign #842, PI 0.504). Above the panel: the
      fitness-distance correlation (FDC; Jones & Forrest 1995), here the
      Spearman rank correlation of the objective with the distance from that
      objective's own best trial over all trials (strongly negative = the
      objective falls off steadily with distance from its optimum = a
      searchable landscape; ~0 = the optimum is a needle).
      Y-SCALING: raw PI spans -5.56 to 0.50 (median -1.47), so on a linear
      axis the 166 trials with PI > 0 would occupy the top 8 % of the panel.
      The axis instead uses the monotone "log-tail" scale the PI campaign
      itself optimized -- PI for PI >= 0, -log(1 - PI) for PI < 0 (C1 at 0,
      so the axis is linear above 0 and only the loss tail is compressed;
      flagged by the axis title "(log scale below 0)") -- with ticks
      labelled in RAW PI. No trial is clipped or dropped (the full range
      maps to [-1.88, 0.50]); `--pi-scale linear` draws the plain linear
      axis instead. The FDC is rank-based, so it does not depend on this
      choice.
  (b) Isobutanol yield of every trial vs its distance from the highest-
      isobutanol-yield trial (ibo_yield #841, 0.381 g/g; its own PI is -0.68;
      it lies 0.40 from #842), same encoding; the trials with PI > 0 are
      highlighted to show where the profitable designs sit on this
      landscape. The 40 ethanol-only PI > 0 designs at the bottom are
      labelled in the panel (see NOTE).
  (c) Normalized semivariance of PERCENTILE RANKS vs pair distance,
      gamma(h) = 0.5 * E[(r_i - r_j)^2 | d_ij in bin] / Var(r), for PI and the
      six process objectives, computed EXACTLY over all n(n-1)/2 = 92.8 million
      trial pairs (no pair sampling; distance bins in DIST_EDGES, each point at
      its bin's mean pair distance). gamma = 1 means a pair at that distance
      is no more alike than two random trials; small gamma at short h = a
      smooth landscape. Ranks make the seven objectives commensurate and PI's
      heavy loss tail harmless. PI is emphasized in the accent colour,
      isobutanol yield in a strong second style, the other five thin
      (isobutanol family solid, ethanol family dashed). The panel shows the
      bins with mean h <= VARIOGRAM_HMAX = 1.6 (49 % of all pairs; every bin
      has >= 700 pairs); ALL bins are printed. Beyond h ~ 1.6 the process
      objectives' gamma rises ABOVE 1 (1.7-3.4 at h 2.1-2.6: distant pairs are
      systematically one high / one low -- the signature of a global trend,
      i.e. a hill) while PI's levels off at ~1.1-1.3; showing it would need a
      y-axis to ~3.5 that flattens the short-range comparison, and it only
      strengthens the figure's message.

CAVEAT (fairness of (a) vs (b)). The pool is not a uniform sample: every GP
campaign exploited around its own optimum, and the IBO-yield campaign
contributed most of the trials near #841 / #842 (and 124 of the 166 PI > 0
trials). The dense neighbourhoods of (a) and (b) are therefore the SAME
cluster, chosen by the IBO-yield GP. Panel (c) answers this: it scores every
objective on the SAME set of pairs, so the sampling design is common to all
seven curves, and PI is the roughest at every distance. It does not remove
the bias entirely: six of the seven campaigns exploited around a PROCESS
objective's optimum, so short-distance pairs are over-represented on
process-objective plateaus. Also, near-zero yields (e.g. isobutanol yield
1e-26 ... 1e-4 in ethanol-only designs) are ranked like any other values,
which adds rank noise to the process objectives -- a bias AGAINST the
figure's message (visible as isobutanol yield's first bin, 0.066 on only 700
pairs).

NOTE on (b): 40 of the 166 trials with PI > 0 are ETHANOL-ONLY designs
(isobutanol yield < 0.013, all with PI <= 0.05; 38 from the PI campaign) and
sit at the bottom of the isobutanol-yield landscape; every trial with
PI > 0.05 is a co-production design (isobutanol yield >= 0.11 g/g; the 20
trials with PI > 0.3 have 0.15-0.37 g/g, 18 of them >= 0.33).

Sim-safe: imports only numpy / pandas / scipy / matplotlib and the data layer
BY FILE PATH; never imports the biorefineries package, never load()s -- safe
alongside any running simulation on any numba-cache state. Run:

    python plots/plot_objective_landscape_fitness_distance.py
        [--out-dir DIR] [--stem STEM] [--dpi 300] [--seed N]
        [--pi-scale {log-tail,linear}]

Writes <out-dir>/<stem>.png (at --dpi) and .pdf, and prints the numbers the
figure shows (FDCs, variogram table). The variogram is exact, so --seed only
seeds the random draw ORDER of the grey cloud (which trial sits on top where
points overlap); the default is fixed so re-runs are reproducible.
"""
import os
import argparse
import importlib.util

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.text
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D, TICKDOWN, TICKLEFT
from matplotlib.ticker import FixedLocator, FixedFormatter, MultipleLocator
from matplotlib.transforms import ScaledTranslation
from scipy.spatial.distance import cdist
from scipy.stats import rankdata, spearmanr

_HERE = os.path.dirname(os.path.abspath(__file__))


def _load_data_layer():
    path = os.path.join(_HERE, '_objective_landscape_data.py')
    spec = importlib.util.spec_from_file_location('_objective_landscape_data',
                                                  path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


old = _load_data_layer()

# -----------------------------------------------------------------------------
# Typeface, font sizes, tick style (from
# across_feed_strat_multipanel_global_with_optima_summary.py)
# -----------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'panel_label': 12, 'legend': 9,
         'annotation': 9, 'fdc': 10}
TICK_LEN = {'major': 4.0, 'minor': 2.0}  # pt; left/bottom ticks extend this far in AND out


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
# Shared encoding (identical across the three objective-landscape figures)
# -----------------------------------------------------------------------------
GREY = '#B8B8B8'
ACCENT = '#C0392B'
INK = '#222222'
MUTED_INK = '#555555'
REF_LINE = dict(color='#8A8A8A', linewidth=0.6, linestyle=':', zorder=2)
ALL_KW = dict(s=4, c=GREY, alpha=0.5, linewidths=0, rasterized=True,
              zorder=1)
POS_KW = dict(s=9, c=ACCENT, alpha=0.9, linewidths=0, rasterized=True,
              zorder=3)
STAR_KW = dict(marker='*', markersize=13, markerfacecolor=ACCENT,
               markeredgecolor='black', markeredgewidth=0.8, linestyle='none',
               zorder=6)
DIAMOND_KW = dict(marker='D', markersize=6.5, markerfacecolor='none',
                  markeredgecolor='black', markeredgewidth=1.2,
                  linestyle='none', zorder=5)

# Variogram series: fixed categorical order (checked with the dataviz palette
# validator: CVD and normal-vision separation pass; the three light hues are
# below 3:1 contrast, relieved by the legend + the console table). PI in the
# accent colour, isobutanol yield in a strong second style, the other five
# thin; isobutanol family solid, ethanol family dashed.
VARIOGRAM_SERIES = (
    # (column, legend label, colour, linewidth, linestyle, marker, zorder)
    ('PI',                'Profitability index',     ACCENT,    1.8, '-',  'o', 6),
    ('IBO yield',         'Isobutanol yield',        '#2a78d6', 1.8, '-',  's', 5),
    ('IBO titer',         'Isobutanol titer',        '#1baf7a', 0.9, '-',  None, 3),
    ('IBO productivity',  'Isobutanol productivity', '#4a3aa7', 0.9, '-',  None, 3),
    ('EtOH yield',        'Ethanol yield',           '#eda100', 0.9, '--', None, 3),
    ('EtOH titer',        'Ethanol titer',           '#008300', 0.9, '--', None, 3),
    ('EtOH productivity', 'Ethanol productivity',    '#e87ba4', 0.9, '--', None, 3),
)
PROCESS_COLUMNS = tuple(s[0] for s in VARIOGRAM_SERIES[1:])

# Distance bins of the variogram (unit-cube Euclidean; the pool's largest
# pair distance is ~3.3), the drawn window and the pair-count floor.
DIST_EDGES = (0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2,
              1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, np.inf)
VARIOGRAM_HMAX = 1.6
MIN_PAIRS = 200
# Coarse bins of the pre-computed facts (console cross-check only).
CHECK_EDGES = (0.1, 0.2, 0.4, 0.8, 1.2, np.inf)

# PI axis ticks (raw PI values; placed on whichever scale is used).
PI_MAJOR = (0.5, 0.0, -0.5, -1.0, -2.0, -5.0)
PI_MINOR = (0.25, -0.25, -0.75, -1.5, -3.0, -4.0)
PI_MAJOR_LINEAR = (0.0, -1.0, -2.0, -3.0, -4.0, -5.0)  # 0.5 would collide with 0

ETHANOL_ONLY_MAX_IBO_YIELD = 0.1  # g/g; PI > 0 designs below it = ethanol-only


# -----------------------------------------------------------------------------
# Computations
# -----------------------------------------------------------------------------
def pi_log_tail(y):
    y = np.asarray(y, float)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(y >= 0, y, -np.log1p(-np.minimum(y, 0.0)))


def pi_log_tail_inverse(t):
    t = np.asarray(t, float)
    with np.errstate(over='ignore'):
        return np.where(t >= 0, t, -np.expm1(-np.minimum(t, 0.0)))


def fdc(values, dist):
    """Fitness-distance correlation: Spearman rho of objective vs distance."""
    ok = np.isfinite(values) & np.isfinite(dist)
    return float(spearmanr(values[ok], dist[ok])[0])


def rank_variogram(U, columns, edges=DIST_EDGES, block=384):
    """Exact all-pairs normalized semivariance of percentile ranks.

    columns -- dict name -> (n,) values. Returns (counts, mean_dist, gamma,
    overall): per-bin pair counts and mean pair distances, gamma = dict name
    -> (n_bins,) 0.5 * mean squared rank difference of the pairs in a bin /
    Var(rank), and the all-pairs value (~1 by construction; a sanity check)."""
    n = len(U)
    names = list(columns)
    R = np.column_stack([rankdata(columns[k]) / n for k in names])
    var = R.var(axis=0)
    edges = np.asarray(edges, float)
    nb = len(edges) - 1
    counts = np.zeros(nb)
    dsum = np.zeros(nb)
    ssum = np.zeros((len(names), nb))
    for s in range(0, n - 1, block):
        e = min(n, s + block)
        D = cdist(U[s:e], U[s:])
        mask = np.arange(n - s)[None, :] > np.arange(e - s)[:, None]  # j > i
        d = D[mask]
        b = np.clip(np.searchsorted(edges, d, side='right') - 1, 0, nb - 1)
        counts += np.bincount(b, minlength=nb)
        dsum += np.bincount(b, weights=d, minlength=nb)
        for q in range(len(names)):
            diff = (R[s:e, q][:, None] - R[s:, q][None, :])[mask]
            ssum[q] += np.bincount(b, weights=diff * diff, minlength=nb)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_dist = dsum / counts
        gamma = {k: 0.5 * ssum[q] / counts / var[q]
                 for q, k in enumerate(names)}
        overall = {k: 0.5 * ssum[q].sum() / counts.sum() / var[q]
                   for q, k in enumerate(names)}
    return counts, mean_dist, gamma, overall


def regroup_variogram(counts, gamma, edges, coarse_edges):
    """Pair-weighted re-aggregation of the fine variogram onto coarse bins."""
    edges = np.asarray(edges, float)
    out_counts, out = [], {k: [] for k in gamma}
    for lo, hi in zip(coarse_edges[:-1], coarse_edges[1:]):
        sel = (edges[:-1] >= lo - 1e-12) & (edges[1:] <= hi + 1e-12)
        c = counts[sel]
        out_counts.append(c.sum())
        for k, g in gamma.items():
            out[k].append(np.nansum(g[sel] * c) / c.sum() if c.sum() else np.nan)
    return np.array(out_counts), {k: np.array(v) for k, v in out.items()}


# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
def _minus(text):
    return text.replace('-', '−')


def _panel_label(fig, ax, letter):
    # left-aligned with the panel's y-axis title, just above the axes
    renderer = fig.canvas.get_renderer()
    x0 = ax.get_tightbbox(renderer).x0 / fig.bbox.width
    fig.text(x0, ax.get_position().y1 + 5 / 72 / fig.get_figheight(), letter,
             fontsize=FONTS['panel_label'], fontweight='bold', color='black',
             ha='left', va='bottom')


def _trials_under_text(fig, ax, text, x, y):
    """Number of plotted trials whose marker centre lies under `text`'s
    bounding box (a legibility check; the texts have no background, so they
    never hide a trial)."""
    # Text's own extent (an Annotation's would include its leader line)
    bb = matplotlib.text.Text.get_window_extent(text,
                                                fig.canvas.get_renderer())
    xy = ax.transData.transform(np.column_stack([x, y]))
    return int(((xy[:, 0] >= bb.x0) & (xy[:, 0] <= bb.x1)
                & (xy[:, 1] >= bb.y0) & (xy[:, 1] <= bb.y1)).sum())


def _fdc_title(fig, ax, rho):
    # above the axes, right-aligned, so it never hides a trial
    ax.text(1, 1, _minus(r'FDC $\rho$' + f' = {rho:.2f}'),
            transform=ax.transAxes + ScaledTranslation(
                0, 4 / 72, fig.dpi_scale_trans),
            ha='right', va='bottom', fontsize=FONTS['fdc'], color=INK)


def make_figure(data, args):
    f = data.frame
    U = data.U
    n = len(f)
    PI = f['PI'].to_numpy(float)
    IY = f['IBO yield'].to_numpy(float)
    a = data.anchor
    b = data.best_of('IBO yield')
    dA = data.distance_from(a)
    dB = data.distance_from(b)
    pos = PI > 0
    eth_only = pos & (IY < ETHANOL_ONLY_MAX_IBO_YIELD)
    rng = np.random.default_rng(args.seed)
    order = rng.permutation(n)  # draw order of the grey cloud only

    # -- numbers ------------------------------------------------------------
    fdcs, bests = {}, {}
    for col in ('PI',) + PROCESS_COLUMNS:
        i = data.best_of(col)
        bests[col] = i
        fdcs[col] = fdc(f[col].to_numpy(float), data.distance_from(i))
    columns = {col: f[col].to_numpy(float) for col in ('PI',) + PROCESS_COLUMNS}
    counts, mean_dist, gamma, overall = rank_variogram(U, columns)
    drawn = (counts >= MIN_PAIRS) & (mean_dist <= VARIOGRAM_HMAX)

    # -- layout -------------------------------------------------------------
    apply_font_rcparams()
    W, H = 7.3, 4.0
    fig = plt.figure(figsize=(W, H))
    gs = fig.add_gridspec(1, 3, width_ratios=(1.0, 1.0, 0.8), left=0.115,
                          right=0.985, top=1 - 0.28 / H, bottom=1.52 / H,
                          wspace=0.52)
    axA, axB, axC = (fig.add_subplot(gs[0, j]) for j in range(3))

    # -- (a) PI vs distance from #842 ----------------------------------------
    logtail = args.pi_scale == 'log-tail'
    if logtail:
        axA.set_yscale('function', functions=(pi_log_tail, pi_log_tail_inverse))
        major = PI_MAJOR
        axA.yaxis.set_minor_locator(FixedLocator(PI_MINOR))
        axA.set_ylim(-5.9, 0.62)
    else:
        major = PI_MAJOR_LINEAR
        axA.yaxis.set_minor_locator(MultipleLocator(0.5))
        axA.set_ylim(PI.min() - 0.25, 0.62)
    axA.yaxis.set_major_locator(FixedLocator(major))
    axA.yaxis.set_major_formatter(FixedFormatter(
        [_minus(f'{v:g}') for v in major]))
    axA.axhline(0.0, **REF_LINE)
    axA.scatter(dA[order], PI[order], **ALL_KW)
    axA.scatter(dA[pos], PI[pos], **POS_KW)
    axA.plot([0.0], [PI[a]], **STAR_KW)
    # the nonlinear axis is flagged in the axis title itself (no free corner
    # of the panel is wide enough for an in-panel note without covering
    # trials)
    axA.set_ylabel('Profitability index' + ('\n(log scale below 0)'
                                            if logtail else ''),
                   fontsize=FONTS['axis_title'], labelpad=2, linespacing=1.1)

    # -- (b) IBO yield vs distance from #841 --------------------------------
    axB.scatter(dB[order], IY[order], **ALL_KW)
    axB.scatter(dB[pos], IY[pos], **POS_KW)
    axB.plot([0.0], [IY[b]], **DIAMOND_KW)
    axB.plot([dB[a]], [IY[a]], **STAR_KW)
    # label the ethanol-only PI > 0 designs lying at the bottom
    xe = float(np.median(dB[eth_only]))
    eth_note = axB.annotate(
                 f'Ethanol-only\ndesigns with\nPI > 0\n(n = {int(eth_only.sum())})',
                 xy=(xe, 0.004), xytext=(3.07, 0.25), textcoords='data',
                 ha='right', va='center', multialignment='left',
                 fontsize=FONTS['annotation'],
                 color=INK, linespacing=1.1, zorder=7,
                 arrowprops=dict(arrowstyle='-', color=MUTED_INK,
                                 linewidth=0.7, shrinkA=2, shrinkB=1))
    axB.set_ylim(-0.012, 0.41)
    axB.yaxis.set_major_locator(MultipleLocator(0.1))
    axB.yaxis.set_minor_locator(MultipleLocator(0.05))
    axB.set_ylabel(r'Isobutanol yield [$\mathrm{g·g}^{-1}$]',
                   fontsize=FONTS['axis_title'], labelpad=2)

    for ax, best in ((axA, a), (axB, b)):
        ax.set_xlim(-0.1, 3.2)
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.set_xlabel(f'Distance from #{int(f["trial_number"][best])}',
                      fontsize=FONTS['axis_title'], labelpad=2)
    _fdc_title(fig, axA, fdcs['PI'])
    _fdc_title(fig, axB, fdcs['IBO yield'])

    # -- (c) rank variogram --------------------------------------------------
    axC.axhline(1.0, **REF_LINE)
    axC.text(0.04, 1.0, 'uncorrelated', transform=axC.get_yaxis_transform(),
             ha='left', va='bottom', fontsize=FONTS['annotation'],
             color=MUTED_INK)
    for col, label, color, lw, ls, marker, z in VARIOGRAM_SERIES:
        axC.plot(mean_dist[drawn], gamma[col][drawn], color=color,
                 linewidth=lw, linestyle=ls, marker=marker,
                 markersize=3.2 if marker else 0,
                 markeredgewidth=0, zorder=z, label=label)
    axC.set_xlim(0, VARIOGRAM_HMAX)
    axC.set_ylim(0, 1.12)
    axC.xaxis.set_major_locator(MultipleLocator(0.5))
    axC.xaxis.set_minor_locator(MultipleLocator(0.1))
    axC.yaxis.set_major_locator(MultipleLocator(0.5))
    axC.yaxis.set_minor_locator(MultipleLocator(0.1))
    axC.set_xlabel('Pair distance', fontsize=FONTS['axis_title'], labelpad=2)
    axC.set_ylabel('Rank semivariance', fontsize=FONTS['axis_title'],
                   labelpad=2)

    for ax in (axA, axB, axC):
        for sp in ax.spines.values():
            sp.set_linewidth(0.8)
    fig.canvas.draw()
    for ax in (axA, axB, axC):
        style_ticks(ax)
    for ax, letter in ((axA, 'a'), (axB, 'b'), (axC, 'c')):
        _panel_label(fig, ax, letter)
    under = {'(b) ethanol-only label': _trials_under_text(fig, axB, eth_note,
                                                          dB, IY)}

    # -- legends (below the panels) ----------------------------------------
    marker_handles = [
        Line2D([], [], marker='o', linestyle='none', markersize=3.2,
               markerfacecolor=GREY, markeredgewidth=0, alpha=0.8),
        Line2D([], [], marker='o', linestyle='none', markersize=4.2,
               markerfacecolor=ACCENT, markeredgewidth=0),
        Line2D([], [], **{**STAR_KW, 'markersize': 11}),
        Line2D([], [], **DIAMOND_KW),
    ]
    marker_labels = [
        f'All trials (n = {n:,})',
        f'PI > 0 (n = {int(pos.sum())})',
        f'Most profitable trial (#{int(f["trial_number"][a])}, '
        f'PI {PI[a]:.2f})',
        f'Highest isobutanol yield (#{int(f["trial_number"][b])})',
    ]
    leg1 = fig.legend(marker_handles, marker_labels, loc='upper center',
                      bbox_to_anchor=(0.5, 0.93 / H), ncol=2,
                      fontsize=FONTS['legend'], frameon=False,
                      handlelength=1.6, handletextpad=0.5, columnspacing=2.0,
                      borderaxespad=0.0)
    line_handles, line_labels = axC.get_legend_handles_labels()
    # row-major order for a 4-column legend (matplotlib fills column-major)
    ncol = 4
    idx = [i for c in range(ncol) for i in range(c, len(line_handles), ncol)]
    leg2 = fig.legend([line_handles[i] for i in idx],
                      [line_labels[i] for i in idx], loc='upper center',
                      bbox_to_anchor=(0.5, 0.50 / H), ncol=ncol,
                      fontsize=FONTS['legend'], frameon=False,
                      handlelength=2.2, handletextpad=0.5, columnspacing=1.2,
                      borderaxespad=0.0)
    for leg in (leg1, leg2):
        for t in leg.get_texts():
            t.set_color(INK)

    summary = dict(n=n, anchor=a, best_ibo=b, fdcs=fdcs, bests=bests,
                   counts=counts, mean_dist=mean_dist, gamma=gamma,
                   overall=overall, drawn=drawn, pos=pos, eth_only=eth_only,
                   PI=PI, IY=IY, dA=dA, dB=dB, under=under)
    return fig, summary


def print_summary(data, s, args):
    f = data.frame
    a, b = s['anchor'], s['best_ibo']
    PI, IY, pos, eo = s['PI'], s['IY'], s['pos'], s['eth_only']
    print(f'{data.n_raw} COMPLETE rows -> {s["n"]} unique trials '
          f'(7 campaigns pooled)')
    print(f'most profitable trial: {data.describe(a)}  PI {PI[a]:.4f}  '
          f'IBO yield {IY[a]:.4f}')
    print(f'highest IBO yield:     {data.describe(b)}  IBO yield {IY[b]:.4f}  '
          f'PI {PI[b]:.4f}  (distance from #842 {s["dA"][b]:.3f})')
    print(f'PI > 0: {int(pos.sum())} trials '
          f'({f.loc[pos, "campaign"].value_counts().to_dict()}); '
          f'PI > 0.3: {int((PI > 0.3).sum())}')
    print(f'  ethanol-only (IBO yield < {ETHANOL_ONLY_MAX_IBO_YIELD}): '
          f'{int(eo.sum())} ({f.loc[eo, "campaign"].value_counts().to_dict()}),'
          f' max PI {PI[eo].max():.3f}, max IBO yield {IY[eo].max():.4f}; '
          f'min IBO yield at PI > 0.05: {IY[PI > 0.05].min():.3f}')
    print(f'PI range {PI.min():.3f} .. {PI.max():.3f} (median '
          f'{np.median(PI):.3f}); panel (a) y-scale: {args.pi_scale}'
          + (f' (maps the range to {pi_log_tail(PI.min()):.3f} .. '
             f'{PI.max():.3f})' if args.pi_scale == 'log-tail' else '')
          + '; trials clipped or dropped: 0')
    print('trials under in-panel text (legibility check): ' + ', '.join(
        f'{k} {v}' for k, v in s['under'].items()))
    print('\nFitness-distance correlation (Spearman rho of objective vs '
          'unit-cube distance from its own best trial, all trials):')
    for col, rho in s['fdcs'].items():
        i = s['bests'][col]
        print(f'  {col:<18s} rho {rho:+.3f}   best {data.describe(i)} '
              f'= {f[col][i]:.4g}')

    cols = list(s['gamma'])
    short = {'PI': 'PI', 'IBO yield': 'IBO_y', 'IBO titer': 'IBO_t',
             'IBO productivity': 'IBO_p', 'EtOH yield': 'EtOH_y',
             'EtOH titer': 'EtOH_t', 'EtOH productivity': 'EtOH_p'}
    print('\nRank variogram, exact over all '
          f'{int(s["counts"].sum()):,} pairs (normalized semivariance of '
          f'percentile ranks; drawn = >= {MIN_PAIRS} pairs and mean h <= '
          f'{VARIOGRAM_HMAX}):')
    print(f'  {"bin":>11s} {"pairs":>10s} {"mean h":>7s} ' + ' '.join(
        f'{short[c]:>7s}' for c in cols) + '  drawn')
    edges = DIST_EDGES
    for k in range(len(edges) - 1):
        hi = 'inf' if np.isinf(edges[k + 1]) else f'{edges[k + 1]:.2f}'
        print(f'  {edges[k]:5.2f}-{hi:>5s} {int(s["counts"][k]):>10,d} '
              f'{s["mean_dist"][k]:7.3f} ' + ' '.join(
                  f'{s["gamma"][c][k]:7.3f}' for c in cols)
              + ('  yes' if s['drawn'][k] else '  no'))
    print(f'  {"all pairs":<30s}' + ' '.join(
        f'{s["overall"][c]:7.3f}' for c in cols) + '  (sanity: ~1)')
    print(f'  drawn bins hold {int(s["counts"][s["drawn"]].sum()):,} pairs '
          f'({s["counts"][s["drawn"]].sum() / s["counts"].sum():.0%})')
    cc, cg = regroup_variogram(s['counts'], s['gamma'], DIST_EDGES,
                               CHECK_EDGES)
    print('\nSame, on the coarse bins of the pre-computed facts:')
    for k in range(len(CHECK_EDGES) - 1):
        hi = ('inf' if np.isinf(CHECK_EDGES[k + 1])
              else f'{CHECK_EDGES[k + 1]:.1f}')
        proc = [cg[c][k] for c in cols[1:]]
        print(f'  {CHECK_EDGES[k]:.1f}-{hi:>3s} pairs {int(cc[k]):>10,d}  '
              f'PI {cg["PI"][k]:.3f}  process {min(proc):.3f}-'
              f'{max(proc):.3f}')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--out-dir', default=os.path.join(
        old.RESULTS_DIR, 'publication', 'Objective-landscape'))
    p.add_argument('--stem', default='objective_landscape_fitness_distance')
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--seed', type=int, default=20260924,
                   help='seed of the grey cloud draw order (the variogram is '
                        'exact and needs no seed)')
    p.add_argument('--pi-scale', choices=('log-tail', 'linear'),
                   default='log-tail',
                   help='y-scale of panel (a) (default: the log-tail scale, '
                        'ticks in raw PI)')
    args = p.parse_args(argv)

    data = old.load_landscape()
    fig, summary = make_figure(data, args)
    print_summary(data, summary, args)
    os.makedirs(args.out_dir, exist_ok=True)
    for ext in ('png', 'pdf'):
        path = os.path.join(args.out_dir, f'{args.stem}.{ext}')
        fig.savefig(path, dpi=args.dpi)
        print(f'wrote {path}')
    plt.close(fig)


if __name__ == '__main__':
    main()

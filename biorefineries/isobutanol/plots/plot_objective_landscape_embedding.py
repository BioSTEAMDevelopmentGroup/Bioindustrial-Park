#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Objective landscape as a 2-D embedding map: profitability index vs
isobutanol yield over the metabolic_split_12d decision space.

Figure. One t-SNE embedding of every unique COMPLETE trial of the seven
2026-09-23 seed-350 `ethanol_isobutanol` x `metabolic_split_12d` GP campaigns
(13,622 decision vectors), drawn twice with identical coordinates:

  (a) coloured by the profitability index PI = NPV(15 % hurdle) / TCI -- a
      rough, salt-and-pepper field: the highest values are isolated dark
      points among much less profitable neighbours (a spiky, "needle"
      landscape);
  (b) coloured by isobutanol yield -- coherent regions and a smooth gradient,
      with the trials of PI > 0 outlined by rings; the co-production ones sit
      on its high-yield arm.

The most profitable trial (ibo_yield campaign #842, PI 0.50) is starred in
both panels. Points are drawn in ascending order of the coloured value in
BOTH panels, so a high value is never hidden under a low one. The PI colour
scale runs linearly from PI = -1 (extend arrow) to the pool maximum: raw PI
reaches about -5.6 and 60 % of the trials lie below -1, and a floor at -1 --
NPV = -TCI, a design that loses more than its entire capital investment at
the 15 % hurdle -- is an economic landmark rather than a tuned quantile. It
leaves the upper colour range to the designs that are near or above
break-even, so the high tail is visible (at a -2 floor PI 0 and PI -0.3 were
nearly the same shade). The isobutanol-yield scale is linear from 0 to the
pool maximum, unclipped. Caveat stated in the figure: 40 of the 166 PI > 0
trials are ethanol-only designs (isobutanol yield <= 0.012 g/g, PI <= 0.05;
38 from the PI campaign, 2 from the ethanol-yield campaign), so those rings
sit on the zero-yield region, not the high-yield one; the 20 trials with
PI > 0.3 all lie at isobutanol yield 0.15-0.37 (median 0.35).

Data. `plots/_objective_landscape_data.py` (loaded BY FILE PATH): the pooled
trajectory CSVs and each trial's decision vector in the campaigns' own unit
cube (log axes log-normalized, linear axes min-max), where Euclidean distance
is the distance the GP itself sees. The embedding is fitted on those unit-cube
coordinates and cached in `<out-dir>/<stem>_tsne_<key>.npz`, keyed by a hash of
the coordinates, the t-SNE settings and the sklearn version (`--recompute`
forces a fresh fit).

Embedding caveat. t-SNE preserves NEIGHBOURHOODS, not distances: map distances
and cluster sizes are not metric, and the axes have no units. How faithfully
the neighbourhoods survive is reported as the trustworthiness T(k) (exact, on
the full set; 1 = every map neighbour is a 12-D neighbour). The smoothness
statistic quoted in each panel -- the Spearman correlation between each
trial's value and that of its nearest neighbour -- is computed in 12-D, NOT
on the map; the map only pictures it.

Smoothness caveat (printed, not drawn). The pooled isobutanol-yield rho (0.89
vs PI 0.75) is partly the split between the two-thirds of the trials that
make no isobutanol (a flat zero floor, the pale right-hand block of panel b)
and the co-production third. Within the co-production trials the 12-D rho is
PI 0.59 vs isobutanol yield 0.67; within the ethanol-only ones 0.75 vs 0.73.
The sharper needle-vs-plateau evidence is at the top: the 50 best trials of
each objective lose on average 0.115 (PI) vs 0.003 (isobutanol yield) in
percentile rank to their 5 nearest 12-D neighbours.

Sim-safe: reads CSVs only; never imports the biorefineries package, never
load()s, so it may run alongside any simulation on any numba-cache state.
"""
import os
import sys
import json
import time
import hashlib
import argparse
import importlib.util

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D, TICKLEFT
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import MultipleLocator, FuncFormatter
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr, rankdata
import sklearn
from sklearn.manifold import TSNE

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_by_path(name, filename):
    spec = importlib.util.spec_from_file_location(name,
                                                  os.path.join(HERE, filename))
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


old = _load_by_path('_objective_landscape_data', '_objective_landscape_data.py')

# -----------------------------------------------------------------------------
# Typeface, font sizes, tick style (after
# across_feed_strat_multipanel_global_with_optima_summary.py)
# -----------------------------------------------------------------------------
FONT_FAMILY = 'Arial'
FONTS = {'tick': 12, 'axis_title': 12, 'panel_label': 12, 'cbar_title': 12,
         'legend': 9, 'note': 9, 'rho': 10}
TICK_LEN = {'major': 4.0, 'minor': 2.0}  # pt; left/bottom ticks extend this far in AND out
INK = '#222222'        # text
NOTE_INK = '#444444'   # in-plot notes

# Shared encoding of the three objective-landscape figures
STAR_COLOR = '#C0392B'
STAR_SIZE = 13         # marker size, pt
RING_COLOR = '#C0392B'

# One sequential single-hue family for both panels (light -> dark); the palest
# 18 % of 'Blues' is dropped so the lowest values stay visible on white.
CMAP = ListedColormap(plt.get_cmap('Blues')(np.linspace(0.18, 1.0, 256)),
                      name='Blues_trunc')

# PI colour floor (extend arrow): PI = -1 <=> NPV = -TCI; see the module
# docstring. Floors tried on the same embedding: -2 (the high tail and the
# ethanol plateau at PI ~ -0.3 share one dark shade), -1 (adopted), -0.5
# (76 % of the trials clipped).
PI_CLIP_LOW = -1.0

# A PI > 0 trial below this isobutanol yield [g/g] is called ethanol-only in
# the figure note (the 40 such trials sit at <= 0.012; the next is at 0.114).
ETHANOL_ONLY_MAX_YIELD = 0.05

TSNE_DEFAULTS = dict(perplexity=40.0, early_exaggeration=12.0,
                     learning_rate='auto', max_iter=1000, init='pca',
                     angle=0.5, metric='euclidean')
TRUST_K = (5, 10, 30)
TRUST_K_SHOWN = 10


def apply_font_rcparams():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = [FONT_FAMILY, 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = FONT_FAMILY
    plt.rcParams['mathtext.it'] = f'{FONT_FAMILY}:italic'
    plt.rcParams['mathtext.bf'] = f'{FONT_FAMILY}:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def style_colorbar_ticks(cax):
    """Vertical colorbar: y ticks on both sides, major and minor; left in and
    out (the same length each way), right inward only; labels on the right.
    Call after fig.canvas.draw() so every tick object exists."""
    for which, L in TICK_LEN.items():
        cax.tick_params(axis='y', which=which, direction='inout', length=2*L,
                        left=True, right=True, labelleft=False,
                        labelright=True, labelsize=FONTS['tick'],
                        labelcolor=INK)
        get = 'get_major_ticks' if which == 'major' else 'get_minor_ticks'
        for tick in getattr(cax.yaxis, get)():
            tick.tick2line.set_marker(TICKLEFT)
            tick.tick2line.set_markersize(L)
    cax.tick_params(axis='x', which='both', bottom=False, top=False,
                    labelbottom=False)


# -----------------------------------------------------------------------------
# Embedding (cached) and neighbourhood statistics
# -----------------------------------------------------------------------------
def _cache_key(U, settings):
    h = hashlib.sha1()
    h.update(np.ascontiguousarray(U, dtype=np.float64).tobytes())
    h.update(json.dumps(settings, sort_keys=True).encode())
    h.update(sklearn.__version__.encode())
    return h.hexdigest()


def exact_trustworthiness(X, Y, k, chunk=256):
    """sklearn.manifold.trustworthiness(X, Y, n_neighbors=k), computed exactly
    in row chunks so the full 13,622-point set fits in memory (sklearn's
    implementation holds two n x n matrices). Rank of j among i's original
    neighbours = 1 + #points closer to i than j (self excluded)."""
    n = len(X)
    nn_Y = cKDTree(Y).query(Y, k=k + 1)[1][:, 1:]
    total = 0.0
    for s in range(0, n, chunk):
        rows = np.arange(s, min(n, s + chunk))
        D = cdist(X[rows], X)
        D[np.arange(len(rows)), rows] = np.inf
        dj = np.take_along_axis(D, nn_Y[rows], axis=1)          # (c, k)
        ranks = (D[:, None, :] < dj[:, :, None]).sum(axis=2) + 1
        excess = ranks - k
        total += excess[excess > 0].sum()
    return 1.0 - total * 2.0 / (n * k * (2.0 * n - 3.0 * k - 1.0))


def embed(U, perplexity, seed, cache_dir, stem, recompute=False):
    settings = dict(TSNE_DEFAULTS, perplexity=float(perplexity),
                    random_state=int(seed))
    key = _cache_key(U, settings)
    path = os.path.join(cache_dir, f'{stem}_tsne_{key[:10]}.npz')
    if not recompute and os.path.isfile(path):
        z = np.load(path, allow_pickle=False)
        if str(z['key']) == key:
            trust = json.loads(str(z['trust']))
            return z['Y'], dict(settings=settings, runtime_s=float(z['runtime_s']),
                                kl=float(z['kl']), n_iter=int(z['n_iter']),
                                trust={int(k): v for k, v in trust.items()},
                                cached=True, path=path)
    print(f't-SNE: fitting {len(U)} points, perplexity {perplexity}, '
          f'seed {seed} (Barnes-Hut) ...', flush=True)
    t0 = time.perf_counter()
    tsne = TSNE(n_components=2, random_state=int(seed), n_jobs=-1,
                **{k: v for k, v in settings.items() if k != 'random_state'})
    Y = tsne.fit_transform(U)
    runtime = time.perf_counter() - t0
    print(f't-SNE: done in {runtime:.1f} s (KL {tsne.kl_divergence_:.4f}, '
          f'{tsne.n_iter_} iterations)', flush=True)
    t0 = time.perf_counter()
    trust = {k: float(exact_trustworthiness(U, Y, k)) for k in TRUST_K}
    print(f'trustworthiness: computed in {time.perf_counter() - t0:.1f} s',
          flush=True)
    os.makedirs(cache_dir, exist_ok=True)
    np.savez(path, Y=Y, key=key, runtime_s=runtime,
             kl=float(tsne.kl_divergence_), n_iter=int(tsne.n_iter_),
             trust=json.dumps(trust), settings=json.dumps(settings))
    return Y, dict(settings=settings, runtime_s=runtime,
                   kl=float(tsne.kl_divergence_), n_iter=int(tsne.n_iter_),
                   trust=trust, cached=False, path=path)


def nn_rank_correlation(X, y):
    """Spearman correlation between each trial's value and that of its nearest
    neighbour in X (the decision vectors are unique, so the self-match is the
    first hit)."""
    nn = cKDTree(X).query(X, k=2)[1][:, 1]
    return float(spearmanr(y, y[nn]).statistic)


def top_isolation(X, y, n_top=50, k=5):
    """Mean percentile-rank drop from each of the n_top best trials to its k
    nearest neighbours in X (0 = the top is a plateau)."""
    pr = rankdata(y) / len(y)
    top = np.argsort(-y, kind='stable')[:n_top]
    nn = cKDTree(X).query(X[top], k=k + 1)[1][:, 1:]
    return float((pr[top][:, None] - pr[nn]).mean())


# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------
def draw_figure(data, Y, info, stats, out_dir, stem, dpi):
    apply_font_rcparams()
    f = data.frame
    a = data.anchor
    pi = f['PI'].to_numpy(float)
    ibo_col, ibo_label, ibo_unit = old.OBJECTIVES[old.PROCESS_EXAMPLE]
    ibo = f[ibo_col].to_numpy(float)
    positive = pi > 0

    # Square, shared map limits
    lo, hi = Y.min(axis=0), Y.max(axis=0)
    centre, half = (lo + hi) / 2, (hi - lo).max() / 2 * 1.04
    xlim = (centre[0] - half, centre[0] + half)
    ylim = (centre[1] - half, centre[1] + half)

    # Layout in inches
    W, H = 7.5, 4.02
    S = 2.6                     # map side
    y0 = 1.12                   # map bottom
    xa, xb = 0.40, 4.02         # map left edges
    cb_gap, cb_w = 0.08, 0.12
    fig = plt.figure(figsize=(W, H))

    def add_ax(x, y, w, h):
        return fig.add_axes([x / W, y / H, w / W, h / H])

    ax_a, ax_b = add_ax(xa, y0, S, S), add_ax(xb, y0, S, S)
    cax_a = add_ax(xa + S + cb_gap, y0, cb_w, S)
    cax_b = add_ax(xb + S + cb_gap, y0, cb_w, S)

    panels = (
        (ax_a, cax_a, pi, Normalize(PI_CLIP_LOW, pi.max(), clip=False), 'min',
         'Profitability index (NPV/TCI)', MultipleLocator(0.5),
         MultipleLocator(0.25), stats['rho_12d']['PI'], 'a'),
        (ax_b, cax_b, ibo, Normalize(0.0, ibo.max()), 'neither',
         f'{ibo_label} ({ibo_unit})', MultipleLocator(0.1),
         MultipleLocator(0.05), stats['rho_12d'][ibo_col], 'b'),
    )
    for ax, cax, v, norm, extend, title, major, minor, rho, label in panels:
        order = np.argsort(v, kind='stable')          # ascending: peaks on top
        sc = ax.scatter(Y[order, 0], Y[order, 1], c=v[order], cmap=CMAP,
                        norm=norm, s=3.0, linewidths=0, rasterized=True,
                        zorder=1)
        ax.plot(Y[a, 0], Y[a, 1], ls='none', marker='*', ms=STAR_SIZE,
                mfc=STAR_COLOR, mec='black', mew=0.6, zorder=5)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.8)
            sp.set_color(INK)
        ax.set_xlabel('t-SNE 1', fontsize=FONTS['axis_title'], color=INK,
                      labelpad=4)
        # smoothness statistic, computed in 12-D (not on the map)
        ax.text(1.0, 1.0 + 0.05 / S, f'12-D nearest-neighbour ρ = {rho:.2f}',
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=FONTS['rho'], color=INK)
        fig.text((ax.get_position().x0 - 0.30 / W), (y0 + S + 0.05) / H,
                 label, fontsize=FONTS['panel_label'], fontweight='bold',
                 color='black', ha='left', va='bottom')
        cbar = fig.colorbar(sc, cax=cax, extend=extend)
        cbar.outline.set_linewidth(0.8)
        cbar.ax.yaxis.set_major_locator(major)
        cbar.ax.yaxis.set_minor_locator(minor)
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: f'{v:g}'.replace('-', '−')))
        cbar.set_label(title, fontsize=FONTS['cbar_title'], color=INK,
                       labelpad=6)
    ax_a.set_ylabel('t-SNE 2', fontsize=FONTS['axis_title'], color=INK,
                    labelpad=4)

    # PI > 0 rings on the isobutanol-yield map
    ax_b.scatter(Y[positive, 0], Y[positive, 1], s=16, facecolors='none',
                 edgecolors=RING_COLOR, linewidths=0.6, rasterized=True,
                 zorder=3)

    handles = [
        Line2D([], [], ls='none', marker='o', ms=4, mfc='none',
               mec=RING_COLOR, mew=0.6),
        Line2D([], [], ls='none', marker='*', ms=STAR_SIZE, mfc=STAR_COLOR,
               mec='black', mew=0.6),
    ]
    labels = [f'PI > 0 (n = {int(positive.sum())})',
              f'Most profitable trial (#{int(f["trial_number"][a])}, '
              f'PI {pi[a]:.2f})']
    fig.legend(handles, labels, loc='center', ncol=2, frameon=False,
               fontsize=FONTS['legend'], bbox_to_anchor=(0.5, 0.72 / H),
               handletextpad=0.4, columnspacing=2.0)

    s = info['settings']
    ethanol_only = positive & (ibo < ETHANOL_ONLY_MAX_YIELD)
    floor = f'{PI_CLIP_LOW:g}'.replace('-', '−')
    note = (
        f't-SNE (perplexity {s["perplexity"]:g}) of {len(f):,} unique trials '
        f'from seven campaigns; trustworthiness '
        f'{info["trust"][TRUST_K_SHOWN]:.2f} (k = {TRUST_K_SHOWN}); '
        'map distances are not metric.\n'
        'ρ: rank correlation of each trial\'s value with that of its nearest '
        'neighbour in the 12-D decision space.\n'
        f'PI colours clipped at {floor} (NPV = −TCI; '
        f'{stats["pi_clipped_frac"] * 100:.0f} % of trials below). '
        f'{int(ethanol_only.sum())} of the PI > 0 trials are ethanol-only '
        f'(PI ≤ {np.ceil(pi[ethanol_only].max() * 100) / 100:g}).')
    fig.text(0.5, 0.06 / H, note, ha='center', va='bottom',
             fontsize=FONTS['note'], color=NOTE_INK, linespacing=1.35)

    fig.canvas.draw()
    style_colorbar_ticks(cax_a)
    style_colorbar_ticks(cax_b)

    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ext in ('png', 'pdf'):
        p = os.path.join(out_dir, f'{stem}.{ext}')
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    plt.close(fig)
    return paths


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--out-dir', default=os.path.join(
        old.RESULTS_DIR, 'publication', 'Objective-landscape'))
    p.add_argument('--stem', default='objective_landscape_embedding')
    p.add_argument('--dpi', type=int, default=300)
    p.add_argument('--perplexity', type=float,
                   default=TSNE_DEFAULTS['perplexity'])
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--recompute', action='store_true',
                   help='ignore a cached embedding and fit a fresh one')
    args = p.parse_args(argv)

    data = old.load_landscape()
    f, U, a = data.frame, data.U, data.anchor
    pi = f['PI'].to_numpy(float)
    print(f'{data.n_raw} COMPLETE rows -> {len(f)} unique trials '
          f'({f["campaign"].nunique()} campaigns)')
    print(f'most profitable trial: {data.describe(a)}  PI {pi[a]:.4f}  '
          f'IBO yield {f["IBO yield"][a]:.4f}')

    Y, info = embed(U, args.perplexity, args.seed, args.out_dir, args.stem,
                    recompute=args.recompute)

    columns = [c for c, _, _ in old.OBJECTIVES.values()]
    columns = ['PI'] + [c for c in columns if c != 'PI (log-tail)']
    stats = dict(
        rho_12d={c: nn_rank_correlation(U, f[c].to_numpy(float))
                 for c in columns},
        rho_map={c: nn_rank_correlation(Y, f[c].to_numpy(float))
                 for c in columns},
        drop_12d={c: top_isolation(U, f[c].to_numpy(float)) for c in columns},
        drop_map={c: top_isolation(Y, f[c].to_numpy(float)) for c in columns},
        pi_clipped_frac=float((pi < PI_CLIP_LOW).mean()),
    )

    s = info['settings']
    print(f"t-SNE: perplexity {s['perplexity']:g}, seed {s['random_state']}, "
          f"init {s['init']}, learning_rate {s['learning_rate']}, "
          f"early_exaggeration {s['early_exaggeration']:g}, "
          f"max_iter {s['max_iter']}, angle {s['angle']}; "
          f"{'cached' if info['cached'] else 'fitted'} "
          f"(fit {info['runtime_s']:.1f} s, KL {info['kl']:.4f}, "
          f"{info['n_iter']} iterations) -> {info['path']}")
    print('trustworthiness (exact, full set): ' + ', '.join(
        f'k={k}: {v:.4f}' for k, v in sorted(info['trust'].items())))
    print(f"{'metric':<20}{'NN rho 12-D':>12}{'NN rho map':>12}"
          f"{'top-50 drop 12-D':>18}{'top-50 drop map':>17}")
    for c in columns:
        print(f"{c:<20}{stats['rho_12d'][c]:>12.3f}{stats['rho_map'][c]:>12.3f}"
              f"{stats['drop_12d'][c]:>18.3f}{stats['drop_map'][c]:>17.3f}")
    print(f'PI colour scale [{PI_CLIP_LOW:g}, {pi.max():.4f}]; '
          f'{stats["pi_clipped_frac"]:.1%} of trials below the floor '
          f'(raw PI min {pi.min():.3f})')
    ibo = f[old.OBJECTIVES[old.PROCESS_EXAMPLE][0]].to_numpy(float)
    positive, top = pi > 0, pi > 0.3
    eo = positive & (ibo < ETHANOL_ONLY_MAX_YIELD)
    print(f'PI > 0: {int(positive.sum())} trials, {int(eo.sum())} of them '
          f'ethanol-only (IBO yield < {ETHANOL_ONLY_MAX_YIELD:g}; PI <= '
          f'{pi[eo].max():.3f}); PI > 0.3: {int(top.sum())} trials at IBO '
          f'yield {ibo[top].min():.3f}-{ibo[top].max():.3f} '
          f'(median {np.median(ibo[top]):.3f})')
    # Pooled IBO-yield smoothness is partly the zero / non-zero split of the
    # pool (two-thirds of the trials make no isobutanol): report it per region
    nn = cKDTree(U).query(U, k=2)[1][:, 1]
    for name, mask in (('co-production (IBO yield >= '
                        f'{ETHANOL_ONLY_MAX_YIELD:g})',
                        ibo >= ETHANOL_ONLY_MAX_YIELD),
                       (f'ethanol-only (< {ETHANOL_ONLY_MAX_YIELD:g})',
                        ibo < ETHANOL_ONLY_MAX_YIELD)):
        print(f'12-D NN rho within {name}, n = {int(mask.sum())}: '
              f'PI {spearmanr(pi[mask], pi[nn][mask]).statistic:.3f}, '
              f'IBO yield {spearmanr(ibo[mask], ibo[nn][mask]).statistic:.3f}')

    paths = draw_figure(data, Y, info, stats, args.out_dir, args.stem,
                        args.dpi)
    for path in paths:
        print('saved', path)


if __name__ == '__main__':
    main()

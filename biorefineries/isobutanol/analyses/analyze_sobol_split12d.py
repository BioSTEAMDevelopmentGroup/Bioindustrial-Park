#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Stage 2 of the Sobol' sensitivity analysis: surrogate per metric ->
closed-subset / first-order / total / Shapley indices over the FEASIBLE
campaign domain -> tables, summary sentence, two figures.

SIM-SAFE: reads <study>_trajectory.csv + <study>_design.json; loads
sensitivity_analysis / kinetic_optimization BY FILE PATH; never imports the
package, never load()s -- safe alongside a running simulation on any numba
cache state (also on a partially complete stage-1 CSV).

    python analyze_sobol_split12d.py --study-name <STUDY_NAME>
    python analyze_sobol_split12d.py --self-test        # synthetic, no data needed

--gp-metrics names the metrics whose surrogate may be a Gaussian process
(default: the two PI metrics); every other metric is fitted with gradient
boosting only, because one GP prediction of the design costs ~0.6 s and every
metric needs 2^d - 2 subsets x n_replicates of them (hours per GP metric),
while the trees predict in milliseconds at a CV Q2 within ~0.02 of the GP's.

Pre-optimization screening additions (spec
docs/superpowers/specs/2026-09-21-sobol-screening-measure-design.md): with a
threshold (--pi-threshold, default = the design record's meta.baseline_PI =
the scenario-A PI) the run also reports TARGET indices on the indicator
'PI > baseline' (what decides whether a design beats A; same base measure)
and, when at least --min-zone-rows simulated rows lie above the threshold,
CONDITIONAL indices on 'PI | PI > baseline' (what sets how much better; the
pick-freeze partners are drawn on feasible AND surrogate-zone rows).
"""
import argparse
import importlib.util
import json
import os
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_DIR = os.path.dirname(HERE)
RESULTS = os.path.join(HERE, 'results')

def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PKG_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

sa = _load('sa', 'sensitivity_analysis.py')
ko = _load('ko', 'kinetic_optimization.py')

HEADLINE = 'PI'
TARGET = 'PI > baseline'            # indicator pseudo-metric (target sensitivity)
CONDITIONAL = 'PI | PI > baseline'  # PI restricted to the zone (conditional sensitivity)
#: Never analysed: IRR is -inf at most points (variance undefined); the
#: convergence diagnostics are solver bookkeeping, not model outputs.
EXCLUDED_METRICS = ('IRR', 'spike_feed_residual', 'n_sims_run', 'final_drift')
#: Metrics a GP may be fitted for by default (--gp-metrics); every other metric
#: uses gradient boosting only -- a GP costs hours of predictions per metric.
DEFAULT_GP_METRICS = ('PI', 'PI (log-tail)')
LABELS = {'k_3': 'k$_3$ (Pdc)', 'k_6': 'k$_6$ (Adh1)', 'k_13': 'k$_{13}$ (ALS)',
          'k_17': 'k$_{17}$ (Adh6)', 'glycolysis': 'glycolysis (k$_{\\mathrm{1e}}$, k$_{\\mathrm{1l}}$, k$_{\\mathrm{1h}}$)',
          'ehrlich_downstream': 'k$_{14}$ (Ilv5), k$_{15}$ (Ilv3), k$_{16}$ (Aro10)',
          'inhib_ethanol': 'ethanol inhibition',
          'inhib_isobutanol': 'isobutanol inhibition',
          'inhib_acetate': 'acetate inhibition',
          'threshold_conc': 'thresh. sugar conc.', 'target_delta': 'target sugar conc.',
          'max_n_spikes': 'max. spikes'}

#%% Data

def load_inputs(study_name, min_rows):
    csv_path = os.path.join(RESULTS, study_name + '_trajectory.csv')
    with open(os.path.join(RESULTS, study_name + '_design.json')) as fh:
        design = json.load(fh)
    df = pd.read_csv(csv_path)
    counts = df['state'].value_counts().to_dict()
    ok = df[df['state'] == 'COMPLETE'].reset_index(drop=True)
    if len(ok) < min_rows:
        raise SystemExit(f'{len(ok)} COMPLETE rows < --min-rows {min_rows}; refusing.')
    space = design['search_space']
    U = np.array([ko.external_to_unit({n: row[n] for n in space}, space)
                  for _, row in ok.iterrows()])
    return design, ok, U, counts

def eligible_metrics(ok, requested):
    metrics, skipped = [], []
    for name in (requested or [m for m in ko.TRACKED_METRICS if m not in EXCLUDED_METRICS]):
        y = pd.to_numeric(ok[name], errors='coerce').to_numpy(dtype=float)
        if not np.isfinite(y).all():
            skipped.append((name, f'{int((~np.isfinite(y)).sum())} non-finite rows'))
        elif np.ptp(y) == 0.0:
            skipped.append((name, 'constant'))
        else:
            metrics.append(name)
    if HEADLINE not in metrics:
        raise SystemExit(f'headline metric {HEADLINE!r} is not analysable: {skipped}')
    return metrics, skipped

#%% Analysis core (shared by the real run, the convergence refits and --self-test)

def resolve_gp_metrics(requested):
    """--gp-metrics -> the collection `analyse` takes: None (the single literal
    'all') = a GP candidate for every metric; omitted = DEFAULT_GP_METRICS; a
    bare --gp-metrics (empty list) = none. Names that are not analysed are
    ignored."""
    if requested is None:
        return set(DEFAULT_GP_METRICS)
    if list(requested) == ['all']:
        return None
    return set(requested)

def analyse(U, Y, vf, *, n_base, n_replicates, seed=0, label='', gp_metrics=None):
    """Y = {metric: y}. `gp_metrics` = the metrics whose surrogate may be a GP
    (None = all of them); every other metric is fitted with gradient boosting
    only, which predicts thousands of times faster.
    Returns (surrogates, S {metric: (R, 2^d)}, fallback)."""
    surrogates = {}
    for m, y in Y.items():
        t0 = time.time()
        # The indicator pseudo-metric is always trees (a 0/1 target; never GP-eligible).
        candidates = (sa.SURROGATE_CANDIDATES
                      if m != TARGET and (gp_metrics is None or m in gp_metrics)
                      else ('hgb',))
        surrogates[m] = sa.fit_surrogates(U, y, seed=seed, candidates=candidates)
        s = surrogates[m]
        q2_text = ', '.join(f'{k} {v:.3f}' for k, v in s.q2.items())
        print(f'  {label}{m}: {s.name} (Q2 {q2_text})'
              f'{"" if s.reliable else "  ** UNRELIABLE **"}  [{time.time() - t0:.0f} s]',
              flush=True)
    t0 = time.time()
    def progress(r, mask, last):
        if mask in (64, 1024) or mask == last:
            rate = (time.time() - t0)/(r*last + mask)
            print(f'    replicate {r + 1}/{n_replicates}, subset {mask}/{last}; '
                  f'~{rate*(n_replicates*last - r*last - mask)/60:.0f} min left', flush=True)
    S, fallback = sa.all_closed_indices(
        {m: s.predict for m, s in surrogates.items()}, vf, vf.d,
        n_base=n_base, n_replicates=n_replicates, seed=seed, progress=progress)
    return surrogates, S, fallback

def conditional_analysis(U, y_pi, threshold, target_surrogate, vf, *, n_base, n_replicates,
                         seed=0, mass_draws=None):
    """Indices of PI RESTRICTED to the zone PI > threshold: the PI surrogate is
    fitted on the zone rows only (trees; CV Q2 on those rows) and the
    pick-freeze base / partner samples are drawn on feasible AND
    surrogate-zone (target_surrogate > 0.5) rows, so the closed indices are of
    Var(PI | zone). Returns (surrogate, S (R, 2^d), fallback, zone_mass) or
    None with a printed reason when the surrogate zone has no mass on
    `mass_draws` (16 n_base) feasible draws."""
    zone = y_pi > threshold
    sur = sa.fit_surrogates(U[zone], y_pi[zone], seed=seed, candidates=('hgb',))
    rng = np.random.default_rng([seed, 99])
    X = sa.sample_feasible(mass_draws or 16*n_base, vf.d, vf, rng)
    mass = float((target_surrogate.predict(X) > 0.5).mean())
    print(f'  {CONDITIONAL}: hgb (Q2 {sur.q2["hgb"]:.3f} on {int(zone.sum())} zone rows)'
          f'{"" if sur.reliable else "  ** UNRELIABLE **"}; surrogate zone mass '
          f'{mass:.2%} of the feasible measure', flush=True)
    if mass == 0.0:
        print('  conditional indices skipped: the surrogate zone has zero mass on '
              f'{len(X)} feasible draws.', flush=True)
        return None
    region = lambda Uq: vf(Uq) & (target_surrogate.predict(Uq) > 0.5)
    t0 = time.time()
    def progress(r, mask, last):
        if mask in (64, 1024) or mask == last:
            rate = (time.time() - t0)/(r*last + mask)
            print(f'    [zone] replicate {r + 1}/{n_replicates}, subset {mask}/{last}; '
                  f'~{rate*(n_replicates*last - r*last - mask)/60:.0f} min left', flush=True)
    S, fallback = sa.all_closed_indices({CONDITIONAL: sur.predict}, region, vf.d,
                                        n_base=n_base, n_replicates=n_replicates,
                                        seed=seed + 1, progress=progress)
    return sur, S[CONDITIONAL], fallback, mass

def _headline_lines(S_m, names, label):
    """'  label: top Shapley ...; smallest set >= 80 %: ...; best triple ...'"""
    d = len(names)
    mean = S_m.mean(axis=0)
    shap = sa.shapley_effects(S_m, d).mean(axis=0)
    top3 = np.argsort(shap)[::-1][:3]
    reach_mask, reach_value = sa.smallest_subset_reaching(mean, d, 0.8)
    best3 = sa.best_subsets(mean, d, sizes=(3,))[3]
    return (f'  {label}: top Shapley ' + ', '.join(f'{names[j]} {shap[j]:.0%}' for j in top3)
            + f' ({shap[top3].sum():.0%} together); smallest set >= 80 %: '
            + ', '.join(sa.mask_names(reach_mask, names)) + f' ({reach_value:.0%}); '
            f'best triple {" + ".join(sa.mask_names(best3[0], names))} = {best3[1]:.0%}.')

def index_table(names, surrogates, S):
    rows = []
    for m, s in surrogates.items():
        d = len(names)
        for kind, fn in (('S1', sa.first_order), ('ST', sa.total_order),
                         ('Shapley', sa.shapley_effects)):
            per_rep = fn(S[m], d)                              # (R, d)
            for j, name in enumerate(names):
                rows.append(dict(metric=m, parameter=name, index=kind,
                                 mean=per_rep[:, j].mean(),
                                 sd=per_rep[:, j].std(ddof=1) if len(per_rep) > 1 else np.nan,
                                 surrogate=s.name, Q2=s.q2[s.name], reliable=s.reliable))
    return pd.DataFrame(rows)

def subset_table(names, S_headline, sizes=(1, 2, 3, 4)):
    d, mean = len(names), S_headline.mean(axis=0)
    sd = S_headline.std(axis=0, ddof=1) if len(S_headline) > 1 else np.full_like(mean, np.nan)
    rows = [dict(size=bin(mask).count('1'), parameters=' + '.join(sa.mask_names(mask, names)),
                 S_closed=mean[mask], sd=sd[mask])
            for mask in range(1, 2**d) if bin(mask).count('1') in sizes]
    return (pd.DataFrame(rows).sort_values(['size', 'S_closed'], ascending=[True, False])
            .reset_index(drop=True))

#%% Figures (project style: Arial, 12/12/10/9 pt, ticks on all sides; PI bars legend inside the panel)

def _style():
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold', 'mathtext.fallback': 'stixsans',
        'axes.labelsize': 12, 'axes.titlesize': 12, 'xtick.labelsize': 10,
        'ytick.labelsize': 10, 'legend.fontsize': 9})

def _ticks(ax):
    """Project tick style: ticks on all four sides, major 4 pt / minor 2 pt;
    left / bottom point in AND out, top / right point inward (half length, so
    every tick protrudes the same distance into the panel)."""
    ax.minorticks_on()
    for which, length in (('major', 4), ('minor', 2)):
        ax.tick_params(which=which, length=length, direction='inout',
                       left=True, bottom=True, top=False, right=False)
    top, right = ax.secondary_xaxis('top'), ax.secondary_yaxis('right')
    top.set_xticks(ax.get_xticks())
    right.minorticks_on()
    for extra in (top, right):
        extra.tick_params(which='major', length=2, direction='in',
                          labeltop=False, labelright=False)
        extra.tick_params(which='minor', length=1, direction='in')
    # categorical x axis: no minor ticks, bottom or top
    ax.xaxis.set_tick_params(which='minor', bottom=False)
    top.xaxis.set_tick_params(which='minor', top=False)

def plot_bars(table, names, q2_text, path):
    _style()
    t = table[table.metric == HEADLINE]
    order = (t[t['index'] == 'Shapley'].set_index('parameter')['mean']
             .sort_values(ascending=False).index.tolist())
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    width, x = 0.27, np.arange(len(order))
    ax.axhline(0.0, color='k', linewidth=0.8, zorder=0.5)
    for k, (kind, colour, label) in enumerate((
            ('S1', '#9ecae1', 'first-order'), ('Shapley', '#08519c', 'Shapley effect'),
            ('ST', '#fdae6b', 'total-order'))):
        sub = t[t['index'] == kind].set_index('parameter').loc[order]
        ax.bar(x + (k - 1)*width, sub['mean'], width, yerr=sub['sd'], capsize=2,
               color=colour, edgecolor='k', linewidth=0.5, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS.get(n, n) for n in order], rotation=40, ha='right')
    ax.set_ylim(-0.04, 0.6)
    ax.set_ylabel('Share of PI variance')
    ax.set_title(f'Sensitivity of the profitability index ({q2_text})', fontweight='bold')
    _ticks(ax)
    # inside the panel: bars are sorted by Shapley effect, so the upper right is empty
    ax.legend(loc='upper right', frameon=False)
    fig.savefig(path + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(path + '.pdf', bbox_inches='tight')
    plt.close(fig)

def plot_heatmap(table, names, path):
    _style()
    t = table[table['index'] == 'Shapley']
    metrics = [HEADLINE] + [m for m in t.metric.unique() if m != HEADLINE]
    grid = t.pivot(index='metric', columns='parameter', values='mean').loc[metrics, names]
    fig, ax = plt.subplots(figsize=(8.5, 0.42*len(metrics) + 2.4))
    im = ax.imshow(grid.to_numpy(), cmap='Blues', vmin=0.0,
                   vmax=max(0.5, float(np.nanmax(grid.to_numpy()))), aspect='auto')
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels([LABELS.get(n, n) for n in names], rotation=40, ha='right')
    unreliable = set(table[~table.reliable].metric)
    ax.set_yticks(range(len(metrics)))
    ax.set_yticklabels([m + (' *' if m in unreliable else '') for m in metrics])
    for i in range(len(metrics)):
        for j in range(len(names)):
            v = grid.iat[i, j]
            if v >= 0.05:
                ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=9,
                        color='w' if v > 0.3 else 'k')
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label('Shapley effect (share of variance)', fontsize=12)
    ax.set_title('Shapley effects across the feasible campaign space'
                 + ('  (* surrogate Q$^2$ < 0.8)' if unreliable else ''), fontweight='bold')
    fig.savefig(path + '.png', dpi=600, bbox_inches='tight')
    fig.savefig(path + '.pdf', bbox_inches='tight')
    plt.close(fig)

#%% Summary text

def write_summary(path, *, study_name, counts, n_rows, irr_finite, skipped, names,
                  surrogates, S, fallback, y_headline, gp_metrics=None, measure='campaign',
                  target=None):
    """`target` (None = no threshold): dict(threshold, source, n_zone, fraction,
    target_skip (str or None), conditional (dict(q2, reliable, mass, fallback)
    or None), conditional_skip (str or None)) -- the TARGET / CONDITIONAL
    block of the summary."""
    d, s = len(names), surrogates[HEADLINE]
    tail_threshold = target['threshold'] if target else 0.0
    mean = S[HEADLINE].mean(axis=0)
    sd = S[HEADLINE].std(axis=0, ddof=1) if len(S[HEADLINE]) > 1 else np.zeros_like(mean)
    best = sa.best_subsets(mean, d, sizes=(1, 2, 3, 4))
    reach_mask, reach_value = sa.smallest_subset_reaching(mean, d, 0.8)
    shap = sa.shapley_effects(S[HEADLINE], d).mean(axis=0)
    top3 = np.argsort(shap)[::-1][:3]
    lines = [
        f"Sobol' / Shapley sensitivity analysis -- campaign {study_name}",
        f'Rows: {n_rows} COMPLETE of {sum(counts.values())} ({counts}); '
        f'IRR finite at {irr_finite:.1%} of COMPLETE rows (not analysed).',
        f'Input measure: {measure}. Headline output: raw PI; surrogate {s.name}, '
        f'CV Q2 {s.q2[s.name]:.3f}'
        + ('' if s.reliable else '  ** BELOW 0.8: treat every number below as unreliable **'),
        f'Share of Var(PI) carried by the PI < {tail_threshold:.4g} tail: '
        f'{sa.tail_variance_share(y_headline, tail_threshold):.1%} (large = the indices mostly '
        'describe what separates working from failing designs).',
        '',
        'HEADLINE',
        f'  The three parameters {", ".join(sa.mask_names(best[3][0], names))} together '
        f'explain {best[3][1]:.0%} (+/- {sd[best[3][0]]:.0%}) of the variability in PI '
        '(closed Sobol\' index: variance explained by knowing them, incl. their interactions).',
        f'  Smallest set reaching 80 %: {", ".join(sa.mask_names(reach_mask, names))} '
        f'({reach_value:.0%}).',
        '  Best subset per size: ' + '; '.join(
            f'{k}: {" + ".join(sa.mask_names(m, names))} = {v:.0%}' for k, (m, v) in best.items()),
        '  Top Shapley effects (sum to 100 % over all 12): ' + ', '.join(
            f'{names[j]} {shap[j]:.0%}' for j in top3)
        + f' = {shap[top3].sum():.0%} together.',
    ]
    if target is not None:
        lines += [
            '',
            'TARGET / CONDITIONAL (pre-optimization screening against the baseline)',
            f'  Threshold: PI > {target["threshold"]:.4f} ({target["source"]}); zone = '
            f'{target["n_zone"]} of {n_rows} rows ({target["fraction"]:.1%}).',
        ]
        if TARGET in S:
            st = surrogates[TARGET]
            lines.append(_headline_lines(
                S[TARGET], names,
                f'Entry (indicator PI > threshold, same base measure; trees Q2 '
                f'{st.q2[st.name]:.3f}{"" if st.reliable else " **UNRELIABLE**"})'))
        else:
            lines.append(f'  Entry indicator not analysed: {target["target_skip"]}')
        if target.get('conditional'):
            c = target['conditional']
            lines.append(_headline_lines(
                S[CONDITIONAL], names,
                f'Within the zone (PI | PI > threshold; trees fitted on the zone rows, Q2 '
                f'{c["q2"]:.3f}{"" if c["reliable"] else " **UNRELIABLE**"}; surrogate zone '
                f'mass {c["mass"]:.2%} of the feasible measure; partner fallback max '
                f'{c["fallback"]:.2%})'))
        else:
            lines.append(f'  Conditional indices skipped: {target["conditional_skip"]}')
        lines.append('  * Entry = what decides whether a design beats the baseline at all; '
                     'within-zone = what sets how much better it is. A near-zero first-order '
                     'index with a large total index means the variable acts only jointly '
                     '(the conjunctive Ehrlich pathway).')
    lines += [
        '',
        'HOW TO READ THIS',
        '  * Domain: the FEASIBLE part of the campaign box (enzyme burden + volume bound) under '
        + ("the campaign's log-uniform measure -- variability across the designs the optimizer "
           'was allowed to propose, not around any one design.' if measure == 'campaign' else
           'the SCREENING measure (an engineering prior, not data: k_13 / k_17 / '
           'ehrlich_downstream linear-uniform from zero, k_3 / k_6 / glycolysis log-uniform on '
           '0.1x-4x of scenario A) -- variability across the designs a pre-optimization screen '
           'considers, not around any one design.'),
        '  * The feasibility constraint makes the inputs dependent: a closed index includes '
        'influence mediated through the constraint; Shapley effects are the attribution that '
        'sums to 100 %; first-order / total indices are reported for reference and need not '
        'bracket the Shapley effect.',
        '  * Indices are of the SURROGATE; Q2 bounds how much of the true variance it carries. '
        'The GP Q2 carries a mild hyper-parameter leak (kernel tuned once on a subsample).',
        '  * GP-eligible metrics (GP vs gradient boosting by CV Q2): '
        + ('all of them.' if gp_metrics is None else
           (', '.join(m for m in surrogates if m in gp_metrics) or 'none')
           + '; every other metric used gradient boosting only (a GP predicts thousands of '
             'times more slowly, and each subset needs one prediction).'),
        '  * max_n_spikes is an integer treated as numeric.',
        f'  * Conditional-sampling fallback (partner kept its own coordinates): max '
        f'{fallback.max():.2%} of base points over all subsets.',
    ]
    if skipped:
        lines += ['', 'Skipped metrics: ' + '; '.join(f'{m} ({why})' for m, why in skipped)]
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))

#%% Entry points

def run(args):
    design, ok, U, counts = load_inputs(args.study_name, args.min_rows)
    names = list(design['search_space'])
    vf = sa.VectorizedFeasibility.from_design(design)
    # Cross-checks of the model-free domain against what stage 1 recorded:
    assert vf(U).all(), f'{int((~vf(U)).sum())} simulated rows are infeasible under VectorizedFeasibility'
    assert np.allclose(vf.phi_M(U), ok['Phi_M'].to_numpy(dtype=float), rtol=1e-6), \
        'VectorizedFeasibility.phi_M disagrees with the recorded Phi_M column'
    metrics, skipped = eligible_metrics(ok, args.metrics)
    gp_metrics = resolve_gp_metrics(args.gp_metrics)
    irr = pd.to_numeric(ok['IRR'], errors='coerce').to_numpy(dtype=float)
    Y = {m: pd.to_numeric(ok[m]).to_numpy(dtype=float) for m in metrics}
    print(f'{len(ok)} COMPLETE rows, {len(metrics)} metrics; d = {vf.d}.')
    # Target / conditional set-up: threshold from the command line or the record.
    meta = design.get('meta', {})
    target = None
    if args.pi_threshold is not None:
        target = dict(threshold=float(args.pi_threshold), source='--pi-threshold')
    elif meta.get('baseline_PI') is not None:
        target = dict(threshold=float(meta['baseline_PI']),
                      source=f"design record: scenario-{meta.get('anchor', 'A')} baseline PI")
    else:
        print('No threshold (no --pi-threshold and the design record predates meta.baseline_PI): '
              'target / conditional indices skipped.', flush=True)
    if target is not None:
        zone = Y[HEADLINE] > target['threshold']
        target.update(n_zone=int(zone.sum()), fraction=float(zone.mean()),
                      target_skip=None, conditional=None, conditional_skip=None)
        if 0 < target['n_zone'] < len(zone):
            Y[TARGET] = zone.astype(float)
        else:
            target['target_skip'] = (f'{target["n_zone"]} of {len(zone)} rows above the '
                                     'threshold (constant indicator)')
        if target['n_zone'] < args.min_zone_rows:
            target['conditional_skip'] = (f'{target["n_zone"]} zone rows < --min-zone-rows '
                                          f'{args.min_zone_rows}')
        print(f'Threshold PI > {target["threshold"]:.4f} ({target["source"]}): '
              f'{target["n_zone"]} zone rows ({target["fraction"]:.1%}).', flush=True)
    surrogates, S, fallback = analyse(U, Y, vf, n_base=args.n_base,
                                      n_replicates=args.n_replicates,
                                      gp_metrics=gp_metrics)
    if target is not None and TARGET in surrogates and target['conditional_skip'] is None:
        cond = conditional_analysis(U, Y[HEADLINE], target['threshold'], surrogates[TARGET], vf,
                                    n_base=args.n_base, n_replicates=args.n_replicates)
        if cond is None:
            target['conditional_skip'] = 'the surrogate zone has zero mass on the feasible measure'
        else:
            sur_c, S_c, fb_c, mass = cond
            surrogates[CONDITIONAL], S[CONDITIONAL] = sur_c, S_c
            target['conditional'] = dict(q2=sur_c.q2['hgb'], reliable=sur_c.reliable,
                                         mass=mass, fallback=float(fb_c.max()))
    elif target is not None and target['conditional_skip'] is None:
        target['conditional_skip'] = 'no entry surrogate'
    out = os.path.join(RESULTS, args.study_name + '_sobol_')
    table = index_table(names, surrogates, S)
    table.to_csv(out + 'indices.csv', index=False)
    subset_table(names, S[HEADLINE]).to_csv(out + 'PI_subsets.csv', index=False)
    if TARGET in S:
        subset_table(names, S[TARGET]).to_csv(out + 'PI_target_subsets.csv', index=False)
    if CONDITIONAL in S:
        subset_table(names, S[CONDITIONAL]).to_csv(out + 'PI_conditional_subsets.csv', index=False)
    # Convergence of the headline with the (nested) design size
    conv = []
    for n in [n for n in args.convergence_sizes if n < len(ok)] + [len(ok)]:
        if n == len(ok):
            s_n, S_n = surrogates[HEADLINE], S[HEADLINE]
        else:
            sur_n, S_all, _ = analyse(U[:n], {HEADLINE: Y[HEADLINE][:n]}, vf,
                                      n_base=args.n_base, n_replicates=args.n_replicates,
                                      label=f'[n={n}] ', gp_metrics=gp_metrics)
            s_n, S_n = sur_n[HEADLINE], S_all[HEADLINE]
        mean = S_n.mean(axis=0)
        best3 = sa.best_subsets(mean, vf.d, sizes=(3,))[3]
        shap = sa.shapley_effects(S_n, vf.d).mean(axis=0)
        conv.append(dict(n_rows=n, surrogate=s_n.name, Q2=s_n.q2[s_n.name],
                         best_triple=' + '.join(sa.mask_names(best3[0], names)),
                         best_triple_S=best3[1],
                         **{f'Shapley_{name}': shap[j] for j, name in enumerate(names)}))
    pd.DataFrame(conv).to_csv(out + 'PI_convergence.csv', index=False)
    s = surrogates[HEADLINE]
    plot_bars(table, names, f'{s.name.upper()} surrogate, Q$^2$ = {s.q2[s.name]:.2f}', out + 'PI_bars')
    plot_heatmap(table, names, out + 'shapley_heatmap')
    write_summary(out + 'summary.txt', study_name=args.study_name, counts=counts,
                  n_rows=len(ok), irr_finite=float(np.isfinite(irr).mean()), skipped=skipped,
                  names=names, surrogates=surrogates, S=S, fallback=fallback,
                  y_headline=Y[HEADLINE], gp_metrics=gp_metrics,
                  measure=meta.get('measure', 'campaign'), target=target)
    print(f'Outputs: {out}*')

def self_test():
    """Synthetic end-to-end run of the analysis core + tables + figures on a
    3-variable box (no stage-1 data needed): y = 4 x0 x1 + 0.2 x2."""
    class Box:
        d, names = 3, ['k_13', 'ehrlich_downstream', 'k_3']
        def __call__(self, U): return np.ones(len(U), dtype=bool)
    rng = np.random.default_rng(0)
    U = rng.random((400, 3))
    Y = {HEADLINE: 4*U[:, 0]*U[:, 1] + 0.2*U[:, 2] - 1.0, 'IBO titer': U[:, 0] + 0.1*U[:, 2]}
    surrogates, S, fallback = analyse(U, Y, Box(), n_base=2048, n_replicates=2)
    table = index_table(Box.names, surrogates, S)
    subsets = subset_table(Box.names, S[HEADLINE], sizes=(1, 2))
    assert subsets.iloc[0]['size'] == 1 and set(subsets.columns) == {'size', 'parameters', 'S_closed', 'sd'}
    top2 = subsets[subsets['size'] == 2].iloc[0]
    assert top2['parameters'] == 'k_13 + ehrlich_downstream' and top2['S_closed'] > 0.9, top2
    shap = table[(table.metric == HEADLINE) & (table['index'] == 'Shapley')]['mean']
    assert abs(shap.sum() - 1.0) < 1e-9
    # --gp-metrics wiring: only a named metric is GP-eligible, the rest are trees
    assert resolve_gp_metrics(None) == set(DEFAULT_GP_METRICS)
    assert resolve_gp_metrics(['all']) is None and resolve_gp_metrics([]) == set()
    sur_trees, _, _ = analyse(U, {'IBO titer': Y['IBO titer']}, Box(), n_base=64,
                              n_replicates=1, gp_metrics={HEADLINE})
    assert set(sur_trees['IBO titer'].q2) == {'hgb'}, sur_trees['IBO titer'].q2
    # Target / conditional path: threshold 0.5 keeps ~1/4 of the rows in the zone
    thr = 0.5
    zone = Y[HEADLINE] > thr
    assert 0.1 < zone.mean() < 0.5, zone.mean()
    Y[TARGET] = zone.astype(float)
    sur_t, S_t, _ = analyse(U, {TARGET: Y[TARGET]}, Box(), n_base=512, n_replicates=1)
    assert set(sur_t[TARGET].q2) == {'hgb'}
    cond = conditional_analysis(U, Y[HEADLINE], thr, sur_t[TARGET], Box(), n_base=512,
                                n_replicates=1)
    assert cond is not None
    sur_c, S_c, fb_c, mass = cond
    assert 0.1 < mass < 0.5 and S_c.shape == (1, 8) and S_c[0, 7] == 1.0, (mass, S_c.shape)
    assert abs(sa.shapley_effects(S_c, 3).sum() - 1.0) < 1e-9
    surrogates[TARGET], S[TARGET] = sur_t[TARGET], S_t[TARGET]
    surrogates[CONDITIONAL], S[CONDITIONAL] = sur_c, S_c
    table = index_table(Box.names, surrogates, S)
    assert set(table.metric) == {HEADLINE, 'IBO titer', TARGET, CONDITIONAL}
    target = dict(threshold=thr, source='self-test', n_zone=int(zone.sum()),
                  fraction=float(zone.mean()), target_skip=None,
                  conditional=dict(q2=sur_c.q2['hgb'], reliable=sur_c.reliable, mass=mass,
                                   fallback=float(fb_c.max())), conditional_skip=None)
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        plot_bars(table, Box.names, 'self-test', os.path.join(tmp, 'bars'))
        plot_heatmap(table, Box.names, os.path.join(tmp, 'heat'))
        write_summary(os.path.join(tmp, 'summary.txt'), study_name='self-test',
                      counts={'COMPLETE': 400}, n_rows=400, irr_finite=0.06, skipped=[],
                      names=Box.names, surrogates=surrogates, S=S, fallback=fallback,
                      y_headline=Y[HEADLINE], measure='screening', target=target)
        made = sorted(os.listdir(tmp))
        assert made == ['bars.pdf', 'bars.png', 'heat.pdf', 'heat.png', 'summary.txt'], made
        with open(os.path.join(tmp, 'summary.txt'), encoding='utf-8') as fh:
            text = fh.read()
        for needle in ('GP-eligible metrics', 'Input measure: screening',
                       'TARGET / CONDITIONAL', 'Entry (indicator', 'Within the zone'):
            assert needle in text, needle
        # the skip wording when the zone is too thin
        thin = {**target, 'conditional': None,
                'conditional_skip': '3 zone rows < --min-zone-rows 200'}
        write_summary(os.path.join(tmp, 'summary.txt'), study_name='self-test',
                      counts={'COMPLETE': 400}, n_rows=400, irr_finite=0.06, skipped=[],
                      names=Box.names,
                      surrogates={k: v for k, v in surrogates.items() if k != CONDITIONAL},
                      S={k: v for k, v in S.items() if k != CONDITIONAL}, fallback=fallback,
                      y_headline=Y[HEADLINE], target=thin)
        with open(os.path.join(tmp, 'summary.txt'), encoding='utf-8') as fh:
            assert 'Conditional indices skipped: 3 zone rows' in fh.read()
    print('SELF-TEST PASSED')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--study-name')
    parser.add_argument('--metrics', nargs='*', default=None)
    parser.add_argument('--gp-metrics', nargs='*', default=None,
                        help="metrics whose surrogate may be a GP (default: PI, "
                             "'PI (log-tail)'); 'all' = every metric, bare flag = none")
    parser.add_argument('--n-base', type=int, default=4096)
    parser.add_argument('--n-replicates', type=int, default=3)
    parser.add_argument('--min-rows', type=int, default=200)
    parser.add_argument('--convergence-sizes', type=int, nargs='*', default=[1024, 2048])
    parser.add_argument('--pi-threshold', type=float, default=None,
                        help='target / conditional threshold on PI (default: the design '
                             "record's meta.baseline_PI = the scenario-A baseline PI)")
    parser.add_argument('--min-zone-rows', type=int, default=200,
                        help='simulated rows above the threshold needed for the conditional '
                             'indices (default 200)')
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif not args.study_name:
        parser.error('--study-name is required (or --self-test)')
    else:
        run(args)

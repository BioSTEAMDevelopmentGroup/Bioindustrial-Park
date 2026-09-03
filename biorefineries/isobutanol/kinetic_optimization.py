#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Bayesian (Optuna TPE) global optimization of the fermentation kinetic
parameters (k_*/K_* on V406's tellurium model; by default the driver
restricts the set to the scenario workbook's rows via include_params) plus
the feeding-strategy
variables (threshold_conc, target_conc via target_delta, spike_conc via
spike_delta -- feasible-by-construction threshold < target < spike --
and the integer max_n_spikes cap), against a named or custom objective --
to prioritize metabolic-engineering / bioprocess research directions.

Import is free of side effects and does not require optuna (imported lazily
inside run_kinetic_optimization); the pure logic here is exercised by
analyses/test_kinetic_optimization_offline.py without a biorefinery load.
The engine itself (run_kinetic_optimization) requires
biorefineries.isobutanol.load() to have been called first; the canonical
entry point is analyses/optimize_kinetics_BO.py.

Every trial runs the FULL system simulation + one side-effect-free TEA
solve, for kinetic-level and system-level objectives alike (the 'level'
tag is metadata). Kinetic parameters and feeding specs are restored to
their scenario baselines in a `finally` after every run.
"""
import csv
import math
import os
import re

import numpy as np

__all__ = ('OBJECTIVE_REGISTRY', 'TRACKED_METRICS',
           'discover_kinetic_parameters', 'build_search_space',
           'parameter_distributions_workbook', 'workbook_kinetic_baselines',
           'kinetic_param_names_from_scenario',
           'baseline_decision_point',
           'trajectory_columns', 'append_trajectory_row', 'load_trajectory',
           'get_handles', 'run_kinetic_optimization', 'restore_baseline',
           'plot_optimization_trajectories', 'plot_parameter_trajectory',
           'plot_best_vs_baseline',
           'pca_decision_matrix', 'plot_pca_projection',
           'StallGuard', 'attempt_outcome')

FEEDING_VARIABLES = ('threshold_conc', 'target_delta', 'spike_delta',
                     'max_n_spikes')
#: Feeding decision variables of studies started before 2026-08-31
#: (target-anchored parameterization); resumable via the legacy bounds
#: kwargs of build_search_space / run_kinetic_optimization.
LEGACY_FEEDING_VARIABLES = ('target_conc', 'threshold_delta', 'spike_conc')

#: Applied-concentration envelope (g/L) for the current parameterization:
#: target_conc = min(TARGET_CONC_MAX, threshold_conc + target_delta);
#: spike_conc = clip(target_conc + spike_delta, SPIKE_CONC_MIN,
#: SPIKE_CONC_MAX).
TARGET_CONC_MAX = 500.0
SPIKE_CONC_MIN = 50.0
SPIKE_CONC_MAX = 600.0

#%% Objective registry and tracked metrics
# Getters are callables over a `handles` dict (see get_handles below):
#   'V406': the fermentation unit (.nsk_results_specific_tau_dict, .tau)
#   'tea': the system TEA (.TCI)
#   'latest_TEA_solution': {'IRR': ..., 'MPSPs': {'ethanol': ..,
#                           'isobutanol': ..}}, refreshed once per trial.
# This indirection keeps every getter testable offline with fakes.

def _nsk(handles):
    return handles['V406'].nsk_results_specific_tau_dict

OBJECTIVE_REGISTRY = {
    'IBO yield': dict(
        getter=lambda h: _nsk(h)['y_IBO_glu_added'],
        direction='maximize', level='kinetic', units='g-IBO/g-sugars'),
    'IBO titer': dict(
        getter=lambda h: _nsk(h)['[s_IBO]'],
        direction='maximize', level='kinetic', units='g-IBO/L-broth'),
    'IBO productivity': dict(
        getter=lambda h: _nsk(h)['[s_IBO]']/_nsk(h)['time'],
        direction='maximize', level='kinetic', units='g-IBO/L-broth/h'),
    'EtOH yield': dict(
        getter=lambda h: _nsk(h)['y_EtOH_glu_added'],
        direction='maximize', level='kinetic', units='g-EtOH/g-sugars'),
    'EtOH titer': dict(
        getter=lambda h: _nsk(h)['[s_EtOH]'],
        direction='maximize', level='kinetic', units='g-EtOH/L-broth'),
    'EtOH productivity': dict(
        getter=lambda h: _nsk(h)['prod_EtOH'],
        direction='maximize', level='kinetic', units='g-EtOH/L-broth/h'),
    'Combined yield': dict(
        getter=lambda h: _nsk(h)['y_EtOH_IBO_glu_added'],
        direction='maximize', level='kinetic', units='g-EtOH-and-IBO/g-sugars'),
    'Cell density': dict(
        getter=lambda h: _nsk(h)['[x]'],
        direction='maximize', level='kinetic', units='g-cell/L-broth'),
    'IRR': dict(
        getter=lambda h: h['latest_TEA_solution']['IRR'],
        direction='maximize', level='system', units=''),
    'EtOH MPSP': dict(
        getter=lambda h: h['latest_TEA_solution']['MPSPs']['ethanol'],
        direction='minimize', level='system', units='$/kg'),
    'IBO MPSP': dict(
        getter=lambda h: h['latest_TEA_solution']['MPSPs']['isobutanol'],
        direction='minimize', level='system', units='$/kg'),
    'TCI': dict(
        getter=lambda h: h['tea'].TCI/1e6,
        direction='minimize', level='system', units='MM$'),
    }

#: Metrics recorded for EVERY trial (spec trajectory (ii)-(vi) + extras).
TRACKED_METRICS = {name: OBJECTIVE_REGISTRY[name]['getter'] for name in
                   ('IBO yield', 'IBO titer', 'IBO productivity',
                    'EtOH yield', 'EtOH titer', 'EtOH productivity',
                    'Cell density', 'IRR', 'TCI')}
TRACKED_METRICS['tau'] = lambda h: h['V406'].tau
TRACKED_METRICS['n_glu_spikes'] = lambda h: _nsk(h)['curr_n_glu_spikes']

#%% Search space

def discover_kinetic_parameters(r_te):
    """{name: baseline} for every kinetic parameter of the tellurium model
    `r_te` -- the k_* rates and K_* constants, selected with the same name
    rule as utils.generate_save_kinetic_parameter_distributions
    (name[:2].lower() == 'k_'). Baselines are the CURRENT values, so call
    this only after the scenario's parameter distributions have been
    loaded (model.load_parameter_distributions + metrics_at_baseline)."""
    return {p: float(getattr(r_te, p))
            for p in r_te.getGlobalParameterIds()
            if p[:2].lower() == 'k_'}

def build_search_space(kinetic_baselines,
                       multiplier_bounds=(0.1, 10.0),
                       param_bounds_override=None,
                       exclude_params=(),
                       include_params=None,
                       threshold_conc_bounds=(0.0, 500.0),
                       target_delta_bounds=(5.0, 500.0),
                       spike_delta_bounds=(0.5, 595.0),
                       max_n_spikes_bounds=(0, 50),
                       target_conc_bounds=None,
                       threshold_delta_bounds=None,
                       spike_conc_bounds=None,
                       ):
    """Build the decision-variable space: {name: {'low', 'high', 'log'}}
    (integer variables additionally carry 'int': True).

    Kinetic parameters default to [m_lo*baseline, m_hi*baseline] sampled
    log-scale (research-opportunity width); entries in
    `param_bounds_override` ({name: (low, high)}) use those absolute
    bounds instead (log-scale only if low > 0). A parameter with a
    nonpositive baseline and no override cannot use the multiplier band
    and is EXCLUDED with a printed warning. `exclude_params` names are
    always excluded (silently). `include_params`
    (None = no restriction) is a WHITELIST: a kinetic parameter is placed
    in the space only if its name is in it, and every other kinetic
    parameter is excluded silently -- the driver passes the kinetic rows
    of a scenario's parameter-distributions workbook here so the search
    set matches the curated uncertainty set. The feeding variables below
    are never filtered by it. The feeding variables (all linear) make
    the required ordering threshold < target < spike feasible BY
    CONSTRUCTION: threshold_conc is sampled absolutely, the target is
    sampled as target_delta above the threshold
    (target_conc = min(TARGET_CONC_MAX, threshold_conc + target_delta)),
    and the spike as spike_delta above the target
    (spike_conc = clip(target_conc + spike_delta, SPIKE_CONC_MIN,
    SPIKE_CONC_MAX)) -- spanning the applied envelope threshold [0, 500],
    target [5, 500], spike [50, 600] g/L at the default bounds.
    max_n_spikes (the glucose-spike cap, fbs_spec.max_n_spikes) is an
    INTEGER variable (0 = forced batch); pass max_n_spikes_bounds=None to
    pin it at the scenario baseline instead.

    Passing any of the LEGACY kwargs (target_conc_bounds,
    threshold_delta_bounds, spike_conc_bounds; unspecified ones fall back
    to the pre-2026-08-31 defaults (180, 300)/(0.5, 30)/(200, 800))
    builds the legacy target-anchored parameterization instead
    (target_conc absolute, threshold_conc = max(0, target_conc -
    threshold_delta), spike_conc absolute) -- required to resume a study
    started under it.

    Returns (space, excluded_parameter_names)."""
    param_bounds_override = dict(param_bounds_override or {})
    m_lo, m_hi = multiplier_bounds
    space, excluded = {}, []
    for name, baseline in kinetic_baselines.items():
        if include_params is not None and name not in include_params:
            excluded.append(name)
        elif name in exclude_params:
            excluded.append(name)
        elif name in param_bounds_override:
            lo, hi = param_bounds_override[name]
            space[name] = dict(low=lo, high=hi, log=lo > 0.0)
        elif baseline <= 0.0:
            print(f'Warning: kinetic parameter {name} has a nonpositive '
                  f'baseline ({baseline}) and no param_bounds_override '
                  'entry; excluding it from the search space.')
            excluded.append(name)
        else:
            space[name] = dict(low=m_lo*baseline, high=m_hi*baseline,
                               log=True)
    if (target_conc_bounds is not None or threshold_delta_bounds is not None
            or spike_conc_bounds is not None):
        # Legacy target-anchored parameterization (resumes of studies
        # started before 2026-08-31).
        tcb = target_conc_bounds or (180.0, 300.0)
        tdb = threshold_delta_bounds or (0.5, 30.0)
        scb = spike_conc_bounds or (200.0, 800.0)
        space['target_conc'] = dict(low=tcb[0], high=tcb[1], log=False)
        space['threshold_delta'] = dict(low=tdb[0], high=tdb[1], log=False)
        space['spike_conc'] = dict(low=scb[0], high=scb[1], log=False)
    else:
        space['threshold_conc'] = dict(low=threshold_conc_bounds[0],
                                       high=threshold_conc_bounds[1],
                                       log=False)
        space['target_delta'] = dict(low=target_delta_bounds[0],
                                     high=target_delta_bounds[1],
                                     log=False)
        space['spike_delta'] = dict(low=spike_delta_bounds[0],
                                    high=spike_delta_bounds[1], log=False)
    if max_n_spikes_bounds is not None:
        space['max_n_spikes'] = dict(low=int(max_n_spikes_bounds[0]),
                                     high=int(max_n_spikes_bounds[1]),
                                     log=False, int=True)
    return space, excluded

def baseline_decision_point(search_space, kinetic_baselines,
                            baseline_model_kwargs,
                            baseline_max_n_spikes=None):
    """The scenario baseline expressed in decision-variable coordinates
    for `search_space` (either feeding parameterization) -- suitable for
    study.enqueue_trial, so a fresh study evaluates the baseline itself
    as trial 0. Only names present in the search space are included.
    Every value is CLIPPED into its [low, high] bounds: an out-of-range
    enqueued value is a hard optuna ValueError on a log-scale variable
    and is silently replaced by a random draw on a linear one, so when a
    baseline lies outside the space (e.g. a zero-baseline kinetic
    parameter under absolute param_bounds_override bounds), trial 0
    evaluates the nearest in-bounds point to the baseline instead."""
    point = {name: kinetic_baselines[name]
             for name in search_space if name in kinetic_baselines}
    thr = baseline_model_kwargs['threshold_conc']
    tgt = baseline_model_kwargs['target_conc']
    spk = baseline_model_kwargs['spike_conc']
    if 'threshold_conc' in search_space:  # threshold-anchored scheme
        point['threshold_conc'] = thr
        point['target_delta'] = tgt - thr
        point['spike_delta'] = spk - tgt
    elif 'target_conc' in search_space:  # legacy target-anchored scheme
        point['target_conc'] = tgt
        point['threshold_delta'] = tgt - thr
        point['spike_conc'] = spk
    if 'max_n_spikes' in search_space and baseline_max_n_spikes is not None:
        point['max_n_spikes'] = int(baseline_max_n_spikes)
    for name, value in point.items():
        sp = search_space[name]
        clipped = min(max(value, sp['low']), sp['high'])
        point[name] = int(clipped) if sp.get('int') else clipped
    return point

#%% Scenario parameter-distribution workbooks

#: Load-statement pattern of a kinetic parameter row in the scenario
#: parameter-distributions workbooks (written by
#: utils.generate_save_kinetic_parameter_distributions).
_TE_LOAD_STATEMENT = re.compile(
    r'V406\.nsk_kinetic_model\._te\.([A-Za-z0-9_]+)\s*=\s*x')

def parameter_distributions_workbook(scenario):
    """Absolute path of the scenario's parameter-distributions workbook
    (analyses/full/parameter_distributions/
    parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx)."""
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'analyses', 'full', 'parameter_distributions',
        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx')

def workbook_kinetic_baselines(scenario):
    """Ordered {te parameter name: workbook Baseline} for every kinetic
    parameter row (load statement `V406.nsk_kinetic_model._te.<name> = x`)
    of the scenario's parameter-distributions workbook -- the curated set
    of kinetic parameters treated as free/uncertain (the physically
    constrained k_6r, k_16r, K_2, K_9 are absent since commit 1e4efee1).
    A plain file read, no simulation, so it can describe a scenario other
    than the loaded one (the driver's start-at-A / set-from-B mode)."""
    import pandas as pd
    baselines = {}
    for _, row in pd.read_excel(
            parameter_distributions_workbook(scenario)).iterrows():
        match = _TE_LOAD_STATEMENT.fullmatch(
            str(row['Load statement']).strip())
        if match:
            baselines[match.group(1)] = float(row['Baseline'])
    return baselines

def kinetic_param_names_from_scenario(scenario):
    """The kinetic parameter names of the scenario's workbook, in workbook
    order (see workbook_kinetic_baselines) -- the include_params of a
    workbook-restricted search space."""
    return list(workbook_kinetic_baselines(scenario))

#%% Trajectory recording

def trajectory_columns(search_space):
    """Column order of the trajectory CSV for a given search space."""
    return ['trial_number', 'state', *search_space.keys(), 'objective',
            *TRACKED_METRICS.keys(), 'error']

def append_trajectory_row(csv_path, columns, record):
    """Append one row (dict; missing keys become '') to `csv_path`,
    writing the header first if the file does not exist. An existing
    file's header must match `columns` exactly -- otherwise the search
    space changed since the trajectory was started, and appending would
    silently misalign the rows. Flushed immediately, so a crash/segfault
    loses at most the in-flight trial."""
    exists = os.path.isfile(csv_path)
    if exists:
        with open(csv_path, newline='') as csvfile:
            existing_header = next(csv.reader(csvfile), None)
        if existing_header != list(columns):
            raise ValueError(
                f'Trajectory CSV {csv_path} has a different column set '
                'than the current search space -- the search space '
                'changed since this trajectory was started. Use a new '
                'study_name (or move the old CSV and .db) to start '
                'fresh.')
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns,
                                extrasaction='ignore')
        if not exists:
            writer.writeheader()
        writer.writerow({k: record.get(k, '') for k in columns})
        csvfile.flush()

def load_trajectory(csv_path):
    """Read a trajectory CSV back as a pandas DataFrame (numeric columns
    parsed; blank cells -> NaN)."""
    import pandas as pd
    return pd.read_csv(csv_path)

#%% Post-run plots
# All three take the trajectory DataFrame (load_trajectory(csv_path)), so
# plots can be regenerated from a saved CSV without re-running anything.
# matplotlib/pandas are imported lazily (headless Agg per the environment).

def _best_so_far(objective_series, direction):
    """Best-so-far series: cumulative max ('maximize') or min
    ('minimize') of the objective, NaNs carried over."""
    return (objective_series.cummax() if direction == 'maximize'
            else objective_series.cummin())

def _completed(df):
    ok = df[df['state'] == 'COMPLETE'].reset_index(drop=True)
    if ok.empty:
        raise ValueError('No completed trials in the trajectory -- '
                         'nothing to plot.')
    return ok

def plot_optimization_trajectories(df, objective_name, direction,
                                   objective_units='', filename=None):
    """Multipanel trajectory figure: the objective (all completed trials
    as scatter + best-so-far line) and every tracked metric vs trial
    number, each metric panel overlaying the metric's value in the
    incumbent best-objective-so-far configuration (NOT that metric's own
    best-so-far). Failed/pruned trials are omitted. Returns (fig, axes)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    best_rows = ok.iloc[_best_row_indices(ok, direction)].reset_index(
        drop=True)
    metric_names = [m for m in TRACKED_METRICS if m in ok.columns]
    n_panels = 1 + len(metric_names)
    ncols = 4
    nrows = math.ceil(n_panels/ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows),
                             squeeze=False)
    flat = axes.ravel()
    ax = flat[0]
    ax.scatter(ok['trial_number'], ok['objective'], s=8, alpha=0.4,
               color='tab:gray', label='trials')
    ax.plot(ok['trial_number'], _best_so_far(ok['objective'], direction),
            color='tab:red', lw=1.5, label='best so far')
    title = f'Objective: {objective_name}'
    if objective_units:
        title += f' ({objective_units})'
    ax.set_title(title)
    ax.set_xlabel('trial')
    ax.legend(fontsize=8)
    for ax, m in zip(flat[1:], metric_names):
        ax.scatter(ok['trial_number'], ok[m], s=8, alpha=0.4,
                   color='tab:blue')
        ax.plot(ok['trial_number'], best_rows[m], color='tab:red', lw=1.5,
                label='at best-so-far objective')
        ax.set_title(m)
        ax.set_xlabel('trial')
        ax.legend(fontsize=7)
    for ax in flat[n_panels:]:
        ax.axis('off')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, axes

def _applied_feeding(rows):
    """(target, threshold, spike) APPLIED concentrations from trajectory
    rows (a DataFrame or a single-row Series), handling both the current
    threshold-anchored parameterization (threshold_conc/target_delta/
    spike_delta) and the legacy target-anchored one (target_conc/
    threshold_delta/spike_conc)."""
    if 'threshold_conc' in rows:
        threshold = rows['threshold_conc']
        target = np.minimum(TARGET_CONC_MAX,
                            threshold + rows['target_delta'])
        spike = np.minimum(SPIKE_CONC_MAX,
                           np.maximum(SPIKE_CONC_MIN,
                                      target + rows['spike_delta']))
    else:
        target = rows['target_conc']
        threshold = np.maximum(0.0, target - rows['threshold_delta'])
        spike = rows['spike_conc']
    return target, threshold, spike

def _best_row_indices(ok, direction):
    """Row index (into `ok`) of the incumbent-best trial as of each
    completed trial."""
    best_idx, cur = [], None
    values = ok['objective'].to_list()
    for i, v in enumerate(values):
        if cur is None or (direction == 'maximize' and v > values[cur]) \
                or (direction == 'minimize' and v < values[cur]):
            cur = i
        best_idx.append(cur)
    return best_idx

def plot_parameter_trajectory(df, kinetic_baselines, direction,
                              filename=None):
    """Best-so-far kinetic CONFIGURATION vs trial number: for each kinetic
    parameter, the incumbent's multiplier (value/baseline, log y) at every
    completed trial, plus a companion panel with the three feeding
    concentrations in absolute units (threshold shown as
    target - threshold_delta). Parameters with nonpositive baselines
    (multiplier undefined; e.g. absolute-bounds overrides of zero-baseline
    params) are skipped here -- they still live in the CSV. Returns
    (fig, axes)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    best_rows = ok.iloc[_best_row_indices(ok, direction)].reset_index(
        drop=True)
    x = ok['trial_number']
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True,
                                   height_ratios=[3, 1])
    for pname, baseline in kinetic_baselines.items():
        if pname not in ok.columns or baseline <= 0:
            continue
        ax1.plot(x, best_rows[pname]/baseline, lw=1, label=pname)
    ax1.set_yscale('log')
    ax1.axhline(1.0, color='k', lw=0.8, ls='--')
    ax1.set_ylabel('best-so-far multiplier vs baseline')
    ax1.legend(fontsize=6, ncol=8, loc='upper center',
               bbox_to_anchor=(0.5, -0.08))
    if 'target_conc' in ok.columns or 'threshold_conc' in ok.columns:
        target, threshold, spike = _applied_feeding(best_rows)
        ax2.plot(x, target, label='target_conc')
        ax2.plot(x, threshold, label='threshold_conc')
        ax2.plot(x, spike, label='spike_conc')
        ax2.legend(fontsize=7)
    ax2.set_ylabel('g/L')
    ax2.set_xlabel('trial')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200, bbox_inches='tight')
    return fig, (ax1, ax2)

def plot_best_vs_baseline(df, kinetic_baselines, direction, filename=None):
    """The research-prioritization headline: horizontal bars of the FINAL
    incumbent's kinetic-parameter multipliers (log x, baseline = 1 dashed
    line), sorted by multiplier, with the optimal feeding concentrations
    and objective value in the title. Nonpositive-baseline parameters are
    skipped (see plot_parameter_trajectory). Returns (fig, ax)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    i_best = (ok['objective'].idxmax() if direction == 'maximize'
              else ok['objective'].idxmin())
    best = ok.loc[i_best]
    names = [p for p, b in kinetic_baselines.items()
             if p in ok.columns and b > 0]
    multipliers = np.array([best[p]/kinetic_baselines[p] for p in names])
    order = np.argsort(multipliers)
    fig, ax = plt.subplots(figsize=(7, 0.28*len(names) + 2))
    ax.barh([names[i] for i in order], multipliers[order],
            color='tab:blue')
    ax.set_xscale('log')
    ax.axvline(1.0, color='k', ls='--', lw=0.8)
    ax.set_xlabel('best/baseline multiplier')
    target, threshold, spike = _applied_feeding(best)
    ax.set_title(f"objective = {best['objective']:.4g} at trial "
                 f"{int(best['trial_number'])}; target = {target:.1f}, "
                 f"threshold = {threshold:.1f}, "
                 f"spike = {spike:.1f} g/L", fontsize=8)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, ax

#%% PCA projection of the sampled decision space

def pca_decision_matrix(df, log_columns=(), decision_columns=None):
    """PCA of the decision vectors of EVERY recorded trial (COMPLETE, FAIL
    and NAN rows alike -- each records its full decision vector), from a
    trajectory DataFrame (load_trajectory(csv_path)).

    `decision_columns` defaults to the columns between 'state' and
    'objective' -- exactly the search-space variables, by
    trajectory_columns construction. Columns named in `log_columns` (the
    log-sampled kinetic multiplier bands) are log10-transformed; all are
    then z-scored, zero-variance columns dropped (e.g. pinned variables),
    and the PCA computed by SVD (no sklearn). Component signs are
    stabilized (largest-|loading| entry positive) so successive mid-run
    figures do not flip axes.

    Returns (coords, explained_var_ratio, loadings, kept_columns, valid):
    `coords` (n_valid x n_components) are the trial scores, `loadings`
    (n_components x n_kept) the variable weights, and `valid` a boolean
    mask aligned with df's rows (False where any decision value is
    missing/non-finite after transformation)."""
    if decision_columns is None:
        cols = list(df.columns)
        decision_columns = cols[cols.index('state') + 1:
                                cols.index('objective')]
    # .copy(): pandas may hand back a read-only zero-copy view here, and
    # the log10 transform below writes in place.
    T = df[list(decision_columns)].to_numpy(dtype=float).copy()
    with np.errstate(divide='ignore', invalid='ignore'):
        for j, name in enumerate(decision_columns):
            if name in log_columns:
                T[:, j] = np.log10(T[:, j])
    valid = np.isfinite(T).all(axis=1)
    if valid.sum() < 3:
        raise ValueError('Fewer than 3 trials with complete decision '
                         'vectors -- nothing to project.')
    Tv = T[valid]
    kept, columns_z = [], []
    for j, name in enumerate(decision_columns):
        col = Tv[:, j]
        std = col.std()
        if std > 0.0 and np.isfinite(std):
            columns_z.append((col - col.mean())/std)
            kept.append(name)
    if len(kept) < 2:
        raise ValueError('Fewer than 2 non-degenerate decision variables '
                         '-- nothing to project.')
    Z = np.column_stack(columns_z)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    signs = np.sign(Vt[np.arange(Vt.shape[0]),
                       np.argmax(np.abs(Vt), axis=1)])
    signs[signs == 0] = 1.0
    Vt = Vt*signs[:, None]
    coords = U*S*signs[None, :]
    explained_var_ratio = S**2/(S**2).sum()
    return coords, explained_var_ratio, Vt, kept, valid

def plot_pca_projection(df, direction, log_columns=(),
                        objective_name='objective', objective_units='',
                        filename=None):
    """Four-panel PCA view of the sampled decision space: (1) the PC1 x
    PC2 landscape -- completed trials colored by objective, FAIL/NAN
    trials as gray/red crosses, the incumbent best-so-far path, the
    enqueued baseline (trial 0) and the current best marked; (2) the
    explained-variance scree of the top 10 PCs; (3) the top-|loading|
    variables on PC1/PC2; (4) PC1 and PC2 of every sampled point vs trial
    number (the sampler-contraction diagnostic). Works with zero
    completed trials (landscape shows only pruned points). Returns
    (fig, axes)."""
    import matplotlib.pyplot as plt
    coords, evr, loadings, kept, valid = pca_decision_matrix(
        df, log_columns=log_columns)
    dfv = df[valid].reset_index(drop=True)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.2, 1.2, 1],
                          height_ratios=[1, 1, 0.7])
    ax_main = fig.add_subplot(gs[0:2, 0:2])
    ax_scree = fig.add_subplot(gs[0, 2])
    ax_load = fig.add_subplot(gs[1:3, 2])
    ax_time = fig.add_subplot(gs[2, 0:2])

    pc1, pc2 = coords[:, 0], coords[:, 1]
    state = dfv['state']
    for st, color, label in (('FAIL', 'tab:gray', 'failed (pruned)'),
                             ('NAN', 'tab:red', 'NaN objective (pruned)')):
        m = (state == st).to_numpy()
        if m.any():
            ax_main.scatter(pc1[m], pc2[m], marker='x', s=18, alpha=0.35,
                            color=color, label=label)
    ok = (state == 'COMPLETE').to_numpy()
    if ok.any():
        sc = ax_main.scatter(pc1[ok], pc2[ok], c=dfv['objective'][ok],
                             cmap='viridis', s=16, alpha=0.85,
                             label='completed')
        clabel = objective_name + (f' ({objective_units})'
                                   if objective_units else '')
        fig.colorbar(sc, ax=ax_main, label=clabel)
        ok_df = dfv[ok].reset_index(drop=True)
        best_path = np.unique(_best_row_indices(ok_df, direction))
        ax_main.plot(pc1[ok][best_path], pc2[ok][best_path], ls='--',
                     color='tab:red', lw=1.2, marker='D', ms=4,
                     label='best-so-far path')
        i_best = (ok_df['objective'].idxmax() if direction == 'maximize'
                  else ok_df['objective'].idxmin())
        ax_main.scatter(pc1[ok][i_best], pc2[ok][i_best], marker='*',
                        s=260, color='gold', edgecolor='k', zorder=5,
                        label='current best')
    m0 = (dfv['trial_number'] == 0).to_numpy()
    if m0.any():
        ax_main.scatter(pc1[m0], pc2[m0], marker='*', s=200, color='k',
                        zorder=5, label='baseline (trial 0)')
    ax_main.set_xlabel(f'PC1 ({100*evr[0]:.1f}% var)')
    ax_main.set_ylabel(f'PC2 ({100*evr[1]:.1f}% var)')
    ax_main.set_title(f'Sampled decision space ({len(kept)} variables), '
                      f'{len(dfv)} trials')
    ax_main.legend(fontsize=7, loc='best')

    n_show = min(10, len(evr))
    ax_scree.bar(np.arange(1, n_show + 1), evr[:n_show],
                 color='tab:blue', alpha=0.8)
    ax_scree.plot(np.arange(1, n_show + 1), np.cumsum(evr[:n_show]),
                  color='tab:red', marker='.', lw=1, label='cumulative')
    ax_scree.set_xlabel('PC')
    ax_scree.set_title('explained variance', fontsize=9)
    ax_scree.legend(fontsize=7)

    strength = np.hypot(loadings[0], loadings[1])
    top = np.argsort(strength)[::-1][:12][::-1]
    y = np.arange(len(top))
    ax_load.barh(y + 0.2, loadings[0][top], height=0.4,
                 color='tab:blue', label='PC1')
    ax_load.barh(y - 0.2, loadings[1][top], height=0.4,
                 color='tab:orange', label='PC2')
    ax_load.set_yticks(y)
    ax_load.set_yticklabels([kept[i] for i in top], fontsize=7)
    ax_load.axvline(0.0, color='k', lw=0.8)
    ax_load.set_title('top loadings', fontsize=9)
    ax_load.legend(fontsize=7)

    x = dfv['trial_number']
    ax_time.scatter(x, pc1, s=6, alpha=0.4, color='tab:blue',
                    label='PC1')
    ax_time.scatter(x, pc2, s=6, alpha=0.4, color='tab:orange',
                    label='PC2')
    ax_time.set_xlabel('trial')
    ax_time.set_ylabel('PC coordinate')
    ax_time.set_title('sampler contraction', fontsize=9)
    ax_time.legend(fontsize=7)

    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, (ax_main, ax_scree, ax_load, ax_time)

#%% Process-level trial-timeout supervision (pure logic)
# A pathological kinetic draw can hang the simulation INSIDE a native
# CVODE/roadrunner integrator call (observed 2026-08-31: one trial at
# 100% of a core for ~2.7 h with no output). No in-process mechanism can
# reliably interrupt that on Windows -- watchdog threads and
# PyThreadState_SetAsyncExc only act at Python bytecode boundaries -- so
# the per-trial wall-clock timeout is implemented one level up: a
# supervisor process (analyses/optimize_kinetics_BO_supervised.py) polls
# the crash-safe trajectory CSV, kills the run when no trial has been
# recorded for the timeout, and relaunches the resumable study (the
# in-flight trial is lost, exactly as in a segfault crash; its RUNNING
# optuna record still counts toward the total budget). The decision
# logic lives here, dependency-free, so the offline test covers it.

class StallGuard:
    """Stall detector over trajectory-CSV row counts. Feed it one
    (rows, now) observation per poll; `update` returns 'stalled' once no
    new row has appeared for `stall_timeout_s` (progress resets the
    clock), else None. Timebase: any monotonic seconds. Call `reset()`
    before each new supervised attempt."""
    def __init__(self, stall_timeout_s=1500.0):
        if stall_timeout_s <= 0:
            raise ValueError('stall_timeout_s must be positive; '
                             f'got {stall_timeout_s!r}')
        self.stall_timeout_s = stall_timeout_s
        self.reset()

    def reset(self):
        self._last_rows = None
        self._last_progress_t = None

    def update(self, rows, now):
        if self._last_rows is None or rows > self._last_rows:
            self._last_rows = rows
            self._last_progress_t = now
            return None
        if now - self._last_progress_t >= self.stall_timeout_s:
            return 'stalled'
        return None

def attempt_outcome(exit_code, rows_before, rows_after,
                    killed_for_stall=False):
    """Supervisor decision after one attempt exits: 'complete' (clean
    exit 0), 'resume' (crash or stall-kill, but the attempt recorded new
    trials -- relaunch the resumable study), or 'abort' (an attempt that
    added no rows; resuming would loop on the same failure, a human
    should look)."""
    if exit_code == 0 and not killed_for_stall:
        return 'complete'
    if rows_after > rows_before:
        return 'resume'
    return 'abort'

#%% Engine

def get_handles():
    """Resolve the live flowsheet objects the getters and engine need.
    Requires biorefineries.isobutanol.load() to have run in this kernel
    (the system module raises an informative error otherwise)."""
    from biorefineries.isobutanol import system as ibo_system
    f = ibo_system.f
    V406 = f.V406
    return {'f': f,
            'V406': V406,
            'r_te': V406.nsk_kinetic_model._te,
            'fbs_spec': V406.fbs_spec,
            'tea': ibo_system.corn_EtOH_IBO_sys_tea,
            'HXN': f.HXN1001,
            'model_specification': ibo_system.model_specification,
            'solve_TEA': ibo_system.solve_TEA,
            'latest_TEA_solution': {
                'IRR': np.nan,
                'MPSPs': {'ethanol': np.nan, 'isobutanol': np.nan}},
            }

def restore_baseline(handles, kinetic_baselines, baseline_model_kwargs,
                     baseline_max_n_spikes=None):
    """Reset every kinetic parameter to its recorded baseline (and, when
    `baseline_max_n_spikes` is given, the glucose-spike cap
    fbs_spec.max_n_spikes), re-simulate at the baseline feeding
    specifications, and refresh the TEA solution -- leaving the process in
    a clean scenario-baseline state. Called in run_kinetic_optimization's
    `finally` (success, exception, or KeyboardInterrupt alike)."""
    r_te = handles['r_te']
    for pname, baseline in kinetic_baselines.items():
        setattr(r_te, pname, baseline)
    if baseline_max_n_spikes is not None:
        handles['fbs_spec'].max_n_spikes = baseline_max_n_spikes
    handles['model_specification'](**baseline_model_kwargs)
    handles['latest_TEA_solution'].update(
        handles['solve_TEA'](stream_IDs=('ethanol', 'isobutanol')))
    print('Restored kinetic parameters and feeding specifications to '
          'baseline. Baseline TEA solution: '
          f"{handles['latest_TEA_solution']}")

def run_kinetic_optimization(objective='IRR',
                             direction=None, level=None,
                             objective_units=None, objective_name=None,
                             scenario_label='B',
                             n_trials=2000, seed=3221,
                             multiplier_bounds=(0.1, 10.0),
                             param_bounds_override=None,
                             exclude_params=(),
                             threshold_conc_bounds=(0.0, 500.0),
                             target_delta_bounds=(5.0, 500.0),
                             spike_delta_bounds=(0.5, 595.0),
                             max_n_spikes_bounds=(0, 50),
                             target_conc_bounds=None,
                             threshold_delta_bounds=None,
                             spike_conc_bounds=None,
                             study_name=None, results_dir=None,
                             handles=None, print_status_every=1,
                             ):
    """Run the Bayesian optimization. `objective` is a name in
    OBJECTIVE_REGISTRY (direction/level/units filled from the entry) or a
    custom callable(handles)->float (then `direction`, and ideally
    `objective_name`/`level`/`objective_units`, must be given).

    `n_trials` is the TOTAL budget of the study: rerunning with the same
    study_name resumes from the on-disk SQLite store and runs only the
    remainder (crash/segfault recovery). A FRESH study evaluates the
    scenario baseline configuration as trial 0 (see
    baseline_decision_point); resumes never re-enqueue it. Trials execute STRICTLY
    sequentially (n_jobs=1; one simulation in flight at a time). Every
    trial appends one row to the trajectory CSV (same stable name as the
    study, '_trajectory.csv' suffix) whether it completes, fails
    (state='FAIL', pruned), or yields a NaN objective (state='NAN',
    pruned). Kinetic parameters and feeding specs are restored to the
    scenario baseline in a `finally`.

    Returns (study, csv_path, kinetic_baselines)."""
    import optuna
    if handles is None:
        handles = get_handles()
    r_te, fbs_spec = handles['r_te'], handles['fbs_spec']

    if isinstance(objective, str):
        entry = OBJECTIVE_REGISTRY[objective]
        objective_getter = entry['getter']
        objective_name = objective_name or objective
        direction = direction or entry['direction']
        level = level or entry['level']
        if objective_units is None:
            objective_units = entry['units']
    else:
        objective_getter = objective
        objective_name = objective_name or 'custom'
        if direction not in ('maximize', 'minimize'):
            raise ValueError("A custom objective callable requires "
                             "direction='maximize' or 'minimize'.")

    kinetic_baselines = discover_kinetic_parameters(r_te)
    search_space, excluded = build_search_space(
        kinetic_baselines,
        multiplier_bounds=multiplier_bounds,
        param_bounds_override=param_bounds_override,
        exclude_params=exclude_params,
        threshold_conc_bounds=threshold_conc_bounds,
        target_delta_bounds=target_delta_bounds,
        spike_delta_bounds=spike_delta_bounds,
        max_n_spikes_bounds=max_n_spikes_bounds,
        target_conc_bounds=target_conc_bounds,
        threshold_delta_bounds=threshold_delta_bounds,
        spike_conc_bounds=spike_conc_bounds)
    n_kinetic = sum(1 for name in search_space if name in kinetic_baselines)
    n_feeding = len(search_space) - n_kinetic
    print(f'Search space: {len(search_space)} decision variables '
          f'({n_kinetic} kinetic + {n_feeding} feeding); '
          f'{len(excluded)} kinetic parameters excluded: {excluded}')

    # Scenario baseline snapshot for restoration (the driver has already
    # baseline-simulated the scenario, so current_specifications IS the
    # scenario baseline).
    baseline_model_kwargs = {
        k: fbs_spec.current_specifications[k]
        for k in ('target_conc', 'threshold_conc', 'spike_conc')}
    baseline_max_n_spikes = fbs_spec.max_n_spikes

    if results_dir is None:
        results_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'analyses', 'results')
    os.makedirs(results_dir, exist_ok=True)
    slug = objective_name.lower().replace(' ', '_')
    if study_name is None:
        study_name = f'kin_opt_{scenario_label}_{slug}'
    csv_path = os.path.join(results_dir, study_name + '_trajectory.csv')
    storage = ('sqlite:///'
               + os.path.join(results_dir, study_name + '.db')
               .replace('\\', '/'))
    columns = trajectory_columns(search_space)

    study = optuna.create_study(study_name=study_name, storage=storage,
                                direction=direction,
                                load_if_exists=True)
    n_done = len(study.trials)
    # Offset the seed by the number of stored trials so a resumed study
    # draws fresh points instead of replaying the original RNG stream.
    study.sampler = optuna.samplers.TPESampler(
        multivariate=True, seed=seed + n_done,
        n_startup_trials=max(10, n_trials//10))
    if n_done == 0:
        # Fresh study: evaluate the scenario baseline itself as trial 0,
        # so the baseline provably participates and TPE learns from it.
        study.enqueue_trial(baseline_decision_point(
            search_space, kinetic_baselines, baseline_model_kwargs,
            baseline_max_n_spikes))
        print('Enqueued the scenario baseline configuration as trial 0.')

    def _objective(trial):
        values = {name: (trial.suggest_int(name, sp['low'], sp['high'])
                         if sp.get('int') else
                         trial.suggest_float(name, sp['low'], sp['high'],
                                             log=sp['log']))
                  for name, sp in search_space.items()}
        if 'threshold_conc' in values:  # current threshold-anchored scheme
            threshold = values['threshold_conc']
            target = min(TARGET_CONC_MAX,
                         threshold + values['target_delta'])
            spike = min(SPIKE_CONC_MAX,
                        max(SPIKE_CONC_MIN,
                            target + values['spike_delta']))
        else:  # legacy target-anchored scheme
            target = values['target_conc']
            threshold = max(0.0, target - values['threshold_delta'])
            spike = values['spike_conc']
        model_kwargs = dict(target_conc=target, threshold_conc=threshold,
                            spike_conc=spike)
        record = {'trial_number': trial.number, **values}
        try:
            for pname in kinetic_baselines:
                if pname in values:
                    setattr(r_te, pname, values[pname])
            if 'max_n_spikes' in values:
                fbs_spec.max_n_spikes = values['max_n_spikes']
            handles['model_specification'](**model_kwargs)
            handles['latest_TEA_solution'].update(
                handles['solve_TEA'](
                    stream_IDs=('ethanol', 'isobutanol')))
            for mname, getter in TRACKED_METRICS.items():
                record[mname] = getter(handles)
            obj = float(objective_getter(handles))
            record['objective'] = obj
        except Exception as e:
            record['state'] = 'FAIL'
            record['error'] = repr(e)[:300]
            append_trajectory_row(csv_path, columns, record)
            print(f'Trial {trial.number}: FAILED ({repr(e)[:120]})')
            raise optuna.TrialPruned() from e
        if not math.isfinite(obj):
            # NaN (not solved) or +-inf (solve_TEA reports an IRR with no
            # real root on the valid domain as -inf since 2026-09-03; the
            # sampler needs finite values, so such trials are pruned like
            # NaN ones -- the CSV keeps the raw value in 'objective')
            record['state'] = 'NAN'
            append_trajectory_row(csv_path, columns, record)
            print(f'Trial {trial.number}: objective is non-finite ({obj}); pruned.')
            raise optuna.TrialPruned()
        record['state'] = 'COMPLETE'
        append_trajectory_row(csv_path, columns, record)
        for mname in TRACKED_METRICS:
            trial.set_user_attr(mname, record[mname])
        if trial.number % print_status_every == 0:
            try:
                best = study.best_value
            except Exception:  # no completed trial stored yet
                best = np.nan
            try:
                print(f'\nTrial {trial.number}/{n_trials}: '
                      f'{objective_name} = {obj:.6g} '
                      f'(best so far {best:.6g})\n'
                      f'integrator: {r_te.integrator.getName()}; '
                      'HXN Qbal error = '
                      f"{handles['HXN'].energy_balance_percent_error:.2f} %")
            except Exception:  # cosmetic only -- never abort the study
                pass
        return obj

    n_remaining = max(0, n_trials - n_done)
    if n_done:
        print(f'Resuming study {study_name}: {n_done} trials stored; '
              f'running {n_remaining} more (budget {n_trials}).')
    try:
        study.optimize(_objective, n_trials=n_remaining, n_jobs=1,
                       gc_after_trial=True)
    finally:
        restore_baseline(handles, kinetic_baselines,
                         baseline_model_kwargs,
                         baseline_max_n_spikes=baseline_max_n_spikes)
    return study, csv_path, kinetic_baselines

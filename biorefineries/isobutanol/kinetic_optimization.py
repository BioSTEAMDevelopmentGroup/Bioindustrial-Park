#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Bayesian (Optuna TPE) global optimization of all fermentation kinetic
parameters (k_*/K_* on V406's tellurium model) plus the feeding-strategy
sugar concentrations (target_conc, threshold_conc via threshold_delta,
spike_conc), against a named or custom objective -- to prioritize
metabolic-engineering / bioprocess research directions.

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

import numpy as np

__all__ = ('OBJECTIVE_REGISTRY', 'TRACKED_METRICS',
           'discover_kinetic_parameters', 'build_search_space',
           'trajectory_columns', 'append_trajectory_row', 'load_trajectory',
           'get_handles', 'run_kinetic_optimization', 'restore_baseline',
           'plot_optimization_trajectories', 'plot_parameter_trajectory',
           'plot_best_vs_baseline')

FEEDING_VARIABLES = ('target_conc', 'threshold_delta', 'spike_conc')

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
                       target_conc_bounds=(180.0, 300.0),
                       threshold_delta_bounds=(0.5, 30.0),
                       spike_conc_bounds=(200.0, 800.0),
                       ):
    """Build the decision-variable space: {name: {'low', 'high', 'log'}}.

    Kinetic parameters default to [m_lo*baseline, m_hi*baseline] sampled
    log-scale (research-opportunity width); entries in
    `param_bounds_override` ({name: (low, high)}) use those absolute
    bounds instead (log-scale only if low > 0). A parameter with a
    nonpositive baseline and no override cannot use the multiplier band
    and is EXCLUDED with a printed warning. `exclude_params` names are
    always excluded (silently). The three feeding variables are appended
    with absolute linear bounds; threshold_conc is optimized as
    threshold_delta = target_conc - threshold_conc, which guarantees
    threshold < target by construction.

    Returns (space, excluded_parameter_names)."""
    param_bounds_override = dict(param_bounds_override or {})
    m_lo, m_hi = multiplier_bounds
    space, excluded = {}, []
    for name, baseline in kinetic_baselines.items():
        if name in exclude_params:
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
    space['target_conc'] = dict(low=target_conc_bounds[0],
                                high=target_conc_bounds[1], log=False)
    space['threshold_delta'] = dict(low=threshold_delta_bounds[0],
                                    high=threshold_delta_bounds[1],
                                    log=False)
    space['spike_conc'] = dict(low=spike_conc_bounds[0],
                               high=spike_conc_bounds[1], log=False)
    return space, excluded

#%% Trajectory recording

def trajectory_columns(search_space):
    """Column order of the trajectory CSV for a given search space."""
    return ['trial_number', 'state', *search_space.keys(), 'objective',
            *TRACKED_METRICS.keys(), 'error']

def append_trajectory_row(csv_path, columns, record):
    """Append one row (dict; missing keys become '') to `csv_path`,
    writing the header first if the file does not exist. Flushed
    immediately, so a crash/segfault loses at most the in-flight trial."""
    exists = os.path.isfile(csv_path)
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
    number. Failed/pruned trials are omitted. Returns (fig, axes)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
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
        ax.set_title(m)
        ax.set_xlabel('trial')
    for ax in flat[n_panels:]:
        ax.axis('off')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, axes

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
    if 'target_conc' in ok.columns:
        ax2.plot(x, best_rows['target_conc'], label='target_conc')
        ax2.plot(x, best_rows['target_conc'] - best_rows['threshold_delta'],
                 label='threshold_conc')
        ax2.plot(x, best_rows['spike_conc'], label='spike_conc')
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
    ax.set_title(f"objective = {best['objective']:.4g} at trial "
                 f"{int(best['trial_number'])}; target = "
                 f"{best['target_conc']:.1f}, threshold = "
                 f"{best['target_conc'] - best['threshold_delta']:.1f}, "
                 f"spike = {best['spike_conc']:.1f} g/L", fontsize=8)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, ax

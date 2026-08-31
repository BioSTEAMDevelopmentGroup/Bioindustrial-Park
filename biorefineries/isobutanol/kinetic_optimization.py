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
variables (target_conc, threshold_conc via threshold_delta, spike_conc,
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

import numpy as np

__all__ = ('OBJECTIVE_REGISTRY', 'TRACKED_METRICS',
           'discover_kinetic_parameters', 'build_search_space',
           'trajectory_columns', 'append_trajectory_row', 'load_trajectory',
           'get_handles', 'run_kinetic_optimization', 'restore_baseline',
           'plot_optimization_trajectories', 'plot_parameter_trajectory',
           'plot_best_vs_baseline')

FEEDING_VARIABLES = ('target_conc', 'threshold_delta', 'spike_conc',
                     'max_n_spikes')

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
                       target_conc_bounds=(5.0, 500.0),
                       threshold_delta_bounds=(0.5, 500.0),
                       spike_conc_bounds=(50.0, 600.0),
                       max_n_spikes_bounds=(0, 20),
                       ):
    """Build the decision-variable space: {name: {'low', 'high', 'log'}}
    (integer variables additionally carry 'int': True).

    Kinetic parameters default to [m_lo*baseline, m_hi*baseline] sampled
    log-scale (research-opportunity width); entries in
    `param_bounds_override` ({name: (low, high)}) use those absolute
    bounds instead (log-scale only if low > 0). A parameter with a
    nonpositive baseline and no override cannot use the multiplier band
    and is EXCLUDED with a printed warning. `exclude_params` names are
    always excluded (silently). The feeding variables are appended with
    absolute linear bounds; threshold_conc is optimized as
    threshold_delta = target_conc - threshold_conc, which guarantees
    threshold < target by construction (the applied threshold is floored
    at 0.0, so wide deltas map to threshold_conc = 0 -- feeding never
    triggered). Note the model itself requires target < spike; draws
    violating it fail the specification check pre-simulation and are
    pruned. max_n_spikes (the glucose-spike cap, fbs_spec.max_n_spikes)
    is an INTEGER variable (0 = forced batch); pass
    max_n_spikes_bounds=None to pin it at the scenario baseline instead
    (pre-2026-08-31 search-space behavior).

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
    if max_n_spikes_bounds is not None:
        space['max_n_spikes'] = dict(low=int(max_n_spikes_bounds[0]),
                                     high=int(max_n_spikes_bounds[1]),
                                     log=False, int=True)
    return space, excluded

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
        ax2.plot(x, (best_rows['target_conc']
                     - best_rows['threshold_delta']).clip(lower=0.0),
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
                             target_conc_bounds=(5.0, 500.0),
                             threshold_delta_bounds=(0.5, 500.0),
                             spike_conc_bounds=(50.0, 600.0),
                             max_n_spikes_bounds=(0, 20),
                             study_name=None, results_dir=None,
                             handles=None, print_status_every=1,
                             ):
    """Run the Bayesian optimization. `objective` is a name in
    OBJECTIVE_REGISTRY (direction/level/units filled from the entry) or a
    custom callable(handles)->float (then `direction`, and ideally
    `objective_name`/`level`/`objective_units`, must be given).

    `n_trials` is the TOTAL budget of the study: rerunning with the same
    study_name resumes from the on-disk SQLite store and runs only the
    remainder (crash/segfault recovery). Trials execute STRICTLY
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
        target_conc_bounds=target_conc_bounds,
        threshold_delta_bounds=threshold_delta_bounds,
        spike_conc_bounds=spike_conc_bounds,
        max_n_spikes_bounds=max_n_spikes_bounds)
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

    def _objective(trial):
        values = {name: (trial.suggest_int(name, sp['low'], sp['high'])
                         if sp.get('int') else
                         trial.suggest_float(name, sp['low'], sp['high'],
                                             log=sp['log']))
                  for name, sp in search_space.items()}
        model_kwargs = dict(
            target_conc=values['target_conc'],
            threshold_conc=max(
                0.0, values['target_conc']-values['threshold_delta']),
            spike_conc=values['spike_conc'])
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
        if math.isnan(obj):
            record['state'] = 'NAN'
            append_trajectory_row(csv_path, columns, record)
            print(f'Trial {trial.number}: objective is NaN; pruned.')
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

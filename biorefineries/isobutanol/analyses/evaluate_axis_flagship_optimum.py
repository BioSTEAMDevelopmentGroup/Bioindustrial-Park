#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""1-D sweeps of any metabolic_split_12d decision variable through the
flagship profitability campaign's optimum.

The generalization of evaluate_k13_flagship_optimum.py: every decision
variable is held at the flagship (relay) campaign's best trial
(`..._pi_log-tail_gp_rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden` #1912: PI
0.808, IRR 0.273, 1 glucose spike) EXCEPT the swept one(s). Two run modes,
chosen by IBO_AXIS_SWEEP (or --sweep; supervise_sweep.py passes no arguments,
so the environment variable is the way to configure a supervised run):

* 'screen' (default) -- a COARSE screen of SCREEN_PARAMS (the 11 variables
  other than k_13, already swept at full resolution): N points per variable
  over its campaign band (log-spaced on a log axis, linear on a linear one,
  every integer on the max_n_spikes axis when there are at most 2N of them),
  plus the trial's own value. Output results/<stem>_screen.csv (one row per
  point, 'param' / 'value' columns). Feeds the slice ranking in
  plots/plot_axis_flagship_sweep.py --screen: which slice makes PI spiky
  while the isobutanol yield stays smooth.
* a variable name, e.g. 'threshold_conc' -- a FULL-resolution sweep of that
  variable alone (IBO_AXIS_N_POINTS, default 160) over its campaign band, or
  over [IBO_AXIS_LOW, IBO_AXIS_HIGH] when given. Output
  results/<stem>_<name>.csv.

A grouped variable (glycolysis, ehrlich_downstream, inhib_*) is re-expanded
into its members with ko.expand_grouped_values (the campaign's own group
references), so the members reaching the model are exactly what the campaign
would have applied; a feeding variable (threshold_conc, target_delta,
max_n_spikes) goes through _resolve_feeding_concs / the spike cap inside
ko._simulate_trial_reproduction (target_conc = min(300, threshold +
target_delta), spike pinned at 600 g/L).

Set-up, anchor check and per-point protocol are those of the k_13 sweep: the
read-only trial reproduction (scenario A, A-referenced burden ON), the trial
simulated first and compared with its recorded row (abort beyond
ANCHOR_PI_TOL), its CONVERGED flowsheet snapshotted, and EVERY point restored
to that snapshot before simulating, so the sweep order can neither create nor
hide a spike. An over-cap point is logged INFEASIBLE, a raising simulation
(incl. the fed-batch volume guard) ERROR.

SIMULATES: ask-first (the screen approved by the user 2026-09-24). Run it in
a fresh process, from this directory, with the hensmith pin used for the
rs350 / flagship campaigns (PYTHONPATH from the session scratchpad's
PINNED_PYTHONPATH.txt), under the supervisor with the matching stem:

    $env:IBO_AXIS_SWEEP = 'screen'
    python supervise_sweep.py evaluate_axis_flagship_optimum.py \
        --stem evaluate_axis_flagship_optimum_screen

The k_3 figure's data (plots/plot_axis_flagship_sweep.py) is the full campaign
band 0.00581-23.24 at 218 log-spaced points (+ the trial's own value; the
density of the earlier 160-point 0.0058-2.5 zoom), approved 2026-09-24:

    $env:IBO_AXIS_SWEEP = 'k_3'; $env:IBO_AXIS_N_POINTS = '218'
    python supervise_sweep.py evaluate_axis_flagship_optimum.py \
        --stem evaluate_axis_flagship_optimum_k_3

Checkpoint + resume (the supervise_sweep.py convention, per run tag):
results/<stem>_<tag>_checkpoint.csv (one flushed row per point) and
results/<stem>_<tag>_inflight.json (written right before each point); a
relaunch logs the point the previous process died in as LOST and resumes
after it (the anchor is re-simulated each attempt to rebuild the snapshot);
the checkpoint's (param, value) pairs are verified against the grid;
IBO_SWEEP_FRESH=1 discards a leftover one. On completion the rows are written
in grid order to results/<stem>_<tag>.csv (plus <stem>_<tag>_anchor.json: the
reproduction check, grid and environment) and the checkpoint is deleted."""
import os
import csv
import json
import math
import time
import argparse

import numpy as np

#%% Settings
STUDY_NAME = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_pi_log-tail_gp_'
              'rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden')
TRIAL_NUMBER = 1912
ANCHOR_SCENARIO = 'A'
ANCHOR_PI_TOL = 0.02    # relative; a larger anchor-PI miss aborts the sweep
SCREEN_PARAMS = ('threshold_conc', 'target_delta', 'max_n_spikes',
                 'glycolysis', 'k_3', 'k_6', 'k_17', 'ehrlich_downstream',
                 'inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate')
SCREEN_N_POINTS = 30
FULL_N_POINTS = 160

_here = os.path.dirname(os.path.abspath(__file__))
_stem = os.path.splitext(os.path.basename(__file__))[0]
RESULTS_DIR = os.path.join(_here, 'results')

def run_paths(tag):
    base = os.path.join(RESULTS_DIR, f'{_stem}_{tag}')
    return dict(checkpoint=base + '_checkpoint.csv',
                inflight=base + '_inflight.json',
                output=base + '.csv', anchor=base + '_anchor.json')

def _env_float(name):
    value = os.environ.get(name, '')
    return float(value) if value else None

#%% Grid
def _round(x):
    return float(f'{x:.12g}')

def axis_grid(entry, n_points, anchor_value, low=None, high=None):
    """Sorted grid over [low, high] (default: the campaign band `entry`) plus
    the trial's own value; log-spaced on a log axis, linear on a linear one;
    an 'int' axis takes every integer when there are at most 2*n_points of
    them, else rounded linear points. Values rounded to 12 significant
    digits so a checkpoint written by an earlier process matches exactly."""
    low = entry['low'] if low is None else low
    high = entry['high'] if high is None else high
    if entry.get('int'):
        lo, hi = int(math.ceil(low)), int(math.floor(high))
        if hi - lo + 1 <= 2*n_points:
            grid = np.arange(lo, hi + 1)
        else:
            grid = np.round(np.linspace(lo, hi, n_points))
        return sorted({int(x) for x in grid} | {int(anchor_value)})
    grid = (np.geomspace(low, high, n_points) if entry['log']
            else np.linspace(low, high, n_points))
    return sorted({_round(x) for x in grid} | {_round(anchor_value)})

#%% Checkpoint helpers
def fieldnames(metric_names, group_names):
    return (['i', 'param', 'value', 'state', 'is_anchor'] + list(metric_names)
            + ['MPSP ethanol', 'MPSP isobutanol', 'threshold', 'target',
               'Phi_M', 'phi_T', 'burden_factor']
            + [installed_column(g) for g in group_names]
            + ['wall_s', 'error'])

def installed_column(group_name):
    """Column of a unit group's installed equipment cost (MM$), so a PI step
    can be traced to the part of the plant whose capital moved."""
    return f'installed MM$ {group_name}'

def _same_value(cell, value):
    return float(cell) == float(value)

def load_checkpoint(path, names, grid):
    done = {}
    if not os.path.exists(path):
        return done
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != names:
            raise RuntimeError(f'checkpoint columns do not match this sweep: '
                               f'{path} (delete it or set IBO_SWEEP_FRESH=1)')
        for row in reader:
            i = int(row['i'])
            if (i >= len(grid) or row['param'] != grid[i][0]
                    or not _same_value(row['value'], grid[i][1])):
                raise RuntimeError(f'checkpoint is of a different grid: {path} '
                                   '(delete it or set IBO_SWEEP_FRESH=1)')
            done[i] = row
    return done

def append_checkpoint(path, names, row):
    is_new = not os.path.exists(path)
    with open(path, 'a', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        if is_new:
            writer.writeheader()
        writer.writerow(row)
        fh.flush()
        os.fsync(fh.fileno())

def write_inflight(path, i, param, value):
    with open(path, 'w') as fh:
        json.dump(dict(i=i, param=param, value=value, started=time.time()), fh)

def clear_inflight(path):
    if os.path.exists(path):
        os.remove(path)

#%% Runner
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--sweep', default=os.environ.get('IBO_AXIS_SWEEP',
                                                          'screen'),
                        help="'screen' or one decision variable name")
    parser.add_argument('--n-points', type=int,
                        default=int(os.environ.get('IBO_AXIS_N_POINTS', 0)) or None)
    parser.add_argument('--low', type=float, default=_env_float('IBO_AXIS_LOW'))
    parser.add_argument('--high', type=float, default=_env_float('IBO_AXIS_HIGH'))
    parser.add_argument('--study-name', default=STUDY_NAME)
    parser.add_argument('--trial-number', type=int, default=TRIAL_NUMBER)
    args = parser.parse_args(argv)
    screen = args.sweep == 'screen'
    if screen and (args.low is not None or args.high is not None):
        raise ValueError('--low / --high apply to a single-variable sweep only')
    tag = args.sweep
    paths = run_paths(tag)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    import biorefineries.isobutanol as isobutanol
    isobutanol.load()   # default both-trains build (S201 split 1.0), as the driver runs
    from biorefineries.isobutanol import kinetic_optimization as ko
    from biorefineries.isobutanol import scenarios, system
    import hensmith

    # --- the trial's decision point (reproduce_split12d_trial's set-up) ---
    csv_path = ko.split12d_trajectory_path(args.study_name)
    row = ko.read_trajectory_row(csv_path, args.trial_number)
    preset = ko.resolve_study_preset(ko.SPLIT12D_STUDY_TARGET_PRODUCTS,
                                     ko.SPLIT12D_STUDY_TYPE)
    scenarios.load_scenario(ANCHOR_SCENARIO, burden=True)
    handles = ko.get_handles()
    r_te, fbs_spec = handles['r_te'], handles['fbs_spec']
    kinetic_baselines = ko.discover_kinetic_parameters(r_te)
    baseline_model_kwargs = {k: fbs_spec.current_specifications[k]
                             for k in ('target_conc', 'threshold_conc',
                                       'spike_conc')}
    search_space, _ = ko.build_search_space(
        kinetic_baselines,
        **{key: value for key, value in preset.items()
           if key not in ('scenario', 'kinetic_bounds_scenario')})
    groups = preset['parameter_groups']
    references = preset['group_references']
    values, applied, cross_check = ko.reconstruct_trial_kinetics(
        row, search_space, groups, kinetic_baselines, references, mode='both')
    burden_model = system.get_active_burden()
    if burden_model is None:
        raise RuntimeError('no active enzyme burden after load_scenario')

    params = SCREEN_PARAMS if screen else (args.sweep,)
    unknown = [p for p in params if p not in search_space]
    if unknown:
        raise KeyError(f'not decision variables of this space: {unknown}; '
                       f'choose from {list(search_space)}')
    n_points = args.n_points or (SCREEN_N_POINTS if screen else FULL_N_POINTS)
    grid, anchor_idx, bands = [], set(), {}
    for p in params:
        entry = search_space[p]
        axis = axis_grid(entry, n_points, values[p], args.low, args.high)
        bands[p] = dict(low=entry['low'], high=entry['high'],
                        log=bool(entry['log']), int=bool(entry.get('int')),
                        swept_low=axis[0], swept_high=axis[-1],
                        anchor=values[p], n=len(axis))
        for x in axis:
            if x == (int(values[p]) if entry.get('int') else _round(values[p])):
                anchor_idx.add(len(grid))
            grid.append((p, x))

    metric_names = list(ko.TRACKED_METRICS)
    unit_groups = list(system.unit_groups)
    names = fieldnames(metric_names, [g.name for g in unit_groups])

    def installed_costs():
        """{column: installed equipment cost (MM$)} of every unit group
        (reporting getters only; no state is touched)."""
        out = {}
        for group in unit_groups:
            metric, = [m for m in group.metrics
                       if m.name == 'Installed equipment cost']
            out[installed_column(group.name)] = float(metric())
        return out

    def simulate(param, x):
        """Simulate the trial point with `param` = x; returns a checkpoint
        row (without 'i')."""
        vals = dict(values, **{param: x})
        app = ko.expand_grouped_values(vals, groups, kinetic_baselines,
                                       references)
        t0 = time.time()
        feeding, reproduced, error = ko._simulate_trial_reproduction(
            handles, kinetic_baselines, vals, app, baseline_model_kwargs)
        capacities = {name: float(getattr(r_te, name))
                      for name in burden_model.required_capacities()}
        burden = burden_model.evaluate(capacities)
        if error is None:
            state = 'OK'
        elif 'BurdenInfeasible' in error:
            state = 'INFEASIBLE'
        else:
            state = 'ERROR'
        out = dict(param=param, value=x, state=state,
                   threshold=feeding['threshold'], target=feeding['target'],
                   Phi_M=burden.Phi_M, phi_T=burden.phi_T,
                   burden_factor=burden.burden_factor,
                   wall_s=round(time.time() - t0, 2), error=error or '')
        metrics = reproduced['metrics'] if error is None else {}
        for name in metric_names:
            out[name] = float(metrics[name]) if name in metrics else math.nan
        mpsps = reproduced['MPSPs'] if error is None else None
        out['MPSP ethanol'] = mpsps['ethanol'] if mpsps else math.nan
        out['MPSP isobutanol'] = mpsps['isobutanol'] if mpsps else math.nan
        costs = installed_costs() if error is None else {}
        for group in unit_groups:
            column = installed_column(group.name)
            out[column] = costs.get(column, math.nan)
        return out, feeding, reproduced

    # --- the anchor: reproduce the trial, then snapshot its converged state ---
    print(f'\nAnchor: {args.study_name} #{args.trial_number}; sweep {tag!r} '
          f'({len(grid)} points over {len(params)} variable(s)); hensmith from '
          f'{hensmith.__file__}')
    anchor, feeding, reproduced = simulate(params[0], values[params[0]])
    if anchor['state'] != 'OK':
        raise RuntimeError(f"anchor simulation failed: {anchor['error']}")
    metric_check, metric_warnings = ko.compare_tracked_metrics(
        reproduced['metrics'], row)
    recorded_PI = float(row['PI'])
    rel_PI = abs(anchor['PI'] - recorded_PI)/abs(recorded_PI)
    print(f"  PI reproduced {anchor['PI']:.5f} vs recorded {recorded_PI:.5f} "
          f"(rel {rel_PI:.2e}); flagged metrics: {metric_warnings or 'none'}")
    if rel_PI > ANCHOR_PI_TOL:
        raise RuntimeError(f'anchor PI off by {rel_PI:.3g} (> {ANCHOR_PI_TOL}); '
                           'the sweep would not pass through the optimum')
    snapshot = handles['snapshot_flowsheet_state']()
    with open(paths['anchor'], 'w') as fh:
        json.dump(dict(
            study_name=args.study_name, trial_number=args.trial_number,
            anchor_scenario=ANCHOR_SCENARIO, burden=True, sweep=tag,
            values=values, feeding=feeding, cross_check=dict(
                max_rel_delta=cross_check['max_rel_delta'],
                n_mismatches=len(cross_check['mismatches'])),
            reproduced_PI=anchor['PI'], recorded_PI=recorded_PI,
            reproduced_IBO_yield=anchor['IBO yield'],
            metric_check=[list(m) for m in metric_check],
            metric_warnings=metric_warnings, bands=bands,
            n_points=n_points, n_total=len(grid),
            hensmith=hensmith.__file__), fh, indent=1, default=float)

    # --- checkpoint / resume ---
    if os.environ.get('IBO_SWEEP_FRESH', '') == '1':
        for path in (paths['checkpoint'], paths['inflight']):
            if os.path.exists(path):
                os.remove(path)
    done = load_checkpoint(paths['checkpoint'], names, grid)
    if os.path.exists(paths['inflight']):
        with open(paths['inflight']) as fh:
            lost = json.load(fh)
        if lost['i'] not in done:
            print(f"  logging point {lost['i']} ({lost['param']} = "
                  f"{lost['value']:g}) LOST: the previous process died in it")
            row_lost = {name: math.nan for name in names}
            row_lost.update(i=lost['i'], param=grid[lost['i']][0],
                            value=grid[lost['i']][1], state='LOST',
                            is_anchor=int(lost['i'] in anchor_idx), error='')
            append_checkpoint(paths['checkpoint'], names, row_lost)
            done[lost['i']] = row_lost
        clear_inflight(paths['inflight'])
    if done:
        print(f'  RESUMING: {len(done)} of {len(grid)} points already logged '
              '(IBO_SWEEP_FRESH=1 for a fresh sweep)')

    # --- the sweep: every point from the anchor's converged state ---
    t_start = time.time()
    todo = [i for i in range(len(grid)) if i not in done]
    for n, i in enumerate(todo, 1):
        param, x = grid[i]
        handles['restore_flowsheet_state'](snapshot)
        write_inflight(paths['inflight'], i, param, x)
        out, _, _ = simulate(param, x)
        out['i'] = i
        out['is_anchor'] = int(i in anchor_idx)
        append_checkpoint(paths['checkpoint'], names, out)
        clear_inflight(paths['inflight'])
        done[i] = out
        eta = (time.time() - t_start)/n*(len(todo) - n)/60
        print(f"  [{len(done):3d}/{len(grid)}] {param:<18} {x:9.4g}  "
              f"{out['state']:<10} PI {out['PI']:9.4f}  IBO yield "
              f"{out['IBO yield']:.4f}  spikes {out['n_glu_spikes']:4.0f}  "
              f"tau {out['tau']:6.2f}  ({out['wall_s']:.1f} s; "
              f"~{eta:.0f} min left)", flush=True)

    # --- final output, in grid order; the checkpoint is then removed ---
    with open(paths['output'], 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        writer.writeheader()
        for i in sorted(done):
            writer.writerow({name: done[i][name] for name in names})
    os.remove(paths['checkpoint'])
    states = [str(done[i]['state']) for i in done]
    print(f"\nWrote {paths['output']} ({len(done)} points: "
          + ', '.join(f'{s} {states.count(s)}' for s in sorted(set(states)))
          + f") and {paths['anchor']}")

if __name__ == '__main__':
    main()

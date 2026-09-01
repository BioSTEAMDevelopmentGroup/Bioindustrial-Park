#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Crash- AND stall-resilient supervisor for the Bayesian kinetic
optimization (analyses/optimize_kinetics_BO.py). Runs the driver in child
processes and, until the study's total trial budget completes cleanly:

- relaunches after a crash (the known nondeterministic native CVODE
  segfault) exactly like the proven bash-loop recipe;
- kills and relaunches a STALLED child -- the per-trial wall-clock
  timeout: a pathological kinetic draw can hang inside a native
  CVODE/roadrunner call, which nothing in-process can interrupt on
  Windows, so the timeout is enforced by polling the crash-safe
  trajectory CSV and killing the child when no trial has been recorded
  for --stall-timeout-min (default 25 min, comfortably above the ~15 min
  worst self-resolving trial observed in production). The in-flight
  trial is lost, as in a segfault; its RUNNING optuna record still
  counts toward the budget;
- aborts (exit 1) after any attempt that recorded no new trials --
  resuming would loop on the same failure.

Decision logic (StallGuard, attempt_outcome) lives in
kinetic_optimization.py and is covered by the offline test. This
supervisor is stdlib-only and loads that module BY FILE PATH -- it must
never import the biorefineries.isobutanol package itself (the package
__init__ pulls biosteam -> eager numba compilation, which corrupts the
shared numba cache when concurrent with the child's own load()).

Long-running, one simulation at a time (children run strictly
sequentially) -- ask-first, like the unsupervised driver. Example:

    python optimize_kinetics_BO_supervised.py --scenario A \\
        --kinetic-bounds-scenario B --objective IRR --n-trials 2000
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import time

ANALYSES_DIR = os.path.dirname(os.path.abspath(__file__))
IBO_DIR = os.path.dirname(ANALYSES_DIR)
DRIVER = os.path.join(ANALYSES_DIR, 'optimize_kinetics_BO.py')
RESULTS_DIR = os.path.join(ANALYSES_DIR, 'results')

_spec = importlib.util.spec_from_file_location(
    '_ko_supervision', os.path.join(IBO_DIR, 'kinetic_optimization.py'))
ko = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ko)


def default_study_name(scenario, objective, kinetic_bounds_scenario):
    """Mirror the driver's stable study naming (resume finds the same
    study)."""
    slug = objective.lower().replace(' ', '_')
    if kinetic_bounds_scenario:
        return f'kin_opt_{scenario}_kb{kinetic_bounds_scenario}_{slug}'
    return f'kin_opt_{scenario}_{slug}'


def row_count(csv_path):
    """Data rows in the trajectory CSV (0 when absent/empty)."""
    try:
        with open(csv_path, 'rb') as f:
            return max(0, sum(1 for _ in f) - 1)
    except OSError:
        return 0


def child_code(scenario, objective, n_trials, kinetic_bounds_scenario,
               make_plots, study_name):
    """The -c program for one supervised attempt of the driver."""
    return (
        'import runpy\n'
        f'ns = runpy.run_path({DRIVER!r})\n'
        f"ns['run'](scenario={scenario!r}, objective={objective!r},\n"
        f'          n_trials={n_trials!r},\n'
        f'          kinetic_bounds_scenario={kinetic_bounds_scenario!r},\n'
        f'          make_plots={make_plots!r},\n'
        f'          study_name={study_name!r})\n')


def supervise(scenario='B', objective='IRR', n_trials=2000,
              kinetic_bounds_scenario=None, make_plots=True,
              study_name=None, stall_timeout_min=25.0, poll_s=30.0,
              settle_s=10.0, python=None, log_path=None):
    """Run attempts until 'complete' or 'abort'; returns the final
    outcome string ('complete' or 'abort')."""
    if study_name is None:
        study_name = default_study_name(scenario, objective,
                                        kinetic_bounds_scenario)
    csv_path = os.path.join(RESULTS_DIR, study_name + '_trajectory.csv')
    if python is None:
        python = sys.executable
    if log_path is None:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        log_path = os.path.join(RESULTS_DIR, study_name + '_run.log')
    code = child_code(scenario, objective, n_trials,
                      kinetic_bounds_scenario, make_plots, study_name)
    guard = ko.StallGuard(stall_timeout_s=60.0*stall_timeout_min)

    def event(msg):
        line = f'SUPERVISOR: {msg}'
        print(line, flush=True)
        with open(log_path, 'a') as f:
            f.write(line + '\n')

    attempt = 0
    while True:
        attempt += 1
        rows_before = row_count(csv_path)
        guard.reset()
        event(f'attempt {attempt} starting ({time.strftime("%c")}); '
              f'{rows_before} trials logged so far; log: {log_path}')
        with open(log_path, 'a') as log:
            child = subprocess.Popen([python, '-u', '-c', code],
                                     stdout=log, stderr=log)
            killed = False
            while child.poll() is None:
                time.sleep(poll_s)
                if guard.update(row_count(csv_path),
                                time.monotonic()) == 'stalled':
                    event(f'attempt {attempt}: no new trial in '
                          f'{stall_timeout_min:g} min -- killing stalled '
                          'child (per-trial timeout; in-flight trial '
                          'lost, study resumes)')
                    child.kill()
                    child.wait()
                    killed = True
                    break
        rows_after = row_count(csv_path)
        outcome = ko.attempt_outcome(child.returncode, rows_before,
                                     rows_after, killed_for_stall=killed)
        if outcome == 'complete':
            event(f'clean exit after attempt {attempt} '
                  f'({rows_after} trials logged)')
            return outcome
        if outcome == 'abort':
            event(f'attempt {attempt} (exit {child.returncode}, '
                  f'stall-killed={killed}) recorded no new trials '
                  f'(rows={rows_after}); aborting -- investigate before '
                  'relaunching')
            return outcome
        event(f'attempt {attempt} ended (exit {child.returncode}, '
              f'stall-killed={killed}) at {rows_after} rows; resuming '
              f'in {settle_s:g} s')
        time.sleep(settle_s)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Supervised (crash- and stall-resilient) kinetic BO '
                    'run; see module docstring.')
    parser.add_argument('--scenario', default='B', choices=('A', 'B'))
    parser.add_argument('--objective', default='IRR',
                        help='name in OBJECTIVE_REGISTRY (callables: use '
                             'the unsupervised driver)')
    parser.add_argument('--n-trials', type=int, default=2000)
    parser.add_argument('--kinetic-bounds-scenario', default=None,
                        choices=('A', 'B'),
                        help="derive kinetic bounds from this scenario's "
                             'workbook (e.g. B for a scenario-A run '
                             'without zero-baseline exclusions)')
    parser.add_argument('--study-name', default=None,
                        help='override the derived stable study name')
    parser.add_argument('--no-plots', action='store_true',
                        help='skip the final driver plots')
    parser.add_argument('--stall-timeout-min', type=float, default=25.0)
    parser.add_argument('--poll-s', type=float, default=30.0)
    parser.add_argument('--settle-s', type=float, default=10.0)
    args = parser.parse_args()
    outcome = supervise(scenario=args.scenario, objective=args.objective,
                        n_trials=args.n_trials,
                        kinetic_bounds_scenario=args.kinetic_bounds_scenario,
                        make_plots=not args.no_plots,
                        study_name=args.study_name,
                        stall_timeout_min=args.stall_timeout_min,
                        poll_s=args.poll_s, settle_s=args.settle_s)
    sys.exit(0 if outcome == 'complete' else 1)

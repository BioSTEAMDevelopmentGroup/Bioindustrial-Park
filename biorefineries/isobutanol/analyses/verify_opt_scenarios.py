#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Verify that ONE opt_* scenario's baseline simulation matches the
reproduction of its campaign's optimal trial (SIMULATES; read-only: writes
no workbook, CSV, sidecar or store).

In one FRESH kernel, on the smoke tests' IBO_EtOH-only build (== both
trains at S201 split 1.0, the campaigns' build):
  1. scenarios.load_scenario(<scenario>) -- the scenario's baseline
     simulation from its workbook + registry feeding strategy, exactly as
     its smoke test runs it -- then solve_TEA and every TRACKED_METRICS
     getter;
  2. resets the kinetics to scenario A's full live kinetics (the A workbook
     alone cannot undo the opt_* k_13-k_17 / isobutanol-inhibition values),
     then kinetic_optimization.reproduce_split12d_trial on the campaign /
     trial build_opt_split12d_workbooks.OPT_SCENARIOS maps the scenario to
     (anchor A, A-referenced burden on, mode 'both');
  3. compares the two MPSPs and every tracked metric (convergence
     diagnostics excepted) at a relative tolerance TOL; the reproduction's
     metric-by-metric comparison with the RECORDED trial is printed too.
Exit 0 + "MATCH" = the scenario's baseline IS the reproduced trial.

Run (one fresh process per scenario; several in parallel only on a warm
numba cache):
  & "C:/Users/saran/anaconda3/envs/IBO_2026/python.exe" analyses/verify_opt_scenarios.py opt_PI_TRY_informed
"""
import os
import sys
import argparse
import importlib.util

#: Max relative delta (MPSPs + tracked metrics) scenario vs reproduction.
TOL = 1e-3

_ANALYSES_DIR = os.path.dirname(os.path.abspath(__file__))


def _load_generator():
    """build_opt_split12d_workbooks by file path (build-free import)."""
    path = os.path.join(_ANALYSES_DIR, 'build_opt_split12d_workbooks.py')
    spec = importlib.util.spec_from_file_location(
        'build_opt_split12d_workbooks', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None):
    gen = _load_generator()
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('scenario', choices=sorted(gen.OPT_SCENARIOS))
    parser.add_argument('--tol', type=float, default=TOL)
    args = parser.parse_args(argv)
    objective_name, study_name, trial_number = gen.OPT_SCENARIOS[args.scenario]

    from biorefineries import isobutanol
    isobutanol.load(separation_processes=('IBO_EtOH',))
    from biorefineries.isobutanol import scenarios
    from biorefineries.isobutanol import kinetic_optimization as ko

    # 1. The scenario's baseline simulation (load_scenario simulates once).
    bundle = scenarios.load_scenario(args.scenario)
    handles = ko.get_handles()
    solution = handles['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'),
                                    IRR_for_MPSP=0.15)
    handles['latest_TEA_solution'].update(solution)
    scenario_values = gen.comparable_values(
        solution['MPSPs'],
        {name: getter(handles) for name, getter in ko.TRACKED_METRICS.items()})
    feeding = bundle['feeding_kwargs']
    n_spikes = bundle['spec'].max_n_spikes

    # 2. The campaign trial's reproduction (restores the anchor kinetics).
    #    The scenario-A workbook has no k_13-k_17 / isobutanol-inhibition
    #    rows, so its load_scenario('A') alone would keep the opt_* values
    #    as the anchor; reset every kinetic parameter to scenario A's full
    #    live kinetics first (load_scenario's per-kernel A reference, taken
    #    on the clean post-load() model before the opt_* workbook).
    r_te = handles['r_te']
    for name, value in scenarios._A_REFERENCE[0].items():
        setattr(r_te, name, value)
    result = ko.reproduce_split12d_trial('A', study_name, trial_number,
                                         mode='both', burden=True,
                                         restore=True)
    if result['error'] is not None:
        print(f'FAIL: the trial reproduction raised {result["error"]}')
        return 1
    reproduced_values = gen.comparable_values(
        result['reproduced']['MPSPs'], result['reproduced']['metrics'])

    # 3. Compare (feeding first: the registry must carry the trial's).
    problems = []
    trial_feeding = result['feeding']
    for key, a, b in (('threshold_conc', feeding['threshold_conc'],
                       trial_feeding['threshold']),
                      ('target_conc', feeding['target_conc'],
                       trial_feeding['target']),
                      ('max_n_spikes', n_spikes,
                       trial_feeding['max_n_spikes'])):
        if abs(float(a) - float(b)) > 1e-9*max(1.0, abs(float(b))):
            problems.append(f'feeding {key}: scenario {a!r} vs trial {b!r}')
    try:
        deltas = gen.assert_reload_matches(reproduced_values,
                                           scenario_values, args.tol)
    except RuntimeError as e:
        problems.append(str(e))
        deltas = {}

    print(f'\n{args.scenario}: scenario baseline vs reproduction of trial '
          f'{trial_number} of {study_name}')
    print(f"  {'value':<24}{'scenario':>16}{'reproduction':>16}{'rel delta':>12}")
    for key, b in reproduced_values.items():
        a = scenario_values.get(key, float('nan'))
        d = deltas.get(key, float('nan'))
        print(f'  {key:<24}{a:>16.8g}{b:>16.8g}{d:>12.3g}')
    if problems:
        print('MISMATCH: ' + '; '.join(problems))
        return 1
    finite = [d for d in deltas.values() if d == d]
    print(f'MATCH: {len(deltas)} values within rel tol {args.tol:g} (max '
          f'{max(finite, default=0.0):.2e}); feeding strategy identical')
    return 0


if __name__ == '__main__':
    sys.exit(main())

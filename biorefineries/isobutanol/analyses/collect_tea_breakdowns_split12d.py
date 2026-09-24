#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Stage 1 of the TEA-breakdown figure of the metabolic_split_12d campaign
optima (spec docs/superpowers/specs/2026-09-24-tea-breakdowns-split12d-design.md;
stage 2, the sim-safe figure, is plots/plot_tea_breakdowns_split12d.py).

SIMULATES (ask-first; one fresh process, never next to another simulation on
a cold cache): one load() (default both-trains build, S201 split 1.0, as the
campaigns ran), the scenario-A baseline, then the objective-optimum trial of
each of the eight most recent ethanol_isobutanol x metabolic_split_12d
campaigns (CAMPAIGNS: the seven seed-350 `_rs350` campaigns and the flagship
PI (log-tail) relay `_rl15c111dc`), each re-simulated with
ko.reproduce_split12d_trial (anchor A, A-referenced burden ON, restore=True;
it leaves the flowsheet at the trial, so the unit groups read the trial).

After every simulation the five metrics of the 19 unit groups are read and
two closures pinned by analyses/test_unit_groups.py are asserted (the
installed-cost column sums to tea.installed_equipment_cost, the
Operating-cost column to tea.AOC / operating_hours); a failed simulation or
closure aborts. A reproduction whose tracked metrics differ from the
recorded row by more than 2 % is recorded and printed, not fatal.

Writes analyses/results/tea_breakdowns_split12d_<stamp>.json: meta (metric
names + units, group order, operating hours, hensmith path) and one record per
scenario (label, campaign, trial, summary metrics, reproduction check,
ABSOLUTE metric values per group; operating cost in USD/hr).

Run (the campaigns ran with hensmith master 2b5b27d pinned via PYTHONPATH;
use the same pin):
    python analyses/collect_tea_breakdowns_split12d.py
"""
import os
import csv
import sys
import json
import math
import argparse
import threading
from datetime import datetime

#%% Settings
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')
OUT_PREFIX = 'tea_breakdowns_split12d'
ANCHOR_SCENARIO = 'A'

_RS350 = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_{}_gp'
          '_rb0.001-4_ib0.75-1.5_aA_rs350_burden')
_RELAY = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_pi_log-tail_gp'
          '_rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden')
#: (key, figure label, study name, objective) in simulation order. A relay's
#: trajectory CSV holds only its SIMULATED trials; its preloaded donor rows
#: (best 0.504) sit below its simulated best, and a donor row could not be
#: reproduced under the relay's name anyway.
CAMPAIGNS = (
    ('profitability', 'Profitability', _RS350.format('pi_log-tail'),
     'PI (log-tail)'),
    ('flagship', 'Profitability (flagship)', _RELAY, 'PI (log-tail)'),
    ('ibo_yield', 'Isobutanol yield', _RS350.format('ibo_yield'), 'IBO yield'),
    ('ibo_titer', 'Isobutanol titer', _RS350.format('ibo_titer'), 'IBO titer'),
    ('ibo_productivity', 'Isobutanol productivity',
     _RS350.format('ibo_productivity'), 'IBO productivity'),
    ('etoh_yield', 'Ethanol yield', _RS350.format('etoh_yield'), 'EtOH yield'),
    ('etoh_titer', 'Ethanol titer', _RS350.format('etoh_titer'), 'EtOH titer'),
    ('etoh_productivity', 'Ethanol productivity',
     _RS350.format('etoh_productivity'), 'EtOH productivity'),
)
BASELINE_KEY, BASELINE_LABEL = 'baseline', 'Baseline'

#: the unit groups' metrics, in biosteam's autofill order (system.load renames
#: 'Material cost' -> 'Operating cost')
METRICS = ('Installed equipment cost', 'Cooling duty', 'Heating duty',
           'Electricity consumption', 'Operating cost')
INSTALLED, OPERATING = METRICS[0], METRICS[-1]
#: closure tolerances of analyses/test_unit_groups.py checks 1 and 5
INSTALLED_RTOL, OPERATING_RTOL = 1e-9, 1e-6
BASELINE_PIN_RTOL = 0.01
WATCHDOG_MIN = 45.0


class BreakdownClosureError(RuntimeError):
    """A unit-group column does not close on the TEA."""


#%% Pure helpers
def objective_optimum(csv_path):
    """(trial_number, objective) of the COMPLETE row with the largest finite
    'objective' (first occurrence on a tie, like DataFrame.idxmax)."""
    best = None
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            if row.get('state') != 'COMPLETE':
                continue
            try:
                value = float(row['objective'])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(value):
                continue
            if best is None or value > best[1]:
                best = (int(float(row['trial_number'])), value)
    if best is None:
        raise ValueError(f'{csv_path}: no COMPLETE row with a finite objective')
    return best


def _close(a, b, rtol):
    return abs(a - b) <= rtol*max(abs(a), abs(b), 1e-300)


def _json_float(x):
    """A metric value as a plain float (numpy scalars and None -> float/NaN)."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


#%% Model readers (need a loaded model)
def read_breakdown(system):
    """({group: {metric: value}}, {metric: units}, tea_totals) of the live
    flowsheet; raises BreakdownClosureError if a column does not close."""
    tea = system.corn_EtOH_IBO_sys_tea
    groups = system.unit_groups
    values, units = {}, {}
    for group in groups:
        names = tuple(m.name for m in group.metrics)
        if names != METRICS:
            raise KeyError(f'unit group {group.name!r} metrics {names} != '
                           f'{METRICS}')
        values[group.name] = {m.name: _json_float(m()) for m in group.metrics}
        units = {m.name: m.units for m in group.metrics}
    installed = sum(v[INSTALLED] for v in values.values())
    installed_ref = tea.installed_equipment_cost/1e6
    if not _close(installed, installed_ref, INSTALLED_RTOL):
        raise BreakdownClosureError(
            f'installed-cost column {installed!r} != '
            f'tea.installed_equipment_cost {installed_ref!r} MM$')
    operating = sum(v[OPERATING] for v in values.values())
    operating_ref = tea.AOC/tea.operating_hours
    if not _close(operating, operating_ref, OPERATING_RTOL):
        raise BreakdownClosureError(
            f'Operating-cost column {operating!r} != tea.AOC/operating_hours '
            f'{operating_ref!r} USD/hr')
    totals = dict(installed_equipment_cost=installed_ref,
                  AOC=tea.AOC/1e6, TCI=tea.TCI/1e6,
                  operating_hours=float(tea.operating_hours))
    return values, units, totals


def baseline_record(system, ko, scenarios):
    """Scenario-A record: load_scenario has just simulated the baseline."""
    handles = ko.get_handles()
    solution = handles['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'))
    handles['latest_TEA_solution'].update(solution)
    metrics = {name: _json_float(getter(handles))
               for name, getter in ko.TRACKED_METRICS.items()}
    mpsps = {k: _json_float(v) for k, v in solution['MPSPs'].items()}
    pin = scenarios.SCENARIOS[ANCHOR_SCENARIO].expected['ethanol']
    pin_ok = _close(mpsps['ethanol'], pin, BASELINE_PIN_RTOL)
    if not pin_ok:
        print(f'  WARNING baseline ethanol MPSP {mpsps["ethanol"]:.5f} is '
              f'not within 1 % of the smoke-test pin {pin}', flush=True)
    values, units, totals = read_breakdown(system)
    return dict(key=BASELINE_KEY, label=BASELINE_LABEL, campaign=None,
                trial_number=None, objective=None, objective_value=None,
                MPSPs=mpsps, IRR=_json_float(solution['IRR']),
                metrics=metrics, reproduction=dict(baseline_pin=pin,
                                                   baseline_pin_ok=pin_ok),
                totals=totals, breakdown=values), units


def campaign_record(system, ko, key, label, study, objective, trial,
                    recorded_objective):
    """Reproduce one trial and read its breakdown (raises on a failed
    simulation; a >2 % tracked-metric deviation is recorded, not raised)."""
    result = ko.reproduce_split12d_trial(ANCHOR_SCENARIO, study, trial,
                                         mode='both', burden=True,
                                         restore=True)
    if result['error'] is not None:
        raise RuntimeError(f'{key}: reproduction of trial {trial} of {study} '
                           f'failed: {result["error"]}')
    reproduced = result['reproduced']
    finite = [rel for name, _, _, rel in result['metric_check']
              if name not in ko.REPRODUCTION_DIAGNOSTIC_METRICS
              and math.isfinite(rel)]
    reproduction = dict(
        metric_warnings=list(result['metric_warnings']),
        max_rel_delta=max(finite) if finite else math.nan,
        cross_check_max_rel_delta=_json_float(
            result['cross_check'].get('max_rel_delta')),
        cross_check_mismatches=len(result['cross_check'].get('mismatches', ())),
        recorded={name: recorded for name, _, recorded, _
                  in result['metric_check']})
    if reproduction['metric_warnings']:
        print(f'  WARNING {key}: reproduced metrics differ from the recorded '
              f'row: {reproduction["metric_warnings"]}', flush=True)
    values, units, totals = read_breakdown(system)
    return dict(key=key, label=label, campaign=study, trial_number=trial,
                objective=objective, objective_value=recorded_objective,
                MPSPs={k: _json_float(v)
                       for k, v in reproduced['MPSPs'].items()},
                IRR=_json_float(reproduced['IRR']),
                metrics={k: _json_float(v)
                         for k, v in reproduced['metrics'].items()},
                reproduction=reproduction, totals=totals,
                breakdown=values), units


#%% Runner
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--results-dir', default=RESULTS_DIR)
    parser.add_argument('--out-dir', default=RESULTS_DIR)
    parser.add_argument('--watchdog-min', type=float, default=WATCHDOG_MIN,
                        help='exit 2 if the whole run exceeds this')
    args = parser.parse_args(argv)

    watchdog = threading.Timer(
        60.0*args.watchdog_min,
        lambda: (print(f'\nWATCHDOG: {args.watchdog_min:g} min exceeded -> '
                       'exit 2', flush=True), os._exit(2)))
    watchdog.daemon = True
    watchdog.start()

    # input errors raise before the ~20 s load
    selected = []
    for key, label, study, objective in CAMPAIGNS:
        path = os.path.join(args.results_dir, f'{study}_trajectory.csv')
        trial, value = objective_optimum(path)
        selected.append((key, label, study, objective, trial, value))
        print(f'{key:<18} optimum trial #{trial} ({objective} = {value:.6g})',
              flush=True)

    import hensmith
    from biorefineries import isobutanol
    from biorefineries.isobutanol import kinetic_optimization as ko
    for *_, objective, _, _ in selected:
        direction = ko.OBJECTIVE_REGISTRY[objective]['direction']
        if direction != 'maximize':
            raise ValueError(f'{objective!r} is a {direction} objective; '
                             'objective_optimum assumes maximize')
    isobutanol.load()
    from biorefineries.isobutanol import scenarios, system

    records = {}
    scenarios.load_scenario(ANCHOR_SCENARIO)
    print('\n=== baseline (scenario A) ===', flush=True)
    records[BASELINE_KEY], units = baseline_record(system, ko, scenarios)
    for key, label, study, objective, trial, value in selected:
        print(f'\n=== {key}: trial #{trial} ===', flush=True)
        records[key], units = campaign_record(system, ko, key, label, study,
                                              objective, trial, value)

    tea = system.corn_EtOH_IBO_sys_tea
    doc = dict(
        meta=dict(created=datetime.now().isoformat(timespec='seconds'),
                  spec='docs/superpowers/specs/'
                       '2026-09-24-tea-breakdowns-split12d-design.md',
                  anchor_scenario=ANCHOR_SCENARIO,
                  metrics=list(METRICS), metric_units=units,
                  groups=[g.name for g in system.unit_groups],
                  operating_hours=float(tea.operating_hours),
                  hensmith=os.path.dirname(hensmith.__file__)),
        order=[BASELINE_KEY] + [s[0] for s in selected],
        scenarios=records)
    stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
    path = os.path.join(args.out_dir, f'{OUT_PREFIX}_{stamp}.json')
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(doc, f, indent=1)
    os.replace(tmp, path)
    watchdog.cancel()

    print(f'\n{"scenario":<18}{"trial":>7}{"IRR":>9}{"PI":>9}{"TCI":>8}'
          f'{"AOC":>8}{"max rel":>10}  warnings')
    for key in doc['order']:
        r = records[key]
        rep = r['reproduction']
        print(f'{key:<18}{str(r["trial_number"] or "-"):>7}{r["IRR"]:>9.4f}'
              f'{r["metrics"]["PI"]:>9.4f}{r["totals"]["TCI"]:>8.1f}'
              f'{r["totals"]["AOC"]:>8.1f}'
              f'{rep.get("max_rel_delta", math.nan):>10.2e}  '
              f'{rep.get("metric_warnings", "")}')
    print(f'\nwrote {path}')
    return path


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
One-shot generator (ASK-FIRST to run; it simulates) for the five opt_*
kinetic-optimum scenario workbooks.

For each of the last five metabolic_minimal_subset studies (IRR, IBO
titer, IBO yield, EtOH titer, EtOH yield), it:
  1. Snapshots scenario A's full live kinetic state + isobutanol price.
  2. Selects the best-objective (argmax) COMPLETE & finite trial from the
     study's trajectory CSV and re-verifies its trial number against the
     spec table.
  3. Reproduces the trial's exact model state (A snapshot + 9 sampled
     rates + 16 applied_* inhibition coefficients + k_7<-k_7_eff,
     k_8<-k_8_eff), sets feeding (threshold, target=min(300,threshold+
     target_delta), max_n_spikes; spike 600 / stage_1_max_x 5.0 pinned),
     simulates (default both-trains build == IBO_EtOH-only at S201 split
     1.0), and cross-checks the reproduced objective against the CSV.
  4. Writes the workbook: copies scenario B, sets every kinetic row's
     Baseline to the reproduced value with a Triangular +/-20%
     distribution, sets the isobutanol-price row to the scenario-A value
     (B's relative spread preserved), keeps all other B rows.
  5. Reloads the written workbook and re-simulates to confirm it
     reproduces (defends against openpyxl save issues).

Finally prints a ScenarioSpec block to paste into scenarios.SCENARIOS.

Run (ask first):
  & "C:/Users/saran/anaconda3/envs/IBO_2026/python.exe" analyses/build_opt_scenario_workbooks.py
"""
import os
import csv
import math
import re

from biorefineries import isobutanol
isobutanol.load()   # default both-trains == IBO_EtOH-only at S201 split 1.0
from biorefineries.isobutanol import scenarios
from biorefineries.isobutanol import kinetic_optimization as ko

IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')
RESULTS_DIR = os.path.join(IBO_filepath, 'analyses', 'results')
WORKBOOK_DIR = os.path.join(IBO_filepath, 'analyses', 'full',
                            'parameter_distributions')
B_WORKBOOK = os.path.join(WORKBOOK_DIR,
                          'parameter-distributions_corn_IBO_EtOH_B.xlsx')

_STEM = 'kin_opt_ethanol_isobutanol_metabolic_minimal_subset_'
_TAIL = '_rb0.001-10_ib0.2-2_burden_trajectory.csv'

# (scenario name, objective registry key, study CSV, expected best trial #)
STUDIES = [
    ('opt_IRR',        'IRR',        _STEM + 'irr'        + _TAIL,  774),
    ('opt_IBO_titer',  'IBO titer',  _STEM + 'ibo_titer'  + _TAIL, 1732),
    ('opt_IBO_yield',  'IBO yield',  _STEM + 'ibo_yield'  + _TAIL, 1538),
    ('opt_EtOH_titer', 'EtOH titer', _STEM + 'etoh_titer' + _TAIL, 1820),
    ('opt_EtOH_yield', 'EtOH yield', _STEM + 'etoh_yield' + _TAIL, 1349),
]

# The 9 sampled rate columns (ko.METABOLIC_MINIMAL_SUBSET_RATES).
RATE_COLS = list(ko.METABOLIC_MINIMAL_SUBSET_RATES)

_LOAD_RE = re.compile(r'_te\.(\w+)\s*=')


def _finite(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def select_best_trial(csv_path, expected_trial):
    """Argmax-objective COMPLETE & finite row; verify trial number."""
    best = None
    with open(csv_path, newline='', encoding='utf-8') as fh:
        for row in csv.DictReader(fh):
            if row.get('state') != 'COMPLETE':
                continue
            obj = _finite(row.get('objective'))
            if obj is None:
                continue
            if best is None or obj > best[0]:
                best = (obj, row)
    if best is None:
        raise RuntimeError(f'No COMPLETE&finite trial in {csv_path}')
    obj, row = best
    trial = int(float(row['trial_number']))
    if trial != expected_trial:
        print(f'  WARNING: best trial #{trial} != spec-table #{expected_trial} '
              f'(objective {obj}); using the argmax #{trial}.')
    return row, obj


def apply_trial(r_te, a_snapshot, row):
    """Set the model's kinetic state to the trial's exact simulated state.
    Returns the applied {name: value} dict."""
    applied = dict(a_snapshot)                 # reset held/non-sampled params
    for name in RATE_COLS:
        applied[name] = float(row[name])
    for col, val in row.items():
        if col.startswith('applied_'):
            applied[col[len('applied_'):]] = float(val)
    applied['k_7'] = float(row['k_7_eff'])     # burden-derated, baked in
    applied['k_8'] = float(row['k_8_eff'])
    for name, value in applied.items():
        setattr(r_te, name, value)
    return applied


def simulate_and_measure(bundle, r_te, a_snapshot, row, objective_name):
    """Apply the trial, set feeding, simulate, return (applied, objective,
    MPSPs dict)."""
    applied = apply_trial(r_te, a_snapshot, row)
    bundle['fbs_spec'].max_n_spikes = int(float(row['max_n_spikes']))
    threshold = float(row['threshold_conc'])
    target = min(ko.TARGET_CONC_MAX, threshold + float(row['target_delta']))
    bundle['model_specification'](threshold_conc=threshold, target_conc=target)
    results = bundle['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'),
                                  IRR_for_MPSP=0.15)
    h = {'V406': bundle['V406'], 'tea': bundle['tea'],
         'latest_TEA_solution': results}
    obj = ko.OBJECTIVE_REGISTRY[objective_name]['getter'](h)
    return applied, obj, results['MPSPs'], threshold, target, \
        int(float(row['max_n_spikes']))


def write_workbook(out_path, applied, ibo_price_A):
    """Copy B; set every kinetic row Baseline to applied[name] with a
    Triangular +/-20% distribution; set the isobutanol-price row to the
    scenario-A value (B's relative spread preserved); keep other rows."""
    from openpyxl import load_workbook
    wb = load_workbook(B_WORKBOOK, data_only=True)
    ws = wb.active
    header = {str(c.value).strip(): i for i, c in enumerate(ws[1]) if c.value}
    col_load = header['Load statement']
    col_base = header['Baseline']
    col_shape = header['Shape']
    col_low = header['Lower']
    col_mid = header['Midpoint']
    col_up = header['Upper']
    n_kin = 0
    for r in ws.iter_rows(min_row=2):
        load_stmt = r[col_load].value
        if not load_stmt:
            continue
        load_stmt = str(load_stmt)
        m = _LOAD_RE.search(load_stmt)
        if m:                                   # kinetic row
            name = m.group(1)
            if name not in applied:
                print(f'  WARNING: kinetic row {name} not in applied dict; '
                      'left unchanged.')
                continue
            base = applied[name]
            r[col_base].value = base
            r[col_shape].value = 'Triangular'
            r[col_low].value = 0.8*base
            r[col_mid].value = 1.0*base
            r[col_up].value = 1.2*base
            n_kin += 1
        elif 'isobutanol_price' in load_stmt:   # price row -> scenario-A value
            b_base = r[col_base].value
            r[col_base].value = ibo_price_A
            if b_base:                          # preserve B's relative spread
                scale = ibo_price_A/float(b_base)
                for c in (col_low, col_mid, col_up):
                    if r[c].value is not None:
                        r[c].value = float(r[c].value)*scale
    wb.save(out_path)
    return n_kin


def reload_verify(name, out_path, objective_name, threshold, target,
                  max_n_spikes, expected_obj, expected_mpsps):
    """Load the just-written workbook fresh and re-simulate; confirm it
    reproduces the recorded objective + MPSPs."""
    model = isobutanol.models.models_EtOH_IBO_corn.model
    ns = isobutanol.models.namespace_dict
    fbs = isobutanol.models.fbs_spec
    ms = model.specification
    solve_TEA = isobutanol.system.solve_TEA
    V406 = model.system.flowsheet.V406
    tea = model.system.TEA
    model.parameters = ()
    model.load_parameter_distributions(out_path, ns)
    model.metrics_at_baseline()
    fbs.max_n_spikes = max_n_spikes
    ms(threshold_conc=threshold, target_conc=target)
    results = solve_TEA(stream_IDs=('ethanol', 'isobutanol'), IRR_for_MPSP=0.15)
    h = {'V406': V406, 'tea': tea, 'latest_TEA_solution': results}
    obj = ko.OBJECTIVE_REGISTRY[objective_name]['getter'](h)
    d_obj = abs(obj - expected_obj)/abs(expected_obj)
    print(f'  reload-verify {name}: objective {obj:.6g} '
          f'(recorded {expected_obj:.6g}, rel delta {d_obj:.2e}); '
          f"MPSPs ethanol {results['MPSPs']['ethanol']:.6g}, "
          f"isobutanol {results['MPSPs']['isobutanol']:.6g}")
    return obj, results['MPSPs']


def main():
    # 1. Snapshot scenario A's live kinetic state + isobutanol price.
    bundle = scenarios.load_scenario('A')
    r_te = ko.get_handles()['r_te']
    a_snapshot = ko.discover_kinetic_parameters(r_te)
    ibo_price_A = float(model_V514_price(bundle))
    print(f'Scenario-A snapshot: {len(a_snapshot)} kinetic params; '
          f'isobutanol price {ibo_price_A:.6g}')

    emitted = []
    for name, objective_name, csv_name, expected_trial in STUDIES:
        csv_path = os.path.join(RESULTS_DIR, csv_name)
        print(f'\n=== {name} ({objective_name}) : {csv_name} ===')
        row, recorded_obj = select_best_trial(csv_path, expected_trial)
        applied, repro_obj, mpsps, threshold, target, max_n_spikes = \
            simulate_and_measure(bundle, r_te, a_snapshot, row, objective_name)
        d = abs(repro_obj - recorded_obj)/abs(recorded_obj)
        print(f'  trial #{int(float(row["trial_number"]))}: reproduced '
              f'{objective_name} {repro_obj:.6g} vs recorded {recorded_obj:.6g} '
              f'(rel delta {d:.2e}); MPSPs ethanol {mpsps["ethanol"]:.6g}, '
              f'isobutanol {mpsps["isobutanol"]:.6g}')
        out_path = os.path.join(
            WORKBOOK_DIR,
            f'parameter-distributions_corn_IBO_EtOH_{name}.xlsx')
        n_kin = write_workbook(out_path, applied, ibo_price_A)
        print(f'  wrote {out_path} ({n_kin} kinetic rows)')
        reload_verify(name, out_path, objective_name, threshold, target,
                      max_n_spikes, recorded_obj, mpsps)
        emitted.append(dict(name=name, workbook=os.path.basename(out_path),
                            max_n_spikes=max_n_spikes, threshold=threshold,
                            target=target,
                            ethanol=mpsps['ethanol'],
                            isobutanol=mpsps['isobutanol'],
                            objective_name=objective_name,
                            objective_value=repro_obj))

    print('\n\n# ---- paste into scenarios.SCENARIOS ----')
    for e in emitted:
        print(
            f"    '{e['name']}': ScenarioSpec(\n"
            f"        name='{e['name']}',\n"
            f"        workbook='{e['workbook']}',\n"
            f"        max_n_spikes={e['max_n_spikes']}, "
            f"threshold_conc={e['threshold']!r}, target_conc={e['target']!r},\n"
            f"        expected={{'ethanol': {e['ethanol']!r}, "
            f"'isobutanol': {e['isobutanol']!r}}},\n"
            f"        objective_name={e['objective_name']!r}, "
            f"objective_value={e['objective_value']!r}),")


def model_V514_price(bundle):
    """The scenario-A isobutanol price (V514.isobutanol_price after
    load_scenario('A'); A's workbook has no price row, so this is the
    system.py default indexed to the price year)."""
    f = bundle['model'].system.flowsheet
    return f.V514.isobutanol_price


if __name__ == '__main__':
    main()

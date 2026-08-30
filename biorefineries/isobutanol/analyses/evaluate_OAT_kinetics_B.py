#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
One-at-a-time (OAT) sensitivity sweep of every non-discrete kinetic parameter
of the V406 fermentation kinetic model under scenario B.

Each global parameter of the nskinetics/tellurium model whose ID starts with
'k_' or 'K_' (rate constants, affinity constants, inhibition coefficients,
reversibility ratios) is swept over np.linspace(0, 10*baseline, N_STEPS) while
all other parameters stay at their scenario-B baselines; at each point the
full biorefinery is simulated and IRR, ethanol MPSP, and isobutanol MPSP are
recorded via the side-effect-free solve_TEA().

Results go to analyses/results/OAT_kinetics_B/, one CSV per parameter,
written row-by-row so the sweep is resumable: on relaunch, parameters (and
rows) already on disk are skipped. Failed simulation points are recorded as
NaN with the error message, matching the convention of the 2-D sweeps.

Plotting lives in plot_OAT_kinetics_B.py (reads the CSVs; no simulation).
"""

import os
import re
import csv
import numpy as np
from datetime import datetime

from warnings import filterwarnings
filterwarnings('ignore')

from biorefineries import isobutanol
isobutanol.load()

model = isobutanol.models.models_EtOH_IBO_corn.model
fbs_spec = isobutanol.models.models_EtOH_IBO_corn.fbs_spec
namespace_dict = isobutanol.models.namespace_dict
solve_TEA = isobutanol.system.solve_TEA
model_specification = model.specification
system = model.system

f = system.flowsheet
ferm_reactor = f.V406
r = ferm_reactor.nsk_kinetic_model._te

N_STEPS = 20
RANGE_MULT = 10.0  # sweep upper bound = RANGE_MULT * baseline

#%% Filepaths
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
results_dir = os.path.join(isobutanol_filepath, 'analyses', 'results',
                           'OAT_kinetics_B')
os.makedirs(results_dir, exist_ok=True)

#%% Load scenario B (parameter distributions -> kinetic baselines, then
#   the scenario-B feeding strategy), mirroring smoke_test_2.load_simulate_baseline
parameter_distributions_filename = isobutanol_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+\
    'parameter-distributions_corn_IBO_EtOH_B.xlsx'

model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
baseline_initial = model.metrics_at_baseline()

model_specification()
fbs_spec.max_n_spikes = 13
model_specification(threshold_conc=216.3, target_conc=226.3)

#%% Enumerate the non-discrete kinetic parameters (same filter as
#   utils.generate_save_kinetic_parameter_distributions) and their
#   scenario-B baselines from the live model
kinetic_params = [i for i in r.getGlobalParameterIds() if i[:2].lower() == 'k_']
baselines = {p: float(getattr(r, p)) for p in kinetic_params}

# Units, parsed from the antimony source ("<param> has <unit>;")
from nskinetics.models import s_cerevisiae_ferm_fb_inhib_mod_ibo as _nsk_mod
antimony_path = os.path.join(
    os.path.dirname(os.path.abspath(_nsk_mod.__file__)),
    's_cerevisiae_ferm_fb_inhib_mod_ibo_antimony.txt')
units = {}
with open(antimony_path) as _af:
    for m in re.finditer(r'^\s*(\w+) has (\w+);', _af.read(), re.M):
        units[m.group(1)] = m.group(2)

with open(os.path.join(results_dir, 'kinetic_params_baseline_B.csv'),
          'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['parameter', 'baseline_B', 'units',
                     'sweep_min', 'sweep_max'])
    for p in kinetic_params:
        writer.writerow([p, baselines[p], units.get(p, ''),
                         0.0, RANGE_MULT*baselines[p]])

print(f'\n{len(kinetic_params)} kinetic parameters: {kinetic_params}\n',
      flush=True)

#%% Sweep
CSV_HEADER = ['index', 'value', 'IRR', 'MPSP_ethanol', 'MPSP_isobutanol',
              'EtOH_titer', 'IBO_titer', 'cell_titer',
              'EtOH_yield', 'IBO_yield',
              'EtOH_productivity', 'IBO_productivity', 'tau', 'error']
N_METRICS = len(CSV_HEADER) - 3  # minus index, value, error

def sweep_csv_name(p):
    """Case-safe per-parameter CSV name: Windows filenames are
    case-insensitive, so OAT_B_k_1l.csv and OAT_B_K_1l.csv would be the
    SAME file -- uppercase (affinity K_*) parameters get a '_cap' suffix."""
    return f'OAT_B_{p}_cap.csv' if p[0] == 'K' else f'OAT_B_{p}.csv'

def n_completed_rows(csv_path):
    """Completed data rows in an existing sweep CSV; 0 if absent or if its
    header does not match CSV_HEADER (stale schema -> redo the sweep)."""
    if not os.path.exists(csv_path): return 0
    with open(csv_path, newline='') as csvfile:
        rows = list(csv.reader(csvfile))
    if not rows or rows[0] != CSV_HEADER: return 0
    return max(0, len(rows) - 1)  # minus header

def evaluate_point():
    """Metrics of the current (already simulated) flowsheet state."""
    TEA_solution = solve_TEA(stream_IDs=('ethanol', 'isobutanol'))
    nsk_res = ferm_reactor.nsk_results_specific_tau_dict
    tau = nsk_res.get('time', np.nan)  # the selected fermentation time, h
    IBO_titer = nsk_res.get('[s_IBO]', np.nan)
    return [TEA_solution['IRR'],
            TEA_solution['MPSPs']['ethanol'],
            TEA_solution['MPSPs']['isobutanol'],
            nsk_res.get('[s_EtOH]', np.nan),
            IBO_titer,
            nsk_res.get('[x]', np.nan),
            nsk_res.get('y_EtOH_glu_added', np.nan),
            nsk_res.get('y_IBO_glu_added', np.nan),
            nsk_res.get('prod_EtOH', np.nan),  # = [s_EtOH]/time
            IBO_titer/tau if tau else np.nan,
            tau]

#%% Baseline metrics at the current (already simulated) scenario-B state
#   (sanity check: ethanol MPSP ~ 0.396)
baseline_metrics = evaluate_point()
print('Scenario B baseline metrics:', flush=True)
print(dict(zip(CSV_HEADER[2:-1], baseline_metrics)), flush=True)
with open(os.path.join(results_dir, 'baseline_TEA_B.csv'),
          'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(CSV_HEADER[2:-1])
    writer.writerow(baseline_metrics)

curr_spec = {k: v for k, v in fbs_spec.current_specifications.items()}
total_no = len(kinetic_params) * N_STEPS
curr_no = 0
errors_dict = {}
t_start = datetime.now()

for i_p, p in enumerate(kinetic_params):
    baseline = baselines[p]
    values = np.linspace(0.0, RANGE_MULT*baseline, N_STEPS)
    csv_path = os.path.join(results_dir, sweep_csv_name(p))
    start_row = n_completed_rows(csv_path)
    curr_no += start_row
    if start_row >= N_STEPS:
        print(f'[{p}] already complete; skipping.', flush=True)
        continue
    if start_row == 0:
        with open(csv_path, 'w', newline='') as csvfile:
            csv.writer(csvfile).writerow(CSV_HEADER)
    else:
        print(f'[{p}] resuming at row {start_row}.', flush=True)

    for i_v in range(start_row, N_STEPS):
        v = values[i_v]
        curr_no += 1
        error_message = None
        row_metrics = [np.nan]*N_METRICS
        try:
            setattr(r, p, v)
            model_specification(**curr_spec)
            row_metrics = evaluate_point()
        except Exception as e:
            error_message = str(e).lower()
            print(f'Error in model spec: {error_message}', flush=True)
            if not 'specifications do not meet required' in error_message:
                errors_dict[(p, v)] = error_message
        with open(csv_path, 'a', newline='') as csvfile:
            csv.writer(csvfile).writerow(
                [i_v, v] + row_metrics + [error_message or ''])
        elapsed = (datetime.now() - t_start).total_seconds()
        print(f'{curr_no}/{total_no}  {p} = {v:.6g}  '
              f'(param {i_p+1}/{len(kinetic_params)}, point {i_v+1}/{N_STEPS})  '
              f'IRR={row_metrics[0]}, MPSP_EtOH={row_metrics[1]}, '
              f'MPSP_IBO={row_metrics[2]}  '
              f'[{elapsed:.0f} s elapsed]', flush=True)

    # Restore this parameter and re-converge the flowsheet at baseline so
    # the next parameter's sweep starts from a clean baseline state
    setattr(r, p, baseline)
    try:
        model_specification(**curr_spec)
        check = solve_TEA(stream_IDs=('ethanol',))
        print(f'[{p}] restored to {baseline:.6g}; baseline re-check '
              f'MPSP_EtOH = {check["MPSPs"]["ethanol"]}', flush=True)
    except Exception as e:
        print(f'[{p}] baseline restore simulation failed: {e}', flush=True)

print(f'\nDone. {len(errors_dict)} hard-errored points.', flush=True)
for k, v in errors_dict.items():
    print(f'  {k}: {v}', flush=True)
print(f'Total wall time: {(datetime.now()-t_start).total_seconds():.0f} s',
      flush=True)

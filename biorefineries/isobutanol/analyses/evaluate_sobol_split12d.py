#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Stage 1 of the Sobol' sensitivity analysis (spec
docs/superpowers/specs/2026-09-20-sobol-sensitivity-split12d-design.md):
simulate N_POINTS FEASIBLE scrambled-Sobol' points of the ethanol_isobutanol
x metabolic_split_12d campaign space -- same 12 variables, bands, groups,
references, scenario-A start, A-referenced enzyme burden ON, volume check ON,
default both-trains build -- through the campaigns' own evaluation site
(ko._prepare_optimization / ko.evaluate_decision_point), so every point is one
flushed trajectory row with all tracked metrics, the in-flight sidecar, LOST
recovery and reset-after-FAIL.

ASK-FIRST (4,096 simulations). Resumable: rerunning skips the rows already in
the CSV (every state; a point is never re-simulated); raising N_POINTS extends
the same nested design. Supervised launch:

    python supervise_sweep.py evaluate_sobol_split12d.py \
        --stem <STUDY_NAME> --checkpoint-suffix _trajectory.csv

Stage 2 (sim-safe): analyze_sobol_split12d.py."""
import json
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

#%% Settings
N_POINTS = int(os.environ.get('IBO_SOBOL_N_POINTS', 4096))
SEED = int(os.environ.get('IBO_SOBOL_SEED', 20260920))
ANCHOR = 'A'
STUDY_TARGET_PRODUCTS = 'ethanol_isobutanol'
STUDY_TYPE = 'metabolic_split_12d'

#%% Load (fresh kernel; one load per process)
from biorefineries import isobutanol
isobutanol.load()

from biorefineries.isobutanol import enzyme_burden as eb
from biorefineries.isobutanol import kinetic_optimization as ko
from biorefineries.isobutanol import scenarios
from biorefineries.isobutanol import sensitivity_analysis as sa
from biorefineries.isobutanol import system as ibo_system

#%% Campaign space, study name
preset, engine_kwargs = sa.campaign_engine_kwargs(ko, STUDY_TARGET_PRODUCTS, STUDY_TYPE)
assert preset['scenario'] == ANCHOR, preset['scenario']
# The campaign name of the same settings, re-badged: kin_opt_ -> kin_sobol_, no
# objective slug, + the design seed. Derived through ko.default_study_name (the
# driver's call) so a changed band shows up in the name.
campaign_name = ko.default_study_name(
    'PI', STUDY_TARGET_PRODUCTS, STUDY_TYPE,
    scenario=ANCHOR, kinetic_bounds_scenario=preset['kinetic_bounds_scenario'],
    burden=True, rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds'],
    inhibition_multiplier_bounds=engine_kwargs['group_multiplier_bounds'],
    exclude_params=engine_kwargs['exclude_params'],
    stage_1_max_x_bounds=engine_kwargs['stage_1_max_x_bounds'],
    n_seeds=0, method='tpe', ibo_pathway_anchoring='scenario_A')
prefix = f'kin_opt_{STUDY_TARGET_PRODUCTS}_{STUDY_TYPE}_pi_'
assert campaign_name.startswith(prefix), campaign_name
STUDY_NAME = (f'kin_sobol_{STUDY_TARGET_PRODUCTS}_{STUDY_TYPE}_'
              + campaign_name[len(prefix):] + f'_seed{SEED}')
print(f'STUDY_NAME={STUDY_NAME}', flush=True)

#%% Scenario + evaluation context (the engines' own set-up)
bundle = scenarios.load_scenario(ANCHOR, burden=True)
ctx = ko._prepare_optimization(
    'PI', direction=None, level=None, objective_units=None, objective_name=None,
    scenario_label=ANCHOR,
    threshold_conc_bounds=(0.0, 300.0), target_delta_bounds=(5.0, 500.0),
    max_n_spikes_bounds=(0, 50),          # run_kinetic_optimization's defaults
    target_conc_bounds=None, threshold_delta_bounds=None, spike_conc_bounds=None,
    study_name=STUDY_NAME, results_dir=None, handles=None,
    burden_model=bundle['burden_model'], volume_feasibility=True, volume_cap=None,
    seed_from=None, method_tag='', **engine_kwargs)
is_feasible = ko.feasibility_predicate(
    burden_on=True, volume_on=True, burden_model=ctx.burden_model,
    parameter_groups=ctx.parameter_groups, kinetic_baselines=ctx.kinetic_baselines,
    baseline_model_kwargs=ctx.baseline_model_kwargs,
    baseline_max_n_spikes=ctx.baseline_max_n_spikes,
    volume_cap=ctx.volume_cap, group_references=ctx.group_references)

#%% Design sidecar (everything stage 2 needs; identical on every resume)
def _git_commit():
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], text=True,
            cwd=os.path.dirname(os.path.abspath(__file__))).strip()
    except Exception:
        return 'unknown'

base = ctx.csv_path[:-len('_trajectory.csv')]
design = sa.design_record(
    search_space=ctx.search_space, parameter_groups=ctx.parameter_groups,
    group_references=ctx.group_references, kinetic_baselines=ctx.kinetic_baselines,
    burden_model=ctx.burden_model, eb=eb,
    baseline_model_kwargs=ctx.baseline_model_kwargs,
    baseline_max_n_spikes=ctx.baseline_max_n_spikes, volume_cap=ctx.volume_cap,
    target_conc_max=ko.TARGET_CONC_MAX,
    meta=dict(study_name=STUDY_NAME, seed=SEED, anchor=ANCHOR,
              study_target_products=STUDY_TARGET_PRODUCTS, study_type=STUDY_TYPE,
              git_commit=_git_commit(), python=sys.version.split()[0],
              numpy=np.__version__, created=datetime.now().isoformat(timespec='seconds')))
design_path = base + '_design.json'
if os.path.isfile(design_path):
    with open(design_path) as fh:
        old = json.load(fh)
    for key in ('search_space', 'parameter_groups', 'group_references', 'burden', 'volume'):
        if json.loads(json.dumps(design[key])) != old[key]:
            raise RuntimeError(f'{design_path}: {key!r} differs from the live model -- '
                               'the design changed since this study started; use a new SEED.')
else:
    with open(design_path, 'w') as fh:
        json.dump(design, fh, indent=1)

#%% Resume state: rows already logged (any state) are never re-simulated
n_done = 0
recorded = None
if os.path.isfile(ctx.csv_path) and os.path.getsize(ctx.csv_path) > 0:
    recorded = ko.load_trajectory(ctx.csv_path)
    n_done = len(recorded)
    assert list(recorded['trial_number']) == list(range(n_done)), \
        'trajectory trial numbers are not 0..n-1'
print(f'{n_done} of {N_POINTS} points already logged.', flush=True)

stream = sa.feasible_sobol_stream(ctx.search_space, is_feasible, SEED, ko.unit_to_external)
candidates = []
for n in range(n_done):                     # regenerate + verify the logged prefix
    candidate_index, values = next(stream)
    candidates.append(candidate_index)
    for name, value in values.items():
        if not np.isclose(float(recorded[name].iloc[n]), float(value), rtol=1e-9, atol=0.0):
            raise RuntimeError(f'row {n}: recorded {name} = {recorded[name].iloc[n]!r} but the '
                               f'regenerated stream gives {value!r} (seed / scipy changed?).')

def _write_candidates():
    with open(base + '_candidates.csv', 'w', newline='') as fh:
        fh.write('trial_number,candidate_index\n')
        fh.writelines(f'{n},{c}\n' for n, c in enumerate(candidates))

#%% Simulate
ibo_system.set_active_burden(ctx.burden_model)
try:
    for n in range(n_done, N_POINTS):
        candidate_index, values = next(stream)
        candidates.append(candidate_index)
        evaluation = ko.evaluate_decision_point(ctx, values, trial_number=n)
        if evaluation.state == 'INFEASIBLE':
            raise RuntimeError(f'point {n} passed the predicate but the evaluation site '
                               f'pruned it: {evaluation.record.get("error")}')
        if n % 50 == 0:
            _write_candidates()
            print(f'[{datetime.now():%H:%M:%S}] point {n}: {evaluation.state}; '
                  f'{candidates[-1] + 1} candidates drawn '
                  f'({(n + 1)/(candidates[-1] + 1):.1%} feasible).', flush=True)
finally:
    try:
        _write_candidates()
        ko.restore_baseline(ctx.handles, ctx.kinetic_baselines,
                            ctx.baseline_model_kwargs,
                            baseline_max_n_spikes=ctx.baseline_max_n_spikes,
                            baseline_stage_1_max_x=ctx.baseline_stage_1_max_x)
    finally:
        ibo_system.set_active_burden(None)

#%% Summary
df = ko.load_trajectory(ctx.csv_path)
print(f'Done: {len(df)} rows -- ' + ', '.join(
    f'{state} {count}' for state, count in df['state'].value_counts().items()))
print(f'Trajectory: {ctx.csv_path}\nDesign:     {design_path}')

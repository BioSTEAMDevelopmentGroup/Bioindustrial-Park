#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline verification of kinetic_optimization's pure logic -- search-space
construction, objective registry, trajectory-CSV round-trip, (from check 6)
plotting on synthetic data, and (checks 16-18) the lost-trial sidecar. No
isobutanol.load(), no simulation; optuna only in check 17, which drives the
real engine with fake handles (SKIP if optuna is absent).
Run in a fresh kernel; exit 0 + 'ALL ... CHECKS PASSED' = clean."""
import os
import tempfile
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

from biorefineries.isobutanol import kinetic_optimization as ko

n_pass = 0
def PASS(msg):
    global n_pass
    n_pass += 1
    print(f'PASS {n_pass}: {msg}')

#%% 1. build_search_space: bounds, zero-baseline exclusion, override, feeding vars
baselines = {'k_1e': 47.1, 'K_1e': 0.12, 'k_13': 0.0, 'k_16': 0.0}
space, excluded = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)})
assert excluded == ['k_16'] and 'k_16' not in space
assert space['k_13'] == dict(low=0.0, high=40.0, log=False)  # lo==0 -> linear
assert space['k_1e'] == dict(low=0.1*47.1, high=10.0*47.1, log=True)
assert space['K_1e'] == dict(low=0.1*0.12, high=10.0*0.12, log=True)
assert space['threshold_conc'] == dict(low=0.0, high=300.0, log=False)
assert space['target_delta'] == dict(low=5.0, high=500.0, log=False)
assert space['spike_delta'] == dict(low=0.5, high=595.0, log=False)
assert space['max_n_spikes'] == dict(low=0, high=50, log=False, int=True)
space_pinned, _ = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)},
    max_n_spikes_bounds=None)
assert 'max_n_spikes' not in space_pinned  # None -> pinned at baseline
# Any legacy kwarg -> legacy target-anchored scheme (unspecified legacy
# bounds fall back to the pre-2026-08-31 defaults).
space_legacy, _ = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)},
    target_conc_bounds=(180.0, 300.0))
assert space_legacy['target_conc'] == dict(low=180.0, high=300.0, log=False)
assert space_legacy['threshold_delta'] == dict(low=0.5, high=30.0, log=False)
assert space_legacy['spike_conc'] == dict(low=200.0, high=800.0, log=False)
assert 'threshold_conc' not in space_legacy
assert 'target_delta' not in space_legacy and 'spike_delta' not in space_legacy
PASS('build_search_space: multiplier band, zero-baseline exclusion, override, feeding vars')

#%% 2. explicit exclusion
space2, excluded2 = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)},
    exclude_params=('k_1e',))
assert 'k_1e' in excluded2 and 'k_1e' not in space2 and 'k_13' in space2
PASS('build_search_space: exclude_params honored')

#%% 3. objective registry completeness
required = {'IBO yield', 'IBO titer', 'IBO productivity',
            'EtOH yield', 'EtOH titer', 'EtOH productivity',
            'Combined yield', 'Cell density',
            'IRR', 'EtOH MPSP', 'IBO MPSP', 'TCI'}
assert required <= set(ko.OBJECTIVE_REGISTRY), required - set(ko.OBJECTIVE_REGISTRY)
for name, entry in ko.OBJECTIVE_REGISTRY.items():
    assert entry['direction'] in ('maximize', 'minimize'), name
    assert entry['level'] in ('kinetic', 'system'), name
    assert callable(entry['getter']) and isinstance(entry['units'], str), name
assert ko.OBJECTIVE_REGISTRY['IRR']['direction'] == 'maximize'
assert ko.OBJECTIVE_REGISTRY['EtOH MPSP']['direction'] == 'minimize'
assert ko.OBJECTIVE_REGISTRY['TCI']['direction'] == 'minimize'
assert ko.OBJECTIVE_REGISTRY['IBO yield']['level'] == 'kinetic'
assert ko.OBJECTIVE_REGISTRY['IRR']['level'] == 'system'
PASS('OBJECTIVE_REGISTRY: names, directions, levels, units')

#%% 4. getters against fake handles
nsk = {'y_IBO_glu_added': 0.1, '[s_IBO]': 20.0, 'time': 40.0,
       'y_EtOH_glu_added': 0.3, '[s_EtOH]': 90.0, 'prod_EtOH': 2.0,
       'y_EtOH_IBO_glu_added': 0.4, '[x]': 30.0, 'curr_n_glu_spikes': 7}
handles = {'V406': SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0),
           'tea': SimpleNamespace(TCI=350e6),
           'latest_TEA_solution': {'IRR': 0.21,
                                   'MPSPs': {'ethanol': 0.4, 'isobutanol': 0.9}}}
assert ko.OBJECTIVE_REGISTRY['IBO yield']['getter'](handles) == 0.1
assert ko.OBJECTIVE_REGISTRY['IBO productivity']['getter'](handles) == 0.5
assert ko.OBJECTIVE_REGISTRY['IBO yield x titer']['getter'](handles) == 0.1*20.0
assert ko.OBJECTIVE_REGISTRY['EtOH productivity']['getter'](handles) == 2.0
assert ko.OBJECTIVE_REGISTRY['Cell density']['getter'](handles) == 30.0
assert ko.OBJECTIVE_REGISTRY['IRR']['getter'](handles) == 0.21
assert ko.OBJECTIVE_REGISTRY['EtOH MPSP']['getter'](handles) == 0.4
assert ko.OBJECTIVE_REGISTRY['IBO MPSP']['getter'](handles) == 0.9
assert ko.OBJECTIVE_REGISTRY['TCI']['getter'](handles) == 350.0
assert ko.TRACKED_METRICS['tau'](handles) == 55.0
assert ko.TRACKED_METRICS['n_glu_spikes'](handles) == 7
assert set(ko.TRACKED_METRICS) == {'IBO yield', 'IBO titer', 'IBO productivity',
                                   'EtOH yield', 'EtOH titer', 'EtOH productivity',
                                   'Cell density', 'IRR', 'TCI',
                                   'tau', 'n_glu_spikes'}
PASS('getters read the handles contract correctly')

#%% 5. trajectory CSV: header, append, round-trip
outdir = tempfile.mkdtemp()
csv_path = os.path.join(outdir, 'traj.csv')
columns = ko.trajectory_columns(space)
assert columns[:2] == ['trial_number', 'state'] and columns[-1] == 'error'
assert 'objective' in columns and 'k_1e' in columns and 'IRR' in columns
rec = {'trial_number': 0, 'state': 'COMPLETE', 'objective': 0.2,
       'k_1e': 50.0, 'K_1e': 0.1, 'k_13': 5.0,
       'target_conc': 220.0, 'threshold_delta': 5.0, 'spike_conc': 600.0,
       'IRR': 0.2}
ko.append_trajectory_row(csv_path, columns, rec)
ko.append_trajectory_row(csv_path, columns,
                         dict(rec, trial_number=1, state='FAIL',
                              objective='', error='boom'))
df = ko.load_trajectory(csv_path)
assert list(df.columns) == columns and len(df) == 2
assert df['state'].tolist() == ['COMPLETE', 'FAIL']
assert df['objective'][0] == 0.2 and np.isnan(df['objective'][1])
PASS('trajectory CSV: columns, header-once, append, round-trip')

#%% 6. plot functions on a synthetic trajectory
rng = np.random.default_rng(0)
n = 60
synth = pd.DataFrame({'trial_number': np.arange(n)})
synth['state'] = ['FAIL' if i % 10 == 3 else 'COMPLETE' for i in range(n)]
synth['objective'] = rng.normal(0.15, 0.05, n)
for name, sp in space.items():
    synth[name] = rng.uniform(sp['low'], sp['high'], n)
for m in ko.TRACKED_METRICS:
    synth[m] = rng.normal(1.0, 0.1, n)
synth['error'] = ''
synth.loc[synth['state'] == 'FAIL', 'objective'] = np.nan

f1 = os.path.join(outdir, 'trajectories.png')
f2 = os.path.join(outdir, 'param_trajectory.png')
f3 = os.path.join(outdir, 'best_vs_baseline.png')
fig1, _ = ko.plot_optimization_trajectories(
    synth, objective_name='IRR', direction='maximize', filename=f1)
fig2, _ = ko.plot_parameter_trajectory(
    synth, baselines, direction='maximize', filename=f2)
fig3, _ = ko.plot_best_vs_baseline(
    synth, baselines, direction='maximize', filename=f3)
for p in (f1, f2, f3):
    assert os.path.isfile(p) and os.path.getsize(p) > 0, p
PASS('three plot functions save non-empty PNGs from synthetic data')

#%% 7. best-so-far
s = pd.Series([1.0, 0.5, 2.0, 1.5])
assert ko._best_so_far(s, 'maximize').tolist() == [1.0, 1.0, 2.0, 2.0]
assert ko._best_so_far(s, 'minimize').tolist() == [1.0, 0.5, 0.5, 0.5]
PASS('_best_so_far: cummax / cummin by direction')

#%% 8. changed search space rejected on append to an existing trajectory
try:
    ko.append_trajectory_row(csv_path, ['trial_number', 'state', 'bogus'],
                             rec)
    raise AssertionError('append with changed columns should have raised')
except ValueError as e:
    assert 'search space' in str(e) and 'study_name' in str(e)
ko.append_trajectory_row(csv_path, columns, rec)  # matching header still OK
assert len(ko.load_trajectory(csv_path)) == 3
PASS('append_trajectory_row: mismatched existing header raises ValueError')

#%% 9. baseline_decision_point: both feeding parameterizations
bmk = dict(target_conc=226.3, threshold_conc=216.3, spike_conc=632.0)
kb = {'k_1e': 47.1, 'K_1e': 0.12}
space9, _ = ko.build_search_space(kb)  # current threshold-anchored scheme
pt = ko.baseline_decision_point(space9, kb, bmk, baseline_max_n_spikes=13)
assert pt['k_1e'] == 47.1 and pt['K_1e'] == 0.12
assert pt['threshold_conc'] == 216.3
assert abs(pt['target_delta'] - 10.0) < 1e-12
assert abs(pt['spike_delta'] - (632.0 - 226.3)) < 1e-12
assert pt['max_n_spikes'] == 13
assert set(pt) == set(space9)  # exactly the decision variables, no extras
space9L, _ = ko.build_search_space(  # legacy target-anchored scheme
    kb, max_n_spikes_bounds=None, target_conc_bounds=(180.0, 300.0))
ptL = ko.baseline_decision_point(space9L, kb, bmk)
assert ptL['target_conc'] == 226.3
assert abs(ptL['threshold_delta'] - 10.0) < 1e-12
assert ptL['spike_conc'] == 632.0
assert 'max_n_spikes' not in ptL and set(ptL) == set(space9L)
PASS('baseline_decision_point: both feeding parameterizations')

#%% 10. baseline_decision_point: out-of-bounds baselines clipped into the space
# The scenario-A-with-B-bounds case: a zero-baseline kinetic parameter given
# absolute (log-scale) override bounds must enqueue at its low bound -- an
# enqueued 0.0 on a log distribution is a hard optuna ValueError at
# suggest_float, and an out-of-range linear value is silently replaced by a
# random draw. Same for a baseline feeding delta below its bound.
kb10 = {'k_13': 0.0, 'k_1e': 47.1}
space10, excl10 = ko.build_search_space(
    kb10, param_bounds_override={'k_13': (0.581, 58.1)})
assert not excl10 and space10['k_13'] == dict(low=0.581, high=58.1, log=True)
bmk10 = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
pt10 = ko.baseline_decision_point(space10, kb10, bmk10,
                                  baseline_max_n_spikes=55)
assert pt10['k_13'] == 0.581           # zero baseline -> clipped to low bound
assert pt10['k_1e'] == 47.1            # in-bounds baseline untouched
assert pt10['threshold_conc'] == 217.125
assert abs(pt10['target_delta'] - 5.0) < 1e-12   # 4.125 -> clipped to low 5.0
assert pt10['max_n_spikes'] == 50      # above high bound -> clipped, stays int
assert isinstance(pt10['max_n_spikes'], int)
PASS('baseline_decision_point: out-of-bounds baselines clipped into bounds')

#%% 11. pca_decision_matrix + plot_pca_projection on synthetic trajectories
# Reorder the synthetic trajectory into the real CSV column order --
# decision columns are inferred as everything between 'state' and
# 'objective' (trajectory_columns construction).
synth11 = synth[ko.trajectory_columns(space)].copy()
log_cols = {'k_1e', 'K_1e'}  # the log-sampled kinetic bands of `space`
coords, evr, load11, kept, valid = ko.pca_decision_matrix(synth11, log_cols)
assert valid.all() and valid.shape == (n,)
assert kept == list(space)                      # all columns vary
assert coords.shape == (n, len(space))
assert load11.shape == (len(space), len(space))
assert abs(evr.sum() - 1.0) < 1e-9
assert all(evr[i] >= evr[i+1] for i in range(len(evr) - 1))
# Zero-variance drop + missing-value mask
synth11b = synth11.copy()
synth11b['max_n_spikes'] = 7                    # pinned -> dropped
synth11b.loc[5, 'k_1e'] = np.nan                # incomplete row -> masked
c2, evr2, l2, kept2, valid2 = ko.pca_decision_matrix(synth11b, log_cols)
assert 'max_n_spikes' not in kept2 and set(kept2) < set(space)
assert not valid2[5] and valid2.sum() == n - 1 and c2.shape[0] == n - 1
# Known 1-D structure: two log-columns perfectly correlated after log10
# -> PC1 carries ~all variance; sign stabilized (max-|loading| positive);
# a nonpositive value in a log column invalidates its row (no crash).
n12 = 30
u = rng.uniform(-1.0, 1.0, n12)
df12 = pd.DataFrame({'trial_number': np.arange(n12),
                     'state': 'COMPLETE',
                     'k_1e': 10.0**u,
                     'K_1e': 10.0**(2.0*u),
                     'objective': u})
c3, evr3, l3, kept3, valid3 = ko.pca_decision_matrix(df12, log_cols)
assert evr3[0] > 0.999
assert l3[0][np.argmax(np.abs(l3[0]))] > 0.0
df12.loc[0, 'k_1e'] = 0.0                       # log10 -> -inf -> masked
_, _, _, _, valid3b = ko.pca_decision_matrix(df12, log_cols)
assert not valid3b[0] and valid3b.sum() == n12 - 1
# Too few rows -> ValueError
try:
    ko.pca_decision_matrix(df12.iloc[:2], log_cols)
    raise AssertionError('pca_decision_matrix on 2 rows should have raised')
except ValueError as e:
    assert 'Fewer than 3' in str(e)
# Plot smoke test (FAIL rows, baseline trial 0, best marker all exercised)
f4 = os.path.join(outdir, 'pca.png')
ko.plot_pca_projection(synth11, 'maximize', log_columns=log_cols,
                       objective_name='IRR', filename=f4)
assert os.path.isfile(f4) and os.path.getsize(f4) > 0
PASS('pca_decision_matrix + plot_pca_projection: transform, mask, EVR, plot')

#%% 12. StallGuard / attempt_outcome + supervised-runner helpers
g = ko.StallGuard(stall_timeout_s=100.0)
assert g.update(10, 0.0) is None                 # first poll initializes
assert g.update(10, 50.0) is None                # quiet but within timeout
assert g.update(11, 99.0) is None                # progress resets the clock
assert g.update(11, 150.0) is None
assert g.update(11, 199.0) == 'stalled'          # 100 s after last progress
g.reset()
assert g.update(11, 500.0) is None               # reset re-initializes
try:
    ko.StallGuard(stall_timeout_s=0.0)
    raise AssertionError('StallGuard(0.0) should have raised')
except ValueError:
    pass
assert ko.attempt_outcome(0, 5, 5) == 'complete'
assert ko.attempt_outcome(0, 5, 9, killed_for_stall=True) == 'resume'
assert ko.attempt_outcome(139, 5, 9) == 'resume'
assert ko.attempt_outcome(139, 5, 5) == 'abort'
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True) == 'abort'
# Supervised-runner helpers (run by file path -- stdlib-only module; its
# study naming MUST mirror the driver's or a resume forks the study).
import runpy as _runpy
sup = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert sup['default_study_name']('A', 'IRR', 'B') == 'kin_opt_A_kbB_irr'
assert sup['default_study_name']('B', 'IRR', None) == 'kin_opt_B_irr'
assert sup['default_study_name']('B', 'IBO titer', None) == 'kin_opt_B_ibo_titer'
assert sup['default_study_name'](None, 'IRR', None) == 'kin_opt_B_irr'   # legacy default scenario
assert sup['row_count'](os.path.join(outdir, 'nonexistent.csv')) == 0
assert sup['row_count'](csv_path) == 3           # check-5/8 trajectory
code12 = sup['child_code']('A', 'IRR', 2000, 'B', False, 'kin_opt_A_kbB_irr')
assert 'runpy' in code12 and "scenario='A'" in code12
assert "study_name='kin_opt_A_kbB_irr'" in code12
import inspect as _inspect12
assert _inspect12.signature(sup['child_code']).parameters['restrict_to_workbook'].default is True
assert _inspect12.signature(sup['supervise']).parameters['restrict_to_workbook'].default is True
assert 'restrict_to_workbook=False' in sup['child_code'](
    'B', 'IRR', 5, None, False, 'x', restrict_to_workbook=False)
assert 'restrict_to_workbook=True' in sup['child_code'](
    'B', 'IRR', 5, None, False, 'x')
# Lost-trial recovery wiring (2026-09-04): the supervisor logs the child's
# in-flight sidecar as a LOST row immediately after every attempt (and
# once at startup); the cause text is a pure helper.
assert sup['lost_cause'](True, 25.0, 1) == 'stall-killed after 25 min (no terminal row)'
assert sup['lost_cause'](True, 2.5, 1) == 'stall-killed after 2.5 min (no terminal row)'
assert sup['lost_cause'](False, 25.0, 3221225477) == \
    'child exited with code 3221225477 (no terminal row)'
assert callable(sup['ko'].recover_inflight) and callable(sup['ko'].inflight_path_for)
_src12 = _inspect12.getsource(sup['supervise'])
assert _src12.count('recover_inflight(') == 2               # startup + per attempt
assert _src12.index('attempt_outcome(') < _src12.rindex('recover_inflight(')  # decided BEFORE recovery
assert 'lost_cause(' in _src12 and 'inflight_path_for(' in _src12
# Abort-guard invariant on the building blocks, in the supervisor's order:
# a sidecar left by an attempt that added no terminal rows still yields
# 'abort' (decided on the pre-recovery count) and is THEN logged, so the
# trial is not silently lost.
outdir12 = tempfile.mkdtemp()
csv12 = os.path.join(outdir12, 's_trajectory.csv')
side12 = sup['ko'].inflight_path_for(outdir12, 's')
ko.append_trajectory_row(csv12, columns, rec)
rows_before12 = sup['row_count'](csv12)
ko.write_inflight(side12, columns, dict(rec, trial_number=1))
rows_after12 = sup['row_count'](csv12)
assert ko.attempt_outcome(1, rows_before12, rows_after12,
                          killed_for_stall=True) == 'abort'
assert ko.recover_inflight(csv12, side12, state='LOST',
                           error=sup['lost_cause'](True, 25.0, 1)) == 1
assert sup['row_count'](csv12) == rows_before12 + 1
assert ko.load_trajectory(csv12)['state'].tolist() == ['COMPLETE', 'LOST']
PASS('StallGuard / attempt_outcome / supervised-runner helpers')

#%% 13. include_params: whitelist restricts the kinetic set; feeding vars untouched
kb13 = {'k_1e': 47.1, 'K_1e': 0.12, 'k_6r': 3.0, 'K_2': 0.5, 'k_13': 0.0}
space13, excl13 = ko.build_search_space(
    kb13, include_params=['k_1e', 'K_1e', 'k_13'],
    param_bounds_override={'k_13': (0.581, 58.1)})
assert set(space13) == {'k_1e', 'K_1e', 'k_13', *ko.FEEDING_VARIABLES}
assert excl13 == ['k_6r', 'K_2']            # silent, model order preserved
assert space13['k_13'] == dict(low=0.581, high=58.1, log=True)  # override honored
assert space13['k_1e'] == dict(low=0.1*47.1, high=10.0*47.1, log=True)
# Composes with exclude_params (within the included set); a whitelisted
# name absent from the model is simply never placed; a non-whitelisted
# zero-baseline parameter is excluded by the whitelist (no warning path).
space13b, excl13b = ko.build_search_space(
    kb13, include_params=('k_1e', 'K_1e', 'k_6r', 'not_on_model'),
    exclude_params=('K_1e',))
assert set(space13b) == {'k_1e', 'k_6r', *ko.FEEDING_VARIABLES}
assert excl13b == ['K_1e', 'K_2', 'k_13']
assert 'not_on_model' not in space13b
# include_params=None -> no restriction (pre-2026-09-03 behaviour)
space13c, _ = ko.build_search_space(
    kb13, include_params=None,
    param_bounds_override={'k_13': (0.581, 58.1)})
assert set(space13c) == set(kb13) | set(ko.FEEDING_VARIABLES)
# baseline_decision_point on the restricted space: only in-space names,
# clipped into bounds (the start-at-A / set-from-B case: k_13 = 0 -> low)
pt13 = ko.baseline_decision_point(
    space13, kb13,
    dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0),
    baseline_max_n_spikes=16)
assert set(pt13) == set(space13)
assert pt13['k_13'] == 0.581 and pt13['k_1e'] == 47.1
assert pt13['max_n_spikes'] == 16 and 'k_6r' not in pt13
PASS('build_search_space: include_params whitelist composes with override/exclude; restricted baseline point')

#%% 14. workbook readers (plain file reads; skipped cleanly if the workbooks are absent)
wb_A = ko.parameter_distributions_workbook('A')
wb_B = ko.parameter_distributions_workbook('B')
assert wb_B.endswith('parameter-distributions_corn_IBO_EtOH_B.xlsx')
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    names_A = ko.kinetic_param_names_from_scenario('A')
    names_B = ko.kinetic_param_names_from_scenario('B')
    assert len(names_B) == 56 and len(names_A) == 40, (len(names_B), len(names_A))
    assert all(n[:2].lower() == 'k_' for n in names_A + names_B)
    assert len(set(names_B)) == len(names_B) and len(set(names_A)) == len(names_A)
    assert names_B[:5] == ['k_1l', 'K_1l', 'k_1h', 'K_1h', 'k_1e']  # workbook order
    for dropped in ('k_6r', 'k_16r', 'K_2', 'K_9'):   # commit 1e4efee1
        assert dropped not in names_B and dropped not in names_A, dropped
    assert 'k_13' in names_B and 'k_13' not in names_A  # IBO pathway: B only
    wb_baselines = ko.workbook_kinetic_baselines('B')
    assert list(wb_baselines) == names_B
    assert all(v > 0.0 for v in wb_baselines.values())
    PASS('workbook readers: 56 B / 40 A kinetic names in workbook order; constrained params absent')
else:
    print('SKIP 14: parameter-distribution workbooks not found')

#%% 15. run_kinetic_optimization: include_params kwarg plumbed (engine not run offline)
import inspect as _inspect
_sig = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'include_params' in _sig and _sig['include_params'].default is None
PASS('run_kinetic_optimization: include_params kwarg present, default None')

#%% 16. in-flight sidecar helpers: path, atomic write, LOST-row recovery, no-op, corrupt
# A hard-killed/segfaulted child cannot write its trial's terminal row; the
# sidecar it wrote beforehand is appended as ONE state='LOST' row carrying
# the trial's OWN trial_number (so the CSV has no gap), decision vector
# intact, objective/metrics blank (NaN on read-back), cause in 'error'.
outdir16 = tempfile.mkdtemp()
csv16 = os.path.join(outdir16, 'kin_opt_x_trajectory.csv')
side16 = ko.inflight_path_for(outdir16, 'kin_opt_x')
assert side16 == os.path.join(outdir16, 'kin_opt_x_inflight.json')
cols16 = ko.trajectory_columns(space)
rec16 = {'trial_number': 7, 'k_1e': 50.0, 'K_1e': 0.1, 'k_13': 5.0,
         'threshold_conc': 210.0, 'target_delta': 10.0, 'spike_delta': 300.0,
         'max_n_spikes': 3}
# Nothing pending: recovery and clearing are no-ops (no CSV created)
assert ko.recover_inflight(csv16, side16) is None
assert not os.path.isfile(csv16)
ko.clear_inflight(side16)
# Round trip
ko.write_inflight(side16, cols16, rec16)
assert os.path.isfile(side16) and not os.path.isfile(side16 + '.tmp')
assert ko.recover_inflight(
    csv16, side16, state='LOST',
    error='stall-killed after 25 min (no terminal row)') == 7
assert not os.path.isfile(side16)
df16 = ko.load_trajectory(csv16)
assert list(df16.columns) == cols16 and len(df16) == 1
row16 = df16.iloc[0]
assert row16['trial_number'] == 7 and row16['state'] == 'LOST'
for k16, v16 in rec16.items():
    assert row16[k16] == v16, (k16, row16[k16], v16)
assert np.isnan(row16['objective']) and np.isnan(row16['IRR'])
assert row16['error'] == 'stall-killed after 25 min (no terminal row)'
# Recovering again with nothing pending appends nothing
assert ko.recover_inflight(csv16, side16) is None
assert len(ko.load_trajectory(csv16)) == 1
# Corrupt / partial sidecars: warned, cleared, nothing appended, no raise
with open(side16, 'w') as f16:
    f16.write('{"columns": ["trial_number", "st')            # truncated JSON
assert ko.recover_inflight(csv16, side16) is None
assert not os.path.isfile(side16) and len(ko.load_trajectory(csv16)) == 1
with open(side16, 'w') as f16:
    f16.write('{"columns": ["trial_number"]}')               # no 'record'
assert ko.recover_inflight(csv16, side16) is None
assert not os.path.isfile(side16) and len(ko.load_trajectory(csv16)) == 1
# A later terminal row appends after the LOST row under the same header
ko.append_trajectory_row(csv16, cols16, dict(rec16, trial_number=8,
                                              state='COMPLETE', objective=0.2))
assert ko.load_trajectory(csv16)['state'].tolist() == ['LOST', 'COMPLETE']
assert set(ko.__all__) >= {'inflight_path_for', 'write_inflight',
                           'clear_inflight', 'recover_inflight'}
PASS('in-flight sidecar: path, atomic write, LOST-row recovery fills the gap, no-op, corrupt handled')

#%% 17. engine bracket: sidecar present DURING every trial, gone after; orphan recovered at start
# The real run_kinetic_optimization driven by FAKE handles (no biorefinery,
# no simulation): a scripted model_specification records whether the
# sidecar exists at call time and captures its content, raises on trial 1
# (FAIL path) and yields a NaN IRR on trial 2 (NAN path). optuna is used
# only here (skipped cleanly if absent).
try:
    import optuna as _optuna
except ImportError:
    _optuna = None
if _optuna is None:
    print('SKIP 17: optuna not installed')
else:
    import json as _json
    outdir17 = tempfile.mkdtemp()
    study17 = 'offline_bracket'
    csv17 = os.path.join(outdir17, study17 + '_trajectory.csv')
    side17 = ko.inflight_path_for(outdir17, study17)

    class _FakeTE:
        k_1e = 47.1
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'K_1e', 'not_kinetic']
    fbs17 = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    seen17 = []      # (sidecar present?, parsed sidecar or None) per call
    st17 = {'n': 0, 'irr': 0.2}
    def _model_specification(**kw):
        present = os.path.isfile(side17)
        with open(side17) if present else open(os.devnull) as f17:
            content = _json.load(f17) if present else None
        seen17.append((present, content))
        st17['n'] += 1
        if st17['n'] == 2:
            raise RuntimeError('boom')                          # trial 1 -> FAIL
        st17['irr'] = float('nan') if st17['n'] == 3 else 0.2  # trial 2 -> NAN
    def _solve_TEA(stream_IDs=None):
        return {'IRR': st17['irr'],
                'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles17 = {
        'r_te': _FakeTE(), 'fbs_spec': fbs17,
        'V406': SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0),
        'tea': SimpleNamespace(TCI=350e6), 'HXN': SimpleNamespace(),
        'model_specification': _model_specification,
        'solve_TEA': _solve_TEA,
        'latest_TEA_solution': {'IRR': np.nan,
                                'MPSPs': {'ethanol': np.nan,
                                          'isobutanol': np.nan}}}
    # An orphan left by a "previous run" that died mid-trial 99 (the
    # unsupervised driver's crash-resume case): must be logged, not
    # overwritten by trial 0's own sidecar.
    space17, _ = ko.build_search_space({'k_1e': 47.1, 'K_1e': 0.12})
    cols17 = ko.trajectory_columns(space17)
    ko.write_inflight(side17, cols17,
                      {'trial_number': 99, 'k_1e': 1.0, 'K_1e': 0.01,
                       'threshold_conc': 100.0, 'target_delta': 50.0,
                       'spike_delta': 100.0, 'max_n_spikes': 2})
    _optuna.logging.set_verbosity(_optuna.logging.WARNING)
    study17_obj, csv17_out, kb17 = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=3, seed=1,
        study_name=study17, results_dir=outdir17, handles=handles17,
        print_status_every=1, burden_model=None)
    assert csv17_out == csv17 and kb17 == {'k_1e': 47.1, 'K_1e': 0.12}
    df17 = ko.load_trajectory(csv17)
    assert df17['trial_number'].tolist() == [99, 0, 1, 2]
    assert df17['state'].tolist() == ['LOST', 'COMPLETE', 'FAIL', 'NAN']
    assert df17['error'][0] == 'recovered at engine startup (no terminal row)'
    assert df17['k_1e'][0] == 1.0 and np.isnan(df17['objective'][0])
    # 3 trials, then the restore_baseline call in the engine's finally
    assert [p for p, _ in seen17] == [True, True, True, False], seen17
    for i17 in range(3):
        side = seen17[i17][1]
        assert side['columns'] == cols17
        assert side['record']['trial_number'] == i17
        assert set(side['record']) == {'trial_number', 'k_1e', 'K_1e',
                                       *ko.FEEDING_VARIABLES}
        assert all(np.isfinite(v) for v in side['record'].values())
        # the sidecar held the SAME vector the CSV row later recorded
        assert np.isclose(side['record']['k_1e'], df17['k_1e'][i17 + 1],
                          rtol=1e-12, atol=0.0)
    assert not os.path.isfile(side17) and not os.path.isfile(side17 + '.tmp')

    # An INTERRUPTED trial (KeyboardInterrupt during the simulation -- a
    # manual abort) writes no terminal row: the sidecar must survive so the
    # next engine start recovers it as a LOST row (otherwise the trial
    # number is consumed and silently missing).
    outdir17b = tempfile.mkdtemp()
    study17b = 'offline_bracket_interrupt'
    csv17b = os.path.join(outdir17b, study17b + '_trajectory.csv')
    side17b = ko.inflight_path_for(outdir17b, study17b)
    def _model_specification_interrupt(**kw):
        raise KeyboardInterrupt
    handles17b = dict(handles17, model_specification=_model_specification_interrupt,
                      latest_TEA_solution={'IRR': np.nan,
                                           'MPSPs': {'ethanol': np.nan,
                                                     'isobutanol': np.nan}})
    try:
        ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=2, seed=1,
            study_name=study17b, results_dir=outdir17b, handles=handles17b,
            print_status_every=1, burden_model=None)
    except KeyboardInterrupt:
        pass
    else:
        raise AssertionError('KeyboardInterrupt did not propagate out of the engine')
    assert os.path.isfile(side17b), 'sidecar was cleared although no terminal row was written'
    assert not os.path.isfile(csv17b) or len(ko.load_trajectory(csv17b)) == 0
    # The next start of the same study recovers it as a LOST row for trial 0
    assert ko.recover_inflight(csv17b, side17b, state='LOST',
                               error='recovered at engine startup (no terminal row)') == 0
    df17b = ko.load_trajectory(csv17b)
    assert df17b['trial_number'].tolist() == [0] and df17b['state'].tolist() == ['LOST']
    assert not os.path.isfile(side17b)

    PASS('engine bracket: sidecar holds the full vector during each trial, cleared on COMPLETE/FAIL/NAN, orphan recovered at start, survives an interrupted trial')

#%% 18. plot safety with LOST rows: excluded from completed-only plots, drawn in the PCA landscape
synth18 = synth11.copy()
lost18 = synth18['trial_number'].isin([7, 22])
synth18.loc[lost18, 'state'] = 'LOST'
synth18.loc[lost18, ['objective', *ko.TRACKED_METRICS]] = np.nan  # blank cells on read-back
synth18.loc[lost18, 'error'] = 'stall-killed after 25 min (no terminal row)'
ok18 = ko._completed(synth18)
assert (ok18['state'] == 'COMPLETE').all()
assert len(ok18) == int((synth18['state'] == 'COMPLETE').sum())
assert not ok18['trial_number'].isin([7, 22]).any()
f18 = [os.path.join(outdir, f'lost_{i}.png') for i in range(4)]
ko.plot_optimization_trajectories(synth18, objective_name='IRR',
                                  direction='maximize', filename=f18[0])
ko.plot_parameter_trajectory(synth18, baselines, direction='maximize',
                             filename=f18[1])
ko.plot_best_vs_baseline(synth18, baselines, direction='maximize',
                         filename=f18[2])
# LOST rows keep a finite decision vector -> they join the PCA fit and are
# drawn as their own crosses (legend entry present).
_, _, _, _, valid18 = ko.pca_decision_matrix(synth18, log_cols)
assert valid18.all()
fig18, axes18 = ko.plot_pca_projection(synth18, 'maximize',
                                       log_columns=log_cols,
                                       objective_name='IRR', filename=f18[3])
labels18 = axes18[0].get_legend_handles_labels()[1]
assert 'lost (stalled/crashed)' in labels18, labels18
assert 'failed (pruned)' in labels18 and 'completed' in labels18
for p in f18:
    assert os.path.isfile(p) and os.path.getsize(p) > 0, p
# A LOST-free trajectory draws no LOST legend entry
_, axes18b = ko.plot_pca_projection(synth11, 'maximize', log_columns=log_cols)
assert 'lost (stalled/crashed)' not in axes18b[0].get_legend_handles_labels()[1]
PASS('LOST rows: excluded by _completed, drawn as crosses in the PCA landscape, other plots unaffected')

#%% 19. rate_multiplier_bounds: k_* band separate from the K_* band; None = legacy single band
kb19 = {'k_1e': 47.1, 'K_1e': 0.12, 'k_7': 0.5, 'K_1i': 2.0, 'k_13': 0.0}
assert ko.DEFAULT_RATE_MULTIPLIER_BOUNDS == (1e-3, 10.0)   # 1e-5 before 2026-09-06 (pm)
assert ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS == (0.1, 10.0)
space19, excl19 = ko.build_search_space(
    kb19, multiplier_bounds=ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS,
    rate_multiplier_bounds=ko.DEFAULT_RATE_MULTIPLIER_BOUNDS)
assert space19['k_1e'] == dict(low=1e-3*47.1, high=10.0*47.1, log=True)
assert space19['k_7'] == dict(low=1e-3*0.5, high=10.0*0.5, log=True)
assert space19['K_1e'] == dict(low=0.1*0.12, high=10.0*0.12, log=True)   # uppercase: saturation band
assert space19['K_1i'] == dict(low=0.1*2.0, high=10.0*2.0, log=True)
assert excl19 == ['k_13']                       # zero baseline still excluded
# Precedence unchanged: override beats the band, whitelist beats everything.
space19b, _ = ko.build_search_space(
    kb19, rate_multiplier_bounds=(1e-5, 10.0),
    param_bounds_override={'k_1e': (1.0, 2.0)}, include_params=['k_1e', 'K_1e'])
assert space19b['k_1e'] == dict(low=1.0, high=2.0, log=True)
assert set(space19b) == {'k_1e', 'K_1e', *ko.FEEDING_VARIABLES}
# rate_multiplier_bounds=None reproduces the previous single-band space EXACTLY.
space19c, excl19c = ko.build_search_space(kb19)
expected19c = {
    'k_1e': dict(low=0.1*47.1, high=10.0*47.1, log=True),
    'K_1e': dict(low=0.1*0.12, high=10.0*0.12, log=True),
    'k_7': dict(low=0.1*0.5, high=10.0*0.5, log=True),
    'K_1i': dict(low=0.1*2.0, high=10.0*2.0, log=True),
    'threshold_conc': dict(low=0.0, high=300.0, log=False),
    'target_delta': dict(low=5.0, high=500.0, log=False),
    'spike_delta': dict(low=0.5, high=595.0, log=False),
    'max_n_spikes': dict(low=0, high=50, log=False, int=True)}
assert space19c == expected19c and excl19c == ['k_13']
assert ko.build_search_space(kb19, rate_multiplier_bounds=None)[0] == expected19c
# Engine kwarg plumbed (default None).
_sig19 = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'rate_multiplier_bounds' in _sig19 and _sig19['rate_multiplier_bounds'].default is None
# workbook_kinetic_bounds: absolute per-prefix bands around the workbook baselines.
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    wbb19 = ko.workbook_kinetic_bounds('B', multiplier_bounds=(0.1, 10.0),
                                       rate_multiplier_bounds=(1e-5, 10.0))
    base19 = ko.workbook_kinetic_baselines('B')
    assert list(wbb19) == list(base19)                    # workbook order, every row (all > 0)
    for n19, b19 in base19.items():
        lo19, hi19 = wbb19[n19]
        assert hi19 == 10.0*b19
        assert lo19 == (1e-5*b19 if n19.startswith('k_') else 0.1*b19), n19
    wbb19_legacy = ko.workbook_kinetic_bounds('B')
    assert all(wbb19_legacy[n] == (0.1*b, 10.0*b) for n, b in base19.items())
else:
    print('SKIP 19b: parameter-distribution workbooks not found')
PASS('rate_multiplier_bounds: per-prefix bands, precedence, None = legacy space; workbook_kinetic_bounds')

#%% 20. kinetic_parameter_roles: nskinetics role table read by FILE PATH (no package import)
import subprocess as _subprocess
import sys as _sys
roles_path20 = ko.kinetic_parameter_roles_path()
assert roles_path20.endswith(os.path.join(
    'models', 's_cerevisiae_ferm_fb_inhib_mod_ibo', 'parameter_categories.py'))
assert os.path.isfile(roles_path20), roles_path20
# The no-heavy-import guarantee is probed in a FRESH interpreter that loads
# the engine by file path (as the stdlib-only supervisor does): this script
# imports ko through the biorefineries.isobutanol package, whose system.py
# imports nskinetics at module top, so sys.modules here proves nothing.
_probe20 = (
    "import importlib.util, sys\n"
    f"spec = importlib.util.spec_from_file_location('ko_probe', {ko.__file__!r})\n"
    "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
    "roles = m.kinetic_parameter_roles()\n"
    "heavy = sorted(k for k in sys.modules if k.split('.')[0] in\n"
    "               ('nskinetics', 'tellurium', 'roadrunner', 'biosteam', 'thermosteam'))\n"
    "print(len(roles), heavy)\n")
_out20 = _subprocess.run([_sys.executable, '-c', _probe20],
                         capture_output=True, text=True)
assert _out20.returncode == 0, _out20.stderr
assert _out20.stdout.strip() == '65 []', _out20.stdout + _out20.stderr
roles20 = ko.kinetic_parameter_roles()
from collections import Counter as _Counter
assert _Counter(roles20.values()) == {
    'capacity': 20, 'affinity': 16, 'product_inhibition': 13,
    'substrate_regulation': 4, 'product_self_inhibition': 4,
    'lethality': 3, 'lethality_threshold': 3, 'initial_state': 2}
assert roles20['k_7'] == 'capacity' and roles20['K_1i'] == 'substrate_regulation'
assert roles20['K_6e'] == 'product_self_inhibition' and roles20['k_10ii'] == 'lethality'
assert roles20['K_13'] == 'affinity' and roles20['k_16ie'] == 'product_inhibition'
assert list(roles20)[:3] == ['k_1h', 'k_1l', 'k_1e']          # table order kept
assert ko.kinetic_parameter_roles() is roles20               # cached
assert ko.kinetic_parameter_roles(path=roles_path20) == roles20   # explicit path: fresh, equal
assert ko.kinetic_parameter_roles(path=roles_path20) is not roles20
PASS('kinetic_parameter_roles: role table loaded by file path, no heavy import in a fresh interpreter, counts pinned, cached')

#%% 21. resolve_study_preset: set sizes by (target products x study type), roles, errors
assert ko.DEFAULT_STUDY_TARGET_PRODUCTS == 'ethanol_isobutanol'
assert ko.DEFAULT_STUDY_TYPE == 'metabolic_protein'
assert set(ko.STUDY_TARGET_PRODUCTS) == {'ethanol_only', 'ethanol_isobutanol'}
assert set(ko.STUDY_TYPE_ROLES) == {'metabolic', 'metabolic_protein', 'metabolic_minimal',
                                    'metabolic_minimal_subset'}
assert set(ko.STUDY_TYPE_ROLES['metabolic_minimal']) == {
    'capacity', 'product_inhibition', 'lethality'}
assert set(ko.STUDY_TYPE_ROLES['metabolic']) == {
    'capacity', 'product_inhibition', 'lethality', 'substrate_regulation'}
assert set(ko.STUDY_TYPE_ROLES['metabolic_protein']) == (
    set(ko.STUDY_TYPE_ROLES['metabolic']) | {'affinity', 'product_self_inhibition'})
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'
assert ko.default_study_name('IBO titer', 'ethanol_only', 'metabolic') \
    == 'kin_opt_ethanol_only_metabolic_ibo_titer'
# An explicit scenario / kinetic_bounds_scenario tags the name ONLY when
# it differs from the preset's own values (Finding 1: an override that
# matches silently used to look identical to a genuine override too).
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             scenario='A', kinetic_bounds_scenario='B') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'   # preset's own values: no tag
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             scenario='B') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_scB'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             kinetic_bounds_scenario='A') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_kbA'
assert ko.default_study_name('IRR', 'ethanol_only', 'metabolic',
                             scenario='B', kinetic_bounds_scenario='B') \
    == 'kin_opt_ethanol_only_metabolic_irr_scB_kbB'
for bad21 in (('ethanol', 'metabolic'), ('ethanol_only', 'protein')):
    try:
        ko.resolve_study_preset(*bad21)
    except ValueError as e21:
        assert 'ethanol_only' in str(e21) or 'metabolic_protein' in str(e21)
    else:
        raise AssertionError(f'unknown preset {bad21} did not raise ValueError')
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    expected21 = {('ethanol_only', 'metabolic'): 29,
                  ('ethanol_only', 'metabolic_protein'): 40,
                  ('ethanol_isobutanol', 'metabolic'): 40,
                  ('ethanol_isobutanol', 'metabolic_protein'): 56}
    roles21 = ko.kinetic_parameter_roles()
    for (stp21, st21), n21 in expected21.items():
        p21 = ko.resolve_study_preset(stp21, st21)
        assert set(p21) == {'scenario', 'kinetic_bounds_scenario', 'include_params',
                            'multiplier_bounds', 'rate_multiplier_bounds',
                            'rate_params', 'parameter_multiplier_bounds',
                            'exclude_params', 'stage_1_max_x_bounds',
                            'parameter_groups', 'group_multiplier_bounds',
                            'spike_delta_bounds'}
        # The three 2026-09-07 keys are inert on every pre-existing preset.
        assert p21['parameter_groups'] is None
        assert p21['group_multiplier_bounds'] == ko.DEFAULT_GROUP_MULTIPLIER_BOUNDS
        assert p21['spike_delta_bounds'] == ko.DEFAULT_SPIKE_DELTA_BOUNDS
        assert p21['scenario'] == 'A'                       # both start at the A baseline
        # k_10 is excluded from every preset by default (2026-09-06 pm): a
        # lower decay rate is a free lunch; it stays in include_params (the
        # workbook set) and is removed by build_search_space.
        assert p21['exclude_params'] == ('k_10',)
        assert p21['exclude_params'] == ko.DEFAULT_EXCLUDED_PARAMETERS
        assert p21['kinetic_bounds_scenario'] == ('A' if stp21 == 'ethanol_only' else 'B')
        assert p21['multiplier_bounds'] == (0.1, 10.0)
        assert p21['rate_multiplier_bounds'] == (1e-3, 10.0)
        # k_10 (active-biomass decay capacity, r10) keeps a 0.1x floor:
        # a near-zero decay rate is not an engineering target.
        assert p21['parameter_multiplier_bounds'] == {'k_10': (0.1, 10.0)}
        assert p21['parameter_multiplier_bounds'] is not ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS
        assert 'k_10' in p21['rate_params'] and 'k_10' in p21['include_params']
        inc21 = p21['include_params']
        assert len(inc21) == n21, (stp21, st21, len(inc21))
        wb21 = ko.kinetic_param_names_from_scenario(p21['kinetic_bounds_scenario'])
        assert inc21 == [n for n in wb21 if n in inc21]     # workbook order, no extras
        assert all(roles21[n] in ko.STUDY_TYPE_ROLES[st21] for n in inc21)
        # rate_params: the RATE CONSTANTS (role capacity) of the set
        # scenario's workbook, workbook order -- the only names the k_*
        # band applies to (inhibition coefficients k_*i* share the K_*
        # band since 2026-09-06). 16 in A's workbook, 20 in B's (+ the
        # four Ehrlich capacities); every one is in every study type.
        rp21 = p21['rate_params']
        assert rp21 == [n for n in wb21 if roles21[n] == 'capacity']
        assert len(rp21) == (16 if stp21 == 'ethanol_only' else 20)
        assert all(n in inc21 for n in rp21)
        assert rp21 == ko.rate_constant_names(wb21, roles=roles21)
        assert not any(n in rp21 for n in ('k_1ie', 'k_1ia', 'k_7ie', 'k_10ie'))
        if st21 == 'metabolic':
            assert [n for n in inc21 if n.startswith('K_')] == ['K_1i', 'K_2i', 'K_5i', 'K_9i']
            assert all(n.startswith('k_') or n.startswith('K_') for n in inc21)
        else:
            assert inc21 == wb21                            # every workbook row
    p21_eo = ko.resolve_study_preset('ethanol_only', 'metabolic_protein')['include_params']
    for ibo21 in ('k_13', 'K_16i', 'k_1ii', 'k_7ii', 'k_10ii', 'k_16ie'):
        assert ibo21 not in p21_eo, ibo21
    # A workbook row missing from the role table must raise, not leak.
    roles_missing21 = dict(roles21); del roles_missing21['k_7']
    try:
        ko.resolve_study_preset('ethanol_only', 'metabolic', roles=roles_missing21)
    except KeyError as e21:
        assert 'k_7' in str(e21)
    else:
        raise AssertionError('missing role-table entry did not raise KeyError')
    PASS('resolve_study_preset: 29/40/40/56 sets, role filter, A start, workbook order, errors')
else:
    print('SKIP 21: parameter-distribution workbooks not found')

#%% 22. trial 0 of (ethanol_isobutanol, metabolic_protein): A baseline clipped into the preset bands
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    p22 = ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_protein')
    base22_A = ko.workbook_kinetic_baselines('A')
    base22_B = ko.workbook_kinetic_baselines('B')
    # What discover_kinetic_parameters returns on the model under A's
    # distributions: A's baselines, and 0 for the Ehrlich capacities that
    # A leaves switched off (the other B-only rows keep the model's value,
    # which equals B's baseline).
    ehrlich22 = ('k_13', 'k_14', 'k_15', 'k_16')
    model22 = {n: (0.0 if n in ehrlich22 else base22_A.get(n, base22_B[n]))
               for n in base22_B}
    override22 = ko.workbook_kinetic_bounds(
        p22['kinetic_bounds_scenario'],
        multiplier_bounds=p22['multiplier_bounds'],
        rate_multiplier_bounds=p22['rate_multiplier_bounds'],
        rate_params=p22['rate_params'],
        parameter_multiplier_bounds=p22['parameter_multiplier_bounds'])
    space22, excl22 = ko.build_search_space(
        model22, multiplier_bounds=p22['multiplier_bounds'],
        rate_multiplier_bounds=p22['rate_multiplier_bounds'],
        rate_params=p22['rate_params'],
        parameter_multiplier_bounds=p22['parameter_multiplier_bounds'],
        param_bounds_override=override22, include_params=p22['include_params'],
        exclude_params=p22['exclude_params'])
    # k_10 (active-biomass decay) is excluded by default (2026-09-06 pm):
    # 55 sampled kinetic parameters, no k_10 column, no k_10 probe.
    assert excl22 == ['k_10'] and len(space22) == 55 + 4
    assert 'k_10' not in space22
    assert all(space22[n]['log'] for n in base22_B if n != 'k_10')
    # Bands by ROLE (2026-09-06): rate constants (capacity) 1e-3x-10x
    # (1e-5x until later that day); inhibition coefficients (k_*i*:
    # product_inhibition, lethality) and the K_* terms (regulation,
    # affinity, self-inhibition) 0.1x-10x.
    for n22, b22 in base22_B.items():
        if n22 == 'k_10':
            continue
        assert space22[n22]['high'] == 10.0*b22
        assert space22[n22]['low'] == (1e-3*b22 if n22 in p22['rate_params']
                                       else 0.1*b22), n22
    for inh22 in ('k_1ie', 'k_1ii', 'k_7ii', 'k_10ie', 'k_10ii', 'k_16ie'):
        assert space22[inh22]['low'] == 0.1*base22_B[inh22], inh22
    for rate22 in ('k_1h', 'k_2', 'k_7', 'k_13'):
        assert space22[rate22]['low'] == 1e-3*base22_B[rate22], rate22
    # The workbook bounds still carry k_10's per-parameter band (in force
    # only when a caller re-includes it with exclude_params=()).
    assert override22['k_10'] == (0.1*0.06, 10.0*0.06)
    space22_k10, excl22_k10 = ko.build_search_space(
        model22, multiplier_bounds=p22['multiplier_bounds'],
        rate_multiplier_bounds=p22['rate_multiplier_bounds'],
        rate_params=p22['rate_params'],
        parameter_multiplier_bounds=p22['parameter_multiplier_bounds'],
        param_bounds_override=override22, include_params=p22['include_params'],
        exclude_params=())
    assert excl22_k10 == [] and len(space22_k10) == 56 + 4
    assert space22_k10['k_10'] == dict(low=0.1*0.06, high=10.0*0.06, log=True)
    assert {n: v for n, v in space22_k10.items() if n != 'k_10'} == space22
    pt22 = ko.baseline_decision_point(
        space22, model22,
        dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0),
        baseline_max_n_spikes=16)
    assert 'k_10' not in pt22
    for n22 in ehrlich22:
        assert pt22[n22] == 1e-3*base22_B[n22], n22       # clipped to the floor, exactly
    for n22, v22 in model22.items():
        if n22 not in ehrlich22 and n22 != 'k_10':
            assert pt22[n22] == v22, n22                  # nonzero baselines untouched
    assert pt22['threshold_conc'] == 217.125 and pt22['max_n_spikes'] == 16
    probes22, at_floor22 = ko.knockout_probe_points(space22, pt22,
                                                    rate_params=p22['rate_params'])
    assert 'k_10' not in probes22 and 'k_10' not in at_floor22
    assert set(at_floor22) == set(ehrlich22)
    assert len(probes22) == 20 - 1 - len(ehrlich22)       # 20 rate constants - k_10 - 4 clipped
    PASS('preset trial 0: role-based bands (capacity 1e-3x, inhibition/K_* 0.1x); k_10 excluded (55 sampled, no probe; re-included on 0.1x-10x with exclude_params=()); Ehrlich rates clipped to exactly 1e-3 x b_B, every other baseline unchanged')
else:
    print('SKIP 22: parameter-distribution workbooks not found')

#%% 23. supervisor: preset flags forwarded to the driver; naming mirrors the engine; legacy switch
sup23 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
# Naming: preset convention when a target-products preset is in force,
# the legacy convention otherwise (what --legacy-flags selects). Every
# preset-derived name carries the EFFECTIVE rate band tag _rb{lo}-{hi}
# (the preset's DEFAULT_RATE_MULTIPLIER_BOUNDS unless overridden; always
# tagged since the band moved 1e-5x -> 1e-3x on 2026-09-06, so a
# new-band study never resumes the 1e-5x study of the same objective,
# e.g. kin_opt_ethanol_isobutanol_metabolic_irr_ib0.1-10_burden) and the
# inhibition-coefficient band tag _ib{lo}-{hi} (the supervisor passes
# the presets' K_* band; since 2026-09-06, when the inhibition
# coefficients k_*i* left the k_* rate band), and the exclusion tag
# _xk10 (the presets' DEFAULT_EXCLUDED_PARAMETERS; since 2026-09-06 pm) and
# the operating-variable band tag _s1x1-50 (DEFAULT_STAGE_1_MAX_X_BOUNDS;
# since 2026-09-06 pm, check 40).
assert sup23['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50'
assert sup23['default_study_name']('A', 'IBO titer', 'B',
                                   study_target_products='ethanol_only',
                                   study_type='metabolic') \
    == 'kin_opt_ethanol_only_metabolic_ibo_titer_kbB_rb0.001-10_ib0.1-10_xk10_s1x1-50'
    # ethanol_only's own kinetic_bounds_scenario is 'A' (STUDY_TARGET_PRODUCTS);
    # the explicit 'B' here differs, so it IS tagged (Finding 1 fix -- this
    # name used to silently drop the override and collide with the default).
assert sup23['default_study_name']('A', 'IRR', 'B', study_target_products=None,
                                   study_type='metabolic') == 'kin_opt_A_kbB_irr'
# Legacy result (study_target_products=None, positional) is untouched.
assert sup23['default_study_name']('A', 'IRR', 'B') == 'kin_opt_A_kbB_irr'
# Mirror of the engine's tag rule (Finding 1): only a differing override
# is tagged; matching the preset's own scenario / kinetic_bounds_scenario
# is silent.
assert sup23['default_study_name']('B', 'IRR', None,
                                   'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_scB_rb0.001-10_ib0.1-10_xk10_s1x1-50'
assert sup23['default_study_name'](None, 'IRR', 'A',
                                   'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_kbA_rb0.001-10_ib0.1-10_xk10_s1x1-50'
assert sup23['default_study_name']('A', 'IRR', 'B',
                                   'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50'  # both match the preset: no sc/kb tag
# child_code forwards both kwargs (and None under --legacy-flags).
code23 = sup23['child_code'](None, 'IRR', 2000, None, False, 'x',
                             study_target_products='ethanol_only',
                             study_type='metabolic')
assert "study_target_products='ethanol_only'" in code23
assert "study_type='metabolic'" in code23 and 'scenario=None' in code23
code23_legacy = sup23['child_code']('A', 'IRR', 2000, 'B', False, 'kin_opt_A_kbB_irr',
                                    study_target_products=None, study_type='metabolic_protein')
assert 'study_target_products=None' in code23_legacy and "scenario='A'" in code23_legacy
# Defaults: supervise() and child_code() agree with the engine's constants;
# the supervisor's scenario default is None (a preset picks it).
_p23 = _inspect.signature(sup23['supervise']).parameters
assert _p23['study_target_products'].default == ko.DEFAULT_STUDY_TARGET_PRODUCTS
assert _p23['study_type'].default == ko.DEFAULT_STUDY_TYPE
assert _p23['scenario'].default is None
# Stall timeout default 3 min (2026-09-06; 25 before) in supervise() and
# the CLI (source-text check, the parser lives under __main__).
assert _p23['stall_timeout_min'].default == 3.0
_c23 = _inspect.signature(sup23['child_code']).parameters
assert 'study_target_products' in _c23 and 'study_type' in _c23
# The argparse choices/defaults are derived from the engine's tables, not
# re-typed as literals (Finding 2), even though the parser itself is built
# only under `if __name__ == '__main__':` (runpy.run_path's default
# run_name is not '__main__', so sup23 has no live parser/args to
# introspect -- these are pinned as source-text checks instead).
src23 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert '--study-target-products' in src23 and '--study-type' in src23
assert '--legacy-flags' in src23
assert 'choices=tuple(ko.STUDY_TARGET_PRODUCTS)' in src23
assert 'choices=tuple(ko.STUDY_TYPE_ROLES)' in src23
assert "'--study-target-products', default=ko.DEFAULT_STUDY_TARGET_PRODUCTS" in src23
assert "'--study-type', default=ko.DEFAULT_STUDY_TYPE" in src23
assert "'--stall-timeout-min', type=float, default=3.0" in src23
assert 'import biorefineries' not in src23 and 'import optuna' not in src23
# Finding 3: --no-restrict-to-workbook without --legacy-flags must fail at
# argparse time (parser.error), not after a full child load.
assert 'args = parser.parse_args()' in src23
assert 'not args.restrict_to_workbook' in src23 and 'not args.legacy_flags' in src23
assert "parser.error('--no-restrict-to-workbook requires --legacy-flags" in src23
# supervise() derives the SAME csv path the child will write: the preset
# name, from the study_name computed before any child launches.
assert 'default_study_name(scenario, objective,' in _inspect.getsource(sup23['supervise'])
assert 'study_target_products=study_target_products' in _inspect.getsource(sup23['supervise'])
# The engine kwarg the driver forwards is accepted by the real engine with
# fake handles: the k_* band reaches optuna's distributions (optuna only).
if _optuna is not None:
    outdir23 = tempfile.mkdtemp()
    st23 = {'irr': 0.2}
    def _model_specification23(**kw):
        pass
    def _solve_TEA23(stream_IDs=None):
        return {'IRR': st23['irr'], 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles23 = dict(handles17, model_specification=_model_specification23,
                     solve_TEA=_solve_TEA23,
                     latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    study23, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=1, seed=1,
        study_name='offline_rate_band', results_dir=outdir23, handles=handles23,
        rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1,
        burden_model=None)
    d23 = study23.trials[0].distributions
    assert np.isclose(d23['k_1e'].low, 1e-5*47.1) and np.isclose(d23['k_1e'].high, 10.0*47.1)
    assert np.isclose(d23['K_1e'].low, 0.1*0.12) and np.isclose(d23['K_1e'].high, 10.0*0.12)
    assert d23['k_1e'].log and d23['K_1e'].log
    # trial 0 = the baseline (in-bounds, so unclipped)
    assert np.isclose(study23.trials[0].params['k_1e'], 47.1)
else:
    print('SKIP 23b: optuna not installed')
# The supervisor's own recover_inflight wrapper never lets the trajectory
# header guard crash the supervisor (a sidecar whose columns mismatch the
# CSV, e.g. a burden sidecar next to a burden-free CSV): it warns with the
# FULL sidecar record (so the lost trial survives in the supervisor log),
# clears the sidecar (so the next start does not hit it again) and
# recovers nothing. supervise() routes both recoveries through it.
import io as _io
import contextlib as _contextlib
outdir23s = tempfile.mkdtemp()
csv23s = os.path.join(outdir23s, 's_trajectory.csv')
side23s = ko.inflight_path_for(outdir23s, 's')
ko.append_trajectory_row(csv23s, columns, rec)
ko.write_inflight(side23s, [*columns, 'Phi_M'],
                  dict(rec, trial_number=7, Phi_M=0.31))
buf23 = _io.StringIO()
with _contextlib.redirect_stdout(buf23):
    got23 = sup23['recover_inflight'](csv23s, side23s, state='LOST', error='x')
assert got23 is None
assert not os.path.isfile(side23s)
out23 = buf23.getvalue()
assert 'different column set' in out23, out23
assert "'trial_number': 7" in out23 and "'Phi_M': 0.31" in out23, out23
assert sup23['row_count'](csv23s) == 1                       # nothing appended
assert ko.load_trajectory(csv23s)['state'].tolist() == ['COMPLETE']
src23_sup = _inspect.getsource(sup23['supervise'])
assert 'ko.recover_inflight(' not in src23_sup and src23_sup.count('recover_inflight(') == 2
# and the wrapper is transparent when the header matches
ko.write_inflight(side23s, columns, dict(rec, trial_number=8))
assert sup23['recover_inflight'](csv23s, side23s, state='LOST', error='y') == 8
assert ko.load_trajectory(csv23s)['state'].tolist() == ['COMPLETE', 'LOST']
PASS('supervisor: preset flags/defaults forwarded, naming mirrors the engine, legacy switch; k_* band reaches optuna; header-mismatch recovery warns, clears, never raises')

#%% 24. burden plumbing: _burden suffix, extra trajectory columns, header guard, INFEASIBLE colour
from biorefineries.isobutanol import enzyme_burden as eb
assert ko.BURDEN_STUDY_SUFFIX == '_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'            # default: unchanged
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_burden'
assert ko.default_study_name('IBO titer', 'ethanol_only', 'metabolic',
                             scenario='B', kinetic_bounds_scenario='B', burden=True) \
    == 'kin_opt_ethanol_only_metabolic_ibo_titer_scB_kbB_burden'    # suffix goes last
cols24 = ko.trajectory_columns(space, extra_columns=eb.BURDEN_COLUMNS)
assert cols24 == ['trial_number', 'state', *space, 'objective',
                  *ko.TRACKED_METRICS, *eb.BURDEN_COLUMNS, 'error']
assert ko.trajectory_columns(space) == ko.trajectory_columns(space, extra_columns=())
# The header guard keeps a burden run out of a burden-free CSV and vice versa.
csv24 = os.path.join(outdir, 'burden_guard.csv')
ko.append_trajectory_row(csv24, cols24, {'trial_number': 0, 'state': 'COMPLETE',
                                         'objective': 0.1, 'Phi_M': 0.0637})
try:
    ko.append_trajectory_row(csv24, ko.trajectory_columns(space),
                             {'trial_number': 1, 'state': 'COMPLETE'})
except ValueError as e:
    assert 'different column set' in str(e)
else:
    raise AssertionError('header guard did not refuse the burden-free column set')
# The same guard as a standalone pre-flight check (run_kinetic_optimization
# calls it before creating the optuna study, so a study name colliding with
# a store of a different column set fails before any sidecar or
# simulation): raises on a mismatch, passes on a match and on an absent CSV.
assert 'check_trajectory_header' in ko.__all__
ko.check_trajectory_header(csv24, cols24)
ko.check_trajectory_header(os.path.join(outdir, 'no_such_trajectory.csv'), cols24)
try:
    ko.check_trajectory_header(csv24, ko.trajectory_columns(space))
except ValueError as e:
    assert 'different column set' in str(e)
else:
    raise AssertionError('check_trajectory_header did not refuse the burden-free column set')
# PCA decision columns are still the search-space variables (burden columns sit after objective)
df24 = ko.load_trajectory(csv24)
assert list(df24.columns) == cols24
assert list(df24.columns)[list(df24.columns).index('state') + 1:
                          list(df24.columns).index('objective')] == list(space)
# INFEASIBLE rows: excluded from completed-only plots, drawn in the PCA landscape
synth24 = synth11.copy()
inf24 = synth24['trial_number'].isin([5, 13, 21])
synth24.loc[inf24, 'state'] = 'INFEASIBLE'
synth24.loc[inf24, ['objective', *ko.TRACKED_METRICS]] = np.nan
synth24.loc[inf24, 'error'] = 'enzyme burden: Phi_M 0.2862 > F_flex 0.2450 g/gDCW'
ok24 = ko._completed(synth24)
assert not ok24['trial_number'].isin([5, 13, 21]).any()
f24 = [os.path.join(outdir, f'infeasible_{i}.png') for i in range(4)]
ko.plot_optimization_trajectories(synth24, objective_name='IRR',
                                  direction='maximize', filename=f24[0])
ko.plot_parameter_trajectory(synth24, baselines, direction='maximize', filename=f24[1])
ko.plot_best_vs_baseline(synth24, baselines, direction='maximize', filename=f24[2])
_, axes24 = ko.plot_pca_projection(synth24, 'maximize', log_columns=log_cols,
                                   objective_name='IRR', filename=f24[3])
labels24 = axes24[0].get_legend_handles_labels()[1]
assert 'infeasible (enzyme burden, pruned)' in labels24, labels24
assert 'completed' in labels24
for p in f24:
    assert os.path.isfile(p) and os.path.getsize(p) > 0, p
_, axes24b = ko.plot_pca_projection(synth11, 'maximize', log_columns=log_cols)
assert 'infeasible (enzyme burden, pruned)' not in axes24b[0].get_legend_handles_labels()[1]
PASS('burden plumbing: _burden suffix, BURDEN_COLUMNS after the metrics, header guard, INFEASIBLE in the PCA legend')

#%% 25. engine hook: INFEASIBLE pruned before the sidecar/simulation; effective k_7 written, sampled k_7 recorded
if _optuna is None:
    print('SKIP 25: optuna not installed')
else:
    outdir25 = tempfile.mkdtemp()
    study25 = 'offline_burden'
    csv25 = os.path.join(outdir25, study25 + '_trajectory.csv')
    side25 = ko.inflight_path_for(outdir25, study25)

    class _FakeTE25:
        # the model's reference capacities (scenario A: Ehrlich off)
        k_1h = 0.584; k_1l = 1.43; k_1e = 47.1; k_2 = 0.501; k_3 = 5.81
        k_4 = 4.8; k_5 = 0.0104; k_5e = 0.775; k_6 = 2.82
        k_7 = 1.203; k_8 = 0.589
        k_13 = 0.0; k_14 = 0.0; k_15 = 0.0; k_16 = 0.0
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1h', 'k_1l', 'k_1e', 'k_2', 'k_3', 'k_4', 'k_5',
                    'k_5e', 'k_6', 'k_7', 'k_8', 'k_13', 'k_14', 'k_15',
                    'k_16', 'K_1e', 'not_kinetic']
    te25 = _FakeTE25()
    seen_k7 = []          # k_7 on the fake model at each model_specification call
    def _model_specification25(**kw):
        seen_k7.append(te25.k_7)
    def _solve_TEA25(stream_IDs=None):
        return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles25 = dict(handles17, r_te=te25,
                     model_specification=_model_specification25,
                     solve_TEA=_solve_TEA25,
                     latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    baselines25 = ko.discover_kinetic_parameters(te25)
    override25 = {'k_13': (0.0, 60.0), 'k_14': (0.0, 50.0),
                  'k_15': (0.0, 50.0), 'k_16': (0.0, 30.0)}
    space25, _ = ko.build_search_space(baselines25, param_bounds_override=override25,
                                       rate_multiplier_bounds=(1e-5, 10.0))
    base25 = ko.baseline_decision_point(
        space25, baselines25, fbs17.current_specifications, fbs17.max_n_spikes)
    # Pre-enqueue three scripted trials in the study the engine will resume
    # (same storage URL as the engine builds): 0 = infeasible (k_13 = 60 ->
    # Phi_M ~ 0.40 > F_flex), 1 = feasible but derated (9x k_7 at wild-type
    # enzymes), 2 = the exact reference (inert). study.trials counts the
    # WAITING trials, so n_trials=6 runs exactly these three, FIFO.
    storage25 = ('sqlite:///' + os.path.join(outdir25, study25 + '.db')
                 .replace('\\', '/'))
    pre25 = _optuna.create_study(study_name=study25, storage=storage25,
                                 direction='maximize')
    pre25.enqueue_trial({**base25, 'k_13': 60.0})
    pre25.enqueue_trial({**base25, 'k_7': 9.0*1.203})
    pre25.enqueue_trial(dict(base25))
    n_sidecar25 = []
    _orig_write_inflight = ko.write_inflight
    def _counting_write_inflight(path, columns, record):
        n_sidecar25.append(record['trial_number'])
        return _orig_write_inflight(path, columns, record)
    ko.write_inflight = _counting_write_inflight
    try:
        study25_obj, csv25_out, kb25 = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=6, seed=1,
            study_name=study25, results_dir=outdir25, handles=handles25,
            param_bounds_override=override25,
            rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1)
    finally:
        ko.write_inflight = _orig_write_inflight
    assert csv25_out == csv25 and kb25 == baselines25
    df25 = ko.load_trajectory(csv25)
    assert list(df25.columns) == ko.trajectory_columns(space25, extra_columns=eb.BURDEN_COLUMNS)
    assert df25['trial_number'].tolist() == [0, 1, 2]
    assert df25['state'].tolist() == ['INFEASIBLE', 'COMPLETE', 'COMPLETE']
    # trial 0: burden columns recorded, no objective, no sidecar, no simulation
    assert df25['k_13'][0] == 60.0 and df25['burden_factor'][0] == 0.0
    assert df25['Phi_M'][0] > df25['F_flex'][0] == eb.F_FLEX == 0.245
    assert np.isnan(df25['objective'][0]) and np.isnan(df25['IRR'][0])
    assert df25['error'][0].startswith('enzyme burden: Phi_M ')
    assert n_sidecar25 == [1, 2], n_sidecar25
    # trial 1: the CSV keeps the SAMPLED k_7; the model now receives the
    # INTENDED (sampled) k_7 too -- the k_7/k_8 derating is delegated to the
    # load_simulate choke point (system._apply_enzyme_burden), which this
    # offline no-load harness does not exercise. The burden columns
    # (k_7_eff etc.) are still recorded by the pre-sim evaluate.
    assert np.isclose(df25['k_7'][1], 9.0*1.203)
    assert np.isclose(df25['burden_factor'][1], 0.17701, rtol=1e-3)   # 0.13276 at TRANSLATION_FRACTION_WT = 0.30
    assert np.isclose(df25['k_7_eff'][1], 1.9165, rtol=1e-3)          # 1.4374 at 0.30
    assert np.isclose(df25['k_8_eff'][1], 0.17701*0.589, rtol=1e-3)
    # Phi_M,wt = the report's 0.0637 (at P = 0.45) x PROTEIN_CONTENT/0.45
    assert np.isclose(df25['Phi_M'][1], 0.0637*eb.PROTEIN_CONTENT/eb.POOL_TABLE_PROTEIN_CONTENT)
    assert np.isclose(df25['phi_T'][1], 9.0*eb.PHI_T_WT)
    assert df25['objective'][1] == 0.2
    # trial 2: the reference is inert
    assert df25['burden_factor'][2] == 1.0 and df25['k_7_eff'][2] == 1.203
    assert np.isclose(df25['pool_r1'][2], eb.NATIVE_STEPS['r1'][0]) and df25['pool_r13'][2] == 0.0
    # model_specification saw trial 1's INTENDED (sampled) k_7 -- the
    # optimizer no longer derates it in-place; that happens at the
    # load_simulate choke point instead -- then trial 2's reference k_7,
    # then restore_baseline's reference k_7 (2 simulations + the finally)
    assert len(seen_k7) == 3, seen_k7
    assert np.isclose(seen_k7[0], 9.0*1.203) and seen_k7[1] == 1.203 and seen_k7[2] == 1.203
    assert te25.k_7 == 1.203 and te25.k_8 == 0.589                  # restored
    assert not os.path.isfile(side25)
    # optuna side: the infeasible trial is PRUNED, its violation reached the
    # sampler constraint (system attr 'constraints'); feasible trials <= 0
    t25 = study25_obj.trials
    TS = _optuna.trial.TrialState
    assert [t.state for t in t25] == [TS.PRUNED, TS.COMPLETE, TS.COMPLETE]
    assert t25[0].user_attrs['burden_violation'] > 0.0
    assert t25[1].user_attrs['burden_violation'] < 0.0 and t25[2].user_attrs['burden_violation'] < 0.0
    assert list(t25[0].system_attrs['constraints'])[0] > 0.0, t25[0].system_attrs
    assert list(t25[1].system_attrs['constraints'])[0] < 0.0
    # Pre-flight header guard (F1): resuming the burden study's name / CSV
    # with the burden OFF (a different column set, e.g. --study-name of a
    # legacy study without --no-burden, or vice versa) fails BEFORE any
    # sidecar is written or any simulation runs -- not at the first row
    # append after a full ~20 s simulation.
    n_sidecar25_before = list(n_sidecar25)
    n_seen25_before = len(seen_k7)
    ko.write_inflight = _counting_write_inflight
    try:
        ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=6, seed=1,
            study_name=study25, results_dir=outdir25, handles=handles25,
            param_bounds_override=override25,
            rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1,
            burden_model=None)
    except ValueError as e:
        assert 'different column set' in str(e), e
    else:
        raise AssertionError('burden-free resume of the burden study was not refused')
    finally:
        ko.write_inflight = _orig_write_inflight
    assert n_sidecar25 == n_sidecar25_before and len(seen_k7) == n_seen25_before
    assert not os.path.isfile(side25)
    assert ko.load_trajectory(csv25)['trial_number'].tolist() == [0, 1, 2]   # untouched
    # A caller-supplied BurdenModel (F2) must be a snapshot of the LIVE
    # baselines: one perturbed reference capacity is named with both values.
    bm25_stale = eb.BurdenModel.from_reference({**baselines25, 'k_3': 6.0})
    try:
        ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=1, seed=1,
            study_name='offline_stale_reference', results_dir=tempfile.mkdtemp(),
            handles=handles25, param_bounds_override=override25,
            rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1,
            burden_model=bm25_stale)
    except ValueError as e:
        assert 'k_3' in str(e) and '6.0' in str(e) and '5.81' in str(e), e
    else:
        raise AssertionError('a stale BurdenModel reference was not refused')
    # ... and a burden_model that is neither 'auto', None nor a model is a TypeError
    try:
        ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=1, seed=1,
            study_name='offline_bad_burden', results_dir=tempfile.mkdtemp(),
            handles=handles25, param_bounds_override=override25,
            rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1,
            burden_model=True)
    except TypeError as e:
        assert 'burden_model' in str(e), e
    else:
        raise AssertionError('burden_model=True was not refused')
    assert len(seen_k7) == n_seen25_before                     # no simulation in either
    # Offline resume of the burden study (F3): same name / sqlite store /
    # CSV, a budget of ONE more trial (the store holds 3, so 4). The store
    # holds only 3 trials -- fewer than TPE's n_startup_trials -- so a
    # freely-sampled resumed trial would just be a random startup draw;
    # enqueue a known-feasible point on the SAME study object (pre25,
    # loaded from the same storage URL) so the resumed trial is forced to
    # that exact point and deterministically exercises the sidecar ->
    # simulate -> COMPLETE path with burden columns on a resumed study:
    # exactly one new row (no header error), the new trial COMPLETE and
    # feasible, carrying the sampler constraint, and the baseline restored
    # afterwards.
    import warnings as _warnings
    n_rows25_before = len(ko.load_trajectory(csv25))
    n_trials25_before = len(study25_obj.trials)
    assert n_trials25_before == 3
    pre25.enqueue_trial(dict(base25))          # exact reference: feasible, inert
    seen_k7.clear()
    # Belt-and-braces guard, not the proof: optuna's "does not have
    # constraint values" warning only fires when the sampler actually
    # samples (rather than dequeues) over stored trials lacking
    # constraints, so it cannot fire here regardless of scenario. The
    # property it would guard -- every stored trial carries a sampler
    # constraint -- is verified directly below via the trials'
    # system_attrs['constraints'].
    with _warnings.catch_warnings(record=True) as w25:
        _warnings.simplefilter('always')
        # +2, not +1: the enqueue above already added a WAITING trial to
        # storage, so the store holds 4 by the time the engine reads it;
        # a budget of n_trials25_before + 2 (5) is what makes the engine
        # run exactly that one already-queued trial (5 - 4 stored = 1).
        study25r, csv25r, kb25r = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=n_trials25_before + 2,
            seed=1, study_name=study25, results_dir=outdir25, handles=handles25,
            param_bounds_override=override25,
            rate_multiplier_bounds=(1e-5, 10.0), print_status_every=1)
    assert csv25r == csv25 and kb25r == baselines25
    df25r = ko.load_trajectory(csv25)
    assert len(df25r) == n_rows25_before + 1, len(df25r)
    assert df25r['trial_number'].tolist() == [0, 1, 2, 3]
    assert df25r['state'].iloc[-1] == 'COMPLETE', df25r['state'].tolist()
    assert list(df25r.columns) == ko.trajectory_columns(space25, extra_columns=eb.BURDEN_COLUMNS)
    t25r = study25r.trials
    assert len(t25r) == n_trials25_before + 1
    assert t25r[-1].state == TS.COMPLETE, t25r[-1].state
    assert 'constraints' in t25r[-1].system_attrs, t25r[-1].system_attrs
    assert list(t25r[-1].system_attrs['constraints'])[0] <= 0.0    # feasible
    assert 'burden_violation' in t25r[-1].user_attrs
    bad25 = [str(x.message) for x in w25 if 'does not have constraint values' in str(x.message)]
    assert not bad25, bad25
    assert te25.k_7 == 1.203 and te25.k_8 == 0.589                  # restored
    assert seen_k7[0] == 1.203                                      # the resumed trial itself
    assert seen_k7[-1] == 1.203                                     # the restore reached the model
    assert not os.path.isfile(side25)
    # burden off: no burden columns, k_7 written as sampled
    outdir25b = tempfile.mkdtemp()
    seen_k7.clear()
    study25b_obj, csv25b, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=1, seed=1,
        study_name='offline_no_burden', results_dir=outdir25b, handles=handles25,
        param_bounds_override=override25, rate_multiplier_bounds=(1e-5, 10.0),
        print_status_every=1, burden_model=None)
    df25b = ko.load_trajectory(csv25b)
    assert 'Phi_M' not in df25b.columns and df25b['state'].tolist() == ['COMPLETE']
    assert 'burden_violation' not in study25b_obj.trials[0].user_attrs
    # the engine's fallback name carries the suffix only when the burden is on
    assert os.path.isfile(os.path.join(outdir25b, 'offline_no_burden_trajectory.csv'))
    outdir25c = tempfile.mkdtemp()
    _, csv25c, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=1, seed=1,
        results_dir=outdir25c, handles=handles25,
        param_bounds_override=override25, rate_multiplier_bounds=(1e-5, 10.0),
        print_status_every=1)
    assert csv25c == os.path.join(outdir25c, 'kin_opt_X_irr_burden_trajectory.csv'), csv25c
    PASS('engine hook: INFEASIBLE pruned pre-sidecar with burden columns; effective k_7 to the model, sampled k_7 in the CSV; constraint reaches optuna; pre-flight header guard, stale/invalid burden_model refused, resume clean; burden off unchanged')

#%% 26. driver / supervisor: burden default on, --no-burden, naming mirrors the engine, reports printed
sup26 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
# Naming: the suffix on both paths, only when asked; legacy defaults untouched.
assert sup26['default_study_name']('A', 'IRR', 'B') == 'kin_opt_A_kbB_irr'
assert sup26['default_study_name']('A', 'IRR', 'B', burden=True) == 'kin_opt_A_kbB_irr_burden'
assert sup26['default_study_name'](None, 'IRR', None, burden=True) == 'kin_opt_B_irr_burden'
assert sup26['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
assert sup26['default_study_name']('B', 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_scB_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
# child_code forwards the flag both ways; supervise()'s default is ON and it
# derives the study name WITH the flag (so resume/stall-kill hit the store
# the child writes).
code26 = sup26['child_code'](None, 'IRR', 2000, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic_protein', burden=False)
assert 'burden=False' in code26
code26b = sup26['child_code'](None, 'IRR', 2000, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic_protein')
assert 'burden=True' in code26b
_p26 = _inspect.signature(sup26['supervise']).parameters
assert _p26['burden'].default is True
_c26 = _inspect.signature(sup26['child_code']).parameters
assert _c26['burden'].default is True
_n26 = _inspect.signature(sup26['default_study_name']).parameters
assert _n26['burden'].default is False
src26_sup = _inspect.getsource(sup26['supervise'])
assert 'burden=burden' in src26_sup                       # into default_study_name and child_code
src26 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--no-burden'" in src26 and "dest='burden'" in src26
assert "action='store_false'" in src26 and 'burden=args.burden' in src26
assert 'import biorefineries' not in src26 and 'import optuna' not in src26
assert 'enzyme_burden' not in src26.replace('enzyme_burden.py', '')   # stdlib-only: never imports it
# The driver cannot be imported offline (it loads the biorefinery): pin its
# plumbing as source text.
drv26 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'burden=True,' in drv26                                    # run() kwarg, default on
assert 'burden_model=burden_model' in drv26                       # forwarded to the engine
assert 'ko.default_study_name(' in drv26 and 'burden=burden' in drv26
assert 'ko.BURDEN_STUDY_SUFFIX if burden else' in drv26           # legacy _kb path
assert 'scenarios.load_scenario(scenario, burden=burden)' in drv26   # A-referenced burden from the bundle
assert "bundle['burden_model']" in drv26                             # consumed, not rebuilt in the driver
assert 'eb.scenario_b_ehrlich()' in drv26 and 'describe_point(' in drv26
assert "'burden_model' in engine_kwargs" in drv26                 # ambiguity guard
PASS('driver/supervisor: burden on by default, --no-burden, _burden naming on both paths, A/B reports wired')

#%% 27. single-knockout probes: pure points; engine enqueues them right after trial 0 (fresh study only); _rb naming tag; driver/supervisor plumbing
# Pure: one probe per log-scale RATE constant (k_*) of the space, at that
# variable's floor, every other decision variable at the baseline point.
# K_*, the feeding variables and a rate whose baseline already sits at its
# floor (the preset's clipped Ehrlich rates: a probe would duplicate trial
# 0) get none.
kb27 = {'k_1e': 47.1, 'k_7': 1.203, 'k_13': 0.0, 'K_1e': 0.12}
override27 = {'k_13': (1e-5*5.81, 10.0*5.81)}   # the preset's Ehrlich band
space27, excl27 = ko.build_search_space(kb27, param_bounds_override=override27,
                                        rate_multiplier_bounds=(0.1, 10.0))
assert excl27 == [] and space27['k_1e']['low'] == 4.71
base27 = ko.baseline_decision_point(space27, kb27,
                                    fbs17.current_specifications, 16)
assert base27['k_13'] == 1e-5*5.81                       # clipped to the floor
probes27, at_floor27 = ko.knockout_probe_points(space27, base27)
assert list(probes27) == ['k_1e', 'k_7'], list(probes27)   # space order; K_1e, feeding, k_13 skipped
assert at_floor27 == ['k_13']
for n27, p27 in probes27.items():
    assert p27 is not base27 and p27[n27] == space27[n27]['low']
    assert {k: v for k, v in p27.items() if k != n27} \
        == {k: v for k, v in base27.items() if k != n27}
assert base27['k_1e'] == 47.1 and base27['k_7'] == 1.203  # inputs untouched
# A space without rate constants yields no probes; a linear (override with
# low <= 0) rate is not a log-scale variable and gets none either.
assert ko.knockout_probe_points({'K_1e': space27['K_1e']}, {'K_1e': 0.12}) == ({}, [])
space27b, _ = ko.build_search_space({'k_13': 0.0},
                                    param_bounds_override={'k_13': (0.0, 58.1)})
assert ko.knockout_probe_points(space27b, {'k_13': 0.0}) == ({}, [])

if _optuna is None:
    print('SKIP 27 (engine part): optuna not installed')
else:
    class _FakeTE27:
        k_1e = 47.1; k_7 = 1.203; k_13 = 0.0; K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'k_7', 'k_13', 'K_1e', 'not_kinetic']
    te27 = _FakeTE27()
    seen27 = []      # (k_1e, k_7, k_13) on the fake model at each simulation
    def _model_specification27(**kw):
        seen27.append((te27.k_1e, te27.k_7, te27.k_13))
    def _solve_TEA27(stream_IDs=None):
        return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    def _handles27():
        return dict(handles17, r_te=te27,
                    model_specification=_model_specification27,
                    solve_TEA=_solve_TEA27,
                    latest_TEA_solution={'IRR': np.nan,
                                         'MPSPs': {'ethanol': np.nan,
                                                   'isobutanol': np.nan}})
    outdir27 = tempfile.mkdtemp()
    study27 = 'offline_knockouts'
    # Opted in (both enqueue defaults are OFF since 2026-09-07): trial 0 =
    # baseline, trials 1-2 = the k_1e / k_7 probes
    # (k_13 sits at its floor -> no probe), trial 3 = the first sampled
    # point. The engine's own baseline point is what the probes copy.
    st27, csv27, kb27_out = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=4, seed=1,
        study_name=study27, results_dir=outdir27, handles=_handles27(),
        param_bounds_override=override27, rate_multiplier_bounds=(0.1, 10.0),
        print_status_every=1, burden_model=None)
    df27 = ko.load_trajectory(csv27)
    assert df27['trial_number'].tolist() == [0, 1, 2, 3]
    assert df27['state'].tolist() == ['COMPLETE']*4
    assert np.isclose(df27['k_1e'][0], 47.1) and np.isclose(df27['k_7'][0], 1.203)
    assert np.isclose(df27['k_1e'][1], 4.71) and np.isclose(df27['k_7'][1], 1.203)
    assert np.isclose(df27['k_1e'][2], 47.1) and np.isclose(df27['k_7'][2], 0.1203)
    assert np.allclose(df27['k_13'][:3], 1e-5*5.81)          # floor-clipped, never probed
    for c27 in ('K_1e', *ko.FEEDING_VARIABLES):               # everything else pinned at baseline
        assert len(set(df27[c27][:3].round(12))) == 1, c27
    # The model actually received the knocked-down capacities.
    assert np.isclose(seen27[1][0], 4.71) and np.isclose(seen27[1][1], 1.203)
    assert np.isclose(seen27[2][0], 47.1) and np.isclose(seen27[2][1], 0.1203)
    # The probes are identifiable in the store (user attr), trial 0 / 3 are not.
    tr27 = sorted(st27.trials, key=lambda t: t.number)
    assert [t.user_attrs.get('knockout_probe') for t in tr27] \
        == [None, 'k_1e', 'k_7', None]
    # Resume: no re-enqueue (nothing WAITING, no duplicate probes).
    st27r, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=5, seed=1,
        study_name=study27, results_dir=outdir27, handles=_handles27(),
        param_bounds_override=override27, rate_multiplier_bounds=(0.1, 10.0),
        print_status_every=1, burden_model=None)
    assert len(st27r.trials) == 5
    assert sum(1 for t in st27r.trials if t.user_attrs.get('knockout_probe')) == 2
    # OFF: trial 1 is a sampled point, not a probe.
    outdir27b = tempfile.mkdtemp()
    st27b, csv27b, _ = ko.run_kinetic_optimization(enqueue_baseline=True,
        objective='IRR', scenario_label='X', n_trials=2, seed=1,
        study_name=study27, results_dir=outdir27b, handles=_handles27(),
        param_bounds_override=override27, rate_multiplier_bounds=(0.1, 10.0),
        enqueue_knockouts=False, print_status_every=1, burden_model=None)
    df27b = ko.load_trajectory(csv27b)
    assert not np.isclose(df27b['k_1e'][1], 4.71) or not np.isclose(df27b['k_7'][1], 1.203)
    assert all(t.user_attrs.get('knockout_probe') is None for t in st27b.trials)
    _e27 = _inspect.signature(ko.run_kinetic_optimization).parameters
    assert _e27['enqueue_knockouts'].default is False   # off by default since 2026-09-07

# Naming: the rate band is tagged _rb{lo}-{hi} WHENEVER it is given (the
# preset's own band included -- since the default moved 1e-5x -> 1e-3x a
# new default-band study must never resume the untagged 1e-5x study of
# the same objective: same columns, so the header guard cannot tell them
# apart); None leaves the name unchanged (older callers).
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic',
                             burden=True, rate_multiplier_bounds=(0.1, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.1-10_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic',
                             burden=True,
                             rate_multiplier_bounds=ko.DEFAULT_RATE_MULTIPLIER_BOUNDS) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic',
                             burden=True, rate_multiplier_bounds=(1e-5, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb1e-05-10_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_burden'
assert ko.default_study_name('IBO titer', 'ethanol_only', 'metabolic_protein',
                             scenario='B', rate_multiplier_bounds=(0.01, 10.0)) \
    == 'kin_opt_ethanol_only_metabolic_protein_ibo_titer_scB_rb0.01-10'
# Supervisor: mirrors the tag, exposes both flags, forwards them to the driver.
sup27 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert sup27['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic', burden=True,
                                   rate_multiplier_bounds=(0.1, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.1-10_ib0.1-10_xk10_s1x1-50_burden'
assert sup27['default_study_name']('A', 'IRR', 'B', rate_multiplier_bounds=(0.1, 10.0)) \
    == 'kin_opt_A_kbB_irr'                                    # legacy path: band not encoded
code27 = sup27['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic',
                             rate_multiplier_bounds=(0.1, 10.0),
                             enqueue_knockouts=True)           # opt in (default off)
assert 'enqueue_knockouts=True' in code27 and 'rate_multiplier_bounds=(0.1, 10.0)' in code27
code27b = sup27['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', enqueue_knockouts=False)
# No explicit band -> the kwarg is OMITTED (passing None would defeat the
# driver's engine_kwargs.setdefault of the preset band).
assert 'enqueue_knockouts=False' in code27b and 'rate_multiplier_bounds' not in code27b
_s27 = _inspect.signature(sup27['supervise']).parameters
assert _s27['enqueue_knockouts'].default is False and _s27['rate_multiplier_bounds'].default is None
src27_sup = _inspect.getsource(sup27['supervise'])
assert 'rate_multiplier_bounds=rate_multiplier_bounds' in src27_sup
assert 'enqueue_knockouts=enqueue_knockouts' in src27_sup
src27 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--enqueue-knockouts'" in src27 and "dest='enqueue_knockouts'" in src27   # opt-in since 2026-09-07
# The flag must be store_true: a regression to store_false with the same
# name/dest would flip the default back ON and pass every check above.
_ekp = src27[src27.index("'--enqueue-knockouts'"):][:300]
assert "action='store_true'" in _ekp and "store_false" not in _ekp
assert "'--rate-multiplier-bounds'" in src27 and 'nargs=2' in src27
assert 'enqueue_knockouts=args.enqueue_knockouts' in src27
drv27 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'enqueue_knockouts=False,' in drv27                          # run() kwarg, default off since 2026-09-07
assert 'enqueue_knockouts=enqueue_knockouts' in drv27              # forwarded to the engine
assert "rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds']" in drv27  # naming sees the EFFECTIVE band
PASS('single-knockout probes: floor points per k_*, enqueued after trial 0 on a fresh study only, identifiable, off switch, _rb naming tag, driver/supervisor flags')

#%% 28. role-based bands (2026-09-06): the k_* rate band applies to RATE
# CONSTANTS (role capacity) only; inhibition coefficients (k_*i*:
# product_inhibition, lethality) share the 0.1x-10x band of the K_* terms.
# `rate_params` (the names the rate band applies to) is threaded through
# build_search_space / workbook_kinetic_bounds / knockout_probe_points /
# run_kinetic_optimization; None = the pre-change lowercase-'k_' prefix
# rule, so every legacy/older call is byte-identical. Preset-derived study
# names carry the inhibition band tag _ib{lo}-{hi}.
assert ko.RATE_CONSTANT_ROLES == ('capacity',)
assert ko.INHIBITION_COEFFICIENT_ROLES == ('product_inhibition', 'lethality')
roles28 = {'k_1e': 'capacity', 'k_1ie': 'product_inhibition',
           'k_10ie': 'lethality', 'K_1i': 'substrate_regulation',
           'K_1e': 'affinity', 'k_6r': 'product_self_inhibition',
           'k_7': 'capacity'}
assert ko.rate_constant_names(['k_1e', 'k_1ie', 'k_10ie', 'K_1i', 'K_1e',
                               'k_6r', 'k_7'], roles=roles28) == ['k_1e', 'k_7']
assert ko.rate_constant_names([], roles=roles28) == []
try:
    ko.rate_constant_names(['k_1e', 'k_99'], roles=roles28)
except KeyError as e28:
    assert 'k_99' in str(e28)
else:
    raise AssertionError('unknown name did not raise KeyError')
# build_search_space: rate_params selects the band; prefix rule when None.
kb28 = {'k_1e': 47.1, 'k_1ie': 0.05, 'k_10ie': 0.02, 'K_1i': 2.0,
        'K_1e': 0.12, 'k_7': 1.203}
rp28 = ko.rate_constant_names(kb28, roles=roles28)
space28, excl28 = ko.build_search_space(
    kb28, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-5, 10.0),
    rate_params=rp28)
assert excl28 == []
assert space28['k_1e'] == dict(low=1e-5*47.1, high=10.0*47.1, log=True)
assert space28['k_7'] == dict(low=1e-5*1.203, high=10.0*1.203, log=True)
assert space28['k_1ie'] == dict(low=0.1*0.05, high=10.0*0.05, log=True)    # inhibition coefficient
assert space28['k_10ie'] == dict(low=0.1*0.02, high=10.0*0.02, log=True)   # lethality coefficient
assert space28['K_1i'] == dict(low=0.1*2.0, high=10.0*2.0, log=True)       # regulation term
assert space28['K_1e'] == dict(low=0.1*0.12, high=10.0*0.12, log=True)
space28_legacy, _ = ko.build_search_space(
    kb28, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-5, 10.0))
assert space28_legacy['k_1ie']['low'] == 1e-5*0.05        # prefix rule: unchanged
assert space28_legacy['k_1e'] == space28['k_1e']
assert space28_legacy['K_1i'] == space28['K_1i']
# An empty rate_params (not None) means NO name gets the rate band.
space28_none, _ = ko.build_search_space(
    kb28, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-5, 10.0),
    rate_params=())
assert space28_none['k_1e']['low'] == 0.1*47.1
# Override and whitelist precedence is untouched.
space28b, _ = ko.build_search_space(
    kb28, rate_multiplier_bounds=(1e-5, 10.0), rate_params=rp28,
    param_bounds_override={'k_1e': (1.0, 2.0)}, include_params=['k_1e', 'k_1ie'])
assert space28b['k_1e'] == dict(low=1.0, high=2.0, log=True)
assert set(space28b) == {'k_1e', 'k_1ie', *ko.FEEDING_VARIABLES}
# knockout_probe_points: rate_params restricts the probes; None = prefix.
base28 = ko.baseline_decision_point(space28, kb28,
                                    fbs17.current_specifications, 16)
probes28, at_floor28 = ko.knockout_probe_points(space28, base28,
                                                rate_params=rp28)
assert list(probes28) == ['k_1e', 'k_7'] and at_floor28 == []
assert probes28['k_1e']['k_1e'] == 1e-5*47.1 and probes28['k_1e']['k_1ie'] == 0.05
assert list(ko.knockout_probe_points(space28, base28)[0]) \
    == ['k_1e', 'k_1ie', 'k_10ie', 'k_7']                   # legacy prefix rule
assert ko.knockout_probe_points(space28, base28, rate_params=())[0] == {}
# workbook_kinetic_bounds: absolute role-based bands (plain file read).
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    roles28_real = ko.kinetic_parameter_roles()
    for sc28 in ('A', 'B'):
        base28_wb = ko.workbook_kinetic_baselines(sc28)
        rp28_wb = ko.rate_constant_names(base28_wb, roles=roles28_real)
        assert len(rp28_wb) == (16 if sc28 == 'A' else 20)
        wbb28 = ko.workbook_kinetic_bounds(
            sc28, multiplier_bounds=(0.1, 10.0),
            rate_multiplier_bounds=(1e-5, 10.0), rate_params=rp28_wb)
        assert list(wbb28) == list(base28_wb)
        for n28, b28 in base28_wb.items():
            assert wbb28[n28] == ((1e-5*b28, 10.0*b28) if n28 in rp28_wb
                                  else (0.1*b28, 10.0*b28)), n28
        inh28 = [n for n in base28_wb
                 if roles28_real[n] in ko.INHIBITION_COEFFICIENT_ROLES]
        assert len(inh28) == (9 if sc28 == 'A' else 16)
        assert all(n.startswith('k_') and wbb28[n][0] == 0.1*base28_wb[n]
                   for n in inh28)
        # Legacy prefix rule (rate_params=None) unchanged.
        wbb28_legacy = ko.workbook_kinetic_bounds(
            sc28, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-5, 10.0))
        assert all(wbb28_legacy[n][0] == 1e-5*base28_wb[n] for n in inh28)
else:
    print('SKIP 28 (workbook part): parameter-distribution workbooks not found')
# Naming: inhibition_multiplier_bounds tags _ib{lo}-{hi} whenever given
# (after _rb, before _burden); None leaves every existing name untouched.
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             inhibition_multiplier_bounds=(0.1, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_ib0.1-10'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             inhibition_multiplier_bounds=(0.1, 10.0), burden=True,
                             rate_multiplier_bounds=(0.1, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.1-10_ib0.1-10_burden'
assert ko.default_study_name('IBO titer', 'ethanol_only', 'metabolic',
                             scenario='B', kinetic_bounds_scenario='B',
                             inhibition_multiplier_bounds=(0.01, 10.0)) \
    == 'kin_opt_ethanol_only_metabolic_ibo_titer_scB_kbB_ib0.01-10'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'   # None: untouched
# Engine kwarg plumbed (default None).
_sig28 = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'rate_params' in _sig28 and _sig28['rate_params'].default is None
# Driver / supervisor plumbing (source text: the driver cannot be
# imported offline).
drv28 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert "'rate_params'" in drv28                                   # preset key forwarded
assert "rate_params=engine_kwargs.get('rate_params')" in drv28    # workbook bounds by role
assert 'inhibition_multiplier_bounds=' in drv28                   # preset names tagged
sup28 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert 'inhibition_multiplier_bounds=' in _inspect.getsource(sup28['default_study_name'])
assert sup28['default_study_name']('A', 'IRR', 'B') == 'kin_opt_A_kbB_irr'   # legacy: no tag
if _optuna is None:
    print('SKIP 28 (engine part): optuna not installed')
else:
    # End to end: with rate_params, the engine probes ONLY the rate
    # constants and samples the inhibition coefficient on its 0.1x band.
    class _FakeTE28:
        k_1e = 47.1; k_1ie = 0.05; k_7 = 1.203; K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'k_1ie', 'k_7', 'K_1e']
    te28 = _FakeTE28()
    def _model_specification28(**kw):
        pass
    def _solve_TEA28(stream_IDs=None):
        return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles28 = dict(handles17, r_te=te28,
                     model_specification=_model_specification28,
                     solve_TEA=_solve_TEA28,
                     latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    outdir28 = tempfile.mkdtemp()
    st28, csv28, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=30, seed=1,
        study_name='offline_role_bands', results_dir=outdir28,
        handles=handles28, multiplier_bounds=(0.1, 10.0),
        rate_multiplier_bounds=(1e-5, 10.0), rate_params=['k_1e', 'k_7'],
        print_status_every=10, burden_model=None)
    df28 = ko.load_trajectory(csv28)
    tr28 = sorted(st28.trials, key=lambda t: t.number)
    assert [t.user_attrs.get('knockout_probe') for t in tr28[:3]] \
        == [None, 'k_1e', 'k_7']                              # no k_1ie probe
    assert tr28[3].user_attrs.get('knockout_probe') is None
    assert np.isclose(df28['k_1e'][1], 1e-5*47.1)
    assert np.allclose(df28['k_1ie'][:3], 0.05)
    assert df28['k_1ie'].min() >= 0.1*0.05 - 1e-12             # 0.1x floor, never below
    assert df28['k_1ie'].max() <= 10.0*0.05 + 1e-12
    assert df28['k_1e'].min() < 0.1*47.1                       # rate band goes below 0.1x
    dist28 = tr28[5].distributions
    assert np.isclose(dist28['k_1ie'].low, 0.1*0.05) and dist28['k_1ie'].log
    assert np.isclose(dist28['k_1e'].low, 1e-5*47.1) and dist28['k_1e'].log
PASS('role-based bands: rate constants on the rate band, inhibition coefficients and K_* 0.1x-10x; rate_params threaded (None = legacy prefix rule); probes only for rate constants; _ib naming tag; driver/supervisor plumbing')

#%% 29. per-parameter multiplier bands (2026-09-06 pm): the rate band is
# 1e-3x-10x and k_10 (active-biomass decay capacity, r10) keeps 0.1x-10x
# via DEFAULT_PARAMETER_MULTIPLIER_BOUNDS; parameter_multiplier_bounds
# threaded through build_search_space / workbook_kinetic_bounds / the
# engine / the driver; precedence absolute override > per-parameter band
# > role band; every preset name tags the effective rate band.
assert ko.DEFAULT_RATE_MULTIPLIER_BOUNDS == (1e-3, 10.0)
assert ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS == {'k_10': (0.1, 10.0)}
assert 'DEFAULT_PARAMETER_MULTIPLIER_BOUNDS' in ko.__all__
kb29 = {'k_1e': 47.1, 'k_10': 0.06, 'k_10ie': 0.04, 'K_1e': 0.12, 'k_7': 1.203}
roles29 = {'k_1e': 'capacity', 'k_10': 'capacity', 'k_10ie': 'lethality',
           'K_1e': 'affinity', 'k_7': 'capacity'}
rp29 = ko.rate_constant_names(kb29, roles=roles29)
assert rp29 == ['k_1e', 'k_10', 'k_7']
space29, excl29 = ko.build_search_space(
    kb29, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-3, 10.0),
    rate_params=rp29, parameter_multiplier_bounds={'k_10': (0.1, 10.0)})
assert excl29 == []
assert space29['k_1e'] == dict(low=1e-3*47.1, high=10.0*47.1, log=True)
assert space29['k_7'] == dict(low=1e-3*1.203, high=10.0*1.203, log=True)
assert space29['k_10'] == dict(low=0.1*0.06, high=10.0*0.06, log=True)     # per-parameter band
assert space29['k_10ie'] == dict(low=0.1*0.04, high=10.0*0.04, log=True)   # untouched
assert space29['K_1e'] == dict(low=0.1*0.12, high=10.0*0.12, log=True)
# Per-parameter band applies regardless of the rate rule (a K_* too), and
# under the legacy prefix rule (rate_params=None).
space29b, _ = ko.build_search_space(
    kb29, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-3, 10.0),
    parameter_multiplier_bounds={'k_10': (0.5, 2.0), 'K_1e': (0.01, 100.0)})
assert space29b['k_10'] == dict(low=0.5*0.06, high=2.0*0.06, log=True)
assert space29b['K_1e'] == dict(low=0.01*0.12, high=100.0*0.12, log=True)
assert space29b['k_10ie']['low'] == 1e-3*0.04                              # prefix rule
# Precedence: absolute override > per-parameter band; whitelist/exclude
# still win over both; a per-parameter entry for an absent/excluded name
# is ignored; None / {} = no per-parameter band (the pre-change space).
space29c, _ = ko.build_search_space(
    kb29, rate_multiplier_bounds=(1e-3, 10.0), rate_params=rp29,
    parameter_multiplier_bounds={'k_10': (0.1, 10.0), 'k_1e': (0.2, 5.0),
                                 'k_99': (0.1, 10.0)},
    param_bounds_override={'k_10': (0.001, 0.002)},
    include_params=['k_10', 'k_1e', 'k_7'], exclude_params=('k_7',))
assert space29c['k_10'] == dict(low=0.001, high=0.002, log=True)
assert space29c['k_1e'] == dict(low=0.2*47.1, high=5.0*47.1, log=True)
assert set(space29c) == {'k_10', 'k_1e', *ko.FEEDING_VARIABLES}
space29d, _ = ko.build_search_space(
    kb29, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-3, 10.0),
    rate_params=rp29, parameter_multiplier_bounds=None)
assert space29d['k_10'] == dict(low=1e-3*0.06, high=10.0*0.06, log=True)
assert ko.build_search_space(
    kb29, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-3, 10.0),
    rate_params=rp29, parameter_multiplier_bounds={})[0] == space29d
# A nonpositive baseline is still excluded (a multiplier band needs one).
space29e, excl29e = ko.build_search_space(
    dict(kb29, k_10=0.0), rate_multiplier_bounds=(1e-3, 10.0), rate_params=rp29,
    parameter_multiplier_bounds={'k_10': (0.1, 10.0)})
assert excl29e == ['k_10'] and 'k_10' not in space29e
# Legacy single-band space (no rate band, no per-parameter band) EXACTLY as before.
assert ko.build_search_space(kb29)[0] == {
    'k_1e': dict(low=0.1*47.1, high=10.0*47.1, log=True),
    'k_10': dict(low=0.1*0.06, high=10.0*0.06, log=True),
    'k_10ie': dict(low=0.1*0.04, high=10.0*0.04, log=True),
    'K_1e': dict(low=0.1*0.12, high=10.0*0.12, log=True),
    'k_7': dict(low=0.1*1.203, high=10.0*1.203, log=True),
    'threshold_conc': dict(low=0.0, high=300.0, log=False),
    'target_delta': dict(low=5.0, high=500.0, log=False),
    'spike_delta': dict(low=0.5, high=595.0, log=False),
    'max_n_spikes': dict(low=0, high=50, log=False, int=True)}
# knockout probes: k_10's probe sits at ITS floor (0.1x, a knock-down).
base29 = ko.baseline_decision_point(space29, kb29,
                                    fbs17.current_specifications, 16)
probes29, at_floor29 = ko.knockout_probe_points(space29, base29, rate_params=rp29)
assert list(probes29) == ['k_1e', 'k_10', 'k_7'] and at_floor29 == []
assert probes29['k_10']['k_10'] == 0.1*0.06 and probes29['k_1e']['k_1e'] == 1e-3*47.1
# Engine kwarg plumbed (default None) and forwarded to build_search_space.
_sig29 = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'parameter_multiplier_bounds' in _sig29 \
    and _sig29['parameter_multiplier_bounds'].default is None
assert 'parameter_multiplier_bounds=parameter_multiplier_bounds' \
    in _inspect.getsource(ko.run_kinetic_optimization)
# workbook_kinetic_bounds: the per-parameter band on the real workbooks.
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    roles29_real = ko.kinetic_parameter_roles()
    assert roles29_real['k_10'] == 'capacity'
    for sc29 in ('A', 'B'):
        base29_wb = ko.workbook_kinetic_baselines(sc29)
        rp29_wb = ko.rate_constant_names(base29_wb, roles=roles29_real)
        assert 'k_10' in rp29_wb
        wbb29 = ko.workbook_kinetic_bounds(
            sc29, multiplier_bounds=(0.1, 10.0),
            rate_multiplier_bounds=(1e-3, 10.0), rate_params=rp29_wb,
            parameter_multiplier_bounds=ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS)
        assert list(wbb29) == list(base29_wb)
        for n29, b29 in base29_wb.items():
            if n29 == 'k_10':
                assert wbb29[n29] == (0.1*b29, 10.0*b29)
            else:
                assert wbb29[n29] == ((1e-3*b29, 10.0*b29) if n29 in rp29_wb
                                      else (0.1*b29, 10.0*b29)), n29
        # Without the table, k_10 is an ordinary rate constant (pre-change).
        wbb29_none = ko.workbook_kinetic_bounds(
            sc29, multiplier_bounds=(0.1, 10.0),
            rate_multiplier_bounds=(1e-3, 10.0), rate_params=rp29_wb)
        assert wbb29_none['k_10'] == (1e-3*base29_wb['k_10'], 10.0*base29_wb['k_10'])
        assert all(wbb29_none[n] == wbb29[n] for n in base29_wb if n != 'k_10')
    # Every preset carries a COPY of the default table.
    for stp29, st29 in (('ethanol_only', 'metabolic'),
                        ('ethanol_isobutanol', 'metabolic_protein')):
        p29 = ko.resolve_study_preset(stp29, st29)
        assert p29['parameter_multiplier_bounds'] == ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS
        assert p29['parameter_multiplier_bounds'] is not ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS
        assert p29['rate_multiplier_bounds'] == (1e-3, 10.0)
else:
    print('SKIP 29 (workbook part): parameter-distribution workbooks not found')
# Naming: the EFFECTIVE rate band is always tagged on preset names (engine
# rule: whenever given), before _ib; the legacy path never encodes it.
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0), burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_burden'
sup29 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert sup29['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
assert sup29['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True,
                                   rate_multiplier_bounds=(1e-5, 10.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb1e-05-10_ib0.1-10_xk10_s1x1-50_burden'
assert sup29['default_study_name']('A', 'IRR', 'B', burden=True) == 'kin_opt_A_kbB_irr_burden'
src29_sup = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'optimize_kinetics_BO_supervised.py')).read()
assert 'ko.DEFAULT_RATE_MULTIPLIER_BOUNDS' in _inspect.getsource(sup29['default_study_name'])
assert "1e-5 10" not in src29_sup                          # stale help text
# Driver: the preset's table is defaulted into engine_kwargs, forwarded to
# the workbook bounds, and the study name sees the EFFECTIVE rate band.
drv29 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert ("'rate_multiplier_bounds', 'rate_params',\n"
        "                    'parameter_multiplier_bounds', 'stage_1_max_x_bounds',\n"
        "                    'parameter_groups', 'group_multiplier_bounds',\n"
        "                    'spike_delta_bounds'):") in drv29
assert "parameter_multiplier_bounds=engine_kwargs.get('parameter_multiplier_bounds')" in drv29
assert "rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds']" in drv29
assert 'explicit_rate_bounds' not in drv29
PASS('per-parameter bands: rate band 1e-3x-10x, k_10 0.1x-10x via DEFAULT_PARAMETER_MULTIPLIER_BOUNDS; precedence override > per-parameter > role; probes at own floor; workbook/preset/engine/driver plumbing; _rb always tagged')

#%% 30. n_startup_trials: the TPE random start-up length is an engine kwarg
# (None = the default rule max(10, n_trials//4)), forwarded by the driver's
# run() and the supervisor's --n-startup-trials; not part of the study name.
_sig30 = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'n_startup_trials' in _sig30 and _sig30['n_startup_trials'].default is None
assert 'n_startup_trials' not in _inspect.signature(ko.default_study_name).parameters
class _FakeTE30:
    k_1e = 47.1; k_7 = 1.203; K_1e = 0.12
    def getGlobalParameterIds(self):
        return ['k_1e', 'k_7', 'K_1e']
def _solve_TEA30(stream_IDs=None):
    return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
handles30 = dict(handles17, r_te=_FakeTE30(),
                 model_specification=lambda **kw: None, solve_TEA=_solve_TEA30,
                 latest_TEA_solution={'IRR': np.nan,
                                      'MPSPs': {'ethanol': np.nan,
                                                'isobutanol': np.nan}})
outdir30 = tempfile.mkdtemp()
common30 = dict(objective='IRR', scenario_label='X', seed=1,
                study_name='offline_startup', results_dir=outdir30,
                handles=handles30, print_status_every=10, burden_model=None,
                enqueue_knockouts=False)
st30, csv30, _ = ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=12, n_startup_trials=3,
                                             **common30)
assert st30.sampler._n_startup_trials == 3
assert len(ko.load_trajectory(csv30)) == 12
# Resume (nothing left to run) with None: the default rule, floor 10.
st30b, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=12, **common30)
assert st30b.sampler._n_startup_trials == 10 == max(10, 12//4)
assert len(ko.load_trajectory(csv30)) == 12                  # no new trials
st30c, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=12, n_startup_trials=0,
                                          **common30)
assert st30c.sampler._n_startup_trials == 0
try:
    ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=12, n_startup_trials=-1, **common30)
except ValueError as e30:
    assert 'n_startup_trials' in str(e30)
else:
    raise AssertionError('negative n_startup_trials did not raise')
# Driver: explicit run() kwarg (default None), forwarded to the engine.
drv30 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'n_startup_trials=None,' in drv30
assert 'n_startup_trials=n_startup_trials' in drv30
# Supervisor: --n-startup-trials, forwarded through supervise() and
# child_code() (the kwarg is OMITTED when None, so the driver's default
# rule applies); no effect on the study name.
sup30 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
_s30 = _inspect.signature(sup30['supervise']).parameters
_c30 = _inspect.signature(sup30['child_code']).parameters
assert _s30['n_startup_trials'].default is None and _c30['n_startup_trials'].default is None
code30 = sup30['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic', n_startup_trials=25)
assert 'n_startup_trials=25' in code30
code30b = sup30['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic')
assert 'n_startup_trials' not in code30b
assert 'n_startup_trials=n_startup_trials' in _inspect.getsource(sup30['supervise'])
assert 'n_startup_trials' not in _inspect.signature(sup30['default_study_name']).parameters
src30 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--n-startup-trials'" in src30 and 'n_startup_trials=args.n_startup_trials' in src30
# Default rule is 25% of the budget floored at 10 (set 2026-09-10; was
# n_trials//10). The n_trials=12 resume above floors to 10 and can't tell
# the two apart, so pin the rule at the source and at a distinguishing n.
_eng_src30 = _inspect.getsource(ko.run_kinetic_optimization)
assert 'n_startup = max(10, n_trials//4)' in _eng_src30      # the active rule
assert max(10, 2000//4) == 500 and max(10, 2000//10) == 200   # 25% vs old 10%
PASS('n_startup_trials: engine kwarg (None = max(10, n_trials//4), 25% floored at 10), validated, resume-safe; driver run() kwarg; supervisor --n-startup-trials; study name untouched')

#%% 31. empty-attempt abort rule (2026-09-06): an attempt that logged no
# new row but whose child had already STARTED a simulation (in-flight
# sidecar present => the reload worked, the first draw hung) resumes; abort
# only when the child never reached a simulation (kill-loop / broken load)
# or after `max_empty_attempts` consecutive empty attempts (safety net).
# Legacy defaults (no inflight flag, max 1) reproduce the old rule exactly.
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True) == 'abort'           # legacy
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True, inflight_lost=True) \
    == 'abort'                                                                 # max 1: still abort
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True, inflight_lost=True,
                          empty_streak=0, max_empty_attempts=5) == 'resume'
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True, inflight_lost=True,
                          empty_streak=3, max_empty_attempts=5) == 'resume'
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True, inflight_lost=True,
                          empty_streak=4, max_empty_attempts=5) == 'abort'     # 5th consecutive
assert ko.attempt_outcome(1, 5, 5, killed_for_stall=True, inflight_lost=False,
                          empty_streak=0, max_empty_attempts=5) == 'abort'     # never simulated
assert ko.attempt_outcome(3221225477, 5, 5, inflight_lost=True,
                          max_empty_attempts=5) == 'resume'                    # crash on first draw
assert ko.attempt_outcome(3221225477, 5, 5, inflight_lost=False,
                          max_empty_attempts=5) == 'abort'
assert ko.attempt_outcome(1, 5, 8, killed_for_stall=True, empty_streak=4,
                          max_empty_attempts=5) == 'resume'                    # progress: streak irrelevant
assert ko.attempt_outcome(0, 5, 5, empty_streak=4, max_empty_attempts=5) == 'complete'
try:
    ko.attempt_outcome(1, 5, 5, max_empty_attempts=0)
except ValueError as e31:
    assert 'max_empty_attempts' in str(e31)
else:
    raise AssertionError('max_empty_attempts=0 did not raise')
# Supervisor: tracks the streak across attempts, checks the sidecar BEFORE
# recovering it (the LOST row still never counts as progress), exposes
# max_empty_attempts (default 5) on supervise() and --max-empty-attempts.
sup31 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert _inspect.signature(sup31['supervise']).parameters['max_empty_attempts'].default == 5
src31 = _inspect.getsource(sup31['supervise'])
assert 'inflight_lost=' in src31 and 'empty_streak' in src31
assert 'max_empty_attempts=max_empty_attempts' in src31
assert src31.index('attempt_outcome(') < src31.rindex('recover_inflight(')     # still decided before recovery
assert src31.index('os.path.isfile(inflight_path)') < src31.index('attempt_outcome(')
file31 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'optimize_kinetics_BO_supervised.py')).read()
assert "'--max-empty-attempts'" in file31 and 'max_empty_attempts=args.max_empty_attempts' in file31
PASS('empty-attempt abort rule: first-draw hang resumes (sidecar present), never-simulated aborts, streak cap; supervisor streak + flag')

#%% 32. feasibility-aware sampling helpers (2026-09-06, spec
# docs/superpowers/specs/2026-09-06-feasible-tpe-sampler-design.md):
# search_space_distributions reproduces the distributions a real suggest_*
# trial records; draw_uniform_feasible is uniform in optuna's internal repr
# (log-uniform for log floats, integer-uniform for ints), in bounds, typed,
# and gives up after exactly max_draws; feasible_candidate_mask converts a
# Parzen batch to external values per candidate.
space32 = {'a': dict(low=0.01, high=100.0, log=True),
           'b': dict(low=0.0, high=5.0, log=False),
           'n': dict(low=0, high=10, log=False, int=True)}
if _optuna is None:
    print('SKIP 32: optuna not installed')
else:
    dists32 = ko.search_space_distributions(space32)
    assert list(dists32) == ['a', 'b', 'n']
    def _obj32(trial):
        trial.suggest_float('a', 0.01, 100.0, log=True)
        trial.suggest_float('b', 0.0, 5.0)
        trial.suggest_int('n', 0, 10)
        return 0.0
    st32 = _optuna.create_study(sampler=_optuna.samplers.RandomSampler(seed=0))
    st32.optimize(_obj32, n_trials=1)
    assert st32.trials[0].distributions == dists32, (st32.trials[0].distributions, dists32)
    # uniform-feasible draw: bounds, types, log-uniformity, draw count
    rng32 = np.random.RandomState(0)
    n_calls32 = []
    def _always32(values):
        n_calls32.append(values)
        return True
    decades32 = np.zeros(4, dtype=int)     # [0.01,0.1) [0.1,1) [1,10) [10,100]
    for _ in range(4000):
        values32, n_draws32, ok32 = ko.draw_uniform_feasible(rng32, dists32, _always32)
        assert ok32 and n_draws32 == 1
        assert set(values32) == {'a', 'b', 'n'}
        assert 0.01 <= values32['a'] <= 100.0 and 0.0 <= values32['b'] <= 5.0
        assert isinstance(values32['a'], float) and isinstance(values32['b'], float)
        assert isinstance(values32['n'], int) and 0 <= values32['n'] <= 10
        decades32[min(3, int(np.floor(np.log10(values32['a'])) + 2))] += 1
    assert len(n_calls32) == 4000
    frac32 = decades32/4000.0
    assert np.all(np.abs(frac32 - 0.25) < 0.03), frac32           # log-uniform in a
    seen_n32 = {v['n'] for v in n_calls32}
    assert seen_n32 == set(range(11)), seen_n32                   # every int reached
    # always-false predicate: exactly max_draws draws, feasible=False
    n_false32 = []
    def _never32(values):
        n_false32.append(values)
        return False
    values32f, n_draws32f, ok32f = ko.draw_uniform_feasible(
        np.random.RandomState(1), dists32, _never32, max_draws=7)
    assert (not ok32f) and n_draws32f == 7 and len(n_false32) == 7
    assert values32f == n_false32[-1]                              # the LAST draw is returned
    try:
        ko.draw_uniform_feasible(np.random.RandomState(1), dists32, _never32, max_draws=0)
    except ValueError as e32:
        assert 'max_draws' in str(e32)
    else:
        raise AssertionError('max_draws=0 did not raise')
    # candidate mask over a Parzen-style batch (internal repr: floats)
    batch32 = {'a': np.array([1.0, 50.0, 2.0]),
               'b': np.array([1.0, 1.0, 4.0]),
               'n': np.array([0.0, 3.0, 10.0])}
    types32 = []
    def _cap32(values):
        types32.append(type(values['n']))
        return values['a'] + values['b'] < 5.0
    mask32 = ko.feasible_candidate_mask(batch32, dists32, _cap32)
    assert mask32.dtype == bool and mask32.tolist() == [True, False, False]
    assert types32 == [int, int, int]                              # ints converted per candidate
    for name in ('search_space_distributions', 'draw_uniform_feasible',
                 'feasible_candidate_mask'):
        assert name in ko.__all__, name
    PASS('feasible-sampling helpers: distributions match suggest_*; uniform-feasible draw in bounds, typed, log-uniform, max_draws; candidate mask')

#%% 33. FeasibleTPESampler end to end on a toy capped problem: zero sampled
# infeasible trials (a plain TPESampler with the same seed samples some),
# enqueued trials bypass the sampler, the objective still improves,
# counters consistent, resume with a second sampler instance, seed
# determinism.
if _optuna is None:
    print('SKIP 33: optuna not installed')
else:
    _optuna.logging.set_verbosity(_optuna.logging.WARNING)
    space33 = {'a': dict(low=0.01, high=100.0, log=True),
               'b': dict(low=0.01, high=100.0, log=True),
               'n': dict(low=0, high=10, log=False, int=True)}
    def feas33(values):                       # the cap: a + b < 5 (~40 % of log-uniform draws)
        return values['a'] + values['b'] < 5.0
    def cons33(frozen_trial):
        return (frozen_trial.user_attrs.get('violation', 0.0),)
    def _toy33(trial):
        a = trial.suggest_float('a', 0.01, 100.0, log=True)
        b = trial.suggest_float('b', 0.01, 100.0, log=True)
        n = trial.suggest_int('n', 0, 10)
        violation = a + b - 5.0
        trial.set_user_attr('violation', violation)
        if violation >= 0.0:                  # the engine's INFEASIBLE guard
            raise _optuna.TrialPruned()
        return -((a - 3.0)**2 + (b - 1.0)**2) - n     # optimum 0 at (3, 1, 0)
    TS33 = _optuna.trial.TrialState
    def _sampled_infeasible33(study):
        return [t.number for t in study.trials
                if t.state == TS33.PRUNED and 'enqueued' not in t.user_attrs]
    def _fresh33(sampler):
        study = _optuna.create_study(direction='maximize', sampler=sampler)
        study.enqueue_trial({'a': 50.0, 'b': 50.0, 'n': 3},
                            user_attrs={'enqueued': True})
        return study
    tpe_kw33 = dict(multivariate=True, n_startup_trials=10, constraints_func=cons33)
    samp33 = ko.feasible_tpe_sampler(space33, feas33, seed=11, **tpe_kw33)
    assert type(samp33).__name__ == 'FeasibleTPESampler'
    assert isinstance(samp33, _optuna.samplers.TPESampler)
    st33 = _fresh33(samp33)
    st33.optimize(_toy33, n_trials=60)
    assert len(st33.trials) == 60
    # the enqueued infeasible trial 0 reached the objective untouched (bypass)
    assert st33.trials[0].state == TS33.PRUNED and st33.trials[0].params == {'a': 50.0, 'b': 50.0, 'n': 3}
    assert _sampled_infeasible33(st33) == [], _sampled_infeasible33(st33)
    assert all(t.state == TS33.COMPLETE for t in st33.trials[1:])
    assert all(feas33(t.params) for t in st33.trials[1:])
    assert st33.best_value > -10.0, st33.best_value          # random feasible points score ~ -30
    assert samp33.n_rejected > 0 and samp33.n_uniform_fallbacks == 0 and samp33.n_unfiltered == 0
    # the same seed under a plain TPESampler proposes infeasible points
    # (expected ~5 of the 9 sampled start-up draws alone)
    plain33 = _optuna.samplers.TPESampler(seed=11, **tpe_kw33)
    stp33 = _fresh33(plain33)
    stp33.optimize(_toy33, n_trials=60)
    assert len(_sampled_infeasible33(stp33)) >= 1, _sampled_infeasible33(stp33)
    # resume: a second sampler instance on the study (which holds a pruned
    # infeasible trial with constraints attrs) keeps sampling feasibly
    assert list(st33.trials[0].system_attrs['constraints'])[0] > 0.0
    samp33b = ko.feasible_tpe_sampler(space33, feas33, seed=12, **tpe_kw33)
    st33.sampler = samp33b
    st33.optimize(_toy33, n_trials=10)
    assert len(st33.trials) == 70 and _sampled_infeasible33(st33) == []
    assert samp33b.n_unfiltered == 0 and all(feas33(t.params) for t in st33.trials[60:])
    # determinism: same seed -> same first sampled trial (start-up path) and
    # same first TPE-phase trial
    def _run33(seed, n):
        s = ko.feasible_tpe_sampler(space33, feas33, seed=seed, **tpe_kw33)
        st = _optuna.create_study(direction='maximize', sampler=s)
        st.optimize(_toy33, n_trials=n)
        return [t.params for t in st.trials]
    assert _run33(5, 11) == _run33(5, 11)
    assert _run33(5, 1) != _run33(6, 1)
    PASS('FeasibleTPESampler: zero sampled infeasible trials on a capped toy (plain TPE samples some), enqueued bypass, improvement, counters, resume, seed determinism')

#%% 34. fallbacks and factory guards: the wrapped "below" estimator falls
# back to one uniform-feasible draw when every Parzen batch is infeasible,
# returns the raw last batch (counted) when the predicate never accepts,
# concatenates partial batches to exactly `size`; the real sampler under an
# always-false predicate never raises and counts exactly; group /
# constant_liar / bad caps / non-callable predicate raise in the factory.
if _optuna is None:
    print('SKIP 34: optuna not installed')
else:
    cls34 = ko._feasible_tpe_sampler_class()
    assert cls34 is ko._feasible_tpe_sampler_class()               # memoized
    dists34 = ko.search_space_distributions(space33)
    class _FakeMPE34:
        def __init__(self, batches):
            self.batches, self.calls = list(batches), 0
        def sample(self, rng, size):
            self.calls += 1
            return self.batches[min(self.calls, len(self.batches)) - 1]
        def log_pdf(self, samples):
            return np.zeros(len(next(iter(samples.values()))))
    def _counters34():
        return SimpleNamespace(max_uniform_draws=200, max_parzen_batches=3,
                               n_rejected=0, n_uniform_fallbacks=0, n_unfiltered=0)
    bad34 = {'a': np.full(4, 50.0), 'b': np.full(4, 50.0), 'n': np.zeros(4)}
    mixed34 = {'a': np.array([1.0, 50.0, 2.0, 60.0]),
               'b': np.array([1.0, 1.0, 1.0, 1.0]),
               'n': np.array([0.0, 1.0, 2.0, 3.0])}
    # (a) every batch infeasible, region reachable -> uniform fallback, size-1 batch
    mpe_a, s_a = _FakeMPE34([bad34]), _counters34()
    out_a = ko._FeasibleParzenEstimator(mpe_a, dists34, feas33, s_a).sample(np.random.RandomState(0), 4)
    assert mpe_a.calls == 3 and s_a.n_uniform_fallbacks == 1 and s_a.n_unfiltered == 0
    assert list(out_a) == ['a', 'b', 'n'] and all(len(v) == 1 for v in out_a.values())
    assert feas33({k: dists34[k].to_external_repr(float(v[0])) for k, v in out_a.items()})
    assert float(out_a['n'][0]).is_integer()
    assert s_a.n_rejected >= 3*4                                     # 12 candidates + the uniform misses
    # (b) predicate never accepts -> raw last batch, unfiltered counted, no raise
    mpe_b, s_b = _FakeMPE34([bad34]), _counters34()
    out_b = ko._FeasibleParzenEstimator(mpe_b, dists34, lambda v: False, s_b).sample(np.random.RandomState(0), 4)
    assert out_b is bad34 and s_b.n_uniform_fallbacks == 1 and s_b.n_unfiltered == 1
    assert s_b.n_rejected == 3*4 + 200
    # (c) half-feasible batches concatenate to exactly `size`, no fallback
    mpe_c, s_c = _FakeMPE34([mixed34]), _counters34()
    w_c = ko._FeasibleParzenEstimator(mpe_c, dists34, feas33, s_c)
    out_c = w_c.sample(np.random.RandomState(0), 4)
    assert mpe_c.calls == 2 and s_c.n_uniform_fallbacks == 0 and s_c.n_unfiltered == 0
    assert out_c['a'].tolist() == [1.0, 2.0, 1.0, 2.0] and out_c['n'].tolist() == [0.0, 2.0, 0.0, 2.0]
    assert s_c.n_rejected == 4
    assert w_c.log_pdf(mixed34).tolist() == [0.0]*4                 # delegates
    # (d) the real sampler under an always-false predicate: 2 raw start-up
    # draws (max_uniform_draws=50 each), then one TPE trial whose 20 Parzen
    # batches of 24 (optuna's n_ei_candidates default) are all rejected,
    # the uniform fallback misses 50 times, and the raw batch is used.
    samp34 = ko.feasible_tpe_sampler(space33, lambda v: False, max_uniform_draws=50,
                                     multivariate=True, seed=3, n_startup_trials=2)
    st34 = _optuna.create_study(direction='maximize', sampler=samp34)
    st34.optimize(lambda t: (t.suggest_float('a', 0.01, 100.0, log=True)
                             + t.suggest_float('b', 0.01, 100.0, log=True)
                             + t.suggest_int('n', 0, 10)), n_trials=3)
    assert len(st34.trials) == 3 and all(t.state == _optuna.trial.TrialState.COMPLETE for t in st34.trials)
    assert samp34.n_unfiltered == 3 and samp34.n_uniform_fallbacks == 1
    assert samp34.n_rejected == 2*50 + 20*24 + 50, samp34.n_rejected
    # factory guards
    for bad_kw in (dict(group=True), dict(constant_liar=True),
                   dict(max_uniform_draws=0), dict(max_parzen_batches=0)):
        try:
            ko.feasible_tpe_sampler(space33, feas33, **bad_kw)
        except ValueError:
            pass
        else:
            raise AssertionError(f'{bad_kw} did not raise')
    try:
        ko.feasible_tpe_sampler(space33, None)
    except TypeError:
        pass
    else:
        raise AssertionError('non-callable predicate did not raise')
    # the relative search space is the FULL engine space (not optuna's
    # intersection of stored trials), single() entries dropped
    st34b = _optuna.create_study(sampler=ko.feasible_tpe_sampler(
        {**space33, 'c': dict(low=1.0, high=1.0, log=False)}, feas33, seed=0))
    rel34 = st34b.sampler.infer_relative_search_space(st34b, None)
    assert rel34 == dists34, rel34
    PASS('feasible sampling fallbacks: uniform fallback, raw-batch last resort, partial-batch concatenation, exact counters on a real sampler; factory guards; full relative space')

#%% 35. engine wiring: feasible_sampling=True (default) with the burden on
# installs a FeasibleTPESampler and the study samples NO infeasible point
# (every sampled trial is COMPLETE with Phi_M < F_flex, n_unfiltered = 0);
# feasible_sampling=False, or burden off, installs the plain TPESampler;
# the printed sampler / summary lines; the kwarg never reaches the study
# name.
if _optuna is None:
    print('SKIP 35: optuna not installed')
else:
    _sig35 = _inspect.signature(ko.run_kinetic_optimization).parameters
    assert 'feasible_sampling' in _sig35 and _sig35['feasible_sampling'].default is True
    assert 'feasible_sampling' not in _inspect.signature(ko.default_study_name).parameters
    class _FakeTE35:
        # the model's reference capacities (scenario A: Ehrlich off), so
        # BurdenModel.from_reference is inert at the baseline
        k_1h = 0.584; k_1l = 1.43; k_1e = 47.1; k_2 = 0.501; k_3 = 5.81
        k_4 = 4.8; k_5 = 0.0104; k_5e = 0.775; k_6 = 2.82
        k_7 = 1.203; k_8 = 0.589
        k_13 = 0.0; k_14 = 0.0; k_15 = 0.0; k_16 = 0.0
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1h', 'k_1l', 'k_1e', 'k_2', 'k_3', 'k_4', 'k_5',
                    'k_5e', 'k_6', 'k_7', 'k_8', 'k_13', 'k_14', 'k_15',
                    'k_16', 'K_1e', 'not_kinetic']
    def _solve_TEA35(stream_IDs=None):
        return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles35 = dict(handles17, r_te=_FakeTE35(),
                     model_specification=lambda **kw: None, solve_TEA=_solve_TEA35,
                     latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    outdir35 = tempfile.mkdtemp()
    override35 = {'k_13': (0.0, 60.0), 'k_14': (0.0, 50.0),
                  'k_15': (0.0, 50.0), 'k_16': (0.0, 30.0)}
    common35 = dict(objective='IRR', scenario_label='X', seed=1,
                    study_name='offline_feasible', results_dir=outdir35,
                    handles=handles35, print_status_every=10,
                    param_bounds_override=override35,
                    rate_multiplier_bounds=(1e-3, 10.0),
                    enqueue_knockouts=False)
    buf35 = _io.StringIO()
    with _contextlib.redirect_stdout(buf35):
        st35, csv35, _ = ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=10, n_startup_trials=4,
                                                     **common35)   # burden 'auto', feasible_sampling default
    out35 = buf35.getvalue()
    assert type(st35.sampler).__name__ == 'FeasibleTPESampler', type(st35.sampler)
    assert st35.sampler._n_startup_trials == 4
    assert 'Sampler: feasibility-aware TPE' in out35, out35
    assert 'feasibility-aware: joint uniform-feasible draws' in out35, out35
    assert 'Feasible sampling: rejected ' in out35 and ' 0 unfiltered draws.' in out35, out35
    df35 = ko.load_trajectory(csv35)
    assert len(df35) == 10 and df35['state'].tolist() == ['COMPLETE']*10, df35['state'].tolist()
    assert (df35['Phi_M'] < df35['F_flex']).all()
    assert df35['trial_number'].tolist() == list(range(10))
    # in this space a 10x k_1e alone (0.44 g/gDCW on the 0.044 r1 pool)
    # breaks the cap, so the 9 sampled trials rejected at least one draw
    assert st35.sampler.n_rejected > 0 and st35.sampler.n_unfiltered == 0
    # the objective's guard never fired: no INFEASIBLE row, no pruned trial
    assert all(t.state == _optuna.trial.TrialState.COMPLETE for t in st35.trials)
    # resume with feasible_sampling=False -> the plain TPESampler, same study
    buf35b = _io.StringIO()
    with _contextlib.redirect_stdout(buf35b):
        st35b, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, n_trials=12, n_startup_trials=4,
                                                  feasible_sampling=False, **common35)
    assert type(st35b.sampler) is _optuna.samplers.TPESampler
    assert 'Sampler: plain TPESampler (feasible_sampling=False)' in buf35b.getvalue()
    assert 'Feasible sampling: rejected' not in buf35b.getvalue()
    assert len(ko.load_trajectory(csv35)) == 12
    # burden off + feasible_sampling=True -> plain sampler (no predicate)
    buf35c = _io.StringIO()
    with _contextlib.redirect_stdout(buf35c):
        st35c, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True,
            n_trials=3, n_startup_trials=2, feasible_sampling=True,
            burden_model=None, **{**common35, 'study_name': 'offline_feasible_noburden'})
    assert type(st35c.sampler) is _optuna.samplers.TPESampler
    assert 'Sampler: plain TPESampler (burden off' in buf35c.getvalue()
    PASS('engine: feasible_sampling=True + burden on installs FeasibleTPESampler (all sampled trials feasible, n_unfiltered 0, summary line); False or burden off = plain TPESampler; not in the study name')

#%% 36. driver run(feasible_sampling=True) forwarded to the engine; supervisor
# --no-feasible-sampling (store_false, default on) forwarded through
# supervise() and child_code() on every attempt, shown in the settings
# event; never part of the study name.
drv36 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'feasible_sampling=True,' in drv36                          # run() kwarg, default on
assert 'feasible_sampling=feasible_sampling' in drv36              # forwarded to the engine
sup36 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert _inspect.signature(sup36['supervise']).parameters['feasible_sampling'].default is True
assert _inspect.signature(sup36['child_code']).parameters['feasible_sampling'].default is True
code36 = sup36['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic')
assert 'feasible_sampling=True,' in code36
code36b = sup36['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', feasible_sampling=False)
assert 'feasible_sampling=False,' in code36b
src36_sup = _inspect.getsource(sup36['supervise'])
assert 'feasible_sampling=feasible_sampling' in src36_sup
assert 'feasible_sampling={feasible_sampling!r}' in src36_sup     # settings event line
assert 'feasible_sampling' not in _inspect.signature(sup36['default_study_name']).parameters
src36 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--no-feasible-sampling'" in src36 and "dest='feasible_sampling'" in src36
assert 'feasible_sampling=args.feasible_sampling' in src36
PASS('feasible_sampling: driver run() kwarg (default on) forwarded; supervisor --no-feasible-sampling through supervise()/child_code(), settings event; study name untouched')

#%% 37. k_10 (active-biomass decay capacity) excluded from every study
# preset by default (2026-09-06 pm): DEFAULT_EXCLUDED_PARAMETERS, the
# preset's exclude_params, the _x{names} study-name tag (an exclusion drops
# a CSV column, so the header guard blocks a resume of the old set and the
# tag keeps the default name off those studies), no knockout probe for an
# excluded rate, driver setdefault + naming, supervisor --exclude-params
# through supervise()/child_code()/default_study_name.
assert ko.DEFAULT_EXCLUDED_PARAMETERS == ('k_10',)
assert 'DEFAULT_EXCLUDED_PARAMETERS' in ko.__all__ and 'excluded_parameters_tag' in ko.__all__
assert ko.excluded_parameters_tag(('k_10',)) == '_xk10'
assert ko.excluded_parameters_tag(['k_10', 'k_7']) == '_xk10+k7'
assert ko.excluded_parameters_tag(()) == '' and ko.excluded_parameters_tag(None) == ''
# Engine naming: tagged after _ib, before _burden; None / () untouched.
_sig37 = _inspect.signature(ko.default_study_name).parameters
assert 'exclude_params' in _sig37 and _sig37['exclude_params'].default is None
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0),
                             exclude_params=('k_10',), burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             exclude_params=('k_10',)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_xk10'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0),
                             exclude_params=(), burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_burden'  # re-included: the old name
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'                             # None: untouched
# Search space: excluded name absent, no probe, everything else identical.
kb37 = {'k_1e': 47.1, 'k_10': 0.06, 'k_10ie': 0.04, 'K_1e': 0.12, 'k_7': 1.203}
roles37 = {'k_1e': 'capacity', 'k_10': 'capacity', 'k_10ie': 'lethality',
           'K_1e': 'affinity', 'k_7': 'capacity'}
rp37 = ko.rate_constant_names(kb37, roles=roles37)
space37, excl37 = ko.build_search_space(
    kb37, multiplier_bounds=(0.1, 10.0), rate_multiplier_bounds=(1e-3, 10.0),
    rate_params=rp37, parameter_multiplier_bounds=dict(ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS),
    exclude_params=ko.DEFAULT_EXCLUDED_PARAMETERS)
assert excl37 == ['k_10'] and 'k_10' not in space37
assert set(space37) == {'k_1e', 'k_10ie', 'K_1e', 'k_7', *ko.FEEDING_VARIABLES}
assert space37['k_1e'] == dict(low=1e-3*47.1, high=10.0*47.1, log=True)
assert space37['k_10ie'] == dict(low=0.1*0.04, high=10.0*0.04, log=True)
base37 = ko.baseline_decision_point(space37, kb37, fbs17.current_specifications, 16)
probes37, at_floor37 = ko.knockout_probe_points(space37, base37, rate_params=rp37)
assert list(probes37) == ['k_1e', 'k_7'] and at_floor37 == []
assert ko.trajectory_columns(space37)[2:2 + len(space37)] == list(space37)  # no k_10 column
# Every preset returns a COPY of the default set (and its band table
# stays, for a re-included k_10).
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    for stp37, st37 in (('ethanol_only', 'metabolic'),
                        ('ethanol_isobutanol', 'metabolic_protein')):
        p37 = ko.resolve_study_preset(stp37, st37)
        assert p37['exclude_params'] == ko.DEFAULT_EXCLUDED_PARAMETERS
        assert isinstance(p37['exclude_params'], tuple)
        assert 'k_10' in p37['include_params'] and 'k_10' in p37['rate_params']
        assert p37['parameter_multiplier_bounds'] == {'k_10': (0.1, 10.0)}
# Driver: exclude_params is a preset key defaulted into engine_kwargs and
# the EFFECTIVE set reaches the study name; the preset print reports the
# effective (include minus exclude) set.
drv37 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert "for key in ('include_params', 'exclude_params', 'multiplier_bounds'," in drv37
assert "exclude_params=engine_kwargs['exclude_params']," in drv37
assert "excluded = tuple(engine_kwargs['exclude_params'] or ())" in drv37
assert 'if n not in excluded and n not in grouped]' in drv37
# Supervisor: --exclude-params (nargs='*', default None = the preset's;
# bare = () re-includes k_10), threaded through supervise() -> child_code()
# (kwarg omitted when None, forwarded as a tuple otherwise) and the naming
# mirror (None -> ko.DEFAULT_EXCLUDED_PARAMETERS), shown in the settings
# event.
sup37 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert _inspect.signature(sup37['supervise']).parameters['exclude_params'].default is None
assert _inspect.signature(sup37['child_code']).parameters['exclude_params'].default is None
assert _inspect.signature(sup37['default_study_name']).parameters['exclude_params'].default is None
code37 = sup37['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic')
assert 'exclude_params' not in code37                                  # None: preset default
code37b = sup37['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', exclude_params=())
assert 'exclude_params=(),' in code37b                                 # re-include k_10
code37c = sup37['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', exclude_params=['k_10', 'k_7'])
assert "exclude_params=('k_10', 'k_7')," in code37c
assert sup37['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic', burden=True,
                                   exclude_params=()) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_ib0.1-10_s1x1-50_burden'   # the production study's name + the _s1x tag (resume it with a bare --stage-1-max-x-bounds)
assert sup37['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic', burden=True,
                                   exclude_params=('k_10', 'k_7')) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_ib0.1-10_xk10+k7_s1x1-50_burden'
assert sup37['default_study_name']('A', 'IRR', 'B', exclude_params=('k_10',)) \
    == 'kin_opt_A_kbB_irr'                                             # legacy path: never tagged
src37_sup = _inspect.getsource(sup37['supervise'])
assert 'exclude_params=exclude_params' in src37_sup
assert 'exclude_params={exclude_params!r}' in src37_sup                # settings event line
assert 'ko.study_type_name_defaults(study_type)' in _inspect.getsource(sup37['default_study_name'])
src37 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--exclude-params', nargs='*', default=None" in src37
assert 'else tuple(args.exclude_params)' in src37
PASS('k_10 excluded by default: DEFAULT_EXCLUDED_PARAMETERS, preset exclude_params, _xk10 tag (None/() untouched), no column/probe; driver setdefault + naming; supervisor --exclude-params plumbing')

#%% 38. stage_1_max_x operating variable: search space, baseline point, no probe, name tag, preset key
assert ko.OPERATING_VARIABLES == ('stage_1_max_x',)
assert ko.DEFAULT_STAGE_1_MAX_X_BOUNDS == (1.0, 50.0)
assert {'OPERATING_VARIABLES', 'DEFAULT_STAGE_1_MAX_X_BOUNDS'} <= set(ko.__all__)
kb38 = {'k_1e': 47.1, 'K_1e': 0.12}
# Default: ABSENT (every study started before it keeps its columns).
space38_off, _ = ko.build_search_space(kb38)
assert 'stage_1_max_x' not in space38_off
assert set(space38_off) == {'k_1e', 'K_1e', *ko.FEEDING_VARIABLES}
# Bounds given: a log-scale float after max_n_spikes; include_params never
# filters it (an operating variable, like the feeding variables).
space38, _ = ko.build_search_space(kb38, include_params=['k_1e'],
                                   stage_1_max_x_bounds=(1.0, 50.0))
assert space38['stage_1_max_x'] == dict(low=1.0, high=50.0, log=True)
assert list(space38) == ['k_1e', *ko.FEEDING_VARIABLES, 'stage_1_max_x']
assert ko.trajectory_columns(space38)[2:2 + len(space38)] == list(space38)
assert ko.trajectory_columns(space38)[2 + len(space38)] == 'objective'
for bad38 in ((0.0, 50.0), (-1.0, 50.0), (50.0, 1.0), (5.0, 5.0)):
    try:
        ko.build_search_space(kb38, stage_1_max_x_bounds=bad38)
    except ValueError as e38:
        assert 'stage_1_max_x_bounds' in str(e38)
    else:
        raise AssertionError(f'stage_1_max_x_bounds={bad38} did not raise')
# Baseline point: the live V406 value when given (clipped into the band);
# absent when not given or when the variable is not in the space.
spec38 = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
base38 = ko.baseline_decision_point(space38, kb38, spec38, 16,
                                    baseline_stage_1_max_x=5.0)
assert base38['stage_1_max_x'] == 5.0
assert ko.baseline_decision_point(space38, kb38, spec38, 16,
                                  baseline_stage_1_max_x=60.0)['stage_1_max_x'] == 50.0
assert 'stage_1_max_x' not in ko.baseline_decision_point(space38, kb38, spec38, 16)
assert 'stage_1_max_x' not in ko.baseline_decision_point(
    space38_off, kb38, spec38, 16, baseline_stage_1_max_x=5.0)
# No knockout probe (not a rate constant) under either rate rule, although
# it is a log-scale variable above its floor.
probes38, at_floor38 = ko.knockout_probe_points(space38, base38)
assert list(probes38) == ['k_1e'] and at_floor38 == []
assert probes38['k_1e']['stage_1_max_x'] == 5.0
probes38r, _ = ko.knockout_probe_points(space38, base38, rate_params=['k_1e'])
assert list(probes38r) == ['k_1e']
# Study-name tag after the exclusion tag, before _burden; None untouched.
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0),
                             exclude_params=('k_10',), burden=True,
                             stage_1_max_x_bounds=(1.0, 50.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             stage_1_max_x_bounds=(0.5, 20.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_s1x0.5-20'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_protein',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0),
                             exclude_params=('k_10',), burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_burden'  # None: untouched
# Every preset carries the band (an immutable tuple equal to the default).
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    for stp38, st38 in (('ethanol_only', 'metabolic'),
                        ('ethanol_isobutanol', 'metabolic_protein')):
        p38 = ko.resolve_study_preset(stp38, st38)
        assert p38['stage_1_max_x_bounds'] == ko.DEFAULT_STAGE_1_MAX_X_BOUNDS
        assert isinstance(p38['stage_1_max_x_bounds'], tuple)
PASS('stage_1_max_x: OPERATING_VARIABLES/DEFAULT_STAGE_1_MAX_X_BOUNDS, opt-in log-scale space entry, ValueError on bad bounds, baseline point (clipped), no probe, _s1x tag, preset key')

#%% 39. engine applies stage_1_max_x via the V406 property per trial and restores the baseline
# Fake handles (no biorefinery, no simulation), the check-17 pattern: the
# scripted model_specification records V406.stage_1_max_x at call time.
if _optuna is None:
    print('SKIP 39: optuna not installed')
else:
    outdir39 = tempfile.mkdtemp()
    study39 = 'offline_stage_1_max_x'
    csv39 = os.path.join(outdir39, study39 + '_trajectory.csv')

    class _FakeTE39:
        k_1e = 47.1
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'K_1e', 'not_kinetic']
    fbs39 = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    V406_39 = SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0,
                              stage_1_max_x=5.0)
    seen39 = []   # V406.stage_1_max_x at every model_specification call
    def _model_specification39(**kw):
        seen39.append(V406_39.stage_1_max_x)
    def _solve_TEA39(stream_IDs=None):
        return {'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles39 = {
        'r_te': _FakeTE39(), 'fbs_spec': fbs39, 'V406': V406_39,
        'tea': SimpleNamespace(TCI=350e6), 'HXN': SimpleNamespace(),
        'model_specification': _model_specification39,
        'solve_TEA': _solve_TEA39,
        'latest_TEA_solution': {'IRR': np.nan,
                                'MPSPs': {'ethanol': np.nan,
                                          'isobutanol': np.nan}}}
    _optuna.logging.set_verbosity(_optuna.logging.WARNING)
    # 4 trials: 0 = baseline, 1 = the k_1e knockout probe (stage_1_max_x
    # stays at the baseline), 2-3 = sampled.
    study39_obj, csv39_out, kb39 = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=4, seed=1,
        study_name=study39, results_dir=outdir39, handles=handles39,
        print_status_every=1, burden_model=None,
        stage_1_max_x_bounds=(1.0, 50.0))
    df39 = ko.load_trajectory(csv39)
    assert 'stage_1_max_x' in df39.columns
    cols39 = list(df39.columns)
    assert cols39.index('state') < cols39.index('stage_1_max_x') < cols39.index('objective')
    assert df39['trial_number'].tolist() == [0, 1, 2, 3]
    assert df39['state'].tolist() == ['COMPLETE']*4
    assert df39['stage_1_max_x'][0] == 5.0                    # trial 0 = baseline
    assert df39['stage_1_max_x'][1] == 5.0                    # probe leaves it at baseline
    assert np.isclose(df39['k_1e'][1], 0.1*47.1)              # ...and knocks k_1e down
    assert ((df39['stage_1_max_x'] >= 1.0) & (df39['stage_1_max_x'] <= 50.0)).all()
    # The property was set BEFORE each simulation to the trial's value, and
    # the finally put the baseline back (one extra call at the baseline).
    assert len(seen39) == 5, seen39
    assert np.allclose(seen39[:4], df39['stage_1_max_x'].to_numpy(dtype=float),
                       rtol=1e-12, atol=0.0)
    assert seen39[4] == 5.0 and V406_39.stage_1_max_x == 5.0
    assert fbs39.max_n_spikes == 16
    # Without bounds the engine never touches the attribute (a fake V406
    # WITHOUT stage_1_max_x, as in check 17, keeps working).
    outdir39b = tempfile.mkdtemp()
    V406_39b = SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0)
    handles39b = dict(handles39, V406=V406_39b,
                      latest_TEA_solution={'IRR': np.nan,
                                           'MPSPs': {'ethanol': np.nan,
                                                     'isobutanol': np.nan}})
    _, csv39b, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=2, seed=1,
        study_name=study39 + '_off', results_dir=outdir39b, handles=handles39b,
        print_status_every=1, burden_model=None)
    assert 'stage_1_max_x' not in ko.load_trajectory(csv39b).columns
    assert not hasattr(V406_39b, 'stage_1_max_x')
    # restore_baseline: sets the property only when a value is given.
    V406_39.stage_1_max_x = 12.0
    ko.restore_baseline(handles39, kb39, fbs39.current_specifications,
                        baseline_max_n_spikes=16, baseline_stage_1_max_x=5.0)
    assert V406_39.stage_1_max_x == 5.0
    V406_39.stage_1_max_x = 12.0
    ko.restore_baseline(handles39, kb39, fbs39.current_specifications,
                        baseline_max_n_spikes=16)
    assert V406_39.stage_1_max_x == 12.0
    PASS('engine: stage_1_max_x sampled in-band, set on V406 before every simulation, baseline in trial 0 and the probe, restored in the finally; absent/untouched without bounds')

#%% 40. driver + supervisor: stage_1_max_x_bounds preset default, naming, --stage-1-max-x-bounds plumbing
drv40 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert ("'parameter_multiplier_bounds', 'stage_1_max_x_bounds',\n"
        "                    'parameter_groups', 'group_multiplier_bounds',\n"
        "                    'spike_delta_bounds'):") in drv40
assert "stage_1_max_x_bounds=engine_kwargs['stage_1_max_x_bounds']," in drv40
sup40 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
UNSET40 = sup40['_UNSET']
for fn40 in ('supervise', 'child_code', 'default_study_name'):
    assert _inspect.signature(sup40[fn40]).parameters['stage_1_max_x_bounds'].default is UNSET40, fn40
# child_code: kwarg omitted when unset (preset default), forwarded as
# None (pin) or a tuple (explicit band).
code40 = sup40['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic')
assert 'stage_1_max_x_bounds' not in code40
code40b = sup40['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', stage_1_max_x_bounds=None)
assert 'stage_1_max_x_bounds=None,' in code40b
code40c = sup40['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic', stage_1_max_x_bounds=[2.0, 30.0])
assert 'stage_1_max_x_bounds=(2.0, 30.0),' in code40c
# Naming mirror: unset -> the preset band tag; None -> no tag; explicit ->
# its own tag; legacy path never tagged.
assert sup40['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
assert sup40['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True,
                                   stage_1_max_x_bounds=None) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_burden'
assert sup40['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True,
                                   stage_1_max_x_bounds=(2.0, 30.0)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x2-30_burden'
assert sup40['default_study_name']('A', 'IRR', 'B', stage_1_max_x_bounds=(1.0, 50.0)) \
    == 'kin_opt_A_kbB_irr'
src40_sup = _inspect.getsource(sup40['supervise'])
assert 'stage_1_max_x_bounds=stage_1_max_x_bounds' in src40_sup
assert 'stage_1_max_x_bounds={stage_1_max_x_bounds!r}' in src40_sup   # settings event line
src40 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--stage-1-max-x-bounds', nargs='*', type=float" in src40
assert "parser.error('--stage-1-max-x-bounds takes" in src40
PASS('stage_1_max_x: driver setdefault + _s1x naming; supervisor _UNSET sentinel, child_code omit/None/tuple, naming mirror, --stage-1-max-x-bounds [LO HI]')

#%% 41. seed points from donor studies (2026-09-07): sim-free reader, clipping,
# column guard, _seed{n} tag; engine enqueues them after the probes of a
# FRESH study only; driver / supervisor plumbing (--seed-from).
kb41 = {'k_1e': 47.1, 'K_1e': 0.12}
space41, _ = ko.build_search_space(kb41, stage_1_max_x_bounds=(1.0, 50.0))
outdir41 = tempfile.mkdtemp()
donor41 = os.path.join(outdir41, 'donor41_trajectory.csv')
cols41 = ko.trajectory_columns(space41)
donor_row41 = dict(trial_number=7, state='COMPLETE', k_1e=100.0, K_1e=5.0,
                   threshold_conc=200.0, target_delta=20.0, spike_delta=300.0,
                   max_n_spikes=12, stage_1_max_x=3.0, objective=0.3)
ko.append_trajectory_row(donor41, cols41, donor_row41)
ko.append_trajectory_row(donor41, cols41, dict(donor_row41, trial_number=8,
                                               state='LOST', objective=''))
# clip_to_search_space: shared clipping rule (int cast, names outside the
# space dropped, moved names reported).
clipped41, moved41 = ko.clip_to_search_space(
    dict(k_1e=1e4, K_1e=0.5, max_n_spikes=60.0, not_a_var=1.0), space41)
assert clipped41 == dict(k_1e=471.0, K_1e=0.5, max_n_spikes=50)
assert moved41 == ['k_1e', 'max_n_spikes']
# Reader: label = the donor study name (file stem), K_1e clipped to its
# 10x ceiling, integer column cast, LOST row accepted with a note.
seeds41, notes41 = ko.seed_points_from_trajectory(donor41, [7, 8], space41)
assert list(seeds41) == ['donor41#7', 'donor41#8']
assert seeds41['donor41#7'] == dict(k_1e=100.0, K_1e=1.2, threshold_conc=200.0,
                                    target_delta=20.0, spike_delta=300.0,
                                    max_n_spikes=12, stage_1_max_x=3.0)
assert isinstance(seeds41['donor41#7']['max_n_spikes'], int)
assert any('donor41#7: clipped' in n and "['K_1e']" in n for n in notes41)
assert any("donor41#8: donor state 'LOST'" in n for n in notes41)
# Guards: absent trial; donor of another column set (no stage_1_max_x).
for bad41, msg41 in (([9], 'not in the trajectory'),):
    try:
        ko.seed_points_from_trajectory(donor41, bad41, space41)
    except ValueError as e41:
        assert msg41 in str(e41), e41
    else:
        raise AssertionError('absent seed trial did not raise')
space41_off, _ = ko.build_search_space(kb41)
donor41_off = os.path.join(outdir41, 'donor41off_trajectory.csv')
ko.append_trajectory_row(donor41_off, ko.trajectory_columns(space41_off),
                         dict(donor_row41))
try:
    ko.seed_points_from_trajectory(donor41_off, [7], space41)
except ValueError as e41:
    assert 'no decision column for' in str(e41) and 'stage_1_max_x' in str(e41)
else:
    raise AssertionError('donor of another column set did not raise')
# A donor with an EXTRA decision column (e.g. it sampled k_10) is fine:
# the extra column is dropped and noted.
seeds41_x, notes41_x = ko.seed_points_from_trajectory(donor41, [7], space41_off)
assert 'stage_1_max_x' not in seeds41_x['donor41#7']
assert any('ignored' in n and 'stage_1_max_x' in n for n in notes41_x)
# Naming tag: after _s1x, before _burden; nothing for 0 / None.
assert ko.seed_points_tag(3) == '_seed3' and ko.seed_points_tag(0) == '' \
    and ko.seed_points_tag(None) == ''
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic',
                             burden=True, rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.1, 10.0),
                             exclude_params=('k_10',),
                             stage_1_max_x_bounds=(1.0, 50.0), n_seeds=3) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_seed3_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic',
                             burden=True, n_seeds=0) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_burden'
# Engine: fresh study of 3 trials = baseline, the k_1e probe, the seed
# (donor given as a STUDY NAME resolved under results_dir); the seed's
# decision vector is the clipped donor point and carries the 'seed' attr.
if _optuna is None:
    print('SKIP 41 (engine part): optuna not installed')
else:
    class _FakeTE41:
        k_1e = 47.1
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'K_1e', 'not_kinetic']
    fbs41 = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    V406_41 = SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0,
                              stage_1_max_x=5.0)
    handles41 = {
        'r_te': _FakeTE41(), 'fbs_spec': fbs41, 'V406': V406_41,
        'tea': SimpleNamespace(TCI=350e6), 'HXN': SimpleNamespace(),
        'model_specification': lambda **kw: None,
        'solve_TEA': lambda stream_IDs=None: {
            'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}},
        'latest_TEA_solution': {'IRR': np.nan,
                                'MPSPs': {'ethanol': np.nan,
                                          'isobutanol': np.nan}}}
    _optuna.logging.set_verbosity(_optuna.logging.WARNING)
    study41, csv41, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=3, seed=1,
        study_name='offline_seeded', results_dir=outdir41, handles=handles41,
        print_status_every=1, burden_model=None,
        stage_1_max_x_bounds=(1.0, 50.0),
        seed_from=[('donor41', [7])])
    df41 = ko.load_trajectory(csv41)
    assert df41['trial_number'].tolist() == [0, 1, 2]
    assert df41['state'].tolist() == ['COMPLETE']*3
    assert df41['k_1e'][0] == 47.1 and np.isclose(df41['k_1e'][1], 4.71)
    assert df41['k_1e'][2] == 100.0 and df41['K_1e'][2] == 1.2
    assert df41['max_n_spikes'][2] == 12 and df41['stage_1_max_x'][2] == 3.0
    assert df41['threshold_conc'][2] == 200.0
    assert study41.trials[2].user_attrs.get('seed') == 'donor41#7'
    assert 'seed' not in study41.trials[1].user_attrs
    # Resume (one more trial): seeds are NOT re-enqueued -- trial 3 is a
    # sampled point without the attr.
    study41b, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=4, seed=1,
        study_name='offline_seeded', results_dir=outdir41, handles=handles41,
        print_status_every=1, burden_model=None,
        stage_1_max_x_bounds=(1.0, 50.0),
        seed_from=[('donor41', [7])])
    assert len(study41b.trials) == 4
    assert 'seed' not in study41b.trials[3].user_attrs
    assert ko.load_trajectory(csv41)['k_1e'][3] != 100.0
    # A bad donor fails BEFORE the store is opened (no .db / CSV created).
    try:
        ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
            objective='IRR', scenario_label='X', n_trials=3, seed=1,
            study_name='offline_seeded_bad', results_dir=outdir41,
            handles=handles41, burden_model=None,
            stage_1_max_x_bounds=(1.0, 50.0),
            seed_from=[('no_such_study', [7])])
    except ValueError as e41:
        assert 'no trajectory CSV' in str(e41)
    else:
        raise AssertionError('missing donor did not raise')
    assert not os.path.exists(os.path.join(outdir41, 'offline_seeded_bad.db'))
# Driver / supervisor plumbing.
drv41 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'seed_from=None,' in drv41 and 'seed_from=seed_from,' in drv41
assert 'n_seeds=n_seeds)' in drv41
assert "+ ko.seed_points_tag(n_seeds)" in drv41          # legacy-path name
sup41 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
seeds41_cli = [('donor_a', (1553, 1914)), ('donor_b', (1162,))]
assert sup41['seed_count'](seeds41_cli) == 3 and sup41['seed_count'](None) == 0
assert sup41['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic', burden=True,
                                   seed_from=seeds41_cli) \
    == 'kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_seed3_burden'
assert sup41['default_study_name']('A', 'IRR', 'B', seed_from=seeds41_cli) \
    == 'kin_opt_A_kbB_irr_seed3'
assert sup41['default_study_name']('A', 'IRR', 'B') == 'kin_opt_A_kbB_irr'
code41 = sup41['child_code'](None, 'IRR', 200, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic')
assert 'seed_from' not in code41
code41b = sup41['child_code'](None, 'IRR', 200, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic',
                              seed_from=[('donor_a', [1553, 1914]), ('donor_b', [1162])])
assert "seed_from=[('donor_a', (1553, 1914)), ('donor_b', (1162,))]," in code41b
src41_sup = _inspect.getsource(sup41['supervise'])
assert 'seed_from=seed_from' in src41_sup
assert 'seed_from={seed_from!r}' in src41_sup   # settings event line
src41 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO_supervised.py')).read()
assert "'--seed-from', nargs='+', action='append'" in src41
assert "parser.error('--seed-from takes STUDY TRIAL [TRIAL ...]')" in src41
PASS('seed points: clip_to_search_space, seed_points_from_trajectory (labels, clipping, int cast, LOST note, column guard, extra column dropped), _seed{n} tag, engine enqueues after the probes on a fresh study only / bad donor fails pre-store, driver + supervisor --seed-from plumbing')

#%% 42. parameter groups (2026-09-07, metabolic_minimal spec): one log-scale
# multiplier per group of kinetic parameters, members removed from the
# individual space; expand_grouped_values; spike_delta_bounds=None pins the
# spike (no column); baseline point 1.0 per group; no knockout probe for a
# group; every validation ValueError; plots tolerate a pinned spike.
assert ko.DEFAULT_GROUP_MULTIPLIER_BOUNDS == (0.2, 2.0)
assert ko.DEFAULT_SPIKE_DELTA_BOUNDS == (0.5, 595.0)
assert {'DEFAULT_GROUP_MULTIPLIER_BOUNDS', 'DEFAULT_SPIKE_DELTA_BOUNDS',
        'expand_grouped_values'} <= set(ko.__all__)
kb42 = {'k_1e': 47.1, 'k_1ie': 0.02, 'k_4ie': 0.04, 'k_1ia': 0.06,
        'k_10': 0.01, 'K_1e': 0.12}
groups42 = {'inhib_ethanol': ['k_1ie', 'k_4ie'], 'inhib_acetate': ['k_1ia']}
# Off by default: byte-identical to the pre-change space.
space42_off, excl42_off = ko.build_search_space(kb42)
assert list(space42_off) == ['k_1e', 'k_1ie', 'k_4ie', 'k_1ia', 'k_10', 'K_1e',
                             *ko.FEEDING_VARIABLES]
assert excl42_off == []
# Groups on: members gone, group entries after the kinetics and before the
# feeding variables (input order), band (0.2, 2.0) log; grouped members are
# NOT listed in `excluded` (they are sampled, through their group).
space42, excl42 = ko.build_search_space(kb42, parameter_groups=groups42,
                                        exclude_params=('k_10',),
                                        spike_delta_bounds=None,
                                        stage_1_max_x_bounds=(1.0, 50.0))
assert list(space42) == ['k_1e', 'K_1e', 'inhib_ethanol', 'inhib_acetate',
                         'threshold_conc', 'target_delta', 'max_n_spikes',
                         'stage_1_max_x'], list(space42)
assert space42['inhib_ethanol'] == dict(low=0.2, high=2.0, log=True)
assert space42['inhib_acetate'] == dict(low=0.2, high=2.0, log=True)
assert excl42 == ['k_10']
assert 'spike_delta' not in space42
# A list of pairs is accepted like a dict; a custom band applies to every group.
space42_pairs, _ = ko.build_search_space(
    kb42, parameter_groups=[('inhib_ethanol', ('k_1ie', 'k_4ie'))],
    group_multiplier_bounds=(0.5, 3.0))
assert space42_pairs['inhib_ethanol'] == dict(low=0.5, high=3.0, log=True)
assert 'k_1ia' in space42_pairs and 'k_1ie' not in space42_pairs
# A grouped member's param_bounds_override entry is IGNORED (the driver's
# preset path passes absolute workbook bounds for every row); an
# include_params whitelist that omits a member still groups it.
space42_ov, _ = ko.build_search_space(
    kb42, parameter_groups=groups42, include_params=['k_1e'],
    param_bounds_override={'k_1ie': (0.001, 0.1), 'k_1e': (1.0, 100.0)})
assert list(space42_ov) == ['k_1e', 'inhib_ethanol', 'inhib_acetate',
                            *ko.FEEDING_VARIABLES]
assert space42_ov['k_1e'] == dict(low=1.0, high=100.0, log=True)
# expand_grouped_values: member = baseline x multiplier, group key dropped,
# everything else passed through; identity copy without groups.
vals42 = {'k_1e': 50.0, 'K_1e': 0.1, 'inhib_ethanol': 0.5, 'inhib_acetate': 2.0,
          'threshold_conc': 100.0, 'target_delta': 50.0, 'max_n_spikes': 3,
          'stage_1_max_x': 5.0}
exp42 = ko.expand_grouped_values(vals42, groups42, kb42)
assert exp42 == {'k_1e': 50.0, 'K_1e': 0.1, 'k_1ie': 0.01, 'k_4ie': 0.02,
                 'k_1ia': 0.12, 'threshold_conc': 100.0, 'target_delta': 50.0,
                 'max_n_spikes': 3, 'stage_1_max_x': 5.0}, exp42
assert 'inhib_ethanol' not in exp42
ident42 = ko.expand_grouped_values(vals42, None, kb42)
assert ident42 == vals42 and ident42 is not vals42
assert ko.expand_grouped_values(vals42, {}, kb42) == vals42
# Baseline point: 1.0 per group, members absent, no spike_delta when pinned.
spec42 = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
base42 = ko.baseline_decision_point(space42, kb42, spec42, 16,
                                    baseline_stage_1_max_x=5.0,
                                    parameter_groups=groups42)
assert base42 == {'k_1e': 47.1, 'K_1e': 0.12, 'inhib_ethanol': 1.0,
                  'inhib_acetate': 1.0, 'threshold_conc': 217.125,
                  'target_delta': 5.0,  # 221.25-217.125=4.125 -> clipped to
                                        # the default target_delta low bound
                                        # 5.0 (same clip as check 10)
                  'max_n_spikes': 16,
                  'stage_1_max_x': 5.0}, base42
# ... and unchanged without groups (spike_delta present as before).
base42_off = ko.baseline_decision_point(space42_off, kb42, spec42, 16)
assert base42_off['spike_delta'] == 600.0 - 221.25 and 'inhib_ethanol' not in base42_off
# Knockout probes: a group is not a rate constant -- no probe under either
# rule; the k_1e probe carries the groups at 1.0.
probes42, at_floor42 = ko.knockout_probe_points(space42, base42)
assert list(probes42) == ['k_1e'] and at_floor42 == []
assert probes42['k_1e']['inhib_ethanol'] == 1.0
probes42r, _ = ko.knockout_probe_points(space42, base42, rate_params=['k_1e', 'inhib_ethanol'])
assert list(probes42r) == ['k_1e', 'inhib_ethanol']   # explicit rate_params wins (never the presets' case)
# Trajectory columns: the group IS a decision column (between state and objective).
cols42 = ko.trajectory_columns(space42)
assert cols42[2:2 + len(space42)] == list(space42)
# Validation ValueErrors, each naming the offender.
def _raises42(msg, **kw):
    try:
        ko.build_search_space(kb42, **kw)
    except ValueError as e:
        assert msg in str(e), (msg, str(e))
    else:
        raise AssertionError(f'no ValueError for {kw} (expected {msg!r})')
_raises42('k_1e', parameter_groups={'k_1e': ['k_1ie']})                    # collides with a kinetic name
_raises42('threshold_conc', parameter_groups={'threshold_conc': ['k_1ie']})  # ... a feeding variable
_raises42('stage_1_max_x', parameter_groups={'stage_1_max_x': ['k_1ie']})   # ... an operating variable
_raises42('k_9', parameter_groups={'g': ['k_9']})                          # unknown member
_raises42('k_1ie', parameter_groups={'g1': ['k_1ie'], 'g2': ['k_1ie']})    # member in two groups
_raises42('k_1ie', parameter_groups={'g': ['k_1ie']}, exclude_params=('k_1ie',))  # excluded member
_raises42('g', parameter_groups={'g': []})                                 # empty group
_raises42('k_z', parameter_groups={'g': ['k_z']})                          # nonpositive baseline (below)
for bad42 in ((0.0, 2.0), (-1.0, 2.0), (2.0, 0.2), (1.0, 1.0)):
    _raises42('group_multiplier_bounds', parameter_groups=groups42,
              group_multiplier_bounds=bad42)
_raises42('spike_delta_bounds', spike_delta_bounds=None, target_conc_bounds=(180.0, 300.0))  # legacy mix
try:
    ko.build_search_space({**kb42, 'k_z': 0.0}, parameter_groups={'g': ['k_z']})
except ValueError as e42:
    assert 'k_z' in str(e42) and 'nonpositive' in str(e42)
else:
    raise AssertionError('nonpositive grouped baseline did not raise')
# Plots on a synthetic trajectory WITHOUT spike_delta (pinned spike): the
# applied spike is NaN, the spike panel is omitted, nothing raises.
df42 = pd.DataFrame([dict(trial_number=i, state='COMPLETE', k_1e=47.1*(1 + i),
                          K_1e=0.12, inhib_ethanol=1.0 + 0.1*i, inhib_acetate=1.0,
                          threshold_conc=200.0, target_delta=20.0, max_n_spikes=5,
                          stage_1_max_x=5.0, objective=0.1*(i + 1),
                          **{m: np.nan for m in ko.TRACKED_METRICS},
                          error='') for i in range(3)])
t42, th42, sp42 = ko._applied_feeding(df42)
assert np.allclose(t42, 220.0) and np.allclose(th42, 200.0) and np.isnan(sp42).all()
t42s, th42s, sp42s = ko._applied_feeding(df42.iloc[0])
assert t42s == 220.0 and np.isnan(sp42s)
fig42, axes42 = ko.plot_optimization_trajectories(df42, 'IRR', 'maximize')
assert not any(ax.get_title().startswith('spike_conc') for ax in axes42.ravel())
assert any(ax.get_title().startswith('target_conc') for ax in axes42.ravel())
fig42b, (ax42b1, ax42b2) = ko.plot_parameter_trajectory(
    df42, {**kb42, 'inhib_ethanol': 1.0, 'inhib_acetate': 1.0}, 'maximize')
labels42 = [line.get_label() for line in ax42b1.get_lines()]
assert 'inhib_ethanol' in labels42 and 'k_1e' in labels42 and 'k_1ie' not in labels42
assert 'spike_conc' not in [line.get_label() for line in ax42b2.get_lines()]
fig42c, ax42c = ko.plot_best_vs_baseline(
    df42, {**kb42, 'inhib_ethanol': 1.0, 'inhib_acetate': 1.0}, 'maximize')
assert 'spike = pinned' in ax42c.get_title(), ax42c.get_title()
import matplotlib.pyplot as _plt42
_plt42.close('all')
PASS('parameter groups: one log multiplier per group after the kinetics, members removed, override ignored, expand_grouped_values, baseline 1.0, no probe, spike_delta_bounds=None pins the spike, every validation ValueError, plots tolerate a pinned spike')

#%% 43. metabolic_minimal preset (2026-09-07): 17 rates + one multiplier per
# inhibition-effector family + 4 feeding/operating variables; effector
# table by file path; study_type_name_defaults shared by driver and
# supervisor; the exact study name (the _ib tag from the GROUP band, so an
# explicit band gets its own study); the ethanol_only 19-variable space;
# existing presets untouched.
assert ko.EFFECTOR_ORDER == ('ethanol', 'isobutanol', 'acetate')
assert ko.STUDY_TYPE_OPTIONS == {
    'metabolic_minimal': dict(exclude_params=('k_10', 'k_7', 'k_8'),
                              group_roles=('product_inhibition', 'lethality'),
                              group_multiplier_bounds=(0.2, 2.0),
                              spike_delta_bounds=None),
    'metabolic_minimal_subset': dict(rate_params=ko.METABOLIC_MINIMAL_SUBSET_RATES,
                                     parameter_groups=ko.METABOLIC_MINIMAL_SUBSET_GROUPS,
                                     group_multiplier_bounds=(0.2, 2.0),
                                     exclude_params=(),
                                     spike_delta_bounds=None,
                                     stage_1_max_x_bounds=None)}
assert {'STUDY_TYPE_OPTIONS', 'EFFECTOR_ORDER', 'kinetic_parameter_effectors',
        'study_type_name_defaults'} <= set(ko.__all__)
# Naming defaults per study type (workbook-free): the minimal type's group
# band stands in for the inhibition band (_ib0.2-2) and its exclusion set
# is k_10 + k_7 + k_8; every other type keeps the module defaults.
assert ko.study_type_name_defaults('metabolic_minimal') == dict(
    inhibition_multiplier_bounds=(0.2, 2.0), exclude_params=('k_10', 'k_7', 'k_8'),
    stage_1_max_x_bounds=(1.0, 50.0))
for st43 in ('metabolic', 'metabolic_protein'):
    assert ko.study_type_name_defaults(st43) == dict(
        inhibition_multiplier_bounds=ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS,
        exclude_params=ko.DEFAULT_EXCLUDED_PARAMETERS,
        stage_1_max_x_bounds=(1.0, 50.0))
try:
    ko.study_type_name_defaults('protein')
except ValueError as e43:
    assert 'protein' in str(e43)
else:
    raise AssertionError('unknown study_type did not raise')
NAME43 = ('kin_opt_ethanol_isobutanol_metabolic_minimal_irr'
          '_rb0.001-10_ib0.2-2_xk10+k7+k8_s1x1-50_burden')
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_minimal',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.2, 2.0),
                             exclude_params=('k_10', 'k_7', 'k_8'),
                             stage_1_max_x_bounds=(1.0, 50.0), burden=True) == NAME43
# Driver: the three new keys are setdefault'ed like the others and the
# preset summary reports the groups / pinned spike (source-text pins, the
# driver cannot be imported offline).
drv43 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert "'parameter_groups', 'group_multiplier_bounds'," in drv43
assert "'spike_delta_bounds'):" in drv43
assert 'engine_kwargs.setdefault(key, preset[key])' in drv43
assert "Parameter groups" in drv43 and 'spike feed pinned at the baseline' in drv43
# The _ib tag follows the band that actually SIZES the inhibition entries:
# group_multiplier_bounds under a grouped study type (the members are not
# sampled individually), multiplier_bounds otherwise.
assert 'inhibition_multiplier_bounds=(' in drv43
assert "engine_kwargs['group_multiplier_bounds']" in drv43
assert "if engine_kwargs['parameter_groups'] else" in drv43
assert "plot_baselines" in drv43                    # group multipliers plotted at baseline 1.0
assert "| set(engine_kwargs.get('parameter_groups') or ())" in drv43   # PCA log columns
assert 'metabolic_minimal' in drv43
# Supervisor: the _ib / _x tags come from ko.study_type_name_defaults, so
# its derived name equals the driver's under metabolic_minimal (it used to
# hardcode DEFAULT_SATURATION_MULTIPLIER_BOUNDS / DEFAULT_EXCLUDED_PARAMETERS).
sup43 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal', burden=True) == NAME43
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal', burden=True,
                                   exclude_params=('k_10',)) \
    == 'kin_opt_ethanol_isobutanol_metabolic_minimal_irr_rb0.001-10_ib0.2-2_xk10_s1x1-50_burden'
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
src43_sup = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'optimize_kinetics_BO_supervised.py')).read()
assert 'ko.study_type_name_defaults(study_type)' in src43_sup
assert 'ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS' not in _inspect.getsource(sup43['default_study_name'])
assert 'metabolic_minimal' in src43_sup
# Effector table read by FILE PATH, cached like the roles; a table row
# without an effector gives None.
eff43 = ko.kinetic_parameter_effectors()
assert eff43 is ko.kinetic_parameter_effectors()            # cached
assert set(eff43) == set(ko.kinetic_parameter_roles())
assert eff43['k_1ie'] == 'ethanol' and eff43['k_1ii'] == 'isobutanol' \
    and eff43['k_16ia'] == 'acetate' and eff43['k_10ie'] == 'ethanol'
assert eff43['k_1e'] is None and eff43['k_7'] is None
assert ko.kinetic_parameter_effectors(ko.kinetic_parameter_roles_path()) == eff43
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    roles43 = ko.kinetic_parameter_roles()
    p43 = ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal')
    assert p43['scenario'] == 'A' and p43['kinetic_bounds_scenario'] == 'B'
    assert p43['exclude_params'] == ('k_10', 'k_7', 'k_8')
    assert p43['multiplier_bounds'] == (0.2, 2.0)             # -> the _ib0.2-2 tag
    assert p43['group_multiplier_bounds'] == (0.2, 2.0)
    assert p43['spike_delta_bounds'] is None
    assert p43['rate_multiplier_bounds'] == ko.DEFAULT_RATE_MULTIPLIER_BOUNDS
    assert p43['parameter_multiplier_bounds'] == {'k_10': (0.1, 10.0)}
    assert p43['stage_1_max_x_bounds'] == (1.0, 50.0)
    assert p43['parameter_groups'] == {
        'inhib_ethanol': ['k_1ie', 'k_4ie', 'k_7ie', 'k_10ie', 'k_16ie'],
        'inhib_isobutanol': ['k_1ii', 'k_4ii', 'k_6ii', 'k_7ii', 'k_10ii'],
        'inhib_acetate': ['k_1ia', 'k_4ia', 'k_6ia', 'k_7ia', 'k_10ia', 'k_16ia']}
    assert list(p43['parameter_groups']) == ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate']
    inc43 = p43['include_params']
    assert len(inc43) == 36                                    # 20 capacities + 16 inhibition rows
    assert all(roles43[n] in ('capacity', 'product_inhibition', 'lethality') for n in inc43)
    assert not any(n.startswith('K_') for n in inc43)         # K_1i etc. are OUT
    grouped43 = {m for ms in p43['parameter_groups'].values() for m in ms}
    individual43 = [n for n in inc43 if n not in grouped43 and n not in p43['exclude_params']]
    assert len(individual43) == 17 and all(roles43[n] == 'capacity' for n in individual43)
    assert not {'k_10', 'k_7', 'k_8'} & set(individual43)
    assert len(p43['rate_params']) == 20                       # the workbook's capacities, as before
    # The resulting space (live baselines = the B workbook values here):
    # 17 + 3 + 4 = 24 decision variables, in the documented order.
    kb43 = ko.workbook_kinetic_baselines('B')
    space43, excl43 = ko.build_search_space(
        kb43, include_params=inc43, exclude_params=p43['exclude_params'],
        rate_multiplier_bounds=p43['rate_multiplier_bounds'],
        rate_params=p43['rate_params'],
        parameter_multiplier_bounds=p43['parameter_multiplier_bounds'],
        parameter_groups=p43['parameter_groups'],
        group_multiplier_bounds=p43['group_multiplier_bounds'],
        spike_delta_bounds=p43['spike_delta_bounds'],
        stage_1_max_x_bounds=p43['stage_1_max_x_bounds'])
    assert len(space43) == 24, len(space43)
    assert list(space43)[:17] == individual43
    assert list(space43)[17:] == ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate',
                                  'threshold_conc', 'target_delta', 'max_n_spikes',
                                  'stage_1_max_x']
    assert set(excl43) == set(kb43) - set(individual43) - grouped43
    # ethanol_only: 13 rates (16 - 3), 2 groups (no isobutanol rows), 4 = 19.
    p43_eo = ko.resolve_study_preset('ethanol_only', 'metabolic_minimal')
    assert p43_eo['parameter_groups'] == {
        'inhib_ethanol': ['k_1ie', 'k_4ie', 'k_7ie', 'k_10ie'],
        'inhib_acetate': ['k_1ia', 'k_4ia', 'k_6ia', 'k_7ia', 'k_10ia']}
    inc43_eo = p43_eo['include_params']
    grouped43_eo = {m for ms in p43_eo['parameter_groups'].values() for m in ms}
    individual43_eo = [n for n in inc43_eo if n not in grouped43_eo
                       and n not in p43_eo['exclude_params']]
    assert len(individual43_eo) == 13
    # ... and the space it builds (live baselines = the A workbook values
    # here): 13 rates + 2 group multipliers + 4 feeding/operating = 19,
    # spike_delta pinned out.
    kb43_eo = ko.workbook_kinetic_baselines('A')
    space43_eo, excl43_eo = ko.build_search_space(
        kb43_eo, include_params=inc43_eo,
        exclude_params=p43_eo['exclude_params'],
        rate_multiplier_bounds=p43_eo['rate_multiplier_bounds'],
        rate_params=p43_eo['rate_params'],
        parameter_multiplier_bounds=p43_eo['parameter_multiplier_bounds'],
        parameter_groups=p43_eo['parameter_groups'],
        group_multiplier_bounds=p43_eo['group_multiplier_bounds'],
        spike_delta_bounds=p43_eo['spike_delta_bounds'],
        stage_1_max_x_bounds=p43_eo['stage_1_max_x_bounds'])
    assert len(space43_eo) == 19, len(space43_eo)
    assert list(space43_eo)[:13] == individual43_eo
    assert list(space43_eo)[13:] == ['inhib_ethanol', 'inhib_acetate',
                                     'threshold_conc', 'target_delta',
                                     'max_n_spikes', 'stage_1_max_x']
    assert 'spike_delta' not in space43_eo
    assert set(excl43_eo) == set(kb43_eo) - set(individual43_eo) - grouped43_eo
    # The four existing presets: parameter_groups None, everything else as
    # in check 21 (their multiplier_bounds / exclude_params come through
    # study_type_name_defaults now, same values).
    for stp43, st43 in (('ethanol_only', 'metabolic'), ('ethanol_only', 'metabolic_protein'),
                        ('ethanol_isobutanol', 'metabolic'),
                        ('ethanol_isobutanol', 'metabolic_protein')):
        q43 = ko.resolve_study_preset(stp43, st43)
        assert q43['parameter_groups'] is None
        assert q43['multiplier_bounds'] == (0.1, 10.0) and q43['exclude_params'] == ('k_10',)
        assert q43['spike_delta_bounds'] == (0.5, 595.0)
    # A grouped row with no effector in the (injected) table must raise.
    eff_missing43 = dict(eff43); eff_missing43['k_1ie'] = None
    try:
        ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal',
                                effectors=eff_missing43)
    except KeyError as e43:
        assert 'k_1ie' in str(e43)
    else:
        raise AssertionError('grouped row without an effector did not raise KeyError')
    # Supervisor name == driver name (the driver's exact default_study_name
    # call on the preset's effective values, _driver_name43 mirroring its
    # inhibition-band choice) for all six presets.
    def _driver_name43(stp, st, **overrides):
        """The driver's derived study name for a preset, with `overrides`
        standing in for explicit run() engine kwargs (the driver
        setdefault's the preset into engine_kwargs, then names from the
        EFFECTIVE values)."""
        kw = dict(ko.resolve_study_preset(stp, st))
        kw.update(overrides)
        return ko.default_study_name(
            'IRR', stp, st, scenario=kw['scenario'],
            kinetic_bounds_scenario=kw['kinetic_bounds_scenario'], burden=True,
            rate_multiplier_bounds=kw['rate_multiplier_bounds'],
            # the band that sizes the inhibition entries: the group band
            # under a grouped preset, multiplier_bounds otherwise
            inhibition_multiplier_bounds=(kw['group_multiplier_bounds']
                                          if kw['parameter_groups'] else
                                          kw['multiplier_bounds']),
            exclude_params=kw['exclude_params'],
            stage_1_max_x_bounds=kw['stage_1_max_x_bounds'], n_seeds=0)
    for stp43, st43 in ((a, b) for a in ko.STUDY_TARGET_PRODUCTS for b in ko.STUDY_TYPE_ROLES):
        driver_name43 = _driver_name43(stp43, st43)
        assert sup43['default_study_name'](None, 'IRR', None, study_target_products=stp43,
                                           study_type=st43, burden=True) == driver_name43, (stp43, st43)
    # An EXPLICIT group band must get its own study: it sizes the same
    # columns, and optuna accepts a changed numeric range on a resume, so
    # only the name keeps a 0.5x-3x run off the preset's 0.2x-2x store.
    assert _driver_name43('ethanol_isobutanol', 'metabolic_minimal') == NAME43
    name43_wide = _driver_name43('ethanol_isobutanol', 'metabolic_minimal',
                                 group_multiplier_bounds=(0.5, 3.0))
    assert '_ib0.5-3' in name43_wide and name43_wide != NAME43, name43_wide
    # ... while an ungrouped preset ignores group_multiplier_bounds entirely.
    assert _driver_name43('ethanol_isobutanol', 'metabolic_protein',
                          group_multiplier_bounds=(0.5, 3.0)) \
        == _driver_name43('ethanol_isobutanol', 'metabolic_protein')
    PASS('metabolic_minimal preset: 3 effector groups (5/5/6) + 17 rates + 4 = 24 and the ethanol_only space built too (2 groups + 13 + 4 = 19, no spike_delta), K_* out, spike pinned, _ib0.2-2 / _xk10+k7+k8 naming via study_type_name_defaults (the _ib tag reads the GROUP band under a grouped preset, so an explicit 0.5-3 renames the study), effector table by file path, existing presets untouched')
else:
    print('SKIP 43 (preset part): parameter-distribution workbooks not found')
    PASS('metabolic_minimal naming + effector table (workbook-free part)')

#%% 44. engine end to end with a parameter group and a pinned spike (the
# check-17 / 39 fake-handle pattern): members set on r_te = baseline x
# multiplier, excluded k_10 untouched, spike_conc = the baseline at every
# call, applied_* columns after the metrics (before 'error'), sidecar
# carries them, baseline restored in the finally; seeding a space that
# samples the members individually resolves them through applied_*, the
# reverse raises; and (burden x groups) a fake burden_model's evaluate/apply
# receive the EXPANDED member values, never the group multiplier, with the
# BURDEN_COLUMNS ahead of the applied_* columns in the header.
if _optuna is None:
    print('SKIP 44: optuna not installed')
else:
    outdir44 = tempfile.mkdtemp()
    study44 = 'offline_grouped'
    csv44 = os.path.join(outdir44, study44 + '_trajectory.csv')

    class _FakeTE44:
        k_1e = 47.1
        k_1ie = 0.02
        k_4ie = 0.04
        k_10 = 0.01
        K_1e = 0.12
        def getGlobalParameterIds(self):
            return ['k_1e', 'k_1ie', 'k_4ie', 'k_10', 'K_1e', 'not_kinetic']
    te44 = _FakeTE44()
    fbs44 = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    seen44 = []   # (spike_conc kwarg, k_1ie, k_4ie, k_10) at every model_specification call
    def _model_specification44(**kw):
        seen44.append((kw['spike_conc'], te44.k_1ie, te44.k_4ie, te44.k_10))
    handles44 = {
        'r_te': te44, 'fbs_spec': fbs44,
        'V406': SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0),
        'tea': SimpleNamespace(TCI=350e6), 'HXN': SimpleNamespace(),
        'model_specification': _model_specification44,
        'solve_TEA': lambda stream_IDs=None: {
            'IRR': 0.2, 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}},
        'latest_TEA_solution': {'IRR': np.nan,
                                'MPSPs': {'ethanol': np.nan,
                                          'isobutanol': np.nan}}}
    groups44 = {'inhib_ethanol': ['k_1ie', 'k_4ie']}
    _optuna.logging.set_verbosity(_optuna.logging.WARNING)
    # 3 trials: 0 = baseline (multiplier 1.0), 1 = the k_1e probe, 2 = sampled.
    study44_obj, csv44_out, kb44 = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=3, seed=1,
        study_name=study44, results_dir=outdir44, handles=handles44,
        print_status_every=1, burden_model=None,
        exclude_params=('k_10',), parameter_groups=groups44,
        spike_delta_bounds=None)
    assert kb44 == {'k_1e': 47.1, 'k_1ie': 0.02, 'k_4ie': 0.04, 'k_10': 0.01, 'K_1e': 0.12}
    df44 = ko.load_trajectory(csv44)
    cols44 = list(df44.columns)
    assert cols44 == ['trial_number', 'state', 'k_1e', 'K_1e', 'inhib_ethanol',
                      'threshold_conc', 'target_delta', 'max_n_spikes', 'objective',
                      *ko.TRACKED_METRICS, 'applied_k_1ie', 'applied_k_4ie', 'error'], cols44
    assert df44['trial_number'].tolist() == [0, 1, 2]
    assert df44['state'].tolist() == ['COMPLETE']*3
    assert df44['inhib_ethanol'][0] == 1.0 and df44['inhib_ethanol'][1] == 1.0
    assert np.isclose(df44['k_1e'][1], 0.1*47.1)              # the k_1e probe at the floor of the default 0.1x band
    assert 0.2 <= df44['inhib_ethanol'][2] <= 2.0
    # The model saw baseline x multiplier for every member, the excluded
    # k_10 untouched, the baseline spike at every call (3 trials + restore).
    assert len(seen44) == 4, seen44
    for i44 in range(3):
        spike44, k1ie44, k4ie44, k10_44 = seen44[i44]
        assert spike44 == 600.0
        assert np.isclose(k1ie44, 0.02*df44['inhib_ethanol'][i44], rtol=1e-12, atol=0.0)
        assert np.isclose(k4ie44, 0.04*df44['inhib_ethanol'][i44], rtol=1e-12, atol=0.0)
        assert k10_44 == 0.01
        assert np.isclose(df44['applied_k_1ie'][i44], k1ie44, rtol=1e-12, atol=0.0)
        assert np.isclose(df44['applied_k_4ie'][i44], k4ie44, rtol=1e-12, atol=0.0)
    assert seen44[3] == (600.0, 0.02, 0.04, 0.01)             # restore_baseline
    assert te44.k_1ie == 0.02 and te44.k_4ie == 0.04 and te44.k_1e == 47.1
    # The sidecar written before each simulation carried the applied_*
    # columns too (so a LOST row is complete): replay one trial's record
    # through write_inflight/recover_inflight with the engine's columns.
    side44 = ko.inflight_path_for(outdir44, study44)
    cols44_engine = ko.trajectory_columns(
        ko.build_search_space(kb44, exclude_params=('k_10',),
                              parameter_groups=groups44, spike_delta_bounds=None)[0],
        extra_columns=['applied_k_1ie', 'applied_k_4ie'])
    assert cols44_engine == cols44
    ko.write_inflight(side44, cols44_engine,
                      {'trial_number': 3, 'k_1e': 47.1, 'K_1e': 0.12, 'inhib_ethanol': 0.5,
                       'threshold_conc': 200.0, 'target_delta': 20.0, 'max_n_spikes': 2,
                       'applied_k_1ie': 0.01, 'applied_k_4ie': 0.02})
    assert ko.recover_inflight(csv44, side44, state='LOST', error='x') == 3
    df44b = ko.load_trajectory(csv44)
    assert df44b['state'].tolist()[-1] == 'LOST' and df44b['applied_k_1ie'].tolist()[-1] == 0.01
    # Seeds: minimal -> individually sampled members (same pinned spike)
    # resolve k_1ie / k_4ie through applied_*; the reverse direction (a
    # donor that sampled the members individually into a grouped space)
    # raises the column error, naming the group.
    space44_full, _ = ko.build_search_space(kb44, exclude_params=('k_10',),
                                            spike_delta_bounds=None)
    assert {'k_1ie', 'k_4ie'} <= set(space44_full)
    seeds44, notes44 = ko.seed_points_from_trajectory(csv44, [2], space44_full)
    pt44 = seeds44[f'{study44}#2']
    assert np.isclose(pt44['k_1ie'], df44['applied_k_1ie'][2]) \
        and np.isclose(pt44['k_4ie'], df44['applied_k_4ie'][2])
    # (isclose, not ==: the reader round-trips the CSV cell exactly through
    # float(), while pandas' default parser can land 1 ULP away.)
    assert np.isclose(pt44['k_1e'], df44['k_1e'][2], rtol=1e-12, atol=0.0) \
        and np.isclose(pt44['threshold_conc'], df44['threshold_conc'][2],
                       rtol=1e-12, atol=0.0)
    assert any('applied_' in n and 'k_1ie' in n for n in notes44), notes44
    donor44_full = os.path.join(outdir44, 'full44_trajectory.csv')
    ko.append_trajectory_row(donor44_full, ko.trajectory_columns(space44_full),
                             dict(trial_number=1, state='COMPLETE', k_1e=47.1, K_1e=0.12,
                                  k_1ie=0.02, k_4ie=0.04, threshold_conc=200.0,
                                  target_delta=20.0, max_n_spikes=2, objective=0.1))
    space44_min, _ = ko.build_search_space(kb44, exclude_params=('k_10',),
                                           parameter_groups=groups44, spike_delta_bounds=None)
    try:
        ko.seed_points_from_trajectory(donor44_full, [1], space44_min)
    except ValueError as e44:
        assert 'inhib_ethanol' in str(e44) and 'no decision column' in str(e44)
        assert 'no unique inverse' in str(e44)
    else:
        raise AssertionError('full -> grouped seeding did not raise')
    # spike_delta has no applied_ column, so the pin only seeds one way:
    # a donor that SAMPLED spike_delta into a pinned-spike space drops the
    # column (noted); a pinned-spike donor into a space that samples it
    # raises (the documented limitation).
    space44_spk_full, _ = ko.build_search_space(kb44, exclude_params=('k_10',))  # members individual, spike_delta sampled
    donor44_spk = os.path.join(outdir44, 'spk44_trajectory.csv')
    ko.append_trajectory_row(donor44_spk, ko.trajectory_columns(space44_spk_full),
                             dict(trial_number=1, state='COMPLETE', k_1e=47.1, K_1e=0.12,
                                  k_1ie=0.02, k_4ie=0.04, threshold_conc=200.0,
                                  target_delta=20.0, spike_delta=300.0, max_n_spikes=2,
                                  objective=0.1))
    seeds44_drop, notes44_drop = ko.seed_points_from_trajectory(
        donor44_spk, [1], space44_full)
    assert 'spike_delta' not in seeds44_drop['spk44#1']
    assert seeds44_drop['spk44#1']['k_1ie'] == 0.02
    assert any('ignored' in n and 'spike_delta' in n for n in notes44_drop)
    space44_spk, _ = ko.build_search_space(kb44, exclude_params=('k_10',),
                                           parameter_groups=groups44)   # spike_delta sampled
    try:
        ko.seed_points_from_trajectory(csv44, [2], space44_spk)
    except ValueError as e44:
        assert 'spike_delta' in str(e44)
        # ... and the message says the DONOR pinned it (a feeding variable
        # is never a column), not that a group multiplier has no inverse.
        assert 'PINNED' in str(e44) and 'no unique inverse' not in str(e44), e44
    else:
        raise AssertionError('pinned-spike donor into a spike-sampling space did not raise')
    # Burden x groups: the burden model must be called with the EXPANDED
    # member values (baseline x sampled multiplier), never the raw group
    # multiplier -- a fake burden_model captures exactly what evaluate()/
    # apply() receive at every trial. Fresh handles/study (te44/fbs44 above
    # were already exercised and restored, but a fresh fixture avoids any
    # cross-run coupling).
    outdir44g = tempfile.mkdtemp()
    study44g = 'offline_grouped_burden'
    csv44g = os.path.join(outdir44g, study44g + '_trajectory.csv')

    class _FakeBurden44:
        # Minimal stand-in for enzyme_burden.BurdenModel: the engine reads
        # F_flex/Phi_M_wt/phi_T_wt unconditionally at start-up (the 'Enzyme
        # burden ON' banner) and `reference` for the stale-snapshot guard
        # (empty -> nothing to compare -> never stale); every trial is kept
        # feasible so evaluate()/apply() run on every one of the 3 trials.
        F_flex = 0.245
        Phi_M_wt = 0.05
        phi_T_wt = 0.05
        reference = {}
        def __init__(self):
            self.seen_evaluate = []   # dicts passed to evaluate(), in order
            self.seen_apply = []      # dicts passed to apply(), in order
        def evaluate(self, values):
            self.seen_evaluate.append(dict(values))
            record = {col: 0.0 for col in eb.BURDEN_COLUMNS}
            return SimpleNamespace(feasible=True, violation=-1.0, Phi_M=0.01,
                                   F_flex=self.F_flex, k_7_eff=0.0, k_8_eff=0.0,
                                   as_record=lambda: record)
        def apply(self, values):
            self.seen_apply.append(dict(values))
            return dict(values)   # no k_7/k_8 sampled here -- pass through
    burden44g = _FakeBurden44()
    te44g = _FakeTE44()
    fbs44g = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    seen44g = []
    def _model_specification44g(**kw):
        seen44g.append((kw['spike_conc'], te44g.k_1ie, te44g.k_4ie, te44g.k_10))
    handles44g = dict(handles44, r_te=te44g, fbs_spec=fbs44g,
                      model_specification=_model_specification44g,
                      latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    # feasible_sampling=False: burden on + feasible_sampling on (the
    # default, check 35) also calls burden_model.evaluate() from the
    # FeasibleTPESampler's own start-up/candidate feasibility checks, which
    # would inflate seen_evaluate with calls unrelated to _objective; off
    # here isolates the one evaluate()/apply() pair per trial this check
    # targets (the objective's own hook, unconditionally exercised either way).
    study44g_obj, csv44g_out, kb44g = ko.run_kinetic_optimization(enqueue_baseline=True, enqueue_knockouts=True,
        objective='IRR', scenario_label='X', n_trials=3, seed=1,
        study_name=study44g, results_dir=outdir44g, handles=handles44g,
        print_status_every=1, burden_model=burden44g, feasible_sampling=False,
        exclude_params=('k_10',), parameter_groups=groups44,
        spike_delta_bounds=None)
    assert csv44g_out == csv44g and kb44g == kb44
    df44g = ko.load_trajectory(csv44g)
    cols44g = list(df44g.columns)
    space44g, _ = ko.build_search_space(kb44, exclude_params=('k_10',),
                                        parameter_groups=groups44,
                                        spike_delta_bounds=None)
    assert cols44g == ko.trajectory_columns(
        space44g, extra_columns=[*eb.BURDEN_COLUMNS, 'applied_k_1ie', 'applied_k_4ie'])
    # BURDEN_COLUMNS precede every applied_* column; 'error' is last.
    burden_idx44g = [cols44g.index(c) for c in eb.BURDEN_COLUMNS]
    applied_idx44g = [cols44g.index(c) for c in ('applied_k_1ie', 'applied_k_4ie')]
    assert max(burden_idx44g) < min(applied_idx44g)
    assert cols44g[-1] == 'error'
    assert df44g['trial_number'].tolist() == [0, 1, 2]
    assert df44g['state'].tolist() == ['COMPLETE']*3
    # The burden model's evaluate() ran on every feasible trial (to record
    # columns + prune); apply() is no longer called -- the k_7/k_8 derating
    # moved to the load_simulate choke point.
    assert len(burden44g.seen_evaluate) == 3 and burden44g.seen_apply == []
    for i44g in range(3):
        seen_e = burden44g.seen_evaluate[i44g]
        assert 'inhib_ethanol' not in seen_e         # the group key never reaches the burden model
        assert 'k_10' not in seen_e                  # excluded param: not sampled, not expanded, not passed
        assert {'k_1ie', 'k_4ie'} <= set(seen_e)      # the expanded members ARE passed
        mult44g = df44g['inhib_ethanol'][i44g]        # the sampled multiplier, read from the CSV
        assert np.isclose(seen_e['k_1ie'], 0.02*mult44g, rtol=1e-12, atol=0.0)
        assert np.isclose(seen_e['k_4ie'], 0.04*mult44g, rtol=1e-12, atol=0.0)
        assert np.isclose(df44g['applied_k_1ie'][i44g], seen_e['k_1ie'], rtol=1e-12, atol=0.0)
        assert np.isclose(df44g['applied_k_4ie'][i44g], seen_e['k_4ie'], rtol=1e-12, atol=0.0)
    # The model itself also received the expanded values (via the setattr
    # loop / `applied` -- apply() is gone), the excluded k_10 untouched, the
    # baseline restored afterwards.
    assert len(seen44g) == 4, seen44g   # 3 trials + the finally's restore
    for i44g in range(3):
        spike44g, k1ie44g, k4ie44g, k10_44g = seen44g[i44g]
        assert spike44g == 600.0
        assert np.isclose(k1ie44g, df44g['applied_k_1ie'][i44g], rtol=1e-12, atol=0.0)
        assert np.isclose(k4ie44g, df44g['applied_k_4ie'][i44g], rtol=1e-12, atol=0.0)
        assert k10_44g == 0.01
    assert seen44g[3] == (600.0, 0.02, 0.04, 0.01)      # restore_baseline
    assert te44g.k_1ie == 0.02 and te44g.k_4ie == 0.04 and te44g.k_1e == 47.1
    # Feasibility-aware sampling x groups: the FeasibleTPESampler's
    # predicate shares this burden model, so it too must be handed the
    # EXPANDED members -- a group multiplier reaching burden_model.evaluate()
    # would be scored as if it were a rate constant. A separate fake/study
    # (the run above is pinned with feasible_sampling=False so its
    # evaluate/apply counts stay exactly one pair per trial).
    outdir44f = tempfile.mkdtemp()
    burden44f = _FakeBurden44()
    te44f = _FakeTE44()
    fbs44f = SimpleNamespace(
        current_specifications=dict(target_conc=221.25,
                                    threshold_conc=217.125,
                                    spike_conc=600.0),
        max_n_spikes=16)
    handles44f = dict(handles44, r_te=te44f, fbs_spec=fbs44f,
                      model_specification=lambda **kw: None,
                      latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    buf44f = _io.StringIO()
    with _contextlib.redirect_stdout(buf44f):
        st44f, _, _ = ko.run_kinetic_optimization(enqueue_baseline=True,
            objective='IRR', scenario_label='X', n_trials=3, seed=1,
            study_name='offline_grouped_feasible', results_dir=outdir44f,
            handles=handles44f, print_status_every=1, burden_model=burden44f,
            feasible_sampling=True, n_startup_trials=1, enqueue_knockouts=False,
            exclude_params=('k_10',), parameter_groups=groups44,
            spike_delta_bounds=None)
    out44f = buf44f.getvalue()
    assert type(st44f.sampler).__name__ == 'FeasibleTPESampler', type(st44f.sampler)
    # more evaluate() calls than the 3 objective ones = the sampler's own
    # feasibility checks ran ...
    assert len(burden44f.seen_evaluate) > 3, len(burden44f.seen_evaluate)
    # ... and EVERY dict the burden model saw (predicate or objective) has
    # the members at baseline x multiplier, never the group key.
    for seen44f in burden44f.seen_evaluate:
        assert 'inhib_ethanol' not in seen44f, seen44f
        assert {'k_1ie', 'k_4ie'} <= set(seen44f), seen44f
        assert np.isclose(seen44f['k_4ie']/0.04, seen44f['k_1ie']/0.02,
                          rtol=1e-9)      # one multiplier, both members
        assert 0.2 <= seen44f['k_1ie']/0.02 <= 2.0
    # The group line records each member's LIVE baseline (the basis of its
    # applied_* column), and the pinned-spike line reads the snapshot.
    assert 'Parameter group inhib_ethanol' in out44f, out44f
    assert 'k_1ie (0.02), k_4ie (0.04)' in out44f, out44f
    assert 'Spike feed pinned at the scenario baseline (600 g/L' in out44f, out44f
    PASS('engine: group multiplier applied to every member before each simulation, excluded k_10 untouched, spike pinned at the baseline, applied_* columns after the metrics (sidecar/LOST complete), baseline restored; seeds resolve grouped members through applied_*, the reverse raises; a fake burden_model receives the EXPANDED member values (not the group multiplier) at evaluate()/apply(), never the excluded k_10, with BURDEN_COLUMNS ahead of applied_* in the header, and the feasible-TPE predicate expanding them too; the group print records each member baseline')

#%% 45. metabolic_minimal_subset preset (2026-09-07): a STANDALONE explicit
# set -- 9 listed rate constants + 3 listed inhibition-effector groups
# (0.2x-2x) + 3 feeding variables, spike feed AND stage_1_max_x pinned;
# intersected with the target's workbook (ethanol_only: 5 rates + 2
# groups = 10); typo guard on the role table BEFORE the intersection;
# stage_1_max_x_bounds is a per-type option key surfaced by
# study_type_name_defaults, so the driver's and the supervisor's names
# agree (no _x tag, no _s1x tag); existing presets and names untouched.
assert ko.METABOLIC_MINIMAL_SUBSET_RATES == (
    'k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6', 'k_13', 'k_14', 'k_15', 'k_16')
assert ko.METABOLIC_MINIMAL_SUBSET_GROUPS == {
    'inhib_ethanol': ('k_1ie', 'k_4ie', 'k_7ie', 'k_10ie', 'k_16ie'),
    'inhib_isobutanol': ('k_1ii', 'k_4ii', 'k_6ii', 'k_7ii', 'k_10ii'),
    'inhib_acetate': ('k_1ia', 'k_4ia', 'k_6ia', 'k_7ia', 'k_10ia', 'k_16ia')}
assert list(ko.METABOLIC_MINIMAL_SUBSET_GROUPS) == [
    'inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate']
assert {'METABOLIC_MINIMAL_SUBSET_RATES',
        'METABOLIC_MINIMAL_SUBSET_GROUPS'} <= set(ko.__all__)
assert ko.STUDY_TYPE_ROLES['metabolic_minimal_subset'] == ()       # no role filter: explicit set
opt45 = ko.STUDY_TYPE_OPTIONS['metabolic_minimal_subset']
assert opt45 == dict(rate_params=ko.METABOLIC_MINIMAL_SUBSET_RATES,
                     parameter_groups=ko.METABOLIC_MINIMAL_SUBSET_GROUPS,
                     group_multiplier_bounds=(0.2, 2.0), exclude_params=(),
                     spike_delta_bounds=None, stage_1_max_x_bounds=None)
assert 'group_roles' not in opt45                                   # explicit, not role-grouped
# Name defaults: the subset pins stage_1_max_x (None -> no _s1x tag), has
# no exclusions (no _x tag) and reports the group band as the inhibition
# band; the three older types keep the default (1, 50) g/L band.
assert ko.study_type_name_defaults('metabolic_minimal_subset') == dict(
    inhibition_multiplier_bounds=(0.2, 2.0), exclude_params=(),
    stage_1_max_x_bounds=None)
for st45 in ('metabolic', 'metabolic_protein', 'metabolic_minimal'):
    assert ko.study_type_name_defaults(st45)['stage_1_max_x_bounds'] == (1.0, 50.0), st45
# The default name: the rate band, the group band as _ib, NO exclusion
# tag (empty set) and NO _s1x tag (pinned); driver (ko.default_study_name
# on the preset's effective values) and supervisor (bounds omitted =
# _UNSET -> the type's name default, no longer the hard-coded (1, 50))
# agree; an EXPLICIT --stage-1-max-x-bounds is still tagged; NAME43 and
# the metabolic_protein default are untouched.
NAME45 = 'kin_opt_ethanol_isobutanol_metabolic_minimal_subset_irr_rb0.001-10_ib0.2-2_burden'
assert ko.default_study_name('IRR', 'ethanol_isobutanol', 'metabolic_minimal_subset',
                             rate_multiplier_bounds=(1e-3, 10.0),
                             inhibition_multiplier_bounds=(0.2, 2.0),
                             exclude_params=(), stage_1_max_x_bounds=None,
                             burden=True) == NAME45
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal_subset',
                                   burden=True) == NAME45
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal_subset',
                                   burden=True, stage_1_max_x_bounds=(2.0, 30.0)) \
    == NAME45.replace('_burden', '_s1x2-30_burden')
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal_subset',
                                   burden=True, stage_1_max_x_bounds=None) == NAME45
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_minimal', burden=True) == NAME43
assert sup43['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_rb0.001-10_ib0.1-10_xk10_s1x1-50_burden'
# Source guard: the omitted-flag default reads the name-defaults table,
# not ko.DEFAULT_STAGE_1_MAX_X_BOUNDS.
src45_sup_name = _inspect.getsource(sup43['default_study_name'])
assert "name_defaults['stage_1_max_x_bounds']" in src45_sup_name
assert 'ko.DEFAULT_STAGE_1_MAX_X_BOUNDS' not in src45_sup_name
assert 'metabolic_minimal_subset' in src43_sup
# The 232-character plot-path budget: the longest plot file name of the
# results dir stays under Windows' 260-character limit (long paths are
# disabled on this machine; the mechanical _x tag of an exclusion-based
# definition would have reached 271).
assert len(NAME45) == 81, len(NAME45)
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    p45 = ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal_subset')
    assert p45['stage_1_max_x_bounds'] is None
    for stp45, st45 in (('ethanol_only', 'metabolic'), ('ethanol_only', 'metabolic_protein'),
                        ('ethanol_isobutanol', 'metabolic'),
                        ('ethanol_isobutanol', 'metabolic_protein'),
                        ('ethanol_isobutanol', 'metabolic_minimal')):
        assert ko.resolve_study_preset(stp45, st45)['stage_1_max_x_bounds'] == (1.0, 50.0), (stp45, st45)
else:
    print('SKIP 45 (preset part): parameter-distribution workbooks not found')
if os.path.isfile(wb_A) and os.path.isfile(wb_B):
    roles45 = ko.kinetic_parameter_roles()
    # ethanol_isobutanol: the whole explicit set is in the B workbook.
    p45 = ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal_subset')
    assert p45['scenario'] == 'A' and p45['kinetic_bounds_scenario'] == 'B'
    assert p45['include_params'] == list(ko.METABOLIC_MINIMAL_SUBSET_RATES)
    assert p45['parameter_groups'] == {
        g: list(ms) for g, ms in ko.METABOLIC_MINIMAL_SUBSET_GROUPS.items()}
    assert list(p45['parameter_groups']) == ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate']
    assert p45['exclude_params'] == ()
    assert p45['multiplier_bounds'] == (0.2, 2.0)              # -> the _ib0.2-2 tag
    assert p45['group_multiplier_bounds'] == (0.2, 2.0)
    assert p45['spike_delta_bounds'] is None
    assert p45['stage_1_max_x_bounds'] is None
    assert p45['rate_multiplier_bounds'] == ko.DEFAULT_RATE_MULTIPLIER_BOUNDS
    assert p45['parameter_multiplier_bounds'] == {'k_10': (0.1, 10.0)}
    assert len(p45['rate_params']) == 20                       # the B workbook's capacities, as for every preset
    assert set(p45['include_params']) <= set(p45['rate_params'])
    assert all(roles45[n] == 'capacity' for n in p45['include_params'])
    assert all(roles45[m] in ('product_inhibition', 'lethality')
               for ms in p45['parameter_groups'].values() for m in ms)
    assert set(p45) == set(ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal'))
    # The space it builds (live baselines = the B workbook values here, so
    # every listed rate sits above its floor): 9 rates (rate band, log)
    # + 3 groups (0.2-2, log) + 3 feeding variables = 15; spike_delta and
    # stage_1_max_x absent; the group members are not individual entries;
    # 9 knockout probes (a live ethanol_isobutanol study starts at A with
    # k_13-k_16 clipped to the floor, so it gets 5 -- not tested here).
    kb45 = ko.workbook_kinetic_baselines('B')
    space45, excl45 = ko.build_search_space(
        kb45, include_params=p45['include_params'],
        exclude_params=p45['exclude_params'],
        rate_multiplier_bounds=p45['rate_multiplier_bounds'],
        rate_params=p45['rate_params'],
        parameter_multiplier_bounds=p45['parameter_multiplier_bounds'],
        parameter_groups=p45['parameter_groups'],
        group_multiplier_bounds=p45['group_multiplier_bounds'],
        spike_delta_bounds=p45['spike_delta_bounds'],
        stage_1_max_x_bounds=p45['stage_1_max_x_bounds'])
    assert len(space45) == 15, list(space45)
    assert set(list(space45)[:9]) == set(ko.METABOLIC_MINIMAL_SUBSET_RATES)   # workbook order
    assert list(space45)[9:] == ['inhib_ethanol', 'inhib_isobutanol', 'inhib_acetate',
                                 'threshold_conc', 'target_delta', 'max_n_spikes']
    assert 'spike_delta' not in space45 and 'stage_1_max_x' not in space45
    for r45 in ko.METABOLIC_MINIMAL_SUBSET_RATES:
        assert space45[r45] == dict(low=1e-3*kb45[r45], high=10.0*kb45[r45], log=True), r45
    for g45 in ko.METABOLIC_MINIMAL_SUBSET_GROUPS:
        assert space45[g45] == dict(low=0.2, high=2.0, log=True), g45
    grouped45 = {m for ms in p45['parameter_groups'].values() for m in ms}
    assert set(excl45) == set(kb45) - set(ko.METABOLIC_MINIMAL_SUBSET_RATES) - grouped45
    assert 'k_10' in excl45 and 'k_7' in excl45 and 'k_2' in excl45     # not in the set: at baseline
    bmk45 = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
    pt45 = ko.baseline_decision_point(space45, kb45, bmk45, baseline_max_n_spikes=16,
                                      parameter_groups=p45['parameter_groups'])
    assert all(pt45[g45] == 1.0 for g45 in ko.METABOLIC_MINIMAL_SUBSET_GROUPS)
    probes45, at_floor45 = ko.knockout_probe_points(space45, pt45,
                                                    rate_params=p45['rate_params'])
    assert set(probes45) == set(ko.METABOLIC_MINIMAL_SUBSET_RATES) and at_floor45 == []
    # ethanol_only: intersected with the A workbook -- no k_13-k_16, no
    # isobutanol coefficients: 5 rates + 2 groups + 3 = 10.
    p45_eo = ko.resolve_study_preset('ethanol_only', 'metabolic_minimal_subset')
    assert p45_eo['scenario'] == 'A' and p45_eo['kinetic_bounds_scenario'] == 'A'
    assert p45_eo['include_params'] == ['k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6']
    assert p45_eo['parameter_groups'] == {
        'inhib_ethanol': ['k_1ie', 'k_4ie', 'k_7ie', 'k_10ie'],
        'inhib_acetate': ['k_1ia', 'k_4ia', 'k_6ia', 'k_7ia', 'k_10ia']}
    assert p45_eo['exclude_params'] == () and p45_eo['stage_1_max_x_bounds'] is None
    assert len(p45_eo['rate_params']) == 16
    kb45_eo = ko.workbook_kinetic_baselines('A')
    space45_eo, excl45_eo = ko.build_search_space(
        kb45_eo, include_params=p45_eo['include_params'],
        exclude_params=p45_eo['exclude_params'],
        rate_multiplier_bounds=p45_eo['rate_multiplier_bounds'],
        rate_params=p45_eo['rate_params'],
        parameter_multiplier_bounds=p45_eo['parameter_multiplier_bounds'],
        parameter_groups=p45_eo['parameter_groups'],
        group_multiplier_bounds=p45_eo['group_multiplier_bounds'],
        spike_delta_bounds=p45_eo['spike_delta_bounds'],
        stage_1_max_x_bounds=p45_eo['stage_1_max_x_bounds'])
    assert len(space45_eo) == 10, list(space45_eo)
    assert set(list(space45_eo)[:5]) == {'k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6'}
    assert list(space45_eo)[5:] == ['inhib_ethanol', 'inhib_acetate',
                                    'threshold_conc', 'target_delta', 'max_n_spikes']
    # Typo guard, BEFORE the intersection: a misspelt rate (k_1x is in no
    # workbook, so an intersection-first design would drop it silently), a
    # rate with a non-capacity role, and a group member with a capacity
    # role each raise KeyError naming the parameter and the preset.
    good45 = ko.STUDY_TYPE_OPTIONS['metabolic_minimal_subset']
    bad45_groups = dict(good45['parameter_groups'])
    bad45_groups['inhib_ethanol'] = bad45_groups['inhib_ethanol'] + ('k_3',)
    for bad45, label45 in (
            (dict(good45, rate_params=good45['rate_params'] + ('k_1x',)), 'k_1x'),
            (dict(good45, rate_params=good45['rate_params'] + ('k_1ie',)), 'k_1ie'),
            (dict(good45, parameter_groups=bad45_groups), 'k_3')):
        ko.STUDY_TYPE_OPTIONS['metabolic_minimal_subset'] = bad45
        try:
            ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal_subset')
        except KeyError as e45:
            assert label45 in str(e45) and 'metabolic_minimal_subset' in str(e45), e45
        else:
            raise AssertionError(f'{label45}: the typo guard did not raise')
        finally:
            ko.STUDY_TYPE_OPTIONS['metabolic_minimal_subset'] = good45
    assert ko.STUDY_TYPE_OPTIONS['metabolic_minimal_subset'] is good45
    # The role-filtered presets are untouched (same values as check 43).
    p45_mm = ko.resolve_study_preset('ethanol_isobutanol', 'metabolic_minimal')
    assert len(p45_mm['include_params']) == 36 and p45_mm['exclude_params'] == ('k_10', 'k_7', 'k_8')
    assert p45_mm['stage_1_max_x_bounds'] == (1.0, 50.0)
    # Driver name == supervisor name for the subset on both targets
    # (_driver_name43 mirrors the driver's inhibition-band choice).
    for stp45 in ko.STUDY_TARGET_PRODUCTS:
        assert sup43['default_study_name'](None, 'IRR', None, study_target_products=stp45,
                                           study_type='metabolic_minimal_subset', burden=True) \
            == _driver_name43(stp45, 'metabolic_minimal_subset'), stp45
    assert _driver_name43('ethanol_isobutanol', 'metabolic_minimal_subset') == NAME45
else:
    print('SKIP 45 (space part): parameter-distribution workbooks not found')
# Driver: no logic change (every returned key is already setdefault'ed,
# the print reports both pins); its docstring and runner example name
# the type. Supervisor: help text names it.
drv45 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert "ns['run'](objective='IRR', study_type='metabolic_minimal_subset')" in drv45
assert "'metabolic_minimal_subset'" in drv45 and 'standalone' in drv45.lower()
assert 'metabolic_minimal_subset = ' in src43_sup          # --study-type help text
PASS('metabolic_minimal_subset preset: explicit constants in __all__, empty role entry, options entry with the new stage_1_max_x_bounds key, name defaults (group band, no exclusions, stage_1_max_x pinned; older types keep (1, 50))')

#%% 46. enqueue_baseline flag: engine kwarg gates the trial-0 baseline
# enqueue. Default FALSE since 2026-09-07 (no enqueued baseline for ALL
# studies): a fresh study enqueues NO baseline point so the sampler draws
# every trial; True evaluates the baseline as trial 0. Threaded engine ->
# driver -> supervisor (opt-in --enqueue-baseline), mirroring
# enqueue_knockouts; not part of any study name; resumes unaffected.
_e46 = _inspect.signature(ko.run_kinetic_optimization).parameters
assert 'enqueue_baseline' in _e46 and _e46['enqueue_baseline'].default is False
assert _e46['enqueue_knockouts'].default is False      # both enqueue flags off by default
_src46_eng = _inspect.getsource(ko.run_kinetic_optimization)
assert 'if enqueue_baseline:' in _src46_eng
assert 'baseline_point = baseline_decision_point(' in _src46_eng
assert 'enqueue_baseline' not in _inspect.getsource(ko.default_study_name)
# Driver: run() accepts it (default False) and forwards it to the engine.
drv46 = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'optimize_kinetics_BO.py')).read()
assert 'enqueue_baseline=False,' in drv46
assert 'enqueue_knockouts=False,' in drv46
assert 'enqueue_baseline=enqueue_baseline,' in drv46
# Supervisor: supervise()/child_code() accept it (default False); child_code
# forwards it both ways into the child program; the opt-in
# --enqueue-baseline (store_true) parses to dest enqueue_baseline; main
# forwards args.enqueue_baseline.
sup46 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
assert _inspect.signature(sup46['supervise']).parameters['enqueue_baseline'].default is False
assert _inspect.signature(sup46['child_code']).parameters['enqueue_baseline'].default is False
assert _inspect.signature(sup46['supervise']).parameters['enqueue_knockouts'].default is False
assert _inspect.signature(sup46['child_code']).parameters['enqueue_knockouts'].default is False
code46 = sup46['child_code'](None, 'IRR', 2000, None, False, 'x',
                             study_target_products='ethanol_isobutanol',
                             study_type='metabolic_minimal_subset',
                             enqueue_baseline=True)
assert 'enqueue_baseline=True' in code46
code46b = sup46['child_code'](None, 'IRR', 2000, None, False, 'x',
                              study_target_products='ethanol_isobutanol',
                              study_type='metabolic_minimal_subset')
assert 'enqueue_baseline=False' in code46b and 'enqueue_knockouts=False' in code46b   # the defaults: no baseline, no probes
src46_sup = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'optimize_kinetics_BO_supervised.py')).read()
assert "'--enqueue-baseline'" in src46_sup and "dest='enqueue_baseline'" in src46_sup
# store_true, not store_false: a regression to store_false with the same
# name/dest would silently restore the enqueued-baseline default.
_ebp = src46_sup[src46_sup.index("'--enqueue-baseline'"):][:300]
assert "action='store_true'" in _ebp and "store_false" not in _ebp
assert "'--no-enqueue-baseline'" not in src46_sup    # the old opt-out is gone
assert "'--enqueue-knockouts'" in src46_sup and "'--no-enqueue-knockouts'" not in src46_sup
assert 'enqueue_baseline=args.enqueue_baseline' in src46_sup
assert 'enqueue_baseline=enqueue_baseline' in _inspect.getsource(sup46['supervise'])
# Engine-level, with the real engine on fake handles (check-23 pattern):
# under the DEFAULT a fresh study enqueues NO point at all (probes off too),
# so trial 0 is a sampled draw -- optuna marks an enqueued trial with
# system_attrs['fixed_params']; a sampled one has none. enqueue_baseline=True
# restores the enqueued baseline (k_1e 47.1). Deterministic under the seed.
if _optuna is not None:
    st46 = {'irr': 0.2}
    def _model_specification46(**kw):
        pass
    def _solve_TEA46(stream_IDs=None):
        return {'IRR': st46['irr'], 'MPSPs': {'ethanol': 0.5, 'isobutanol': 1.0}}
    handles46 = dict(handles17, model_specification=_model_specification46,
                     solve_TEA=_solve_TEA46,
                     latest_TEA_solution={'IRR': np.nan,
                                          'MPSPs': {'ethanol': np.nan,
                                                    'isobutanol': np.nan}})
    study46, _, _ = ko.run_kinetic_optimization(
        objective='IRR', scenario_label='X', n_trials=1, seed=1,
        study_name='offline_no_baseline_default',
        results_dir=tempfile.mkdtemp(), handles=handles46,
        print_status_every=1,
        burden_model=None)    # PURE defaults: enqueue_baseline=False, enqueue_knockouts=False
    assert 'fixed_params' not in study46.trials[0].system_attrs   # sampled, not enqueued
    assert not np.isclose(study46.trials[0].params['k_1e'], 47.1)
    study46b, _, _ = ko.run_kinetic_optimization(
        objective='IRR', scenario_label='X', n_trials=1, seed=1,
        study_name='offline_baseline_opt_in',
        results_dir=tempfile.mkdtemp(), handles=handles46,
        enqueue_baseline=True, enqueue_knockouts=False,
        print_status_every=1, burden_model=None)
    assert np.isclose(study46b.trials[0].system_attrs['fixed_params']['k_1e'], 47.1)
    assert np.isclose(study46b.trials[0].params['k_1e'], 47.1)
else:
    print('SKIP 46 (engine part): optuna not installed')
PASS('enqueue_baseline + enqueue_knockouts: BOTH default FALSE for all studies (a default fresh study enqueues NO point: pure-defaults engine run has no fixed_params on trial 0; opt-in restores the enqueued baseline); driver run() and supervisor supervise()/child_code() thread both; opt-in --enqueue-baseline / --enqueue-knockouts, the --no-* opt-outs gone; off study name')

#%% 47. default seed from the study's launch datetime (set 2026-09-10,
# replacing the fixed 3221):
# seed = int((year/day**2)*month*(hour+1)*(minute+1)), computed by the engine
# (ko.default_seed_from_datetime) whenever seed is None -- the new default
# across the engine, the driver run() and the supervisor.
import datetime as _dt47
assert 'default_seed_from_datetime' in ko.__all__
# The exact equation on fixed datetimes.
assert ko.default_seed_from_datetime(_dt47.datetime(2026, 9, 10, 14, 30)) == 84788
assert ko.default_seed_from_datetime(_dt47.datetime(2024, 3, 6, 20, 50)) == 180642
# +1 on hour/minute keeps the hour==0/minute==0 midnight corner positive ...
assert ko.default_seed_from_datetime(_dt47.datetime(2020, 1, 1, 0, 0)) == 2020
# ... and year/day**2 (>= ~2.1 at day 31) keeps even the former zero case >= 1.
assert ko.default_seed_from_datetime(_dt47.datetime(2025, 1, 31, 0, 0)) == 2
_s47 = ko.default_seed_from_datetime()                # when=None -> now()
assert isinstance(_s47, int) and _s47 >= 1
# seed default is the None sentinel in the engine, driver and supervisor.
assert _inspect.signature(ko.run_kinetic_optimization).parameters['seed'].default is None
assert 'seed=None,' in drv30                          # driver run() default
assert _inspect.signature(sup30['supervise']).parameters['seed'].default is None
assert _inspect.signature(sup30['child_code']).parameters['seed'].default is None
assert "'--seed', type=int, default=None" in src30    # supervisor CLI default
PASS('default seed from launch datetime: int((year/day**2)*month*(hour+1)*(minute+1)) '
     'via ko.default_seed_from_datetime (exported), computed by the engine when seed '
     'is None; +1 keeps the midnight corner positive; seed default None across engine '
     '+ driver run() + supervisor supervise()/child_code()/--seed')

#%% 48. LHSDesign: shape, determinism, one-per-stratum (continuous), int column,
# empty design, and the no-log-IntDistribution guard (2026-09-10, LHS start-up).
if _optuna is None:
    print('SKIP 48: optuna not installed')
else:
    space48 = {'a': dict(low=0.01, high=100.0, log=True),
               'b': dict(low=2.0, high=8.0, log=False),
               'n': dict(low=0, high=50, log=False, int=True)}
    d48 = ko.LHSDesign(space48, n_startup=500, seed=7)
    assert d48.size == 500 == len(d48)
    assert 'a' in d48 and 'z' not in d48
    # determinism: same seed -> identical rows; different seed -> different
    d48b = ko.LHSDesign(space48, 500, 7)
    d48c = ko.LHSDesign(space48, 500, 8)
    assert d48.external_point(0) == d48b.external_point(0)
    assert d48.external_point(0) != d48c.external_point(0)
    # every external point carries exactly the space's keys, in-bounds
    for k in (0, 123, 499):
        pt = d48.external_point(k)
        assert set(pt) == {'a', 'b', 'n'}
        assert 0.01 <= pt['a'] <= 100.0 and 2.0 <= pt['b'] <= 8.0
        assert isinstance(pt['n'], int) and 0 <= pt['n'] <= 50
    # one-per-stratum for the LINEAR float column: the 500 internal values of
    # 'b' fall one into each of 500 equal-width bins on [2, 8]
    import numpy as _np48
    b_int = _np48.array([d48.internal_point(k)['b'] for k in range(500)])
    b_bins = _np48.floor((b_int - 2.0) / (8.0 - 2.0) * 500).astype(int)
    b_bins = _np48.clip(b_bins, 0, 499)
    assert sorted(b_bins.tolist()) == list(range(500))
    # one-per-stratum for the LOG float column: log10('a') stratified evenly
    a_int = _np48.array([d48.internal_point(k)['a'] for k in range(500)])
    a_log = _np48.log(a_int)
    a_bins = _np48.floor((a_log - _np48.log(0.01))
                         / (_np48.log(100.0) - _np48.log(0.01)) * 500).astype(int)
    a_bins = _np48.clip(a_bins, 0, 499)
    assert sorted(a_bins.tolist()) == list(range(500))
    # the int column: integer-valued, in range, approximately uniform (floor of
    # a stratified unit-cube column; 500 draws over 51 integers)
    n_vals = _np48.array([d48.internal_point(k)['n'] for k in range(500)])
    assert _np48.all(n_vals == _np48.floor(n_vals)) and n_vals.min() >= 0 and n_vals.max() <= 50
    counts48 = _np48.bincount(n_vals.astype(int), minlength=51)
    assert counts48.min() >= 3 and counts48.max() <= 20   # ~500/51 ≈ 9.8 per value
    # small n_startup <= (high-low+1): the int floor map is monotone/correct
    d48s = ko.LHSDesign({'n': dict(low=0, high=9, log=False, int=True)}, 10, 3)
    small = sorted(d48s.external_point(k)['n'] for k in range(10))
    assert small == list(range(10))
    # empty design: size 0, indexing raises
    d48e = ko.LHSDesign(space48, 0, 1)
    assert d48e.size == 0 and len(d48e) == 0
    try:
        d48e.external_point(0)
    except IndexError:
        pass
    else:
        raise AssertionError('empty LHSDesign did not raise on indexing')
    # a log-scale IntDistribution would be silently linear-floored: guarded
    try:
        ko.LHSDesign({'m': dict(low=1, high=100, log=True, int=True)}, 5, 0)
    except AssertionError:
        pass
    else:
        raise AssertionError('log IntDistribution was not rejected')
    PASS('LHSDesign: shape/determinism, one-per-stratum for linear & log floats, '
         'integer-valued approx-uniform int column, monotone small-n floor, empty '
         'design raises on index, log-IntDistribution guarded')

#%% 49. Trial-count helpers: _n_startup_finished counts ALL COMPLETE|PRUNED
# (enqueued included); _n_sampler_drawn_finished excludes fixed_params trials
# (enqueued baseline/probe/seed), giving the LHS row index k (2026-09-10).
if _optuna is None:
    print('SKIP 49: optuna not installed')
else:
    TS49 = _optuna.trial.TrialState
    st49 = _optuna.create_study()
    # two enqueued (fixed_params) trials + three plain sampler-drawn trials
    st49.enqueue_trial({'x': 0.1})
    st49.enqueue_trial({'x': 0.2})
    st49.optimize(lambda t: t.suggest_float('x', 0.0, 1.0), n_trials=5)
    assert len(st49.trials) == 5
    n_enqueued49 = sum(1 for t in st49.trials if 'fixed_params' in t.system_attrs)
    assert n_enqueued49 == 2
    assert ko._n_startup_finished(st49) == 5
    assert ko._n_sampler_drawn_finished(st49) == 3
    assert (ko._n_startup_finished(st49)
            - ko._n_sampler_drawn_finished(st49)) == n_enqueued49
    # nothing enqueued -> the two counts coincide
    st49b = _optuna.create_study()
    st49b.optimize(lambda t: t.suggest_float('x', 0.0, 1.0), n_trials=4)
    assert ko._n_startup_finished(st49b) == ko._n_sampler_drawn_finished(st49b) == 4
    PASS('trial-count helpers: _n_startup_finished counts enqueued trials, '
         '_n_sampler_drawn_finished excludes fixed_params (row index k); '
         'they differ by n_enqueued and coincide when nothing is enqueued')

#%% 50. FeasibleTPESampler with an LHSDesign: start-up rows come from the design
# when feasible, infeasible rows fall back to draw_uniform_feasible, zero sampled
# INFEASIBLE trials, n_lhs_infeasible_fallbacks matches the infeasible rows, and
# lhs_design=None reproduces today's uniform-feasible start-up (2026-09-10).
if _optuna is None:
    print('SKIP 50: optuna not installed')
else:
    d50 = ko.LHSDesign(space33, n_startup=10, seed=4)
    samp50 = ko.feasible_tpe_sampler(space33, feas33, seed=11,
                                     lhs_design=d50, **tpe_kw33)
    assert type(samp50).__name__ == 'FeasibleTPESampler'
    assert samp50.n_lhs_infeasible_fallbacks == 0
    st50 = _optuna.create_study(direction='maximize', sampler=samp50)
    st50.optimize(_toy33, n_trials=40)          # no enqueued trials: k == trial number
    # zero sampled INFEASIBLE trials and every trial feasible
    assert [t.number for t in st50.trials if t.state == TS33.PRUNED] == []
    assert all(feas33(t.params) for t in st50.trials)
    # the fallback counter equals the number of infeasible design rows over the
    # start-up phase (the first 10 sampler-drawn trials)
    n_infeasible_rows = sum(0 if feas33(d50.external_point(k)) else 1
                            for k in range(10))
    assert samp50.n_lhs_infeasible_fallbacks == n_infeasible_rows
    # a FEASIBLE early row was used verbatim (find the first feasible design row)
    first_feasible = next(k for k in range(10) if feas33(d50.external_point(k)))
    assert st50.trials[first_feasible].params == d50.external_point(first_feasible)
    # lhs_design=None -> byte-for-byte today's uniform-feasible start-up
    import copy as _copy50
    samp50n = ko.feasible_tpe_sampler(space33, feas33, seed=11, **tpe_kw33)
    assert not hasattr(samp50n, '_lhs_design') or samp50n._lhs_design is None
    assert samp50n.n_lhs_infeasible_fallbacks == 0
    st50n = _fresh33(samp50n)
    st50n.optimize(_toy33, n_trials=40)
    samp50n2 = ko.feasible_tpe_sampler(space33, feas33, seed=11, lhs_design=None,
                                       **tpe_kw33)
    st50n2 = _fresh33(samp50n2)
    st50n2.optimize(_toy33, n_trials=40)
    assert ([t.params for t in st50n.trials]
            == [t.params for t in st50n2.trials])
    PASS('FeasibleTPESampler + LHSDesign: LHS start-up rows used when feasible, '
         'infeasible rows fall back to uniform-feasible (counter matches), zero '
         'sampled INFEASIBLE; lhs_design=None reproduces the uniform start-up')

#%% 51. LHSStartupTPESampler (plain path): start-up trials reproduce the design
# rows in order; TPE phase falls through to super().sample_independent; B1
# regression -- with enqueued trials, TPE begins at total COMPLETE|PRUNED ==
# n_startup (identical to plain TPESampler) and only design-row prefix
# 0..n_startup-n_enqueued-1 is used (2026-09-10).
if _optuna is None:
    print('SKIP 51: optuna not installed')
else:
    space51 = {'a': dict(low=0.01, high=100.0, log=True),
               'b': dict(low=2.0, high=8.0, log=False)}
    def _toy51(trial):
        a = trial.suggest_float('a', 0.01, 100.0, log=True)
        b = trial.suggest_float('b', 2.0, 8.0)
        return -((a - 3.0)**2 + (b - 5.0)**2)
    d51 = ko.LHSDesign(space51, n_startup=8, seed=5)
    samp51 = ko.lhs_startup_tpe_sampler(d51, multivariate=True, seed=1,
                                        n_startup_trials=8)
    assert type(samp51).__name__ == 'LHSStartupTPESampler'
    assert isinstance(samp51, _optuna.samplers.TPESampler)
    st51 = _optuna.create_study(direction='maximize', sampler=samp51)
    st51.optimize(_toy51, n_trials=8)           # nothing enqueued: k == trial number
    # the first n_startup sampler-drawn trials reproduce the design rows in order
    for k in range(8):
        assert st51.trials[k].params == d51.external_point(k), k
    # TPE phase: run more trials; they no longer equal design rows (guidance on)
    st51.optimize(_toy51, n_trials=6)
    assert len(st51.trials) == 14
    d51_big = ko.LHSDesign(space51, 8, 5)        # same design, rows 8+ don't exist
    assert st51.trials[8].params != (d51.external_point(7))  # not stuck on last row
    # B1 regression: two enqueued (fixed_params) trials. TPE must still begin at
    # total COMPLETE|PRUNED == n_startup (== 8), using only design rows 0..5.
    d51b = ko.LHSDesign(space51, 8, 5)
    samp51b = ko.lhs_startup_tpe_sampler(d51b, multivariate=True, seed=1,
                                         n_startup_trials=8)
    st51b = _optuna.create_study(direction='maximize', sampler=samp51b)
    st51b.enqueue_trial({'a': 50.0, 'b': 7.0})
    st51b.enqueue_trial({'a': 40.0, 'b': 6.0})
    st51b.optimize(_toy51, n_trials=8)
    # trials 0,1 are the enqueued fixed_params points (bypass the sampler)
    assert st51b.trials[0].params == {'a': 50.0, 'b': 7.0}
    assert st51b.trials[1].params == {'a': 40.0, 'b': 6.0}
    # the 6 sampler-drawn start-up trials (2..7) reproduce design rows 0..5
    for k in range(6):
        assert st51b.trials[2 + k].params == d51b.external_point(k), k
    # one more trial: total finished is now 8 == n_startup, so TPE takes over
    st51b.optimize(_toy51, n_trials=1)
    assert st51b.trials[8].params != d51b.external_point(6)   # NOT design row 6
    PASS('LHSStartupTPESampler: start-up trials reproduce design rows in order, '
         'TPE phase falls through to super(); B1 -- with enqueued trials TPE '
         'begins at total finished == n_startup, only the row prefix is used')

print(f'\nALL {n_pass} CHECKS PASSED')

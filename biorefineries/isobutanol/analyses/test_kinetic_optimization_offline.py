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
    study17_obj, csv17_out, kb17 = ko.run_kinetic_optimization(
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
        ko.run_kinetic_optimization(
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
assert ko.DEFAULT_RATE_MULTIPLIER_BOUNDS == (1e-5, 10.0)
assert ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS == (0.1, 10.0)
space19, excl19 = ko.build_search_space(
    kb19, multiplier_bounds=ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS,
    rate_multiplier_bounds=ko.DEFAULT_RATE_MULTIPLIER_BOUNDS)
assert space19['k_1e'] == dict(low=1e-5*47.1, high=10.0*47.1, log=True)
assert space19['k_7'] == dict(low=1e-5*0.5, high=10.0*0.5, log=True)
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
assert set(ko.STUDY_TYPE_ROLES) == {'metabolic', 'metabolic_protein'}
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
                            'multiplier_bounds', 'rate_multiplier_bounds'}
        assert p21['scenario'] == 'A'                       # both start at the A baseline
        assert p21['kinetic_bounds_scenario'] == ('A' if stp21 == 'ethanol_only' else 'B')
        assert p21['multiplier_bounds'] == (0.1, 10.0)
        assert p21['rate_multiplier_bounds'] == (1e-5, 10.0)
        inc21 = p21['include_params']
        assert len(inc21) == n21, (stp21, st21, len(inc21))
        wb21 = ko.kinetic_param_names_from_scenario(p21['kinetic_bounds_scenario'])
        assert inc21 == [n for n in wb21 if n in inc21]     # workbook order, no extras
        assert all(roles21[n] in ko.STUDY_TYPE_ROLES[st21] for n in inc21)
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
        rate_multiplier_bounds=p22['rate_multiplier_bounds'])
    space22, excl22 = ko.build_search_space(
        model22, multiplier_bounds=p22['multiplier_bounds'],
        rate_multiplier_bounds=p22['rate_multiplier_bounds'],
        param_bounds_override=override22, include_params=p22['include_params'])
    assert excl22 == [] and len(space22) == 56 + 4
    assert all(space22[n]['log'] for n in base22_B)
    for n22, b22 in base22_B.items():
        assert space22[n22]['high'] == 10.0*b22
        assert space22[n22]['low'] == (1e-5*b22 if n22.startswith('k_') else 0.1*b22)
    pt22 = ko.baseline_decision_point(
        space22, model22,
        dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0),
        baseline_max_n_spikes=16)
    for n22 in ehrlich22:
        assert pt22[n22] == 1e-5*base22_B[n22], n22       # clipped to the floor, exactly
    for n22, v22 in model22.items():
        if n22 not in ehrlich22:
            assert pt22[n22] == v22, n22                  # nonzero baselines untouched
    assert pt22['threshold_conc'] == 217.125 and pt22['max_n_spikes'] == 16
    PASS('preset trial 0: Ehrlich rates clipped to exactly 1e-5 x b_B, every other baseline unchanged')
else:
    print('SKIP 22: parameter-distribution workbooks not found')

#%% 23. supervisor: preset flags forwarded to the driver; naming mirrors the engine; legacy switch
sup23 = _runpy.run_path(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    'optimize_kinetics_BO_supervised.py'))
# Naming: preset convention when a target-products preset is in force,
# the legacy convention otherwise (what --legacy-flags selects).
assert sup23['default_study_name'](None, 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'
assert sup23['default_study_name']('A', 'IBO titer', 'B',
                                   study_target_products='ethanol_only',
                                   study_type='metabolic') \
    == 'kin_opt_ethanol_only_metabolic_ibo_titer_kbB'
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
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_scB'
assert sup23['default_study_name'](None, 'IRR', 'A',
                                   'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_kbA'
assert sup23['default_study_name']('A', 'IRR', 'B',
                                   'ethanol_isobutanol', 'metabolic_protein') \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr'  # both match the preset: no tag
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
    study23, _, _ = ko.run_kinetic_optimization(
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
synth24.loc[inf24, 'error'] = 'enzyme burden: Phi_M 0.2806 > F_flex 0.2250 g/gDCW'
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
        study25_obj, csv25_out, kb25 = ko.run_kinetic_optimization(
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
    assert df25['Phi_M'][0] > df25['F_flex'][0] == 0.225
    assert np.isnan(df25['objective'][0]) and np.isnan(df25['IRR'][0])
    assert df25['error'][0].startswith('enzyme burden: Phi_M ')
    assert n_sidecar25 == [1, 2], n_sidecar25
    # trial 1: the CSV keeps the SAMPLED k_7, the model received the EFFECTIVE k_7
    assert np.isclose(df25['k_7'][1], 9.0*1.203)
    assert np.isclose(df25['burden_factor'][1], 0.13276, rtol=1e-3)
    assert np.isclose(df25['k_7_eff'][1], 1.4374, rtol=1e-3)
    assert np.isclose(df25['k_8_eff'][1], 0.13276*0.589, rtol=1e-3)
    assert np.isclose(df25['Phi_M'][1], 0.0637) and np.isclose(df25['phi_T'][1], 9.0*0.135)
    assert df25['objective'][1] == 0.2
    # trial 2: the reference is inert
    assert df25['burden_factor'][2] == 1.0 and df25['k_7_eff'][2] == 1.203
    assert df25['pool_r1'][2] == 0.044 and df25['pool_r13'][2] == 0.0
    # model_specification saw trial 1's effective k_7, trial 2's reference k_7,
    # then restore_baseline's reference k_7 (2 simulations + the finally)
    assert len(seen_k7) == 3, seen_k7
    assert np.isclose(seen_k7[0], 1.4374, rtol=1e-3) and seen_k7[1] == 1.203 and seen_k7[2] == 1.203
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
        ko.run_kinetic_optimization(
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
        ko.run_kinetic_optimization(
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
        ko.run_kinetic_optimization(
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
        study25r, csv25r, kb25r = ko.run_kinetic_optimization(
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
    study25b_obj, csv25b, _ = ko.run_kinetic_optimization(
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
    _, csv25c, _ = ko.run_kinetic_optimization(
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
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_burden'
assert sup26['default_study_name']('B', 'IRR', None,
                                   study_target_products='ethanol_isobutanol',
                                   study_type='metabolic_protein', burden=True) \
    == 'kin_opt_ethanol_isobutanol_metabolic_protein_irr_scB_burden'
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
assert 'eb.BurdenModel.from_reference(' in drv26
assert 'eb.scenario_b_ehrlich()' in drv26 and 'describe_point(' in drv26
assert "'burden_model' in engine_kwargs" in drv26                 # ambiguity guard
PASS('driver/supervisor: burden on by default, --no-burden, _burden naming on both paths, A/B reports wired')

print(f'\nALL {n_pass} CHECKS PASSED')

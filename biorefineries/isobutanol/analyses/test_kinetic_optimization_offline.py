#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline verification of kinetic_optimization's pure logic -- search-space
construction, objective registry, trajectory-CSV round-trip, and (from check 6)
plotting on synthetic data. No isobutanol.load(), no simulation, no optuna.
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
assert space['threshold_conc'] == dict(low=0.0, high=500.0, log=False)
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

print(f'\nALL {n_pass} CHECKS PASSED')

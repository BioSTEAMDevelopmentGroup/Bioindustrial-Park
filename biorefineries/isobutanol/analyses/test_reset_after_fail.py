#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression: a FAILED kinetic-BO trial must leave NO trace in
the flowsheet state the next trial starts from. Motivation (2026-09-14):
in the two 09-14 metabolic_14d GP studies a FAIL followed a FAIL 4x (PI
study: 21 % vs a 5 % base rate) and 3x (titer study: 65 % vs 22 %) more
often than chance, and the two garbage PI rows 1401 / 1768 were both
produced by biosteam's TEA NaN fallback firing on a flowsheet a preceding
SYS14 non-convergence FAIL had left dirty (diverged recycles, the barrage's
converge_method switch left in place, a collapsed S301 split). Mechanism:
`system.snapshot_flowsheet_state()` (every stream's StreamData + the
feed/spike splitter split + every (sub)system's convergence method) /
`system.restore_flowsheet_state(snapshot)` (set_data on every stream, the
split and methods written back, unit caches reset); the BO engine's ONE
evaluation site `ko.evaluate_decision_point` refreshes the context's
snapshot after every COMPLETE / NAN trial (a converged flowsheet) and
restores it after every FAIL trial. One load(), scenario A. Checks: (1) a
round-trip snapshot -> dirty (a different feeding point simulated, recycles
emptied, converge method switched, split collapsed) -> restore reproduces
every stream, the split and the method EXACTLY; (2) re-simulating the
baseline from the restored state reproduces the pinned MPSP and converges
in the same number of sweeps as from the untouched converged baseline (the
restore is harmless); (3) a fresh optimization context carries the current
state as its snapshot; (4) a FAIL trial (the fed-batch volume guard, < 1 s,
with the state dirtied beforehand) restores the snapshot; (5) a COMPLETE
trial (the scenario baseline) reproduces the baseline objective from the
restored state and refreshes the snapshot. The first sweep after a restore
may register drift on the two UNIT-level drift terms (BT801.utility_cost,
V406's spike-feed residual -- unit results, not stream state), so checks 2
and 5 allow at most one extra sweep to the same fixed point. Exit 0 + ALL
CHECKS PASSED = clean. Fresh kernel; 600 s watchdog."""
import os, threading, math, tempfile
import numpy as np
from biorefineries import isobutanol
isobutanol.load()
from biorefineries.isobutanol import scenarios, system
from biorefineries.isobutanol import kinetic_optimization as ko

_watchdog = threading.Timer(600, lambda: (print('WATCHDOG: 600 s exceeded'), os._exit(2)))
_watchdog.daemon = True   # a traceback must end the process, not wait out the timer
_watchdog.start()

n_pass = 0; failures = []
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}', flush=True)
def check(cond, m):
    if cond: PASS(m)
    else: failures.append(m); print(f'FAIL: {m}', flush=True)

spec = scenarios.SCENARIOS['A']
b = scenarios.load_scenario('A')
V406, fbs, r_te, sysm, tea = b['V406'], b['fbs_spec'], b['V406'].nsk_kinetic_model._te, system.corn_EtOH_IBO_sys, system.corn_EtOH_IBO_sys_tea
BASELINE = dict(threshold_conc=spec.threshold_conc, target_conc=spec.target_conc)
PI_get = ko.OBJECTIVE_REGISTRY['PI']['getter']
method_0 = sysm.converge_method

def state_vector():
    """Independent capture of the state the snapshot must reproduce."""
    streams = {}
    for s in sysm.streams:
        streams[s] = (np.array(s.imol.data, copy=True), float(s.T), float(s.P), tuple(s.phases))
    return dict(streams=streams, split=np.array(fbs.splitter.split, copy=True),
                method=sysm.converge_method)

def same_state(a, b):
    if a['method'] != b['method']: return f'method {a["method"]} != {b["method"]}'
    if not np.array_equal(a['split'], b['split']): return 'splitter split differs'
    for s, (mol, T, P, ph) in b['streams'].items():
        mol_a, T_a, P_a, ph_a = a['streams'][s]
        if ph_a != ph or T_a != T or P_a != P or not np.array_equal(mol_a, mol):
            return f'stream {s.ID!r} differs'
    return None

def baseline_solution():
    system.model_specification(**BASELINE)
    res = system.solve_TEA(stream_IDs=('ethanol', 'isobutanol'))
    return dict(mpsp=res['MPSPs']['ethanol'], n_sims=system.last_convergence['n_sims_run'],
                PI=PI_get({'tea': tea}))

# Reference: the converged scenario-A baseline, re-simulated from itself once
ref = baseline_solution()
print(f'baseline re-sim from converged state: MPSP {ref["mpsp"]:.5f}, PI {ref["PI"]:.5f}, {ref["n_sims"]} sweeps', flush=True)
assert abs(ref['mpsp'] - spec.expected['ethanol']) < 1e-2*spec.expected['ethanol']

# (1) round-trip: snapshot -> dirty -> restore == snapshot, exactly
before = state_vector()
snap = system.snapshot_flowsheet_state()
fbs.max_n_spikes = 0
system.model_specification(threshold_conc=34.25, target_conc=140.0)   # a different (batch) feeding point
fbs.max_n_spikes = spec.max_n_spikes
sysm.empty_recycles()
sysm.converge_method = 'fixedpoint' if method_0 != 'fixedpoint' else 'aitken'
fbs.splitter.split = 3e-7                                               # the collapsed-split pathology
dirty = state_vector()
assert same_state(dirty, before) is not None, 'the dirtying did not change the state'
system.restore_flowsheet_state(snap)
diff = same_state(state_vector(), before)
check(diff is None, f'restore reproduces every stream, the S301 split and the converge method exactly ({diff or "no difference"})')

# (2) the restore is harmless: baseline from the restored state == baseline from the converged state
after = baseline_solution()
# (the first sweep after a restore may register drift on the two unit-level
# drift terms -- BT801.utility_cost and V406's spike-feed residual -- which
# are unit results, not stream state, and are re-established by that sweep;
# so at most ONE extra sweep, and the same fixed point)
check(abs(after['mpsp'] - ref['mpsp']) < 1e-6*ref['mpsp'] and after['n_sims'] <= ref['n_sims'] + 1,
      f'baseline from the restored state reproduces MPSP {after["mpsp"]:.6f} (ref {ref["mpsp"]:.6f}) in {after["n_sims"]} sweeps (ref {ref["n_sims"]}, at most +1)')

# (3) the engine: a fresh context snapshots the current (converged) state
tmp = tempfile.mkdtemp(prefix='kin_opt_reset_test_')
ctx = ko._prepare_optimization(
    'PI', direction=None, level=None, objective_units=None, objective_name=None,
    scenario_label='A', multiplier_bounds=(0.1, 10.0), param_bounds_override=None,
    exclude_params=(), include_params=('k_3',), rate_multiplier_bounds=None,
    rate_params=None, parameter_multiplier_bounds=None,
    threshold_conc_bounds=(0.0, 300.0), target_delta_bounds=(5.0, 500.0),
    spike_delta_bounds=ko.DEFAULT_SPIKE_DELTA_BOUNDS, max_n_spikes_bounds=(0, 50),
    stage_1_max_x_bounds=None, target_conc_bounds=None, threshold_delta_bounds=None,
    spike_conc_bounds=None, study_name='test_reset_after_fail_tmp', results_dir=tmp,
    handles=None, burden_model=None, volume_feasibility=False, volume_cap=None,
    seed_from=None, parameter_groups=None,
    group_multiplier_bounds=ko.DEFAULT_GROUP_MULTIPLIER_BOUNDS)
converged = state_vector()
system.restore_flowsheet_state(ctx.state_snapshot)   # a no-op iff the snapshot IS the current state
diff = same_state(state_vector(), converged)
check(ctx.state_snapshot is not None and diff is None,
      f'a fresh optimization context snapshots the converged state it was built on ({diff or "no difference"})')

# (4) a FAIL trial restores the snapshot (state dirtied beforehand)
k3 = ctx.kinetic_baselines['k_3']
sysm.empty_recycles(); sysm.converge_method = 'fixedpoint' if method_0 != 'fixedpoint' else 'aitken'; fbs.splitter.split = 3e-7
# test_fed_batch_volume_guard's blow-up point: spike 10 g/L above the target -> x21 per spike, the
# runtime guard raises FeedingStrategyError after the kinetic run and BEFORE the flowsheet simulate
blowup = dict(k_3=k3, threshold_conc=100.0, target_delta=200.0, spike_delta=10.0, max_n_spikes=16)
ev = ko.evaluate_decision_point(ctx, blowup, 1)
diff = same_state(state_vector(), converged)
check(ev.state == 'FAIL' and 'max_fed_batch_volume_ratio' in str(ev.exception) and diff is None,
      f'a FAIL trial ({ev.state}: {str(ev.exception)[:60]!r}) restores the converged snapshot ({diff or "no difference"})')

# (5) a COMPLETE trial from the restored state reproduces the baseline and refreshes the snapshot
base = dict(k_3=k3, threshold_conc=spec.threshold_conc,
            target_delta=spec.target_conc - spec.threshold_conc,
            spike_delta=600.0 - spec.target_conc, max_n_spikes=spec.max_n_spikes)
old_snapshot = ctx.state_snapshot
ev = ko.evaluate_decision_point(ctx, base, 2)
completed = state_vector()
system.restore_flowsheet_state(ctx.state_snapshot)
diff = same_state(state_vector(), completed)
check(ev.state == 'COMPLETE' and abs(ev.objective - ref['PI']) < 1e-5*abs(ref['PI'])
      and ev.record['n_sims_run'] <= ref['n_sims'] + 1 and ctx.state_snapshot is not old_snapshot and diff is None,
      f'the baseline trial after the FAIL is COMPLETE with PI {ev.objective if ev.objective is not None else float("nan"):.5f} '
      f'(ref {ref["PI"]:.5f}) in {ev.record.get("n_sims_run")} sweeps (ref {ref["n_sims"]}, at most +1) and refreshes the snapshot')

if failures:
    print(f'\n{len(failures)} CHECK(S) FAILED, {n_pass} passed'); os._exit(1)
print(f'\nALL {n_pass} CHECKS PASSED'); os._exit(0)

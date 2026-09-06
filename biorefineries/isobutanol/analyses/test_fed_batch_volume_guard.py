#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression test of the fed-batch working-volume guard in
system.load_simulate (added 2026-09-06 after the kinetic-BO stall debugging).

Background: nskinetics' FeedSpike adds env*(target - s)/(spike - target) of
volume per glucose spike, so the working volume multiplies by
(spike_conc - threshold)/(spike_conc - target) per spike with no cap. A
feeding point with the spike concentration barely above the target compounds
to a 1e6-1e12x volume growth; NSKBatchReactor scales its effluent by that
ratio on the first flowsheet pass, the WWT membrane bioreactor's tank-count
loop (O(F_vol)) then runs for hours ("stall"), and the feeding spec leaves
S301.split ~ 1e-7, which poisons the NEXT simulation (the study's "could not
reach desired concentration ... F301" actuator error). The guard raises a
FeedingStrategyError BEFORE the flowsheet simulate when the batch's
final/initial volume ratio exceeds parameters['max_fed_batch_volume_ratio']
(20 = initial charge >= 5 % of the final working volume), restoring the
splitter so the process state stays clean; model_specification re-raises it
at once instead of running the convergence-recovery barrage.

Sequential #%% checks; each prints 'PASS n:'; exit 0 + 'ALL n CHECKS PASSED'
= clean. Runs in a fresh kernel (one load()).
"""
import os, sys, time, math, threading
import numpy as np

HANG_TIMEOUT_S = 180   # a guard-less run floods the flowsheet and hangs in WWT


def _hard_exit_on_hang(label):
    def _die():
        sys.stdout.flush()
        sys.stderr.write(f'\nFAIL: {label} did not return within {HANG_TIMEOUT_S} s '
                         '(the fed-batch volume blow-up hang)\n')
        sys.stderr.flush()
        os._exit(1)
    t = threading.Timer(HANG_TIMEOUT_S, _die)
    t.daemon = True
    t.start()
    return t


from biorefineries import isobutanol
isobutanol.load()
from biorefineries.isobutanol import system as ibo_system
from nskinetics.exceptions import FeedingStrategyError

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
model_specification = model.specification
solve_TEA = ibo_system.solve_TEA
parameters = ibo_system.parameters
f = model.system.flowsheet
V406, S301 = f.V406, f.S301
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')

# Scenario-A baseline (smoke-test / BO-driver recipe)
model.parameters = ()
model.load_parameter_distributions(
    IBO_filepath + '\\analyses\\full\\parameter_distributions\\'
    'parameter-distributions_corn_IBO_EtOH_A.xlsx', namespace_dict)
model.metrics_at_baseline()
model_specification()
fbs_spec.max_n_spikes = 16
A_KWARGS = dict(threshold_conc=217.125, target_conc=221.25)
A_ETHANOL_MPSP = 0.86604
model_specification(**A_KWARGS)
baseline_TEA = solve_TEA()
assert abs(baseline_TEA['MPSPs']['ethanol'] - A_ETHANOL_MPSP)/A_ETHANOL_MPSP < 0.01, baseline_TEA
print('baseline ready:', baseline_TEA)

cv = fbs_spec.control_variables
vol_col = cv.resolve_volume_col(V406)
add_col = cv.resolve_feed_volume_added_col(V406)

def volume_ratio():
    d = V406.nsk_results_specific_tau_dict
    env, added = d[vol_col], d[add_col]
    return env/(env - added)

def split_value():
    return float(np.asarray(S301.split).ravel()[0])

# A blow-up feeding point: spike only 10 g/L above the target, so every spike
# adds (300 - 100)/(310 - 300) = 20 working volumes (x21 per spike).
BLOWUP = dict(threshold_conc=100.0, target_conc=300.0, spike_conc=310.0)
n_checks = 0

#%% 1. The cap is a named process parameter: initial charge >= 5 % of the final volume
assert parameters['max_fed_batch_volume_ratio'] == 20.0, parameters.get('max_fed_batch_volume_ratio')
n_checks += 1; print(f'PASS {n_checks}: max_fed_batch_volume_ratio = 20 (initial charge >= 5 %)')

#%% 2. Baseline fed batch is far below the cap (guard inert at baseline)
r_base = volume_ratio()
assert 1.0 <= r_base < 2.0, r_base
n_checks += 1; print(f'PASS {n_checks}: scenario-A baseline volume ratio {r_base:.4f} < cap')

#%% 3. load_simulate raises FeedingStrategyError on a blow-up point BEFORE the flowsheet runs
split_before = split_value()
n_log_before = len(ibo_system.convergence_log)
timer = _hard_exit_on_hang('load_simulate(blow-up point)')
t0 = time.time()
try:
    ibo_system.load_simulate(**BLOWUP)
except FeedingStrategyError as e:
    err = e
else:
    raise AssertionError('load_simulate accepted a fed batch whose volume grows > 20x')
finally:
    timer.cancel()
elapsed = time.time() - t0
msg = str(err)
assert 'working volume' in msg.lower() and 'max_fed_batch_volume_ratio' in msg, msg
assert len(ibo_system.convergence_log) == n_log_before, 'the flowsheet was simulated before the guard fired'
assert elapsed < 60, elapsed
n_checks += 1; print(f'PASS {n_checks}: guard raised in {elapsed:.1f} s: {msg[:160]}')

#%% 4. The splitter is restored on violation (no stale ~0 split poisoning the next run)
assert math.isclose(split_value(), split_before, rel_tol=1e-9), (split_value(), split_before)
n_checks += 1; print(f'PASS {n_checks}: S301.split restored to {split_before:.5f}')

#%% 5. model_specification re-raises immediately (no convergence-recovery barrage)
n_log_before = len(ibo_system.convergence_log)
timer = _hard_exit_on_hang('model_specification(blow-up point)')
try:
    model_specification(**BLOWUP)
except FeedingStrategyError:
    pass
else:
    raise AssertionError('model_specification accepted the blow-up point')
finally:
    timer.cancel()
assert len(ibo_system.convergence_log) == n_log_before, \
    'model_specification ran the recovery barrage (reset_and_reload) on a volume-guard error'
n_checks += 1; print(f'PASS {n_checks}: model_specification re-raised without the barrage')

#%% 6. The process state is clean: the baseline reproduces its pinned MPSP afterwards
model_specification(**A_KWARGS)
after = solve_TEA()
assert abs(after['MPSPs']['ethanol'] - A_ETHANOL_MPSP)/A_ETHANOL_MPSP < 0.01, after
assert abs(after['MPSPs']['ethanol'] - baseline_TEA['MPSPs']['ethanol'])/A_ETHANOL_MPSP < 5e-3, (after, baseline_TEA)
n_checks += 1; print(f'PASS {n_checks}: baseline after the guarded failures: {after}')

#%% 7. A batch (max_n_spikes = 0) has ratio exactly 1 and passes the guard
ibo_system.load_simulate(threshold_conc=34.25, target_conc=140.0, spike_conc=600.0, max_n_spikes=0)
assert volume_ratio() == 1.0, volume_ratio()
n_checks += 1; print(f'PASS {n_checks}: batch (no spikes) volume ratio 1.0 passes')
fbs_spec.max_n_spikes = 16
model_specification(**A_KWARGS)

print(f'ALL {n_checks} CHECKS PASSED')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression: the enzyme burden must derate V406's kinetics on
EVERY path that runs them, not only inside load_simulate. Motivation
(2026-09-14): biosteam's TEA re-simulates the system BARE on a NaN cashflow
(_tea.py, `self.system.simulate()` before raising 'nan encountered in
cashflow array'), outside system._apply_enzyme_burden(); a burden-ON point
then re-ran V406 at the INTENDED k_7 and the kinetic-BO PI study logged the
un-derated fermentation with a garbage TEA (trials 1401 / 1768: PI -4.8e6 /
-1.0e8, TCI 14 / 6.9 MM$; true values PI -3.7, TCI ~105). One load(),
scenario A, then PI-study trial 1768's point (a stalled culture whose
derated run is tau = tau_max with 0 spikes, and whose UN-derated run is a
17.5 h, 39-spike, 108 g/L-cells batch -- an unmistakable fingerprint).
Checks: (1) _apply_enzyme_burden nests without derating twice; (2) the
engine path (model_specification) gives the derated fingerprint; (3) a bare
V406.simulate() and (4) a bare corn_EtOH_IBO_sys.simulate() leave that
fingerprint unchanged and restore the intended k_7; (5) with no active
burden the same bare V406.simulate() runs UN-derated (the guard defers to
the active burden); (6) re-activating the burden derates again.
Exit 0 + ALL CHECKS PASSED = clean. Fresh kernel; 600 s watchdog."""
import os, threading, math
from biorefineries import isobutanol
isobutanol.load()
from biorefineries.isobutanol import scenarios, system
from biorefineries.isobutanol import kinetic_optimization as ko

threading.Timer(600, lambda: (print('WATCHDOG: 600 s exceeded'), os._exit(2))).start()

n_pass = 0; failures = []
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}', flush=True)
def check(cond, m):
    if cond: PASS(m)
    else: failures.append(m); print(f'FAIL: {m}', flush=True)

b = scenarios.load_scenario('A')
V406, fbs, r, sysm = b['V406'], b['fbs_spec'], b['V406'].nsk_kinetic_model._te, system.corn_EtOH_IBO_sys
burden = system.get_active_burden(); assert burden is not None

# PI-study trial 1768 (kin_opt_ethanol_isobutanol_metabolic_14d_pi_gp_..._20260914):
# individual rates + the glycolysis-group / inhibition-family members as applied.
POINT = dict(k_3=0.015510430089350047, k_6=0.13226846931687894,
             k_13=22.297579387234045, k_14=3.202785013306183,
             k_15=0.05127249358413282, k_16=0.008436450173016873,
             k_1l=1.8139941226168712, k_1h=0.7408199773484285, k_1e=59.747638584094155,
             k_1ie=0.017610668881296303, k_4ie=0.035221337762592606,
             k_7ie=0.035221337762592606, k_10ie=0.035221337762592606,
             k_16ie=0.017610668881296303,   # inert since the 2026-09-15 r16/r17 split
             k_17ie=0.017610668881296303,   # the live r17 twin (same inhib_ethanol multiplier)
             k_1ii=0.07124378993369837, k_4ii=0.14248757986739674,
             k_6ii=0.07124378993369837, k_7ii=0.14248757986739674,
             k_10ii=0.14248757986739674,
             k_1ia=0.015607831370166643, k_4ia=0.031215662740333286,
             k_6ia=0.015607831370166643, k_7ia=0.031215662740333286,
             k_10ia=0.031215662740333286,
             k_16ia=0.015607831370166643,   # inert since the 2026-09-15 r16/r17 split
             k_17ia=0.015607831370166643)   # the live r17 twin (same inhib_acetate multiplier)
for k, v in POINT.items(): setattr(r, k, v)
fbs.max_n_spikes = 39
V406.stage_1_max_x = 1.9600132013337246
FEED = dict(threshold_conc=279.7508210875094,
            target_conc=min(ko.TARGET_CONC_MAX, 279.7508210875094 + 37.42941116448492))
k7_intended = float(r.k_7)

def fingerprint():
    d = V406.nsk_results_specific_tau_dict
    return dict(tau=float(V406.tau), spikes=float(d['curr_n_glu_spikes']), x=float(d['[x]']))
def same(a, b):
    return (abs(a['tau'] - b['tau']) <= 1e-6*max(1.0, abs(b['tau'])) and a['spikes'] == b['spikes']
            and abs(a['x'] - b['x']) <= 1e-6*max(1e-9, abs(b['x'])))

# (1) the burden context must be re-entrant: nesting derates ONCE, never twice
with system._apply_enzyme_burden():
    k7_outer = float(r.k_7)
    with system._apply_enzyme_burden():
        k7_inner = float(r.k_7)
check(k7_outer < k7_intended and math.isclose(k7_inner, k7_outer, rel_tol=1e-12) and math.isclose(float(r.k_7), k7_intended, rel_tol=1e-12),
      f'_apply_enzyme_burden nests without a double derate (outer {k7_outer:.5g} == inner {k7_inner:.5g} < intended {k7_intended:.5g}; restored after)')

# (2) the engine path: derated run = stalled culture (tau_max, no spikes)
b['model_specification'](**FEED)
ref = fingerprint()
check(ref['tau'] >= 100.0 and ref['spikes'] == 0 and ref['x'] < 0.1,
      f'model_specification gives the DERATED fingerprint (tau {ref["tau"]:.4g} h, {ref["spikes"]:g} spikes, x {ref["x"]:.3g} g/L)')

# (3) a bare V406.simulate() must give the same derated result and leave k_7 intended
V406.simulate(); fp = fingerprint()
check(same(fp, ref) and math.isclose(float(r.k_7), k7_intended, rel_tol=1e-12),
      f'bare V406.simulate() stays derated (tau {fp["tau"]:.4g} h, {fp["spikes"]:g} spikes, x {fp["x"]:.3g} g/L; k_7 restored)')

# (4) a bare system simulate (what biosteam's TEA NaN fallback runs) likewise
sysm.simulate(); fp = fingerprint()
check(same(fp, ref) and math.isclose(float(r.k_7), k7_intended, rel_tol=1e-12),
      f'bare corn_EtOH_IBO_sys.simulate() stays derated (tau {fp["tau"]:.4g} h, {fp["spikes"]:g} spikes, x {fp["x"]:.3g} g/L)')

# (5) control: with NO active burden the same bare call runs UN-derated
system.set_active_burden(None)
try:
    V406.simulate(); fp = fingerprint()
    check(fp['spikes'] == 39 and 15.0 < fp['tau'] < 20.0 and fp['x'] > 50.0,
          f'no active burden -> bare V406.simulate() runs UN-derated (tau {fp["tau"]:.4g} h, {fp["spikes"]:g} spikes, x {fp["x"]:.3g} g/L)')
finally:
    system.set_active_burden(burden)

# (6) re-activated burden derates again
V406.simulate(); fp = fingerprint()
check(same(fp, ref), f're-activated burden derates again (tau {fp["tau"]:.4g} h, {fp["spikes"]:g} spikes)')

if failures:
    print(f'\n{len(failures)} CHECK(S) FAILED, {n_pass} passed'); os._exit(1)
print(f'\nALL {n_pass} CHECKS PASSED'); os._exit(0)

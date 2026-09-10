#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based regression of scenario-level enzyme-burden enforcement
(scenarios.load_scenario + the system.py choke point). One load(),
scenario A. Verifies: (1) A is burden-ON by default yet inert (baseline
MPSP unchanged, active burden installed); (2) an over-cap kinetic point
raises EnzymeBurdenInfeasibleError through model_specification while a
feasible neighbour simulates; (3) burden=False installs no active burden.
Exit 0 + ALL CHECKS PASSED = clean. Fresh kernel."""
import math
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios, system
from biorefineries.isobutanol import enzyme_burden as eb

n_pass = 0
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}')

# (1) A default -> burden ON but inert
b = scenarios.load_scenario('A')
assert b['burden_on'] is True and b['burden_model'] is not None
assert system.get_active_burden() is b['burden_model']
res = b['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'), IRR_for_MPSP=0.15)
assert abs(res['MPSPs']['ethanol'] - 0.86604)/0.86604 < 0.01
PASS('scenario A burden ON by default, active, and inert (ethanol MPSP 0.86604)')

# (2) an over-cap point -> EnzymeBurdenInfeasibleError; a feasible one simulates
r = b['V406'].nsk_kinetic_model._te
base_k1e = float(getattr(r, 'k_1e'))
try:
    setattr(r, 'k_1e', base_k1e*1000.0)      # blow up the r1 native pool -> over cap
    raised = False
    try:
        b['model_specification'](**b['feeding_kwargs'])
    except system.EnzymeBurdenInfeasibleError:
        raised = True
    assert raised, 'over-cap point did not raise EnzymeBurdenInfeasibleError'
    PASS('over-cap kinetic point raises EnzymeBurdenInfeasibleError via model_specification')
finally:
    setattr(r, 'k_1e', base_k1e)
b['model_specification'](**b['feeding_kwargs'])   # feasible again
res2 = b['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'), IRR_for_MPSP=0.15)
assert abs(res2['MPSPs']['ethanol'] - 0.86604)/0.86604 < 0.01
PASS('feasible neighbour simulates cleanly after the infeasible point')

# (3) burden=False installs no active burden
b3 = scenarios.load_scenario('A', burden=False)
assert b3['burden_on'] is False and b3['burden_model'] is None
assert system.get_active_burden() is None
PASS('burden=False installs no active burden')

print(f'\nALL {n_pass} CHECKS PASSED')

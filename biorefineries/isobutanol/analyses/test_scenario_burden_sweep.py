#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Sim-based check that a kinetic sweep inherits the enzyme burden through
the load_simulate choke point exactly as the evaluate_* sweeps do (set a
capacity on r, call model_specification inside try/except -> NaN). One
load(), scenario A. Verifies: an over-cap point -> NaN (not a crash); a
feasible point -> a finite MPSP; and burden-ON == burden-OFF for the
scenario-A baseline (A is inert). Exit 0 + ALL CHECKS PASSED = clean."""
import math
import numpy as np
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios, system

n_pass = 0
def PASS(m):
    global n_pass; n_pass += 1; print(f'PASS {n_pass}: {m}')

def swept_MPSP(bundle, r, k1e_value):
    """The evaluate_* pattern: set a capacity, simulate, read MPSP, NaN on
    any raise (incl. EnzymeBurdenInfeasibleError)."""
    base = float(getattr(r, 'k_1e'))
    try:
        setattr(r, 'k_1e', k1e_value)
        bundle['model_specification'](**bundle['feeding_kwargs'])
        res = bundle['solve_TEA'](stream_IDs=('ethanol', 'isobutanol'),
                                  IRR_for_MPSP=0.15)
        return res['MPSPs']['ethanol']
    except Exception:
        return np.nan
    finally:
        setattr(r, 'k_1e', base)

# Burden ON (A default): over-cap point -> NaN; feasible baseline -> finite
b = scenarios.load_scenario('A')
r = b['V406'].nsk_kinetic_model._te
base_k1e = float(getattr(r, 'k_1e'))
over = swept_MPSP(b, r, base_k1e*1000.0)
assert math.isnan(over), f'over-cap point returned {over}, expected NaN'
PASS('burden ON: an over-cap swept k_1e yields NaN (not a crash)')
b['model_specification'](**b['feeding_kwargs'])              # restore baseline sim
feasible_on = swept_MPSP(b, r, base_k1e)                     # baseline value
assert math.isfinite(feasible_on) and abs(feasible_on - 0.86604)/0.86604 < 0.01
PASS('burden ON: the feasible baseline k_1e yields a finite MPSP (0.86604)')

# Burden OFF: same feasible baseline point is identical (A is inert either way)
b_off = scenarios.load_scenario('A', burden=False)
r_off = b_off['V406'].nsk_kinetic_model._te
feasible_off = swept_MPSP(b_off, r_off, float(getattr(r_off, 'k_1e')))
assert math.isfinite(feasible_off)
assert abs(feasible_on - feasible_off)/abs(feasible_off) < 5e-3
PASS('burden ON == burden OFF for the scenario-A baseline (inert)')

print(f'\nALL {n_pass} CHECKS PASSED')

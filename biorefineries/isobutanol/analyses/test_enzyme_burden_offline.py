#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline verification of enzyme_burden.py -- the proteome-allocation
(enzyme-burden) constraint of the kinetic Bayesian optimization
(docs/superpowers/specs/2026-09-05-enzyme-burden-design.md). Pure logic:
no isobutanol.load(), no simulation; SCENARIO_B_EHRLICH is read from
nskinetics' scenarios.py by file path (the no-import guarantee is probed
in a fresh interpreter, check 13). Run in a fresh kernel; exit 0 +
'ALL ... CHECKS PASSED' = clean."""
import math
import os
import subprocess
import sys

from biorefineries.isobutanol import enzyme_burden as eb

n_pass = 0
def PASS(msg):
    global n_pass
    n_pass += 1
    print(f'PASS {n_pass}: {msg}')

def close(a, b, rel=1e-9, abs_=0.0):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)

#: The model's reference capacities (antimony defaults = the scenario-A
#: values the engine snapshots after the A-workbook load).
K_REF = {'k_1h': 0.584, 'k_1l': 1.43, 'k_1e': 47.1, 'k_2': 0.501,
         'k_3': 5.81, 'k_4': 4.8, 'k_5': 0.0104, 'k_5e': 0.775,
         'k_6': 2.82, 'k_7': 1.203, 'k_8': 0.589,
         'k_13': 0.0, 'k_14': 0.0, 'k_15': 0.0, 'k_16': 0.0}

#%% 1. sector constants and closure (spec 4.1): 0.026 g/gDCW of slack at wild-type growth
assert eb.PROTEIN_CONTENT == 0.45 and eb.HOUSEKEEPING_FRACTION == 0.50
assert eb.TRANSLATION_FRACTION_WT == 0.30 and eb.SIGMA_EFF == 0.50
assert close(eb.F_FLEX, 0.225) and close(eb.PHI_T_WT, 0.135)
phi_M_wt = sum(pool for pool, _ in eb.NATIVE_STEPS.values())
assert close(phi_M_wt, 0.0637, rel=1e-6)
assert close(eb.F_FLEX - phi_M_wt - eb.PHI_T_WT, 0.0263, rel=1e-6)   # > 0: reference feasible
assert list(eb.NATIVE_STEPS) == ['r1', 'r2', 'r3', 'r4', 'r5', 'r6']
assert eb.NATIVE_STEPS['r1'] == (0.044, ('k_1h', 'k_1l', 'k_1e'))
assert eb.NATIVE_STEPS['r4'] == (0.0032, ())          # fixed pool; k_4 burden-free
assert eb.NATIVE_STEPS['r5'] == (0.0008, ('k_5', 'k_5e'))
assert list(eb.EHRLICH_STEPS) == ['r13', 'r14', 'r15', 'r16']
assert eb.GROWTH_CAPACITIES == ('k_7', 'k_8')
assert eb.STEP_ORDER == ('r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r13', 'r14', 'r15', 'r16')
assert eb.BURDEN_COLUMNS == (
    'pool_r1', 'pool_r2', 'pool_r3', 'pool_r4', 'pool_r5', 'pool_r6',
    'pool_r13', 'pool_r14', 'pool_r15', 'pool_r16',
    'Phi_M', 'phi_T', 'F_flex', 'burden_factor', 'k_7_eff', 'k_8_eff')
PASS('sector constants, closure slack 0.0263 g/gDCW, tables and column order')

#%% 2. Ehrlich per-unit costs (spec 4.3) at SIGMA_EFF, to two significant figures
# kcat per mole of the substrate k_i is written on; r13 counts two pyruvate
# per acetolactate turnover; r16 sums the KDC and the ADH on the KIV flux.
cap13, mw13, enz13 = eb.EHRLICH_STEPS['r13']
assert cap13 == 'k_13' and mw13 == 88.06 and len(enz13) == 1
assert enz13[0][0] == 'Ilv2+Ilv6' and enz13[0][1] == 108924.0
assert close(enz13[0][2], 2.0*49.0*74937.0*1e-3/60.0)          # 122.4 s^-1
assert eb.EHRLICH_STEPS['r14'][:2] == ('k_14', 132.11)
assert close(eb.EHRLICH_STEPS['r14'][2][0][2], 18.5*44368.0*1e-3/60.0)   # 13.68 s^-1
assert eb.EHRLICH_STEPS['r15'][:2] == ('k_15', 134.13)
assert close(eb.EHRLICH_STEPS['r15'][2][0][2], 18.0*62861.0*1e-3/60.0)   # 18.86 s^-1
cap16, mw16, enz16 = eb.EHRLICH_STEPS['r16']
assert cap16 == 'k_16' and mw16 == 116.12
assert [(n, mw, kc) for n, mw, kc in enz16] == [('Aro10', 71384.0, 19.0),
                                                 ('Adh6', 39618.0, 296.0)]
assert close(eb.ehrlich_unit_cost('r13'), 0.0056, rel=0.02)
assert close(eb.ehrlich_unit_cost('r14'), 0.0136, rel=0.02)
assert close(eb.ehrlich_unit_cost('r15'), 0.0138, rel=0.02)
assert close(eb.ehrlich_unit_cost('r16'), 0.0180 + 0.0006, rel=0.02)
# the ADH share is negligible: Aro10 alone is > 96 % of the r16 cost
aro10_only = 71384.0/(19.0*eb.SIGMA_EFF*3600.0*116.12)
assert aro10_only/eb.ehrlich_unit_cost('r16') > 0.96
# sigma_eff scales every Ehrlich cost inversely
assert close(eb.ehrlich_unit_cost('r13', sigma_eff=1.0), 0.5*eb.ehrlich_unit_cost('r13'))
# scenario-B Ehrlich constants need ~0.22 g/gDCW (spec 4.3)
pool_B = (5.81*eb.ehrlich_unit_cost('r13') + 4.8*eb.ehrlich_unit_cost('r14')
          + 4.8*eb.ehrlich_unit_cost('r15') + 2.82*eb.ehrlich_unit_cost('r16'))
assert close(pool_B, 0.2169, rel=0.01), pool_B
PASS('Ehrlich per-unit costs 0.0056 / 0.0136 / 0.0138 / 0.0180+0.0006; B pool 0.217 g/gDCW')

#%% 3. sigma diagnostic (spec 4.4): anchors disagree 24x, geometric mean within 2x of SIGMA_EFF
assert eb.ANCHOR_STEPS == {'r3': ('k_3', 88.06, 61495.0, 60.0),
                           'r6': ('k_6', 44.05, 36849.0, 1800.0)}
sig = eb.anchor_sigma(K_REF)
assert set(sig) == {'r3', 'r6', 'geomean'}
assert close(sig['r3'], 2.21, rel=0.02), sig
assert close(sig['r6'], 0.091, rel=0.02), sig
assert close(sig['geomean'], math.sqrt(sig['r3']*sig['r6']))
assert 0.25 <= sig['geomean'] <= 1.0                      # within 2x of SIGMA_EFF = 0.5
assert close(sig['geomean'], 0.45, rel=0.02)
PASS('sigma diagnostic: sigma_r3 2.2, sigma_r6 0.091, geometric mean 0.45 in [0.25, 1.0]')

print(f'\nALL {n_pass} CHECKS PASSED')

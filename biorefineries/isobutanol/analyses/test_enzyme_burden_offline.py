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

#%% 1. sector constants and closure (spec 4.1, amended 2026-09-06): 0.029 g/gDCW of slack at wild-type growth
assert eb.PROTEIN_CONTENT == 0.49 and eb.HOUSEKEEPING_FRACTION == 0.50   # Jach et al. 2022, Table 1
assert eb.POOL_TABLE_PROTEIN_CONTENT == 0.45                            # the pool table's basis
assert eb.TRANSLATION_FRACTION_WT == 0.30 and eb.SIGMA_EFF == 0.50
assert close(eb.F_FLEX, 0.245) and close(eb.PHI_T_WT, 0.147)
SCALE = eb.PROTEIN_CONTENT/eb.POOL_TABLE_PROTEIN_CONTENT                # 1.0889
phi_M_wt = sum(pool for pool, _ in eb.NATIVE_STEPS.values())
assert close(phi_M_wt, 0.0637*SCALE, rel=1e-9) and close(phi_M_wt, 0.06936, rel=1e-4)
assert close(eb.F_FLEX - phi_M_wt - eb.PHI_T_WT, 0.028638, rel=1e-4)   # > 0: reference feasible
assert list(eb.NATIVE_STEPS) == ['r1', 'r2', 'r3', 'r4', 'r5', 'r6']
# the pools are the report's 0.45-basis values x SCALE (report 3: linear in the protein content)
R1, R4, R5 = (eb.NATIVE_STEPS[s][0] for s in ('r1', 'r4', 'r5'))
assert close(R1, 0.044*SCALE) and eb.NATIVE_STEPS['r1'][1] == ('k_1h', 'k_1l', 'k_1e')
assert close(R4, 0.0032*SCALE) and eb.NATIVE_STEPS['r4'][1] == ()      # fixed pool; k_4 burden-free
assert close(R5, 0.0008*SCALE) and eb.NATIVE_STEPS['r5'][1] == ('k_5', 'k_5e')
assert [round(p/SCALE, 4) for p, _ in eb.NATIVE_STEPS.values()] == [0.044, 0.0032, 0.0085, 0.0032, 0.0008, 0.0040]
assert list(eb.EHRLICH_STEPS) == ['r13', 'r14', 'r15', 'r16']
assert eb.GROWTH_CAPACITIES == ('k_7', 'k_8')
assert eb.STEP_ORDER == ('r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r13', 'r14', 'r15', 'r16')
assert eb.BURDEN_COLUMNS == (
    'pool_r1', 'pool_r2', 'pool_r3', 'pool_r4', 'pool_r5', 'pool_r6',
    'pool_r13', 'pool_r14', 'pool_r15', 'pool_r16',
    'Phi_M', 'phi_T', 'F_flex', 'burden_factor', 'k_7_eff', 'k_8_eff')
PASS('sector constants (P = 0.49), closure slack 0.0286 g/gDCW, tables (0.45-basis x 1.089) and column order')

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
assert close(sig['r3'], 2.03, rel=0.02), sig                # 2.21 at P = 0.45
assert close(sig['r6'], 0.0836, rel=0.02), sig              # 0.091 at P = 0.45
assert close(sig['geomean'], math.sqrt(sig['r3']*sig['r6']))
assert 0.25 <= sig['geomean'] <= 1.0                      # within 2x of SIGMA_EFF = 0.5
assert close(sig['geomean'], 0.412, rel=0.02)               # 0.45 at P = 0.45
PASS('sigma diagnostic: sigma_r3 2.0, sigma_r6 0.084, geometric mean 0.41 in [0.25, 1.0]')

#%% 4. inert at the reference (spec 4.5): d = 1, Phi_M = Phi_M,wt, apply is the identity on k_7/k_8
bm = eb.BurdenModel.from_reference(K_REF)
assert isinstance(bm, eb.BurdenModel) and bm.sigma_eff == eb.SIGMA_EFF
assert bm.required_capacities() == ('k_1h', 'k_1l', 'k_1e', 'k_2', 'k_3',
                                    'k_5', 'k_5e', 'k_6',
                                    'k_13', 'k_14', 'k_15', 'k_16',
                                    'k_7', 'k_8')
assert 'k_4' not in bm.reference and bm.reference['k_7'] == 1.203
assert close(bm.Phi_M_wt, phi_M_wt) and bm.F_flex == eb.F_FLEX and bm.phi_T_wt == eb.PHI_T_WT
ref = bm.evaluate(K_REF)
assert isinstance(ref, eb.BurdenResult)
assert ref.burden_factor == 1.0 and ref.feasible
assert close(ref.Phi_M, phi_M_wt) and close(ref.phi_T, eb.PHI_T_WT) and ref.F_flex == eb.F_FLEX
assert ref.k_7_eff == 1.203 and ref.k_8_eff == 0.589          # exact: d == 1.0
assert close(ref.violation, phi_M_wt - eb.F_FLEX)
assert all(close(ref.pools[s], eb.NATIVE_STEPS[s][0]) for s in eb.NATIVE_STEPS)
assert all(ref.pools[s] == 0.0 for s in eb.EHRLICH_STEPS)     # A: Ehrlich off
assert bm.reference_result.burden_factor == 1.0
applied = bm.apply(K_REF)
assert applied == K_REF and applied is not K_REF               # identity, a copy
rec = ref.as_record()
assert tuple(rec) == eb.BURDEN_COLUMNS
assert rec['pool_r1'] == ref.pools['r1'] and rec['burden_factor'] == 1.0
assert rec['k_7_eff'] == 1.203 and rec['F_flex'] == eb.F_FLEX
assert bm.sigma_diagnostic == eb.anchor_sigma(K_REF)
PASS('inert at the reference: d = 1, pools = wild type, apply identity, record keyed by BURDEN_COLUMNS')

#%% 5. Phi_M is monotone non-decreasing in every capacity; k_7/k_8 do not enter Phi_M
base = bm.evaluate(K_REF).Phi_M
for name in bm.required_capacities():
    if name in eb.GROWTH_CAPACITIES:
        continue
    up = {**K_REF, name: 2.0*K_REF[name] if K_REF[name] > 0.0 else 1.0}
    down = {**K_REF, name: 0.5*K_REF[name]}
    assert bm.evaluate(up).Phi_M > base + 1e-12, name
    assert bm.evaluate(down).Phi_M <= base + 1e-15, name
assert bm.evaluate({**K_REF, 'k_7': 12.03, 'k_8': 5.89}).Phi_M == base
PASS('Phi_M monotone in every native and Ehrlich capacity; growth capacities do not enter it')

#%% 6. max-multiplier rule for the multi-capacity steps r1 and r5 (Q9)
r1 = bm.evaluate({**K_REF, 'k_1h': 3.0*0.584, 'k_1l': 2.0*1.43}).pools['r1']
assert close(r1, 3.0*R1)                                      # sized by the most-raised term
assert close(bm.evaluate({**K_REF, 'k_1h': 0.1*0.584}).pools['r1'], R1)      # one term down: block unchanged
assert close(bm.evaluate({**K_REF, 'k_1h': 0.1*0.584, 'k_1l': 0.1*1.43,
                          'k_1e': 0.1*47.1}).pools['r1'], 0.1*R1)             # all down: block shrinks
assert close(bm.evaluate({**K_REF, 'k_5e': 5.0*0.775}).pools['r5'], 5.0*R5)
assert close(bm.evaluate({**K_REF, 'k_5': 5.0*0.0104}).pools['r5'], 5.0*R5)
PASS('max-multiplier rule: r1 and r5 sized by their most-raised capacity')

#%% 7. r4 is fixed (Q10): 10x k_4 costs nothing
r4 = bm.evaluate({**K_REF, 'k_4': 48.0})
assert r4.pools['r4'] == R4 and r4.Phi_M == base
PASS('r4 pool fixed at wild type; k_4 burden-free')

#%% 8. r16 sums the KDC and the ADH, both on the KIV molar flux
r16 = bm.evaluate({**K_REF, 'k_16': 2.82}).pools['r16']
expected16 = 2.82*(71384.0/(19.0*eb.SIGMA_EFF*3600.0*116.12)
                   + 39618.0/(296.0*eb.SIGMA_EFF*3600.0*116.12))
assert close(r16, expected16)
assert close(r16, 2.82*eb.ehrlich_unit_cost('r16'))
PASS('r16 = Aro10 + Adh6 pools on the KIV flux')

#%% 9. linear squeeze (Q4/Q5d): d = 1 inside the slack, 0.5 halfway, exactly 0 at Phi_M = F_flex
cost13 = eb.ehrlich_unit_cost('r13')
slack = eb.F_FLEX - bm.Phi_M_wt - eb.PHI_T_WT                 # 0.0286
assert bm.evaluate({**K_REF, 'k_13': 0.9*slack/cost13}).burden_factor == 1.0
half = bm.evaluate({**K_REF, 'k_13': (slack + 0.5*eb.PHI_T_WT)/cost13})
assert close(half.burden_factor, 0.5, rel=1e-9) and half.feasible
assert close(half.k_7_eff, 0.5*1.203) and close(half.k_8_eff, 0.5*0.589)
k13_star = (eb.F_FLEX - bm.Phi_M_wt)/cost13                   # Phi_M == F_flex
at_cap = bm.evaluate({**K_REF, 'k_13': k13_star})
assert abs(at_cap.violation) < 1e-12 and abs(at_cap.burden_factor) < 1e-12
# (feasibility is not asserted AT the cap: rounding may land d at +-1e-16)
just_below = bm.evaluate({**K_REF, 'k_13': k13_star*(1.0 - 1e-6)})
assert just_below.feasible and 0.0 < just_below.burden_factor < 1e-4
just_over = bm.evaluate({**K_REF, 'k_13': k13_star*(1.0 + 1e-6)})
assert not just_over.feasible and just_over.burden_factor == 0.0
assert just_over.k_7_eff == 0.0 and just_over.k_8_eff == 0.0
over = bm.evaluate({**K_REF, 'k_13': 2.0*k13_star})
assert not over.feasible and over.burden_factor == 0.0 and over.violation > 0.0
PASS('d = 1 inside the slack, linear to 0, zero exactly at Phi_M = F_flex and flagged infeasible')

#%% 10. k_8 alone inflates phi_T (Q8); one machinery sized by the larger demand
k8 = bm.evaluate({**K_REF, 'k_8': 4.0*0.589})
assert close(k8.phi_T, 4.0*eb.PHI_T_WT)
assert close(k8.burden_factor, (eb.F_FLEX - phi_M_wt)/(4.0*eb.PHI_T_WT))   # 0.2987
assert close(k8.k_7_eff, k8.burden_factor*1.203) and close(k8.k_8_eff, k8.burden_factor*4.0*0.589)
both = bm.evaluate({**K_REF, 'k_7': 2.0*1.203, 'k_8': 4.0*0.589})
assert close(both.phi_T, 4.0*eb.PHI_T_WT)                     # max, not sum
assert close(bm.evaluate({**K_REF, 'k_7': 0.5*1.203, 'k_8': 0.5*0.589}).burden_factor, 1.0)
PASS('phi_T = phi_T,wt * max(k_7/k_7,ref, k_8/k_8,ref); both growth capacities derated by d')

#%% 11. 10x k_7 at wild-type enzymes is derated by the burden itself (spec 1-close)
ten = bm.evaluate({**K_REF, 'k_7': 10.0*1.203})
assert ten.feasible and close(ten.burden_factor, (eb.F_FLEX - phi_M_wt)/(10.0*eb.PHI_T_WT))
assert close(ten.burden_factor, 0.1195, rel=1e-3)
assert close(ten.k_7_eff, 1.4374, rel=1e-3)                   # = (F_flex - Phi_M)/phi_T,wt * k_7,ref
applied10 = bm.apply({**K_REF, 'k_7': 10.0*1.203})
assert close(applied10['k_7'], ten.k_7_eff) and close(applied10['k_8'], ten.k_8_eff)
assert applied10['k_3'] == 5.81 and applied10['k_13'] == 0.0     # nothing else touched
PASS('10x k_7 with wild-type enzymes: d = 0.12, k_7_eff = 1.44 (growth capped by the burden)')

#%% 12. missing keys fall back to the reference; unknown keys are ignored and passed through
empty = bm.evaluate({})
assert empty.as_record() == ref.as_record()
assert bm.apply({}) == {'k_7': 1.203, 'k_8': 0.589}
partial = bm.evaluate({'k_13': 5.0, 'K_1e': 0.12, 'k_16r': 0.0125})
assert close(partial.pools['r13'], 5.0*cost13) and close(partial.pools['r1'], R1)
assert bm.apply({'k_16r': 0.0125}) == {'k_16r': 0.0125, 'k_7': 1.203, 'k_8': 0.589}
# construction guards
try:
    eb.BurdenModel.from_reference({k: v for k, v in K_REF.items() if k != 'k_6'})
except KeyError as e:
    assert 'k_6' in str(e)
else:
    raise AssertionError('missing capacity did not raise')
try:
    eb.BurdenModel.from_reference({**K_REF, 'k_3': 0.0})
except ValueError as e:
    assert 'k_3' in str(e)
else:
    raise AssertionError('non-positive native reference did not raise')
# sigma_eff override scales the Ehrlich pools inversely
bm1 = eb.BurdenModel.from_reference(K_REF, sigma_eff=1.0)
assert close(bm1.evaluate({'k_13': 5.0}).pools['r13'], 0.5*partial.pools['r13'])
PASS('missing keys -> reference, unknown keys ignored/passed through, construction guards')

#%% 13. scenario_b_ehrlich: SCENARIO_B_EHRLICH read from nskinetics BY FILE PATH (no package import)
path13 = eb.scenarios_path()
assert path13.endswith(os.path.join('models', 's_cerevisiae_ferm_fb_inhib_mod_ibo',
                                    'scenarios.py')) and os.path.isfile(path13), path13
B = eb.scenario_b_ehrlich()
assert B == {'k_13': 5.81, 'k_14': 4.8, 'k_15': 4.8, 'k_16': 2.82, 'k_16r': 0.0125}, B
assert eb.scenario_b_ehrlich(path=path13) == B and eb.scenario_b_ehrlich() is not B
# The no-heavy-import guarantee is probed in a FRESH interpreter that loads
# enzyme_burden.py by file path (as the stdlib-only supervisor would): this
# script imports it through the biorefineries.isobutanol package, whose
# system.py imports nskinetics at module top, so sys.modules here proves
# nothing.
probe13 = (
    "import importlib.util, sys\n"
    f"spec = importlib.util.spec_from_file_location('eb_probe', {eb.__file__!r})\n"
    "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
    "b = m.scenario_b_ehrlich()\n"
    "heavy = sorted(k for k in sys.modules if k.split('.')[0] in\n"
    "               ('nskinetics', 'tellurium', 'roadrunner', 'biosteam',\n"
    "                'thermosteam', 'numpy', 'optuna'))\n"
    "print(b['k_13'], heavy)\n")
out13 = subprocess.run([sys.executable, '-c', probe13], capture_output=True, text=True)
assert out13.returncode == 0, out13.stderr
assert out13.stdout.strip() == '5.81 []', out13.stdout + out13.stderr
PASS('scenario_b_ehrlich read by file path; enzyme_burden loads with no heavy import in a fresh interpreter')

#%% 14. reports (spec 4.5 / Q11): the A reference is feasible, the B point is infeasible
report_A = bm.describe_point(bm.reference, label='scenario-A reference')
print(report_A)
assert 'scenario-A reference' in report_A and 'status: FEASIBLE' in report_A
assert 'INFEASIBLE' not in report_A
assert 'sigma_r3 = 2.03' in report_A and 'sigma_r6 = 0.084' in report_A
assert 'Phi_M = 0.0694' in report_A and 'F_flex = 0.2450' in report_A
assert 'd = 1.0000' in report_A
point_B = {**bm.reference, **B}
res_B = bm.evaluate(point_B)
ehrlich_B = sum(res_B.pools[s] for s in eb.EHRLICH_STEPS)
assert close(ehrlich_B, pool_B) and close(ehrlich_B, 0.2169, rel=0.01)
assert close(res_B.Phi_M, 0.2862, rel=0.01) and res_B.violation > 0.04   # 0.2806 / +0.056 at P = 0.45
assert not res_B.feasible and res_B.burden_factor == 0.0
assert res_B.k_7_eff == 0.0 and res_B.k_8_eff == 0.0
report_B = bm.describe_point(point_B, label='scenario-B Ehrlich constants')
print(report_B)
assert 'status: INFEASIBLE' in report_B and 'Phi_M = 0.2862' in report_B
assert 'pruned' in report_B
# k_16r is not a burden capacity: it is ignored (and passed through by apply)
assert bm.apply(point_B)['k_16r'] == 0.0125
# A B-start reference cannot build a burden model (Q11: reported, not repaired)
try:
    eb.BurdenModel.from_reference(point_B)
except ValueError as e:
    assert 'infeasible' in str(e) and '--no-burden' in str(e)
else:
    raise AssertionError('an infeasible reference did not raise')
# a tenth of B is free, two-thirds is derated (d ~ 0.24), the cap sits at
# 0.81 x B ((F_flex - Phi_M,wt)/pool_B; ~0.67 x B at P = 0.45, where 0.75 x B
# was already pruned), so 0.85 x B is pruned (spec 4.5, amended 2026-09-06)
def frac_B(f):
    return bm.evaluate({**bm.reference, **{k: f*v for k, v in B.items()}})
assert frac_B(0.10).burden_factor == 1.0
assert 0.2 < frac_B(0.65).burden_factor < 0.3
assert close((eb.F_FLEX - bm.Phi_M_wt)/pool_B, 0.81, rel=0.01)
assert frac_B(0.75).feasible and not frac_B(0.85).feasible
PASS('describe_point: A reference FEASIBLE, B point INFEASIBLE (Phi_M 0.286 > 0.245); B-start reference raises')

print(f'\nALL {n_pass} CHECKS PASSED')

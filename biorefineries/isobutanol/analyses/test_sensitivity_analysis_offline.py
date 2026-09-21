#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline verification of sensitivity_analysis (the Sobol' / Shapley
analysis of the metabolic_split_12d campaign space). No isobutanol.load(),
no simulation: sa / ko / eb are loaded BY FILE PATH (never through the
package), so this is sim-safe on any numba-cache state.
Run in a fresh kernel; exit 0 + 'ALL n CHECKS PASSED' = clean."""
import importlib.util
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PKG_DIR = os.path.dirname(HERE)

def _load(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(PKG_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

sa = _load('sa', 'sensitivity_analysis.py')
ko = _load('ko', 'kinetic_optimization.py')
eb = _load('eb', 'enzyme_burden.py')

n_pass = 0
def PASS(msg):
    global n_pass
    n_pass += 1
    print(f'PASS {n_pass}: {msg}')

#%% Shared fixture: the real split_12d space, built offline as the driver builds it
# Scenario-A-like live baselines: the B workbook's rows with the Ehrlich rates
# zeroed and k_17 at the corrected antimony value (what A carries live).
kb = {**ko.workbook_kinetic_baselines('B'),
      'k_13': 0.0, 'k_14': 0.0, 'k_15': 0.0, 'k_16': 0.0, 'k_17': 0.1077}
preset, engine_kwargs = sa.campaign_engine_kwargs(
    ko, 'ethanol_isobutanol', 'metabolic_split_12d')
space, _ = ko.build_search_space(kb, **engine_kwargs)
bm = eb.BurdenModel.from_reference(kb)
BASE_FEED = dict(target_conc=221.25, threshold_conc=217.125, spike_conc=600.0)
scalar_feasible = ko.feasibility_predicate(
    burden_on=True, volume_on=True, burden_model=bm,
    parameter_groups=engine_kwargs['parameter_groups'], kinetic_baselines=kb,
    baseline_model_kwargs=BASE_FEED, baseline_max_n_spikes=16,
    volume_cap=20.0, group_references=engine_kwargs['group_references'])
design = sa.design_record(
    search_space=space, parameter_groups=engine_kwargs['parameter_groups'],
    group_references=engine_kwargs['group_references'], kinetic_baselines=kb,
    burden_model=bm, eb=eb, baseline_model_kwargs=BASE_FEED,
    baseline_max_n_spikes=16, volume_cap=20.0,
    target_conc_max=ko.TARGET_CONC_MAX, meta={'seed': 7})

#%% 1. campaign space + JSON-able design record
assert list(space) == ['k_3', 'k_6', 'k_13', 'k_17', 'glycolysis',
                       'ehrlich_downstream', 'inhib_ethanol',
                       'inhib_isobutanol', 'inhib_acetate', 'threshold_conc',
                       'target_delta', 'max_n_spikes'], list(space)
assert preset['scenario'] == 'A'
# the preset's A-anchored override must win over the B-workbook-derived band
for name, band in (preset['param_bounds_override'] or {}).items():
    assert engine_kwargs['param_bounds_override'][name] == band, name
import json
design_rt = json.loads(json.dumps(design))
assert list(design_rt['search_space']) == list(space)
assert design_rt['meta']['seed'] == 7
PASS('campaign_engine_kwargs reproduces the 12-variable split_12d space; '
     'design_record is JSON round-trippable')

#%% 2. unit_to_values == ko.unit_to_external, elementwise
rng = np.random.default_rng(0)
U = rng.random((500, len(space)))
vals = sa.unit_to_values(U, space)
for r in range(len(U)):
    ext = ko.unit_to_external(U[r], space)
    for name in space:
        assert np.isclose(vals[name][r], ext[name], rtol=1e-12, atol=0.0), (r, name)
assert vals['max_n_spikes'].dtype.kind == 'i'
PASS('unit_to_values matches ko.unit_to_external on 500 random points (log, linear, int)')

#%% 3. feasible Sobol' stream: deterministic, nested, feasible; vectorized == scalar
def take(n, seed):
    return list(itertools.islice(sa.feasible_sobol_stream(
        space, scalar_feasible, seed, ko.unit_to_external), n))
a, b = take(40, 11), take(40, 11)
assert a == b                                   # deterministic
assert take(15, 11) == a[:15]                   # nested prefix
assert take(40, 12) != a                        # seed matters
idx = [c for c, _ in a]
assert idx == sorted(idx) and len(set(idx)) == len(idx)
assert all(scalar_feasible(v) for _, v in a)
assert all(isinstance(v['max_n_spikes'], int) for _, v in a)
vf = sa.VectorizedFeasibility.from_design(design_rt)     # from the JSON round trip
assert vf.d == 12 and vf.names == list(space)
U = np.random.default_rng(1).random((3000, vf.d))
vec = vf(U)
sca = np.array([scalar_feasible(ko.unit_to_external(u, space)) for u in U])
assert (vec == sca).all(), f'{(vec != sca).sum()} disagreements'
assert 0.05 < vec.mean() < 0.95, vec.mean()     # a real, non-trivial constraint
phi = vf.phi_M(U[:50])
for r in range(50):
    ext = ko.unit_to_external(U[r], space)
    res = bm.evaluate(ko.expand_grouped_values(
        ext, engine_kwargs['parameter_groups'], kb, engine_kwargs['group_references']))
    assert np.isclose(phi[r], res.Phi_M, rtol=1e-10), r
print(f'   feasible fraction of the campaign box: {vec.mean():.3f}')
PASS('feasible_sobol_stream deterministic / nested / feasible; VectorizedFeasibility '
     '== ko.feasibility_predicate on 3000 points; phi_M == BurdenModel.evaluate')

#%% Done
print(f'ALL {n_pass} CHECKS PASSED')

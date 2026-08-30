#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 8 -- scenario B baseline with BOTH separation trains built
(``isobutanol.load(separation_processes=('IBO_EtOH', 'ethanol'))``, the
default) but the gating splitter re-gated to send ALL broth to the
ethanol-primary train: ``sep_udct['S201'].split = 0.0``.

The scenario-B analog of smoke_test_6, and the first re-gated integrated
run with NONZERO broth isobutanol: the idle zero-flow IBO/EtOH train must
be economically equivalent to an absent one (smoke_test_4's ethanol-only
build) while HXN/WWT carry branch-2 duties for the full B-feed IBO, which
leaves via the rectifier (D303) bottoms to WWT unrecovered
(sub-decantable), so the isobutanol product is empty.

Gates (asserted inside ``load_simulate_baseline`` before it returns, so a
violation exits non-zero exactly like a traceback):

- purity-adjusted ethanol MPSP within 1% of 1.41371 (the smoke_test_4
  reference; the 2026-08-30 verification run gave 1.4137094435, matching
  the ethanol-only build to ~6 sig figs)
- isobutanol MPSP is nan (empty product)

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and each scenario mutates global V406
state. Running the file directly prints nothing -- a runner must call
``load_simulate_baseline()`` and print the returned dict.
"""
import math
import biosteam as bst
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH', 'ethanol'))

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
solve_TEA = isobutanol.system.solve_TEA
sep_udct = isobutanol.system.sep_udct
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')

def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'), # products whose MPSPs are solved
                           IRR_for_MPSP=0.15, # fixed IRR at which MPSPs are solved
                           ):
    # Re-gate BEFORE any simulation in this function: all broth to the
    # ethanol-primary train (branch 2). load() above already baseline-
    # simulated at the default split = 1.0; every pass below runs re-gated.
    sep_udct['S201'].split = 0.0

    # Scenario-B body, mirroring smoke_test_2.py.
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        'parameter-distributions_corn_IBO_EtOH_B.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
    baseline_initial = model.metrics_at_baseline()

    model_specification()

    fbs_spec.max_n_spikes = 13
    model_specification(threshold_conc=216.3, target_conc=226.3)

    results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
    MPSP_ethanol = results['MPSPs']['ethanol']
    MPSP_isobutanol = results['MPSPs']['isobutanol']
    assert abs(MPSP_ethanol - 1.41371)/1.41371 < 0.01, \
        f'ethanol MPSP {MPSP_ethanol} not within 1% of 1.41371'
    assert math.isnan(MPSP_isobutanol), \
        f'isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
    return results

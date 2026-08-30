#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 7 -- scenario A baseline with BOTH separation trains built
(``isobutanol.load(separation_processes=('IBO_EtOH', 'ethanol'))``, the
default) but the gating splitter re-gated to send ALL broth to the
ethanol-primary train: ``sep_udct['S201'].split = 0.0``.

This is the first re-gated (split < 1.0) INTEGRATED run -- the standalone
wrapper tests covered splits 0.0/0.5/1.0 without the HXN/WWT/facilities
coupling, so this is the configuration CLAUDE.md flags as a new
verification event. With all flow on branch 2, the IBO/EtOH train idles at
zero flow (design/cost skipped), D103 bottoms is empty (ProcessWaterCenter
draws more makeup water), and the isobutanol product is empty. Results are
expected to track smoke_test_5 (ethanol-only build, same scenario-A flow
routing) closely -- the 2026-08-30 verification runs agreed to ~6
significant figures (0.7979571 here vs 0.7979563 there).

Gates (asserted inside ``load_simulate_baseline`` before it returns, so a
violation exits non-zero exactly like a traceback):

- purity-adjusted ethanol MPSP within 1% of 0.79796 (matching the
  smoke_test_5 reference: an idle zero-flow IBO/EtOH branch must be
  economically equivalent to an absent one)
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

    # Scenario-A body, mirroring smoke_test_1.py.
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        'parameter-distributions_corn_IBO_EtOH_A.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
    baseline_initial = model.metrics_at_baseline()

    model_specification()

    fbs_spec.max_n_spikes = 16
    model_specification(threshold_conc=217.125, target_conc=221.25)

    results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
    MPSP_ethanol = results['MPSPs']['ethanol']
    MPSP_isobutanol = results['MPSPs']['isobutanol']
    assert abs(MPSP_ethanol - 0.79796)/0.79796 < 0.01, \
        f'ethanol MPSP {MPSP_ethanol} not within 1% of 0.79796'
    assert math.isnan(MPSP_isobutanol), \
        f'isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
    return results

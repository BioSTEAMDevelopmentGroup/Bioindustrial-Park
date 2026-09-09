#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 3 -- scenario A baseline with ONLY the IBO/EtOH separation
train: ``isobutanol.load(separation_processes=('IBO_EtOH',))``.

The scenario-A twin of smoke_test_4: dropping the pass-through gating
splitter S201 from the both-trains build (smoke_test_1's configuration)
must not move results beyond simulation tolerance, so the gate is the
both-trains scenario-A baseline. Scenario A's broth carries zero
isobutanol, so the (built) IBO/EtOH train runs in its zero-IBO
beer-column mode and the isobutanol product is empty.

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 0.86604 (the both-trains
  scenario-A baseline, reproduced to full precision by the 2026-09-02
  re-pin run, 0.866036453655084. Re-pinned 2026-09-02 for
  lang_factor=None (per-unit bare-module capital; pin under Lang 3.0
  with 2023$ prices was 0.98872), earlier that day for the 2023
  price year (stream/utility prices indexed to 2023$ with the BLS
  chemicals PPI; pin under the 2023 CEPCI with unindexed prices was
  0.85293), then for the 2023 CEPCI 797.9 (bst.CE had been
  biosteam's 567.5 default; pin under Lang 3.0 alone was 0.76073) and, earlier that day, the Lang factor 3.0 TEA
  (Huang et al. 2016; previous pin 0.84057 under corn's uncited Lang
  factor 4), and 2026-09-01 after (a) the
  kinetics-synced parameter xlsx (d467f0aa / 3720fa21) and (b) the
  vent-scrubber-bottoms recycle to the separation feed (MX8) with molar
  L/G = 2.0 wash water; the pin before that was 0.81733)
- isobutanol MPSP is nan (empty product)
- every MPSP stable against the first simulation's (relative drift
  < 5e-3, ~3 significant figures; nan stays nan)

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and the separation configuration is fixed
for the kernel's lifetime. Running the file directly prints nothing -- a
runner must call ``load_simulate_baseline()`` and print the returned list
(one solve_TEA dict per simulation).
"""
import math
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios


def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'),
                           IRR_for_MPSP=0.15,
                           n_sims=3,
                           ):
    # IBO_EtOH-only build reproduces the both-trains scenario-A baseline
    # (0.86604) to full precision; pins from the registry.
    bundle = scenarios.load_scenario('A')
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = bundle['expected']

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

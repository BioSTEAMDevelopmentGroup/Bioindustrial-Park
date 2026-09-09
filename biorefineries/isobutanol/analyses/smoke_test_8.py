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

The scenario-B twin of smoke_test_7, and the first re-gated integrated
run with NONZERO broth isobutanol: the idle zero-flow IBO/EtOH train must
be economically equivalent to an absent one (smoke_test_6's ethanol-only
build) while HXN/WWT carry branch-2 duties for the full B-feed IBO, which
leaves via the rectifier (D303) bottoms to WWT unrecovered
(sub-decantable), so the isobutanol product is empty.

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 1.4868 (the smoke_test_6
  reference; the 2026-09-03 re-pin run gave 1.4868429277, matching
  the ethanol-only build to ~6 sig figs. Re-pinned 2026-09-03 for the
  IRR-optimal batch feeding strategy from the 25x25 feeding-strategy
  sweep (target 140 g/L, threshold 34.25 g/L, max_n_spikes 0; pin under
  the 216.3/226.3/13-spike fed-batch strategy was 1.8820), 2026-09-02 for
  lang_factor=None (per-unit bare-module capital; pin under Lang 3.0
  with 2023$ prices was 2.3817), earlier that day for the 2023 price
  year (stream/utility prices indexed to 2023$ with the BLS chemicals
  PPI; pin under the 2023 CEPCI with unindexed prices was 2.1569; the
  IRR at default product prices has been unsolvable since then -- reported
  as -inf by solve_TEA from 2026-09-03, nan before),
  then for the 2023 CEPCI 797.9 (bst.CE had been biosteam's
  567.5 default; pin under Lang 3.0 alone was 1.8149) and, earlier that day, the Lang factor
  3.0 TEA (Huang et al. 2016; previous pin 2.1278 under corn's uncited
  Lang factor 4), and 2026-09-01 after
  (a) the kinetics-synced parameter xlsx (d467f0aa / 3720fa21) and (b) the
  vent-scrubber-bottoms recycle to the separation feed (MX8) with molar
  L/G = 2.0 wash water; the previous pin before that was 1.41371)
- isobutanol MPSP is nan (empty product)
- every MPSP stable against the first simulation's (relative drift
  < 5e-3, ~3 significant figures; nan stays nan)

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and each scenario mutates global V406
state. Running the file directly prints nothing -- a runner must call
``load_simulate_baseline()`` and print the returned list (one solve_TEA
dict per simulation).
"""
import math
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH', 'ethanol'))
from biorefineries.isobutanol import scenarios

sep_udct = isobutanol.system.sep_udct


def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'),
                           IRR_for_MPSP=0.15,
                           n_sims=3,
                           ):
    # Re-gate all broth to the ethanol-primary train (branch 2) before any
    # simulation in this function; matches smoke_test_6 to ~6 sig figs.
    sep_udct['S201'].split = 0.0

    bundle = scenarios.load_scenario('B')
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = {'ethanol': 1.4868, 'isobutanol': math.nan}

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

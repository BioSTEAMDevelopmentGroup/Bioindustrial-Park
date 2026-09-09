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
routing) closely -- the 2026-09-02 re-pin runs agreed to ~6
significant figures (0.8445818 here vs 0.8445811 there).

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 0.84458 (matching the
  smoke_test_5 reference: an idle zero-flow IBO/EtOH branch must be
  economically equivalent to an absent one. Re-pinned 2026-09-02 for
  lang_factor=None (per-unit bare-module capital; pin under Lang 3.0
  with 2023$ prices was 0.96339), earlier that day for the 2023 price
  year (stream/utility prices indexed to 2023$ with the BLS chemicals
  PPI; pin under the 2023 CEPCI with unindexed prices was 0.83131),
  then for the 2023 CEPCI 797.9 (bst.CE had been
  biosteam's 567.5 default; pin under Lang 3.0 alone was 0.74200) and, earlier that day, the Lang factor
  3.0 TEA (Huang et al. 2016; previous pin 0.81943 under corn's uncited
  Lang factor 4), and 2026-09-01 after
  (a) the kinetics-synced parameter xlsx (d467f0aa / 3720fa21) and (b) the
  vent-scrubber-bottoms recycle to the separation feed (MX8) with molar
  L/G = 2.0 wash water; the pin before that was 0.79796)
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
    # simulation in this function; matches smoke_test_5 to ~6 sig figs.
    sep_udct['S201'].split = 0.0

    bundle = scenarios.load_scenario('A')
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = {'ethanol': 0.84458, 'isobutanol': math.nan}

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

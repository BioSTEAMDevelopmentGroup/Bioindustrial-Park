#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 4 -- scenario B baseline with ONLY the IBO/EtOH separation train:
``isobutanol.load(separation_processes=('IBO_EtOH',))``.

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 0.59698 (was 0.61216 before the
  2026-09-13 P508 detach re-pin -- corn's orphaned rectifier-bottoms pump P508
  no longer feeds its stale build-time outlet (~14,160 kg/hr, essentially
  water) into MX5 -> M501 / WWT on every simulation -- 0.60727 before
  the 2026-09-13 r14-r16 re-pin (nskinetics b61360e / 623ff6f: r16 irreversible
  Michaelis-Menten in KIV with no k_16r / K_16i terms, r16 redox 0.363 -> 0.138
  $Red per g KIV, r14 charged one NADPH, K_14/K_15/K_16 anchored 2.64e-4/2.64e-4/
  0.034 -> 0.017/0.080/0.27 g/L, K_16i workbook row dropped), 0.67182 before the
  2026-09-13 r13 rate-law re-pin -- acetolactate synthase Michaelis-Menten in
  pyruvate, K_13 = 0.10 g/L, nskinetics a2b001c -- and 0.66834 before the
  2026-09-12 DDGS-dryer acid-split re-pin -- broth acetic acid routed to the
  dryer exhaust / X611 instead of the DDGS product -- and 0.6642 before the
  2026-09-12 nskinetics anaerobic_growth_mult 1.0 -> 0.75 re-pin; the both-trains
  scenario-B baseline; dropping the pass-through gating splitter S201 must
  not move results beyond simulation tolerance)
- purity-adjusted isobutanol MPSP within 1% of 1.2808 (was 1.2948 before the
  2026-09-13 P508 detach re-pin, 1.2864 before the
  2026-09-13 r14-r16 re-pin, 1.3018 before the
  2026-09-13 r13 rate-law re-pin, 1.2976 before the
  2026-09-12 DDGS-dryer acid-split re-pin, 1.2924 before the 2026-09-12
  anaerobic_growth_mult re-pin)
  (Both re-pinned 2026-09-03 for the IRR-optimal batch feeding strategy
  from the 25x25 feeding-strategy sweep (target 140 g/L, threshold 34.25
  g/L, max_n_spikes 0; pins under the 216.3/226.3/13-spike fed-batch
  strategy were 1.0784 / 1.7960), earlier 2026-09-02 for lang_factor=None (per-unit bare-module
  capital; pins under Lang 3.0 with 2023$ prices were 1.4479 / 2.2448),
  earlier that day for the 2023 price year (stream/utility prices
  indexed to 2023$ with the BLS chemicals PPI; pins under the 2023 CEPCI
  with unindexed prices were 1.0870 / 1.7961), then for the 2023 CEPCI 797.9 (bst.CE had been biosteam's 567.5 default; pins
  under Lang 3.0 alone were 0.81176 / 1.4618) and, earlier that day, the Lang factor 3.0 TEA (Huang et al.
  2016; previous pins 1.0521 / 1.7538 under corn's uncited Lang factor
  4), and 2026-09-01 after (a) the kinetics-synced parameter xlsx
  (d467f0aa / 3720fa21) and (b) the vent-scrubber-bottoms recycle to the
  separation feed (MX8) with molar L/G = 2.0 wash water; the previous
  pins before that were 0.39536 / 0.93321.)
- every MPSP stable against the first simulation's (relative drift
  < 5e-3, ~3 significant figures; nan stays nan)

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and the separation configuration is fixed
for the kernel's lifetime. Running the file directly prints nothing -- a
runner must call ``load_simulate_baseline()`` and print the returned list
(one solve_TEA dict per simulation).
"""
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios


def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'),
                           IRR_for_MPSP=0.15,
                           n_sims=3,
                           ):
    # IBO_EtOH-only build reproduces the both-trains scenario-B baseline
    # (0.59698 / 1.2808) to full precision; pins from the registry.
    bundle = scenarios.load_scenario('B')
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

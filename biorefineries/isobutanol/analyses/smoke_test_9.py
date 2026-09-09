#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 9 -- scenario opt_IRR (best-IRR trial #774 of the 2026-09-07
metabolic_minimal_subset IRR study) with ONLY the IBO/EtOH separation
train: ``isobutanol.load(separation_processes=('IBO_EtOH',))``.

Reproduces the study's best-objective trial: the workbook's Baseline
column carries the trial's exact simulated kinetic state (scenario-A
held/non-sampled params + the sampled rates + expanded inhibition
coefficients + burden-derated k_7/k_8), and the registry carries the
trial's feeding strategy and reproduced pins. IBO_EtOH-only == both-trains
at S201 split 1.0 (the studies' build), so the pins transfer exactly.

Loads once (at import) and simulates n_sims (default 3) times, verifying
after EACH simulation: MPSPs pinned (1%) + stable vs the first sim (5e-3),
and the study objective (IRR) reproduced within the registry tolerance.
Exit 0 = pass. Must run in a FRESH kernel.
"""
import math
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))
from biorefineries.isobutanol import scenarios

SCENARIO = 'opt_IRR'


def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'),
                           IRR_for_MPSP=0.15,
                           n_sims=3,
                           ):
    bundle = scenarios.load_scenario(SCENARIO)
    model_specification = bundle['model_specification']
    solve_TEA = bundle['solve_TEA']
    feeding_kwargs = bundle['feeding_kwargs']
    expected_MPSPs = bundle['expected']
    spec = bundle['spec']

    all_results = []
    for i in range(n_sims):
        model_specification(**feeding_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        scenarios.assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        scenarios.assert_objective_reproduced(spec, results, bundle, i+1)
        if all_results:
            scenarios.assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

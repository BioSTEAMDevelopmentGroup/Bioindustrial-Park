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

- purity-adjusted ethanol MPSP within 1% of 2.1569 (the smoke_test_6
  reference; the 2026-09-02 re-pin run gave 2.1569206141, matching
  the ethanol-only build to ~6 sig figs. Re-pinned 2026-09-02 for the
  2023 CEPCI 797.9 (bst.CE had been biosteam's 567.5 default; pin under
  Lang 3.0 alone was 1.8149) and, earlier that day, the Lang factor
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

def _assert_MPSPs_stable(reference, current, sim_number, rel_tol=5e-3):
    """Verify the current simulation's MPSPs match the first simulation's to
    ~3 significant figures (relative drift < rel_tol; nan stays nan)."""
    for ID, ref in reference['MPSPs'].items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref) or math.isnan(cur):
            assert math.isnan(ref) and math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} vs first-sim {ref} (nan mismatch)'
        else:
            assert abs(cur - ref)/abs(ref) < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} drifted from first-sim '
                 f'value {ref} beyond rel tol {rel_tol}')

def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'), # products whose MPSPs are solved
                           IRR_for_MPSP=0.15, # fixed IRR at which MPSPs are solved
                           n_sims=3, # simulations in this kernel (in-process stability check)
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
    all_results = []
    for i in range(n_sims):
        model_specification(threshold_conc=216.3, target_conc=226.3)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        MPSP_ethanol = results['MPSPs']['ethanol']
        MPSP_isobutanol = results['MPSPs']['isobutanol']
        assert abs(MPSP_ethanol - 2.1569)/2.1569 < 0.01, \
            f'sim {i+1}: ethanol MPSP {MPSP_ethanol} not within 1% of 2.1569'
        assert math.isnan(MPSP_isobutanol), \
            f'sim {i+1}: isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
        if all_results:
            _assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

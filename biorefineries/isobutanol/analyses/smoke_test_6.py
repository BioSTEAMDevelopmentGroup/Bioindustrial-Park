#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 6 -- scenario B baseline with ONLY the ethanol-primary separation
train: ``isobutanol.load(separation_processes=('ethanol',))``.

This is the first full-envelope run of the ethanol-primary train: full
scenario-B isobutanol in the broth (the standalone branch tests covered IBO
up to 0.3x of the B feed) with HXN/WWT/facilities coupled to its duties at
full flow. All broth IBO leaves via the rectifier bottoms to WWT
(sub-decantable; recovery infeasible), so the isobutanol product is empty.

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 1.8149 (the reference baseline
  for this configuration, set 3-sim-stable by the 2026-09-02 re-pin run
  and recorded in CLAUDE.md; well above the both-trains B baseline 0.81176
  because the IBO revenue is lost and its burden shifts to ethanol.
  Re-pinned 2026-09-02 for the Lang factor 3.0 TEA (Huang et al. 2016;
  previous pin 2.1278 under corn's uncited Lang factor 4), and 2026-09-01
  after (a) the kinetics-synced parameter xlsx
  (d467f0aa / 3720fa21) and (b) the vent-scrubber-bottoms recycle to the
  separation feed (MX8) with molar L/G = 2.0 wash water; the previous
  pin before that was 1.41371)
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
import biosteam as bst
from biorefineries import isobutanol
isobutanol.load(separation_processes=('ethanol',))

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
solve_TEA = isobutanol.system.solve_TEA
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
        assert abs(MPSP_ethanol - 1.8149)/1.8149 < 0.01, \
            f'sim {i+1}: ethanol MPSP {MPSP_ethanol} not within 1% of 1.8149'
        assert math.isnan(MPSP_isobutanol), \
            f'sim {i+1}: isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
        if all_results:
            _assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 5 -- scenario A baseline with ONLY the ethanol-primary separation
train: ``isobutanol.load(separation_processes=('ethanol',))``.

Scenario A's broth carries zero isobutanol at baseline, so this exercises
the ethanol-primary train (stock corn-ethanol purification wrapped with
feed-adaptive specs) as a plain ethanol-water beer-column/rectifier/sieve
train at full scenario-A flow, with HXN/WWT/facilities coupled to its
duties. The IBO/EtOH train is absent; V514 is fed a permanently-empty
placeholder, so the isobutanol product is empty.

``load_simulate_baseline`` loads once (at import) and simulates
``n_sims`` (default 3) times via ``model_specification``, verifying after
EACH simulation. Gates (a violation exits non-zero exactly like a
traceback):

- purity-adjusted ethanol MPSP within 1% of 0.81943 (the reference baseline
  for this configuration, set 3-sim-stable by the 2026-09-01 re-pin run;
  slightly below the both-trains A baseline 0.84057 -- the stock
  corn-ethanol purification train is marginally cheaper than the IBO/EtOH
  train for a zero-IBO broth. Re-pinned 2026-09-01 after (a) the
  kinetics-synced parameter xlsx (d467f0aa / 3720fa21) and (b) the
  vent-scrubber-bottoms recycle to the separation feed (MX8) with molar
  L/G = 2.0 wash water; the previous pinned value was 0.79796)
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
    # Scenario-A body, mirroring smoke_test_1.py.
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        'parameter-distributions_corn_IBO_EtOH_A.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
    baseline_initial = model.metrics_at_baseline()

    model_specification()

    fbs_spec.max_n_spikes = 16
    all_results = []
    for i in range(n_sims):
        model_specification(threshold_conc=217.125, target_conc=221.25)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        MPSP_ethanol = results['MPSPs']['ethanol']
        MPSP_isobutanol = results['MPSPs']['isobutanol']
        assert abs(MPSP_ethanol - 0.81943)/0.81943 < 0.01, \
            f'sim {i+1}: ethanol MPSP {MPSP_ethanol} not within 1% of 0.81943'
        assert math.isnan(MPSP_isobutanol), \
            f'sim {i+1}: isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
        if all_results:
            _assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)
    return all_results

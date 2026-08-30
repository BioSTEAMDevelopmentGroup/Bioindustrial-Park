#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Smoke test 4 -- scenario B baseline with ONLY the ethanol-primary separation
train: ``isobutanol.load(separation_processes=('ethanol',))``.

This is the first full-envelope run of the ethanol-primary train: full
scenario-B isobutanol in the broth (the standalone branch tests covered IBO
up to 0.3x of the B feed) with HXN/WWT/facilities coupled to its duties at
full flow. All broth IBO leaves via the rectifier bottoms to WWT
(sub-decantable; recovery infeasible), so the isobutanol product is empty.

Gates (asserted inside ``load_simulate_baseline`` before it returns, so a
violation exits non-zero exactly like a traceback):

- purity-adjusted ethanol MPSP > 0.396 (IBO revenue lost; its burden
  shifts to ethanol, so the both-trains B baseline is a floor)
- purity-adjusted ethanol MPSP < 1.5 (generous pathology cap)
- isobutanol MPSP is nan (empty product)

No pinned reference value yet: once 3-run-stable, the observed ethanol MPSP
is recorded in CLAUDE.md as this configuration's reference baseline.

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and the separation configuration is fixed
for the kernel's lifetime. Running the file directly prints nothing -- a
runner must call ``load_simulate_baseline()`` and print the returned dict.
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

def load_simulate_baseline(stream_IDs=('ethanol', 'isobutanol'), # products whose MPSPs are solved
                           IRR_for_MPSP=0.15, # fixed IRR at which MPSPs are solved
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
    model_specification(threshold_conc=216.3, target_conc=226.3)

    results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
    MPSP_ethanol = results['MPSPs']['ethanol']
    MPSP_isobutanol = results['MPSPs']['isobutanol']
    assert MPSP_ethanol > 0.396, \
        f'ethanol MPSP {MPSP_ethanol} <= 0.396 (should exceed the both-trains B baseline: IBO revenue is lost)'
    assert MPSP_ethanol < 1.5, \
        f'ethanol MPSP {MPSP_ethanol} >= 1.5 (pathology cap)'
    assert math.isnan(MPSP_isobutanol), \
        f'isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
    return results

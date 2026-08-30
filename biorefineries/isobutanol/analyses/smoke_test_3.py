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

Gates (asserted inside ``load_simulate_baseline`` before it returns, so a
violation exits non-zero exactly like a traceback):

- purity-adjusted ethanol MPSP within 1% of 0.81733 (the both-trains
  scenario-A baseline; the 2026-08-30 verification run reproduced it to
  full precision, 0.8173294638658914)
- isobutanol MPSP is nan (empty product)

Must run in a FRESH kernel/process: ``isobutanol.load(...)`` runs at import
below, rebuilds are unsupported, and the separation configuration is fixed
for the kernel's lifetime. Running the file directly prints nothing -- a
runner must call ``load_simulate_baseline()`` and print the returned dict.
"""
import math
import biosteam as bst
from biorefineries import isobutanol
isobutanol.load(separation_processes=('IBO_EtOH',))

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
    assert abs(MPSP_ethanol - 0.81733)/0.81733 < 0.01, \
        f'ethanol MPSP {MPSP_ethanol} not within 1% of 0.81733'
    assert math.isnan(MPSP_isobutanol), \
        f'isobutanol MPSP {MPSP_isobutanol} expected nan (empty product)'
    return results

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Consolidated scenario loading -- the single source of truth for every
scenario's baseline kinetics + uncertainty distributions (from its
parameter-distributions workbook) and baseline feeding strategy.

A "scenario" = (workbook: baseline kinetics + distributions)
             + (feeding strategy: max_n_spikes, threshold_conc, target_conc)
             + (expected pins for the smoke tests).

`load_scenario(scenario)` runs the recipe formerly duplicated across
smoke tests 1-8 and the BO driver against the ALREADY-LOADED model
(the caller calls `isobutanol.load(...)` first, choosing the separation
configuration). It does NOT call `isobutanol.load()`.

Registered as a submodule in `__init__.py`; model/system objects and
`kinetic_optimization` are imported lazily (only valid after load()), so
importing this module is build-free.
"""
import os
import math
from dataclasses import dataclass, field

__all__ = ('ScenarioSpec', 'SCENARIOS', 'load_scenario',
           'assert_MPSPs_pinned', 'assert_MPSPs_stable',
           'assert_objective_reproduced')

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
_WORKBOOK_DIR = os.path.join(_PKG_DIR, 'analyses', 'full',
                             'parameter_distributions')


@dataclass(frozen=True)
class ScenarioSpec:
    """One scenario. `workbook` is a bare filename under
    analyses/full/parameter_distributions/. `expected` pins the smoke-test
    MPSPs, {'ethanol': float, 'isobutanol': float} (nan for an empty
    product). For opt_* scenarios `objective_name` names an
    OBJECTIVE_REGISTRY key and `objective_value` its reproduced value.
    `spike_conc` / `stage_1_max_x` are None for every current scenario
    (load() defaults 600 g/L / 5.0 g/L); a non-None value is not yet
    supported (see load_scenario)."""
    name: str
    workbook: str
    max_n_spikes: int
    threshold_conc: float
    target_conc: float
    expected: dict
    objective_name: str = None
    objective_value: float = None
    objective_tol: float = 0.02
    spike_conc: float = None
    stage_1_max_x: float = None

    @property
    def workbook_path(self):
        return os.path.join(_WORKBOOK_DIR, self.workbook)


SCENARIOS = {
    'A': ScenarioSpec(
        name='A',
        workbook='parameter-distributions_corn_IBO_EtOH_A.xlsx',
        max_n_spikes=16, threshold_conc=217.125, target_conc=221.25,
        expected={'ethanol': 0.86604, 'isobutanol': math.nan}),
    'B': ScenarioSpec(
        name='B',
        workbook='parameter-distributions_corn_IBO_EtOH_B.xlsx',
        max_n_spikes=0, threshold_conc=34.25, target_conc=140.0,
        expected={'ethanol': 0.6642, 'isobutanol': 1.2924}),
    # opt_* scenarios are appended by analyses/build_opt_scenario_workbooks.py
    # (ask-first generator) -- see the implementation plan, Task 5.
}


def _handles():
    """Live model/system handles. Valid only after isobutanol.load()."""
    from biorefineries import isobutanol
    m = isobutanol.models
    model = m.models_EtOH_IBO_corn.model
    return dict(model=model,
                namespace_dict=m.namespace_dict,
                fbs_spec=m.fbs_spec,
                model_specification=model.specification,
                solve_TEA=isobutanol.system.solve_TEA,
                V406=model.system.flowsheet.V406,
                tea=model.system.TEA)


def load_scenario(scenario, *, apply=True):
    """Load `scenario`'s baseline kinetics + distributions from its
    workbook, set the baseline feeding strategy, and (when apply=True) run
    one baseline `model_specification`. Returns a bundle the caller reuses
    for its stability loop. Requires the model to be built already
    (isobutanol.load(...))."""
    spec = SCENARIOS[scenario]
    h = _handles()
    model = h['model']

    model.parameters = ()
    model.load_parameter_distributions(spec.workbook_path, h['namespace_dict'])
    model.metrics_at_baseline()          # workbook Baseline column SETS kinetics

    h['fbs_spec'].max_n_spikes = spec.max_n_spikes
    if spec.spike_conc is not None or spec.stage_1_max_x is not None:
        raise NotImplementedError(
            'load_scenario does not yet override spike_conc / stage_1_max_x; '
            'every registered scenario pins them at the load() defaults '
            '(600 g/L / 5.0 g/L).')
    feeding_kwargs = dict(threshold_conc=spec.threshold_conc,
                          target_conc=spec.target_conc)
    if apply:
        h['model_specification'](**feeding_kwargs)

    return dict(spec=spec, feeding_kwargs=feeding_kwargs,
                expected=spec.expected, model=model,
                model_specification=h['model_specification'],
                solve_TEA=h['solve_TEA'], fbs_spec=h['fbs_spec'],
                namespace_dict=h['namespace_dict'], V406=h['V406'],
                tea=h['tea'])


def assert_MPSPs_pinned(expected, current, sim_number, rel_tol=0.01):
    """Verify current MPSPs against pinned baselines (within rel_tol; nan
    stays nan). Pins live in the SCENARIOS registry."""
    for ID, ref in expected.items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref):
            assert math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} expected nan (empty product)'
        else:
            assert abs(cur - ref)/ref < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} not within rel tol '
                 f'{rel_tol} of pinned {ref}')


def assert_MPSPs_stable(reference, current, sim_number, rel_tol=5e-3):
    """Verify current MPSPs match the first simulation's to ~3 sig figs
    (relative drift < rel_tol; nan stays nan)."""
    for ID, ref in reference['MPSPs'].items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref) or math.isnan(cur):
            assert math.isnan(ref) and math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} vs first-sim {ref} (nan mismatch)'
        else:
            assert abs(cur - ref)/abs(ref) < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} drifted from first-sim '
                 f'value {ref} beyond rel tol {rel_tol}')


def assert_objective_reproduced(spec, results, bundle, sim_number, rel_tol=None):
    """For an opt_* scenario, recompute the study's objective from the
    current simulated model state (via kinetic_optimization.OBJECTIVE_REGISTRY,
    the exact study definition) and verify it reproduces the best trial's
    recorded value. No-op for scenarios without an objective (A/B)."""
    if spec.objective_name is None:
        return
    tol = spec.objective_tol if rel_tol is None else rel_tol
    from biorefineries.isobutanol import kinetic_optimization as ko
    getter = ko.OBJECTIVE_REGISTRY[spec.objective_name]['getter']
    h = {'V406': bundle['V406'], 'tea': bundle['tea'],
         'latest_TEA_solution': results}
    value = getter(h)
    ref = spec.objective_value
    assert abs(value - ref)/abs(ref) < tol, \
        (f'sim {sim_number}: {spec.name} objective {spec.objective_name} '
         f'{value} not within rel tol {tol} of recorded best-trial {ref}')

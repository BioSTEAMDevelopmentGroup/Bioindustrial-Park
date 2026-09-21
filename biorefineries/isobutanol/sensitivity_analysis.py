#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Variance-based global sensitivity analysis over a kinetic-optimization
campaign space (spec docs/superpowers/specs/2026-09-20-sobol-sensitivity-
split12d-design.md). SIM-FREE and package-free: stdlib + numpy / scipy /
sklearn only, loadable by file path; callers pass the kinetic_optimization
module (`ko`) and the enzyme_burden module (`eb`) in where needed.

The campaign space is a box with a closed-form feasibility constraint
(enzyme burden Phi_M < F_flex; fed-batch volume-ratio bound <= cap), so the
inputs are DEPENDENT on the feasible domain. The quantities estimated here
are well defined under dependence: the closed subset index
S_u = Var(E[Y | X_u]) / Var(Y) (pick-freeze with conditional rejection
sampling) and the exact Shapley effects built from all 2^d of them."""
import math

import numpy as np

__all__ = ('campaign_engine_kwargs', 'design_record', 'unit_to_values',
           'feasible_sobol_stream', 'VectorizedFeasibility')

#%% Campaign space

_PRESET_ENGINE_KEYS = ('include_params', 'exclude_params', 'multiplier_bounds',
                       'rate_multiplier_bounds', 'rate_params',
                       'parameter_multiplier_bounds', 'stage_1_max_x_bounds',
                       'parameter_groups', 'group_multiplier_bounds',
                       'spike_delta_bounds', 'group_references',
                       'param_bounds_override')

def campaign_engine_kwargs(ko, study_target_products, study_type):
    """(preset, engine_kwargs) exactly as analyses/optimize_kinetics_BO.run
    derives them for a named preset: the preset's keys, with
    param_bounds_override = the bounds workbook's bands
    (ko.workbook_kinetic_bounds of preset['kinetic_bounds_scenario'])
    UPDATED by the preset's own override (the scenario-A-anchored IBO
    pathway bands win). engine_kwargs feeds ko.build_search_space and
    ko._prepare_optimization unchanged."""
    preset = ko.resolve_study_preset(study_target_products, study_type)
    engine_kwargs = {key: preset[key] for key in _PRESET_ENGINE_KEYS}
    derived = ko.workbook_kinetic_bounds(
        preset['kinetic_bounds_scenario'],
        multiplier_bounds=engine_kwargs['multiplier_bounds'],
        rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds'],
        rate_params=engine_kwargs['rate_params'],
        parameter_multiplier_bounds=engine_kwargs['parameter_multiplier_bounds'])
    derived.update(engine_kwargs['param_bounds_override'] or {})
    engine_kwargs['param_bounds_override'] = derived
    return preset, engine_kwargs

def design_record(*, search_space, parameter_groups, group_references,
                  kinetic_baselines, burden_model, eb, baseline_model_kwargs,
                  baseline_max_n_spikes, volume_cap, target_conc_max,
                  meta=None):
    """Everything stage 2 needs to rebuild the feasible domain WITHOUT the
    model, as a JSON-able dict (written by stage 1 as <study>_design.json).
    The burden block is the closed form BurdenModel.pools evaluates: native
    steps pool_wt x max_c(value_c / reference_c), Ehrlich steps value x
    unit cost; feasible <=> Phi_M < F_flex."""
    for legacy in ('spike_delta', 'target_conc', 'threshold_delta', 'spike_conc'):
        if legacy in search_space:
            raise NotImplementedError(
                f'{legacy!r} is sampled: VectorizedFeasibility covers only a '
                'pinned spike feed on the threshold-anchored scheme.')
    return dict(
        search_space={name: {k: (bool(v) if k in ('log', 'int') else float(v))
                             for k, v in sp.items()}
                      for name, sp in search_space.items()},
        parameter_groups={g: list(m) for g, m in (parameter_groups or {}).items()},
        group_references={g: {m: float(r) for m, r in refs.items()}
                          for g, refs in (group_references or {}).items()},
        kinetic_baselines={n: float(v) for n, v in kinetic_baselines.items()},
        burden=dict(
            F_flex=float(burden_model.F_flex),
            reference={n: float(v) for n, v in burden_model.reference.items()},
            native_steps={step: [float(pool_wt), list(caps)]
                          for step, (pool_wt, caps) in eb.NATIVE_STEPS.items()},
            ehrlich_costs={cap: float(eb.ehrlich_unit_cost(step, burden_model.sigma_eff))
                           for step, (cap, _, _) in eb.EHRLICH_STEPS.items()}),
        volume=dict(spike_conc=float(baseline_model_kwargs['spike_conc']),
                    baseline_max_n_spikes=int(baseline_max_n_spikes),
                    volume_cap=float(volume_cap),
                    target_conc_max=float(target_conc_max)),
        meta=dict(meta or {}))

#%% Unit cube <-> decision values (vectorized twin of ko.unit_to_external)

def unit_to_values(U, search_space):
    """{name: ndarray} for unit-cube rows U (n, d), columns in search-space
    order: log floats exp(ln lo + u (ln hi - ln lo)) clipped into [lo, hi],
    ints lo + floor(u (hi - lo + 1)) clipped to hi (int64), linear floats
    lo + u (hi - lo) -- ko.unit_to_internal's measure, elementwise."""
    U = np.clip(np.asarray(U, dtype=float), 0.0, 1.0)
    values = {}
    for j, (name, sp) in enumerate(search_space.items()):
        u, lo, hi = U[:, j], sp['low'], sp['high']
        if sp.get('int'):
            lo, hi = int(lo), int(hi)
            values[name] = np.minimum(hi, lo + np.floor(u*(hi - lo + 1))).astype(np.int64)
        elif sp['log']:
            x = np.exp(math.log(lo) + u*(math.log(hi) - math.log(lo)))
            values[name] = np.clip(x, lo, hi)
        else:
            values[name] = lo + u*(hi - lo)
    return values

def feasible_sobol_stream(search_space, is_feasible, seed, unit_to_external,
                          block=1024):
    """Endless iterator of (candidate_index, values): scrambled-Sobol'
    candidates of the unit cube in blocks of `block` (a power of 2),
    mapped by `unit_to_external` (pass ko.unit_to_external -- the ONE
    definition of the campaign measure) and kept when `is_feasible(values)`.
    Rejection of a uniform stream = uniform on the feasible set. Fully
    deterministic in `seed`, so the accepted sequence is NESTED (its first n
    points are a valid design for every n) and resuming = regenerating and
    skipping."""
    from scipy.stats import qmc
    engine = qmc.Sobol(d=len(search_space), scramble=True, seed=seed)
    candidate_index = 0
    while True:
        for u in engine.random(block):
            values = unit_to_external(u, search_space)
            if is_feasible(values):
                yield candidate_index, values
            candidate_index += 1

#%% Vectorized feasibility (stage 2 needs ~1e7-1e8 checks)

class VectorizedFeasibility:
    """Batch twin of ko.feasibility_predicate(burden_on=True, volume_on=True)
    over unit-cube arrays, rebuilt from a design_record dict. The scalar
    predicate stays the authority: the offline test and stage 2 both
    cross-check this class against it / against the recorded Phi_M."""

    def __init__(self, design):
        self.search_space = design['search_space']
        self.names = list(self.search_space)
        self.d = len(self.names)
        self.parameter_groups = design['parameter_groups']
        self.group_references = design['group_references']
        self.kinetic_baselines = design['kinetic_baselines']
        self.burden = design['burden']
        self.volume = design['volume']

    @classmethod
    def from_design(cls, design):
        return cls(design)

    def _applied(self, values):
        """ko.expand_grouped_values, vectorized: member = basis x multiplier
        (basis = the group's reference when it has one, else the live
        baseline); individually sampled kinetics pass through."""
        applied = {}
        for name, x in values.items():
            if name in self.parameter_groups:
                refs = self.group_references.get(name)
                for member in self.parameter_groups[name]:
                    basis = (self.kinetic_baselines[member] if refs is None
                             else refs[member])
                    applied[member] = basis*x
            elif name in self.kinetic_baselines:
                applied[name] = x
        return applied

    def _phi_M(self, values):
        applied, reference = self._applied(values), self.burden['reference']
        n = len(next(iter(values.values())))
        def val(c):
            return applied[c] if c in applied else np.full(n, reference[c])
        phi = np.zeros(n)
        for pool_wt, caps in self.burden['native_steps'].values():
            # A native step with no capacities (r4, the fixed Ald6 pool) has a
            # constant unit multiplier -- BurdenModel.pools uses default=1.0 --
            # so guard the empty reduction (np.max([]) has no identity).
            multiplier = (np.max([val(c)/reference[c] for c in caps], axis=0)
                          if caps else 1.0)
            phi += pool_wt*multiplier
        for cap, cost in self.burden['ehrlich_costs'].items():
            phi += val(cap)*cost
        return phi

    def phi_M(self, U):
        return self._phi_M(unit_to_values(U, self.search_space))

    def __call__(self, U):
        values = unit_to_values(U, self.search_space)
        ok = self._phi_M(values) < self.burden['F_flex']
        vol = self.volume
        thr = values['threshold_conc']
        tgt = np.minimum(vol['target_conc_max'], thr + values['target_delta'])
        spk = vol['spike_conc']
        n_spikes = (values['max_n_spikes'] if 'max_n_spikes' in values
                    else np.full(len(thr), vol['baseline_max_n_spikes']))
        with np.errstate(divide='ignore', invalid='ignore'):
            log_ratio = n_spikes*np.log((spk - thr)/(spk - tgt))
        log_ratio = np.where(spk <= tgt, np.inf, log_ratio)
        log_ratio = np.where(n_spikes <= 0, 0.0, log_ratio)
        return ok & (log_ratio <= math.log(vol['volume_cap']))

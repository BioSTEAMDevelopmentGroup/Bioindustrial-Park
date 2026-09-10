#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Bayesian (Optuna TPE) global optimization of the fermentation kinetic
parameters (k_*/K_* on V406's tellurium model; by default the driver
restricts the set to the scenario workbook's rows via include_params) plus
the feeding-strategy
variables (threshold_conc, target_conc via target_delta, spike_conc via
spike_delta -- feasible-by-construction threshold < target < spike --
and the integer max_n_spikes cap), against a named or custom objective --
to prioritize metabolic-engineering / bioprocess research directions.

Named study presets (resolve_study_preset) pick the search set and bands
on two axes: study_target_products ('ethanol_only' = scenario-A workbook
rows; 'ethanol_isobutanol' = scenario-B rows, i.e. plus the Ehrlich block
and the isobutanol-inhibition coefficients; both start at the A baseline)
and study_type ('metabolic' = capacity, product-inhibition, lethality and
substrate-regulation roles; 'metabolic_protein' = plus affinity and
product self-inhibition), with the rate constants on [1e-3x, 10x] and the
inhibition coefficients / K_* terms on [0.1x, 10x] log bands; k_10 (the
active-biomass decay capacity) is excluded from every preset by default
(DEFAULT_EXCLUDED_PARAMETERS -- a lower decay rate is a free lunch, not an
engineering target; study-name tag `_xk10`). Roles come from nskinetics'
parameter_categories table, read by file path.

Import is free of side effects and does not require optuna (imported lazily
inside run_kinetic_optimization); the pure logic here is exercised by
analyses/test_kinetic_optimization_offline.py without a biorefinery load.
The engine itself (run_kinetic_optimization) requires
biorefineries.isobutanol.load() to have been called first; the canonical
entry point is analyses/optimize_kinetics_BO.py.

Every trial runs the FULL system simulation + one side-effect-free TEA
solve, for kinetic-level and system-level objectives alike (the 'level'
tag is metadata). Kinetic parameters and feeding specs are restored to
their scenario baselines in a `finally` after every run.
"""
import csv
import importlib.util
import json
import math
import os
import re

import numpy as np

__all__ = ('OBJECTIVE_REGISTRY', 'TRACKED_METRICS',
           'discover_kinetic_parameters', 'build_search_space',
           'parameter_distributions_workbook', 'workbook_kinetic_baselines',
           'kinetic_param_names_from_scenario', 'workbook_kinetic_bounds',
           'DEFAULT_RATE_MULTIPLIER_BOUNDS',
           'DEFAULT_PARAMETER_MULTIPLIER_BOUNDS',
           'DEFAULT_SATURATION_MULTIPLIER_BOUNDS',
           'DEFAULT_EXCLUDED_PARAMETERS', 'excluded_parameters_tag',
           'OPERATING_VARIABLES', 'DEFAULT_STAGE_1_MAX_X_BOUNDS',
           'DEFAULT_SPIKE_DELTA_BOUNDS', 'DEFAULT_GROUP_MULTIPLIER_BOUNDS',
           'expand_grouped_values',
           'RATE_CONSTANT_ROLES', 'INHIBITION_COEFFICIENT_ROLES',
           'kinetic_parameter_roles_path', 'kinetic_parameter_roles',
           'rate_constant_names',
           'STUDY_TARGET_PRODUCTS', 'STUDY_TYPE_ROLES',
           'STUDY_TYPE_OPTIONS', 'EFFECTOR_ORDER',
           'METABOLIC_MINIMAL_SUBSET_RATES', 'METABOLIC_MINIMAL_SUBSET_GROUPS',
           'kinetic_parameter_effectors', 'study_type_name_defaults',
           'DEFAULT_STUDY_TARGET_PRODUCTS', 'DEFAULT_STUDY_TYPE',
           'resolve_study_preset', 'default_study_name',
           'BURDEN_STUDY_SUFFIX',
           'baseline_decision_point', 'knockout_probe_points',
           'clip_to_search_space', 'seed_points_from_trajectory',
           'seed_points_tag',
           'trajectory_columns', 'check_trajectory_header',
           'append_trajectory_row', 'load_trajectory',
           'inflight_path_for', 'write_inflight', 'clear_inflight',
           'recover_inflight',
           'get_handles', 'run_kinetic_optimization', 'restore_baseline',
           'plot_optimization_trajectories', 'plot_parameter_trajectory',
           'plot_best_vs_baseline',
           'pca_decision_matrix', 'plot_pca_projection',
           'StallGuard', 'attempt_outcome',
           'search_space_distributions', 'draw_uniform_feasible',
           'feasible_candidate_mask', 'feasible_tpe_sampler')

FEEDING_VARIABLES = ('threshold_conc', 'target_delta', 'spike_delta',
                     'max_n_spikes')
#: Feeding decision variables of studies started before 2026-08-31
#: (target-anchored parameterization); resumable via the legacy bounds
#: kwargs of build_search_space / run_kinetic_optimization.
LEGACY_FEEDING_VARIABLES = ('target_conc', 'threshold_delta', 'spike_conc')

#: Non-feeding OPERATING decision variables of the study presets (since
#: 2026-09-06 pm): stage_1_max_x, the fermentor's aerobic stage-1 biomass
#: cutoff [g/L] -- the nskinetics event `x >= stage_1_max_x` sets
#: is_aerobic = 0 (whichever of it and `time >= stage_1_max_time`, 25 h,
#: fires first ends stage 1). It is an nskinetics OPERATION parameter,
#: not a kinetic role, so it is never in kinetic_baselines / the burden
#: model and gets no knockout probe. It is applied through the
#: V406.stage_1_max_x PROPERTY -- whose setter mirrors the value onto the
#: kinetic model (r_te.stage_1_max_x) AND the AerationSpec that sizes
#: the air supply (K330/V330) -- never by setattr on r_te, which would
#: leave aeration sizing on the stale cutoff. Baseline = the nskinetics
#: factory default 5.0 g/L (system.py never passes it). Sampled
#: log-scale on DEFAULT_STAGE_1_MAX_X_BOUNDS; in the space only when
#: build_search_space is given stage_1_max_x_bounds (None = absent, so
#: every study started before it keeps its columns; the presets pass the
#: default); tagged `_s1x{lo}-{hi}` into preset study names by
#: default_study_name.
OPERATING_VARIABLES = ('stage_1_max_x',)
DEFAULT_STAGE_1_MAX_X_BOUNDS = (1.0, 50.0)

#: Applied-concentration envelope (g/L) for the current parameterization:
#: target_conc = min(TARGET_CONC_MAX, threshold_conc + target_delta);
#: spike_conc = clip(target_conc + spike_delta, SPIKE_CONC_MIN,
#: SPIKE_CONC_MAX).
TARGET_CONC_MAX = 300.0
SPIKE_CONC_MIN = 50.0
SPIKE_CONC_MAX = 600.0
#: Default band of the spike_delta feeding variable (build_search_space /
#: run_kinetic_optimization `spike_delta_bounds`); the presets return it
#: too, except metabolic_minimal, which pins the spike (None).
DEFAULT_SPIKE_DELTA_BOUNDS = (0.5, 595.0)

#: Default log-scale band of a PARAMETER GROUP multiplier (x every
#: member's baseline; build_search_space `group_multiplier_bounds`): the
#: metabolic_minimal preset samples one multiplier per inhibition-effector
#: family on it (0.2x-2x, baseline 1.0), so a family's intra-family
#: ratios are preserved while its overall strength varies.
DEFAULT_GROUP_MULTIPLIER_BOUNDS = (0.2, 2.0)

def expand_grouped_values(values, parameter_groups, kinetic_baselines):
    """Copy of the decision dict `values` with every PARAMETER-GROUP
    multiplier replaced by its members' applied values (member baseline x
    multiplier, from `kinetic_baselines`), the group key dropped and every
    other entry (individual kinetics, feeding and operating variables)
    passed through unchanged. `parameter_groups` is the
    {group_name: [member names]} mapping given to build_search_space;
    None / {} gives a plain copy. Pure: used by the engine's objective
    (what reaches the model and the burden), by the `applied_<member>`
    trajectory columns and by the offline test. The inverse for the
    baseline point is trivial (every group multiplier = 1.0)."""
    if not parameter_groups:
        return dict(values)
    groups = dict(parameter_groups)
    out = {}
    for name, value in values.items():
        if name in groups:
            for member in groups[name]:
                out[member] = kinetic_baselines[member]*value
        else:
            out[name] = value
    return out

#: Default log-scale multiplier bands (× the workbook baseline) of the
#: named study presets (resolve_study_preset), assigned by nskinetics
#: ROLE (kinetic_parameter_roles) since 2026-09-06:
#:   1. RATE CONSTANTS (role capacity: k_1h, k_2, ..., k_13-k_16) --
#:      DEFAULT_RATE_MULTIPLIER_BOUNDS, [1e-3×, 10×] (1e-5× until later
#:      on 2026-09-06): an effective
#:      knock-out is reachable (still log-uniform below and above the
#:      baseline).
#:   2. INHIBITION COEFFICIENTS (roles product_inhibition and lethality:
#:      k_1ie, k_1ii, k_7ii, k_10ie, ...) -- DEFAULT_SATURATION_
#:      MULTIPLIER_BOUNDS, [0.1×, 10×].
#:   3. REGULATION, AFFINITY and SELF-INHIBITION TERMS (K_1i, K_2i, ...,
#:      K_1e, ..., K_6e, K_16i) -- [0.1×, 10×] (a zero saturation constant
#:      has no engineering meaning).
#: Before 2026-09-06 the rate band applied to every lowercase 'k_' name
#: (inhibition coefficients included); that prefix rule is what
#: build_search_space still applies when `rate_params` is None (legacy
#: and older-study resumes). The legacy single band of
#: build_search_space's multiplier_bounds default, (0.1, 10.0), is
#: unchanged.
DEFAULT_RATE_MULTIPLIER_BOUNDS = (1e-3, 10.0)
DEFAULT_SATURATION_MULTIPLIER_BOUNDS = (0.1, 10.0)

#: Per-parameter multiplier bands of the study presets, {name: (m_lo,
#: m_hi)} x baseline, taking precedence over the ROLE band of that name
#: (build_search_space / workbook_kinetic_bounds `parameter_multiplier_
#: bounds`; an absolute param_bounds_override entry still wins). k_10 is
#: the capacity of r10, the ACTIVE-BIOMASS DECAY step: it is a rate
#: constant by role, but a near-zero decay rate (the 1e-3x floor of the
#: rate band, an effectively immortal culture) is not an engineering
#: target, so it keeps the saturation band [0.1×, 10×] (its knockout
#: probe therefore sits at 0.1x, a knock-down). The presets hand each
#: study a COPY of this table; changing an entry changes the sampled
#: space of every preset without changing its columns, so a changed
#: table needs a fresh study name (the header guard cannot tell).
DEFAULT_PARAMETER_MULTIPLIER_BOUNDS = {'k_10': (0.1, 10.0)}

#: Kinetic parameters every study preset keeps OUT of the search space
#: (resolve_study_preset `exclude_params`; since 2026-09-06 pm). k_10, the
#: active-biomass decay capacity, is a free lunch for the optimizer -- a
#: lower decay rate is a longer-lived culture with no engineering lever
#: behind it (the best trial of the aborted 2026-09-06 role-band study was
#: its 0.1x knock-down probe) -- so it stays at the scenario baseline and
#: gets no knockout probe (the probes follow the search space). The
#: workbook set of include_params is unchanged (29/40/40/56); the
#: effective sampled set is include minus exclude (28/39/39/55). A caller
#: re-includes k_10 with exclude_params=() (driver run(exclude_params=()),
#: supervisor bare --exclude-params); it then lands on the
#: DEFAULT_PARAMETER_MULTIPLIER_BOUNDS band above, which is kept for that
#: case. An exclusion REMOVES a trajectory-CSV column, so a study of the
#: old set cannot be resumed silently (the header guard raises); the
#: effective exclusion set is tagged into every preset-derived study name
#: (excluded_parameters_tag, `_xk10`) so the default name can launch next
#: to the older studies at all.
DEFAULT_EXCLUDED_PARAMETERS = ('k_10',)

def excluded_parameters_tag(names):
    """Study-name tag of an exclusion set: '_x' + the names joined by '+'
    with their underscores dropped, in input order (('k_10',) -> '_xk10';
    ('k_10', 'k_7') -> '_xk10+k7'); '' for None or an empty set."""
    names = tuple(names or ())
    if not names:
        return ''
    return '_x' + '+'.join(name.replace('_', '') for name in names)

#: nskinetics roles (parameter_categories.ROLES) that make a parameter a
#: RATE CONSTANT -- the only names the k_* rate band applies to under the
#: study presets (rate_constant_names) -- and the roles of the
#: INHIBITION COEFFICIENTS (documentation / reporting; they simply take
#: the saturation band). Every other role (affinity,
#: substrate_regulation, product_self_inhibition) is a K_* term on the
#: saturation band too.
RATE_CONSTANT_ROLES = ('capacity',)
INHIBITION_COEFFICIENT_ROLES = ('product_inhibition', 'lethality')

#%% Objective registry and tracked metrics
# Getters are callables over a `handles` dict (see get_handles below):
#   'V406': the fermentation unit (.nsk_results_specific_tau_dict, .tau)
#   'tea': the system TEA (.TCI)
#   'latest_TEA_solution': {'IRR': ..., 'MPSPs': {'ethanol': ..,
#                           'isobutanol': ..}}, refreshed once per trial.
# This indirection keeps every getter testable offline with fakes.

def _nsk(handles):
    return handles['V406'].nsk_results_specific_tau_dict

# 'level' ('kinetic' vs 'system') is metadata only: every trial runs the full
# system simulation AND one TEA solve regardless (the system-level metrics
# IRR/TCI/... in TRACKED_METRICS are always recorded), so a 'kinetic' objective
# costs the same per trial as a 'system' one -- it does not skip the TEA solve.
OBJECTIVE_REGISTRY = {
    'IBO yield': dict(
        getter=lambda h: _nsk(h)['y_IBO_glu_added'],
        direction='maximize', level='kinetic', units='g-IBO/g-sugars'),
    'IBO titer': dict(
        getter=lambda h: _nsk(h)['[s_IBO]'],
        direction='maximize', level='kinetic', units='g-IBO/L-broth'),
    'IBO productivity': dict(
        getter=lambda h: _nsk(h)['[s_IBO]']/_nsk(h)['time'],
        direction='maximize', level='kinetic', units='g-IBO/L-broth/h'),
    'IBO yield x titer': dict(
        getter=lambda h: _nsk(h)['y_IBO_glu_added']*_nsk(h)['[s_IBO]'],
        direction='maximize', level='kinetic',
        units='(g-IBO/g-sugars)(g-IBO/L-broth)'),
    'EtOH yield': dict(
        getter=lambda h: _nsk(h)['y_EtOH_glu_added'],
        direction='maximize', level='kinetic', units='g-EtOH/g-sugars'),
    'EtOH titer': dict(
        getter=lambda h: _nsk(h)['[s_EtOH]'],
        direction='maximize', level='kinetic', units='g-EtOH/L-broth'),
    'EtOH productivity': dict(
        getter=lambda h: _nsk(h)['prod_EtOH'],
        direction='maximize', level='kinetic', units='g-EtOH/L-broth/h'),
    'Combined yield': dict(
        getter=lambda h: _nsk(h)['y_EtOH_IBO_glu_added'],
        direction='maximize', level='kinetic', units='g-EtOH-and-IBO/g-sugars'),
    'Cell density': dict(
        getter=lambda h: _nsk(h)['[x]'],
        direction='maximize', level='kinetic', units='g-cell/L-broth'),
    'IRR': dict(
        getter=lambda h: h['latest_TEA_solution']['IRR'],
        direction='maximize', level='system', units=''),
    'EtOH MPSP': dict(
        getter=lambda h: h['latest_TEA_solution']['MPSPs']['ethanol'],
        direction='minimize', level='system', units='$/kg'),
    'IBO MPSP': dict(
        getter=lambda h: h['latest_TEA_solution']['MPSPs']['isobutanol'],
        direction='minimize', level='system', units='$/kg'),
    'TCI': dict(
        getter=lambda h: h['tea'].TCI/1e6,
        direction='minimize', level='system', units='MM$'),
    }

#: Metrics recorded for EVERY trial (spec trajectory (ii)-(vi) + extras).
TRACKED_METRICS = {name: OBJECTIVE_REGISTRY[name]['getter'] for name in
                   ('IBO yield', 'IBO titer', 'IBO productivity',
                    'EtOH yield', 'EtOH titer', 'EtOH productivity',
                    'Cell density', 'IRR', 'TCI')}
TRACKED_METRICS['tau'] = lambda h: h['V406'].tau
TRACKED_METRICS['n_glu_spikes'] = lambda h: _nsk(h)['curr_n_glu_spikes']

#%% Search space

def discover_kinetic_parameters(r_te):
    """{name: baseline} for every kinetic parameter of the tellurium model
    `r_te` -- the k_* rates and K_* constants, selected with the same name
    rule as utils.generate_save_kinetic_parameter_distributions
    (name[:2].lower() == 'k_'). Baselines are the CURRENT values, so call
    this only after the scenario's parameter distributions have been
    loaded (model.load_parameter_distributions + metrics_at_baseline)."""
    return {p: float(getattr(r_te, p))
            for p in r_te.getGlobalParameterIds()
            if p[:2].lower() == 'k_'}

def _rate_predicate(rate_params, rate_prefix='k_'):
    """name -> bool: does the k_* RATE band apply to `name`? An explicit
    `rate_params` (any iterable of names, possibly empty) is the exact
    set; None is the pre-2026-09-06 lowercase-prefix rule."""
    if rate_params is None:
        return lambda name: name.startswith(rate_prefix)
    rate_set = frozenset(rate_params)
    return lambda name: name in rate_set

def build_search_space(kinetic_baselines,
                       multiplier_bounds=(0.1, 10.0),
                       param_bounds_override=None,
                       exclude_params=(),
                       include_params=None,
                       threshold_conc_bounds=(0.0, 300.0),
                       target_delta_bounds=(5.0, 500.0),
                       spike_delta_bounds=DEFAULT_SPIKE_DELTA_BOUNDS,
                       max_n_spikes_bounds=(0, 50),
                       target_conc_bounds=None,
                       threshold_delta_bounds=None,
                       spike_conc_bounds=None,
                       rate_multiplier_bounds=None,
                       rate_params=None,
                       parameter_multiplier_bounds=None,
                       stage_1_max_x_bounds=None,
                       parameter_groups=None,
                       group_multiplier_bounds=DEFAULT_GROUP_MULTIPLIER_BOUNDS,
                       ):
    """Build the decision-variable space: {name: {'low', 'high', 'log'}}
    (integer variables additionally carry 'int': True).

    Kinetic parameters default to [m_lo*baseline, m_hi*baseline] sampled
    log-scale (research-opportunity width); entries in
    `param_bounds_override` ({name: (low, high)}) use those absolute
    bounds instead (log-scale only if low > 0). A parameter with a
    nonpositive baseline and no override cannot use the multiplier band
    and is EXCLUDED with a printed warning.

    `rate_multiplier_bounds` (None = use `multiplier_bounds` for every
    name, the single-band behaviour of studies started before
    2026-09-04) is a separate (m_lo, m_hi) band for the RATE constants,
    so the study presets can let a rate reach an effective knock-out
    (DEFAULT_RATE_MULTIPLIER_BOUNDS) while every other kinetic
    parameter keeps `multiplier_bounds`. WHICH names are rate constants
    is `rate_params` (an iterable of names; the presets pass the
    capacity-role rows of the workbook, rate_constant_names, so the
    inhibition coefficients k_*i* stay on `multiplier_bounds`); None
    keeps the pre-2026-09-06 rule -- every name starting with lowercase
    'k_' -- required to resume a study started under it. Either way the
    rate band applies only where a band applies at all (after
    include/exclude/override).

    `parameter_multiplier_bounds` ({name: (m_lo, m_hi)}, None/{} = none)
    is a PER-PARAMETER multiplier band that takes precedence over the
    role band of that name (rate or saturation, under either rate
    rule) -- the presets pass DEFAULT_PARAMETER_MULTIPLIER_BOUNDS, which
    keeps k_10 (active-biomass decay) on [0.1x, 10x] while the other
    rate constants take DEFAULT_RATE_MULTIPLIER_BOUNDS. Precedence:
    include/exclude > absolute `param_bounds_override` > per-parameter
    band > role band; an entry for a name that is absent, excluded or
    overridden is ignored, and a nonpositive baseline is still excluded
    (a multiplier band needs one).

    `exclude_params` names are
    always excluded (silently). `include_params`
    (None = no restriction) is a WHITELIST: a kinetic parameter is placed
    in the space only if its name is in it, and every other kinetic
    parameter is excluded silently -- the driver passes the kinetic rows
    of a scenario's parameter-distributions workbook here so the search
    set matches the curated uncertainty set. The feeding variables below
    are never filtered by it. The feeding variables (all linear) make
    the required ordering threshold < target < spike feasible BY
    CONSTRUCTION: threshold_conc is sampled absolutely, the target is
    sampled as target_delta above the threshold
    (target_conc = min(TARGET_CONC_MAX, threshold_conc + target_delta)),
    and the spike as spike_delta above the target
    (spike_conc = clip(target_conc + spike_delta, SPIKE_CONC_MIN,
    SPIKE_CONC_MAX)) -- spanning the applied envelope threshold [0, 300],
    target [5, 300], spike [50, 600] g/L at the default bounds.
    max_n_spikes (the glucose-spike cap, fbs_spec.max_n_spikes) is an
    INTEGER variable (0 = forced batch); pass max_n_spikes_bounds=None to
    pin it at the scenario baseline instead.
    `stage_1_max_x_bounds` (None, the default = NOT in the space) adds the
    OPERATING variable stage_1_max_x (see OPERATING_VARIABLES), sampled
    LOG-scale on the given (lo, hi) g/L (0 < lo < hi, ValueError
    otherwise), after max_n_spikes; the presets pass
    DEFAULT_STAGE_1_MAX_X_BOUNDS. Like the feeding variables it is never
    filtered by include_params / exclude_params. Absent by default so
    every study started before 2026-09-06 pm keeps its columns.

    Passing any of the LEGACY kwargs (target_conc_bounds,
    threshold_delta_bounds, spike_conc_bounds; unspecified ones fall back
    to the pre-2026-08-31 defaults (180, 300)/(0.5, 30)/(200, 800))
    builds the legacy target-anchored parameterization instead
    (target_conc absolute, threshold_conc = max(0, target_conc -
    threshold_delta), spike_conc absolute) -- required to resume a study
    started under it.

    `spike_delta_bounds=None` (since 2026-09-07; a tuple by default,
    DEFAULT_SPIKE_DELTA_BOUNDS) REMOVES spike_delta from the space: the
    spike concentration is then pinned at the engine's scenario-baseline
    snapshot (600 g/L on both scenarios), like max_n_spikes_bounds=None /
    stage_1_max_x_bounds=None pin theirs. Threshold-anchored scheme only
    (ValueError with a legacy kwarg). A missing column, so the header
    guard refuses to resume a study that sampled it.

    `parameter_groups` (None = none; since 2026-09-07) is a
    {group_name: [member kinetic names]} mapping (or a list of pairs).
    Each group is ONE decision variable, a log-scale MULTIPLIER on
    `group_multiplier_bounds` (0 < lo < hi; DEFAULT_GROUP_MULTIPLIER_
    BOUNDS) applied to every member's baseline (expand_grouped_values;
    baseline value 1.0), so the members move together and keep their
    ratios. Members are REMOVED from the individual kinetic space (never
    sampled on their own; not listed in `excluded`); groups are appended
    after the individual kinetic entries and before the feeding
    variables, in input order. ValueError: a group name colliding with a
    kinetic parameter, feeding, legacy-feeding or operating variable; an
    empty group; a member not in `kinetic_baselines`, listed in two
    groups, in `exclude_params`, or with a nonpositive baseline. A
    member absent from `include_params` is still grouped (the whitelist
    governs INDIVIDUAL sampling; the group is the explicit instruction),
    and a member's `param_bounds_override` entry is ignored (the
    driver's preset path passes absolute workbook bounds for every row).
    The metabolic_minimal preset groups the inhibition coefficients by
    effector (resolve_study_preset).

    Returns (space, excluded_parameter_names)."""
    param_bounds_override = dict(param_bounds_override or {})
    parameter_multiplier_bounds = dict(parameter_multiplier_bounds or {})
    m_lo, m_hi = multiplier_bounds
    r_lo, r_hi = (multiplier_bounds if rate_multiplier_bounds is None
                  else rate_multiplier_bounds)
    is_rate = _rate_predicate(rate_params)
    parameter_groups = {str(group): list(members)
                        for group, members in dict(parameter_groups or {}).items()}
    grouped = {}  # member name -> group name
    if parameter_groups:
        g_lo, g_hi = group_multiplier_bounds
        if not (0.0 < g_lo < g_hi):
            raise ValueError('group_multiplier_bounds must satisfy 0 < lo < hi '
                             '(a log-scale multiplier band); got '
                             f'{tuple(group_multiplier_bounds)!r}')
        reserved = (set(kinetic_baselines) | set(FEEDING_VARIABLES)
                    | set(LEGACY_FEEDING_VARIABLES) | set(OPERATING_VARIABLES))
        for group, members in parameter_groups.items():
            if group in reserved:
                raise ValueError(f'parameter group name {group!r} collides '
                                 'with a kinetic parameter, feeding or '
                                 'operating variable name')
            if not members:
                raise ValueError(f'parameter group {group!r} is empty')
            for member in members:
                if member not in kinetic_baselines:
                    raise ValueError(f'parameter group {group!r}: member '
                                     f'{member!r} is not a kinetic parameter '
                                     'of the model')
                if member in grouped:
                    raise ValueError(f'kinetic parameter {member!r} is listed '
                                     f'in two parameter groups ({grouped[member]!r} '
                                     f'and {group!r})')
                if member in exclude_params:
                    raise ValueError(f'kinetic parameter {member!r} is both in '
                                     f'parameter group {group!r} and in '
                                     'exclude_params')
                if kinetic_baselines[member] <= 0.0:
                    raise ValueError(f'parameter group {group!r}: member '
                                     f'{member!r} has a nonpositive baseline '
                                     f'({kinetic_baselines[member]}); a group '
                                     'multiplier needs a positive one')
                grouped[member] = group
    space, excluded = {}, []
    for name, baseline in kinetic_baselines.items():
        if name in grouped:
            continue  # sampled through its group, never individually
        if include_params is not None and name not in include_params:
            excluded.append(name)
        elif name in exclude_params:
            excluded.append(name)
        elif name in param_bounds_override:
            lo, hi = param_bounds_override[name]
            space[name] = dict(low=lo, high=hi, log=lo > 0.0)
        elif baseline <= 0.0:
            print(f'Warning: kinetic parameter {name} has a nonpositive '
                  f'baseline ({baseline}) and no param_bounds_override '
                  'entry; excluding it from the search space.')
            excluded.append(name)
        else:
            if name in parameter_multiplier_bounds:
                lo_m, hi_m = parameter_multiplier_bounds[name]
            else:
                lo_m, hi_m = (r_lo, r_hi) if is_rate(name) else (m_lo, m_hi)
            space[name] = dict(low=lo_m*baseline, high=hi_m*baseline,
                               log=True)
    for group in parameter_groups:
        space[group] = dict(low=float(g_lo), high=float(g_hi), log=True)
    if (target_conc_bounds is not None or threshold_delta_bounds is not None
            or spike_conc_bounds is not None):
        # Legacy target-anchored parameterization (resumes of studies
        # started before 2026-08-31).
        if spike_delta_bounds is None:
            raise ValueError('spike_delta_bounds=None (spike pinned at the '
                             'baseline) is only supported by the '
                             'threshold-anchored scheme; do not combine it '
                             'with target_conc_bounds / threshold_delta_bounds '
                             '/ spike_conc_bounds')
        tcb = target_conc_bounds or (180.0, 300.0)
        tdb = threshold_delta_bounds or (0.5, 30.0)
        scb = spike_conc_bounds or (200.0, 800.0)
        space['target_conc'] = dict(low=tcb[0], high=tcb[1], log=False)
        space['threshold_delta'] = dict(low=tdb[0], high=tdb[1], log=False)
        space['spike_conc'] = dict(low=scb[0], high=scb[1], log=False)
    else:
        space['threshold_conc'] = dict(low=threshold_conc_bounds[0],
                                       high=threshold_conc_bounds[1],
                                       log=False)
        space['target_delta'] = dict(low=target_delta_bounds[0],
                                     high=target_delta_bounds[1],
                                     log=False)
        if spike_delta_bounds is not None:  # None = spike pinned at the baseline
            space['spike_delta'] = dict(low=spike_delta_bounds[0],
                                        high=spike_delta_bounds[1], log=False)
    if max_n_spikes_bounds is not None:
        space['max_n_spikes'] = dict(low=int(max_n_spikes_bounds[0]),
                                     high=int(max_n_spikes_bounds[1]),
                                     log=False, int=True)
    if stage_1_max_x_bounds is not None:
        lo, hi = stage_1_max_x_bounds
        if not (0.0 < lo < hi):
            raise ValueError('stage_1_max_x_bounds must satisfy 0 < lo < hi '
                             f'(a log-scale band, g/L); got '
                             f'{tuple(stage_1_max_x_bounds)!r}')
        space['stage_1_max_x'] = dict(low=float(lo), high=float(hi), log=True)
    return space, excluded

def baseline_decision_point(search_space, kinetic_baselines,
                            baseline_model_kwargs,
                            baseline_max_n_spikes=None,
                            baseline_stage_1_max_x=None,
                            parameter_groups=None):
    """The scenario baseline expressed in decision-variable coordinates
    for `search_space` (either feeding parameterization) -- suitable for
    study.enqueue_trial, so a fresh study evaluates the baseline itself
    as trial 0. Only names present in the search space are included.
    Every value is CLIPPED into its [low, high] bounds: an out-of-range
    enqueued value is a hard optuna ValueError on a log-scale variable
    and is silently replaced by a random draw on a linear one, so when a
    baseline lies outside the space (e.g. a zero-baseline kinetic
    parameter under absolute param_bounds_override bounds), trial 0
    evaluates the nearest in-bounds point to the baseline instead.
    `baseline_stage_1_max_x` (the live V406.stage_1_max_x, 5.0 g/L at
    the factory default) fills point['stage_1_max_x'] when that
    operating variable is in the space and a value is given (clipped
    into its band like every other entry).

    `parameter_groups` (the build_search_space mapping; None = none)
    gives every group in the space its baseline multiplier 1.0. A space
    without spike_delta (spike_delta_bounds=None) gets no such entry."""
    point = {name: kinetic_baselines[name]
             for name in search_space if name in kinetic_baselines}
    for group in dict(parameter_groups or {}):
        if group in search_space:
            point[group] = 1.0
    thr = baseline_model_kwargs['threshold_conc']
    tgt = baseline_model_kwargs['target_conc']
    spk = baseline_model_kwargs['spike_conc']
    if 'threshold_conc' in search_space:  # threshold-anchored scheme
        point['threshold_conc'] = thr
        point['target_delta'] = tgt - thr
        if 'spike_delta' in search_space:  # absent = spike pinned at the baseline
            point['spike_delta'] = spk - tgt
    elif 'target_conc' in search_space:  # legacy target-anchored scheme
        point['target_conc'] = tgt
        point['threshold_delta'] = tgt - thr
        point['spike_conc'] = spk
    if 'max_n_spikes' in search_space and baseline_max_n_spikes is not None:
        point['max_n_spikes'] = int(baseline_max_n_spikes)
    if 'stage_1_max_x' in search_space and baseline_stage_1_max_x is not None:
        point['stage_1_max_x'] = float(baseline_stage_1_max_x)
    point, _ = clip_to_search_space(point, search_space)
    return point

def clip_to_search_space(point, search_space):
    """Clip every entry of `point` (a {name: value} decision point; only
    names in `search_space` are kept) into its [low, high] bounds, integer
    variables cast to int. Returns (clipped point, [names that moved]).
    The clipping rule of baseline_decision_point, shared with the seed
    points: an out-of-range enqueued value is a hard optuna ValueError on
    a log-scale variable and is silently replaced by a random draw on a
    linear one."""
    clipped, moved = {}, []
    for name, value in point.items():
        if name not in search_space:
            continue
        sp = search_space[name]
        new = min(max(value, sp['low']), sp['high'])
        new = int(new) if sp.get('int') else new
        if new != value:
            moved.append(name)
        clipped[name] = new
    return clipped, moved

def seed_points_from_trajectory(csv_path, trial_numbers, search_space,
                                label=None, parameter_groups=None):
    """Decision points of the trials `trial_numbers` of a DONOR study's
    trajectory CSV, in `search_space` coordinates -- the SEED points a
    fresh study enqueues right after its knockout probes (see
    run_kinetic_optimization(seed_from=...)), so TPE starts with a
    foothold in a basin another study found instead of having to hit it
    by chance (the 2026-09-06 IRR study's 575 completed trials had none
    with an isobutanol titer above 10 g/L, while the IBO-yield-x-titer
    study's trial 1553 in that basin scores a higher IRR than the IRR
    study's own optimum: the ethanol-only optimum was a local one).

    Read sim-free with the stdlib csv module (the driver builds the seeds
    inside the child, the offline test on a temp file). The decision
    columns are the CSV columns between 'state' and 'objective'
    (trajectory_columns). Every name of `search_space` must be a donor
    decision column (ValueError otherwise -- the donor sampled a
    different space; a donor column NOT in the space, e.g. a k_10 the new
    study pins, is dropped and noted); every value is clipped into the
    new study's bounds by clip_to_search_space (a donor under wider bands
    seeds the nearest in-bounds point; noted). A trial number absent from
    the CSV is a ValueError; a seed of any state is accepted (a LOST/FAIL
    row still has a decision vector) but a non-COMPLETE state is noted.
    Since 2026-09-07 a space name absent from the donor's decision
    columns is read from the donor's derived `applied_<name>` column when
    there is one (the donor grouped it; noted), so a grouped
    (metabolic_minimal) study can seed a study that samples those members
    individually; the reverse raises. A donor without spike_delta (pinned
    spike) can only seed a space without it.

    The ValueError names WHY each missing column is missing: a feeding /
    operating variable of this space (spike_delta, stage_1_max_x) means
    the DONOR pinned it (its value is never a column), while a group
    multiplier means the donor sampled the members individually (no
    unique inverse). `parameter_groups` ({group: [members]} of this
    space, optional) is what tells the two apart; without it a missing
    non-feeding name keeps the group wording.

    Returns ({label: point}, [notes]): labels are
    '{label or the CSV's study name}#{trial_number}' (the donor study
    name = the file name minus '_trajectory.csv'), stored as the optuna
    user attr 'seed' of the enqueued trial; notes are human-readable
    lines for the engine to print."""
    stem = os.path.basename(csv_path)
    if stem.endswith('_trajectory.csv'):
        stem = stem[:-len('_trajectory.csv')]
    label = label or stem
    wanted = {int(n) for n in trial_numbers}
    if not wanted:
        raise ValueError(f'seed_from {csv_path}: no trial numbers given')
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader, None)
        if header is None or 'state' not in header or 'objective' not in header:
            raise ValueError(f'{csv_path} is not a trajectory CSV (no '
                             "'state'/'objective' columns)")
        decision = header[header.index('state') + 1:header.index('objective')]
        # A name the donor did not sample individually may still be
        # recorded as a DERIVED applied_<name> column (a member of one of
        # the donor's parameter groups): read it from there. The reverse
        # (a group multiplier of this space, sampled individually by the
        # donor) has no unique inverse and stays an error.
        column_for, via_applied, missing = {}, [], []
        for name in search_space:
            if name in decision:
                column_for[name] = name
            elif f'applied_{name}' in header:
                column_for[name] = f'applied_{name}'
                via_applied.append(name)
            else:
                missing.append(name)
        if missing:
            # Why the donor has no column, per missing name: a FEEDING /
            # OPERATING variable this space samples means the donor PINNED
            # it (e.g. spike_delta under spike_delta_bounds=None: pinned at
            # the donor's own scenario baseline, which the CSV never
            # records); a group multiplier of this space means the donor
            # sampled its members individually, and a group multiplier has
            # no unique inverse. `parameter_groups` (when given) names this
            # space's groups; without it a missing non-feeding name keeps
            # the group wording.
            pinned = [name for name in missing
                      if name in FEEDING_VARIABLES
                      or name in LEGACY_FEEDING_VARIABLES
                      or name in OPERATING_VARIABLES]
            groups_here = dict(parameter_groups or {})
            grouped_missing = [name for name in missing
                               if name not in pinned
                               and (name in groups_here
                                    or parameter_groups is None)]
            reasons = []
            if pinned:
                reasons.append(
                    f'{pinned}: the donor PINNED these (no column -- e.g. '
                    'spike_delta under spike_delta_bounds=None, pinned at '
                    'the donor scenario baseline), so it cannot seed a '
                    'study that samples them')
            if grouped_missing:
                reasons.append(
                    f'{grouped_missing}: a group multiplier has no unique '
                    'inverse in a donor that sampled its members '
                    'individually')
            other = [name for name in missing
                     if name not in pinned and name not in grouped_missing]
            if other:
                reasons.append(f'{other}: the donor sampled a different '
                               'search space')
            raise ValueError(
                f'seed_from {csv_path}: the donor study has no decision '
                f'column for {missing} -- ' + '; '.join(reasons)
                + '. Seeds must come from a study of the same columns (or '
                'one that recorded the name as an applied_<name> column) '
                f'(donor decision columns: {decision})')
        rows = {}
        i_trial, i_state = header.index('trial_number'), header.index('state')
        for row in reader:
            if not row:
                continue
            try:
                n = int(float(row[i_trial]))
            except ValueError:
                continue
            if n in wanted and n not in rows:
                rows[n] = row
    absent = sorted(wanted - set(rows))
    if absent:
        raise ValueError(f'seed_from {csv_path}: trial_number {absent} not '
                         'in the trajectory')
    dropped = [name for name in decision if name not in search_space]
    points, notes = {}, []
    if dropped:
        notes.append(f'{label}: donor decision columns not in this '
                     f'search space are ignored: {dropped}')
    if via_applied:
        notes.append(f'{label}: read from the donor\'s derived applied_* '
                     f'columns (members of its parameter groups): {via_applied}')
    for n in sorted(rows):
        row = rows[n]
        raw = {}
        for name in search_space:
            cell = row[header.index(column_for[name])]
            if cell == '':
                raise ValueError(f'seed_from {csv_path}: trial {n} has no '
                                 f'value for {name}')
            raw[name] = float(cell)
        point, moved = clip_to_search_space(raw, search_space)
        key = f'{label}#{n}'
        points[key] = point
        state = row[i_state]
        if state != 'COMPLETE':
            notes.append(f'{key}: donor state {state!r} (not COMPLETE)')
        if moved:
            notes.append(f'{key}: clipped into this study\'s bounds: {moved}')
    return points, notes

def seed_points_tag(n_seeds):
    """Study-name tag of a seeded study: `_seed{n}` for n > 0 seed points
    (default_study_name / the driver), nothing for 0 / None. Seeds change
    a study's trajectory but not its columns, so the header guard cannot
    tell a seeded study from the unseeded one of the same name -- the tag
    is what launches it next to that study. It carries only the COUNT: a
    different panel of the same size needs an explicit study_name."""
    return f'_seed{int(n_seeds)}' if n_seeds else ''

def knockout_probe_points(search_space, baseline_point, rate_prefix='k_',
                          rate_params=None):
    """The single-knockout probes of a study: for every log-scale RATE
    constant of `search_space` -- the names in `rate_params` when given
    (the presets' capacity-role rows, so inhibition coefficients k_*i*
    get no probe: their 0.1x floor is a knock-down, not a knock-out),
    else every name starting with `rate_prefix` (the pre-2026-09-06
    rule) -- a copy
    of `baseline_point` with that one variable at its band FLOOR
    (search_space[name]['low']) and every other decision variable at the
    baseline. Enqueued right after trial 0 by run_kinetic_optimization
    (enqueue_knockouts=True; default False since 2026-09-07), they teach the TPE sampler
    the single-parameter lethality map explicitly -- which rate can be
    knocked down alone and which cannot -- before it starts drawing
    multi-parameter points, instead of relying on random startup draws
    that (under a wide band) knock out many rates at once.

    The probe value is the floor of the study's own band, so it is always
    in range (an out-of-range enqueued value is a hard optuna error on a
    log variable): under the preset band (DEFAULT_RATE_MULTIPLIER_BOUNDS,
    1e-3x) it is an effective knock-out, under a 0.1x band a 10x
    knock-down (so the probe of a rate on a per-parameter band --
    k_10 under DEFAULT_PARAMETER_MULTIPLIER_BOUNDS -- is a knock-down). Saturation constants K_*, the feeding variables, integer
    and linear-scale variables get no probe. A rate whose baseline already
    sits at (or below) its floor -- e.g. the Ehrlich rates of the
    ethanol_isobutanol preset, clipped up to the floor from scenario A's
    zeros by baseline_decision_point -- is skipped too (its probe would
    duplicate the baseline point) and reported in the second return value.

    Returns ({name: point}, [names skipped as already at the floor]), both
    in search-space order; `baseline_point` is not modified."""
    is_rate = _rate_predicate(rate_params, rate_prefix)
    probes, at_floor = {}, []
    for name, sp in search_space.items():
        if (not is_rate(name) or not sp.get('log')
                or sp.get('int') or name not in baseline_point):
            continue
        if baseline_point[name] <= sp['low']:
            at_floor.append(name)
            continue
        point = dict(baseline_point)
        point[name] = sp['low']
        probes[name] = point
    return probes, at_floor

#%% Scenario parameter-distribution workbooks

#: Load-statement pattern of a kinetic parameter row in the scenario
#: parameter-distributions workbooks (written by
#: utils.generate_save_kinetic_parameter_distributions).
_TE_LOAD_STATEMENT = re.compile(
    r'V406\.nsk_kinetic_model\._te\.([A-Za-z0-9_]+)\s*=\s*x')

def parameter_distributions_workbook(scenario):
    """Absolute path of the scenario's parameter-distributions workbook
    (analyses/full/parameter_distributions/
    parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx)."""
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'analyses', 'full', 'parameter_distributions',
        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx')

def workbook_kinetic_baselines(scenario):
    """Ordered {te parameter name: workbook Baseline} for every kinetic
    parameter row (load statement `V406.nsk_kinetic_model._te.<name> = x`)
    of the scenario's parameter-distributions workbook -- the curated set
    of kinetic parameters treated as free/uncertain (the physically
    constrained k_6r, k_16r, K_2, K_9 are absent since commit 1e4efee1).
    A plain file read, no simulation, so it can describe a scenario other
    than the loaded one (the driver's start-at-A / set-from-B mode)."""
    import pandas as pd
    baselines = {}
    for _, row in pd.read_excel(
            parameter_distributions_workbook(scenario)).iterrows():
        match = _TE_LOAD_STATEMENT.fullmatch(
            str(row['Load statement']).strip())
        if match:
            baselines[match.group(1)] = float(row['Baseline'])
    return baselines

def kinetic_param_names_from_scenario(scenario):
    """The kinetic parameter names of the scenario's workbook, in workbook
    order (see workbook_kinetic_baselines) -- the include_params of a
    workbook-restricted search space."""
    return list(workbook_kinetic_baselines(scenario))

def workbook_kinetic_bounds(scenario, multiplier_bounds=(0.1, 10.0),
                            rate_multiplier_bounds=None, rate_params=None,
                            parameter_multiplier_bounds=None):
    """Absolute (low, high) bounds -- the multiplier band around the
    scenario workbook's baseline -- for every positive-baseline kinetic
    row of that workbook (workbook_kinetic_baselines), keyed by te name
    and in workbook order. Same rule as build_search_space: the RATE
    constants -- the names in `rate_params` when given (the presets pass
    rate_constant_names of the workbook rows), else every name starting
    with lowercase 'k_' (the pre-2026-09-06 rule) -- use
    `rate_multiplier_bounds` when it is given, every other name
    (inhibition coefficients k_*i*, K_*) uses `multiplier_bounds`; a
    name in `parameter_multiplier_bounds` ({name: (m_lo, m_hi)}, the
    presets' DEFAULT_PARAMETER_MULTIPLIER_BOUNDS -- k_10 on 0.1x-10x)
    uses that band instead of either. A plain
    file read (no simulation), so it can parameterize a run of a
    DIFFERENT scenario: passed as param_bounds_override it gives the
    IBO-pathway rates zeroed on the model under scenario A their
    scenario-B search bands instead of degenerate zero-baseline
    exclusion (the driver's start-at-A / set-from-B mode and every study
    preset)."""
    m_lo, m_hi = multiplier_bounds
    r_lo, r_hi = (multiplier_bounds if rate_multiplier_bounds is None
                  else rate_multiplier_bounds)
    is_rate = _rate_predicate(rate_params)
    parameter_multiplier_bounds = dict(parameter_multiplier_bounds or {})
    bounds = {}
    for name, baseline in workbook_kinetic_baselines(scenario).items():
        if baseline > 0.0:
            if name in parameter_multiplier_bounds:
                lo_m, hi_m = parameter_multiplier_bounds[name]
            else:
                lo_m, hi_m = (r_lo, r_hi) if is_rate(name) else (m_lo, m_hi)
            bounds[name] = (lo_m*baseline, hi_m*baseline)
    return bounds

#%% nskinetics parameter roles (read by file path)

#: Path of the role table relative to the nskinetics package directory.
_ROLE_TABLE_RELPATH = ('models', 's_cerevisiae_ferm_fb_inhib_mod_ibo',
                       'parameter_categories.py')

def kinetic_parameter_roles_path():
    """Absolute path of nskinetics' parameter_categories.py for the shipped
    S. cerevisiae ethanol/isobutanol model (the source of truth for each
    kinetic parameter's role: capacity, affinity, substrate_regulation,
    product_inhibition, product_self_inhibition, lethality,
    lethality_threshold, initial_state). Located via the package's spec
    WITHOUT importing it (importlib.util.find_spec only resolves the
    file)."""
    spec = importlib.util.find_spec('nskinetics')
    if spec is None or not spec.origin:
        raise ImportError('nskinetics is not installed (needed for the '
                          'kinetic-parameter role table).')
    return os.path.join(os.path.dirname(spec.origin), *_ROLE_TABLE_RELPATH)

_kinetic_parameter_table_cache = None

def _kinetic_parameter_table(path=None):
    """nskinetics' parameter_categories.KINETIC_PARAMETERS ({name:
    ParameterInfo(role, reactions, modules, effector)}) executed BY PATH
    (spec_from_file_location + exec_module): the file is pure Python
    (stdlib math + dataclasses only) and documented as readable without
    loading the model, whereas importing it through the package pulls
    tellurium/roadrunner/biosteam (~15 s) -- which the stdlib-only
    supervisor and the offline test must never do. The default-path
    table is cached after the first call; an explicit `path` (tests) is
    always read afresh and never cached."""
    global _kinetic_parameter_table_cache
    if path is None and _kinetic_parameter_table_cache is not None:
        return _kinetic_parameter_table_cache
    table_path = kinetic_parameter_roles_path() if path is None else path
    spec = importlib.util.spec_from_file_location(
        '_nskinetics_parameter_categories', table_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    table = dict(module.KINETIC_PARAMETERS)
    if path is None:
        _kinetic_parameter_table_cache = table
    return table

_kinetic_parameter_roles_cache = None

def kinetic_parameter_roles(path=None):
    """{kinetic parameter name: role} from nskinetics'
    parameter_categories.KINETIC_PARAMETERS, in table order, read by
    file path (_kinetic_parameter_table; cached for the default path,
    same pattern as kinetic_parameter_effectors)."""
    global _kinetic_parameter_roles_cache
    if path is None and _kinetic_parameter_roles_cache is not None:
        return _kinetic_parameter_roles_cache
    roles = {name: info.role
             for name, info in _kinetic_parameter_table(path).items()}
    if path is None:
        _kinetic_parameter_roles_cache = roles
    return roles

_kinetic_parameter_effectors_cache = None

def kinetic_parameter_effectors(path=None):
    """{kinetic parameter name: effector or None} from the same table
    (the `effector` attribute: 'ethanol' / 'isobutanol' / 'acetate' on
    the inhibition, lethality and self-inhibition rows, 'glucose' /
    'acetaldehyde' on the regulation rows, None on capacities and
    affinities), in table order; cached for the default path like the
    roles. The metabolic_minimal preset groups the inhibition
    coefficients by it (resolve_study_preset)."""
    global _kinetic_parameter_effectors_cache
    if path is None and _kinetic_parameter_effectors_cache is not None:
        return _kinetic_parameter_effectors_cache
    effectors = {name: getattr(info, 'effector', None)
                 for name, info in _kinetic_parameter_table(path).items()}
    if path is None:
        _kinetic_parameter_effectors_cache = effectors
    return effectors

def rate_constant_names(names, roles=None):
    """The RATE CONSTANTS among `names` (role in RATE_CONSTANT_ROLES, i.e.
    the capacities k_1h, k_2, ..., k_13-k_16), in input order -- the names
    the presets' k_* band (DEFAULT_RATE_MULTIPLIER_BOUNDS) applies to;
    inhibition coefficients (k_1ie, k_1ii, k_10ie, ...) and every K_*
    term are left on the saturation band. `roles` (default
    kinetic_parameter_roles()) is the {name: role} table; a name absent
    from it raises KeyError(name) so an unclassified parameter can never
    silently land on the wrong band."""
    if roles is None:
        roles = kinetic_parameter_roles()
    out = []
    for name in names:
        if name not in roles:
            raise KeyError(f'{name!r} has no entry in the nskinetics '
                           'kinetic-parameter role table; cannot assign '
                           'its multiplier band.')
        if roles[name] in RATE_CONSTANT_ROLES:
            out.append(name)
    return out

#%% Study presets

#: `study_target_products` axis: which products the strain is engineered
#: for. Both presets START at the scenario-A baseline (the current
#: ethanol strain); the parameter SET comes from the named workbook
#: (A: native pathways only -- no Ehrlich block, no isobutanol-inhibition
#: coefficients; B: plus both). Exclusions (k_6r, k_16r, K_2, K_9, P_10*,
#: docs/reports/kinetic-parameter-exclusions.md) are inherited from the
#: workbooks and never re-enter.
STUDY_TARGET_PRODUCTS = {
    'ethanol_only': {'scenario': 'A', 'parameter_set_scenario': 'A'},
    'ethanol_isobutanol': {'scenario': 'A', 'parameter_set_scenario': 'B'},
}

#: `study_type` axis: which kinetic-parameter ROLES (nskinetics
#: parameter_categories) are decision variables. 'metabolic' = expression
#: and tolerance engineering: enzyme abundance (capacities incl. the
#: growth rate), whole-cell product tolerance (exponential inhibition and
#: lethality coefficients) and regulatory wiring (glucose/acetaldehyde
#: repression constants K_1i, K_2i, K_5i, K_9i) -- no change to any
#: enzyme's intrinsic kinetics. 'metabolic_protein' additionally allows
#: enzyme redesign: substrate affinities and competitive product
#: self-inhibition constants. Membership is by ROLE, not name prefix.
STUDY_TYPE_ROLES = {
    'metabolic': ('capacity', 'product_inhibition', 'lethality',
                  'substrate_regulation'),
    'metabolic_protein': ('capacity', 'product_inhibition', 'lethality',
                          'substrate_regulation', 'affinity',
                          'product_self_inhibition'),
    # 'metabolic_minimal' (2026-09-07): the 'metabolic' roles minus the
    # substrate-regulation terms, with the inhibition coefficients
    # sampled as ONE multiplier per effector family (STUDY_TYPE_OPTIONS).
    'metabolic_minimal': ('capacity', 'product_inhibition', 'lethality'),
    # 'metabolic_minimal_subset' (2026-09-07): NO role filter -- the set
    # is listed outright (METABOLIC_MINIMAL_SUBSET_RATES / _GROUPS via
    # STUDY_TYPE_OPTIONS; the empty tuple would otherwise select no
    # rows, so it cannot be mistaken for a role-filtered type). The
    # table stays the registry of valid study_type values.
    'metabolic_minimal_subset': (),
}

#: The STANDALONE metabolic_minimal_subset preset (2026-09-07): its
#: variables are listed outright here, NOT derived from metabolic_minimal
#: (every one of them also appears there, but nothing in the code relates
#: the two). The rate constants, in this order, each on the rate band
#: (DEFAULT_RATE_MULTIPLIER_BOUNDS); resolve_study_preset intersects the
#: list with the target's workbook (ethanol_only lacks k_13-k_16).
METABOLIC_MINIMAL_SUBSET_RATES = ('k_1l', 'k_1h', 'k_1e', 'k_3', 'k_6',
                                  'k_13', 'k_14', 'k_15', 'k_16')
#: Its inhibition-effector groups: ONE log multiplier per group on
#: (0.2, 2.0) x every member's LIVE baseline (expand_grouped_values),
#: group and member order as written; a member absent from the target's
#: workbook is dropped and an emptied group omitted (ethanol_only has no
#: isobutanol coefficients). The rest of the definition: the three
#: feeding variables threshold_conc / target_delta / max_n_spikes on the
#: engine's default bands; the spike feed (spike_delta_bounds=None,
#: 600 g/L) and the aerobic stage-1 cutoff (stage_1_max_x_bounds=None,
#: 5.0 g/L) PINNED at the scenario baseline (k_2, k_4,
#: k_5, k_5e, k_7, k_8, k_9, k_9e, k_9c, k_10, k_11 and every K_* term
#: are simply not in the set -- baseline values, no knockout probe).
METABOLIC_MINIMAL_SUBSET_GROUPS = {
    'inhib_ethanol':    ('k_1ie', 'k_4ie', 'k_7ie', 'k_10ie', 'k_16ie'),
    'inhib_isobutanol': ('k_1ii', 'k_4ii', 'k_6ii', 'k_7ii', 'k_10ii'),
    'inhib_acetate':    ('k_1ia', 'k_4ia', 'k_6ia', 'k_7ia', 'k_10ia', 'k_16ia'),
}

#: Per-study_type options beyond the role filter (a type absent here
#: takes the defaults: DEFAULT_EXCLUDED_PARAMETERS, no groups,
#: DEFAULT_GROUP_MULTIPLIER_BOUNDS, DEFAULT_SPIKE_DELTA_BOUNDS,
#: DEFAULT_STAGE_1_MAX_X_BOUNDS). Keys: exclude_params, group_roles (the
#: roles grouped by effector), group_multiplier_bounds, spike_delta_bounds
#: (None = pinned), stage_1_max_x_bounds (None = pinned; since the
#: metabolic_minimal_subset preset -- honoured for every type), and the
#: EXPLICIT-definition pair rate_params + parameter_groups (present =
#: the set is listed outright and intersected with the workbook;
#: group_roles absent).
#: 'metabolic_minimal' = the compact, interpretable space (24 variables
#: for ethanol_isobutanol, 19 for ethanol_only): exclude_params = k_10
#: (decay, as everywhere) + k_7 and k_8 (the growth capacities, so the
#: burden's phi_T stays at wild type); group_roles = the inhibition
#: coefficients, grouped by the role table's effector into
#: inhib_ethanol / inhib_isobutanol / inhib_acetate (EFFECTOR_ORDER),
#: each ONE log multiplier on group_multiplier_bounds (0.2x-2x); and
#: spike_delta_bounds = None, the spike feed pinned at the scenario
#: baseline (600 g/L). The driver tags the group band as the inhibition
#: band (`_ib0.2-2`) and the exclusion set as `_xk10+k7+k8`
#: (study_type_name_defaults); the group columns and the missing
#: spike_delta column keep the header guard from any cross-resume.
STUDY_TYPE_OPTIONS = {
    'metabolic_minimal': dict(
        exclude_params=('k_10', 'k_7', 'k_8'),
        group_roles=('product_inhibition', 'lethality'),
        group_multiplier_bounds=(0.2, 2.0),
        spike_delta_bounds=None,
    ),
    # The standalone explicit set (METABOLIC_MINIMAL_SUBSET_*): 9 rates +
    # 3 groups + 3 feeding variables = 15 for ethanol_isobutanol (10 for
    # ethanol_only after the workbook intersection); no exclusions, spike
    # and stage_1_max_x pinned; default name
    # kin_opt_ethanol_isobutanol_metabolic_minimal_subset_irr_rb0.001-10_ib0.2-2_burden
    # (no _x / _s1x tag). Its column set differs from every other study's,
    # so the CSV header guard refuses any cross-resume regardless.
    'metabolic_minimal_subset': dict(
        rate_params=METABOLIC_MINIMAL_SUBSET_RATES,
        parameter_groups=METABOLIC_MINIMAL_SUBSET_GROUPS,
        group_multiplier_bounds=(0.2, 2.0),
        exclude_params=(),
        spike_delta_bounds=None,
        stage_1_max_x_bounds=None,
    ),
}
#: Order of the effector groups of a grouped preset (group name
#: `inhib_{effector}`); effectors with no rows in the workbook set are
#: omitted (ethanol_only has no isobutanol coefficients).
EFFECTOR_ORDER = ('ethanol', 'isobutanol', 'acetate')

def study_type_name_defaults(study_type):
    """The per-study_type defaults that ENTER THE STUDY NAME, as
    dict(inhibition_multiplier_bounds=(lo, hi), exclude_params=(names),
    stage_1_max_x_bounds=(lo, hi) or None):
    the STUDY_TYPE_OPTIONS entry's group_multiplier_bounds, exclude_params
    and stage_1_max_x_bounds when the type has them (metabolic_minimal:
    (0.2, 2.0) and ('k_10', 'k_7', 'k_8'); metabolic_minimal_subset:
    (0.2, 2.0), () and None = pinned), else
    DEFAULT_SATURATION_MULTIPLIER_BOUNDS, DEFAULT_EXCLUDED_PARAMETERS and
    DEFAULT_STAGE_1_MAX_X_BOUNDS. resolve_study_preset builds its
    multiplier_bounds / exclude_params / stage_1_max_x_bounds from here
    and the supervisor's default_study_name reads its `_ib` / `_x` /
    `_s1x` tags from here, so the two names can never drift (the
    supervisor's stall watchdog counts rows of the name IT derives).
    Unknown study_type: ValueError."""
    if study_type not in STUDY_TYPE_ROLES:
        raise ValueError(f'Unknown study_type {study_type!r}; expected one '
                         f'of {sorted(STUDY_TYPE_ROLES)}.')
    options = STUDY_TYPE_OPTIONS.get(study_type, {})
    stage_1_max_x_bounds = options.get('stage_1_max_x_bounds',
                                       DEFAULT_STAGE_1_MAX_X_BOUNDS)
    return dict(
        inhibition_multiplier_bounds=tuple(options.get(
            'group_multiplier_bounds', DEFAULT_SATURATION_MULTIPLIER_BOUNDS)),
        exclude_params=tuple(options.get('exclude_params',
                                         DEFAULT_EXCLUDED_PARAMETERS)),
        stage_1_max_x_bounds=(None if stage_1_max_x_bounds is None
                              else tuple(stage_1_max_x_bounds)))

DEFAULT_STUDY_TARGET_PRODUCTS = 'ethanol_isobutanol'
DEFAULT_STUDY_TYPE = 'metabolic_protein'

def resolve_study_preset(study_target_products, study_type, roles=None,
                         effectors=None):
    """The driver run() kwargs of a named study preset (sim-free):
    dict(scenario='A', kinetic_bounds_scenario=<A|B>,
    include_params=[workbook rows of that scenario whose role is in
    STUDY_TYPE_ROLES[study_type], in workbook order],
    multiplier_bounds=DEFAULT_SATURATION_MULTIPLIER_BOUNDS,
    rate_multiplier_bounds=DEFAULT_RATE_MULTIPLIER_BOUNDS,
    rate_params=[the workbook's RATE CONSTANTS -- rate_constant_names,
    role capacity; 16 in A's workbook, 20 in B's -- the only names the
    rate band applies to, so the inhibition coefficients k_*i* sample
    the saturation band like the K_* terms],
    parameter_multiplier_bounds=a COPY of
    DEFAULT_PARAMETER_MULTIPLIER_BOUNDS [k_10, the active-biomass decay
    capacity, on 0.1x-10x instead of the rate band -- in force only when
    k_10 is re-included],
    exclude_params=a COPY of DEFAULT_EXCLUDED_PARAMETERS [('k_10',): the
    decay capacity is NOT a decision variable by default -- it stays at
    the scenario baseline and gets no knockout probe; the driver tags the
    effective exclusion into the study name, `_xk10`]).
    Also stage_1_max_x_bounds=study_type_name_defaults(study_type)
    ['stage_1_max_x_bounds'] [a tuple COPY of DEFAULT_STAGE_1_MAX_X_BOUNDS,
    (1.0, 50.0) g/L, log-scale, for every type but
    metabolic_minimal_subset, which pins it (None): the fermentor's
    aerobic stage-1 biomass cutoff, an OPERATING variable
    (OPERATING_VARIABLES) applied via V406.stage_1_max_x; the driver
    tags a band `_s1x1-50` into the study name, nothing when pinned].
    Set sizes (include_params, the workbook rows): ethanol_only 29
    (metabolic) / 40 (metabolic_protein); ethanol_isobutanol 40 / 56 --
    one fewer each in the sampled space after the exclusion. `roles` (default
    kinetic_parameter_roles()) is the {name: role} table; a workbook row
    absent from it raises KeyError(name) so a future workbook/model change
    can never leak a parameter into a set silently. Unknown axis values
    raise ValueError.

    Since 2026-09-07 every preset also returns `parameter_groups`
    ({group: [members]} for a grouped study_type, None otherwise),
    `group_multiplier_bounds` and `spike_delta_bounds` (a tuple, or None
    = spike pinned at the baseline), and its `multiplier_bounds` /
    `exclude_params` come from study_type_name_defaults(study_type).
    'metabolic_minimal' (STUDY_TYPE_OPTIONS): include_params = the
    capacity + product_inhibition + lethality rows (36 for the B
    workbook, 25 for A's); exclude_params ('k_10', 'k_7', 'k_8'); the
    inhibition rows grouped by the role table's effector (`effectors`,
    default kinetic_parameter_effectors(); a grouped row whose effector
    is None or not in EFFECTOR_ORDER raises KeyError) into
    inhib_ethanol / inhib_isobutanol / inhib_acetate in EFFECTOR_ORDER,
    members in workbook order, effectors without rows omitted; so the
    sampled space is 17 rates + 3 multipliers (+ 4 feeding/operating)
    for ethanol_isobutanol and 13 + 2 (+ 4) for ethanol_only.
    Every group member is scaled from its LIVE baseline -- the model value
    at study start (for a row absent from the scenario-A workbook, the
    nskinetics model default), NOT the workbook value; the engine prints
    each member's baseline next to its name and records the products as
    the applied_<member> columns.

    'metabolic_minimal_subset' (2026-09-07) is the STANDALONE explicit
    preset: STUDY_TYPE_ROLES has no roles for it (the empty tuple) and
    its STUDY_TYPE_OPTIONS entry lists the set outright -- rate_params
    = METABOLIC_MINIMAL_SUBSET_RATES (9), parameter_groups =
    METABOLIC_MINIMAL_SUBSET_GROUPS (3 effector groups), no exclusions,
    spike AND stage_1_max_x pinned. Every listed rate must be a
    capacity row and every member a product_inhibition / lethality row
    of the role table (KeyError naming the parameter and the preset,
    checked BEFORE the intersection so a typo can never be dropped
    silently); then include_params = the listed rates that are rows of
    the target's workbook (list order) and parameter_groups = the
    members that are (empty groups omitted; None if all are), so
    ethanol_isobutanol samples 9 rates + 3 multipliers + 3 feeding
    variables = 15 and ethanol_only (no k_13-k_16, no isobutanol
    coefficients in the A workbook) 5 + 2 + 3 = 10. The returned
    rate_params is still the workbook's every capacity row (the listed
    rates are a subset); `effectors` is not consulted.
    """
    if study_target_products not in STUDY_TARGET_PRODUCTS:
        raise ValueError(
            f'Unknown study_target_products {study_target_products!r}; '
            f'expected one of {sorted(STUDY_TARGET_PRODUCTS)}.')
    if study_type not in STUDY_TYPE_ROLES:
        raise ValueError(f'Unknown study_type {study_type!r}; expected one '
                         f'of {sorted(STUDY_TYPE_ROLES)}.')
    target = STUDY_TARGET_PRODUCTS[study_target_products]
    allowed_roles = set(STUDY_TYPE_ROLES[study_type])
    if roles is None:
        roles = kinetic_parameter_roles()
    set_scenario = target['parameter_set_scenario']
    include_params = []
    workbook_rows = list(workbook_kinetic_baselines(set_scenario))
    for name in workbook_rows:
        if name not in roles:
            raise KeyError(
                f'{name!r} (scenario-{set_scenario} workbook row) has no '
                'entry in the nskinetics kinetic-parameter role table; '
                'refusing to build the study preset.')
        if roles[name] in allowed_roles:
            include_params.append(name)
    options = STUDY_TYPE_OPTIONS.get(study_type, {})
    name_defaults = study_type_name_defaults(study_type)
    group_roles = set(options.get('group_roles', ()))
    parameter_groups = None
    if options.get('rate_params') is not None:
        # EXPLICIT preset (metabolic_minimal_subset): the set is listed
        # outright. (1) Typo guard against the role table FIRST -- a
        # misspelt name must never be dropped silently by the
        # intersection below; (2) the listed rates present in the
        # target's workbook, list order; (3) the groups' members present
        # in it, empty groups omitted (None if every group is).
        explicit_rates = tuple(options['rate_params'])
        explicit_groups = dict(options['parameter_groups'])
        for name in explicit_rates:
            if roles.get(name) not in RATE_CONSTANT_ROLES:
                raise KeyError(
                    f'{name!r} (listed rate constant of the {study_type!r} '
                    'preset) is not a capacity-role row of the nskinetics '
                    'kinetic-parameter role table (role '
                    f'{roles.get(name)!r}); refusing to build the preset.')
        for group, members in explicit_groups.items():
            for name in members:
                if roles.get(name) not in INHIBITION_COEFFICIENT_ROLES:
                    raise KeyError(
                        f'{name!r} (member of group {group!r} of the '
                        f'{study_type!r} preset) is not a product_inhibition '
                        '/ lethality row of the nskinetics kinetic-parameter '
                        f'role table (role {roles.get(name)!r}); refusing '
                        'to build the preset.')
        workbook_set = set(workbook_rows)
        include_params = [name for name in explicit_rates
                          if name in workbook_set]
        parameter_groups = {
            group: [name for name in members if name in workbook_set]
            for group, members in explicit_groups.items()}
        parameter_groups = {group: members
                            for group, members in parameter_groups.items()
                            if members} or None
    elif group_roles:
        if effectors is None:
            effectors = kinetic_parameter_effectors()
        families = {effector: [] for effector in EFFECTOR_ORDER}
        for name in include_params:
            if roles[name] not in group_roles:
                continue
            effector = effectors.get(name)
            if effector is None:
                raise KeyError(
                    f'{name!r} (role {roles[name]!r}, to be grouped by '
                    'effector) has no effector in the nskinetics '
                    'kinetic-parameter table; refusing to build the '
                    f'{study_type!r} preset.')
            if effector not in families:
                raise KeyError(
                    f'{name!r} has effector {effector!r}, not in '
                    f'EFFECTOR_ORDER {EFFECTOR_ORDER}; refusing to build '
                    f'the {study_type!r} preset.')
            families[effector].append(name)
        parameter_groups = {f'inhib_{effector}': members
                            for effector, members in families.items()
                            if members}
    return dict(scenario=target['scenario'],
                kinetic_bounds_scenario=set_scenario,
                include_params=include_params,
                multiplier_bounds=name_defaults['inhibition_multiplier_bounds'],
                rate_multiplier_bounds=DEFAULT_RATE_MULTIPLIER_BOUNDS,
                rate_params=rate_constant_names(workbook_rows, roles=roles),
                parameter_multiplier_bounds=dict(
                    DEFAULT_PARAMETER_MULTIPLIER_BOUNDS),
                exclude_params=name_defaults['exclude_params'],
                stage_1_max_x_bounds=name_defaults['stage_1_max_x_bounds'],
                parameter_groups=parameter_groups,
                group_multiplier_bounds=tuple(options.get(
                    'group_multiplier_bounds', DEFAULT_GROUP_MULTIPLIER_BOUNDS)),
                spike_delta_bounds=options.get('spike_delta_bounds',
                                               DEFAULT_SPIKE_DELTA_BOUNDS))

#: Study-name suffix of a burden-enabled study (enzyme_burden.py): it
#: records extra columns and a different physiology, so it must never
#: share a trajectory CSV or SQLite store with a burden-free study. Both
#: naming paths append it (default_study_name; the driver's and the
#: engine's legacy fallbacks); the CSV header guard is the second line
#: of defence.
BURDEN_STUDY_SUFFIX = '_burden'

def default_study_name(objective, study_target_products, study_type,
                       scenario=None, kinetic_bounds_scenario=None,
                       burden=False, rate_multiplier_bounds=None,
                       inhibition_multiplier_bounds=None,
                       exclude_params=None, stage_1_max_x_bounds=None,
                       n_seeds=None):
    """Stable study name of a preset study:
    kin_opt_{study_target_products}_{study_type}_{objective slug}
    (slug = lower-cased, spaces -> '_'), e.g.
    kin_opt_ethanol_isobutanol_metabolic_protein_irr. New names, so a
    preset study can never collide with the CSV/SQLite of a legacy
    kin_opt_{scenario}[_kb{X}]_{slug} study. The driver and the supervisor
    both derive it from here.

    `scenario` / `kinetic_bounds_scenario` tag the name ONLY when an
    explicit override differs from the preset's own values -- every
    preset starts at scenario 'A' and its kinetic_bounds_scenario is
    STUDY_TARGET_PRODUCTS[study_target_products]['parameter_set_scenario']
    -- so the tags appear only for a genuine override: `_sc{scenario}`
    when scenario is given and != 'A', `_kb{kinetic_bounds_scenario}` when
    it is given and differs from the preset's own workbook scenario.
    Passing back the preset's own values (or leaving both None) reproduces
    the base name above exactly, so default preset names are unchanged
    and still cannot collide with the legacy family. Without this tag, an
    explicit scenario override under a preset (e.g. run(scenario='B'))
    silently resumed the default-scenario study instead (same search-space
    columns, so the CSV header guard could not catch the mix).

    `rate_multiplier_bounds` (the k_* band, (m_lo, m_hi) x baseline) tags
    the name `_rb{m_lo:g}-{m_hi:g}` (e.g. `_rb0.001-10` at the presets'
    own DEFAULT_RATE_MULTIPLIER_BOUNDS, `_rb0.1-10`) WHENEVER it is
    given: a different band samples a different space over the SAME
    columns, so without the tag a run would silently resume the study
    of the same objective under another band (the CSV header guard
    cannot tell them apart). Until the band moved 1e-5x -> 1e-3x later
    on 2026-09-06 only a band DIFFERING from the default was tagged, so
    the 1e-5x studies carry no `_rb` tag (e.g.
    kin_opt_ethanol_isobutanol_metabolic_irr_ib0.1-10_burden); the driver
    and supervisor now pass the EFFECTIVE band on every preset-derived
    name, and only an explicit --study-name / study_name can resume
    those. None (older callers) leaves the name unchanged. The presets'
    per-parameter bands (DEFAULT_PARAMETER_MULTIPLIER_BOUNDS, k_10) are
    part of the preset's identity and are NOT tagged.

    `inhibition_multiplier_bounds` (the band of the inhibition
    coefficients k_*i*, (m_lo, m_hi) x baseline) tags the name
    `_ib{m_lo:g}-{m_hi:g}` (e.g. `_ib0.1-10`) WHENEVER it is given,
    after `_rb`: since 2026-09-06 the presets assign bands by role
    (rate constants alone on the k_* band; inhibition coefficients on
    the saturation band, `multiplier_bounds`), and the driver and
    supervisor pass that band here on every preset-derived name, so a
    role-band study never resumes the CSV/SQLite of a study started
    under the old prefix rule (same columns, so the header guard cannot
    tell them apart). None (older callers) leaves the name unchanged.

    `exclude_params` (the kinetic parameters kept OUT of the search
    space; the presets' DEFAULT_EXCLUDED_PARAMETERS, ('k_10',), since
    2026-09-06 pm) tags the name with excluded_parameters_tag after
    `_ib` (`_xk10`) whenever it is non-empty. An exclusion removes a
    trajectory-CSV column, so the header guard already refuses to resume
    a study of the old set -- the tag is what lets the default name
    launch next to those studies at all (e.g. the 2026-09-06 production
    study kin_opt_ethanol_isobutanol_metabolic_irr_rb0.001-10_ib0.1-10_
    burden, which sampled k_10). None or an empty set (older callers; a
    run that re-includes k_10 with exclude_params=()) leaves the name
    unchanged -- same columns, so resuming those studies is legitimate.

    `stage_1_max_x_bounds` (the (lo, hi) g/L band of the operating
    variable stage_1_max_x; the presets' DEFAULT_STAGE_1_MAX_X_BOUNDS,
    since 2026-09-06 pm) tags the name `_s1x{lo:g}-{hi:g}` (`_s1x1-50`)
    after the exclusion tag whenever it is given. The variable is a new
    decision COLUMN, so the header guard already refuses to resume a
    study without it -- the tag is what lets the default name launch
    next to the existing `..._xk10_burden` studies at all. None (older
    callers; a run that pins the variable with stage_1_max_x_bounds=None)
    leaves the name unchanged -- same columns, so resuming those studies
    is legitimate.

    `n_seeds` (the number of seed points enqueued from donor studies,
    run_kinetic_optimization(seed_from=...); since 2026-09-07) tags the
    name `_seed{n}` (seed_points_tag) after `_s1x` whenever it is > 0.
    Seeds change the trajectory but not the columns, so without the tag
    a seeded run would silently resume the unseeded study of the same
    name; the count alone is tagged, so a different panel of the same
    size needs an explicit study_name. None / 0 leaves the name
    unchanged.

    `burden=True` appends BURDEN_STUDY_SUFFIX ('_burden') after every
    other tag: a burden study (enzyme_burden.py; the driver's default)
    can never resume a burden-free study's CSV/SQLite, or vice versa.
    """
    slug = objective.lower().replace(' ', '_')
    name = f'kin_opt_{study_target_products}_{study_type}_{slug}'
    if scenario is not None and scenario != 'A':
        name += f'_sc{scenario}'
    preset_kinetic_bounds_scenario = STUDY_TARGET_PRODUCTS.get(
        study_target_products, {}).get('parameter_set_scenario')
    if (kinetic_bounds_scenario is not None
            and kinetic_bounds_scenario != preset_kinetic_bounds_scenario):
        name += f'_kb{kinetic_bounds_scenario}'
    if rate_multiplier_bounds is not None:
        lo, hi = rate_multiplier_bounds
        name += f'_rb{lo:g}-{hi:g}'
    if inhibition_multiplier_bounds is not None:
        lo, hi = inhibition_multiplier_bounds
        name += f'_ib{lo:g}-{hi:g}'
    name += excluded_parameters_tag(exclude_params)
    if stage_1_max_x_bounds is not None:
        lo, hi = stage_1_max_x_bounds
        name += f'_s1x{lo:g}-{hi:g}'
    name += seed_points_tag(n_seeds)
    if burden:
        name += BURDEN_STUDY_SUFFIX
    return name

#%% Trajectory recording

def trajectory_columns(search_space, extra_columns=()):
    """Column order of the trajectory CSV for a given search space.
    `extra_columns` (a burden study passes enzyme_burden.BURDEN_COLUMNS)
    go after the tracked metrics and before 'error', so the decision
    columns (between 'state' and 'objective', pca_decision_matrix) are
    unchanged and the header guard of append_trajectory_row separates
    burden from burden-free trajectories."""
    return ['trial_number', 'state', *search_space.keys(), 'objective',
            *TRACKED_METRICS.keys(), *extra_columns, 'error']

def check_trajectory_header(csv_path, columns):
    """Raise ValueError if an existing trajectory CSV at `csv_path` has a
    header other than `columns` (the search space -- or the burden
    column set -- changed since the trajectory was started, so appending
    would silently misalign the rows). A missing or empty file passes.
    The guard of append_trajectory_row; run_kinetic_optimization also
    calls it up front, before the optuna study is opened, so a study
    name colliding with a store of a different column set fails before
    any sidecar or simulation."""
    if not os.path.isfile(csv_path):
        return
    with open(csv_path, newline='') as csvfile:
        existing_header = next(csv.reader(csvfile), None)
    if existing_header is not None and existing_header != list(columns):
        raise ValueError(
            f'Trajectory CSV {csv_path} has a different column set '
            'than the current search space -- the search space '
            'changed since this trajectory was started. Use a new '
            'study_name (or move the old CSV and .db) to start '
            'fresh.')

def append_trajectory_row(csv_path, columns, record):
    """Append one row (dict; missing keys become '') to `csv_path`,
    writing the header first if the file does not exist. An existing
    file's header must match `columns` exactly (check_trajectory_header)
    -- otherwise the search space changed since the trajectory was
    started, and appending would silently misalign the rows. Flushed
    immediately, so a crash/segfault loses at most the in-flight
    trial."""
    # An empty file passes the header check; treat it as absent so the
    # header is written before its first row.
    exists = os.path.isfile(csv_path) and os.path.getsize(csv_path) > 0
    check_trajectory_header(csv_path, columns)
    with open(csv_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns,
                                extrasaction='ignore')
        if not exists:
            writer.writeheader()
        writer.writerow({k: record.get(k, '') for k in columns})
        csvfile.flush()

def load_trajectory(csv_path):
    """Read a trajectory CSV back as a pandas DataFrame (numeric columns
    parsed; blank cells -> NaN)."""
    import pandas as pd
    return pd.read_csv(csv_path)

#%% In-flight sidecar (lost-trial recovery)
# A trial that hangs inside the native integrator is hard-killed by the
# supervisor (analyses/optimize_kinetics_BO_supervised.py), or dies in the
# known CVODE segfault, and can never write its own terminal CSV row. Its
# trial_number is consumed anyway (optuna leaves the trial RUNNING and the
# resume numbers the next one afresh), so the trajectory CSV ends up with a
# gap -- 33 of 500 trials on kin_opt_A_kbB_irr_20260904 (2026-09-04). The
# engine therefore records each trial's decision vector to a per-study
# sidecar JSON the instant it is known (after every suggest_*, before the
# simulation) and deletes it on every in-process outcome; a sidecar still
# present after the child has ended IS the lost trial, and recover_inflight
# appends it as a state='LOST' row -- with its own trial_number, filling
# the gap. stdlib-only, so the supervisor can call it.

def inflight_path_for(results_dir, study_name):
    """Path of the per-study in-flight sidecar (next to the trajectory
    CSV); the one place the name is defined, shared by the engine (writer)
    and the supervisor (reader)."""
    return os.path.join(results_dir, study_name + '_inflight.json')

def write_inflight(inflight_path, columns, record):
    """Atomically write the in-flight trial as {'columns': [...],
    'record': {...}} (tmp file + os.replace: a kill mid-write leaves either
    the previous sidecar or the new one, never a partial file). `columns`
    is the trajectory column order, `record` the trial's
    {'trial_number': ..., <decision variables>} dict."""
    tmp = inflight_path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump({'columns': list(columns), 'record': record}, f)
    os.replace(tmp, inflight_path)

def clear_inflight(inflight_path):
    """Remove the sidecar if present (no error when already gone)."""
    try:
        os.remove(inflight_path)
    except FileNotFoundError:
        pass

def recover_inflight(csv_path, inflight_path, state='LOST', error=''):
    """If a sidecar is present (the child ended without writing that
    trial's terminal row), append it to the trajectory CSV as one row with
    the given `state`/`error` -- decision vector as recorded, objective and
    tracked metrics blank (NaN on read-back) -- clear the sidecar and
    return the recovered trial_number. Returns None when there is nothing
    to recover. A sidecar that does not parse or lacks the expected keys
    is reported, cleared and treated as nothing to recover, so a corrupt
    file never aborts the supervisor."""
    if not os.path.isfile(inflight_path):
        return None
    try:
        with open(inflight_path) as f:
            sidecar = json.load(f)
        columns, record = sidecar['columns'], sidecar['record']
        trial_number = record['trial_number']
    except (OSError, ValueError, KeyError, TypeError) as e:
        print(f'Warning: unreadable in-flight sidecar {inflight_path} '
              f'({repr(e)[:120]}); discarded, nothing recovered.')
        clear_inflight(inflight_path)
        return None
    append_trajectory_row(csv_path, columns,
                          {**record, 'state': state, 'error': error})
    clear_inflight(inflight_path)
    return trial_number

#%% Post-run plots
# All three take the trajectory DataFrame (load_trajectory(csv_path)), so
# plots can be regenerated from a saved CSV without re-running anything.
# matplotlib/pandas are imported lazily (headless Agg per the environment).

def _best_so_far(objective_series, direction):
    """Best-so-far series: cumulative max ('maximize') or min
    ('minimize') of the objective, NaNs carried over."""
    return (objective_series.cummax() if direction == 'maximize'
            else objective_series.cummin())

def _completed(df):
    ok = df[df['state'] == 'COMPLETE'].reset_index(drop=True)
    if ok.empty:
        raise ValueError('No completed trials in the trajectory -- '
                         'nothing to plot.')
    return ok

def plot_optimization_trajectories(df, objective_name, direction,
                                   objective_units='', filename=None):
    """Multipanel trajectory figure: the objective (all completed trials
    as scatter + best-so-far line), then every tracked metric and the
    feeding/operating decision variables -- the applied target, threshold
    and spike sugar concentrations (from _applied_feeding) and the aerobic
    stage-1 biomass cutoff stage_1_max_x -- vs trial number, each panel
    overlaying that quantity's value in the incumbent best-objective-so-far
    configuration (NOT its own best-so-far). Failed/pruned trials are
    omitted; the concentration panels appear whenever the feeding columns
    are present and the stage_1_max_x panel whenever that column is (older
    trajectories without it are unaffected). Returns (fig, axes)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    best_rows = ok.iloc[_best_row_indices(ok, direction)].reset_index(
        drop=True)
    metric_names = [m for m in TRACKED_METRICS if m in ok.columns]
    # Non-objective panels, each (title, sampled series, incumbent series):
    # every tracked metric, then the applied feeding concentrations
    # (target/threshold/spike, derived from the feeding decision variables
    # via _applied_feeding) and the operating variable stage_1_max_x. All
    # share the metric-panel style (sampled scatter + best-so-far overlay).
    panels = [(m, ok[m], best_rows[m]) for m in metric_names]
    if 'threshold_conc' in ok.columns or 'target_conc' in ok.columns:
        t_s, th_s, sp_s = _applied_feeding(ok)
        t_b, th_b, sp_b = _applied_feeding(best_rows)
        panels += [('target_conc (g/L)', t_s, t_b),
                   ('threshold_conc (g/L)', th_s, th_b)]
        if 'spike_delta' in ok.columns or 'spike_conc' in ok.columns:
            panels.append(('spike_conc (g/L)', sp_s, sp_b))  # absent = pinned
    if 'stage_1_max_x' in ok.columns:
        panels.append(('stage_1_max_x (g/L)',
                       ok['stage_1_max_x'], best_rows['stage_1_max_x']))
    n_panels = 1 + len(panels)
    ncols = 4
    nrows = math.ceil(n_panels/ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows),
                             squeeze=False)
    flat = axes.ravel()
    ax = flat[0]
    ax.scatter(ok['trial_number'], ok['objective'], s=8, alpha=0.4,
               color='tab:gray', label='trials')
    ax.plot(ok['trial_number'], _best_so_far(ok['objective'], direction),
            color='tab:red', lw=1.5, label='best so far')
    title = f'Objective: {objective_name}'
    if objective_units:
        title += f' ({objective_units})'
    ax.set_title(title)
    ax.set_xlabel('trial')
    ax.legend(fontsize=8)
    for ax, (title, y_sampled, y_incumbent) in zip(flat[1:], panels):
        ax.scatter(ok['trial_number'], y_sampled, s=8, alpha=0.4,
                   color='tab:blue')
        ax.plot(ok['trial_number'], y_incumbent, color='tab:red', lw=1.5,
                label='at best-so-far objective')
        ax.set_title(title)
        ax.set_xlabel('trial')
        ax.legend(fontsize=7)
    for ax in flat[n_panels:]:
        ax.axis('off')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, axes

def _applied_feeding(rows):
    """(target, threshold, spike) APPLIED concentrations from trajectory
    rows (a DataFrame or a single-row Series), handling both the current
    threshold-anchored parameterization (threshold_conc/target_delta/
    spike_delta) and the legacy target-anchored one (target_conc/
    threshold_delta/spike_conc). A trajectory without spike_delta (the
    spike pinned at the baseline, spike_delta_bounds=None) gets NaN for
    the spike."""
    if 'threshold_conc' in rows:
        threshold = rows['threshold_conc']
        target = np.minimum(TARGET_CONC_MAX,
                            threshold + rows['target_delta'])
        if 'spike_delta' in rows:
            spike = np.minimum(SPIKE_CONC_MAX,
                               np.maximum(SPIKE_CONC_MIN,
                                          target + rows['spike_delta']))
        else:
            # spike pinned at the scenario baseline (spike_delta_bounds=
            # None): not a column, so its value is unknown here -> NaN.
            spike = target*np.nan
    else:
        target = rows['target_conc']
        threshold = np.maximum(0.0, target - rows['threshold_delta'])
        spike = rows['spike_conc']
    return target, threshold, spike

def _best_row_indices(ok, direction):
    """Row index (into `ok`) of the incumbent-best trial as of each
    completed trial."""
    best_idx, cur = [], None
    values = ok['objective'].to_list()
    for i, v in enumerate(values):
        if cur is None or (direction == 'maximize' and v > values[cur]) \
                or (direction == 'minimize' and v < values[cur]):
            cur = i
        best_idx.append(cur)
    return best_idx

def plot_parameter_trajectory(df, kinetic_baselines, direction,
                              filename=None):
    """Best-so-far kinetic CONFIGURATION vs trial number: for each kinetic
    parameter, the incumbent's multiplier (value/baseline, log y) at every
    completed trial, plus a companion panel with the three feeding
    concentrations in absolute units (threshold shown as
    target - threshold_delta). Parameters with nonpositive baselines
    (multiplier undefined; e.g. absolute-bounds overrides of zero-baseline
    params) are skipped here -- they still live in the CSV. A
    parameter-group multiplier column plots directly when the caller
    passes it in `kinetic_baselines` with baseline 1.0 (the driver does).
    Returns (fig, axes)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    best_rows = ok.iloc[_best_row_indices(ok, direction)].reset_index(
        drop=True)
    x = ok['trial_number']
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), sharex=True,
                                   height_ratios=[3, 1])
    for pname, baseline in kinetic_baselines.items():
        if pname not in ok.columns or baseline <= 0:
            continue
        ax1.plot(x, best_rows[pname]/baseline, lw=1, label=pname)
    ax1.set_yscale('log')
    ax1.axhline(1.0, color='k', lw=0.8, ls='--')
    ax1.set_ylabel('best-so-far multiplier vs baseline')
    ax1.legend(fontsize=6, ncol=8, loc='upper center',
               bbox_to_anchor=(0.5, -0.08))
    if 'target_conc' in ok.columns or 'threshold_conc' in ok.columns:
        target, threshold, spike = _applied_feeding(best_rows)
        ax2.plot(x, target, label='target_conc')
        ax2.plot(x, threshold, label='threshold_conc')
        if np.isfinite(np.asarray(spike, dtype=float)).any():  # NaN = pinned
            ax2.plot(x, spike, label='spike_conc')
        ax2.legend(fontsize=7)
    ax2.set_ylabel('g/L')
    ax2.set_xlabel('trial')
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200, bbox_inches='tight')
    return fig, (ax1, ax2)

def plot_best_vs_baseline(df, kinetic_baselines, direction, filename=None):
    """The research-prioritization headline: horizontal bars of the FINAL
    incumbent's kinetic-parameter multipliers (log x, baseline = 1 dashed
    line), sorted by multiplier, with the optimal feeding concentrations
    and objective value in the title. Nonpositive-baseline parameters are
    skipped (see plot_parameter_trajectory). A parameter-group multiplier
    column plots directly when the caller passes it in `kinetic_baselines`
    with baseline 1.0 (the driver does). Returns (fig, ax)."""
    import matplotlib.pyplot as plt
    ok = _completed(df)
    i_best = (ok['objective'].idxmax() if direction == 'maximize'
              else ok['objective'].idxmin())
    best = ok.loc[i_best]
    names = [p for p, b in kinetic_baselines.items()
             if p in ok.columns and b > 0]
    multipliers = np.array([best[p]/kinetic_baselines[p] for p in names])
    order = np.argsort(multipliers)
    fig, ax = plt.subplots(figsize=(7, 0.28*len(names) + 2))
    ax.barh([names[i] for i in order], multipliers[order],
            color='tab:blue')
    ax.set_xscale('log')
    ax.axvline(1.0, color='k', ls='--', lw=0.8)
    ax.set_xlabel('best/baseline multiplier')
    target, threshold, spike = _applied_feeding(best)
    spike_text = (f'{spike:.1f} g/L' if np.isfinite(spike)
                  else 'pinned')  # NaN = spike_delta not sampled
    ax.set_title(f"objective = {best['objective']:.4g} at trial "
                 f"{int(best['trial_number'])}; target = {target:.1f}, "
                 f"threshold = {threshold:.1f}, "
                 f"spike = {spike_text}", fontsize=8)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, ax

#%% PCA projection of the sampled decision space

def pca_decision_matrix(df, log_columns=(), decision_columns=None):
    """PCA of the decision vectors of EVERY recorded trial (COMPLETE, FAIL,
    NAN and LOST rows alike -- each records its full decision vector), from a
    trajectory DataFrame (load_trajectory(csv_path)).

    `decision_columns` defaults to the columns between 'state' and
    'objective' -- exactly the search-space variables, by
    trajectory_columns construction. Columns named in `log_columns` (the
    log-sampled kinetic multiplier bands) are log10-transformed; all are
    then z-scored, zero-variance columns dropped (e.g. pinned variables),
    and the PCA computed by SVD (no sklearn). Component signs are
    stabilized (largest-|loading| entry positive) so successive mid-run
    figures do not flip axes.

    Returns (coords, explained_var_ratio, loadings, kept_columns, valid):
    `coords` (n_valid x n_components) are the trial scores, `loadings`
    (n_components x n_kept) the variable weights, and `valid` a boolean
    mask aligned with df's rows (False where any decision value is
    missing/non-finite after transformation)."""
    if decision_columns is None:
        cols = list(df.columns)
        decision_columns = cols[cols.index('state') + 1:
                                cols.index('objective')]
    # .copy(): pandas may hand back a read-only zero-copy view here, and
    # the log10 transform below writes in place.
    T = df[list(decision_columns)].to_numpy(dtype=float).copy()
    with np.errstate(divide='ignore', invalid='ignore'):
        for j, name in enumerate(decision_columns):
            if name in log_columns:
                T[:, j] = np.log10(T[:, j])
    valid = np.isfinite(T).all(axis=1)
    if valid.sum() < 3:
        raise ValueError('Fewer than 3 trials with complete decision '
                         'vectors -- nothing to project.')
    Tv = T[valid]
    kept, columns_z = [], []
    for j, name in enumerate(decision_columns):
        col = Tv[:, j]
        std = col.std()
        if std > 0.0 and np.isfinite(std):
            columns_z.append((col - col.mean())/std)
            kept.append(name)
    if len(kept) < 2:
        raise ValueError('Fewer than 2 non-degenerate decision variables '
                         '-- nothing to project.')
    Z = np.column_stack(columns_z)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    signs = np.sign(Vt[np.arange(Vt.shape[0]),
                       np.argmax(np.abs(Vt), axis=1)])
    signs[signs == 0] = 1.0
    Vt = Vt*signs[:, None]
    coords = U*S*signs[None, :]
    explained_var_ratio = S**2/(S**2).sum()
    return coords, explained_var_ratio, Vt, kept, valid

def plot_pca_projection(df, direction, log_columns=(),
                        objective_name='objective', objective_units='',
                        filename=None, baseline_trial=None):
    """Four-panel PCA view of the sampled decision space: (1) the PC1 x
    PC2 landscape -- completed trials colored by objective, FAIL/NAN/LOST/INFEASIBLE
    trials as gray/red/brown/purple crosses (LOST = stall-killed or crashed
    before writing its row, recovered from the in-flight sidecar;
    INFEASIBLE = over the enzyme-burden cap, pruned before simulating), the incumbent best-so-far path, the
    enqueued baseline (only when `baseline_trial` is given -- its
    trial_number; None, the default, marks no baseline, since a study run
    with enqueue_baseline=False has no baseline trial and its trial 0 is
    just a sampled draw) and the current best marked; (2) the
    explained-variance scree of the top 10 PCs; (3) the top-|loading|
    variables on PC1/PC2; (4) PC1 and PC2 of every sampled point vs trial
    number (the sampler-contraction diagnostic). Works with zero
    completed trials (landscape shows only pruned points). Returns
    (fig, axes)."""
    import matplotlib.pyplot as plt
    coords, evr, loadings, kept, valid = pca_decision_matrix(
        df, log_columns=log_columns)
    dfv = df[valid].reset_index(drop=True)
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.2, 1.2, 1],
                          height_ratios=[1, 1, 0.7])
    ax_main = fig.add_subplot(gs[0:2, 0:2])
    ax_scree = fig.add_subplot(gs[0, 2])
    ax_load = fig.add_subplot(gs[1:3, 2])
    ax_time = fig.add_subplot(gs[2, 0:2])

    pc1, pc2 = coords[:, 0], coords[:, 1]
    state = dfv['state']
    for st, color, label in (('FAIL', 'tab:gray', 'failed (pruned)'),
                             ('NAN', 'tab:red', 'NaN objective (pruned)'),
                             ('LOST', 'tab:brown', 'lost (stalled/crashed)'),
                             ('INFEASIBLE', 'tab:purple',
                              'infeasible (enzyme burden, pruned)')):
        m = (state == st).to_numpy()
        if m.any():
            ax_main.scatter(pc1[m], pc2[m], marker='x', s=18, alpha=0.35,
                            color=color, label=label)
    ok = (state == 'COMPLETE').to_numpy()
    if ok.any():
        sc = ax_main.scatter(pc1[ok], pc2[ok], c=dfv['objective'][ok],
                             cmap='viridis', s=16, alpha=0.85,
                             label='completed')
        clabel = objective_name + (f' ({objective_units})'
                                   if objective_units else '')
        fig.colorbar(sc, ax=ax_main, label=clabel)
        ok_df = dfv[ok].reset_index(drop=True)
        best_path = np.unique(_best_row_indices(ok_df, direction))
        ax_main.plot(pc1[ok][best_path], pc2[ok][best_path], ls='--',
                     color='tab:red', lw=1.2, marker='D', ms=4,
                     label='best-so-far path')
        i_best = (ok_df['objective'].idxmax() if direction == 'maximize'
                  else ok_df['objective'].idxmin())
        ax_main.scatter(pc1[ok][i_best], pc2[ok][i_best], marker='*',
                        s=260, color='gold', edgecolor='k', zorder=5,
                        label='current best')
    if baseline_trial is not None:
        m0 = (dfv['trial_number'] == baseline_trial).to_numpy()
        if m0.any():
            ax_main.scatter(pc1[m0], pc2[m0], marker='*', s=200, color='k',
                            zorder=5,
                            label=f'baseline (trial {baseline_trial})')
    ax_main.set_xlabel(f'PC1 ({100*evr[0]:.1f}% var)')
    ax_main.set_ylabel(f'PC2 ({100*evr[1]:.1f}% var)')
    ax_main.set_title(f'Sampled decision space ({len(kept)} variables), '
                      f'{len(dfv)} trials')
    ax_main.legend(fontsize=7, loc='best')

    n_show = min(10, len(evr))
    ax_scree.bar(np.arange(1, n_show + 1), evr[:n_show],
                 color='tab:blue', alpha=0.8)
    ax_scree.plot(np.arange(1, n_show + 1), np.cumsum(evr[:n_show]),
                  color='tab:red', marker='.', lw=1, label='cumulative')
    ax_scree.set_xlabel('PC')
    ax_scree.set_title('explained variance', fontsize=9)
    ax_scree.legend(fontsize=7)

    strength = np.hypot(loadings[0], loadings[1])
    top = np.argsort(strength)[::-1][:12][::-1]
    y = np.arange(len(top))
    ax_load.barh(y + 0.2, loadings[0][top], height=0.4,
                 color='tab:blue', label='PC1')
    ax_load.barh(y - 0.2, loadings[1][top], height=0.4,
                 color='tab:orange', label='PC2')
    ax_load.set_yticks(y)
    ax_load.set_yticklabels([kept[i] for i in top], fontsize=7)
    ax_load.axvline(0.0, color='k', lw=0.8)
    ax_load.set_title('top loadings', fontsize=9)
    ax_load.legend(fontsize=7)

    x = dfv['trial_number']
    ax_time.scatter(x, pc1, s=6, alpha=0.4, color='tab:blue',
                    label='PC1')
    ax_time.scatter(x, pc2, s=6, alpha=0.4, color='tab:orange',
                    label='PC2')
    ax_time.set_xlabel('trial')
    ax_time.set_ylabel('PC coordinate')
    ax_time.set_title('sampler contraction', fontsize=9)
    ax_time.legend(fontsize=7)

    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=200)
    return fig, (ax_main, ax_scree, ax_load, ax_time)

#%% Process-level trial-timeout supervision (pure logic)
# A pathological kinetic draw can hang the simulation INSIDE a native
# CVODE/roadrunner integrator call (observed 2026-08-31: one trial at
# 100% of a core for ~2.7 h with no output). No in-process mechanism can
# reliably interrupt that on Windows -- watchdog threads and
# PyThreadState_SetAsyncExc only act at Python bytecode boundaries -- so
# the per-trial wall-clock timeout is implemented one level up: a
# supervisor process (analyses/optimize_kinetics_BO_supervised.py) polls
# the crash-safe trajectory CSV, kills the run when no trial has been
# recorded for the timeout, and relaunches the resumable study (the
# in-flight trial is lost, exactly as in a segfault crash; its RUNNING
# optuna record still counts toward the total budget). The decision
# logic lives here, dependency-free, so the offline test covers it.

class StallGuard:
    """Stall detector over trajectory-CSV row counts. Feed it one
    (rows, now) observation per poll; `update` returns 'stalled' once no
    new row has appeared for `stall_timeout_s` (progress resets the
    clock), else None. Timebase: any monotonic seconds. Call `reset()`
    before each new supervised attempt."""
    def __init__(self, stall_timeout_s=1500.0):
        if stall_timeout_s <= 0:
            raise ValueError('stall_timeout_s must be positive; '
                             f'got {stall_timeout_s!r}')
        self.stall_timeout_s = stall_timeout_s
        self.reset()

    def reset(self):
        self._last_rows = None
        self._last_progress_t = None

    def update(self, rows, now):
        if self._last_rows is None or rows > self._last_rows:
            self._last_rows = rows
            self._last_progress_t = now
            return None
        if now - self._last_progress_t >= self.stall_timeout_s:
            return 'stalled'
        return None

def attempt_outcome(exit_code, rows_before, rows_after,
                    killed_for_stall=False, inflight_lost=False,
                    empty_streak=0, max_empty_attempts=1):
    """Supervisor decision after one attempt exits: 'complete' (clean
    exit 0), 'resume' (crash or stall-kill, but the attempt recorded new
    trials -- relaunch the resumable study), or 'abort'.

    An attempt that added NO rows (rows_after == rows_before) is judged
    by what its child got to do (2026-09-06; before that every empty
    attempt aborted, which is what the legacy defaults still do):
    `inflight_lost=True` means the child's in-flight sidecar existed when
    it ended -- it reloaded, sampled and STARTED simulating a trial, so
    the failure was a pathological first draw (a hung native integrator
    call, or a crash) and a reseeded resume will draw a different point:
    'resume'. `inflight_lost=False` means the child never reached a
    simulation (a stall timeout shorter than the ~18 s reload, a broken
    load): 'abort'. `empty_streak` is the number of consecutive empty
    attempts BEFORE this one; once this one makes it
    `max_empty_attempts` (a positive int; the supervisor's default is 5,
    --max-empty-attempts) the study aborts regardless -- the safety net
    against a systematically hanging state (the 2026-09-06 production
    study aborted at trial 763 under the legacy rule after ONE first-draw
    hang under the new 3-min timeout; simulated trials there take 2 s
    median, < 16 s at the 99th percentile, so such a hang is genuine).
    Progress resets the streak (the supervisor tracks it)."""
    if not isinstance(max_empty_attempts, int) or max_empty_attempts < 1:
        raise ValueError('max_empty_attempts must be a positive int; '
                         f'got {max_empty_attempts!r}')
    if exit_code == 0 and not killed_for_stall:
        return 'complete'
    if rows_after > rows_before:
        return 'resume'
    if empty_streak + 1 >= max_empty_attempts:
        return 'abort'
    return 'resume' if inflight_lost else 'abort'

#%% Feasibility-aware TPE sampling (2026-09-06)
# The enzyme-burden cap (enzyme_burden.BurdenModel.evaluate(values).feasible,
# i.e. Phi_M < F_flex) is a closed-form function of the sampled values, so
# the sampler can reject over-cap proposals before they cost a trial
# number. The helpers below are optuna-free at import (optuna is imported
# inside); the sampler subclass is built lazily by
# _feasible_tpe_sampler_class and instantiated by feasible_tpe_sampler.

def search_space_distributions(search_space):
    """{name: optuna distribution} for an engine search space
    ({name: {'low', 'high', 'log'[, 'int']}}, see build_search_space):
    FloatDistribution(low, high, log=) for float entries and
    IntDistribution(low, high) for 'int': True entries -- exactly what the
    objective's suggest_float(..., log=) / suggest_int(...) calls record,
    so optuna accepts the sampler's values as relative parameters."""
    from optuna.distributions import FloatDistribution, IntDistribution
    distributions = {}
    for name, sp in search_space.items():
        if sp.get('int'):
            distributions[name] = IntDistribution(int(sp['low']), int(sp['high']))
        else:
            distributions[name] = FloatDistribution(float(sp['low']),
                                                    float(sp['high']),
                                                    log=bool(sp['log']))
    return distributions

def _uniform_internal_draw(rng, dist):
    """One value of `dist` (an optuna Float/IntDistribution), uniform in
    optuna's internal representation: log-uniform for log floats,
    uniform for linear floats, integer-uniform for ints. Returned as the
    internal-repr float (ints are integer-valued floats)."""
    from optuna.distributions import IntDistribution
    if isinstance(dist, IntDistribution):
        return float(rng.randint(dist.low, dist.high + 1))
    if dist.log:
        x = math.exp(rng.uniform(math.log(dist.low), math.log(dist.high)))
        return float(min(dist.high, max(dist.low, x)))   # exp() round-off
    return float(rng.uniform(dist.low, dist.high))

def _to_external(internal, distributions):
    """{name: external value} of an internal-repr point; a value outside
    its distribution is a programming error (RuntimeError)."""
    values = {}
    for name, dist in distributions.items():
        x = internal[name]
        if not dist._contains(x):
            raise RuntimeError(f'sampled value {x!r} for {name} is outside '
                               f'its distribution {dist}')
        values[name] = dist.to_external_repr(x)
    return values

def draw_uniform_feasible(rng, distributions, is_feasible, max_draws=10_000):
    """Joint uniform draw over `distributions` (search_space_distributions),
    redrawn until `is_feasible(values)` (values in EXTERNAL repr: floats,
    ints for IntDistribution) is true. Returns (values, n_draws,
    feasible): the first feasible draw, or after `max_draws` draws the
    LAST draw with feasible=False. Uniform in optuna's internal repr, so
    the accepted points are exactly uniform on the feasible set under the
    measure the study already uses for its random start-up."""
    max_draws = int(max_draws)
    if max_draws < 1:
        raise ValueError(f'max_draws must be >= 1; got {max_draws!r}')
    values = None
    for n_draws in range(1, max_draws + 1):
        internal = {name: _uniform_internal_draw(rng, dist)
                    for name, dist in distributions.items()}
        values = _to_external(internal, distributions)
        if is_feasible(values):
            return values, n_draws, True
    return values, max_draws, False

def feasible_candidate_mask(samples, distributions, is_feasible):
    """Boolean array over a Parzen candidate batch (`samples` = {name:
    ndarray} in internal repr, as returned by the estimator's sample()):
    candidate i is converted to external values with to_external_repr
    and passed to `is_feasible`."""
    names = list(samples)
    size = len(samples[names[0]])
    mask = np.zeros(size, dtype=bool)
    for i in range(size):
        values = {name: distributions[name].to_external_repr(float(samples[name][i]))
                  for name in names}
        mask[i] = bool(is_feasible(values))
    return mask

class _FeasibleParzenEstimator:
    """Wraps TPE's 'below' Parzen estimator so its candidate batch holds
    only feasible points: sample(rng, size) draws batches of `size` from
    the wrapped estimator, keeps the candidates `is_feasible` accepts, and
    redraws until `size` feasible candidates are collected or
    `sampler.max_parzen_batches` batches have been drawn; then it falls
    back to ONE uniform-feasible draw (draw_uniform_feasible, at most
    `sampler.max_uniform_draws` draws; returned as a size-1 batch in
    internal repr), and if that fails too, returns the last raw batch so
    the engine's in-objective INFEASIBLE guard prunes the trial instead
    of the study crashing. log_pdf delegates unchanged, so optuna's
    expected-improvement scoring and _compare pick among feasible
    candidates only. Counters on `sampler`: n_rejected (candidates and
    uniform draws rejected), n_uniform_fallbacks (batches exhausted),
    n_unfiltered (raw batch returned). Relies on optuna 4.9.0's
    TPESampler._sample using only sample() and log_pdf() of the object
    _build_parzen_estimator returns."""

    def __init__(self, mpe, distributions, is_feasible, sampler):
        self._mpe = mpe
        self._distributions = distributions
        self._is_feasible = is_feasible
        self._sampler = sampler

    def log_pdf(self, samples_dict):
        return self._mpe.log_pdf(samples_dict)

    def sample(self, rng, size):
        kept = {name: [] for name in self._distributions}
        n_kept, last = 0, None
        for _ in range(self._sampler.max_parzen_batches):
            batch = self._mpe.sample(rng, size)
            last = batch
            mask = feasible_candidate_mask(batch, self._distributions,
                                           self._is_feasible)
            self._sampler.n_rejected += int((~mask).sum())
            for name in kept:
                kept[name].append(np.asarray(batch[name])[mask])
            n_kept += int(mask.sum())
            if n_kept >= size:
                break
        if n_kept:
            return {name: np.concatenate(parts)[:size]
                    for name, parts in kept.items()}
        self._sampler.n_uniform_fallbacks += 1
        values, n_draws, feasible = draw_uniform_feasible(
            rng, self._distributions, self._is_feasible,
            self._sampler.max_uniform_draws)
        self._sampler.n_rejected += n_draws - (1 if feasible else 0)
        if feasible:
            return {name: np.array([dist.to_internal_repr(values[name])])
                    for name, dist in self._distributions.items()}
        self._sampler.n_unfiltered += 1
        return last

_FEASIBLE_TPE_CLASS = {}

def _feasible_tpe_sampler_class():
    """The FeasibleTPESampler class (optuna imported here; memoized)."""
    if 'cls' in _FEASIBLE_TPE_CLASS:
        return _FEASIBLE_TPE_CLASS['cls']
    import optuna
    from optuna.trial import TrialState

    class FeasibleTPESampler(optuna.samplers.TPESampler):
        """TPESampler that never proposes a point `is_feasible` rejects
        (values in external repr, the objective's suggest_* values).
        Start-up phase (fewer finished COMPLETE+PRUNED trials than
        n_startup_trials): a JOINT uniform-feasible draw over the full
        engine search space (draw_uniform_feasible; at most
        max_uniform_draws draws, the last raw draw otherwise). TPE phase:
        the base _sample with the 'below' estimator wrapped in
        _FeasibleParzenEstimator (candidates filtered; at most
        max_parzen_batches batches). infer_relative_search_space returns
        the full search space, so every parameter is sampled jointly from
        the first sampled trial on (optuna's intersection space is empty
        until a trial with every parameter completes). Enqueued trials
        never reach the sampler. Counters: n_rejected,
        n_uniform_fallbacks, n_unfiltered (see _FeasibleParzenEstimator).
        group=True and constant_liar=True are refused (both alter
        sample_relative's control flow)."""

        def __init__(self, search_space, is_feasible, *,
                     max_uniform_draws=10_000, max_parzen_batches=20,
                     **tpe_kwargs):
            if tpe_kwargs.get('group'):
                raise ValueError('FeasibleTPESampler does not support group=True')
            if tpe_kwargs.get('constant_liar'):
                raise ValueError('FeasibleTPESampler does not support '
                                 'constant_liar=True')
            if not callable(is_feasible):
                raise TypeError('is_feasible must be callable(values) -> bool; '
                                f'got {is_feasible!r}')
            if int(max_uniform_draws) < 1:
                raise ValueError('max_uniform_draws must be >= 1; got '
                                 f'{max_uniform_draws!r}')
            if int(max_parzen_batches) < 1:
                raise ValueError('max_parzen_batches must be >= 1; got '
                                 f'{max_parzen_batches!r}')
            tpe_kwargs.setdefault('multivariate', True)
            super().__init__(**tpe_kwargs)
            self._feasible_distributions = search_space_distributions(search_space)
            self._is_feasible = is_feasible
            self.max_uniform_draws = int(max_uniform_draws)
            self.max_parzen_batches = int(max_parzen_batches)
            self.n_rejected = 0
            self.n_uniform_fallbacks = 0
            self.n_unfiltered = 0

        def infer_relative_search_space(self, study, trial):
            return {name: dist
                    for name, dist in self._feasible_distributions.items()
                    if not dist.single()}

        def _sample_relative(self, study, trial, search_space):
            if search_space == {}:
                return {}
            states = (TrialState.COMPLETE, TrialState.PRUNED)
            trials = study._get_trials(deepcopy=False, states=states,
                                       use_cache=True)
            if len(trials) < self._n_startup_trials:
                values, n_draws, feasible = draw_uniform_feasible(
                    self._rng.rng, search_space, self._is_feasible,
                    self.max_uniform_draws)
                self.n_rejected += n_draws - (1 if feasible else 0)
                if not feasible:
                    self.n_unfiltered += 1
                return values
            return self._sample(study, trial, search_space)

        def _build_parzen_estimator(self, study, search_space, trials,
                                    handle_below):
            mpe = super()._build_parzen_estimator(study, search_space,
                                                  trials, handle_below)
            if handle_below:
                return _FeasibleParzenEstimator(mpe, search_space,
                                                self._is_feasible, self)
            return mpe

    _FEASIBLE_TPE_CLASS['cls'] = FeasibleTPESampler
    return FeasibleTPESampler

def feasible_tpe_sampler(search_space, is_feasible, **kwargs):
    """A FeasibleTPESampler over the engine `search_space`
    (build_search_space format) with the predicate `is_feasible(values)
    -> bool` (external values). `kwargs`: max_uniform_draws (10_000),
    max_parzen_batches (20), and any optuna TPESampler keyword (seed,
    n_startup_trials, constraints_func, ...; multivariate defaults to
    True; group / constant_liar refused)."""
    return _feasible_tpe_sampler_class()(search_space, is_feasible, **kwargs)

#%% Engine

def get_handles():
    """Resolve the live flowsheet objects the getters and engine need.
    Requires biorefineries.isobutanol.load() to have run in this kernel
    (the system module raises an informative error otherwise)."""
    from biorefineries.isobutanol import system as ibo_system
    f = ibo_system.f
    V406 = f.V406
    return {'f': f,
            'V406': V406,
            'r_te': V406.nsk_kinetic_model._te,
            'fbs_spec': V406.fbs_spec,
            'tea': ibo_system.corn_EtOH_IBO_sys_tea,
            'HXN': f.HXN1001,
            'model_specification': ibo_system.model_specification,
            'solve_TEA': ibo_system.solve_TEA,
            'latest_TEA_solution': {
                'IRR': np.nan,
                'MPSPs': {'ethanol': np.nan, 'isobutanol': np.nan}},
            }

def restore_baseline(handles, kinetic_baselines, baseline_model_kwargs,
                     baseline_max_n_spikes=None, baseline_stage_1_max_x=None):
    """Reset every kinetic parameter to its recorded baseline (and, when
    `baseline_max_n_spikes` is given, the glucose-spike cap
    fbs_spec.max_n_spikes; when `baseline_stage_1_max_x` is given, the
    fermentor's aerobic stage-1 biomass cutoff through the
    V406.stage_1_max_x property, which mirrors it onto the kinetic model
    and the aeration spec), re-simulate at the baseline feeding
    specifications, and refresh the TEA solution -- leaving the process in
    a clean scenario-baseline state. Called in run_kinetic_optimization's
    `finally` (success, exception, or KeyboardInterrupt alike)."""
    r_te = handles['r_te']
    for pname, baseline in kinetic_baselines.items():
        setattr(r_te, pname, baseline)
    if baseline_max_n_spikes is not None:
        handles['fbs_spec'].max_n_spikes = baseline_max_n_spikes
    if baseline_stage_1_max_x is not None:
        handles['V406'].stage_1_max_x = baseline_stage_1_max_x
    handles['model_specification'](**baseline_model_kwargs)
    handles['latest_TEA_solution'].update(
        handles['solve_TEA'](stream_IDs=('ethanol', 'isobutanol')))
    print('Restored kinetic parameters and feeding specifications to '
          'baseline. Baseline TEA solution: '
          f"{handles['latest_TEA_solution']}")

def run_kinetic_optimization(objective='IRR',
                             direction=None, level=None,
                             objective_units=None, objective_name=None,
                             scenario_label='B',
                             n_trials=2000, seed=3221,
                             multiplier_bounds=(0.1, 10.0),
                             param_bounds_override=None,
                             exclude_params=(),
                             include_params=None,
                             rate_multiplier_bounds=None,
                             rate_params=None,
                             parameter_multiplier_bounds=None,
                             threshold_conc_bounds=(0.0, 300.0),
                             target_delta_bounds=(5.0, 500.0),
                             spike_delta_bounds=DEFAULT_SPIKE_DELTA_BOUNDS,
                             max_n_spikes_bounds=(0, 50),
                             stage_1_max_x_bounds=None,
                             target_conc_bounds=None,
                             threshold_delta_bounds=None,
                             spike_conc_bounds=None,
                             study_name=None, results_dir=None,
                             handles=None, print_status_every=1,
                             burden_model='auto',
                             enqueue_baseline=False,
                             enqueue_knockouts=False,
                             n_startup_trials=None,
                             feasible_sampling=True,
                             seed_from=None,
                             parameter_groups=None,
                             group_multiplier_bounds=DEFAULT_GROUP_MULTIPLIER_BOUNDS,
                             ):
    """Run the Bayesian optimization. `objective` is a name in
    OBJECTIVE_REGISTRY (direction/level/units filled from the entry) or a
    custom callable(handles)->float (then `direction`, and ideally
    `objective_name`/`level`/`objective_units`, must be given).

    `n_trials` is the TOTAL budget of the study: rerunning with the same
    study_name resumes from the on-disk SQLite store and runs only the
    remainder (crash/segfault recovery). A FRESH study enqueues NO
    baseline point by default (`enqueue_baseline=False`, the default since
    2026-09-07): the sampler draws every trial from trial 0 -- an enqueued
    scenario-A baseline anchored TPE in the pure-ethanol basin in every
    earlier IRR study. `enqueue_baseline=True` evaluates the scenario
    baseline configuration as trial 0 instead (see
    baseline_decision_point; baseline_point is computed either way, since
    the probes are derived from it). Then -- `enqueue_knockouts=True`
    (default False since 2026-09-07, so by default NO point at all is
    enqueued) -- the single-knockout probes of knockout_probe_points
    (one log-scale rate constant k_* at its band floor, all else at the
    baseline; a rate already at its floor gets none), enqueued right
    after the baseline (if any), each
    tagged with the optuna user attr 'knockout_probe' = its parameter
    name; optuna stores enqueued trials as WAITING, so a resume finishes
    any probes a crash interrupted without re-enqueueing (resumes never
    re-enqueue anything). The probes count toward TPE's random
    start-up phase, so TPE guidance starts as soon as the lethality map
    is in.

    `n_startup_trials` (default None) is the length of that random
    start-up phase (optuna TPESampler n_startup_trials: trials drawn
    uniformly at random, enqueued trials included, before TPE's density
    guidance begins). None applies the rule the engine always used,
    max(10, n_trials//10) -- 200 for a 2000-trial study, which in the
    burden-constrained preset space completed 0 of 183 random draws on
    2026-09-06 -- so an explicit value (e.g. 20-30) is the way to shorten
    it; a
    non-negative integer, ValueError otherwise. It is compared with the
    number of trials already stored, so a resumed study past the
    start-up count starts in TPE mode at once, and a resume may change
    it freely: it is not part of the study's identity (no study-name
    tag, same columns). Trials execute STRICTLY
    sequentially (n_jobs=1; one simulation in flight at a time). Every
    trial appends one row to the trajectory CSV (same stable name as the
    study, '_trajectory.csv' suffix) whether it completes, fails
    (state='FAIL', pruned), yields a NaN objective (state='NAN',
    pruned), or exceeds the enzyme-burden cap (state='INFEASIBLE',
    pruned before simulating). A trial the PROCESS never finishes (hard-killed by the
    supervisor on a stall, or a native segfault) cannot write its row;
    its decision vector is kept in a per-study sidecar
    ('_inflight.json', written before the simulation starts and removed
    once the trial's terminal row has been written; an interrupted trial
    that wrote no row -- e.g. a KeyboardInterrupt during the simulation --
    leaves it for recovery) and appended as a state='LOST' row --
    with its own trial_number, so the CSV has no gaps -- by the
    supervisor as soon as the child ends, or here at the next start of
    the study (recover_inflight). Kinetic parameters and feeding specs
    are restored to the scenario baseline in a `finally`.

    `include_params` (None = every k_*/K_* on the model) restricts the
    kinetic decision variables to the named parameters (see
    build_search_space); the driver passes the scenario workbook's rows
    (kinetic_param_names_from_scenario). Names not on the model are
    ignored with a printed warning. A restricted study has a different
    search space than an unrestricted one of the same study_name -- the
    trajectory-CSV header guard raises rather than misaligning columns,
    so use a fresh study_name.

    `rate_multiplier_bounds` (None = single band) is the separate k_*
    band of build_search_space; the study presets pass
    DEFAULT_RATE_MULTIPLIER_BOUNDS (their absolute bands actually arrive
    via param_bounds_override, see workbook_kinetic_bounds).
    `rate_params` (None = every lowercase 'k_' name, the pre-2026-09-06
    rule) names the RATE CONSTANTS that band -- and the single-knockout
    probes -- apply to; the presets pass their workbook's capacity-role
    rows (rate_constant_names), leaving the inhibition coefficients
    k_*i* on `multiplier_bounds` with no probe.
    `parameter_multiplier_bounds` ({name: (m_lo, m_hi)}, None = none) is
    the per-parameter band of build_search_space, taking precedence over
    the role band; the presets pass DEFAULT_PARAMETER_MULTIPLIER_BOUNDS
    (k_10 on 0.1x-10x; its probe is then a knock-down at 0.1x) -- in
    force only when k_10 is in the space: the presets also pass
    `exclude_params` = DEFAULT_EXCLUDED_PARAMETERS (('k_10',), since
    2026-09-06 pm), so by default k_10 is not sampled and has no probe.

    `stage_1_max_x_bounds` (None = not sampled; the driver passes the
    preset's DEFAULT_STAGE_1_MAX_X_BOUNDS, (1.0, 50.0) g/L) adds the
    OPERATING variable stage_1_max_x -- the fermentor's aerobic stage-1
    biomass cutoff, log-scale -- to the space (build_search_space). Its
    baseline is read off handles['V406'].stage_1_max_x at study start
    (5.0 g/L at the nskinetics factory default) and enters trial 0; each
    trial sets handles['V406'].stage_1_max_x (the property mirrors onto
    the kinetic model and the AerationSpec) right before
    model_specification, and restore_baseline puts it back in the
    `finally`. Not a rate constant: no knockout probe; not a kinetic
    name: ignored by the burden model. A new decision column, so a study
    with it can never resume one without it (header guard); the driver
    tags `_s1x{lo}-{hi}` into the study name.

    `burden_model` (default 'auto') is the enzyme-burden (proteome-
    allocation) constraint of enzyme_burden.py: 'auto' builds
    BurdenModel.from_reference(kinetic_baselines) -- the live values
    after the scenario workbook load, so the burden is exactly inert at
    trial 0 -- None disables it (the legacy burden-free study; older
    studies), and a BurdenModel instance is used as given (the driver
    builds one to print its reports first). With the burden on, every
    trial's sampled vector is evaluated BEFORE the sidecar and the
    simulation: its burden quantities (enzyme_burden.BURDEN_COLUMNS) are
    recorded in the CSV, `burden_violation` (= Phi_M - F_flex, <= 0
    feasible) is set as a trial user attr and fed to the TPE sampler's
    constraints_func, an over-cap point is logged as state='INFEASIBLE'
    and pruned without simulating, and a feasible point reaches the
    model with the DERATED growth capacities k_7_eff/k_8_eff while the
    CSV keeps the sampled k_7/k_8 as the decision. The default study
    name gains BURDEN_STUDY_SUFFIX; restore_baseline is unchanged
    (sampled-space baselines are written back).

    `feasible_sampling` (default True; since 2026-09-06) makes the
    sampler feasibility-aware when the burden is on: feasible_tpe_sampler
    binds FeasibleTPESampler to `burden_model.evaluate(values).feasible`
    (exactly the cap Phi_M < F_flex), so the random start-up draws are
    joint uniform-feasible vectors and every TPE candidate is filtered
    before it is scored -- no sampled trial should ever be INFEASIBLE
    (the 2026-09-06 production study proposed 118 of 200 start-up and
    261 of 1060 TPE trials over the cap). The in-objective INFEASIBLE
    guard and the constraints_func stay as the safety net (a resumed
    study already holds INFEASIBLE rows), and the counters are printed
    at the end of the run (`unfiltered draws` > 0 means the feasible
    region was too small for max_uniform_draws to find; expected 0).
    False (or burden_model=None, where there is no predicate) keeps the
    plain TPESampler. A sampler setting like n_startup_trials: no
    column or study-name change, free to change on a resume.

    `seed_from` (default None; since 2026-09-07) is a sequence of
    (donor, trial_numbers) pairs -- `donor` a trajectory-CSV path or a
    study name resolved to '{results_dir}/{donor}_trajectory.csv' -- whose
    decision points a FRESH study enqueues right after the knockout
    probes (seed_points_from_trajectory: same columns required, values
    clipped into this study's bounds; optuna user attr 'seed' = the
    '{donor study}#{trial}' label; stored WAITING like the probes, so a
    crash mid-seeds resumes them and a resume never re-enqueues). Seeds
    give TPE a foothold in a basin another study found -- e.g. a
    high-isobutanol cell from an IBO-objective study in an IRR study
    whose own sampling never reached that basin (a local optimum) --
    and count toward the random start-up phase. The seed points are
    resolved BEFORE the optuna store is opened, so a bad donor fails
    before any simulation. Same columns as the unseeded study: the
    driver tags `_seed{n}` into the study name (seed_points_tag).

    `parameter_groups` / `group_multiplier_bounds` (since 2026-09-07;
    None = none) are build_search_space's PARAMETER GROUPS: each group
    is one log-scale multiplier decision variable whose members (removed
    from the individual space) receive baseline x multiplier
    (expand_grouped_values) right after sampling -- what the burden
    model evaluates and what is set on the kinetic model. The CSV keeps
    the group multiplier as the DECISION column and adds one derived
    `applied_<member>` column per member after the burden columns
    (before 'error'), so pca_decision_matrix still sees exactly the
    search space and a seed reader can recover the members. Trial 0 has
    every multiplier at 1.0; a group is no rate constant, so it gets no
    knockout probe. `spike_delta_bounds=None` pins the spike
    concentration at the scenario-baseline snapshot (fbs_spec.spike_conc
    at study start; no spike_delta column). The metabolic_minimal preset
    passes all three (resolve_study_preset).

    Returns (study, csv_path, kinetic_baselines)."""
    import optuna
    if handles is None:
        handles = get_handles()
    r_te, fbs_spec = handles['r_te'], handles['fbs_spec']

    if isinstance(objective, str):
        entry = OBJECTIVE_REGISTRY[objective]
        objective_getter = entry['getter']
        objective_name = objective_name or objective
        direction = direction or entry['direction']
        level = level or entry['level']
        if objective_units is None:
            objective_units = entry['units']
    else:
        objective_getter = objective
        objective_name = objective_name or 'custom'
        if direction not in ('maximize', 'minimize'):
            raise ValueError("A custom objective callable requires "
                             "direction='maximize' or 'minimize'.")

    kinetic_baselines = discover_kinetic_parameters(r_te)
    if isinstance(burden_model, str):
        if burden_model != 'auto':
            raise ValueError("burden_model must be 'auto', None or a "
                             f"BurdenModel; got {burden_model!r}.")
        from biorefineries.isobutanol.enzyme_burden import BurdenModel
        burden_model = BurdenModel.from_reference(kinetic_baselines)
    elif burden_model is not None:
        if not all(hasattr(burden_model, a)
                   for a in ('evaluate', 'apply', 'reference')):
            raise TypeError("burden_model must be 'auto', None or a "
                            'BurdenModel (an object with evaluate, apply '
                            f'and reference); got {burden_model!r}.')
        # A caller-built model must be a snapshot of the LIVE baselines
        # (the values on the model now, after the scenario workbook load),
        # or the ratio route is not inert at trial 0 and every native
        # step is mis-charged. Exact equality: both sides are
        # float(getattr(r_te, name)).
        stale = [(name, ref, kinetic_baselines.get(name))
                 for name, ref in burden_model.reference.items()
                 if kinetic_baselines.get(name) != ref]
        if stale:
            raise ValueError(
                'burden_model.reference differs from the live kinetic '
                'baselines of the model for '
                + ', '.join(f'{name} (model reference {ref!r}, live '
                            f'baseline {live!r})'
                            for name, ref, live in stale)
                + '; build the BurdenModel from the same kinetic_baselines '
                "(discover_kinetic_parameters(r_te)) or pass 'auto'.")
    burden_on = burden_model is not None
    if burden_on:
        from biorefineries.isobutanol.enzyme_burden import BURDEN_COLUMNS
        extra_columns = BURDEN_COLUMNS
        print('Enzyme burden ON (enzyme_burden.py): F_flex = '
              f'{burden_model.F_flex:.4f}, Phi_M,wt = {burden_model.Phi_M_wt:.4f}, '
              f'phi_T,wt = {burden_model.phi_T_wt:.4f} g/gDCW; over-cap trials '
              'are logged INFEASIBLE and pruned before simulating.')
    else:
        extra_columns = ()
        print('Enzyme burden OFF (burden_model=None): legacy burden-free study.')
    if include_params is not None:
        missing = [p for p in include_params if p not in kinetic_baselines]
        if missing:
            print(f'Warning: {len(missing)} include_params names are not '
                  f'kinetic parameters of the model and are ignored: '
                  f'{missing}')
    search_space, excluded = build_search_space(
        kinetic_baselines,
        multiplier_bounds=multiplier_bounds,
        param_bounds_override=param_bounds_override,
        exclude_params=exclude_params,
        include_params=include_params,
        rate_multiplier_bounds=rate_multiplier_bounds,
        rate_params=rate_params,
        parameter_multiplier_bounds=parameter_multiplier_bounds,
        threshold_conc_bounds=threshold_conc_bounds,
        target_delta_bounds=target_delta_bounds,
        spike_delta_bounds=spike_delta_bounds,
        max_n_spikes_bounds=max_n_spikes_bounds,
        target_conc_bounds=target_conc_bounds,
        threshold_delta_bounds=threshold_delta_bounds,
        spike_conc_bounds=spike_conc_bounds,
        stage_1_max_x_bounds=stage_1_max_x_bounds,
        parameter_groups=parameter_groups,
        group_multiplier_bounds=group_multiplier_bounds)
    parameter_groups = {str(group): list(members)
                        for group, members in dict(parameter_groups or {}).items()}
    applied_columns = [f'applied_{member}'
                       for members in parameter_groups.values()
                       for member in members]
    n_kinetic = sum(1 for name in search_space if name in kinetic_baselines)
    n_group = sum(1 for name in search_space if name in parameter_groups)
    n_operating = sum(1 for name in search_space if name in OPERATING_VARIABLES)
    n_feeding = len(search_space) - n_kinetic - n_group - n_operating
    restriction = ('' if include_params is None else
                   f', restricted to {len(include_params)} named parameters')
    print(f'Search space: {len(search_space)} decision variables '
          f'({n_kinetic} kinetic + {n_group} group multipliers + '
          f'{n_feeding} feeding + {n_operating} operating{restriction}); '
          f'{len(excluded)} kinetic parameters excluded: {excluded}')
    for group, members in parameter_groups.items():
        sp = search_space[group]
        # Each member's LIVE baseline is printed next to its name: the
        # multiplier is applied to it, so this line is the log's record of
        # the basis of every applied_<member> column.
        listed = ', '.join(f'{member} ({kinetic_baselines[member]:g})'
                           for member in members)
        print(f"Parameter group {group}: one log-scale multiplier on "
              f"[{sp['low']:g}, {sp['high']:g}] x baseline applied to "
              f"{len(members)} members {listed} (baseline 1.0; recorded as "
              f"applied_<member> columns).")
    if 'stage_1_max_x' in search_space:
        sp = search_space['stage_1_max_x']
        print(f"Operating variable stage_1_max_x (aerobic stage-1 biomass "
              f"cutoff) sampled log-scale on [{sp['low']:g}, {sp['high']:g}] "
              f"g/L via V406.stage_1_max_x (baseline "
              f"{handles['V406'].stage_1_max_x:g} g/L).")
    if rate_multiplier_bounds is not None:
        _is_rate = _rate_predicate(rate_params)
        n_rate = sum(1 for name in search_space
                     if name in kinetic_baselines and _is_rate(name))
        rule = ('by role (rate_params)' if rate_params is not None
                else "by the lowercase 'k_' prefix (legacy rule)")
        print(f'Rate band {tuple(rate_multiplier_bounds)} x baseline on '
              f'{n_rate} rate constants {rule}; the other '
              f'{n_kinetic - n_rate} kinetic parameters (inhibition '
              f'coefficients, K_* terms) on {tuple(multiplier_bounds)} x '
              'baseline (where no override applies).')
    if parameter_multiplier_bounds:
        applied = {name: tuple(band)
                   for name, band in parameter_multiplier_bounds.items()
                   if name in search_space and name in kinetic_baselines
                   and name not in (param_bounds_override or {})}
        print(f'Per-parameter bands (x baseline) in force: {applied}')

    # Scenario baseline snapshot for restoration (the driver has already
    # baseline-simulated the scenario, so current_specifications IS the
    # scenario baseline).
    baseline_model_kwargs = {
        k: fbs_spec.current_specifications[k]
        for k in ('target_conc', 'threshold_conc', 'spike_conc')}
    baseline_max_n_spikes = fbs_spec.max_n_spikes
    if 'threshold_conc' in search_space and 'spike_delta' not in search_space:
        # Read the snapshot, not the live spec: the snapshot is what
        # _objective actually applies as the pinned spike concentration.
        print('Spike feed pinned at the scenario baseline '
              f"({baseline_model_kwargs['spike_conc']:g} g/L; "
              'spike_delta_bounds=None).')
    # The live cutoff IS the scenario baseline (never set by the build);
    # read only when the variable is sampled, so handles without the
    # attribute (older callers, offline fakes) keep working.
    baseline_stage_1_max_x = (float(handles['V406'].stage_1_max_x)
                              if 'stage_1_max_x' in search_space else None)

    if results_dir is None:
        results_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'analyses', 'results')
    os.makedirs(results_dir, exist_ok=True)
    slug = objective_name.lower().replace(' ', '_')
    if study_name is None:
        study_name = f'kin_opt_{scenario_label}_{slug}'
        if burden_on:
            study_name += BURDEN_STUDY_SUFFIX
    csv_path = os.path.join(results_dir, study_name + '_trajectory.csv')
    inflight_path = inflight_path_for(results_dir, study_name)
    storage = ('sqlite:///'
               + os.path.join(results_dir, study_name + '.db')
               .replace('\\', '/'))
    columns = trajectory_columns(search_space,
                                 extra_columns=[*extra_columns, *applied_columns])
    # Seed points from donor studies: resolved now (sim-free), so a donor
    # of another column set / a missing trial fails before the store and
    # the first ~20 s simulation; enqueued below on a fresh study only.
    seed_points, seed_notes = {}, []
    for donor, trial_numbers in (seed_from or ()):
        donor_csv = (donor if os.path.isfile(donor) else
                     os.path.join(results_dir, donor + '_trajectory.csv'))
        if not os.path.isfile(donor_csv):
            raise ValueError(f'seed_from donor {donor!r}: no trajectory CSV '
                             f'at {donor_csv}')
        points, notes = seed_points_from_trajectory(
            donor_csv, trial_numbers, search_space,
            parameter_groups=parameter_groups)
        seed_points.update(points)
        seed_notes.extend(notes)
    # Pre-flight: a study name colliding with a trajectory of a different
    # column set (search space or burden on/off changed, e.g. a legacy
    # study_name resumed without burden_model=None) must fail HERE --
    # before the sidecar recovery, the optuna store and any ~20 s
    # simulation -- not at the first row append.
    check_trajectory_header(csv_path, columns)
    # An orphaned sidecar means a previous run of this study died
    # mid-trial without writing that trial's row (unsupervised
    # crash-resume; under the supervisor it has already been recovered).
    # Log it before trial 0's own sidecar could overwrite it.
    lost = recover_inflight(csv_path, inflight_path, state='LOST',
                            error='recovered at engine startup (no terminal row)')
    if lost is not None:
        print(f'Recovered orphaned in-flight trial {lost} from a previous '
              'run as a LOST row (no terminal row had been written).')

    study = optuna.create_study(study_name=study_name, storage=storage,
                                direction=direction,
                                load_if_exists=True)
    n_done = len(study.trials)
    # Offset the seed by the number of stored trials so a resumed study
    # draws fresh points instead of replaying the original RNG stream.
    def _burden_constraints(frozen_trial):
        # <= 0 feasible. optuna evaluates this for COMPLETE and PRUNED
        # trials (samplers/_base._process_constraints_after_trial), so
        # the pruned INFEASIBLE trials populate the sampler's infeasible
        # set and steer sampling toward the feasible region. The attr is
        # set on every trial right after sampling, before any prune.
        return (frozen_trial.user_attrs.get('burden_violation', 0.0),)
    if n_startup_trials is None:
        n_startup = max(10, n_trials//10)
        startup_rule = 'default rule max(10, n_trials//10)'
    else:
        if (isinstance(n_startup_trials, bool)
                or int(n_startup_trials) != n_startup_trials
                or n_startup_trials < 0):
            raise ValueError('n_startup_trials must be a non-negative '
                             f'integer or None; got {n_startup_trials!r}')
        n_startup = int(n_startup_trials)
        startup_rule = 'explicit'
    feasible_on = bool(burden_on and feasible_sampling)
    print(f'TPE random start-up: {n_startup} trials ({startup_rule}); '
          f'{n_done} trials already stored, so guidance begins '
          f'{"now" if n_done >= n_startup else f"after trial {n_startup - 1}"}'
          + (' (feasibility-aware: joint uniform-feasible draws)'
             if feasible_on else '') + '.')
    if feasible_on:
        study.sampler = feasible_tpe_sampler(
            search_space,
            # The predicate must see what the burden model will actually be
            # given in _objective: the EXPANDED member values (baseline x
            # the sampled group multiplier), never the group key itself.
            # Without groups expand_grouped_values is an identity copy.
            lambda values: burden_model.evaluate(
                expand_grouped_values(values, parameter_groups,
                                      kinetic_baselines)).feasible,
            multivariate=True, seed=seed + n_done,
            n_startup_trials=n_startup,
            constraints_func=_burden_constraints)
        print('Sampler: feasibility-aware TPE (FeasibleTPESampler): every '
              'start-up draw and TPE candidate is checked against the '
              'enzyme-burden cap before it is proposed.')
    else:
        study.sampler = optuna.samplers.TPESampler(
            multivariate=True, seed=seed + n_done,
            n_startup_trials=n_startup,
            constraints_func=_burden_constraints if burden_on else None)
        print('Sampler: plain TPESampler '
              + ('(feasible_sampling=False).' if burden_on
                 else '(burden off: no feasibility predicate).'))
    if n_done == 0:
        # Fresh study. By default (enqueue_baseline=False, since
        # 2026-09-07) NO baseline point is enqueued, so the sampler draws
        # every trial from trial 0 -- an enqueued scenario-A baseline
        # anchored TPE in the pure-ethanol basin in every earlier IRR
        # study. enqueue_baseline=True evaluates the baseline as trial 0
        # so it provably participates. baseline_point is computed either
        # way -- the knockout probes below are derived from it.
        baseline_point = baseline_decision_point(
            search_space, kinetic_baselines, baseline_model_kwargs,
            baseline_max_n_spikes,
            baseline_stage_1_max_x=baseline_stage_1_max_x,
            parameter_groups=parameter_groups)
        if enqueue_baseline:
            study.enqueue_trial(baseline_point)
            print('Enqueued the scenario baseline configuration as trial 0.')
        else:
            print('Baseline NOT enqueued (enqueue_baseline=False): no trial '
                  'is pre-seeded; the sampler draws every trial from trial 0.')
        if enqueue_knockouts:
            # Then the single-knockout probes (one k_* at its floor, all
            # else at the baseline), FIFO in search-space order, so the
            # sampler learns which rates tolerate a lone knock-down
            # before it draws multi-parameter points. Stored WAITING in
            # the optuna DB: a crash mid-probes resumes them without any
            # re-enqueue (this block runs on a fresh study only).
            probes, at_floor = knockout_probe_points(search_space,
                                                    baseline_point,
                                                    rate_params=rate_params)
            for pname, point in probes.items():
                study.enqueue_trial(point,
                                    user_attrs={'knockout_probe': pname})
            skipped = (f'; {len(at_floor)} already at the floor, no probe: '
                       f'{at_floor}' if at_floor else '')
            # The probes follow the baseline only when it was enqueued, so
            # they start at trial 1 with enqueue_baseline and at trial 0
            # without it.
            first = 1 if enqueue_baseline else 0
            span = (f'as trials {first}-{first + len(probes) - 1} '
                    if probes else '')
            print(f'Enqueued {len(probes)} single-knockout probes {span}'
                  f'(each k_* alone at its band floor, the rest at the '
                  f'baseline){skipped}.')
        if seed_points:
            # Then the seed points of the donor studies (seed_from), after
            # the probes: a foothold in a basin another study found. Same
            # WAITING-store semantics as the probes.
            for label, point in seed_points.items():
                study.enqueue_trial(point, user_attrs={'seed': label})
            for note in seed_notes:
                print(f'Seed note: {note}')
            print(f'Enqueued {len(seed_points)} seed points from donor '
                  f'studies as the next trials: {list(seed_points)}.')
    elif seed_points:
        print(f'Resumed study: the {len(seed_points)} seed points are NOT '
              're-enqueued (a fresh study enqueues them once).')

    def _objective(trial):
        values = {name: (trial.suggest_int(name, sp['low'], sp['high'])
                         if sp.get('int') else
                         trial.suggest_float(name, sp['low'], sp['high'],
                                             log=sp['log']))
                  for name, sp in search_space.items()}
        # Group multipliers -> individual member values (baseline x m):
        # what the burden model evaluates and what reaches the model. The
        # CSV keeps the multipliers as the decision columns and records
        # the members as applied_<member>.
        applied_kinetics = expand_grouped_values(values, parameter_groups,
                                                 kinetic_baselines)
        if 'threshold_conc' in values:  # current threshold-anchored scheme
            threshold = values['threshold_conc']
            target = min(TARGET_CONC_MAX,
                         threshold + values['target_delta'])
            if 'spike_delta' in values:
                spike = min(SPIKE_CONC_MAX,
                            max(SPIKE_CONC_MIN,
                                target + values['spike_delta']))
            else:  # spike pinned at the scenario baseline
                spike = baseline_model_kwargs['spike_conc']
        else:  # legacy target-anchored scheme
            target = values['target_conc']
            threshold = max(0.0, target - values['threshold_delta'])
            spike = values['spike_conc']
        model_kwargs = dict(target_conc=target, threshold_conc=threshold,
                            spike_conc=spike)
        record = {'trial_number': trial.number, **values}
        # applied_<member> is the SAMPLED member value (its live baseline x
        # the group multiplier), recorded BEFORE any burden derating -- the
        # derated capacities are the burden columns (k_7_eff / k_8_eff).
        for members in parameter_groups.values():
            for member in members:
                record[f'applied_{member}'] = applied_kinetics[member]
        if burden_on:
            # Proteome-allocation burden (enzyme_burden.py): known from the
            # sampled values alone, so an over-cap point is logged and
            # pruned BEFORE the sidecar and the ~20 s simulation, and no
            # baseline state is disturbed. The CSV keeps the sampled
            # k_7/k_8 as the decision; the model receives the INTENDED
            # k_7/k_8 through `applied` below (see the comment further
            # down for where the derating actually happens).
            burden = burden_model.evaluate(applied_kinetics)
            record.update(burden.as_record())
            trial.set_user_attr('burden_violation', burden.violation)
            if not burden.feasible:
                record['state'] = 'INFEASIBLE'
                record['error'] = (f'enzyme burden: Phi_M {burden.Phi_M:.4f} '
                                   f'> F_flex {burden.F_flex:.4f} g/gDCW')
                append_trajectory_row(csv_path, columns, record)
                print(f'Trial {trial.number}: INFEASIBLE under the enzyme '
                      f'burden (Phi_M {burden.Phi_M:.4f} > F_flex '
                      f'{burden.F_flex:.4f} g/gDCW); pruned before simulating.')
                raise optuna.TrialPruned()
            # Feasible: the k_7/k_8 derating is applied by the shared
            # load_simulate choke point (system._apply_enzyme_burden) using
            # the active burden installed below, so the model receives the
            # INTENDED capacities here and the choke point derates growth.
            # (The pre-sim evaluate above still records the burden columns
            # and prunes over-cap points before the ~20 s simulation.)
            applied = applied_kinetics
        else:
            applied = applied_kinetics
        # The decision vector is complete here (every suggest_* has run)
        # and the hang-prone simulation has not started: record it, so a
        # hard kill / segfault during this trial leaves the sidecar for
        # recover_inflight (the supervisor, or the next engine start) to
        # log as a LOST row. Cleared in the finally only once this trial's
        # terminal row has been written (COMPLETE return, FAIL/NAN prune);
        # an interrupted trial (e.g. KeyboardInterrupt) that wrote no row
        # leaves the sidecar in place for recovery.
        write_inflight(inflight_path, columns, record)
        row_written = False
        try:
            try:
                for pname in kinetic_baselines:
                    if pname in applied:
                        setattr(r_te, pname, applied[pname])
                if 'max_n_spikes' in values:
                    fbs_spec.max_n_spikes = values['max_n_spikes']
                if 'stage_1_max_x' in values:
                    # The V406 property mirrors onto r_te AND the
                    # AerationSpec (air-supply sizing); never setattr r_te.
                    handles['V406'].stage_1_max_x = values['stage_1_max_x']
                handles['model_specification'](**model_kwargs)
                handles['latest_TEA_solution'].update(
                    handles['solve_TEA'](
                        stream_IDs=('ethanol', 'isobutanol')))
                for mname, getter in TRACKED_METRICS.items():
                    record[mname] = getter(handles)
                obj = float(objective_getter(handles))
                record['objective'] = obj
            except Exception as e:
                record['state'] = 'FAIL'
                record['error'] = repr(e)[:300]
                append_trajectory_row(csv_path, columns, record)
                row_written = True
                print(f'Trial {trial.number}: FAILED ({repr(e)[:120]})')
                raise optuna.TrialPruned() from e
            if not math.isfinite(obj):
                # NaN (not solved) or +-inf (solve_TEA reports an IRR with
                # no real root on the valid domain as -inf since
                # 2026-09-03; the sampler needs finite values, so such
                # trials are pruned like NaN ones -- the CSV keeps the raw
                # value in 'objective')
                record['state'] = 'NAN'
                append_trajectory_row(csv_path, columns, record)
                row_written = True
                print(f'Trial {trial.number}: objective is non-finite '
                      f'({obj}); pruned.')
                raise optuna.TrialPruned()
            record['state'] = 'COMPLETE'
            append_trajectory_row(csv_path, columns, record)
            row_written = True
            for mname in TRACKED_METRICS:
                trial.set_user_attr(mname, record[mname])
            if trial.number % print_status_every == 0:
                try:
                    best = study.best_value
                except Exception:  # no completed trial stored yet
                    best = np.nan
                try:
                    print(f'\nTrial {trial.number}/{n_trials}: '
                          f'{objective_name} = {obj:.6g} '
                          f'(best so far {best:.6g})\n'
                          f'integrator: {r_te.integrator.getName()}; '
                          'HXN Qbal error = '
                          f"{handles['HXN'].energy_balance_percent_error:.2f} %")
                except Exception:  # cosmetic only -- never abort the study
                    pass
            return obj
        finally:
            # Clear the sidecar only once this trial's terminal row exists.
            # A KeyboardInterrupt (manual abort) during the simulation, or
            # an exception escaping the FAIL branch's own row append, hits
            # this finally with row_written still False -- leave the
            # sidecar in place so recover_inflight logs it as a LOST row
            # at the next engine start (or the supervisor), rather than
            # silently losing the trial number.
            if row_written:
                clear_inflight(inflight_path)

    n_remaining = max(0, n_trials - n_done)
    if n_done:
        print(f'Resuming study {study_name}: {n_done} trials stored; '
              f'running {n_remaining} more (budget {n_trials}).')
    from biorefineries.isobutanol import system as _system
    if burden_on:
        _system.set_active_burden(burden_model)
    try:
        study.optimize(_objective, n_trials=n_remaining, n_jobs=1,
                       gc_after_trial=True)
    finally:
        # restore_baseline re-simulates the scenario baseline; keep the
        # burden active for that (inert at the scenario-A reference), then
        # clear it so a later burden-free caller in the same kernel is not
        # silently constrained. The nested finally guarantees the clear even
        # if restore_baseline or the sampler-print raises.
        try:
            restore_baseline(handles, kinetic_baselines,
                             baseline_model_kwargs,
                             baseline_max_n_spikes=baseline_max_n_spikes,
                             baseline_stage_1_max_x=baseline_stage_1_max_x)
            if feasible_on:
                s = study.sampler
                print(f'Feasible sampling: rejected {s.n_rejected} '
                      f'draws/candidates, {s.n_uniform_fallbacks} uniform '
                      f'fallbacks, {s.n_unfiltered} unfiltered draws.')
        finally:
            _system.set_active_burden(None)
    return study, csv_path, kinetic_baselines

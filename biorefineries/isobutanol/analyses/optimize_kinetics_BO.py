#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Bayesian global optimization of the fermentation kinetic parameters +
feeding sugar concentrations and glucose-spike cap (see
kinetic_optimization.py), run against a scenario baseline. Smoke-test
convention: importing this file loads and baseline-simulates the
biorefinery; the optimization itself runs only when a runner calls
run(...). Long-running at the default budget (2000 trials) -- NOT on the
run-without-asking list; rerunning with the same study name resumes a
crashed/interrupted study and runs only the remaining trials.

Search set and bands come from a named STUDY PRESET (default
study_target_products='ethanol_isobutanol', study_type='metabolic_protein':
start at the scenario-A baseline, the B workbook's 56 kinetic rows; log
bands by nskinetics ROLE since 2026-09-06 -- rate constants (capacity:
k_1h, k_2, ..., k_13-k_16) on [1e-3x, 10x], inhibition coefficients
(k_1ie, k_1ii, k_10ie, ...) and the regulation / affinity / self-
inhibition terms K_* on [0.1x, 10x]; k_10 (active-biomass decay) is
EXCLUDED by default -- a free lunch, not an engineering target -- so the
sampled set is 55; see ko.resolve_study_preset and run()'s docstring; the
derived study name carries the band tags `_rb0.001-10_ib0.1-10` and the
exclusion tag `_xk10`). study_target_products=None is the legacy flag
path (scenario / kinetic_bounds_scenario / single band, k_10 sampled)
for resuming older studies.

study_type='metabolic_minimal' (2026-09-07) is the compact preset: the
capacities minus k_10, k_7 and k_8 (17 for ethanol_isobutanol, 13 for
ethanol_only) on the rate band, ONE 0.2x-2x log multiplier per
inhibition-effector family (inhib_ethanol / inhib_isobutanol /
inhib_acetate, scaling every coefficient of that effector together;
recorded as applied_<member> CSV columns), no K_* terms, and the four
feeding/operating variables with the spike feed pinned at the baseline
600 g/L (no spike_delta column) -- 24 / 19 decision variables; name
kin_opt_ethanol_isobutanol_metabolic_minimal_irr_rb0.001-10_ib0.2-2_xk10+k7+k8_s1x1-50_burden.

study_type='metabolic_minimal_subset' (2026-09-07) is a STANDALONE
explicit set, not derived from metabolic_minimal: 9 listed rate
constants (k_1l, k_1h, k_1e, k_3, k_6, k_13, k_14, k_15, k_16;
ko.METABOLIC_MINIMAL_SUBSET_RATES) on the rate band, the three
inhibition-effector multipliers (ko.METABOLIC_MINIMAL_SUBSET_GROUPS,
0.2x-2x) and the three feeding variables threshold_conc / target_delta
/ max_n_spikes, with BOTH the spike feed and stage_1_max_x pinned at
the baseline (no spike_delta / stage_1_max_x column) -- 15 decision
variables for ethanol_isobutanol, 10 for ethanol_only (the listed set
intersected with the A workbook: no k_13-k_16, no isobutanol
coefficients); nothing excluded (the other rates and every K_* stay at
the baseline with no probe); name
kin_opt_ethanol_isobutanol_metabolic_minimal_subset_irr_rb0.001-10_ib0.2-2_burden
(no _x / _s1x tag).

The enzyme-burden (proteome-allocation) constraint of enzyme_burden.py
is ON by default (burden=True): sampled capacities are charged to the
cell's flexible protein sector, growth (k_7/k_8) is derated linearly as
it fills, and over-cap trials are logged INFEASIBLE and pruned before
simulating; the study name gains a '_burden' suffix, and run() prints
the burden report of the scenario-A reference and of the scenario-B
Ehrlich constants (which exceed the cap -- a finding, not repaired)
before the study starts. burden=False is the legacy burden-free study
(required to resume any study started before 2026-09-05).

Runner pattern (fresh kernel, one process):
    import runpy
    ns = runpy.run_path(r'<this file>')
    # default preset: kin_opt_ethanol_isobutanol_metabolic_protein_irr
    #   _rb0.001-10_ib0.1-10_xk10_burden
    study, csv_path = ns['run'](objective='IRR')
    # re-include k_10 (its 0.1x-10x band; no _xk10 tag)
    study, csv_path = ns['run'](objective='IRR', exclude_params=())
    # ethanol-only strain, expression/tolerance engineering only (29 params)
    study, csv_path = ns['run'](study_target_products='ethanol_only',
                                study_type='metabolic')
    # compact 24-variable space (rates minus k_10/k_7/k_8 + 3 effector
    # multipliers + 4 feeding/operating; spike pinned)
    study, csv_path = ns['run'](objective='IRR', study_type='metabolic_minimal')
    # standalone 15-variable set (9 listed rates + 3 effector multipliers
    # + 3 feeding; spike AND stage_1_max_x pinned; no _x / _s1x tag)
    study, csv_path = ns['run'](objective='IRR', study_type='metabolic_minimal_subset')
    # legacy: resume a pre-2026-09-04 study under its old name/space
    study, csv_path = ns['run'](scenario='A', kinetic_bounds_scenario='B',
                                study_target_products=None)
    # burden-free legacy study (older names have no _burden suffix)
    study, csv_path = ns['run'](objective='IRR', burden=False)
    # seeded: enqueue donor trials after the probes (name gains _seed3)
    study, csv_path = ns['run'](objective='IRR', study_type='metabolic',
                                seed_from=[('<donor study name>', [1553, 1914]),
                                           ('<other donor>', [1162])])
"""
from datetime import datetime

from biorefineries import isobutanol
isobutanol.load()

from biorefineries.isobutanol import kinetic_optimization as ko
from biorefineries.isobutanol import scenarios

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
model_specification = model.specification
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')


def kinetic_bounds_from_scenario(bounds_scenario,
                                 multiplier_bounds=(0.1, 10.0),
                                 rate_multiplier_bounds=None,
                                 rate_params=None,
                                 parameter_multiplier_bounds=None):
    """Absolute (low, high) bounds -- the multiplier band(s) around
    `bounds_scenario`'s workbook baseline -- for every positive-baseline
    kinetic row of that scenario's parameter-distributions workbook,
    keyed by te parameter name (ko.workbook_kinetic_bounds; the rate
    constants -- `rate_params` when given (the presets' capacity-role
    rows), else every lowercase k_* name -- use `rate_multiplier_bounds`
    when given, every other row uses `multiplier_bounds`; a row named in
    `parameter_multiplier_bounds` -- the presets' k_10 on 0.1x-10x --
    uses that band instead of either).
    Read directly from the workbook (no simulation), so it can
    parameterize a run of a DIFFERENT scenario: passed as
    param_bounds_override it gives the IBO-pathway rates zeroed in
    scenario A their scenario-B search bands instead of degenerate
    zero-baseline exclusion."""
    return ko.workbook_kinetic_bounds(
        bounds_scenario, multiplier_bounds=multiplier_bounds,
        rate_multiplier_bounds=rate_multiplier_bounds,
        rate_params=rate_params,
        parameter_multiplier_bounds=parameter_multiplier_bounds)


def run(scenario=None,  # 'A' or 'B'; None = the preset's start scenario
        # (or 'B' on the legacy path, study_target_products=None)
        objective='IRR',  # name in ko.OBJECTIVE_REGISTRY, or a callable
        n_trials=2000,  # TOTAL study budget (resume-aware)
        seed=None,  # None = the engine's launch-datetime default seed,
        # ko.default_seed_from_datetime: (year/day**2)*month*((hour+1)/10)*((minute+1)/10);
        # give parallel studies distinct explicit seeds (same-minute launches
        # collide). See ko.run_kinetic_optimization.
        make_plots=True,
        study_name=None,  # default: preset convention
        # kin_opt_{study_target_products}_{study_type}_{objective slug}
        # _rb{lo}-{hi}_ib{lo}-{hi}[_x{excluded}][_burden] (the effective
        # rate and inhibition bands and exclusion set, `_xk10` at the
        # default), plus _sc{scenario} / _kb{X} only when an
        # explicit override differs from the preset; legacy path:
        # kin_opt_{scenario}[_kb{X}]_{objective slug}
        kinetic_bounds_scenario=None,  # e.g. 'B': the kinetic search SET
        # (restrict_to_workbook) AND its bounds come from THAT scenario's
        # workbook (see kinetic_bounds_from_scenario); explicit
        # param_bounds_override entries win over the derived ones; on the
        # legacy path the default study_name gains a _kb{scenario} tag.
        restrict_to_workbook=True,  # kinetic decision variables = the
        # kinetic rows of the (kinetic_bounds_scenario or scenario)
        # workbook; False = every k_*/K_* on the model (the pre-2026-09-03
        # behaviour, for reproducing older studies; legacy path only). An
        # explicit include_params in engine_kwargs wins.
        study_target_products=ko.DEFAULT_STUDY_TARGET_PRODUCTS,
        # 'ethanol_only' | 'ethanol_isobutanol' | None (= legacy path)
        study_type=ko.DEFAULT_STUDY_TYPE,
        # 'metabolic' | 'metabolic_protein' | 'metabolic_minimal'
        # | 'metabolic_minimal_subset'
        burden=True,  # enzyme-burden (proteome-allocation) constraint,
        # enzyme_burden.py: default ON (study name + '_burden'); False =
        # legacy burden-free study (older studies). Do not pass the
        # engine's burden_model in engine_kwargs -- this flag owns it.
        enqueue_baseline=False,  # default OFF since 2026-09-07: a FRESH
        # study enqueues NO baseline point, so the sampler draws every
        # trial from trial 0 (an enqueued scenario-A baseline anchored TPE
        # in the pure-ethanol basin in every earlier IRR study). True =
        # evaluate the scenario baseline configuration as trial 0
        # (ko.baseline_decision_point). The knockout probes, if on, are
        # still derived from the baseline either way. Not part of the
        # study name (like enqueue_knockouts); a resume is unaffected
        # (only fresh studies enqueue).
        enqueue_knockouts=False,  # default OFF since 2026-09-07 (together
        # with enqueue_baseline, so a default fresh study enqueues NO point
        # at all). True = a FRESH study evaluates, right after the baseline
        # (if enqueued), one single-knockout probe per rate constant k_* of
        # the search space (that rate alone at its band floor, the rest at
        # the baseline; ko.knockout_probe_points), so TPE learns the
        # lethality map before sampling. Off the study name.
        n_startup_trials=None,  # TPE random start-up length (trials drawn
        # at random, probes included, before TPE guidance); None = the
        # engine's rule max(10, n_trials//4) (25 % of the budget floored at
        # 10, so 500 for 2000 trials; was n_trials//10 = 200, of which 0 of
        # 183 random draws completed in the burden-constrained preset
        # space on 2026-09-06); e.g. 20-30 shortens it. Not part of the
        # study name; a resume may change it.
        feasible_sampling=True,  # with the burden on, sample under the
        # feasibility-aware TPE (ko.feasible_tpe_sampler): start-up draws
        # and TPE candidates are checked against the burden cap BEFORE
        # they are proposed, so no sampled trial is INFEASIBLE; False =
        # the plain TPESampler (pre-2026-09-06 behaviour). Not part of
        # the study name; a resume may change it.
        seed_from=None,  # [(donor study name or trajectory-CSV path,
        # [trial numbers]), ...]: decision points of those donor trials
        # are enqueued after the knockout probes of a FRESH study
        # (ko.seed_points_from_trajectory; same decision columns required,
        # values clipped into this study's bands), giving TPE a foothold
        # in a basin another study found; the derived study name gains
        # `_seed{n}` (n = total seed count). None = no seeds.
        **engine_kwargs,  # bounds/overrides/etc. -> run_kinetic_optimization
        ):
    """Set up the scenario baseline (same recipe as the smoke tests), run
    the Bayesian optimization, and (optionally) save the trajectory plots
    next to the trajectory CSV. Returns (study, csv_path).

    STUDY PRESETS (default). `study_target_products` x `study_type` name
    the search set and bands (ko.resolve_study_preset): both target
    presets start at the scenario-A baseline; 'ethanol_only' samples the
    A workbook's rows, 'ethanol_isobutanol' the B workbook's (plus the
    Ehrlich block and isobutanol-inhibition coefficients; trial 0 = the
    A baseline with the zero Ehrlich rates clipped up to 1e-3 x their B
    baseline); 'metabolic' keeps the capacity / product-inhibition /
    lethality / substrate-regulation roles (all k_* + K_1i, K_2i, K_5i,
    K_9i), 'metabolic_protein' every row. 'metabolic_minimal' the
    capacities minus k_10/k_7/k_8 plus one 0.2x-2x multiplier per
    inhibition-effector family (parameter_groups; applied_<member>
    columns) and no K_* terms, spike pinned (spike_delta_bounds=None);
    run(parameter_groups=..., group_multiplier_bounds=...,
    spike_delta_bounds=...) override the preset's like the other keys.
    'metabolic_minimal_subset' the standalone explicit 15-variable set
    (ko.METABOLIC_MINIMAL_SUBSET_RATES / _GROUPS intersected with the
    workbook; 10 for ethanol_only), no exclusions, spike and
    stage_1_max_x pinned (the preset's stage_1_max_x_bounds=None; an
    explicit run(stage_1_max_x_bounds=(lo, hi)) re-samples it and tags
    `_s1x`).
    Bands (log-scale, x baseline)
    by nskinetics ROLE since 2026-09-06: rate constants (role capacity;
    the preset's `rate_params`, ko.rate_constant_names) [1e-3x, 10x]
    (1e-5x until later that day) EXCEPT k_10, the active-biomass decay
    capacity, on [0.1x, 10x] (the preset's `parameter_multiplier_bounds`,
    a copy of ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS -- a near-zero
    decay rate is not an engineering target); inhibition coefficients
    (product_inhibition, lethality: k_1ie, k_1ii, k_7ii, k_10ie, ...)
    [0.1x, 10x]; regulation / affinity / self-inhibition terms K_*
    [0.1x, 10x]. k_10 itself is EXCLUDED from every preset by default
    (the preset's `exclude_params`, a copy of
    ko.DEFAULT_EXCLUDED_PARAMETERS = ('k_10',), since 2026-09-06 pm): a
    lower decay rate is a free lunch for the optimizer, not an
    engineering target, so it stays at the scenario baseline and gets no
    knockout probe (the metabolic_minimal_subset preset excludes
    nothing); the sampled set is the workbook rows minus it
    (28/39/39/55). Each preset entry is applied ONLY where the caller
    passed nothing: an explicit `scenario`, `kinetic_bounds_scenario`,
    `include_params`, `exclude_params` (pass () to re-include k_10 -- it
    then samples its 0.1x-10x per-parameter band), `multiplier_bounds`,
    `rate_multiplier_bounds`, `rate_params` or
    `parameter_multiplier_bounds` (pass {} to put a re-included k_10 on
    the rate band) wins over the preset. A preset is itself a
    workbook restriction, so restrict_to_workbook=False raises. Every
    preset-derived study name carries the EFFECTIVE rate band tag
    `_rb{lo}-{hi}` (`_rb0.001-10` at the default), the inhibition band
    tag `_ib{lo}-{hi}` (`_ib0.1-10`) and the exclusion tag of the
    effective `exclude_params` (`_xk10` at the default, nothing for an
    empty set; ko.default_study_name / ko.excluded_parameters_tag), so a
    study under the current preset never resumes one started under the
    1e-5x rate band (e.g. kin_opt_ethanol_isobutanol_metabolic_irr_ib0.1-
    10_burden, untagged because only a differing band was tagged then),
    under the pre-role prefix rule (e.g. the aborted 2026-09-05
    kin_opt_ethanol_isobutanol_metabolic_irr_burden), or with k_10 still
    sampled (the 2026-09-06 production study kin_opt_ethanol_isobutanol_
    metabolic_irr_rb0.001-10_ib0.1-10_burden -- one more CSV column, so
    the header guard would refuse it anyway); resume those with an
    explicit study_name (and exclude_params=() for the k_10 studies).
    The per-parameter k_10 band is part of the preset's identity and is
    not tagged.

    LEGACY PATH (`study_target_products=None`): exactly the pre-preset
    behaviour, for resuming older studies. `scenario` (default 'B') is
    the STARTING state (its workbook distributions are loaded, its
    feeding baseline set, enqueued as trial 0 only when
    enqueue_baseline=True, and restored in the finally);
    `kinetic_bounds_scenario or scenario` is the
    workbook that defines WHICH kinetic parameters are decision
    variables (when restrict_to_workbook, the default) and centers their
    single multiplier band. run(scenario='A', kinetic_bounds_scenario='B',
    study_target_products=None) therefore starts at the A baseline over
    B's 56-parameter set. A restricted study has a different search
    space than a pre-2026-09-03 full-set study of the same name: reusing
    its trajectory CSV raises in append_trajectory_row -- use a fresh
    study_name.

    ENZYME BURDEN (`burden`, default True). The engine's
    run_kinetic_optimization(burden_model=...) receives a BurdenModel
    snapshotted from the live kinetic parameters right after the
    scenario baseline is set (so it is exactly inert at the scenario
    baseline); the
    study name carries ko.BURDEN_STUDY_SUFFIX on both naming paths.
    Before the study starts, the burden reports of (a) the scenario
    reference and (b) the scenario-B Ehrlich constants
    (eb.scenario_b_ehrlich()) on that reference are printed -- (b) is
    the Q11 sanity report: B's constants need ~0.22 g/gDCW of Ehrlich
    enzyme and exceed F_flex, so a burden study explores lower Ehrlich
    capacities and lower growth than the unburdened kin_opt_B_irr did.
    A burden study cannot START from an infeasible reference (e.g.
    scenario='B'): BurdenModel.from_reference raises; pass burden=False.

    BASELINE ENQUEUE (`enqueue_baseline`, default False since
    2026-09-07). By default a fresh study enqueues NO baseline point, so
    the sampler draws every trial from trial 0 (an enqueued scenario-A
    baseline anchored TPE in the pure-ethanol basin in every earlier IRR
    study). enqueue_baseline=True evaluates the scenario baseline as
    trial 0 so it provably participates (forwarded to
    ko.run_kinetic_optimization; supervisor --enqueue-baseline). Not part
    of the study name; only a fresh study enqueues.

    SINGLE-KNOCKOUT PROBES (`enqueue_knockouts`, default False since
    2026-09-07; probes added 2026-09-06). A fresh study enqueues, after
    the baseline (if any), one probe per
    RATE CONSTANT of the search space (the preset's `rate_params`, role
    capacity; inhibition coefficients k_*i* get none -- their 0.1x floor
    is a knock-down, not a knock-out) -- that rate alone at the floor
    of its band, everything else at the scenario baseline
    (ko.knockout_probe_points; a rate already at its floor, such as the
    ethanol_isobutanol preset's clipped Ehrlich rates, gets none, and
    neither does an excluded rate -- k_10 by default). Under
    the preset band the floor is an effective knock-out (1e-3x); under a
    narrower `rate_multiplier_bounds`, or for a re-included k_10 on its
    0.1x per-parameter band, a knock-down. The effective
    `rate_multiplier_bounds` (explicit or the preset's) is always tagged
    into the derived study name `_rb{lo}-{hi}` (ko.default_study_name),
    so a study never resumes one of the same name under another band.

    TPE START-UP (`n_startup_trials`, default None). The number of
    trials optuna draws uniformly at random (the baseline and the
    probes count) before TPE's density guidance starts; None keeps the
    engine's rule max(10, n_trials//4). Forwarded to
    ko.run_kinetic_optimization; compared with the trials already
    stored, so it can be changed on a resume and never tags the study
    name (supervisor --n-startup-trials).

    FEASIBLE SAMPLING (`feasible_sampling`, default True; since
    2026-09-06). With the burden on, the engine samples under
    ko.feasible_tpe_sampler: the random start-up draws are joint
    uniform-feasible vectors and every TPE candidate is filtered by the
    cap Phi_M < F_flex before it is scored, so no sampled trial should
    ever be logged INFEASIBLE (the 2026-09-06 production study proposed
    118 of 200 start-up and 261 of 1060 TPE trials over the cap). The
    in-objective guard stays as the safety net and the engine prints the
    rejection counters at the end. Meaningless with burden=False. A
    sampler setting: same columns, no study-name tag, so the production
    study resumes under it (supervisor --no-feasible-sampling).

    STAGE-1 CUTOFF (`stage_1_max_x_bounds`, an engine kwarg defaulted
    from the preset since 2026-09-06 pm). Every preset but
    metabolic_minimal_subset samples the fermentor's aerobic stage-1
    biomass cutoff V406.stage_1_max_x
    (ko.OPERATING_VARIABLES) log-scale on ko.DEFAULT_STAGE_1_MAX_X_BOUNDS
    = (1.0, 50.0) g/L (baseline 5.0 g/L, the nskinetics factory default;
    stage 1 also ends at stage_1_max_time = 25 h, whichever fires first).
    The engine sets the V406 property -- which mirrors onto the kinetic
    model and the aeration spec -- before every simulation and restores
    it in its finally; no knockout probe (not a rate constant). Pass
    run(stage_1_max_x_bounds=None) to pin it at the baseline (the
    variable is then absent and the name carries no tag, so an older
    study resumes as before) or an explicit (lo, hi). The effective band
    is tagged `_s1x{lo}-{hi}` into the derived study name (a new decision
    COLUMN: the header guard refuses to resume a study without it, and
    the tag is what lets the default name launch next to the existing
    `..._xk10_burden` studies). Legacy path (study_target_products=None):
    not sampled unless passed explicitly.

    SEED POINTS (`seed_from`, default None; since 2026-09-07). A fresh
    study enqueues, right after its knockout probes, the decision points
    of the named trials of DONOR studies: seed_from = [(donor, [trial,
    ...]), ...], donor = a study name (resolved to
    analyses/results/{donor}_trajectory.csv) or a CSV path
    (ko.seed_points_from_trajectory). The donor must have sampled the
    SAME decision columns (ValueError otherwise); values are clipped into
    this study's bands; each seed trial carries the optuna user attr
    'seed' = '{donor}#{trial}'. Motivation: the 2026-09-06 IRR study of
    the metabolic preset converged on a pure-ethanol cell (IRR 0.174)
    without ever sampling an isobutanol titer above 10 g/L, while the
    IBO-yield-x-titer study's trial 1553 (IBO 72 g/L) scores IRR 0.186 on
    the same model -- a local optimum TPE could not leave because it had
    no foothold in the narrow burden-feasible IBO basin. Pick seeds by
    re-scoring the donor's rows under the NEW objective (its tracked
    metrics are in the CSV), plus the donor's own optimum and the other
    product-mix corner; a seed that scores below the incumbent basin only
    teaches TPE to avoid its region. Same columns as the unseeded study,
    so the derived name is tagged `_seed{n}` (n = the total seed count;
    ko.seed_points_tag) -- a different panel of the same size needs an
    explicit study_name. Resumes never re-enqueue."""
    if 'burden_model' in engine_kwargs:
        raise ValueError("pass burden=True/False to run(), not the engine's "
                         'burden_model (run() builds it so the reports can '
                         'be printed first).')
    slug = (objective if isinstance(objective, str)
            else engine_kwargs.get('objective_name', 'custom')
            ).lower().replace(' ', '_')
    seed_from = [(donor, tuple(int(n) for n in trials))
                 for donor, trials in (seed_from or ())]
    n_seeds = sum(len(trials) for _, trials in seed_from)
    if study_target_products is not None:
        if not restrict_to_workbook:
            raise ValueError(
                'restrict_to_workbook=False cannot be combined with a study '
                'preset (the preset IS a workbook restriction); pass '
                'study_target_products=None for the legacy all-model set.')
        preset = ko.resolve_study_preset(study_target_products, study_type)
        if scenario is None:
            scenario = preset['scenario']
        if kinetic_bounds_scenario is None:
            kinetic_bounds_scenario = preset['kinetic_bounds_scenario']
        for key in ('include_params', 'exclude_params', 'multiplier_bounds',
                    'rate_multiplier_bounds', 'rate_params',
                    'parameter_multiplier_bounds', 'stage_1_max_x_bounds',
                    'parameter_groups', 'group_multiplier_bounds',
                    'spike_delta_bounds'):
            engine_kwargs.setdefault(key, preset[key])
        if study_name is None:
            study_name = ko.default_study_name(
                objective if isinstance(objective, str)
                else engine_kwargs.get('objective_name', 'custom'),
                study_target_products, study_type,
                scenario=scenario, kinetic_bounds_scenario=kinetic_bounds_scenario,
                burden=burden,
                # The EFFECTIVE bands (explicit or the preset's) are tagged
                # on EVERY preset name: the rate band (so a 1e-3x study
                # never resumes the untagged 1e-5x study of the same
                # objective) and the inhibition coefficients' saturation
                # band (role-band scheme).
                rate_multiplier_bounds=engine_kwargs['rate_multiplier_bounds'],
                # Under a GROUPED study type the band that actually sizes
                # the inhibition entries is group_multiplier_bounds (the
                # members are not sampled individually), so an explicit
                # group band must get its own study: same columns, same
                # numeric-range-tolerant optuna store, so only the name
                # keeps it off the preset-band study.
                inhibition_multiplier_bounds=(
                    engine_kwargs['group_multiplier_bounds']
                    if engine_kwargs['parameter_groups'] else
                    engine_kwargs['multiplier_bounds']),
                # The effective exclusion set (explicit or the preset's
                # DEFAULT_EXCLUDED_PARAMETERS): an excluded name is a
                # missing CSV column, so the tag keeps the default name
                # off the studies that still sampled it.
                exclude_params=engine_kwargs['exclude_params'],
                # The operating variable's band (the preset's, an explicit
                # one, or None = pinned -> no tag): a new column, tagged
                # so the default name launches next to the older studies.
                stage_1_max_x_bounds=engine_kwargs['stage_1_max_x_bounds'],
                # The seed count: same columns as the unseeded study, so
                # the tag is what keeps a seeded run off its store.
                n_seeds=n_seeds)
        excluded = tuple(engine_kwargs['exclude_params'] or ())
        groups = dict(engine_kwargs['parameter_groups'] or {})
        grouped = {m for members in groups.values() for m in members}
        effective = [n for n in engine_kwargs['include_params']
                     if n not in excluded and n not in grouped]
        n_rate = sum(1 for n in effective if n in engine_kwargs['rate_params'])
        print(f'Study preset: study_target_products={study_target_products!r}, '
              f'study_type={study_type!r} -> start at scenario {scenario}, '
              f'{len(effective)} individually sampled kinetic parameters '
              f'from the scenario-{kinetic_bounds_scenario} workbook '
              f'({len(engine_kwargs["include_params"])} rows minus the '
              f'excluded {list(excluded) or "none"}, which stay at the '
              f'baseline, minus the {len(grouped)} grouped ones below); bands '
              f'(x baseline, log) by role: {n_rate} rate constants '
              f'{engine_kwargs["rate_multiplier_bounds"]}, the other '
              f'{len(effective) - n_rate} (inhibition '
              f'coefficients, K_* terms) {engine_kwargs["multiplier_bounds"]}; '
              'per-parameter bands '
              f'{engine_kwargs["parameter_multiplier_bounds"] or "none"}; '
              'operating variable stage_1_max_x (aerobic stage-1 biomass '
              'cutoff, log-scale, g/L) '
              + ('pinned at the baseline'
                 if engine_kwargs['stage_1_max_x_bounds'] is None else
                 f'on {tuple(engine_kwargs["stage_1_max_x_bounds"])}')
              + '; spike feed '
              + ('pinned at the baseline (spike_delta_bounds=None)'
                 if engine_kwargs['spike_delta_bounds'] is None else
                 f'sampled as spike_delta on '
                 f'{tuple(engine_kwargs["spike_delta_bounds"])}')
              + '.')
        if groups:
            print('Parameter groups (one log-scale multiplier each on '
                  f'{tuple(engine_kwargs["group_multiplier_bounds"])} x '
                  'baseline, preserving intra-group ratios): '
                  + '; '.join(f'{g}[{len(m)}]: {", ".join(m)}'
                              for g, m in groups.items())
                  + '.')
    elif scenario is None:
        scenario = 'B'  # legacy default
    param_set_scenario = kinetic_bounds_scenario or scenario
    if restrict_to_workbook:
        if 'include_params' in engine_kwargs:
            if study_target_products is None:
                print('Kinetic search set: the '
                      f'{len(engine_kwargs["include_params"])} explicitly '
                      'passed include_params (the workbook set is not '
                      'derived).')
        else:
            names = ko.kinetic_param_names_from_scenario(param_set_scenario)
            engine_kwargs['include_params'] = names
            print(f'Kinetic search set: the {len(names)} kinetic rows of the '
                  f'scenario-{param_set_scenario} parameter-distributions '
                  'workbook.')
    if kinetic_bounds_scenario is not None:
        derived = kinetic_bounds_from_scenario(
            kinetic_bounds_scenario,
            multiplier_bounds=engine_kwargs.get('multiplier_bounds',
                                                (0.1, 10.0)),
            rate_multiplier_bounds=engine_kwargs.get('rate_multiplier_bounds'),
            rate_params=engine_kwargs.get('rate_params'),
            parameter_multiplier_bounds=engine_kwargs.get('parameter_multiplier_bounds'))
        derived.update(engine_kwargs.get('param_bounds_override') or {})
        engine_kwargs['param_bounds_override'] = derived
        if study_name is None:  # legacy path only (presets set it above)
            study_name = (f'kin_opt_{scenario}_'
                          f'kb{kinetic_bounds_scenario}_{slug}'
                          + ko.seed_points_tag(n_seeds)
                          + (ko.BURDEN_STUDY_SUFFIX if burden else ''))
    if seed_from:
        print(f'Seed points: {n_seeds} donor trials enqueued after the '
              f'probes of a fresh study: {seed_from}')
    # Consolidated scenario baseline: workbook kinetics + distributions +
    # feeding strategy + one baseline model_specification (single source of
    # truth in scenarios.SCENARIOS). The BO samples on top of this baseline.
    bundle = scenarios.load_scenario(scenario, burden=burden)

    if burden:
        from biorefineries.isobutanol import enzyme_burden as eb
        burden_model = bundle['burden_model']   # A-referenced, installed active
        print(burden_model.describe_point(
            burden_model.reference,
            label=f'scenario-{scenario} reference (A-calibrated)'
                  + (' (trial 0)' if enqueue_baseline else '')))
        print(burden_model.describe_point(
            {**burden_model.reference, **eb.scenario_b_ehrlich()},
            label='scenario-B Ehrlich constants on this reference '
                  '(Q11 sanity report; infeasible by design, not repaired)'))
    else:
        burden_model = None
        print('Enzyme burden OFF (burden=False): legacy burden-free study.')

    study, csv_path, kinetic_baselines = ko.run_kinetic_optimization(
        objective=objective,
        scenario_label=scenario,
        n_trials=n_trials,
        seed=seed,
        study_name=study_name,
        burden_model=burden_model,
        enqueue_baseline=enqueue_baseline,
        enqueue_knockouts=enqueue_knockouts,
        n_startup_trials=n_startup_trials,
        feasible_sampling=feasible_sampling,
        seed_from=seed_from,
        **engine_kwargs)

    if make_plots:
        try:
            import optuna
            direction = ('maximize'
                         if study.direction == optuna.study.StudyDirection.MAXIMIZE
                         else 'minimize')
            objective_name = (objective if isinstance(objective, str)
                              else engine_kwargs.get('objective_name', 'custom'))
            objective_units = (ko.OBJECTIVE_REGISTRY[objective]['units']
                               if isinstance(objective, str)
                               else engine_kwargs.get('objective_units', ''))
            df = ko.load_trajectory(csv_path)
            stamp = datetime.now().strftime('%Y.%m.%d-%H.%M')
            base = csv_path[:-len('_trajectory.csv')]
            ko.plot_optimization_trajectories(
                df, objective_name=objective_name, direction=direction,
                objective_units=objective_units,
                filename=base + f'_trajectories_{stamp}.png')
            # Group multipliers are plotted directly (baseline 1.0); the
            # grouped members are not decision columns and are skipped.
            plot_baselines = {**kinetic_baselines,
                              **{g: 1.0 for g in
                                 (engine_kwargs.get('parameter_groups') or {})}}
            ko.plot_parameter_trajectory(
                df, plot_baselines, direction=direction,
                filename=base + f'_param_trajectory_{stamp}.png')
            ko.plot_best_vs_baseline(
                df, plot_baselines, direction=direction,
                filename=base + f'_best_vs_baseline_{stamp}.png')
            # Log-sampled columns mirror build_search_space: kinetic
            # params are log-scale unless a param_bounds_override entry
            # has low <= 0.
            override = engine_kwargs.get('param_bounds_override') or {}
            log_columns = {
                p for p in kinetic_baselines if p in df.columns
                and (override[p][0] > 0.0 if p in override else True)
            } | set(engine_kwargs.get('parameter_groups') or ())  # log multipliers
            ko.plot_pca_projection(
                df, direction=direction, log_columns=log_columns,
                objective_name=objective_name,
                objective_units=objective_units,
                filename=base + f'_pca_{stamp}.png',
                # Trial 0 is the baseline ONLY when it was enqueued;
                # otherwise it is a sampled draw and gets no marker.
                baseline_trial=(0 if enqueue_baseline else None))
            print(f'Plots saved next to {csv_path}')
        except Exception as e:
            print('Plotting failed (the trajectory CSV and study are '
                  f'intact on disk): {repr(e)[:300]}')

    return study, csv_path

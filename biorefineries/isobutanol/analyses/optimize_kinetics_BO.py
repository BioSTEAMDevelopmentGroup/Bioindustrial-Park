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
k_1h, k_2, ..., k_13-k_16) on [1e-5x, 10x], inhibition coefficients
(k_1ie, k_1ii, k_10ie, ...) and the regulation / affinity / self-
inhibition terms K_* on [0.1x, 10x]; see ko.resolve_study_preset and
run()'s docstring; the derived study name carries the inhibition band
tag `_ib0.1-10`). study_target_products=None is the legacy flag path
(scenario / kinetic_bounds_scenario / single band) for resuming older
studies.

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
    study, csv_path = ns['run'](objective='IRR')
    # ethanol-only strain, expression/tolerance engineering only (29 params)
    study, csv_path = ns['run'](study_target_products='ethanol_only',
                                study_type='metabolic')
    # legacy: resume a pre-2026-09-04 study under its old name/space
    study, csv_path = ns['run'](scenario='A', kinetic_bounds_scenario='B',
                                study_target_products=None)
    # burden-free legacy study (older names have no _burden suffix)
    study, csv_path = ns['run'](objective='IRR', burden=False)
"""
from datetime import datetime

from biorefineries import isobutanol
isobutanol.load()

from biorefineries.isobutanol import kinetic_optimization as ko

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
        seed=3221,
        make_plots=True,
        study_name=None,  # default: preset convention
        # kin_opt_{study_target_products}_{study_type}_{objective slug}
        # _rb{lo}-{hi}_ib{lo}-{hi}[_burden] (the effective rate and
        # inhibition bands), plus _sc{scenario} / _kb{X} only when an
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
        # 'metabolic' | 'metabolic_protein'
        burden=True,  # enzyme-burden (proteome-allocation) constraint,
        # enzyme_burden.py: default ON (study name + '_burden'); False =
        # legacy burden-free study (older studies). Do not pass the
        # engine's burden_model in engine_kwargs -- this flag owns it.
        enqueue_knockouts=True,  # a FRESH study evaluates, right after
        # the baseline (trial 0), one single-knockout probe per rate
        # constant k_* of the search space (that rate alone at its band
        # floor; ko.knockout_probe_points), so TPE learns the lethality
        # map before sampling; False = baseline only (pre-2026-09-06).
        n_startup_trials=None,  # TPE random start-up length (trials drawn
        # at random, probes included, before TPE guidance); None = the
        # engine's rule max(10, n_trials//10) (200 for 2000 trials -- 0 of
        # 183 random draws completed in the burden-constrained preset
        # space on 2026-09-06); e.g. 20-30 shortens it. Not part of the
        # study name; a resume may change it.
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
    K_9i), 'metabolic_protein' every row. Bands (log-scale, x baseline)
    by nskinetics ROLE since 2026-09-06: rate constants (role capacity;
    the preset's `rate_params`, ko.rate_constant_names) [1e-3x, 10x]
    (1e-5x until later that day) EXCEPT k_10, the active-biomass decay
    capacity, on [0.1x, 10x] (the preset's `parameter_multiplier_bounds`,
    a copy of ko.DEFAULT_PARAMETER_MULTIPLIER_BOUNDS -- a near-zero
    decay rate is not an engineering target); inhibition coefficients
    (product_inhibition, lethality: k_1ie, k_1ii, k_7ii, k_10ie, ...)
    [0.1x, 10x]; regulation / affinity / self-inhibition terms K_*
    [0.1x, 10x]. Each preset entry is applied ONLY where the caller
    passed nothing: an explicit `scenario`, `kinetic_bounds_scenario`,
    `include_params`, `multiplier_bounds`, `rate_multiplier_bounds`,
    `rate_params` or `parameter_multiplier_bounds` (pass {} to put k_10
    back on the rate band) wins over the preset. A preset is itself a
    workbook restriction, so restrict_to_workbook=False raises. Every
    preset-derived study name carries the EFFECTIVE rate band tag
    `_rb{lo}-{hi}` (`_rb0.001-10` at the default) and the inhibition band
    tag `_ib{lo}-{hi}` (`_ib0.1-10`; ko.default_study_name), so a study
    under the current bands never resumes one started under the 1e-5x
    rate band (e.g. kin_opt_ethanol_isobutanol_metabolic_irr_ib0.1-10_
    burden, untagged because only a differing band was tagged then) or
    under the pre-role prefix rule (e.g. the aborted 2026-09-05
    kin_opt_ethanol_isobutanol_metabolic_irr_burden); resume those with
    an explicit study_name. The per-parameter k_10 band is part of the
    preset's identity and is not tagged.

    LEGACY PATH (`study_target_products=None`): exactly the pre-preset
    behaviour, for resuming older studies. `scenario` (default 'B') is
    the STARTING state (its workbook distributions are loaded, its
    feeding baseline set, it is enqueued as trial 0, and it is restored
    in the finally); `kinetic_bounds_scenario or scenario` is the
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
    scenario baseline is set (so it is exactly inert at trial 0); the
    study name carries ko.BURDEN_STUDY_SUFFIX on both naming paths.
    Before the study starts, the burden reports of (a) the scenario
    reference and (b) the scenario-B Ehrlich constants
    (eb.scenario_b_ehrlich()) on that reference are printed -- (b) is
    the Q11 sanity report: B's constants need ~0.22 g/gDCW of Ehrlich
    enzyme and exceed F_flex, so a burden study explores lower Ehrlich
    capacities and lower growth than the unburdened kin_opt_B_irr did.
    A burden study cannot START from an infeasible reference (e.g.
    scenario='B'): BurdenModel.from_reference raises; pass burden=False.

    SINGLE-KNOCKOUT PROBES (`enqueue_knockouts`, default True; since
    2026-09-06). A fresh study enqueues, after trial 0, one probe per
    RATE CONSTANT of the search space (the preset's `rate_params`, role
    capacity; inhibition coefficients k_*i* get none -- their 0.1x floor
    is a knock-down, not a knock-out) -- that rate alone at the floor
    of its band, everything else at the scenario baseline
    (ko.knockout_probe_points; a rate already at its floor, such as the
    ethanol_isobutanol preset's clipped Ehrlich rates, gets none). Under
    the preset band the floor is an effective knock-out (1e-3x); under a
    narrower `rate_multiplier_bounds`, or for k_10 on its 0.1x
    per-parameter band, a knock-down. The effective
    `rate_multiplier_bounds` (explicit or the preset's) is always tagged
    into the derived study name `_rb{lo}-{hi}` (ko.default_study_name),
    so a study never resumes one of the same name under another band.

    TPE START-UP (`n_startup_trials`, default None). The number of
    trials optuna draws uniformly at random (the baseline and the
    probes count) before TPE's density guidance starts; None keeps the
    engine's rule max(10, n_trials//10). Forwarded to
    ko.run_kinetic_optimization; compared with the trials already
    stored, so it can be changed on a resume and never tags the study
    name (supervisor --n-startup-trials)."""
    if 'burden_model' in engine_kwargs:
        raise ValueError("pass burden=True/False to run(), not the engine's "
                         'burden_model (run() builds it so the reports can '
                         'be printed first).')
    slug = (objective if isinstance(objective, str)
            else engine_kwargs.get('objective_name', 'custom')
            ).lower().replace(' ', '_')
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
        for key in ('include_params', 'multiplier_bounds',
                    'rate_multiplier_bounds', 'rate_params',
                    'parameter_multiplier_bounds'):
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
                inhibition_multiplier_bounds=engine_kwargs['multiplier_bounds'])
        n_rate = sum(1 for n in engine_kwargs['include_params']
                     if n in engine_kwargs['rate_params'])
        print(f'Study preset: study_target_products={study_target_products!r}, '
              f'study_type={study_type!r} -> start at scenario {scenario}, '
              f'{len(engine_kwargs["include_params"])} kinetic parameters '
              f'from the scenario-{kinetic_bounds_scenario} workbook; bands '
              f'(x baseline, log) by role: {n_rate} rate constants '
              f'{engine_kwargs["rate_multiplier_bounds"]}, the other '
              f'{len(engine_kwargs["include_params"]) - n_rate} (inhibition '
              f'coefficients, K_* terms) {engine_kwargs["multiplier_bounds"]}; '
              'per-parameter bands '
              f'{engine_kwargs["parameter_multiplier_bounds"] or "none"}.')
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
                          + (ko.BURDEN_STUDY_SUFFIX if burden else ''))
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename,
                                       namespace_dict)
    model.metrics_at_baseline()

    if scenario == 'A':
        fbs_spec.max_n_spikes = 16
        baseline_kwargs = dict(threshold_conc=217.125, target_conc=221.25)
    elif scenario == 'B':
        fbs_spec.max_n_spikes = 0  # batch: no glucose spikes
        baseline_kwargs = dict(threshold_conc=34.25, target_conc=140.0)
    else:
        raise ValueError(f'Scenario {scenario} not found.')
    model_specification(**baseline_kwargs)

    if burden:
        from biorefineries.isobutanol import enzyme_burden as eb
        k_ref = ko.discover_kinetic_parameters(ko.get_handles()['r_te'])
        burden_model = eb.BurdenModel.from_reference(k_ref)
        print(burden_model.describe_point(
            burden_model.reference,
            label=f'scenario-{scenario} reference (trial 0)'))
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
        enqueue_knockouts=enqueue_knockouts,
        n_startup_trials=n_startup_trials,
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
            ko.plot_parameter_trajectory(
                df, kinetic_baselines, direction=direction,
                filename=base + f'_param_trajectory_{stamp}.png')
            ko.plot_best_vs_baseline(
                df, kinetic_baselines, direction=direction,
                filename=base + f'_best_vs_baseline_{stamp}.png')
            # Log-sampled columns mirror build_search_space: kinetic
            # params are log-scale unless a param_bounds_override entry
            # has low <= 0.
            override = engine_kwargs.get('param_bounds_override') or {}
            log_columns = {
                p for p in kinetic_baselines if p in df.columns
                and (override[p][0] > 0.0 if p in override else True)}
            ko.plot_pca_projection(
                df, direction=direction, log_columns=log_columns,
                objective_name=objective_name,
                objective_units=objective_units,
                filename=base + f'_pca_{stamp}.png')
            print(f'Plots saved next to {csv_path}')
        except Exception as e:
            print('Plotting failed (the trajectory CSV and study are '
                  f'intact on disk): {repr(e)[:300]}')

    return study, csv_path

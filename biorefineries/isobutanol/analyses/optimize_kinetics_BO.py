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
start at the scenario-A baseline, the B workbook's 56 kinetic rows, k_* on
[1e-5x, 10x] and K_* on [0.1x, 10x] log bands; see ko.resolve_study_preset
and run()'s docstring). study_target_products=None is the legacy flag path
(scenario / kinetic_bounds_scenario / single band) for resuming older
studies.

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
                                 rate_multiplier_bounds=None):
    """Absolute (low, high) bounds -- the multiplier band(s) around
    `bounds_scenario`'s workbook baseline -- for every positive-baseline
    kinetic row of that scenario's parameter-distributions workbook,
    keyed by te parameter name (ko.workbook_kinetic_bounds; k_* names use
    `rate_multiplier_bounds` when given, K_* use `multiplier_bounds`).
    Read directly from the workbook (no simulation), so it can
    parameterize a run of a DIFFERENT scenario: passed as
    param_bounds_override it gives the IBO-pathway rates zeroed in
    scenario A their scenario-B search bands instead of degenerate
    zero-baseline exclusion."""
    return ko.workbook_kinetic_bounds(
        bounds_scenario, multiplier_bounds=multiplier_bounds,
        rate_multiplier_bounds=rate_multiplier_bounds)


def run(scenario=None,  # 'A' or 'B'; None = the preset's start scenario
        # (or 'B' on the legacy path, study_target_products=None)
        objective='IRR',  # name in ko.OBJECTIVE_REGISTRY, or a callable
        n_trials=2000,  # TOTAL study budget (resume-aware)
        seed=3221,
        make_plots=True,
        study_name=None,  # default: preset convention
        # kin_opt_{study_target_products}_{study_type}_{objective slug},
        # plus _sc{scenario} / _kb{X} only when an explicit override
        # differs from the preset; legacy path:
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
    A baseline with the zero Ehrlich rates clipped up to 1e-5 x their B
    baseline); 'metabolic' keeps the capacity / product-inhibition /
    lethality / substrate-regulation roles (all k_* + K_1i, K_2i, K_5i,
    K_9i), 'metabolic_protein' every row. Bands: k_* [1e-5x, 10x], K_*
    [0.1x, 10x], log-scale. Each preset entry is applied ONLY where the
    caller passed nothing: an explicit `scenario`,
    `kinetic_bounds_scenario`, `include_params`, `multiplier_bounds` or
    `rate_multiplier_bounds` wins over the preset. A preset is itself a
    workbook restriction, so restrict_to_workbook=False raises.

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
    study_name."""
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
                    'rate_multiplier_bounds'):
            engine_kwargs.setdefault(key, preset[key])
        if study_name is None:
            study_name = ko.default_study_name(
                objective if isinstance(objective, str)
                else engine_kwargs.get('objective_name', 'custom'),
                study_target_products, study_type,
                scenario=scenario, kinetic_bounds_scenario=kinetic_bounds_scenario)
        print(f'Study preset: study_target_products={study_target_products!r}, '
              f'study_type={study_type!r} -> start at scenario {scenario}, '
              f'{len(engine_kwargs["include_params"])} kinetic parameters '
              f'from the scenario-{kinetic_bounds_scenario} workbook, bands '
              f'k_* {engine_kwargs["rate_multiplier_bounds"]} / K_* '
              f'{engine_kwargs["multiplier_bounds"]} (x baseline, log).')
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
            rate_multiplier_bounds=engine_kwargs.get('rate_multiplier_bounds'))
        derived.update(engine_kwargs.get('param_bounds_override') or {})
        engine_kwargs['param_bounds_override'] = derived
        if study_name is None:  # legacy path only (presets set it above)
            study_name = (f'kin_opt_{scenario}_'
                          f'kb{kinetic_bounds_scenario}_{slug}')
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

    study, csv_path, kinetic_baselines = ko.run_kinetic_optimization(
        objective=objective,
        scenario_label=scenario,
        n_trials=n_trials,
        seed=seed,
        study_name=study_name,
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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Bayesian global optimization of all fermentation kinetic parameters +
feeding sugar concentrations (see kinetic_optimization.py), run against a
scenario baseline. Smoke-test convention: importing this file loads and
baseline-simulates the biorefinery; the optimization itself runs only when
a runner calls run(...). Long-running at the default budget (2000 trials)
-- NOT on the run-without-asking list; rerunning with the same study name
resumes a crashed/interrupted study and runs only the remaining trials.

Runner pattern (fresh kernel, one process):
    import runpy
    ns = runpy.run_path(r'<this file>')
    study, csv_path = ns['run'](scenario='B', objective='IRR')
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


def run(scenario='B',  # 'A' or 'B'
        objective='IRR',  # name in ko.OBJECTIVE_REGISTRY, or a callable
        n_trials=2000,  # TOTAL study budget (resume-aware)
        seed=3221,
        make_plots=True,
        study_name=None,  # default: kin_opt_{scenario}_{objective slug}
        **engine_kwargs,  # bounds/overrides/etc. -> run_kinetic_optimization
        ):
    """Set up the scenario baseline (same recipe as the smoke tests), run
    the Bayesian optimization, and (optionally) save the three trajectory
    plots next to the trajectory CSV. Returns (study, csv_path)."""
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
        fbs_spec.max_n_spikes = 13
        baseline_kwargs = dict(threshold_conc=216.3, target_conc=226.3)
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
        import optuna
        direction = ('maximize'
                     if study.direction == optuna.study.StudyDirection.MAXIMIZE
                     else 'minimize')
        objective_name = (objective if isinstance(objective, str)
                          else engine_kwargs.get('objective_name', 'custom'))
        objective_units = (ko.OBJECTIVE_REGISTRY[objective]['units']
                           if isinstance(objective, str) else '')
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
        print(f'Plots saved next to {csv_path}')

    return study, csv_path

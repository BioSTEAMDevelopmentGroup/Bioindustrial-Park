# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 18:10:51 2026

@author: sarangbhagwat
"""
import biosteam as bst
from biorefineries import isobutanol
from matplotlib import pyplot as plt

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
optimize_1D_feeding_strategy_for_MPSP = isobutanol.models.optimize_1D_feeding_strategy_for_MPSP
optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP =  isobutanol.models.optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP
optimize_max_n_glu_spikes_for_MPSP = isobutanol.models.optimize_max_n_glu_spikes_for_MPSP
plot_kinetic_results = isobutanol.models.plot_kinetic_results
unit_groups_dict = isobutanol.models.unit_groups_dict
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406
r = V406.nsk_kinetic_model
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')


def load_simulate_baseline(scenario='A', # 'A' or 'B'
                           plot=False,
                           **tea_kwargs, # TEA-solve kwargs forwarded to
                                         # load_simulate_solve_TEA via
                                         # model_specification: solve_for
                                         # ('MPSP'|'IRR'|'NPV'),
                                         # product_stream, IRR, NPV,
                                         # n_tea_solves
                           ):
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
    baseline_initial = model.metrics_at_baseline()

    model_specification(**tea_kwargs)

    # for forced batch mode:
    # fbs_spec.max_n_spikes = 0
    if scenario=='A':
        fbs_spec.max_n_spikes = 16
        result = model_specification(threshold_conc=217.125, target_conc=221.25, **tea_kwargs)
    elif scenario=='B':
        fbs_spec.max_n_spikes = 13
        result = model_specification(threshold_conc=216.3, target_conc=226.3, **tea_kwargs)
    else:
        raise ValueError(f'Scenario {scenario} not found.')


    if plot:
        # Plot conc v time
        fig, ax = plot_kinetic_results(xlim=(0,80), ylim=(0,250))

    return result # purity-adjusted MPSP, IRR, or NPV per solve_for (default: ethanol MPSP)

# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 18:10:51 2026

@author: sarangbhagwat
"""
import math
import biosteam as bst
from biorefineries import isobutanol
isobutanol.load()
from matplotlib import pyplot as plt

model = isobutanol.models.models_EtOH_IBO_corn.model
namespace_dict = isobutanol.models.namespace_dict
fbs_spec = isobutanol.models.fbs_spec
optimize_1D_feeding_strategy_for_MPSP = isobutanol.models.optimize_1D_feeding_strategy_for_MPSP
optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP =  isobutanol.models.optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP
optimize_max_n_glu_spikes_for_MPSP = isobutanol.models.optimize_max_n_glu_spikes_for_MPSP
plot_kinetic_results = isobutanol.models.plot_kinetic_results
solve_TEA = isobutanol.system.solve_TEA
unit_groups_dict = isobutanol.models.unit_groups_dict
model_specification = model.specification
f = model.system.flowsheet
V406 = f.V406
r = V406.nsk_kinetic_model
IBO_filepath = isobutanol.__file__.replace('\\__init__.py', '')

def _assert_MPSPs_stable(reference, current, sim_number, rel_tol=5e-3):
    """Verify the current simulation's MPSPs match the first simulation's to
    ~3 significant figures (relative drift < rel_tol; nan stays nan)."""
    for ID, ref in reference['MPSPs'].items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref) or math.isnan(cur):
            assert math.isnan(ref) and math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} vs first-sim {ref} (nan mismatch)'
        else:
            assert abs(cur - ref)/abs(ref) < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} drifted from first-sim '
                 f'value {ref} beyond rel tol {rel_tol}')

def load_simulate_baseline(scenario='B', # 'A' or 'B'
                           plot=False,
                           stream_IDs=('ethanol', 'isobutanol'), # products whose MPSPs are solved
                           IRR_for_MPSP=0.15, # fixed IRR at which MPSPs are solved
                           n_sims=3, # simulations in this kernel (in-process stability check)
                           ):
    parameter_distributions_filename = IBO_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+\
        f'parameter-distributions_corn_IBO_EtOH_{scenario}.xlsx'

    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
    baseline_initial = model.metrics_at_baseline()

    model_specification()

    # for forced batch mode:
    # fbs_spec.max_n_spikes = 0
    if scenario=='A':
        fbs_spec.max_n_spikes = 16
        sim_kwargs = dict(threshold_conc=217.125, target_conc=221.25)
    elif scenario=='B':
        fbs_spec.max_n_spikes = 13
        sim_kwargs = dict(threshold_conc=216.3, target_conc=226.3)
    else:
        raise ValueError(f'Scenario {scenario} not found.')

    # Load once (above), simulate n_sims times; each solve_TEA dict is
    # {'IRR': <solved IRR at default product prices>,
    #  'MPSPs': {ID: <purity-adjusted MPSP at IRR_for_MPSP>, ...}} (NaN for an
    # empty product, e.g. isobutanol in scenario A); no simulation inside
    # solve_TEA, prices restored. After each simulation the MPSPs are
    # verified stable against the first simulation's.
    all_results = []
    for i in range(n_sims):
        model_specification(**sim_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        if all_results:
            _assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)

    if plot:
        # Plot conc v time
        fig, ax = plot_kinetic_results(xlim=(0,80), ylim=(0,250))

    return all_results

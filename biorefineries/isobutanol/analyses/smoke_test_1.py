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

def _assert_MPSPs_pinned(expected, current, sim_number, rel_tol=0.01):
    """Verify the current simulation's MPSPs against the pinned baseline
    values (within rel_tol, 1% by default; nan stays nan). Re-pinned
    2026-09-02 for lang_factor=None (capital from each unit's own
    bare-module factors, no separate indirects; pins under Lang 3.0
    with 2023$ prices were A 0.98872 / B 1.4479 + 2.2448), earlier
    that day for the 2023 price year (every stream and utility price
    indexed to 2023$ with the BLS chemicals PPI, process_settings
    .index_prices_to_price_year, and the workbook price rows converted
    to 2023$; pins under the 2023 CEPCI with unindexed prices were
    A 0.85293 / B 1.0870 + 1.7961), then for the 2023
    CEPCI 797.9 (bst.CE had been biosteam's 567.5 default; pins under
    Lang 3.0 alone were A 0.76073 / B 0.81176 + 1.4618) and, before
    that, the Lang factor 3.0 TEA (Huang et al. 2016; the previous pins
    under corn's uncited Lang factor 4 were A 0.84057 / B 1.0521 +
    1.7538). Those 2026-09-01 pins followed (a) the
    kinetics-synced parameter xlsx (d467f0aa /
    3720fa21) and (b) the vent-scrubber-bottoms recycle to the separation
    feed (MX8) with molar L/G = 2.0 wash water; before that these tests
    asserted in-process stability only (A ~0.818 / B ~0.396 pre-sync)."""
    for ID, ref in expected.items():
        cur = current['MPSPs'][ID]
        if math.isnan(ref):
            assert math.isnan(cur), \
                f'sim {sim_number}: {ID} MPSP {cur} expected nan (empty product)'
        else:
            assert abs(cur - ref)/ref < rel_tol, \
                (f'sim {sim_number}: {ID} MPSP {cur} not within rel tol '
                 f'{rel_tol} of pinned {ref}')

def load_simulate_baseline(scenario='A', # 'A' or 'B'
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
        # both-trains scenario-A baseline (isobutanol product empty)
        expected_MPSPs = {'ethanol': 0.86604, 'isobutanol': math.nan}
    elif scenario=='B':
        fbs_spec.max_n_spikes = 13
        sim_kwargs = dict(threshold_conc=216.3, target_conc=226.3)
        # both-trains scenario-B baseline
        expected_MPSPs = {'ethanol': 1.0784, 'isobutanol': 1.7960}
    else:
        raise ValueError(f'Scenario {scenario} not found.')

    # Load once (above), simulate n_sims times; each solve_TEA dict is
    # {'IRR': <solved IRR at default product prices>,
    #  'MPSPs': {ID: <purity-adjusted MPSP at IRR_for_MPSP>, ...}} (NaN for an
    # empty product, e.g. isobutanol in scenario A); no simulation inside
    # solve_TEA, prices restored. After each simulation the MPSPs are
    # verified against the pinned baseline (1%) and stable against the
    # first simulation's.
    all_results = []
    for i in range(n_sims):
        model_specification(**sim_kwargs)
        results = solve_TEA(stream_IDs=stream_IDs, IRR_for_MPSP=IRR_for_MPSP)
        _assert_MPSPs_pinned(expected_MPSPs, results, i+1)
        if all_results:
            _assert_MPSPs_stable(all_results[0], results, i+1)
        all_results.append(results)

    if plot:
        # Plot conc v time
        fig, ax = plot_kinetic_results(xlim=(0,80), ylim=(0,250))

    return all_results

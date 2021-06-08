#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-2021, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.utils import TicToc
from biosteam.plots import plot_montecarlo_across_coordinate
from biorefineries.lactic import create_model, \
    SSCF_flowsheet, SSCF_funcs, SHF_flowsheet, SHF_funcs


# %%

def evaluate_uncertainties(kind='SSCF', seed=None, N_simulation=1000,
                           sampling_rule='L',
                           percentiles = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1],
                           if_plot=True, report_name='1_full_evaluation.xlsx'):
    if 'SSCF' in str(kind).upper():
        flowsheet = SSCF_flowsheet
        funcs = SSCF_funcs
    elif 'SHF' in str(kind).upper():
        flowsheet = SHF_flowsheet
        funcs = SHF_funcs
    else:
        raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')

    simulate_get_MPSP = funcs['simulate_get_MPSP']
    simulate_get_MPSP()
    model_dct = create_model(kind)
    model = model_dct['model']

    # =============================================================================
    # Evaluate and organize results for Monte Carlo analysis
    # =============================================================================

    # Initiate a timer
    timer = TicToc('timer')
    timer.tic()

    # Set seed to make sure each time the same set of random numbers will be used
    if seed:
        np.random.seed(seed)
    samples = model.sample(N=N_simulation, rule=sampling_rule)
    model.load_samples(samples)

    baseline_initial = model.metrics_at_baseline()
    baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]),
                            columns=baseline_initial.keys())

    model.evaluate()

    # Baseline results
    baseline_end = model.metrics_at_baseline()
    baseline = baseline.append(baseline_end, ignore_index=True)
    baseline.index = ('initial', 'end')
    # baseline.to_excel('0_baseline.xlsx')

    # Parameters
    parameters = model.get_parameters()
    index_parameters = len(model.get_baseline_sample())
    parameter_values = model.table.iloc[:, :index_parameters].copy()

    # TEA results
    index_TEA = index_parameters + model_dct['index_TEA']
    TEA_results = \
        model.table.iloc[:, index_parameters:index_TEA].copy()
    TEA_percentiles = TEA_results.quantile(q=percentiles)

    # MPSP across IRRs
    index_IRR = index_parameters + model_dct['index_IRR']
    IRR_results = \
        model.table.iloc[:, index_TEA:index_IRR].copy()
    IRR_percentiles = IRR_results.quantile(q=percentiles)

    # LCA_results
    LCA_results = \
        model.table.iloc[:, index_IRR::].copy()
    LCA_percentiles = LCA_results.quantile(q=percentiles)

    # Spearman's rank correlation
    spearman_metrics = model.metrics[0:2] + model.metrics[6:8] + \
        model.metrics[index_IRR:index_IRR+2]

    model.table = model.table.dropna()
    spearman_parameters = parameters
    spearman_results = model.spearman_r(spearman_parameters, spearman_metrics)[0]
    spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])

    # Calculate the cumulative probabilitie of each parameter
    probabilities = {}
    for i in range(index_parameters):
        p = parameters[i]
        p_values = parameter_values.iloc[:, 2*i]
        probabilities[p.name] = p.distribution.cdf(p_values)
        parameter_values.insert(loc=2*i+1,
                          column=(parameter_values.iloc[:, 2*i].name[0], 'Probability'),
                          value=probabilities[p.name],
                          allow_duplicates=True)

    run_number = samples.shape[0]

    # =============================================================================
    # Evaluate the min/max of one parameter each time to ensure the parameter can
    # independently affect the system
    # =============================================================================

    p_values = [[], [], []]
    MPSPs = [[], [], []]
    GWPs = [[], [], []]
    FECs = [[], [], []]

    bst.speed_up()
    for p in parameters:
        p_min = p.distribution.lower.item()
        p_max = p.distribution.upper.item()
        p_baseline = p.baseline
        p_value = (p_min, p_max, p_baseline)
        p.system = flowsheet.system.lactic_sys
        for i in range(len(p_value)):
            p.setter(p_value[i])
            p_values[i].append(p_value[i])
            MPSP = simulate_get_MPSP()
            MPSPs[i].append(MPSP)
            GWPs[i].append(funcs['get_GWP']())
            FECs[i].append(funcs['get_FEC']())
            run_number += 1

    MPSP_baseline = np.asarray(MPSPs[2])
    MPSP_min_diff = np.asarray(MPSPs[0]) - MPSP_baseline
    MPSP_max_diff = np.asarray(MPSPs[1]) - MPSP_baseline

    GWP_baseline = np.asarray(GWPs[2])
    GWP_min_diff = np.asarray(GWPs[0]) - GWP_baseline
    GWP_max_diff = np.asarray(GWPs[1]) - GWP_baseline

    FEC_baseline = np.asarray(FECs[2])
    FEC_min_diff = np.asarray(FECs[0]) - FEC_baseline
    FEC_max_diff = np.asarray(FECs[1]) - FEC_baseline

    one_p_df = pd.DataFrame({
        ('Parameter', 'Name'): [i.name_with_units for i in parameters],
        ('Parameter', 'Baseline'): p_values[2],
        ('Parameter', 'Min'): p_values[0],
        ('Parameter', 'Max'): p_values[1],
        ('MPSP [$/kg]', 'MPSP baseline'): MPSP_baseline,
        ('MPSP [$/kg]', 'MPSP min'): MPSPs[0],
        ('MPSP [$/kg]', 'MPSP min diff'): MPSP_min_diff,
        ('MPSP [$/kg]', 'MPSP max'): MPSPs[1],
        ('MPSP [$/kg]', 'MPSP max diff'): MPSP_max_diff,
        ('GWP [kg CO2-eq/kg]', 'GWP baseline'): GWP_baseline,
        ('GWP [kg CO2-eq/kg]', 'GWP min'): GWPs[0],
        ('GWP [kg CO2-eq/kg]', 'GWP min diff'): GWP_min_diff,
        ('GWP [kg CO2-eq/kg]', 'GWP max'): GWPs[1],
        ('GWP [kg CO2-eq/kg]', 'GWP max diff'): GWP_max_diff,
        ('FEC [MJ/kg]', 'FEC baseline'): FEC_baseline,
        ('FEC [MJ/kg]', 'FEC min'): FECs[0],
        ('FEC [MJ/kg]', 'FEC min diff'): FEC_min_diff,
        ('FEC [MJ/kg]', 'FEC max'): FECs[1],
        ('FEC [MJ/kg]', 'FEC max diff'): FEC_max_diff,
        })

    time = timer.elapsed_time / 60
    print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


    if if_plot:
        IRR_plot_indices = [metric.index for metric in model.metrics
                            if 'IRR' in metric.index[0] and 'MPSP' in metric.index[1]]
        IRR_plot_data = IRR_results[IRR_plot_indices].copy()
        IRR_plot_data.columns = model_dct['IRRs'].copy()
        IRR_plot_y = IRR_plot_data.sort_index(axis=1)
        IRR_plot_y = IRR_plot_y.dropna()
        IRR_plot_x = model_dct['IRRs'].copy()
        IRR_plot_x.sort()
        plot_montecarlo_across_coordinate(IRR_plot_x, IRR_plot_y)

    if report_name:
        with pd.ExcelWriter(report_name) as writer:
            parameter_values.to_excel(writer, sheet_name='Parameters')
            TEA_results.to_excel(writer, sheet_name='TEA results')
            TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
            IRR_results.to_excel(writer, sheet_name='IRR results')
            IRR_percentiles.to_excel(writer, sheet_name='IRR percentiles')
            LCA_results.to_excel(writer, sheet_name='LCA results')
            LCA_percentiles.to_excel(writer, sheet_name='LCA percentiles')
            spearman_results.to_excel(writer, sheet_name='Spearman')
            one_p_df.to_excel(writer, sheet_name='One-parameter')
            model.table.to_excel(writer, sheet_name='Raw data')

evaluate_uncertainties(kind='SSCF', seed=None, N_simulation=10, sampling_rule='L',
                        percentiles = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1],
                        if_plot=True, report_name='1_full_evaluation.xlsx')
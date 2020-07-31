#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon Apr 13 10:24:42 2020

Modified from the biorefineries constructed in [1] and [2] for the production of
lactic acid from lignocellulosic feedstocks

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted.
    July, 2020.

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.utils import TicToc
from biosteam.plots import plot_montecarlo_across_coordinate
from lactic.system import lactic_sys, simulate_get_MPSP
from lactic.analyses import models

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


# %%

# =============================================================================
# Evaluate and organize results for Monte Carlo analysis
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

model = models.model_full

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221)
N_simulation = 1000 # 1000
samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)

baseline_initial = model.metrics_at_baseline()
baseline_df = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                           columns=baseline_initial.keys())

model.evaluate()

# Baseline results
# model._system.TEA.IRR = 0.1
baseline_end = model.metrics_at_baseline()
baseline_df = baseline_df.append(baseline_end, ignore_index=True)
baseline_df.index = ('initial', 'end')
baseline_df.to_excel('0_baseline.xlsx')

# Parameters
parameters = model.get_parameters()
parameter_len = len(model.get_baseline_sample())
parameter_values = model.table.iloc[:, :parameter_len].copy()

# Add baseline values to the end
parameter_values.loc['baseline'] = model.get_baseline_sample()

# All results
all_results = model.table.iloc[:, parameter_len::].copy()
all_results.loc['baseline initial'] = baseline_initial
all_results.loc['baseline end'] = baseline_end

# Monte Carlo results
metrics_no_IRRs_len = len(models.metrics_no_IRRs)
Monte_Carlo_results = \
    model.table.iloc[:, parameter_len:parameter_len+metrics_no_IRRs_len].copy()
Monte_Carlo_percentiles = Monte_Carlo_results.quantile(q=percentiles)

# MPSP across IRRs
IRR_results = \
    model.table.iloc[:, parameter_len+metrics_no_IRRs_len::].copy()
IRR_percentiles = IRR_results.quantile(q=percentiles)

# Spearman's rank correlation
spearman_metrics = model.metrics[0:3] + model.metrics[6:9]
spearman_parameters = parameters
spearman_results = model.spearman(spearman_parameters, spearman_metrics)
spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])

# Calculate the cumulative probabilitie of each parameter
probabilities = {}
for i in range(parameter_len):
    p = parameters[i]
    p_values = parameter_values.iloc[:, 2*i]
    probabilities[p.name] = p.distribution.cdf(p_values)
    parameter_values.insert(loc=2*i+1, 
                      column=(parameter_values.iloc[:, 2*i].name[0], 'Probability'), 
                      value=probabilities[p.name],
                      allow_duplicates=True)

run_number = samples.shape[0]


# %%

# =============================================================================
# Evaluate the min/max of one parameter each time to ensure the parameter can
# independently affect the system
# =============================================================================

p_values = [[], [], []]
MPSPs = [[], [], []]

bst.speed_up()
for p in parameters:
    [p_min], [p_max] = p.distribution.range().tolist()
    p_baseline = p.baseline
    p_value = (p_min, p_max, p_baseline)
    p.system = lactic_sys
    for i in range(len(p_value)):
        p.setter(p_value[i])
        p_values[i].append(p_value[i])
        MPSP = simulate_get_MPSP()
        MPSPs[i].append(MPSP)
        run_number += 1

MPSP_baseline = np.asarray(MPSPs[2])
MPSP_min_diff = np.asarray(MPSPs[0])-MPSP_baseline
MPSP_max_diff = np.asarray(MPSPs[1])-MPSP_baseline

one_p_df = pd.DataFrame({
    'Parameter': [i.name_with_units for i in parameters],
    'Minimum': p_values[0],
    'Maximum': p_values[1],
    'Baseline': p_values[2],
    'Min-value MPSP [$/kg]': MPSPs[0],
    'Max-value MPSP [$/kg]': MPSPs[1],
    'Baseline MPSP [$/kg]': MPSPs[2],
    'Min-value MPSP change [$/kg]': MPSP_min_diff,
    'Max-value MPSP change [$/kg]': MPSP_max_diff,
    'Min-value MPSP change [%]': MPSP_min_diff/MPSP_baseline,
    'Max-value MPSP change [%]': MPSP_max_diff/MPSP_baseline
    })

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %%

'''Get a quick plot'''
IRR_plot_indices = [metric.index for metric in model.metrics
                    if 'IRR' in metric.index[0] and 'MPSP' in metric.index[1]]
IRR_plot_data = IRR_results[IRR_plot_indices].copy()
IRR_plot_data.columns = models.IRRs.copy()
IRR_plot_y = IRR_plot_data.sort_index(axis=1)
IRR_plot_x = models.IRRs.copy()
IRR_plot_x.sort()
plot_montecarlo_across_coordinate(IRR_plot_x, IRR_plot_y.values)

'''Output to Excel'''
with pd.ExcelWriter('1_full_evaluation.xlsx') as writer:
    parameter_values.to_excel(writer, sheet_name='Parameters')
    Monte_Carlo_results.to_excel(writer, sheet_name='Monte Carlo results')
    Monte_Carlo_percentiles.to_excel(writer, sheet_name='Monte Carlo percentiles')
    IRR_results.to_excel(writer, sheet_name='IRR results')
    IRR_percentiles.to_excel(writer, sheet_name='IRR percentiles')
    spearman_results.to_excel(writer, sheet_name='Spearman')
    one_p_df.to_excel(writer, sheet_name='One-parameter')
    model.table.to_excel(writer, sheet_name='Raw data')


# %%

# =============================================================================
# Temporary codes for debugging
# =============================================================================

# import numpy as np
# from lactic.analyses import models

# model = models.model_full
# np.random.seed(3221)

# try: parameter = [parameters[-2]]
# except: parameter = parameters
# model.set_parameters(parameter)

# N_simulation = 20 # 1000
# # samples = model.sample(N=N_simulation, rule='L')
# samples = np.ones(N_simulation)
# samples *= 0.070

# model.load_samples(samples)
# model.evaluate()

# model.table.to_clipboard()









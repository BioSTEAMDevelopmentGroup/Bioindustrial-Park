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
from biosteam.utils import TicToc
from lactic import models

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


# %% 

# =============================================================================
# Evaluate and calculate probabilities
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

model = models.model_full

'''Quick look at baseline values'''
# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221)
N_simulation = 10 # 1000
samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)
baseline = model.metrics_at_baseline()
baseline_df = pd.DataFrame(data=np.array([[i for i in baseline.values],]), 
                            index=('baseline',), columns=baseline.keys())
# baseline_df.to_excel('baseline.xlsx')

'''Full evaluation'''
model.evaluate()
# Parameters and probabilities
parameter_len = len(model.get_baseline_sample())
parameters = model.table.iloc[:, :parameter_len].copy()
# Add baseline values to the end
parameters.loc['baseline'] = model.get_baseline_sample()

Monte_Carlo_results = model.table.iloc[:, parameter_len::].copy()
Monte_Carlo_percentiles = Monte_Carlo_results.quantile(q=percentiles)
# Add baseline values to the end
Monte_Carlo_results.loc['baseline'] = model.metrics_at_baseline()

# Note that if only one metric is used, then need to make sure it's a tuple
spearman_metrics = model.metrics[0:4]
spearman_results = model.spearman(spearman_metrics)

# Calculate the probabilities of each parameter and the overall scenario
probabilities = {}
for i in range(parameter_len):
    p = model.get_parameters()[i]
    p_values = parameters.iloc[:, 2*i]
    #!!! cdf vs. 1-cdf, needs to be mannually reviewed
    if spearman_results.iloc[:, 0][i]>0:
        probabilities[p.name] = p.distribution.cdf(p_values)
    else:
        probabilities[p.name] = np.ones(len(p_values))-p.distribution.cdf(p_values)   
    parameters.insert(loc=2*i+1, 
                      column=(parameters.iloc[:, 2*i].name[0], 'Probability'), 
                      value=probabilities[p.name],
                      allow_duplicates=True)

# '''Output to Excel'''
# with pd.ExcelWriter('1_full_evaluation.xlsx') as writer:
#     parameters.to_excel(writer, sheet_name='Parameters')
#     Monte_Carlo_results.to_excel(writer, sheet_name='Monte Carlo')
#     Monte_Carlo_percentiles.to_excel(writer, sheet_name='Monte Carlo Percentiles')
#     spearman_results.to_excel(writer, sheet_name='Spearman')
#     model.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation
time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %%

# =============================================================================
# Evaluate across feedstock succinic acid content
# =============================================================================

# # Initiate a timer
# timer_across_succinic = TicToc('timer_across_succinic')
# timer_across_succinic.tic()

# from chaospy import distributions as shape
# model_succinic = models.model_succinic
# set_succinic = models.set_succinic
# R301.set_titer_limit = False
# R301.neutralization = True
# R401.bypass = False
# S402.bypass = False

# '''Evaluate'''
# np.random.seed(3221)
# N_simulation = 2 # 1000
# succinic_samples = model_succinic.sample(N=N_simulation, rule='L')
# model_succinic.load_samples(succinic_samples)

# succnic_coordinate = np.linspace(0, 0.05, 5+1) # step = 0.01
# succnic_distribution = shape.Uniform(0, 0.05)
# succnic_probability = np.ones(len(succnic_coordinate)) - \
#     succnic_distribution.cdf(succnic_coordinate)

# succinic_data = model_succinic.evaluate_across_coordinate(
#     'Feedstock succinic acid content', set_succinic, succnic_coordinate, notify=True)

# columns = pd.MultiIndex.from_arrays([succnic_coordinate, succnic_probability],
#                                     names=['Feedstock succinic acid content', 
#                                            'Probability'])

# succinic_MPSP = pd.DataFrame(
#     data=succinic_data[('Biorefinery', 'Minimum product selling price [$/kg]')],
#     columns=columns)

# succinic_MPSP_percentiles = succinic_MPSP.quantile(q=percentiles)

# # Frequency of simulation that has MPSP < target price
# succinic_MPSP_f_upper = succinic_MPSP[succinic_MPSP<MPSP_target_upper].count()/N_simulation
# succinic_MPSP_f_upper.name = 'frequency_upper'
# succinic_MPSP_f_lower = succinic_MPSP[succinic_MPSP<MPSP_target_lower].count()/N_simulation
# succinic_MPSP_f_lower.name = 'frequency_lower'
# succinic_MPSP = succinic_MPSP.append(succinic_MPSP_f_upper)
# succinic_MPSP = succinic_MPSP.append(succinic_MPSP_f_lower)

# '''To get a quick plot'''
# #!!! Plot not looking good, consider removing
# from biosteam.plots import plot_montecarlo_across_coordinate
# succinic_plot_index = ('Biorefinery', 'Minimum product selling price [$/kg]')
# succinic_plot_data = succinic_data[succinic_plot_index]
# plot_montecarlo_across_coordinate(succnic_coordinate, succinic_plot_data)

# # '''Output to Excel'''
# # with pd.ExcelWriter('5_across_feedstock_succinic_acid_content.xlsx') as writer:
# #     succinic_MPSP.to_excel(writer, sheet_name='Succnic')
# #     succinic_MPSP_percentiles.to_excel(writer, sheet_name='Succnic percentiles')
# #     model_succinic.table.to_excel(writer, sheet_name='Raw data')

# run_number = N_simulation * len(succnic_coordinate)
# time = timer_across_succinic.elapsed_time / 60
# print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')









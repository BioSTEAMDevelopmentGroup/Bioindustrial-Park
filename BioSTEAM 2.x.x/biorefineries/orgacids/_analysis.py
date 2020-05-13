#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:24:42 2020

@author: yalinli_cabbi
"""

'''
TODO:
    NPV not converging... not sure if it's because of the not converging
    
    Consider evaluating across the coordinate of:
        Heavy organic acid content in feedstock/fermentation broth
        Purity of lactic acid product
    (should use evaluate_across_coordinate if the variable affects system simulation)
'''


# %% Setup

import numpy as np
import pandas as pd
from orgacids.model import orgacids_model, orgacids_model_IRR

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]



# np.random.seed(3221)
# samples = orgacids_model_IRR.sample(N=10, rule='L')
# parameter_len = len(orgacids_model_IRR.get_baseline_sample())

# orgacids_model_IRR.load_samples(samples)
# orgacids_model_IRR.evaluate()

# IRR_results = orgacids_model_IRR.table.iloc[:, parameter_len::].copy()
# IRR_percentiles = IRR_results.quantile(q=percentiles)


# print(IRR_results)



# %% Perform Monte Carlo and calculate scenario probabilities

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221)
samples = orgacids_model.sample(N=10, rule='L')
orgacids_model.load_samples(samples)
orgacids_model.evaluate()

# Parameters and probabilities
parameter_len = len(orgacids_model.get_baseline_sample())
parameters = orgacids_model.table.iloc[:, :parameter_len].copy()
# Add baseline values to the end
parameters.loc['baseline'] = orgacids_model.get_baseline_sample()

Monte_Carlo_results = orgacids_model.table.iloc[:, parameter_len::].copy()
Monte_Carlo_percentiles = Monte_Carlo_results.quantile(q=percentiles)
# Add baseline values to the end
Monte_Carlo_results.loc['baseline'] = orgacids_model.metrics_at_baseline()

# Note that if only one metric is used, then need to make sure it's a tuple
spearman_metrics = orgacids_model.metrics[0:4]
spearman_results = orgacids_model.spearman(spearman_metrics)

# Calculate the probabilities of each parameter and the overall scenario
probabilities = {}
for i in range(parameter_len):
    p = orgacids_model.get_parameters()[i]
    p_values = parameters.iloc[:, 2*i]
    #!!! cdf vs. 1-cdf, might need being mannually reviewed
    if spearman_results.iloc[:, 0][i]>0:
        probabilities[p.name] = p.distribution.cdf(p_values)
    else:
        probabilities[p.name] = np.ones(len(p_values))-p.distribution.cdf(p_values)
    parameters.insert(loc=2*i+1, 
                      column=(parameters.iloc[:, 2*i].name[0], 'Probability'), 
                      value=probabilities[p.name],
                      allow_duplicates=True)

# Calculate the overall probabilities of scenario
probability_list = list(i for i in probabilities.values())
scenario_probability = probability_list[0].copy()
for i in range(parameter_len-1):
    scenario_probability *= list(i for i in probabilities.values())[i+1]
probabilities['scenario_abs'] = scenario_probability

# Normalize scenario probabilities by baseline probability
scenario_probability_normalzied = probabilities['scenario_abs'].copy()
scenario_probability_normalzied /= \
    np.ones(len(probabilities['scenario_abs']))*probabilities['scenario_abs'][-1]
probabilities['scenario_normalized'] = scenario_probability_normalzied
parameters.insert(loc=parameters.shape[1],
                  column=('Scenario', 'Absolute probability'), 
                  value=probabilities['scenario_abs'])
parameters.insert(loc=parameters.shape[1],
                  column=('Scenario', 'Normalized probability'), 
                  value=probabilities['scenario_normalized'])

# Add scenario probabilities to Monte Carlo results
Monte_Carlo_results.insert(loc=0, column=('Scenario', 'Normalized probability'),
                            value=parameters.iloc[:, parameters.shape[1]-1])

# with pd.ExcelWriter('Monte Carlo scenarios and probabilities.xlsx') as writer:
#     parameters.to_excel(writer, sheet_name='Parameters')
#     Monte_Carlo_results.to_excel(writer, sheet_name='Monte Carlo')
#     Monte_Carlo_percentiles.to_excel(writer, sheet_name='Monte Carlo Percentiles')
#     spearman_results.to_excel(writer, sheet_name='Spearman')


# %% Evaluated resutls across the IRR coordinate

# Use the same set of samples as orgacids_model
orgacids_model_IRR.load_samples(samples)
orgacids_model_IRR.evaluate()

IRR_results = orgacids_model_IRR.table.iloc[:, parameter_len::].copy()
IRR_percentiles = IRR_results.quantile(q=percentiles)

# # To get a quick plot
# from biosteam.evaluation.evaluation_tools import plot_montecarlo_across_coordinate
# from orgacids.model import IRRs
# MSP_indices = [metric.index for metric in orgacids_model_IRR.metrics][0:-1]
# MSP_data = orgacids_model_IRR.table[MSP_indices]
# plot_montecarlo_across_coordinate(IRRs, MSP_data.values)

#!!! If having multiple coordinates, consider summarizing all percentiles into
# one df and make percentiles as the columns, make multi-level column labels

# with pd.ExcelWriter('Monte Carlo across coordinates.xlsx') as writer:
#     IRR_results.to_excel(writer, sheet_name='IRR')
#     IRR_percentiles.to_excel(writer, sheet_name='IRR Percentiles')


# %% Reference codes for using evaluate_across_coordinate

# from orgacids.model import orgacids_model_IRR
# orgacids_model.load_samples(samples)
# orgacids_model.evaluate()

# def set_IRR(IRR):
#     orgacids_tea.IRR = IRR

# orgacids_model_IRR.load_samples(samples)
# np.random.seed(3221)
# coordinate = np.linspace(0.01, 0.4, 40)

# orgacids_model_IRR.evaluate_across_coordinate('IRR', set_IRR, coordinate,
#                                               xlfile='Monte Carlo across IRR3.xlsx',
#                                               notify=False)














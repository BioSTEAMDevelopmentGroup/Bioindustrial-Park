#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 10:24:42 2020

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup (required by all models)
# =============================================================================

import numpy as np
import pandas as pd
from biosteam.utils import TicToc

# Based on commercial lactic acid (88%) price in IHS Markit report
MSP_target_upper = 2.1
MSP_target_lower = 1.7

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


# %% 

# =============================================================================
# Evaluate and calculate probabilities
# =============================================================================

# Initiate a timer
timer_full_evaluation = TicToc('timer_full_evaluation')
timer_full_evaluation.tic()

from TAL.models import TAL_model

'''Quick look at baseline values'''
# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221)
N_simulation = 100 # 1000
samples = TAL_model.sample(N=N_simulation, rule='L')
TAL_model.load_samples(samples)
baseline = TAL_model.metrics_at_baseline()
baseline_df = pd.DataFrame(data=np.array([[i for i in baseline.values()],]), 
                            index=('baseline',), columns=baseline.keys())
baseline_df.to_excel('baseline_88.xlsx')

'''Full evaluation'''
TAL_model.evaluate()
# Parameters and probabilities
parameter_len = len(TAL_model.get_baseline_sample())
parameters = TAL_model.table.iloc[:, :parameter_len].copy()
# Add baseline values to the end
parameters.loc['baseline'] = TAL_model.get_baseline_sample()

Monte_Carlo_results = TAL_model.table.iloc[:, parameter_len::].copy()
Monte_Carlo_percentiles = Monte_Carlo_results.quantile(q=percentiles)
# Add baseline values to the end
Monte_Carlo_results.loc['baseline'] = TAL_model.metrics_at_baseline()

# Note that if only one metric is used, then need to make sure it's a tuple
spearman_metrics = TAL_model.metrics[0:6]
spearman_results = TAL_model.spearman(spearman_metrics)

# Calculate the probabilities of each parameter and the overall scenario
probabilities = {}
for i in range(parameter_len):
    p = TAL_model.get_parameters()[i]
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

'''Output to Excel'''
with pd.ExcelWriter('E1_baseline_evaluation_88.xlsx') as writer:
    parameters.to_excel(writer, sheet_name='Parameters')
    Monte_Carlo_results.to_excel(writer, sheet_name='Monte Carlo')
    Monte_Carlo_percentiles.to_excel(writer, sheet_name='Monte Carlo Percentiles')
    spearman_results.to_excel(writer, sheet_name='Spearman')
    TAL_model.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation
time = timer_full_evaluation.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %% 

# =============================================================================
# Evaluating across internal rate of return
# =============================================================================

# Initiate a timer
timer_IRR = TicToc('timer_IRR')
timer_IRR.tic()

from TAL.models import TAL_model_IRR

'''
Note:
    Not using `evaluate_across_coordinate` function as changing IRR
    does not affect the system, using IRR as the metrics for evaluation 
    will save considerable time.
'''

'''Evaluate'''
np.random.seed(3221)
N_simulation = 100 # 1000
samples = TAL_model_IRR.sample(N=N_simulation, rule='L')
TAL_model_IRR.load_samples(samples)
TAL_model_IRR.evaluate()

parameter_len = len(TAL_model_IRR.get_baseline_sample())
IRR_results = TAL_model_IRR.table.iloc[:, parameter_len::].copy()
IRR_percentiles = IRR_results.quantile(q=percentiles)

'''To get a quick plot'''
from biosteam.evaluation.evaluation_tools import plot_montecarlo_across_coordinate
from TAL.models import IRRs
MSP_indices = [metric.index for metric in TAL_model_IRR.metrics
                if 'Net present value' not in metric.index[1]]
MSP_data = TAL_model_IRR.table[MSP_indices]
plot_montecarlo_across_coordinate(IRRs, MSP_data.values)

'''Output to Excel'''
with pd.ExcelWriter('E2_across_IRR_pure.xlsx') as writer:
    IRR_results.to_excel(writer, sheet_name='IRR')
    IRR_percentiles.to_excel(writer, sheet_name='IRR Percentiles')
    TAL_model_IRR.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation
time = timer_IRR.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %% 

# =============================================================================
# Evaluate across lactic acid yield
# =============================================================================

# Initiate a timer
timer_across_LA = TicToc('timer_across_LA')
timer_across_LA.tic()

from chaospy import distributions as shape
from TAL.models import TAL_model_LA_yield, set_LA_yield

'''Evaluate'''
np.random.seed(3221)
N_simulation = 20 # 1000
LA_yield_samples = TAL_model_LA_yield.sample(N=N_simulation, rule='L')
TAL_model_LA_yield.load_samples(LA_yield_samples)

LA_yield_coordinate = np.linspace(0.5, 1, 11) # step = 0.05
LA_yield_distribution = shape.Triangle(0.55, 0.76, 0.93)
LA_yield_probability = np.ones(len(LA_yield_coordinate)) - \
    LA_yield_distribution.cdf(LA_yield_coordinate)

LA_yield_data = TAL_model_LA_yield.evaluate_across_coordinate(
    'Lactic acid yield', set_LA_yield, LA_yield_coordinate, notify=True)

columns = pd.MultiIndex.from_arrays([LA_yield_coordinate, LA_yield_probability],
                                    names=['Lactic acid yield', 'Probability'])

LA_yield_MSP = pd.DataFrame(
    data=LA_yield_data[('Biorefinery', 'Minimum selling price [$/kg]')],
    columns=columns)

LA_yield_MSP_percentiles = LA_yield_MSP.quantile(q=percentiles)

# Frequency of simulation that has MSP < target price
LA_yield_MSP_f_upper = LA_yield_MSP[LA_yield_MSP<MSP_target_upper].count()/N_simulation
LA_yield_MSP_f_upper.name = 'frequency_upper'
LA_yield_MSP_f_lower = LA_yield_MSP[LA_yield_MSP<MSP_target_lower].count()/N_simulation
LA_yield_MSP_f_lower.name = 'frequency_lower'
LA_yield_MSP = LA_yield_MSP.append(LA_yield_MSP_f_upper)
LA_yield_MSP = LA_yield_MSP.append(LA_yield_MSP_f_lower)

'''Output to Excel'''
with pd.ExcelWriter('E3_across_lactic_acid_yield_88.xlsx') as writer:
    LA_yield_MSP.to_excel(writer, sheet_name='LA yield')
    LA_yield_MSP_percentiles.to_excel(writer, sheet_name='LA yield percentiles')
    TAL_model_LA_yield.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation * len(LA_yield_coordinate)
time = timer_across_LA.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %% 

# =============================================================================
# Evaluate across feedstock price and carbohydrate content
# =============================================================================

# Initiate a timer
timer_across_carbs = TicToc('timer_across_carbs')
timer_across_carbs.tic()

from TAL.models import TAL_model_feedstock, set_feedstock_carbs, prices

'''Evaluate'''
np.random.seed(3221)
# This is not a Monte Carlo simulation, this evaluation uses the baseline parameters
# to see the impacts of feedstock carbohydrate content
# The parameter is a fake one to enable the evaluation
N_simulation = 1
feedstock_samples_1d = TAL_model_feedstock.sample(N=N_simulation, rule='L')
feedstock_samples = feedstock_samples_1d[:, np.newaxis]
TAL_model_feedstock.load_samples(feedstock_samples)

feedstock_carbs_coordinate = np.linspace(0.4, 0.7, 31) # step = 0.01

feedstock_data = TAL_model_feedstock.evaluate_across_coordinate(
    'Feedstock carbohydate content', set_feedstock_carbs, 
    feedstock_carbs_coordinate, notify=True)

feedstock_MSP = pd.DataFrame({
    ('Parameter','Carbohydrate content [dry mass %]'): feedstock_carbs_coordinate
    })

for (i, j) in zip(feedstock_data.keys(), feedstock_data.values()):
    feedstock_MSP[i] = j[0]

'''Organize data for easy plotting'''
x_axis = [f'{i:.3f}' for i in feedstock_carbs_coordinate]
x_axis *= len(prices)
y_axis = sum(([f'{i:.0f}']*len(feedstock_carbs_coordinate) for i in prices), [])

MSP = []
NPV = []
for i in range(feedstock_MSP.columns.shape[0]):
    if 'Minimum selling price' in feedstock_MSP.columns[i][1]:
        MSP +=  feedstock_MSP[feedstock_MSP.columns[i]].to_list()
    if 'Net present value' in feedstock_MSP.columns[i][1]:
        NPV +=  feedstock_MSP[feedstock_MSP.columns[i]].to_list()

feedstock_MSP_plot = pd.DataFrame()
feedstock_MSP_plot['Carbohydrate content [dry mass %]'] = x_axis
feedstock_MSP_plot['Price [$/dry-ton]'] = y_axis
feedstock_MSP_plot['Minimum selling price [$/kg]'] = MSP
feedstock_MSP_plot['Net present value [$]'] = NPV

'''Output to Excel'''
name = 'E4_across_feedstock_cost_and_carbohydrate_content_88.xlsx'
with pd.ExcelWriter(name) as writer:
    feedstock_MSP.to_excel(writer, sheet_name='Evaluation data')
    feedstock_MSP_plot.to_excel(writer, sheet_name='For plotting')
    TAL_model_feedstock.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation * len(feedstock_carbs_coordinate)
time = timer_across_carbs.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %%

# =============================================================================
# Evaluate across feedstock succinic acid content
# =============================================================================

# Initiate a timer
timer_across_SA = TicToc('timer_across_SA')
timer_across_SA.tic()

from chaospy import distributions as shape
from TAL.models import TAL_model_SA_content, set_feedstock_succinic_acid_content

'''Evaluate'''
np.random.seed(3221)
N_simulation = 30 # 1000
SA_content_samples = TAL_model_SA_content.sample(N=N_simulation, rule='L')
TAL_model_SA_content.load_samples(SA_content_samples)

SA_content_coordinate = np.linspace(0, 0.1, 5+1) # step = 0.01
SA_content_distribution = shape.Uniform(0, 0.05)
SA_content_probability = np.ones(len(SA_content_coordinate)) - \
    SA_content_distribution.cdf(SA_content_coordinate)

SA_content_data = TAL_model_SA_content.evaluate_across_coordinate(
    'Feedstock succinic acid content', set_feedstock_succinic_acid_content, 
    SA_content_coordinate, notify=True)

columns = pd.MultiIndex.from_arrays([SA_content_coordinate, SA_content_probability],
                                    names=['Feedstock succinic acid content', 
                                           'Probability'])

SA_content_MSP = pd.DataFrame(
    data=SA_content_data[('Biorefinery', 'Minimum selling price [$/kg]')],
    columns=columns)

SA_content_MSP_percentiles = SA_content_MSP.quantile(q=percentiles)

# Frequency of simulation that has MSP < target price
SA_content_MSP_f_upper = SA_content_MSP[SA_content_MSP<MSP_target_upper].count()/N_simulation
SA_content_MSP_f_upper.name = 'frequency_upper'
SA_content_MSP_f_lower = SA_content_MSP[SA_content_MSP<MSP_target_lower].count()/N_simulation
SA_content_MSP_f_lower.name = 'frequency_lower'
SA_content_MSP = SA_content_MSP.append(SA_content_MSP_f_upper)
SA_content_MSP = SA_content_MSP.append(SA_content_MSP_f_lower)

'''To get a quick plot'''
#!!! The plot doesn't look great for SA though
from biosteam.evaluation.evaluation_tools import plot_montecarlo_across_coordinate
MSP_indx = ('Biorefinery', 'Minimum selling price [$/kg]')
MSP_data = SA_content_data[MSP_indx]
plot_montecarlo_across_coordinate(SA_content_coordinate, MSP_data)

'''Output to Excel'''
with pd.ExcelWriter('E5_across_feedstock_succinic_acid_content_88.xlsx') as writer:
    SA_content_MSP.to_excel(writer, sheet_name='SA content')
    SA_content_MSP_percentiles.to_excel(writer, sheet_name='SA content percentiles')
    TAL_model_SA_content.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation * len(SA_content_coordinate)
time = timer_across_SA.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')




# %%

# =============================================================================
# Evaluate across HXN minimum approach temperature
# =============================================================================

# Initiate a timer
timer_across_T_min_app = TicToc('timer_across_T_min_app')
timer_across_T_min_app.tic()

from chaospy import distributions as shape
from TAL.models import TAL_model_HXN_T_min_app, set_HXN_T_min_app

'''Evaluate'''
np.random.seed(3221)
N_simulation = 60 # 1000
T_min_app_samples = TAL_model_HXN_T_min_app.sample(N=N_simulation, rule='L')
TAL_model_HXN_T_min_app.load_samples(T_min_app_samples)

T_min_app_coordinate = np.linspace(1, 10, 60) 
T_min_app_distribution = shape.Uniform(1, 10)
T_min_app_probability = np.ones(len(T_min_app_coordinate)) - \
    T_min_app_distribution.cdf(T_min_app_coordinate)

T_min_app_data = TAL_model_HXN_T_min_app.evaluate_across_coordinate(
    'HXN minimum approach temperature', set_HXN_T_min_app, 
    T_min_app_coordinate, notify=True)

columns = pd.MultiIndex.from_arrays([T_min_app_coordinate, T_min_app_probability],
                                    names=['HXN minimum approach temperature', 
                                           'Probability'])

T_min_app_MSP = pd.DataFrame(
    data=T_min_app_data[('Biorefinery', 'Minimum selling price [$/kg]')],
    columns=columns)

T_min_app_MSP_percentiles = T_min_app_MSP.quantile(q=percentiles)

# Frequency of simulation that has MSP < target price
T_min_app_MSP_f_upper = T_min_app_MSP[T_min_app_MSP<MSP_target_upper].count()/N_simulation
T_min_app_MSP_f_upper.name = 'frequency_upper'
T_min_app_MSP_f_lower = T_min_app_MSP[T_min_app_MSP<MSP_target_lower].count()/N_simulation
T_min_app_MSP_f_lower.name = 'frequency_lower'
T_min_app_MSP = T_min_app_MSP.append(T_min_app_MSP_f_upper)
T_min_app_MSP = T_min_app_MSP.append(T_min_app_MSP_f_lower)

'''To get a quick plot'''

from biosteam.evaluation.evaluation_tools import plot_montecarlo_across_coordinate
# MSP_indx = ('Biorefinery', 'Minimum selling price [$/kg]')
# MSP_data = T_min_app_data[MSP_indx]

sep_util_indx = ('Utility cost', 'separation_sys [10^6 $/yr]')
HXN_util_indx = ('Utility cost', 'HXN [10^6 $/yr]')
sep_util_data = T_min_app_data[sep_util_indx]
HXN_util_data = T_min_app_data[HXN_util_indx]
util_savings_data = - 100* (HXN_util_data)/(sep_util_data - HXN_util_data)

plot_montecarlo_across_coordinate(T_min_app_coordinate, util_savings_data)

'''Output to Excel'''
with pd.ExcelWriter('E6_across_T_min_app_88.xlsx') as writer:
    T_min_app_MSP.to_excel(writer, sheet_name='T_min_app')
    T_min_app_MSP_percentiles.to_excel(writer, sheet_name='T_min_app percentiles')
    TAL_model_HXN_T_min_app.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation * len(T_min_app_coordinate)
time = timer_across_T_min_app.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')








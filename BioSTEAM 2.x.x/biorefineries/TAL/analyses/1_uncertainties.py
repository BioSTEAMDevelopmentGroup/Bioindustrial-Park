#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Sarang Bhagwat <sarangb2@illinois.edu>
# Yoel Cortes-Pena <yoelcortes@gmail.com>, and Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon Apr 13 10:24:42 2020

Modified from the biorefineries constructed in [1], [2], and [3] for the production of
[1] 3-hydroxypropionic acid, [2] lactic acid, and [3] ethanol from lignocellulosic feedstocks

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
    https://doi.org/10.1021/acssuschemeng.9b07040
    

@author: sarangbhagwat
"""


# %% 

# =============================================================================
# Setup
# =============================================================================
from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
# import biosteam as bst
from biosteam.utils import TicToc
from biorefineries.TAL.system_TAL_adsorption_glucose import (
    spec, TAL_sys,
)
from biorefineries.TAL.analyses import models
from datetime import datetime
import os

# R301 = flowsheet('R301')
percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


# %%

# =============================================================================
# Evaluate and organize results for Monte Carlo analysis
# =============================================================================
# spec.load_specifications(spec_1=0.49, spec_2=54.8, spec_3=0.76)
# Initiate a timer
timer = TicToc('timer')
timer.tic()

# model = models.model_full
model = models.TAL_model
# R301.set_titer_limit = True

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3221) # 3221
N_simulation = 2000 # 2000

samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)


###############################
# Bugfix barrage
###############################

system = TAL_sys
baseline_spec = {'spec_1': spec.baseline_yield,
                 'spec_2': spec.baseline_titer,
                 'spec_3': spec.baseline_productivity,}

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_specifications(**baseline_spec)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    system.simulate()
    
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                raise e
###############################

spec.load_spec_1 = spec.load_yield
# spec.load_spec_2 = spec.load_titer # defined in system
spec.load_spec_3 = spec.load_productivity

full_path = TAL_sys.path
fermenter_index = full_path.index(spec.titer_inhibitor_specification.reactor)
pre_fermenter_units_path = full_path[0:fermenter_index]

def model_specification():
    # TODO: bugfix barrage was removed for speed up, no failed evaluations found for now
    # try:
    for i in pre_fermenter_units_path: i._run()
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    model._system.simulate()   
    # except Exception as e:
    #     str_e = str(e).lower()
    #     print('Error in model spec: %s'%str_e)
    #     # raise e
    #     if 'sugar concentration' in str_e:
    #         # flowsheet('AcrylicAcid').F_mass /= 1000.
    #         raise e
    #     else:
    #         run_bugfix_barrage()
            
model.specification = model_specification

model.exception_hook = 'warn'

baseline_initial = model.metrics_at_baseline()
baseline = pd.DataFrame(data=np.array([[i for i in baseline_initial.values],]), 
                        columns=baseline_initial.keys())

model.evaluate(notify=50, autoload=True, autosave=20, file='unfinished_evaluation')
model.table.to_excel('all_results.xlsx')

# Baseline results
baseline_end = model.metrics_at_baseline()
dateTimeObj = datetime.now()
file_to_save = 'TAL_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)\
    + '_' + str(N_simulation) + 'sims'

baseline = baseline.append(baseline_end, ignore_index=True)
baseline.index = ('initial', 'end')
baseline.to_excel(file_to_save+'_0_baseline.xlsx')

# Parameters
parameters = model.get_parameters()
index_parameters = len(model.get_baseline_sample())
parameter_values = model.table.iloc[:, :index_parameters].copy()

#%%
# TEA results
for index_TEA, i in enumerate(models.metrics):
    if i.element == 'LCA': break
index_TEA = index_parameters + index_TEA
TEA_results = model.table.iloc[:, index_parameters:index_TEA].copy()
TEA_percentiles = TEA_results.quantile(q=percentiles)

# LCA_results
LCA_results = \
    model.table.iloc[:, index_TEA::].copy()
LCA_percentiles = LCA_results.quantile(q=percentiles)

# # Spearman's rank correlation

table = model.table

model.table = model.table.dropna()

spearman_results = model.spearman()
spearman_results.columns = pd.Index([i.name_with_units for i in model.metrics])

model.table = table

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


#%%
'''Output to Excel'''
with pd.ExcelWriter(file_to_save+'_1_full_evaluation.xlsx') as writer:
    parameter_values.to_excel(writer, sheet_name='Parameters')
    TEA_results.to_excel(writer, sheet_name='TEA results')
    TEA_percentiles.to_excel(writer, sheet_name='TEA percentiles')
    LCA_results.to_excel(writer, sheet_name='LCA results')
    LCA_percentiles.to_excel(writer, sheet_name='LCA percentiles')
    spearman_results.to_excel(writer, sheet_name='Spearman')
    # one_p_df.to_excel(writer, sheet_name='One-parameter')
    model.table.to_excel(writer, sheet_name='Raw data')


os.remove("unfinished_evaluation") 
# %%

# =============================================================================
# Temporary codes for debugging
# =============================================================================

# import numpy as np
# from TAL.analyses import models

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






#%% Evaluate MC across glycerol yield

# import numpy as np
# import os 
# from biosteam.plots import plot_montecarlo_across_coordinate


# # model = models.TAL_model
# # spec = models.spec

# prefix = os.path.dirname(__file__)
# yield_gly_setter = spec.load_yield_glycerol


# ygs = np.linspace(0.0001, 0.25, 20)
# spec.TRY_analysis = False
# model.evaluate_across_coordinate('Yield glycerol [% theo]', yield_gly_setter, ygs, 
#                                   xlfile = os.path.join(prefix, 'mc_across_glycerol_yield_2.xlsx'),
#                                   notify = 50)



# #%% Plot MC across glycerol yield
# import pandas as pd
# from biosteam.plots import plot_montecarlo_across_coordinate
# import os
# import numpy as np

# # ygs = np.linspace(0.0001, 0.25, 5)
# prefix = os.path.dirname(__file__)
# MC_across_gly_yield_df = pd.read_excel(os.path.join(prefix, 'mc_across_glycerol_yield_2.xlsx'),
#                                                 'Bior. Mini. sell. price')


# MC_across_gly_yield_df_array = np.array(MC_across_gly_yield_df)
# y_list = []

# for i in MC_across_gly_yield_df_array:
#     y_list.append(i[1:])

# MC_across_gly_yield = np.array(y_list)

# MC_across_gly_yield = MC_across_gly_yield[~np.isnan(MC_across_gly_yield).any(axis=1), :] # filter out nans
 
# mcac = plot_montecarlo_across_coordinate(ygs, MC_across_gly_yield)







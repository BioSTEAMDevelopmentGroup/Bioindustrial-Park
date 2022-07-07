#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Sarang Bhagwat <sarangb2@illinois.edu>
# Yoel Cortes-Pena <yoelcortes@gmail.com>, and Yalin Li <mailto.yalin.li@gmail.com> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon Apr 13 10:24:42 2020

Modified from the biorefineries constructed in [1] and [2] for the production of
BDO from lignocellulosic feedstocks

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted,
    2020.

@author: yalinli_cabbi
"""


# %% 

# =============================================================================
# Setup
# =============================================================================
from warnings import filterwarnings
filterwarnings('ignore')
import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.utils import TicToc
from biorefineries.BDO.system_MS3 import (
    spec, BDO_sys, get_MEK_MPSP, get_GWP, get_FEC, flowsheet, BDO_tea
)
from biorefineries.BDO.analyses import models
from datetime import datetime


R301 = flowsheet('R301')
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
model = models.BDO_model
# R301.set_titer_limit = True

# Set seed to make sure each time the same set of random numbers will be used
np.random.seed(3222) # 3221
N_simulation = 500 # 1000

samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)


###############################
# Bugfix barrage
###############################

system = BDO_sys
baseline_spec = {'spec_1': spec.spec_1,
                 'spec_2': spec.spec_2,
                 'spec_3': spec.spec_3,}

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
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

full_path = BDO_sys.path
evaporator_index = full_path.index(spec.titer_inhibitor_specification.evaporator)
pre_evaporator_units_path = full_path[0:evaporator_index]

def model_specification():
    # TODO: bugfix barrage was removed for speed up, no failed evaluations found for now
    # try:
    for i in pre_evaporator_units_path: i._run()
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
file_to_save = 'BDO_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)

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


# %%

# # =============================================================================
# # Evaluate the min/max of one parameter each time to ensure the parameter can
# # independently affect the system
# # =============================================================================

# p_values = [[], [], []]
# MPSPs = [[], [], []]
# GWPs = [[], [], []]
# FECs = [[], [], []]

# # import pdb
# for p in parameters:
#     # pdb.set_trace()
#     # [p_min], [p_max] = p.distribution.range().tolist()
#     p_dist = p.distribution
#     # import pdb
#     # pdb.set_trace()
#     [p_min], [p_max] = p_dist.lower.tolist(), p_dist.upper.tolist()
#     # [p_min], [p_max] = p_dist.range()[0], p_dist.range()[1]
#     p_baseline = p.baseline
#     p_value = (p_min, p_max, p_baseline)
#     p.system = BDO_sys
#     for i in range(len(p_value)):
#         p.setter(p_value[i])
#         p_values[i].append(p_value[i])
#         model._system._converge()
#         spec.titer_inhibitor_specification.run_units()
#         spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
#         MPSP = get_MEK_MPSP()
#         MPSPs[i].append(MPSP)
#         GWPs[i].append(get_GWP())
#         FECs[i].append(get_FEC())
#         run_number += 1

# MPSP_baseline = np.asarray(MPSPs[2])
# MPSP_min_diff = np.asarray(MPSPs[0]) - MPSP_baseline
# MPSP_max_diff = np.asarray(MPSPs[1]) - MPSP_baseline

# GWP_baseline = np.asarray(GWPs[2])
# GWP_min_diff = np.asarray(GWPs[0]) - GWP_baseline
# GWP_max_diff = np.asarray(GWPs[1]) - GWP_baseline

# FEC_baseline = np.asarray(FECs[2])
# FEC_min_diff = np.asarray(FECs[0]) - FEC_baseline
# FEC_max_diff = np.asarray(FECs[1]) - FEC_baseline

# one_p_df = pd.DataFrame({
#     ('Parameter', 'Name'): [i.name_with_units for i in parameters],
#     ('Parameter', 'Baseline'): p_values[2],
#     ('Parameter', 'Min'): p_values[0],
#     ('Parameter', 'Max'): p_values[1],
#     ('MPSP [$/kg]', 'MPSP baseline'): MPSP_baseline,
#     ('MPSP [$/kg]', 'MPSP min'): MPSPs[0],
#     ('MPSP [$/kg]', 'MPSP min diff'): MPSP_min_diff,
#     ('MPSP [$/kg]', 'MPSP max'): MPSPs[1],
#     ('MPSP [$/kg]', 'MPSP max diff'): MPSP_max_diff,
#     ('GWP [kg CO2-eq/kg]', 'GWP baseline'): GWP_baseline,
#     ('GWP [kg CO2-eq/kg]', 'GWP min'): GWPs[0],
#     ('GWP [kg CO2-eq/kg]', 'GWP min diff'): GWP_min_diff,
#     ('GWP [kg CO2-eq/kg]', 'GWP max'): GWPs[1],
#     ('GWP [kg CO2-eq/kg]', 'GWP max diff'): GWP_max_diff,
#     ('FEC [MJ/kg]', 'FEC baseline'): FEC_baseline,
#     ('FEC [MJ/kg]', 'FEC min'): FECs[0],
#     ('FEC [MJ/kg]', 'FEC min diff'): FEC_min_diff,
#     ('FEC [MJ/kg]', 'FEC max'): FECs[1],
#     ('FEC [MJ/kg]', 'FEC max diff'): FEC_max_diff,
#     })

# time = timer.elapsed_time / 60.
# print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %%

# '''Get a quick plot'''
# IRR_plot_indices = [metric.index for metric in model.metrics
#                     if 'IRR' in metric.index[0] and 'MPSP' in metric.index[1]]
# IRR_plot_data = IRR_results[IRR_plot_indices].copy()
# IRR_plot_data.columns = models.IRRs.copy()
# IRR_plot_y = IRR_plot_data.sort_index(axis=1)
# IRR_plot_y = IRR_plot_y.dropna()
# IRR_plot_x = models.IRRs.copy()
# IRR_plot_x.sort()
# plot_montecarlo_across_coordinate(IRR_plot_x, IRR_plot_y)

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


# %%

# =============================================================================
# Temporary codes for debugging
# =============================================================================

# import numpy as np
# from BDO.analyses import models

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
import numpy as np
import os 
from biosteam.plots import plot_montecarlo_across_coordinate


# model = models.BDO_model
# spec = models.spec

prefix = os.path.dirname(__file__)
yield_gly_setter = spec.load_yield_glycerol


ygs = np.linspace(0.0001, 0.25, 20)
spec.TRY_analysis = False
model.evaluate_across_coordinate('Yield glycerol [% theo]', yield_gly_setter, ygs, 
                                  xlfile = os.path.join(prefix, 'mc_across_glycerol_yield_2.xlsx'),
                                  notify = 50)

#%% Plot MC across glycerol yield
import pandas as pd
from biosteam.plots import plot_montecarlo_across_coordinate
import os
import numpy as np

# ygs = np.linspace(0.0001, 0.25, 5)
prefix = os.path.dirname(__file__)
MC_across_gly_yield_df = pd.read_excel(os.path.join(prefix, 'mc_across_glycerol_yield_2.xlsx'),
                                                'Bior. Mini. sell. price')


MC_across_gly_yield_df_array = np.array(MC_across_gly_yield_df)
y_list = []

for i in MC_across_gly_yield_df_array:
    y_list.append(i[1:])

MC_across_gly_yield = np.array(y_list)

MC_across_gly_yield = MC_across_gly_yield[~np.isnan(MC_across_gly_yield).any(axis=1), :] # filter out nans
 
mcac = plot_montecarlo_across_coordinate(ygs, MC_across_gly_yield)







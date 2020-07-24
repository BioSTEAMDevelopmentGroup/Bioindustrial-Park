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
Created on Wed Jul 22 19:45:17 2020

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
from biosteam.plots import plot_montecarlo_across_coordinate
from lactic import models

percentiles = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1]


# %% 

# =============================================================================
# Evaluating across internal rate of return
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

model = models.model_IRR

'''
Note:
    Not using `evaluate_across_coordinate` function as changing IRR
    does not affect the system, using IRR as the metrics for evaluation 
    will save considerable time.
'''

'''Evaluate'''
np.random.seed(3221)
N_simulation = 100 # 1000
samples = model.sample(N=N_simulation, rule='L')

model.load_samples(samples)
model.evaluate()

parameter_len = len(model.get_baseline_sample())
results = model.table.iloc[:, parameter_len::].copy()
percentiles = results.quantile(q=percentiles)

'''To get a quick plot'''
plot_indices = [metric.index for metric in model.metrics
                if 'Net present value' not in metric.index[1]]
plot_data = model.table[plot_indices]
plot_montecarlo_across_coordinate(models.IRRs, plot_data.values)

'''Output to Excel'''
with pd.ExcelWriter('2_IRR.xlsx') as writer:
    results.to_excel(writer, sheet_name='IRR')
    percentiles.to_excel(writer, sheet_name='IRR percentiles')
    model.table.to_excel(writer, sheet_name='Raw data')

run_number = N_simulation
time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')





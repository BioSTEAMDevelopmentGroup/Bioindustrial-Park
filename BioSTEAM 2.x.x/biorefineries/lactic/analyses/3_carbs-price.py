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
Created on Wed Jul 22 19:48:14 2020

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
from lactic.analyses import models


# %%

# =============================================================================
# Evaluate across feedstock price and carbohydrate content
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

model = models.model_carbs_price
set_carbs = models.set_carbs
prices = models.prices

'''Evaluate'''
np.random.seed(3221)
# This is not a Monte Carlo simulation, this evaluation uses the baseline parameters
# to see the impacts of feedstock carbohydrate content
# The parameter is a fake one to enable the evaluation
N_simulation = 1
samples_1d = model.sample(N=N_simulation, rule='L')
samples = samples_1d[:, np.newaxis]
model.load_samples(samples)

carb_contents = np.arange(0.4, 0.701, 0.01)

data = model.evaluate_across_coordinate(
    'Carbohydate content', set_carbs, carb_contents, notify=True)

results = pd.DataFrame({
    ('Parameter', 'Carbohydrate content [dw%]'): carb_contents})

for (i, j) in zip(data.keys(), data.values()):
    results[i] = j[0]

'''Organize TEA data for easy plotting'''
TEA_x = [f'{i:.3f}' for i in carb_contents]
TEA_x *= len(prices)
TEA_y = sum(([f'{i:.0f}']*len(carb_contents) for i in prices), [])

MPSPs = []
GWPs = []
freshwater = []
for i in range(results.columns.shape[0]):
    if 'MPSP' in results.columns[i][1]:
        MPSPs +=  results[results.columns[i]].to_list()
    if 'GWP' in results.columns[i][1]:
        GWPs +=  results[results.columns[i]].to_list()
    if 'Freshwater' in results.columns[i][1]:
        freshwater +=  results[results.columns[i]].to_list()

TEA_plot_data = pd.DataFrame({
    'Carbohydrate content [dw%]': TEA_x,
    'Price [$/dry-ton]': TEA_y,
    'MPSP [$/kg]': MPSPs
    })

LCA_plot_data = pd.DataFrame({
    'Carbohydrate content [dw%]': carb_contents,    
    'GWP [kg CO2-eq/kg lactic acid]': GWPs,
    'Freshwater consumption [kg H2O/kg lactic acid]': freshwater
    })

'''Output to Excel'''
with pd.ExcelWriter('3_carbs-price.xlsx') as writer:
    TEA_plot_data.to_excel(writer, sheet_name='TEA plotting')
    LCA_plot_data.to_excel(writer, sheet_name='LCA plotting')
    results.to_excel(writer, sheet_name='Raw data')

run_number = samples.shape[0] * len(carb_contents)
time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')













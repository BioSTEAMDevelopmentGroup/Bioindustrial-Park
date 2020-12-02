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
Created on Tue Sep  1 17:43:13 2020

Modified from the biorefineries constructed in [1] and [2] for the production of
lactic acid from lignocellulosic feedstocks

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

import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.utils import TicToc
from lactic import system_diluted
from lactic._utils import set_yield


# %% 

# =============================================================================
# Evaluate system at different lactic acid titer and yield (conversion),
# using either regular strain (need lime addition during fermentation to neutralize
# the produced lactic acid) or acid-resistant strain (no neutralization need)
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()
run_number = 0
limits = [{}, {}] # regular, acid-resistant
lactics = {'yield': [],
           'titer': []}
MPSPs = {0.89: [],
         0.18: [],
         1.92: []}
NPVs = {0.89: [],
        0.18: [],
        1.92: []}
GWPs = {0.89: [],
        0.18: [],
        1.92: []}
FECs = {0.89: [],
        0.18: [],
        1.92: []}

yield_range = np.arange(0.3, 1.01, 0.025) - 1e-6
# yield_range = np.arange(0.3, 1.01, 0.5) - 1e-6

R301 = system_diluted.R301
R302 = system_diluted.R302
R401 = system_diluted.R401
S402 = system_diluted.S402

lactic_acid = system_diluted.lactic_acid
lactic_sys = system_diluted.lactic_sys
lactic_tea = system_diluted.lactic_tea

def solve_TEA():
    lactic_acid.price = 0
    for i in range(3):
        MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
    return MPSP

def update_productivity(productivity):
    R301.productivity = productivity
    R302.productivity = productivity * R302.ferm_ratio
    for unit in (R301, R302):
        unit._design()
        unit._cost()

def simulate_log_results(return_limit=False):
    try:
        lactic_sys.simulate()
        limit = R301.effluent_titer
        lactics['yield'].append(R301.cofermentation_rxns.X[0])
        lactics['titer'].append(R301.effluent_titer)
        for productivity in (0.89, 0.18, 1.92):
            update_productivity(productivity)
            MPSPs[productivity].append(solve_TEA())
            NPVs[productivity].append(lactic_tea.NPV)
            GWPs[productivity].append(system_diluted.get_GWP())
            FECs[productivity].append(system_diluted.get_FEC())
    except:
        limit = np.nan
        lactic_sys.empty_recycles()
        lactic_sys.reset_cache()
        for m in (lactics, MPSPs, NPVs, GWPs, FECs):
            for n in m.keys():
                m[n].append(np.nan)
    global run_number
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')
    if return_limit:
        return limit

def save_data_clear():
    df = pd.DataFrame({
        ('Lactic acid', 'Yield [g/g]'): lactics['yield'],
        ('Lactic acid', 'Titer [g/L]'): lactics['titer'],
        ('0.89 [g/L/hr]', 'MPSP [$/kg]'): MPSPs[0.89],
        ('0.89 [g/L/hr]', 'NPV [$]'): NPVs[0.89],
        ('0.89 [g/L/hr]', 'GWP [kg CO2-eq/kg]'): GWPs[0.89],
        ('0.89 [g/L/hr]', 'FEC [MJ/kg]'): FECs[0.89],
        ('0.18 [g/L/hr]', 'MPSP [$/kg]'): MPSPs[0.18],
        ('0.18 [g/L/hr]', 'NPV [$]'): NPVs[0.18],
        ('0.18 [g/L/hr]', 'GWP [kg CO2-eq/kg]'): GWPs[0.18],
        ('0.18 [g/L/hr]', 'FEC [MJ/kg]'): FECs[0.18],
        ('1.92 [g/L/hr]', 'MPSP [$/kg]'): MPSPs[1.92],
        ('1.92 [g/L/hr]', 'NPV [$]'): NPVs[1.92],
        ('1.92 [g/L/hr]', 'GWP [kg CO2-eq/kg]'): GWPs[1.92],
        ('1.92 [g/L/hr]', 'FEC [MJ/kg]'): FECs[1.92]
        })
    for i in (lactics, MPSPs, NPVs, GWPs, FECs):
        for j in i.keys():
            i[j] = []
    return df


# %%

# =============================================================================
# Regular strain
# =============================================================================

bst.speed_up()
R301.neutralization = True
R301.allow_concentration = False
R401.bypass = False
S402.bypass = False

print('\n---------- Regular Strain ----------')
# First determine the maximum achievable titer at a given yield
R301.allow_dilution = False
yield_with_baseline = [0.76] + yield_range.tolist()
for i in yield_with_baseline:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[0][i] = simulate_log_results(return_limit=True)

regular_limit = save_data_clear()

# Dilute the saccharified stream to achieve lower titers
R301.allow_dilution = True
for i in yield_range:
    titer_range = np.arange(40, limits[0][i], 2.5)
    titer_range = titer_range.tolist() + [limits[0][i]]
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

regular_data = save_data_clear()

with pd.ExcelWriter('regular1.xlsx') as writer:
    regular_limit.to_excel(writer, sheet_name='Regular limit')
    regular_data.to_excel(writer, sheet_name='Regular data')


# %%

# =============================================================================
# Acid-resistant strain    
# =============================================================================

bst.speed_up()
R301.neutralization = False
R301.allow_concentration = False
R401.bypass = True
S402.bypass = True

print('\n---------- Acid-resistant Strain ----------')
# First determine the maximum achievable titer at a given yield
R301.allow_dilution = False
for i in yield_with_baseline:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[1][i] = simulate_log_results(return_limit=True)

resistant_limit = save_data_clear()

# Only simulate for achievable titers
R301.allow_dilution = True
for i in yield_range:
    titer_range = np.arange(40, limits[0][i], 2.5)
    titer_range = titer_range.tolist() + [limits[0][i]]
    for j in titer_range:        
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

resistant_data = save_data_clear()


# %%

'''Output to Excel'''
with pd.ExcelWriter('2-1_ferm_diluted.xlsx') as writer:
    regular_limit.to_excel(writer, sheet_name='Regular limit')
    regular_data.to_excel(writer, sheet_name='Regular data')
    resistant_limit.to_excel(writer, sheet_name='Acid-resistant limit')
    resistant_data.to_excel(writer, sheet_name='Acid-resistant data')

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')



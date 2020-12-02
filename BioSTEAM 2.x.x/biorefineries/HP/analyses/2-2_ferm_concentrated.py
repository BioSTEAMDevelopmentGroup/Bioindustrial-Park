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
from lactic import system_concentrated
from lactic._chemicals import sugars
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
sugar_concs = [[], []]
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

yield_range = np.arange(0.3, 1.01, 0.05) - 1e-6
# yield_range = np.arange(0.3, 1.01, 0.025) - 1e-6
# yield_range = np.arange(0.3, 1.01, 0.5) - 1e-6

# For the concentrated, batch scenario (regular and acid-resistant scenarios).
# Below the lower limit is the diluted scenario; between lower and mid limits
# is the concentrated, continuous reactor scenario; the upper limit is the
# simulation cap
limits = [[{}, {}], [{}, {}], [{}, {}]]

E301 = system_concentrated.E301
R301 = system_concentrated.R301
R302 = system_concentrated.R302
R401 = system_concentrated.R401
S402 = system_concentrated.S402

lactic_acid = system_concentrated.lactic_acid
lactic_sys = system_concentrated.lactic_sys
lactic_tea = system_concentrated.lactic_tea

def compute_sugar_conc(stream):
    return stream.imass[sugars].sum()/stream.F_vol

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
        lactic_yield = R301.cofermentation_rxns.X[0]        

        for productivity in (0.89, 0.18, 1.92):
            update_productivity(productivity)
            MPSP = solve_TEA()
            GWP = system_concentrated.get_GWP()
            FEC = system_concentrated.get_FEC()
            if lactic_yield < 0.3:
                MPSPs[productivity].append(MPSP)
                GWPs[productivity].append(GWP)
                FECs[productivity].append(FEC)
            elif lactic_yield==lactics['yield'][-1] and MPSP>MPSPs[productivity][-1]:
                MPSPs[productivity].append(np.nan)
                GWPs[productivity].append(np.nan)
                FECs[productivity].append(np.nan)
            else:
                MPSPs[productivity].append(MPSP)
                GWPs[productivity].append(GWP)
                FECs[productivity].append(FEC)
            NPVs[productivity].append(lactic_tea.NPV)
        lactics['yield'].append(lactic_yield)
        lactics['titer'].append(R301.effluent_titer)
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
R301.allow_dilution = False
R301.neutralization = True
R401.bypass = False
S402.bypass = False

R301.mode = 'Batch'
# This titer cannot be achieved, the system will find the maximum achievable titer
R301.titer_limit = 600

print('\n---------- Regular Strain Batch Limits ----------')
# Find the maximum achievable titer without concentration
E301.bypass = True
R301.allow_concentration = False
bst.speed_up()
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[0][0][i] = simulate_log_results(return_limit=True)

regular_limit1 = save_data_clear()

# Find the maximum achievable titer with concentration in batch mode
E301.bypass = False
R301.allow_concentration = True
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[1][0][i] = simulate_log_results(return_limit=True)

regular_limit2 = save_data_clear()

print('\n---------- Regular Strain Batch Mode ----------')
# Concentrate the saccharified stream to achieve higher titers in batch mode
for i in yield_range:
    titer_range = np.linspace(limits[0][0][i], limits[1][0][i], 5)
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

regular_batch = save_data_clear()
regular_batch.to_excel('regular_batch.xlsx')

print('\n---------- Regular Strain Continuous Limits ----------')
# Find the maximum achievable titer with concentration in continuous mode
R301.mode = 'Continuous'
R301.titer_limit = 600
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[2][0][i] = simulate_log_results(return_limit=True)

regular_limit3 = save_data_clear()
regular_limits = pd.concat((regular_limit1, regular_limit2, regular_limit3))

print('\n---------- Regular Strain Continuous Mode ----------')
# Concentrate the saccharified stream to achieve higher titers in continuous mode
for i in yield_range:
    if limits[2][0][i] is np.nan:
        limit = 205
    else:
        limit = min(limits[2][0][i], 205)
    titer_range = np.arange(limits[1][0][i], limit, 2.5)
    titer_range = titer_range.tolist() + [limit]
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

regular_continuous = save_data_clear()

with pd.ExcelWriter('regular2.xlsx') as writer:
    regular_limits.to_excel(writer, sheet_name='Regular limits')
    # regular_batch.to_excel(writer, sheet_name='Regular batch')
    regular_continuous.to_excel(writer, sheet_name='Regular continuous')


# %%

# =============================================================================
# Acid-resistant strain    
# =============================================================================

R301.neutralization = False
R401.bypass = True
S402.bypass = True

R301.mode = 'Batch'
R301.titer_limit = 600

print('\n---------- Acid-resistant Strain Batch Limits ----------')
# Find the maximum achievable titer without concentration
E301.bypass = True
R301.allow_concentration = False
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[0][1][i] = simulate_log_results(return_limit=True)

resistant_limit1 = save_data_clear()

# Find the maximum achievable titer with concentration in batch mode
E301.bypass = False
R301.allow_concentration = True
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[1][1][i] = simulate_log_results(return_limit=True)

resistant_limit2 = save_data_clear()

print('\n---------- Acid-resistant Strain Batch Mode ----------')
# Concentrate the saccharified stream to achieve higher titers in batch mode
for i in yield_range:
    titer_range = np.arange(limits[0][1][i], limits[1][1][i], 2.5)
    titer_range = titer_range.tolist() + [limit]
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

resistant_batch = save_data_clear()

print('\n---------- Acid-resistant Strain Continuous Limits ----------')
# Find the maximum achievable titer with concentration in continuous mode
R301.mode = 'Continuous'
R301.titer_limit = 600
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    limits[2][1][i] = simulate_log_results(return_limit=True)

resistant_limit3 = save_data_clear()
resistant_limits = pd.concat((resistant_limit1, resistant_limit2, resistant_limit3))

print('\n---------- Acid-resistant Strain Continuous Mode ----------')
# Concentrate the saccharified stream to achieve higher titers in continuous mode
for i in yield_range:
    if limits[2][1][i] is np.nan:
        limit = 205
    else:
        limit = min(limits[2][1][i], 205)
    titer_range = np.arange(limits[1][1][i], limit, 2.5)
    titer_range = titer_range.tolist() + [limit]
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        set_yield(i, R301, R302)
        simulate_log_results(return_limit=False)

resistant_continuous = save_data_clear()


# %%

'''Output to Excel'''
with pd.ExcelWriter('2-2_ferm_concentrated.xlsx') as writer:
    regular_limits.to_excel(writer, sheet_name='Regular limits')
    resistant_limits.to_excel(writer, sheet_name='Acid-resistant limits')
    regular_batch.to_excel(writer, sheet_name='Regular batch')
    resistant_batch.to_excel(writer, sheet_name='Acid-resistant batch')
    regular_continuous.to_excel(writer, sheet_name='Regular continuous')
    resistant_continuous.to_excel(writer, sheet_name='Acid-resistant continuous')

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')



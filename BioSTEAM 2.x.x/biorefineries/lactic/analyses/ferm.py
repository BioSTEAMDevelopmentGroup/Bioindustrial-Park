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


# %% 

# =============================================================================
# Setup
# =============================================================================

import numpy as np
import pandas as pd
import biosteam as bst
from biosteam.utils import TicToc
from biorefineries.lactic import (system_SSCF as SSCF, system_SHF as SHF)
from biorefineries.lactic._chemicals import sugars
from biorefineries.lactic._utils import set_yield

SSCF.simulate_and_print()
SHF.simulate_and_print()


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
lactics = {'target yield': [],
           'target titer': [],
           'actual yield': [],
           'actual titer': []}
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

def compute_sugar_conc(stream):
    return stream.imass[sugars].sum()/stream.F_vol

def solve_TEA(lactic_acid, lactic_tea):
    lactic_acid.price = 0
    for i in range(3):
        MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
    return MPSP

def update_productivity(R301, R302, productivity):
    R301.productivity = productivity
    R302.productivity = productivity * R302.ferm_ratio
    for unit in (R301, R302):
        unit._design()
        unit._cost()

def simulate_log_results(system, return_limit=False):
    R301 = system.R301
    R302 = system.R302
    lactic_acid = system.lactic_acid
    lactic_sys = system.lactic_sys
    lactic_tea = system.lactic_tea
    
    try:
        lactic_sys.simulate()
        # limit = R301.effluent_titer
    
        for productivity in (0.89, 0.18, 1.92):
            update_productivity(R301, R302, productivity)
            MPSP = solve_TEA(lactic_acid, lactic_tea)
            GWP = system.get_GWP()
            FEC = system.get_FEC()
            MPSPs[productivity].append(MPSP)
            GWPs[productivity].append(GWP)
            FECs[productivity].append(FEC)
            NPVs[productivity].append(lactic_tea.NPV)
        lactics['actual yield'].append(R301.lactic_yield)
        lactics['actual titer'].append(R301.effluent_titer)
    except:
        breakpoint()
    # except:
    #     limit = np.nan
    #     lactic_sys.empty_recycles()
    #     lactic_sys.reset_cache()
    #     for m in (lactics, MPSPs, NPVs, GWPs, FECs):
    #         for n in m.keys():
    #             m[n].append(np.nan)
    global run_number
    run_number += 1
    print(f'Target yield is {round(R301.target_yield,2)}, actual yield is {round(R301.cofermentation_rxns.X[0], 2)}; ' \
          f'target titer is {R301.target_titer}, actual titer is {round(R301.effluent_titer)}')
    # print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')
    # if return_limit:
    #     return limit

def save_data_clear():
    df = pd.DataFrame({
        ('Lactic acid', 'Target Yield [g/g]'): lactics['target yield'],
        ('Lactic acid', 'Target Titer [g/L]'): lactics['target titer'],
        ('Lactic acid', 'Actual Yield [g/g]'): lactics['actual yield'],
        ('Lactic acid', 'Actual Titer [g/L]'): lactics['actual titer'],
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

def run_TRY(yield_range, system, mode, feed_freq, if_resistant, titer_range):
    bst.speed_up()
    R301 = system.R301
    R302 = system.R302
    R401 = system.R401
    S402 = system.S402
    R301.mode = mode
    R301.feed_freq = feed_freq
    R301.allow_dilution = True
    if system is SHF:
        R301.allow_concentration = True
    R401.bypass = if_resistant
    S402.bypass = if_resistant
    for i in yield_range:
        for j in titer_range:
            lactics['target yield'].append(i)
            lactics['target titer'].append(j)
            R301.target_yield = i
            R301.target_titer = j
            set_yield(i, R301, R302)
            simulate_log_results(system, return_limit=False)
            if R301.effluent_titer+2<j:
                n_remained = len(titer_range) - titer_range.index(j) - 1
                lactics['target yield'].extend([i]*n_remained)
                lactics['target titer'].extend(titer_range[titer_range.index(j)+1:])
                lactics['actual yield'].extend([i]*n_remained)
                lactics['actual titer'].extend(['']*n_remained)
                for p in (0.89, 0.18, 1.92):
                    for m in (MPSPs, NPVs,GWPs, FECs):
                        m[p].extend(['']*n_remained)
                break


# %%

# =============================================================================
# Separate hydrolysis and fermentation, regular strain
# =============================================================================

# yield_range = np.arange(0.3, 1.01, 0.025) - 1e-6
yield_range = np.arange(0.3, 1.01, 0.5) - 1e-6
yield_range = yield_range.tolist()

# titer_range = np.arange(50, 250, 1)
titer_range = np.arange(50, 220, 50)
titer_range = titer_range.tolist()

print('\n---------- SSCF Regular Strain Batch Mode ----------')
run_TRY(yield_range=yield_range, system=SSCF, mode='batch', feed_freq=1,
        if_resistant=False, titer_range=titer_range)
SSCF_reg_b = save_data_clear()
SSCF_reg_b.to_excel('SSCF_reg_batch.xlsx')

for i in range(3):
    print(f'\n---------- SHF Regular Strain Batch Feed {i+1} Times ----------')
    run_TRY(yield_range=yield_range, system=SHF, mode='batch', feed_freq=i+1,
            if_resistant=False, titer_range=titer_range)
    SHF_reg_b = save_data_clear()
    SHF_reg_b.to_excel(f'SHF_reg_batch{i+1}.xlsx')

print('\n---------- SHF Regular Strain Continuous ----------')
run_TRY(yield_range=yield_range, system=SHF, mode='continuous', feed_freq=1,
        if_resistant=False, titer_range=titer_range)
SHF_reg_c = save_data_clear()
SHF_reg_c.to_excel('SHF_reg_continuous.xlsx')


# %%

# print('\n---------- Regular Strain Continuous Limits ----------')
# # Find the maximum achievable titer with concentration in continuous mode
# R301.mode = 'continuous'
# R301.target_titer = 600
# for i in yield_range:
#     R301.target_yield = i
#     set_yield(i, R301, R302)
#     limits[2][0][i] = simulate_log_results(return_limit=True)

# regular_limit3 = save_data_clear()
# regular_limits = pd.concat((regular_limit1, regular_limit2, regular_limit3))

# print('\n---------- Regular Strain Continuous Mode ----------')
# # Concentrate the saccharified stream to achieve higher titers in continuous mode
# for i in yield_range:
#     if limits[2][0][i] is np.nan:
#         limit = 205
#     else:
#         limit = min(limits[2][0][i], 205)
#     # titer_range = np.arange(limits[1][0][i], limit, 2.5)
#     titer_range = np.arange(limits[1][0][i], limit, 10)
#     titer_range = titer_range.tolist() + [limit]
#     for j in titer_range:
#         R301.target_yield = i
#         R301.target_titer = j
#         set_yield(i, R301, R302)
#         simulate_log_results(return_limit=False)

# regular_continuous = save_data_clear()

# # with pd.ExcelWriter('regular-2.xlsx') as writer:
# #     regular_limits.to_excel(writer, sheet_name='Regular limits')
# #     regular_batch.to_excel(writer, sheet_name='Regular batch')
# #     regular_continuous.to_excel(writer, sheet_name='Regular continuous')


# # %%

# # =============================================================================
# # Separate hydrolysis and fermentation, acid-resistant strain    
# # =============================================================================

# R301.neutralization = False
# R401.bypass = True
# S402.bypass = True

# R301.mode = 'batch'
# R301.target_titer = 600

# print('\n---------- Acid-resistant Strain Batch Limits ----------')
# # Find the maximum achievable titer without concentration
# R301.allow_concentration = False
# for i in yield_range:
#     R301.target_yield = i
#     set_yield(i, R301, R302)
#     limits[0][1][i] = simulate_log_results(return_limit=True)

# resistant_limit1 = save_data_clear()

# # Find the maximum achievable titer with concentration in batch mode
# R301.allow_concentration = True
# for i in yield_range:
#     R301.target_yield = i
#     set_yield(i, R301, R302)
#     limits[1][1][i] = simulate_log_results(return_limit=True)

# resistant_limit2 = save_data_clear()

# print('\n---------- Acid-resistant Strain Batch Mode ----------')
# # Concentrate the saccharified stream to achieve higher titers in batch mode
# for i in yield_range:
#     limit = limits[1][1][i]
#     titer_range = np.arange(limits[0][1][i], limits[1][1][i], 2.5)
#     titer_range = titer_range.tolist() + [limit]
#     for j in titer_range:
#         R301.target_yield = i
#         R301.target_titer = j
#         set_yield(i, R301, R302)
#         simulate_log_results(return_limit=False)

# resistant_batch = save_data_clear()

# print('\n---------- Acid-resistant Strain Continuous Limits ----------')
# # Find the maximum achievable titer with concentration in continuous mode
# R301.mode = 'continuous'
# R301.target_titer = 600
# for i in yield_range:
#     R301.target_yield = i
#     set_yield(i, R301, R302)
#     limits[2][1][i] = simulate_log_results(return_limit=True)

# resistant_limit3 = save_data_clear()
# resistant_limits = pd.concat((resistant_limit1, resistant_limit2, resistant_limit3))

# print('\n---------- Acid-resistant Strain Continuous Mode ----------')
# # Concentrate the saccharified stream to achieve higher titers in continuous mode
# for i in yield_range:
#     if limits[2][1][i] is np.nan:
#         limit = 205
#     else:
#         limit = min(limits[2][1][i], 205)
#     titer_range = np.arange(limits[1][1][i], limit, 2.5)
#     titer_range = titer_range.tolist() + [limit]
#     for j in titer_range:
#         R301.target_yield = i
#         R301.target_titer = j
#         set_yield(i, R301, R302)
#         simulate_log_results(return_limit=False)

# resistant_continuous = save_data_clear()


# # %%

# '''Output to Excel'''
# with pd.ExcelWriter('2-2_ferm_concentrated.xlsx') as writer:
#     regular_limits.to_excel(writer, sheet_name='Regular limits')
#     resistant_limits.to_excel(writer, sheet_name='Acid-resistant limits')
#     regular_batch.to_excel(writer, sheet_name='Regular batch')
#     resistant_batch.to_excel(writer, sheet_name='Acid-resistant batch')
#     regular_continuous.to_excel(writer, sheet_name='Regular continuous')
#     resistant_continuous.to_excel(writer, sheet_name='Acid-resistant continuous')

# time = timer.elapsed_time / 60
# print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-2021, Yalin Li <yalinli2@illinois.edu>,
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
from warnings import warn
from biosteam.utils import TicToc
from biorefineries.lactic import load_system, \
    SSCF_flowsheet, SSCF_funcs, SHF_flowsheet, SHF_funcs
from biorefineries.lactic._chemicals import sugars
from biorefineries.lactic._utils import set_yield

load_system('SSCF')
load_system('SHF')


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

def simulate_log_results(kind):
    kind = kind.upper()
    if 'SSCF' in kind:
        flowsheet = SSCF_flowsheet
        funcs = SSCF_funcs
    elif 'SHF' in kind:
        flowsheet = SHF_flowsheet
        funcs = SHF_funcs
    bst.main_flowsheet.set_flowsheet(flowsheet)
    load_system(kind)
    
    R301 = flowsheet.unit.R301
    R302 = flowsheet.unit.R302
    lactic_acid = flowsheet.stream.lactic_acid
    lactic_sys = flowsheet.system.lactic_sys
    lactic_tea = lactic_sys.TEA
    
    try:
        lactic_sys.simulate()
        for productivity in (0.89, 0.18, 1.92):
            update_productivity(R301, R302, productivity)
            MPSP = solve_TEA(lactic_acid, lactic_tea)
            GWP = funcs['get_GWP']()
            FEC = funcs['get_FEC']()
            MPSPs[productivity].append(MPSP)
            GWPs[productivity].append(GWP)
            FECs[productivity].append(FEC)
            NPVs[productivity].append(lactic_tea.NPV)
        lactics['actual yield'].append(R301.lactic_yield)
        lactics['actual titer'].append(R301.effluent_titer)
    except:
        warn(f'Simulation failed at target yield {round(R301.target_yield,2)}, ' \
              f'target titer {R301.target_titer}.')
        lactic_sys.empty_recycles()
        lactic_sys.reset_cache()
        lactics['actual yield'].append(np.nan)
        lactics['actual titer'].append(np.nan)
        for m in (MPSPs, NPVs, GWPs, FECs):
            for n in m.keys():
                m[n].append(np.nan)

    global run_number
    run_number += 1
    print(f'Target yield is {round(R301.target_yield,2)}, actual yield is {round(R301.cofermentation_rxns.X[0], 2)}; ' \
          f'target titer is {R301.target_titer}, actual titer is {round(R301.effluent_titer)}.')


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

def run_TRY(yield_range, kind, mode, feed_freq, if_resistant, titer_range):
    bst.speed_up()
    if 'SSCF' in str(kind).upper():
        flowsheet = SSCF_flowsheet
    elif 'SHF' in str(kind).upper():
        flowsheet = SHF_flowsheet
    bst.main_flowsheet.set_flowsheet(flowsheet)

    u = flowsheet.unit
    R301 = u.R301
    R302 = u.R302
    R401 = u.R401
    S402 = u.S402
    R301.mode = mode
    R301.feed_freq = feed_freq
    R301.allow_dilution = True
    if hasattr(R301, 'allow_concentration'):
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
            simulate_log_results(kind)
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
    
    print(f'\nSimulation time is {timer.elapsed_time/60:.1f} min for {run_number} runs.')


# %%

# =============================================================================
# Regular strain
# =============================================================================

# Change this to increase resolution
yield_range = np.arange(0.3, 1.01, 0.2) - 1e-6
# yield_range = np.arange(0.3, 1.01, 0.025) - 1e-6
yield_range = yield_range.tolist()

# Change this to increase resolution
titer_range = np.arange(50, 220, 50)
# titer_range = np.arange(50, 250, 1)
titer_range = titer_range.tolist()

print('\n---------- SSCF Regular Strain Batch Mode ----------')
run_TRY(yield_range=yield_range, kind='SSCF', mode='batch', feed_freq=1,
        if_resistant=False, titer_range=titer_range)
SSCF_reg_b = save_data_clear()
SSCF_reg_b.to_excel('SSCF_reg_batch.xlsx')

# Change this for fed-batch mode
# for i in (1, 3, 5, 10):
for i in (1,):
    print(f'\n---------- SHF Regular Strain Batch Feed {i+1} Times ----------')
    run_TRY(yield_range=yield_range, kind='SHF', mode='batch', feed_freq=i+1,
            if_resistant=False, titer_range=titer_range)
    SHF_reg_b = save_data_clear()
    SHF_reg_b.to_excel(f'SHF_reg_batch{i+1}.xlsx')

print('\n---------- SHF Regular Strain Continuous ----------')
run_TRY(yield_range=yield_range, kind='SHF', mode='continuous', feed_freq=1,
        if_resistant=False, titer_range=titer_range)
SHF_reg_c = save_data_clear()
SHF_reg_c.to_excel('SHF_reg_continuous.xlsx')


# %%

# =============================================================================
# Acid-resistant strain
# =============================================================================

print('\n---------- SSCF Regular Strain Batch Mode ----------')
run_TRY(yield_range=yield_range, kind='SSCF', mode='batch', feed_freq=1,
        if_resistant=True, titer_range=titer_range)
SSCF_reg_b = save_data_clear()
SSCF_reg_b.to_excel('SSCF_reg_batch.xlsx')

# Change this for fed-batch mode
# for i in (1, 3, 5, 10):
for i in (1,):
    print(f'\n---------- SHF Regular Strain Batch Feed {i} Times ----------')
    run_TRY(yield_range=yield_range, kind='SHF', mode='batch', feed_freq=i,
            if_resistant=True, titer_range=titer_range)
    SHF_reg_b = save_data_clear()
    SHF_reg_b.to_excel(f'SHF_reg_batch{i+1}.xlsx')

print('\n---------- SHF Regular Strain Continuous ----------')
run_TRY(yield_range=yield_range, kind='SHF', mode='continuous', feed_freq=1,
        if_resistant=True, titer_range=titer_range)
SHF_reg_c = save_data_clear()
SHF_reg_c.to_excel('SHF_reg_continuous.xlsx')



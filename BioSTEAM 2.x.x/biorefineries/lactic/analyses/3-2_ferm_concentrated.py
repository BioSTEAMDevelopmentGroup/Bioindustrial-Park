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
from lactic.chemicals import sugars
from lactic.utils import set_yield


# %% 

# =============================================================================
# Evaluate system at different lactic acid titer and yield (conversion),
# using either regular strain (need lime addition during fermentation to neutralize
# produced lactic acid) or acid-resistant strain (no neutralization need)
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

yield_range = np.arange(0.3, 1.01, 0.1) - 1e-6
# yield_range = np.arange(0.3, 1.01, 0.5) - 1e-6

# 0.76 is the baseline
yield_range = [0.76] + yield_range.tolist()
limits = [{}, {}]

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


# %%

# =============================================================================
# Regular strain
# =============================================================================

bst.speed_up()
R301.neutralization = True
R301.mode = 'Continuous'
R301.allow_dilution = False
R401.bypass = False
S402.bypass = False

run_number = 0
sugars_regular = [[], []]
lactic_regular = [[], []]
MPSPs_regular = [[], [], []]
NPVs_regular = [[], [], []]
LCA_regular = [[], []]

print('\n---------- Regular Strain ----------')
# First determine the maximum achievable titer without concentration
R301.allow_concentration = False
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    lactic_sys.simulate()
    limits[0][i] = R301.effluent_titer
    sugars_regular[0].append(compute_sugar_conc(E301.outs[0]))
    sugars_regular[1].append(compute_sugar_conc(R301.outs[0]))
    lactic_regular[0].append(R301.cofermentation_rxns.X[0])
    lactic_regular[1].append(R301.effluent_titer)
    LCA_regular[0].append(system_concentrated.get_GWP())
    LCA_regular[1].append(system_concentrated.get_FEC())
    
    update_productivity(0.89)
    MPSPs_regular[0].append(solve_TEA())
    NPVs_regular[0].append(lactic_tea.NPV)
    # Alternative productivities, productivity only affects reactor sizing
    # thus no need to simulate the system again
    update_productivity(0.18)
    MPSPs_regular[1].append(solve_TEA())
    NPVs_regular[1].append(lactic_tea.NPV) 
    update_productivity(1.92)
    MPSPs_regular[2].append(solve_TEA())
    NPVs_regular[2].append(lactic_tea.NPV)
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

# Concentrate the saccharified stream to achieve higher titers
R301.allow_concentration = True 
for i in yield_range:
    titer_range = np.arange(limits[0][i], 200, 10)
    for j in titer_range:
        R301.yield_limit = i
        R301.titer_limit = j
        # Baseline productivity
        set_yield(i, R301, R302)
        lactic_sys.simulate()
        sugars_regular[0].append(compute_sugar_conc(E301.outs[0]))
        sugars_regular[1].append(compute_sugar_conc(R301.outs[0]))
        lactic_regular[0].append(R301.cofermentation_rxns.X[0])
        lactic_regular[1].append(R301.effluent_titer)
        LCA_regular[0].append(system_concentrated.get_GWP())
        LCA_regular[1].append(system_concentrated.get_FEC())
        
        update_productivity(0.89)
        MPSPs_regular[0].append(solve_TEA())
        NPVs_regular[0].append(lactic_tea.NPV)
        # Alternative productivities, productivity only affects reactor sizing
        # thus no need to simulate the system again
        update_productivity(0.18)
        MPSPs_regular[1].append(solve_TEA())
        NPVs_regular[1].append(lactic_tea.NPV) 
        update_productivity(1.92)
        MPSPs_regular[2].append(solve_TEA())
        NPVs_regular[2].append(lactic_tea.NPV)
        run_number += 1
        print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

regular_data = pd.DataFrame({
    ('Concentrated', 'Sugar [g/L]'): sugars_regular[0],
    ('Fermenter effluent', 'Sugar [g/L]'): sugars_regular[1],
    ('Lactic', 'Yield [g/g]'): lactic_regular[0],
    ('Lactic', 'Titer [g/L]'): lactic_regular[1],
    ('Productivity=0.89 [g/L/hr] (baseline)', 'MPSP [$/kg]'):
        MPSPs_regular[0],
    ('Productivity=0.89 [g/L/hr] (baseline)', 'NPV [$]'):
        NPVs_regular[0],
    ('Productivity=0.18 [g/L/hr] (min)', 'MPSP [$/kg]'):
        MPSPs_regular[1],
    ('Productivity=0.18 [g/L/hr] (min)', 'NPV [$]'):
        NPVs_regular[1],
    ('Productivity=1.92 [g/L/hr] (max)', 'MPSP [$/kg]'):
        MPSPs_regular[2],
    ('Productivity=1.92 [g/L/hr] (max)', 'NPV [$]'):
        NPVs_regular[2],
    ('LCA', 'GWP [kg CO2-eq/kg]'): LCA_regular[0],
    ('LCA', 'FEC [MJ/kg]'): LCA_regular[1]
    })


# %%

# =============================================================================
# Acid-resistant strain    
# =============================================================================

bst.speed_up()
R301.neutralization = False
R301.allow_dilution = False
R401.bypass = True
S402.bypass = True

sugars_resistant = [[], []]
lactic_resistant = [[], []]
MPSPs_resistant = [[], [], []]
NPVs_resistant = [[], [], []]
LCA_resistant = [[], []]

print('\n---------- Acid-resistant Strain ----------')
# First determine the maximum achievable titer without concentration
R301.allow_concentration = False
for i in yield_range:
    R301.yield_limit = i
    set_yield(i, R301, R302)
    lactic_sys.simulate()
    limits[1][i] = R301.effluent_titer
    sugars_resistant[0].append(compute_sugar_conc(E301.outs[0]))
    sugars_resistant[1].append(compute_sugar_conc(R301.outs[0]))
    lactic_resistant[0].append(R301.cofermentation_rxns.X[0])
    lactic_resistant[1].append(R301.effluent_titer)
    LCA_resistant[0].append(system_concentrated.get_GWP())
    LCA_resistant[1].append(system_concentrated.get_FEC())
    
    update_productivity(0.89)
    MPSPs_resistant[0].append(solve_TEA())
    NPVs_resistant[0].append(lactic_tea.NPV)
    # Alternative productivities, productivity only affects reactor sizing
    # thus no need to simulate the system again
    update_productivity(0.18)
    MPSPs_resistant[1].append(solve_TEA())
    NPVs_resistant[1].append(lactic_tea.NPV) 
    update_productivity(1.92)
    MPSPs_resistant[2].append(solve_TEA())
    NPVs_resistant[2].append(lactic_tea.NPV)
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

# Concentrate the saccharified stream to achieve higher titers
R301.allow_concentration = True
for i in yield_range:
    titer_range = np.arange(limits[1][i], 200, 10)
    for j in titer_range:        
        R301.yield_limit = i
        R301.titer_limit = j
        # Baseline productivity
        set_yield(i, R301, R302)
        lactic_sys.simulate()
        sugars_resistant[0].append(compute_sugar_conc(E301.outs[0]))
        sugars_resistant[1].append(compute_sugar_conc(R301.outs[0]))
        lactic_resistant[0].append(R301.cofermentation_rxns.X[0])
        lactic_resistant[1].append(R301.effluent_titer)
        LCA_resistant[0].append(system_concentrated.get_GWP())
        LCA_resistant[1].append(system_concentrated.get_FEC())
        
        update_productivity(0.89)
        MPSPs_resistant[0].append(solve_TEA())
        NPVs_resistant[0].append(lactic_tea.NPV)
        # Alternative productivities, productivity only affects reactor sizing
        # thus no need to simulate the system again
        update_productivity(0.18)
        MPSPs_resistant[1].append(solve_TEA())
        NPVs_resistant[1].append(lactic_tea.NPV) 
        update_productivity(1.92)
        MPSPs_resistant[2].append(solve_TEA())
        NPVs_resistant[2].append(lactic_tea.NPV)
        run_number += 1
        print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

acid_resistant_data = pd.DataFrame({
    ('Concentrated', 'Sugar [g/L]'): sugars_resistant[0],
    ('Fermenter effluent', 'Sugar [g/L]'): sugars_resistant[1],
    ('Lactic', 'Yield [g/g]'): lactic_resistant[0],
    ('Lactic', 'Titer [g/L]'): lactic_resistant[1],
    ('Productivity=0.89 [g/L/hr] (baseline)', 'MPSP [$/kg]'):
        MPSPs_resistant[0],
    ('Productivity=0.89 [g/L/hr] (baseline)', 'NPV [$]'):
        NPVs_resistant[0],
    ('Productivity=0.18 [g/L/hr] (min)', 'MPSP [$/kg]'):
        MPSPs_resistant[1],
    ('Productivity=0.18 [g/L/hr] (min)', 'NPV [$]'):
        NPVs_resistant[1],
    ('Productivity=1.92 [g/L/hr] (max)', 'MPSP [$/kg]'):
        MPSPs_resistant[2],
    ('Productivity=1.92 [g/L/hr] (max)', 'NPV [$]'):
        NPVs_resistant[2],
    ('LCA', 'GWP [kg CO2-eq/kg]'): LCA_resistant[0],
    ('LCA', 'FEC [MJ/kg]'): LCA_resistant[1]
    })

    
# %%

'''Output to Excel'''
with pd.ExcelWriter(f'3-2_ferm_concentrated_{R301.mode}.xlsx') as writer:
    regular_data.to_excel(writer, sheet_name='Regular')
    acid_resistant_data.to_excel(writer, sheet_name='Acid-resistant')

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')



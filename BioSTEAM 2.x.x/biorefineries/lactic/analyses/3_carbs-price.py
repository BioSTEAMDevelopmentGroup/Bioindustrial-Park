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

_kg_per_ton = 907.18474
_feedstock_factor = _kg_per_ton / 0.8


# %%

# =============================================================================
# Model to evalute system across feedstock carbohydate content at different
# feedstock price metrics
# =============================================================================

def set_carbs(carbs_content, feedstock):
    carbs = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')
    dry_mass = feedstock.F_mass.copy() - feedstock.imass[('H2O',)].copy()
    old_carbs_mass_total = feedstock.imass[carbs].sum().copy()
    ratio = feedstock.get_normalized_mass(carbs)
    new_carbs_mass = dry_mass * carbs_content * ratio
    feedstock.set_flow(new_carbs_mass, 'kg/hr', carbs)
    mass_diff = new_carbs_mass.sum() - old_carbs_mass_total
    feedstock.imass['Extractives'] -= mass_diff
    if any(feedstock.mass < 0):
        raise ValueError(f'Carbohydrate content of {carbs_content*100:.1f} dw% is infeasible')

prices = np.arange(0, 210, 10)


# %%

# =============================================================================
# Evaluate across feedstock price and carbohydrate content
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()
run_number = 0

TEA_carbs = []
TEA_prices = []
titers = []
MPSPs = []
NPVs = []
GWPs = []
FECs = []

# Configuration 1
from biorefineries.lactic import system_concentrated as concentrated
carb_contents1 = np.arange(0.25, 0.59, 0.01)
carb_contents1 = carb_contents1.tolist() + [0.589]
concentrated.R301.allow_dilution = False
concentrated.R301.allow_concentration = True
concentrated.R301.mode = 'Batch'
concentrated.R301.titer_limit = 97.5
concentrated.E301.bypass = False

# Using two loops are not optimal, can potentially use Model and Metric to speed up
bst.speed_up()
for i in carb_contents1:
    set_carbs(i, concentrated.feedstock)
    concentrated.lactic_sys.simulate()
    titers.append(concentrated.R301.effluent_titer)
    GWPs.append(concentrated.get_GWP())
    FECs.append(concentrated.get_FEC())
    for j in prices:
        TEA_carbs.append(i)
        TEA_prices.append(j)
        concentrated.feedstock.price = j / _feedstock_factor
        concentrated.lactic_acid.price = 0
        for m in range(3):
            MPSP = concentrated.lactic_acid.price = \
                concentrated.lactic_tea.solve_price(concentrated.lactic_acid)
        MPSPs.append(MPSP)
        NPVs.append(concentrated.lactic_tea.NPV)
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

# Then concentration needed to get to the baseline titer
from biorefineries.lactic import system_diluted as diluted
carb_contents2 = np.arange(0.59, 0.701, 0.01).tolist()
diluted.R301.allow_dilution = True
diluted.R301.allow_concentration = False
diluted.R301.titer_limit = 97.5

bst.speed_up()
for i in carb_contents2:
    set_carbs(i, diluted.feedstock)
    diluted.lactic_sys.simulate()
    titers.append(diluted.R301.effluent_titer)
    GWPs.append(diluted.get_GWP())
    FECs.append(diluted.get_FEC())
    for j in prices:
        TEA_carbs.append(i)
        TEA_prices.append(j)
        diluted.feedstock.price = j / _feedstock_factor
        diluted.lactic_acid.price = 0
        for m in range(3):
            MPSP = diluted.lactic_acid.price = \
                diluted.lactic_tea.solve_price(diluted.lactic_acid)
        MPSPs.append(MPSP)
        NPVs.append(diluted.lactic_tea.NPV)
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

TEA_plot_data = pd.DataFrame({
    'Carbohydrate content [dw%]': TEA_carbs,
    'Price [$/dry-ton]': TEA_prices,
    'MPSP [$/kg]': MPSPs,
    'NPVs [$]': NPVs
    })

LCA_plot_data = pd.DataFrame({
    'Carbohydrate content [dw%]': carb_contents1+carb_contents2,
    'Titers [g/L]': titers,
    'GWP [kg CO2-eq/kg]': GWPs,
    'FEC [MJ/kg]': FECs
    })

'''Output to Excel'''
with pd.ExcelWriter('3_carbs-price.xlsx') as writer:
    TEA_plot_data.to_excel(writer, sheet_name='TEA')
    LCA_plot_data.to_excel(writer, sheet_name='LCA')

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')













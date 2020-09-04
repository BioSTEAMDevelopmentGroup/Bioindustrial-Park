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
Created on Wed Jul 22 20:02:49 2020

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
from lactic import system

U101 = system.U101


# %% 

# =============================================================================
# Diverge a portion of the feedstock directly to CHP to reduce reliance on
# natural gas
# =============================================================================

# Initiate a timer
timer = TicToc('timer')
timer.tic()

ratios = np.arange(0, 1, 0.05)

bst.speed_up()
run_number = 0
actuals_regular = [[], []]
MPSPs = []
NPVs = []
GWPs = []
FECs = []

for i in ratios:
    system.U101.diverge_to_CHP = i
    MPSPs.append(system.simulate_get_MPSP())
    NPVs.append(system.lactic_tea.NPV)
    GWPs.append(system.get_GWP())
    FECs.append(system.get_FEC())
    run_number += 1
    print(f'Run #{run_number}: {timer.elapsed_time:.0f} sec')

data = pd.DataFrame({
    'Diverge ratio': ratios,
    'MPSP [$/kg]': MPSPs,
    'NPV [$]': NPVs,
    'GWP [kg CO2-eq/kg]': GWPs,
    'FEC [MJ/kg]': FECs
    })


# %%

'''Get quick a quick plot'''
import matplotlib.pyplot as plt
from biosteam.utils import colors

fig, axs = plt.subplots(3, 1, figsize=(4, 5), sharex=True)
fig.subplots_adjust(hspace=0)

axs[0].plot(ratios, MPSPs, '-', color=colors.brown_shade.RGBn, linewidth=1.5)
axs[0].set_ylabel('MPSP [$/kg]')
axs[0].set_yticks(np.arange(1, 2.2, 0.2))
axs[0].set_ylim(1, 2)

axs[1].plot(ratios, GWPs, '-', color=colors.green_shade.RGBn, linewidth=1.5)
axs[1].set_ylabel('GWP [kg CO2-eq/kg]')
axs[1].set_yticks(np.arange(-20, 10, 5))
axs[1].set_ylim(-20, 10)

axs[2].plot(ratios, FECs, '-', color=colors.red_tint.RGBn, linewidth=1.5)
axs[2].set_ylabel('FEC [MJ/kg]')
axs[2].set_yticks(np.arange(-200, 100, 50))
axs[2].set_ylim(-200, 100)

axs[2].set_xlabel('Diverge ratio')
axs[2].set_xticks(np.arange(0, 1.2, 0.2))
axs[2].set_xlim(0, 1)

plt.show()

'''Output to Excel'''
with pd.ExcelWriter('2_diversion_to_CHP.xlsx') as writer:
    data.to_excel(writer, sheet_name='Results')

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')



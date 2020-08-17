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
Created on Mon Jul 27 14:58:24 2020

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
from lactic import system_succinic

system_succinic.R301.set_titer_limit = True


# %%

# =============================================================================
# Evaluate across feedstock succinic acid content, the baseline separation process 
# is not sufficient to get 88% lactic acid for feedstock with high succinic acid
# contents, therefore needs to add a supplementary distillation column D405
# =============================================================================

# Function to set feedstock succinic acid content
def set_succinic(feedstock, content):
    dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    feedstock.imass['SuccinicAcid'] = content * dry_mass
    # Use Extractives to close mass balance
    feedstock.imass['Extractives'] -= (feedstock.F_mass-feedstock.imass['H2O']) - dry_mass
    if any(feedstock.mass<0):
        raise ValueError(f'Succinic acid content of {content*100:.1f} dw% is infeasible')

# Initiate a timer
timer = TicToc('timer')
timer.tic()

succinic_contents = np.arange(0.07, 0.15, 0.01)
# succinic_contents = np.arange(0.14, 0.06, -0.01)
succinic_contents = [0] + succinic_contents.tolist() + [0.145]
run_number = 0


# %%

# =============================================================================
# With supplementary purification
# =============================================================================

print('\n-------- With D405 --------')
D405 = system_succinic.D405
D405_Lrs = []
MPSPs = []
purities = []
NPVs = []
GWPs = []
freshwater = []

bst.speed_up()

def simulate_log_results():
    system_succinic.lactic_acid.price = 0
    MPSP = system_succinic.simulate_get_MPSP()
    MPSPs.append(MPSP)
    purity = system_succinic.lactic_acid.get_mass_composition('LacticAcid')
    purities.append(purity)
    NPVs.append(system_succinic.lactic_tea.NPV)
    GWPs.append(system_succinic.get_functional_GWP())
    freshwater.append(system_succinic.get_functional_H2O())
    D405_Lrs.append(D405.Lr)
    print(f'{i:.1%} succinic acid:')
    print(f'D405 Lr: {D405.Lr:.3f}, MPSP: ${MPSP:.3f}/kg, purity: {purity:.1%}\n')

for i in succinic_contents:
    j = round(i, 3)
    if j == 0:
        system_succinic.F402.V = 0.9
    elif j < 0.09:
        system_succinic.F402.V = 0.92
    elif j < 0.13:
        system_succinic.F402.V = 0.94
    elif j == 0.13:
        system_succinic.F402.V = 0.95
        system_succinic.D405.specification = None
        system_succinic.D405.Lr = 0.8
    set_succinic(system_succinic.feedstock, j)
    simulate_log_results()
    run_number += 1

w_D405_data = pd.DataFrame({
    'Succinic acid content [%]': succinic_contents,
    'D405 Lr': D405_Lrs,
    'Lactic acid purity [%]': purities,
    'MPSP [$/kg]': MPSPs,
    'NPV [$]': NPVs,
    'GWP [kg CO2-eq/kg]': GWPs,
    'Freshwater [kg H2O/kg]': freshwater
    })

time = timer.elapsed_time / 60
print(f'\nSimulation time for {run_number} runs is: {time:.1f} min')


# %%

'''Output to Excel'''
with pd.ExcelWriter('4-2_succinic_w_D405.xlsx') as writer:
    w_D405_data.to_excel(writer, sheet_name='With D405')







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
from lactic import system

system.R301.set_titer_limit = True


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

succinic_contents = np.arange(0, 0.11, 0.01)
# succinic_contents = np.arange(0.07, -0.01, -0.01)
# succinic_contents = succinic_contents.tolist() + [0.08, 0.09, 0.1]
run_number = 0


# %%

# =============================================================================
# Without supplementary purification
# =============================================================================

print('\n-------- Without D405 --------')
succinic_axis = []
F402 = system.F402
F402_Vs = []
MPSPs = []
purities = []
NPVs = []
GWPs = []
freshwater = []

bst.speed_up()

def simulate_log_results():
    system.lactic_acid.price = 0
    MPSP = system.simulate_get_MPSP()
    MPSPs.append(MPSP)
    purity = system.lactic_acid.get_mass_composition('LacticAcid')
    purities.append(purity)
    NPVs.append(system.lactic_tea.NPV)
    GWPs.append(system.get_functional_GWP())
    freshwater.append(system.get_functional_H2O())
    F402_Vs.append(F402.V)
    print(f'{i:.1%} succinic acid:')
    print(f'F402 V: {F402.V:.3f}, MPSP: ${MPSP:.3f}/kg, purity: {purity:.1%}\n')

for i in succinic_contents:
    j = round(i, 3)
    set_succinic(system.feedstock, j)
    if j < 0.08:
        simulate_log_results()
        succinic_axis.append(i)
        run_number += 1
    elif j == 0.08:
        F402.specification = None
        F402.V = 0.941
        simulate_log_results()
        succinic_axis.append(i)
        run_number += 1
    elif j == 0.09:
        for V in (0.97, 0.98, 0.99):
            F402.V = V
            simulate_log_results()
            succinic_axis.append(i)
            run_number += 1
    elif j == 0.1:
        for V in (0.98, 0.99):
            F402.V = V
            simulate_log_results()
            succinic_axis.append(i)
            run_number += 1
    
wo_D405_data = pd.DataFrame({
    'Succinic acid content [%]': succinic_axis,
    'F402 V': F402_Vs,
    'Lactic acid purity [%]': purities,
    'MPSP [$/kg]': MPSPs,
    'Net present value [$]': NPVs,
    'GWP [kg CO2-eq/kg]': GWPs,
    'Freshwater [kg H2O/kg]': freshwater
    })


# %%

'''Output to Excel'''
with pd.ExcelWriter('4-1_succinic_wo_D405.xlsx') as writer:
    wo_D405_data.to_excel(writer, sheet_name='Without D405')







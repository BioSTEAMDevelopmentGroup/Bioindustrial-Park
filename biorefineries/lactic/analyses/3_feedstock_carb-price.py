#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Setup
# =============================================================================

import numpy as np, pandas as pd, biosteam as bst
from biosteam.utils import Timer
from biorefineries.lactic import load, create_funcs
from biorefineries.lactic.utils import feedstock_factor


# %%

# =============================================================================
# Model to evaluate system across feedstock carbohydrate content at different
# feedstock price metrics
# =============================================================================

carbs = ('Glucan', 'Xylan', 'Arabinan', 'Galactan', 'Mannan')

def set_carbs(carbs_content, feedstock):
    dry_mass = feedstock.F_mass.copy() - feedstock.imass[('H2O',)].copy()
    old_carbs_mass_total = feedstock.imass[carbs].sum().copy()
    ratio = feedstock.get_normalized_mass(carbs)
    new_carbs_mass = dry_mass * carbs_content * ratio
    feedstock.set_flow(new_carbs_mass, 'kg/hr', carbs)
    mass_diff = new_carbs_mass.sum() - old_carbs_mass_total
    feedstock.imass['Extractives'] -= mass_diff
    if any(feedstock.mass < 0):
        raise ValueError(f'Carbohydrate content of {carbs_content*100:.1f} dw% is infeasible')

get_feedstock_carbs = \
    lambda: feedstock.imass[carbs].sum()/(feedstock.F_mass-feedstock.imass['H2O'])

prices = np.arange(0, 210, 10)


# %%

# =============================================================================
# Evaluate across feedstock price and carbohydrate content
# =============================================================================

# Initiate a timer
timer = Timer('timer')
timer.start()

TEA_set_carbs = []
TEA_actual_carbs = []
TEA_prices = []
TEA_titers = []
TEA_yields = []
LCA_actual_carbs = []
LCA_titers = []
LCA_yields = []
MPSPs = []
NPVs = []
GWPs = []
FECs = []

# Configuration 2
load('SHF')
carb_contents1 = np.arange(0.25, 0.59, 0.01)
carb_contents1 = carb_contents1.tolist() + [0.589]
SHF_flowsheet = bst.main_flowsheet
R301 = SHF_flowsheet.unit.R301
R301.allow_dilution = False
R301.allow_concentration = True
R301.mode = 'batch'
R301.feed_freq = 1
R301.target_titer = 97.5

feedstock = SHF_flowsheet.stream.feedstock
lactic_acid = SHF_flowsheet.stream.lactic_acid
lactic_sys = SHF_flowsheet.system.lactic_sys
lactic_tea = lactic_sys.TEA

# Using two loops are not optimal, can potentially use Model and Metric to speed up
for i in carb_contents1:
    set_carbs(i, feedstock)
    LCA_actual_carbs.append(get_feedstock_carbs())
    lactic_sys.simulate()
    LCA_titers.append(R301.effluent_titer)
    LCA_yields.append(R301.lactic_yield)
    SHF_funcs = create_funcs(flowsheet=SHF_flowsheet)
    GWPs.append(SHF_funcs['get_GWP']())
    FECs.append(SHF_funcs['get_FEC']())
    for j in prices:
        TEA_set_carbs.append(i)
        TEA_actual_carbs.append(get_feedstock_carbs())
        TEA_prices.append(j)
        TEA_titers.append(R301.effluent_titer)
        TEA_yields.append(R301.lactic_yield)
        feedstock.price = j / feedstock_factor
        lactic_acid.price = 0
        for m in range(3):
            MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
        MPSPs.append(MPSP)
        NPVs.append(lactic_tea.NPV)

# Then concentration needed to get to the baseline titer
load('SSCF')
carb_contents2 = np.arange(0.59, 0.701, 0.01).tolist()
SSCF_flowsheet = bst.main_flowsheet
R301 = SSCF_flowsheet.unit.R301
R301.allow_dilution = True
R301.allow_concentration = False
R301.mode = 'batch'
R301.feed_freq = 1
R301.target_titer = 97.5

feedstock = SSCF_flowsheet.stream.feedstock
lactic_acid = SSCF_flowsheet.stream.lactic_acid
lactic_sys = SSCF_flowsheet.system.lactic_sys
lactic_tea = lactic_sys.TEA

for i in carb_contents2:
    set_carbs(i, feedstock)
    LCA_actual_carbs.append(get_feedstock_carbs())
    lactic_sys.simulate()
    LCA_titers.append(R301.effluent_titer)
    LCA_yields.append(R301.lactic_yield)
    SSCF_funcs = create_funcs(flowsheet=SSCF_flowsheet)
    GWPs.append(SSCF_funcs['get_GWP']())
    FECs.append(SSCF_funcs['get_FEC']())
    for j in prices:
        TEA_set_carbs.append(i)
        TEA_actual_carbs.append(get_feedstock_carbs())
        TEA_prices.append(j)
        TEA_titers.append(R301.effluent_titer)
        TEA_yields.append(R301.lactic_yield)
        feedstock.price = j / feedstock_factor
        lactic_acid.price = 0
        for m in range(3):
            MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
        MPSPs.append(MPSP)
        NPVs.append(lactic_tea.NPV)

TEA_data = pd.DataFrame({
    'Target carb content [dw%]': TEA_set_carbs,
    'Actual carb content [dw%]': TEA_actual_carbs,
    'Price [$/dry-ton]': TEA_prices,
    'Titer [g/L]': TEA_titers,
    'Yield [dw%]': TEA_yields,
    'MPSP [$/kg]': MPSPs,
    'NPVs [$]': NPVs
    })

LCA_data = pd.DataFrame({
    'Target carb content [dw%]': carb_contents1+carb_contents2,
    'Actual carb content [dw%]': LCA_actual_carbs,
    'Titer [g/L]': LCA_titers,
    'Yield [dw%]': LCA_yields,
    'GWP [kg CO2-eq/kg]': GWPs,
    'FEC [MJ/kg]': FECs
    })

'''Output to Excel'''
with pd.ExcelWriter('3_feedstock_carb-price.xlsx') as writer:
    TEA_data.to_excel(writer, sheet_name='TEA')
    LCA_data.to_excel(writer, sheet_name='LCA')

print(f'\nSimulation time is {timer.elapsed_time/60:.1f} min')
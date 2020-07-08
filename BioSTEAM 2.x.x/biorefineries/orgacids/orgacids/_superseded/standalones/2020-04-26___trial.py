#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.12.6
thermosteam v0.12.4
'''

import biosteam as bst
from orgacids.system import *

find = bst.main_flowsheet
broth_for_separation = find.stream.broth_for_separation

# # Right now minimum selling price changes a bit each time the system simulates 
# # due to separation_sys starting F401, the outs change each simulation
# print(F401.outs[1].F_mass)

# for i in range(3): orgacids_sys.simulate()
# for i in range(3): 
#     lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_BT_tea)
    
print(lactic_acid.price)

# orgacids_sys.save_report('orgacids_sys.xlsx')


# %% Sugar concentration

ferm_in = H301.outs[0]

total_sugar = 0
for sugar in ['Glucan', 'Glucose', 'Xylan', 'Xylose', 'Arabinose', 'Galactan']:
    total_sugar += ferm_in.imass[sugar]

total_sugar_conc = total_sugar/ferm_in.F_vol    
print(total_sugar_conc)


# %% Feedstock handling and pretreatment systems
    
from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()


# %% Fermentation system

# fermentation_sys_units = (H301, M301, M302, R301, R302, T301, M303)

# for unit in fermentation_sys_units:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

fermentation_sys.simulate()


# %% Separation system

# update_stripping_water()

# separation_sys_units_1 = (U401, M401, S401, R401)

# for unit in separation_sys_units_1:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

# update_separation_sulfuric_acid()

# separation_sys_units_2 = (T401, M402, S402, 
#                           F401, H401,
#                           R402, S403,
#                           H402, R403, S404, 
#                           H403, S405, H404)

# for unit in separation_sys_units_2:
#     try:
#         unit.simulate()
#     except: 
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

separation_sys.simulate()


# %% Wastewtaer treatment

# wastewater_sys_units_1 = (M501, WWT_cost, R501)

# for unit in wastewater_sys_units_1:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break
                          
# aerobic_digestion_sys_units = (M502, R502, S501, S502, M503, M504, S503)

# for unit in aerobic_digestion_sys_units:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

# aerobic_digestion_sys.simulate()

# wastewater_sys_units_2 = (S504, M505)

# for unit in wastewater_sys_units_2:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

wastewater_sys.simulate()


# %% Facilities

# facilities_sys_units = (CIP, ADP, FWT, BT,
#                         S601, T601, T602, T603, S602, T604, T605,
#                         P601, T606,
#                         CT, PWC)

# for unit in facilities_sys_units:
#     try:
#         unit.simulate()
#     except:
#         print(unit)
#         unit.show(N=100)
#         unit.simulate()
#         break

# update_fresh_streams()

facilities_sys.simulate()

# %%

import thermosteam as tmo

LacticAcid = tmo.Chemical('LacticAcid', Hfus=11.34e3)
print(LacticAcid.Hfus)




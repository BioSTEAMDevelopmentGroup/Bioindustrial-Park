#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Set up

from orgacids.system import *

# orgacids_sys.simulate()
# MPSP = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)
# print(MPSP)
# orgacids_sys.save_report('1.xlsx')

for i in range(0, 10):
    orgacids_sys.simulate()
    MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    print(MPSP)
    
# orgacids_sys.save_report('11.xlsx')


# %%
    
from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()
fermentation_sys.simulate()



# %% Separation system

# update_stripping_water()

# separation_sys_units_1 = (U401, M401, S401, R401)

# for unit in separation_sys_units_1:
#     try:
#         unit.simulate()
#         unit.show(N=100)
#     except: 
#         print(unit)
#         break

# update_separation_sulfuric_acid()

# separation_sys_units_2 = (T401, M402, S402, 
#                           F401, H401,
#                           F402, H402, R402,
#                           S403,
#                           H403, R403, S404, 
#                           H404, S405, H405, M403)

# for unit in separation_sys_units_2:
#     try:
#         unit.simulate()
#         unit.show(N=100)
#     except: 
#         print(unit)
#         break

separation_sys.simulate()


# %% Wastewtaer treatment

# wastewater_sys_units_1 = (M501, WWT_cost, R501)

# for unit in wastewater_sys_units_1:
#     try:
#         unit.simulate()
#         unit.show(N=100)
#     except: 
#         print(unit)
#         break

# update_aerobic_input_streams()
                          
# aerobic_digestion_sys_units = (M502, R502, S501, S502, M503, M504, S503)

# for unit in aerobic_digestion_sys_units:
#     try:
#         unit.simulate()
#         unit.show(N=100)
#     except: 
#         print(unit)
#         break

# aerobic_digestion_sys.simulate()

# wastewater_sys_units_2 = (aerobic_digestion_sys, S504, M505)

# for unit in wastewater_sys_units_2:
#     try:
#         unit.simulate()
#         unit.show(N=100)
#     except: 
#         print(unit)
#         break

wastewater_sys.simulate()


# %% Facilities

facilities_sys_units_1 = (CIP, ADP, FWT, BT, J601,
                          S601, T601, T602, T603, T604, S602, T605,
                          P601, T606, T607, S603, M601, T608, P602, 
                          CT, J602)


for unit in facilities_sys_units_1:
    try:
        unit.simulate()
        unit.show(N=100)
    except: 
        print(unit)
        break

update_process_water()

CWP.simulate()
S604.simulate()
PWC.simulate()
update_discharged_water()
update_fresh_streams()

facilities_sys.simulate()


# %%

from biosteam.units import HXutility
from thermosteam import Stream, settings
settings.set_thermo(['Water', 'Ethanol'])
feed = Stream('feed', Water=200, Ethanol=200)
hx = HXutility('hx', ins=feed, outs='product', T=50+273.15,
              rigorous=True) # Ignore VLE
hx.simulate()
hx.show()


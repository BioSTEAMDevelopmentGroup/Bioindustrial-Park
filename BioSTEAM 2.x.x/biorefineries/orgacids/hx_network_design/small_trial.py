# -*- coding: utf-8 -*-
"""
Created on Sun May  3 21:58:36 2020

@author: sarangbhagwat
"""
from hx_network_design import Temperature_Interval_Method, Design_Network

min_app_T = 20

Flow_vector = [1,1,1,1,1]

# T_in_vector = [180, 150, 30, 80] 

T_in_vector = [240, 180, 280, 50, 110]
VF_in_vector = [0,0,0,0,0]


# T_out_vector = [60, 30, 135, 140]
T_out_vector = [40, 120, 80, 230, 150]
VF_out_vector = [0,0,0,0,0]

# Cp_vector = [3, 1, 2, 5]
Cp_vector = [4, 10, 2.5, 5, 20]
LH_vector = [0,0,0,0,0]

t_pinches, h, c, hst, cst, hi, ci= Temperature_Interval_Method(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T)

print('\n\nH, C')
print(h, c)

print(hst)
print(t_pinches)
print(cst)

print(hi)
print(ci)

matches_hs, matches_cs, Q_hot_side, Q_cold_side, unavailables, hot_util_load, cold_util_load = Design_Network(Flow_vector, T_in_vector, T_out_vector, VF_in_vector, VF_out_vector, Cp_vector, LH_vector, min_app_T)
print('\n\nFinal hot side matches')
print(matches_hs)
print(Q_hot_side)
print('\n\nFinal cold side matches')
print(matches_cs)

print(Q_cold_side)
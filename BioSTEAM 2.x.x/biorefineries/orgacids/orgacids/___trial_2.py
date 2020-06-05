#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 08:25:48 2020

@author: sarangbhagwat
"""


# %% Whole-system simulation

'''
biosteam v2.12.8
thermosteam v0.12.16
'''
from biosteam.units import HXprocess, HXutility
import biosteam as bst
import copy
# from Sarang.system import *
from orgacids.system import *
from hx_network_design.hx_network_design import Temperature_Interval_Method, Design_Network

find = bst.main_flowsheet
# broth_for_separation = find.stream.broth_for_separation
MPSP = 0
System = bst.System
System.use_least_squares_solution = False
System.maxiter = 1500
System.molar_tolerance = 0.006

MPSPs = []

trial_iters = 3
# bst.HeatUtility.heating_agents[-1].show()

for i in range(trial_iters):
    orgacids_sys.simulate()   
# for i in range(5):
    lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
print(f'after {trial_iters} iterations, MPSP = {MPSP}')
# orgacids_sys.save_report(f'after_{trial_iters}.xlsx')


# %% 

import biosteam as bst
from orgacids.system import *


















# %%
for i in range(trial_iters):
    print('\n%s / %s'%(i+1, trial_iters))
    # print(i)
    # try:
    # System.converge_method = 'fixed-point' 
    System.converge_method = 'aitken' 
    
    # orgacids_sys0.simulate()
    orgacids_sys.simulate()
    
    lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    
    MPSPs.append(MPSP)
    # orgacids_sys.reset_flows()
    print('\nMPSP = %s $/kg\n----------\n'%round(MPSP,3))
    # except:
        # System.converge_method = 'fixed-point' 
        
        # # orgacids_sys0.simulate()
        # orgacids_sys.simulate()
        
        # lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
        # MPSPs.append(MPSP)
    # if i%10==0:
    #     print('MPSP = %s $/kg'%MPSP)
        # print(M501.outs[0].F_mass)
        # print(M501.outs[0].show(N=100))
        # print([(a.ID, a.imol['H2O']) for a in M501.ins])
    # orgacids_chemicals.LacticAcid
    # except: 
    #     System.converge_method = 'Wegstein'
    #     orgacids_sys.simulate()
    # print(F401.outs[1].F_mass)
    # Right now MPSP changes a bit each time the system simulates due to separation_sys
    # starting F401, the outs change each simulation
    
    # for i in range(5):    

orgacids_sys.diagram('thorough')

separation_sys = System('separation_sys', list(orgacids_sub_sys['separation_sys']))

separation_sys.diagram('thorough')


# Unit-level mass balance
print('\n------ Unit-level mass balance ------\n')
problem_units = []
for i in orgacids_sys.units:
    try: 
        in_by_out = sum(i.mass_in)/sum(i.mass_out)
        # print('%s\t%s:\t%s'%(i.ID, type(i), in_by_out))
        print('%s:\t\t%s'%(i.ID, in_by_out))
        if abs(in_by_out - 1)>0.02:
            problem_units.append(i)
    except:
        print('%s:\t\t%s, %s'%(i.ID, sum(i.mass_in),sum(i.mass_out)))        


S404.show(N=100, flow = 'kg/hr', composition = True)
# tot_mass = 0

# for species in S404.outs[1].chemicals:
#     if species.phase_ref=='l' and S404.outs[1].imass[species.ID]>0:
        
#         print(species.ID)
#         print(species.Tb)
#         if not species.ID == 'H2SO4':
#             tot_mass +=S404.outs[1].imass[species.ID]
# LA_purity = round(100*S404.outs[1].imass['LacticAcid']/tot_mass)

# Product stream
LA_purity = 100*S404.outs[1].imass['LacticAcid']/S404.outs[1].F_mass
LA_yield = 100*S404.outs[1].imol['LacticAcid']/(M401.outs[0].imol['LacticAcid'] + 2*M401.outs[0].imol['CalciumLactate'])
print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(round(LA_purity, 2), round(LA_yield, 2), round(MPSP,3)))
#

if problem_units:
    print('Problem units:\n')
    print(problem_units)

# Installation, purchase, and annual utility costs
sep_inst_wo_HXN, sep_pur_wo_HXN, sep_ut_wo_HXN = 0, 0, 0
for unit in separation_sys.units:
    sep_inst_wo_HXN += unit.installation_cost
    sep_pur_wo_HXN += unit.purchase_cost
    sep_ut_wo_HXN += unit.utility_cost

sep_ut_wo_HXN *= 365*24
tot_inst, tot_pur, tot_ut = \
    orgacids_tea.installation_cost, orgacids_tea.purchase_cost, orgacids_tea.utility_cost
HXN_inst, HXN_pur, HXN_ut = \
    HXN.installation_cost, HXN.purchase_cost, 365*24*HXN.utility_cost

tot_inst_wo_HXN, tot_pur_wo_HXN, tot_ut_wo_HXN = \
    tot_inst - HXN_inst, tot_pur - HXN_pur, tot_ut - HXN_ut
#

print('\n------ Cost contributions ------\n')

print('\nI. Savings from HXN implementation:\n')

print('\tInstallation: %s p.c.'%round(100*(1-tot_inst/tot_inst_wo_HXN),2))
print('\tPurchase: %s p.c.'%round(100*(1-tot_pur/tot_pur_wo_HXN),2))
print('\tUtility: %s p.c.'%round(100*(1-tot_ut/tot_ut_wo_HXN),2))


# print('\nII. Contribution of separation to orgacids_sys w/o HXN:\n')

# print('\tInstallation: %s p.c.'%round(100*sep_inst_wo_HXN/tot_inst_wo_HXN,2))
# print('\tPurchase: %s p.c.'%round(100*sep_pur_wo_HXN/tot_pur_wo_HXN,2))
# print('\tUtility: %s p.c.'%round(100*sep_ut_wo_HXN/tot_ut_wo_HXN,2))

# print('\nIII. Contribution of separation to orgacids_sys w/ HXN:\n')

# print('\tInstallation: %s p.c.'%round(100*sep_inst/tot_inst,2))
# print('\tPurchase: %s p.c.'%round(100*sep_pur/tot_pur,2))
# print('\tUtility: %s p.c.'%round(100*sep_ut/tot_ut,2))

print('\n\n--------------------------------')

# sys_outs = []
# sys_ins = []

# sum_mass_sys_outs = 0
# sum_mass_sys_ins = 0

# for i in find.stream:
#     if i.sink and not i.source:
#         sys_ins.append(i)
#         sum_mass_sys_ins += i.F_mass
#     elif i.source and not i.sink:
#         sys_outs.append(i)
#         sum_mass_sys_outs += i.F_mass
        

################################


# print('--- Heat Exchanger Network ---')

# print('\n\nFinal hot side matches')
# print(matches_hs)
# print(Q_hot_side)
# print('\n\nFinal cold side matches')
# print(matches_cs)

# print(Q_cold_side)

# print('\n\n')
# print('Duties for selected streams w/o HXN: %s kJ/hr heat, %s kJ/hr cooling'%(init_hot_util, init_cold_util))
# print('\nDuties for selected streams w/ HXN: %s kJ/hr heat, %s kJ/hr cooling'%(hot_util_load, cold_util_load))

# print('\n\nNew duty : Original duty ratio is %s for heat, %s for cooling'%(hot_util_load/init_hot_util, cold_util_load/init_cold_util))

#################################







#!!! M203.heat_utilities raises error
#!!! save_report not saving diagrams
# orgacids_sys.save_report('orgacids_sys.xlsx')


lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
print(f'after {trial_iters+trial_iters} iterations, MPSP = {MPSP}')
# orgacids_sys.save_report(f'after_{trial_iters+trial_iters}.xlsx')


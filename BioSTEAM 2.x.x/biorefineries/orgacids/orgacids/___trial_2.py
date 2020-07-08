#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 08:25:48 2020

@author: sarangbhagwat
"""


# %% Whole-system simulation

'''
biosteam v2.19.4
thermosteam v0.19.2
'''

import biosteam as bst
from biosteam.utils import TicToc

def set_feedstock_succinic_acid_content(feedstock, SA_content):
    dry_mass = feedstock.F_mass - feedstock.imass['H2O']
    feedstock.imass['SuccinicAcid'] = SA_content * dry_mass
    feedstock.imass['Extract'] -= (feedstock.F_mass-feedstock.imass['H2O']) - dry_mass
    if any(feedstock.mass<0):
        raise ValueError(f'Succinic acid content of {SA_content*100:.0f}% dry weight is infeasible')


# %%

# print('\n-------- orgacids_sys_wo_HXN --------')

# from orgacids.system_wo_HXN import *
# bst.main_flowsheet.set_flowsheet(bst.Flowsheet('orgacids_wo_HXN'))

# System = bst.System
# System.converge_method = 'fixed-point'
# # System.converge_method = 'aitken' 
# System.maxiter = 1500
# System.molar_tolerance = 0.1

# MPSPs_wo_HXN = []
# timer_wo_HXN = TicToc('timer_wo_HXN')
# trial_iters = 10

# from orgacids.system_wo_HXN import feedstock as feedstock_wo_HXN
# set_feedstock_succinic_acid_content(feedstock_wo_HXN, 0.05)
# feedstock.show(N=100)

# for i in range(trial_iters):
#     print('\n%s / %s'%(i+1, trial_iters))
    
#     timer_wo_HXN.tic()
#     orgacids_sys_wo_HXN.simulate()
#     lactic_acid.price = MPSP_wo_HXN = orgacids_tea_wo_HXN.solve_price(lactic_acid, orgacids_no_BT_tea_wo_HXN)    
#     MPSPs_wo_HXN.append(MPSP_wo_HXN)
#     # orgacids_sys_wo_HXN.reset_flows()
    
#     print('\nMPSP = %s $/kg\n----------\n'%round(MPSP_wo_HXN,3))
#     print(f'\nSimulation time: {timer_wo_HXN.elapsed_time:.3f} sec')
    
# print('\n-------- end orgacids_sys_wo_HXN --------')    
# operating_days = orgacids_tea_wo_HXN.operating_days

# tot_inst_wo_HXN, tot_pur_wo_HXN, tot_ut_wo_HXN = \
#     orgacids_tea_wo_HXN.installed_cost, orgacids_tea_wo_HXN.purchase_cost, orgacids_tea_wo_HXN.utility_cost

# separation_sys_wo_HXN = System('separation_sys', list(orgacids_sub_sys_wo_HXN['separation_sys']))
# sep_inst_wo_HXN, sep_pur_wo_HXN, sep_ut_wo_HXN = 0, 0, 0
# for unit in separation_sys_wo_HXN.units:
#     sep_inst_wo_HXN += unit.installed_cost
#     sep_pur_wo_HXN += unit.purchase_cost
#     sep_ut_wo_HXN += unit.utility_cost
# sep_ut_wo_HXN *= operating_days*24


# %%

print('\n-------- orgacids_sys_with_HXN --------')

from orgacids.system import *
bst.main_flowsheet.set_flowsheet(bst.Flowsheet('orgacids'))
# find = bst.main_flowsheet

System = bst.System
System.converge_method = 'fixed-point'
# System.converge_method = 'aitken' 
System.maxiter = 1500
System.molar_tolerance = 0.1

MPSPs = []
timer = TicToc('timer')
trial_iters = 3

from orgacids.system import feedstock as feedstock
set_feedstock_succinic_acid_content(feedstock, 0.05)
feedstock.show(N=100)

for i in range(trial_iters):
    print('\n%s / %s'%(i+1, trial_iters))
    
    timer.tic()
    orgacids_sys.simulate()
    lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_no_BT_tea)    
    MPSPs.append(MPSP)
    # orgacids_sys.reset_flows()
    
    print('\nMPSP = %s $/kg\n----------\n'%round(MPSP,3))
    print(f'\nSimulation time: {timer.elapsed_time:.3f} sec')
    
separation_sys = System('separation_sys', list(orgacids_sub_sys['separation_sys']))

# separation_sys.diagram('thorough')


# Unit-level mass balance
print('\n------ Unit-level mass balance ------\n')
problem_units = []
for i in orgacids_sys.units:
    try: 
        out_by_in = sum(i.mass_out)/sum(i.mass_in)
        # print('%s\t%s:\t%s'%(i.ID, type(i), out_by_in))
        print('%s:\t\t%s'%(i.ID, out_by_in))
        if abs(out_by_in - 1)>0.02:
            problem_units.append(i)
    except:
        print('%s:\t\t%s, %s'%(i.ID, sum(i.mass_in),sum(i.mass_out)))        


#!!! S404 was renamed to D405
D405.show(N=100, flow = 'kg/hr', composition = True)
# tot_mass = 0

# for species in D405.outs[1].chemicals:
#     if species.phase_ref=='l' and D405.outs[1].imass[species.ID]>0:
        
#         print(species.ID)
#         print(species.Tb)
#         if not species.ID == 'H2SO4':
#             tot_mass +=D405.outs[1].imass[species.ID]
# LA_purity = round(100*D405.outs[1].imass['LacticAcid']/tot_mass)

# Product stream
LA_purity = 100*D405.outs[1].imass['LacticAcid']/D405.outs[1].F_mass
LA_yield = 100*D405.outs[1].imol['LacticAcid']/(R301.outs[0].imol['LacticAcid'] + 2*R301.outs[0].imol['CalciumLactate'])
print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(round(LA_purity, 2), round(LA_yield, 2), round(MPSP,3)))
#

if problem_units:
    print('Problem units:\n')
    print(problem_units)

operating_days = orgacids_tea.operating_days
# Installation, purchase, and annual utility costs
sep_inst, sep_pur, sep_ut = 0, 0, 0
for unit in separation_sys.units:
    sep_inst += unit.installed_cost
    sep_pur += unit.purchase_cost
    sep_ut += unit.utility_cost
sep_ut *= operating_days*24

tot_inst, tot_pur, tot_ut = \
    orgacids_tea.installed_cost, orgacids_tea.purchase_cost, orgacids_tea.utility_cost
HXN_inst, HXN_pur, HXN_ut = \
    HXN.installed_cost, HXN.purchase_cost, operating_days*24*HXN.utility_cost


lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_no_BT_tea)
print(f'\nafter {trial_iters} iterations, MPSP = {MPSP}')

print('\n-------- end orgacids_sys_with_HXN --------')   
# tot_inst_wo_HXN, tot_pur_wo_HXN, tot_ut_wo_HXN = \
#     tot_inst - HXN_inst, tot_pur - HXN_pur, tot_ut - HXN_ut
#

# Update this using costs of a system w/o HXN rather than subtracting HXN costs, since 
# CT uses less elecricity since it has to regenerate less cooling water; however, 
# # BT generates less steam.
# print('\n\n------ Cost contributions ------\n')

# print('\nI. Savings from HXN implementation:\n')

# print('\tInstallation: %s p.c.'%round(100*(1-tot_inst/tot_inst_wo_HXN),2))
# print('\tPurchase: %s p.c.'%round(100*(1-tot_pur/tot_pur_wo_HXN),2))
# print('\tUtility: %s p.c.'%round(100*(1-tot_ut/tot_ut_wo_HXN),2))

# format(number, '%.0')

# print('\nII. Contribution of separation to orgacids_sys w/o HXN:\n')

# print('\tInstallation: %s p.c.'%round(100*sep_inst_wo_HXN/tot_inst_wo_HXN,2))
# print('\tPurchase: %s p.c.'%round(100*sep_pur_wo_HXN/tot_pur_wo_HXN,2))
# print('\tUtility: %s p.c.'%round(100*sep_ut_wo_HXN/tot_ut_wo_HXN,2))

# print('\nIII. Contribution of separation to orgacids_sys w/ HXN:\n')

# print('\tInstallation: %s p.c.'%round(100*sep_inst/tot_inst,2))
# print('\tPurchase: %s p.c.'%round(100*sep_pur/tot_pur,2))
# print('\tUtility: %s p.c.'%round(100*sep_ut/tot_ut,2))

# print('\n------ end Cost contributions ------\n')

# sys_outs = []
# sys_ins = []

# sum_mass_sys_outs = 0
# sum_mass_sys_ins = 0

# for i in bst.main_flowsheet.stream:
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

# print('\n\nNew duty: Original duty ratio is %s for heat, %s for cooling'%(hot_util_load/init_hot_util, cold_util_load/init_cold_util))

#################################





# orgacids_sys.save_report(f'after_{trial_iters+trial_iters}.xlsx')

# =============================================================================
# Name_change (both old and new names should work)
# =============================================================================

# name_dict = {
#     'H401': F401_H,
#     'S4ex': D401,
#     'pre_S403': D402,
#     'H404': D402_H,
#     'S403': D403,
#     'H402': D403_H,
#     'Split_S403': S403,
#     'pre_S404': D404,
#     'H403': D404_H,
#     'F_pre_S404': F402,
#     'H_pre_S404': F402_H,
#     'S404': D405,
#     'H_S404': D405_H,
#     'boiler_sys': BT_sys,
#     'orgacids_no_boiler_sys_tea': orgacids_no_BT_tea,
#     'boiler_sys_tea': BT_tea
#     }


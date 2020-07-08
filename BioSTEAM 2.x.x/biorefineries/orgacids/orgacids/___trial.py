#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 09:44:46 2020

@author: yalinli_cabbi
"""

import biosteam as bst
from biosteam.utils import TicToc
from orgacids.models import set_feedstock_carbs, set_feedstock_succinic_acid_content

trial_iters = 5


# %%

print('\n-------- START --------')   

from orgacids.system import *
bst.main_flowsheet.set_flowsheet(bst.Flowsheet('orgacids'))

System = bst.System
System.converge_method = 'fixed-point'
# System.converge_method = 'aitken' 
System.maxiter = 1500
System.molar_tolerance = 0.1

MPSPs = []
timer = TicToc('timer')

from orgacids.system import feedstock as feedstock
# set_feedstock_succinic_acid_content(feedstock, 0.05)
# set_feedstock_carbs(0.5899)

feedstock.show(N=100,flow='kg/hr', composition=True)

# %%

bst.speed_up()

for i in range(trial_iters):
    print('\n%s / %s'%(i+1, trial_iters))
    
    timer.tic()
    orgacids_sys.simulate()
    lactic_acid.price = MPSP = orgacids_tea.solve_price(lactic_acid, orgacids_no_BT_tea)
    MPSPs.append(MPSP)
    
    print('\nMPSP = %s $/kg\n----------\n'%round(MPSP,3))
    print(f'\nSimulation time: {timer.elapsed_time:.3f} sec')

LA_purity = 100*D405.outs[1].imass['LacticAcid']/D405.outs[1].F_mass
LA_yield = 100*D405.outs[1].imol['LacticAcid']/(R301.outs[0].imol['LacticAcid'] + 2*R301.outs[0].imol['CalciumLactate'])
# print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'% \
#       (round(LA_purity, 2), round(LA_yield, 2), round(MPSP,3)))

if 'HXN' in (i.ID for i in orgacids_sys.facilities):
    scenario = 'With HXN'
else:
    scenario = 'Without HXN'
print(f'\n{scenario}:')
print(f'\nPurity = {LA_purity/100:.2%} mass')
print(f'\nSeparation yield = {LA_yield/100:.2%} mass')
print(f'\nMPSP = ${MPSP:.3f}/kg')
print(f'\nAdjusted MPSP = ${MPSP/LA_purity*100:.3f}/kg')

print('\n-------- END --------')  







#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Whole-system simulation

'''
biosteam v2.17.7
thermosteam v0.17.15
'''

import biosteam as bst
# from Sarang.system import *
from orgacids.system import *


# %% MPSP

System.converge_method = 'Fixed-point'
System.maxiter = 1000
System.molar_tolerance = 1

# orgacids_sys.diagram()

def simulate_and_print():
    orgacids_sys.simulate()
    MPSP = lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)
    LA_purity = round(100*S404.outs[1].imass['LacticAcid']/S404.outs[1].F_mass)
    LA_yield = round(100*S404.outs[1].imol['LacticAcid']/(M401.outs[0].imol['LacticAcid'] + 2*M401.outs[0].imol['CalciumLactate']))
    # print(f'for lactic yield = {R301.cofermentation_rxns.X[0]}')
    # print('\nR301 cofermentation_rxns')
    # R301.cofermentation_rxns.show()
    # print('\nR302 cofermentation_rxns')
    # R302.cofermentation_rxns.show()
    print('\n\nPurity = %s percent mass \t Separation yield = %s percent mass \t MPSP = %s $/kg\n\n'%(LA_purity, LA_yield, round(MPSP,3)))

simulate_and_print()
# orgacids_sys.save_report('orgacids_sys.xlsx')


# %% Sugar concentration

from orgacids.utils import get_sugar_conc
sugar_conc = get_sugar_conc(R301.saccharified_stream, ('Glucose','Xylose'))
print(f'Sugar concentration of R301 after saccharification is {sugar_conc}')



# %%

import thermosteam as tmo
from orgacids import units
from orgacids.chemicals import orgacids_chemicals

tmo.settings.set_thermo(orgacids_chemicals)
a = T202.ins[0].copy()

# T1 = bst.units.MixTank('T1',tau=)


F1 = units.OrganicAcidsFlash('F1', ins=a, Q=0, P=101325)

























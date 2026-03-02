# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 23:27:48 2025

@author: saran
"""
import numpy as np

from biorefineries.isobutanol.system import corn_EtOH_IBO_sys as sys

from biorefineries.isobutanol.lumped_yeast_glucose_ethanol_isobutanol import MW_C5H7O2N_by_5, MW_Yeast

from matplotlib import pyplot as plt

f = sys.flowsheet

#%%

S404, V405 = f.S404, f.V405

steps = 100
# E_per_Cs = np.linspace(494_000/2, 494_000*2, steps) * MW_C5H7O2N_by_5 / MW_Yeast
umaxes = np.linspace(2e-6, 0.333e-3, steps)
MPSPs = []
TCIs = []

rxn_sys = V405.nsk_reaction_sys
# for E_per_C in E_per_Cs:
for umax in umaxes:
    print(umax)
    rxn_sys.reactions[0].reaction.rate_params['umax'] = umax
    
    # print(E_per_C)
    # rxn_sys.reactions[1].E_per_C = E_per_C
    
    S404.split = 0.99
    sys.simulate()
    f.isobutanol.price = 1.43 * f.isobutanol.imass['Isobutanol']/f.isobutanol.F_mass
    res_A = sys.TEA.solve_price(f.ethanol)*f.ethanol.F_mass/f.ethanol.imass['Ethanol'], sys.TEA.TCI
    
    # S404.split = 0.01
    sys.simulate()
    f.isobutanol.price = 1.43 * f.isobutanol.imass['Isobutanol']/f.isobutanol.F_mass
    res_B = sys.TEA.solve_price(f.ethanol)*f.ethanol.F_mass/f.ethanol.imass['Ethanol'], sys.TEA.TCI
    
    if res_A[0]<res_B[0]:
        MPSPs.append(res_A[0])
        TCIs.append(res_A[1])
    else:
        MPSPs.append(res_B[0])
        TCIs.append(res_B[1])
  
        
#%%
# plt.plot(E_per_Cs[:len(MPSPs)], MPSPs)
plt.plot(umaxes*3600.0, MPSPs)
plt.xlabel('Max. specific growth rate [M/h]')
plt.ylabel('Ethanol MPSP [$/kg]')

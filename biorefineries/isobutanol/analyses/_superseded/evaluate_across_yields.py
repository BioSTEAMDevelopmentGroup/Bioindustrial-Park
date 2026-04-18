# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 23:27:48 2025

@author: saran
"""
import numpy as np

from biorefineries.isobutanol.system import corn_EtOH_IBO_sys as sys

from matplotlib import pyplot as plt

f = sys.flowsheet

#%%

S404, V405 = f.S404, f.V405

steps = 10
y_min = 0.01
y_total = 0.95
ys_IBO = np.linspace(y_min, y_total-y_min, steps)
MPSPs = []
TCIs = []
for y in ys_IBO:
    V405.fermentation_reactions[0].X = y_total-y
    V405.fermentation_reactions[1].X = y
    
    S404.split = 0.99
    sys.simulate()
    f.isobutanol.price = 1.43 * f.isobutanol.imass['Isobutanol']/f.isobutanol.F_mass
    res_A = sys.TEA.solve_price(f.ethanol)*f.ethanol.F_mass/f.ethanol.imass['Ethanol'], sys.TEA.TCI
    
    S404.split = 0.01
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
plt.plot(ys_IBO[:8], MPSPs[:8])

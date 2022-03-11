# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:52:52 2022

@author: sarangbhagwat
"""

from biorefineries.TAL.system_TAL_adsorption_glucose import *
from matplotlib import pyplot as plt 
import numpy as np

#%%
def fun(v, t):
    AC1.mean_velocity=v
    AC2.mean_velocity=v
    
    AC1.cycle_time = t
    AC2.cycle_time = t
    return get_SA_MPSP()

mean_vels = np.linspace(1., 20., 20)
cycle_times = np.linspace(0.5, 4., 20)
MPSPs = []

for i in mean_vels:
    MPSPs.append([])
    for j in cycle_times:
        MPSPs[-1].append(fun(i, j))

plt.contourf(mean_vels, cycle_times, MPSPs)
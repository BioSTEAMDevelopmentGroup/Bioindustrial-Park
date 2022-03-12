# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:52:52 2022

@author: sarangbhagwat
"""

from biorefineries.TAL.system_TAL_adsorption_glucose import *
from matplotlib import pyplot as plt 
import numpy as np

#%% Across adsorption design space
def MPSP_at_adsorption_design(v, t):
    AC1.regeneration_velocity = v
    AC1.cycle_time = t
    return get_SA_MPSP()

regen_vels = np.linspace(1., 20., 10)
cycle_times = np.linspace(0.5, 4., 10)
MPSPs_ads_ds = []
#%%
for i in regen_vels:
    MPSPs_ads_ds.append([])
    for j in cycle_times:
        MPSPs_ads_ds[-1].append(MPSP_at_adsorption_design(i, j))
#%% Set parameters to optimal
min_MPSP = np.min(MPSPs_ads_ds)
opt_indices = np.where(MPSPs_ads_ds==min_MPSP)
AC1.regeneration_velocity = regen_vels[opt_indices[0][0]]
AC1.cycle_time = cycle_times[opt_indices[1][0]]
print(min_MPSP, get_SA_MPSP())
#%% Plot
plt.contourf(regen_vels, cycle_times, MPSPs_ads_ds)

#%% Across titer without adsorption design optimization
def MPSP_at_titer(t):
    spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
    AC1.regeneration_velocity = 3. + (17./25.)*t
    return get_SA_MPSP()

titers = np.linspace(2., 25., 10)

MPSPs_titer = []

#%%
for i in titers:
    MPSPs_titer.append(MPSP_at_titer(i))

#%% Plot
plt.plot(titers, MPSPs_titer)


#%% Across titer with adsorption design optimization

regen_vels = np.linspace(1., 20., 20)
cycle_times = np.linspace(0.5, 4., 20)

opt_regen_vels = []
opt_cycle_times = []
def MPSP_at_titer(t):
    spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
    MPSPs_ads_ds = []

    for i in regen_vels:
        MPSPs_ads_ds.append([])
        for j in cycle_times:
            MPSPs_ads_ds[-1].append(MPSP_at_adsorption_design(i, j))
    min_MPSP = np.min(MPSPs_ads_ds)
    opt_indices = np.where(MPSPs_ads_ds==min_MPSP)
    
    opt_regen_vels.append(regen_vels[opt_indices[0][0]])
    opt_cycle_times.append(cycle_times[opt_indices[1][0]])
    
    AC1.regeneration_velocity = opt_regen_vels[-1]
    AC1.cycle_time = opt_cycle_times[-1]
    print('titer =', t)
    print(min_MPSP, get_SA_MPSP(), AC1.ins[1].F_mass, AC1.regeneration_velocity, AC1.cycle_time)
    print('\n')
    return get_SA_MPSP()

titers = np.linspace(2., 25., 20)

MPSPs_titer = []

#%%
for i in titers:
    MPSPs_titer.append(MPSP_at_titer(i))

#%% Plot
plt.plot(titers, MPSPs_titer)



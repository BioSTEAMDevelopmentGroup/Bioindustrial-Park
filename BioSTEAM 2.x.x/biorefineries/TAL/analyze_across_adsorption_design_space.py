# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:52:52 2022

@author: sarangbhagwat
"""

from biorefineries.TAL.system_TAL_adsorption_glucose import *
from matplotlib import pyplot as plt 
import numpy as np

column = AC401

#%% Across adsorption design space
def MPSP_at_adsorption_design(v, t):
    column.regeneration_velocity = v
    column.cycle_time = t
    return get_SA_MPSP(), AC401.installed_cost/1e3

regen_vels = np.linspace(1., 20., 10)
cycle_times = np.linspace(0.5, 4., 10)
MPSPs_ads_ds = []
column_costs_ads_r_t = []
#%%
for i in regen_vels:
    MPSPs_ads_ds.append([])
    column_costs_ads_r_t.append([])
    for j in cycle_times:
        MPSP, cost = MPSP_at_adsorption_design(i, j)
        MPSPs_ads_ds[-1].append(MPSP)
        column_costs_ads_r_t[-1].append(cost)
#%% Set parameters to optimal
min_MPSP = np.min(MPSPs_ads_ds)
opt_indices = np.where(MPSPs_ads_ds==min_MPSP)
column.regeneration_velocity = regen_vels[opt_indices[0][0]]
column.cycle_time = cycle_times[opt_indices[1][0]]
print(min_MPSP, get_SA_MPSP())
#%% Plot MPSP
fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(regen_vels, cycle_times, MPSPs_ads_ds)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_xlabel('Regeneration solvent velocity [m/s]')
ax2.set_ylabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('MPSP [$/kg]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)
#%% Plot column cost
fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(regen_vels, cycle_times, column_costs_ads_r_t)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_xlabel('Regeneration solvent velocity [m/s]')
ax2.set_ylabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('Column installed cost [10^3 USD]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)


#%% Across titer without adsorption design optimization
def MPSP_at_titer(t):
    spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
    column.regeneration_velocity = 3. + (17./25.)*t
    return get_SA_MPSP()

titers = np.linspace(2., 25., 10)

# MPSPs_titer = []

#%%
MPSPs_titer = []
for i in titers:
    MPSPs_titer.append(MPSP_at_titer(i))

#%% Plot
plt.plot(titers, MPSPs_titer)


#%% Across titer with adsorption design optimization

# regen_vels = np.linspace(1., 20., 20)
# cycle_times = np.linspace(0.5, 4., 20)

# opt_regen_vels = []
# opt_cycle_times = []
# def MPSP_at_titer(t):
#     spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
#     MPSPs_ads_ds = []

#     for i in regen_vels:
#         MPSPs_ads_ds.append([])
#         for j in cycle_times:
#             MPSPs_ads_ds[-1].append(MPSP_at_adsorption_design(i, j))
#     min_MPSP = np.min(MPSPs_ads_ds)
#     opt_indices = np.where(MPSPs_ads_ds==min_MPSP)
    
#     opt_regen_vels.append(regen_vels[opt_indices[0][0]])
#     opt_cycle_times.append(cycle_times[opt_indices[1][0]])
    
#     column.regeneration_velocity = opt_regen_vels[-1]
#     column.cycle_time = opt_cycle_times[-1]
#     print('titer =', t)
#     print(min_MPSP, get_SA_MPSP(), column.ins[1].F_mass, column.regeneration_velocity, column.cycle_time)
#     print('\n')
#     return get_SA_MPSP()

# titers = np.linspace(2., 25., 20)

# MPSPs_titer = []

# #%%
# for i in titers:
#     MPSPs_titer.append(MPSP_at_titer(i))

# #%% Plot
# plt.plot(titers, MPSPs_titer)


#%%
AC401.target_recovery=None
mean_velocities = np.linspace(4., 15., 10)
cycle_times = np.linspace(0.1, 2, 10)
MPSPs = []
column_costs = []
for m in mean_velocities:
    AC401.mean_velocity = m
    MPSPs.append([])
    column_costs.append([])
    for t in cycle_times:
        AC401.cycle_time = t
        MPSPs[-1].append(get_SA_MPSP())
        column_costs[-1].append(AC401.installed_cost/1e3)

#%% Plot column cost
# plt.contourf(mean_velocities, cycle_times, MPSPs)

fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(mean_velocities, cycle_times, column_costs)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_xlabel('Mean feed velocity [m/s]')
ax2.set_ylabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('Column installed cost [10^3 USD]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

#%% Plot MPSP
# plt.contourf(mean_velocities, cycle_times, MPSPs)

fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(mean_velocities, cycle_times, MPSPs)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_xlabel('Mean feed velocity [m/s]')
ax2.set_ylabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('MPSP [$/kg]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)
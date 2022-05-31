# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 13:52:52 2022

@author: sarangbhagwat
"""

from biorefineries.TAL.system_TAL_adsorption_glucose import *
from matplotlib import pyplot as plt 
import numpy as np

column = AC401

#%% Across regeneration fluid velocity and cycle time
def MPSP_at_adsorption_design(v, t):
    column.regeneration_velocity = v
    column.cycle_time = t
    return get_SA_MPSP(), AC401.installed_cost/1e6

regen_vels = np.linspace(3., 20., 40)
cycle_times = np.linspace(1., 4., 40)
MPSPs_ads_ds = []
column_costs_ads_r_t = []
#%%
for i in regen_vels:
    MPSPs_ads_ds.append([])
    column_costs_ads_r_t.append([])
    for j in cycle_times:
        MPSP, cost = None, None
        try:
            MPSP, cost = MPSP_at_adsorption_design(i, j)
        except:
            print(i, j)
            MPSP, cost = np.nan, np.nan
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
CS = ax2.contourf(cycle_times, regen_vels, MPSPs_ads_ds, levels=[4., 4.5, 5, 5.5,  6., 6.5, 7.])

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_ylabel('Regeneration solvent velocity [m/s]')
ax2.set_xlabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('MPSP [$/kg]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)
#%% Plot column cost
fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(cycle_times, regen_vels, column_costs_ads_r_t, 
                   levels=[0, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5],
                  )

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_ylabel('Regeneration solvent velocity [m/s]')
ax2.set_xlabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('Column installed cost [10^6 USD]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

#%%
AC401.regeneration_velocity = 14.4
AC401.target_recovery=None
superficial_velocities = np.linspace(4., 15., 9)
cycle_times = np.linspace(1., 4., 10)
MPSPs = []
column_costs = []
for m in superficial_velocities:
    AC401.superficial_velocity = m
    MPSPs.append([])
    column_costs.append([])
    for t in cycle_times:
        AC401.cycle_time = t
        MPSPs[-1].append(get_SA_MPSP())
        column_costs[-1].append(AC401.installed_cost/1e6)

#%% Plot column cost
# plt.contourf(superficial_velocities, cycle_times, MPSPs)

fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(cycle_times, superficial_velocities, column_costs)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_ylabel('Superficial feed velocity [m/s]')
ax2.set_xlabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('Column installed cost [10^6 USD]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

#%% Plot MPSP
# plt.contourf(superficial_velocities, cycle_times, MPSPs)

fig1, ax2 = plt.subplots(constrained_layout=True)
CS = ax2.contourf(cycle_times, superficial_velocities, MPSPs)

CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# ax2.set_title('Nonsense (3 masked regions)')
ax2.set_ylabel('Superficial feed velocity [m/s]')
ax2.set_xlabel('Cycle time [h]')

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel('MPSP [$/kg]')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)


#%% Across titer
AC401.regeneration_velocity = 14.4
AC401.target_recovery = 0.99
AC401.cycle_time = 2.

titers = np.linspace(2., 25., 10)
MPSPs_titer_only = []
costs_titer_only = []

for t in titers:
    spec.load_specifications(spec_1=spec.baseline_yield, spec_2=t, spec_3=spec.baseline_productivity)
    MPSPs.append(get_SA_MPSP())
    costs_titer_only.append(AC401.installed_cost)

spec.load_specifications(spec_1=spec.baseline_yield, spec_2=spec.baseline_titer, spec_3=spec.baseline_productivity)

#%% Plot MPSP
plt.plot(titers, MPSPs_titer_only)

#%% Plot column cost
plt.plot(titers, costs_titer_only)
#%% Across titer and target recovery


# AC401.regeneration_velocity = 14.4
# AC401.target_recovery = 0.99

# # def MPSP_at_titer(t):
# #     spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
# #     column.regeneration_velocity = 3. + (17./25.)*t
# #     return get_SA_MPSP()

# titers = np.linspace(2., 25., 10)
# recoveries = np.linspace(0.5, 0.99, 10)
# # MPSPs_titer = []

# #%%
# MPSPs_titer = []
# costs_titer = []
# for t in titers:
#     MPSPs_titer.append([])
#     costs_titer.append([])
#     for r in recoveries:
#         spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
#         AC401.target_recovery = r
#         MPSPs_titer[-1].append(get_SA_MPSP())
#         costs_titer[-1].append(AC401.installed_cost)

# spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)


# #%% Plot MPSP
# fig1, ax2 = plt.subplots(constrained_layout=True)
# CS = ax2.contourf(recoveries, titers, MPSPs_titer, 
#                   # levels=[0., 2.5, 5., 7.5, 10, 12.5, 15, 17.5, 20],
#                   )

# CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# # ax2.set_title('Nonsense (3 masked regions)')
# ax2.set_ylabel('Fermentation titer [g/L]')
# ax2.set_xlabel('Target adsorbate recovery [% of influent]')

# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = fig1.colorbar(CS)
# cbar.ax.set_ylabel('MPSP [$/kg]')
# # Add the contour line levels to the colorbar
# cbar.add_lines(CS2)

# #%% Plot column cost
# fig1, ax2 = plt.subplots(constrained_layout=True)
# CS = ax2.contourf(recoveries, titers, costs_titer, 
#                   # levels=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18 ,20],
#                   )

# CS2 = ax2.contour(CS, levels=CS.levels[::2], colors='black', origin='lower')

# # ax2.set_title('Nonsense (3 masked regions)')
# ax2.set_ylabel('Regeneration solvent velocity [m/s]')
# ax2.set_xlabel('Cycle time [h]')

# # Make a colorbar for the ContourSet returned by the contourf call.
# cbar = fig1.colorbar(CS)
# cbar.ax.set_ylabel('Column installed cost [10^6 USD]')
# # Add the contour line levels to the colorbar
# cbar.add_lines(CS2)

#%% Across titer with rigorous adsorption design optimization

AC401.regeneration_velocity = 14.4
AC401.target_recovery = 0.99
AC401.cycle_time = 2.

regen_vels = np.linspace(1., 14.4, 20)
# cycle_times = np.linspace(0.5, 4., 20)

opt_regen_vels = []
opt_cycle_times = []

def MPSP_and_cost_at_regen_vel(v):
    column.regeneration_velocity = v
    return get_SA_MPSP(), AC401.installed_cost/1e6

def MPSP_at_titer(t):
    spec.load_specifications(spec_1=spec.spec_1, spec_2=t, spec_3=spec.spec_3)
    MPSPs_ads_ds = []
    costs_ads_ds = []

    for i in regen_vels:
        m, c = MPSP_and_cost_at_regen_vel(i)
        MPSPs_ads_ds.append(m)
        costs_ads_ds.append(c)
        
    min_MPSP = np.min(MPSPs_ads_ds)
    opt_indices = np.where(MPSPs_ads_ds==min_MPSP)
    
    opt_regen_vels.append(regen_vels[opt_indices[0][0]])
    # opt_cycle_times.append(cycle_times[opt_indices[1][0]])
    
    column.regeneration_velocity = opt_regen_vels[-1]
    # column.cycle_time = opt_cycle_times[-1]
    print('titer =', t)
    print(min_MPSP, column.ins[1].F_mass, column.regeneration_velocity, column.cycle_time)
    print('\n')
    return min_MPSP

titers = np.linspace(3., 30, 20)


#%%
MPSPs_titer = []
for i in titers:
    MPSPs_titer.append(MPSP_at_titer(i))
    
spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)

#%% Plot MPSP
plt.plot(titers, MPSPs_titer)

#%% Plot optimum regeneration velocity
plt.plot(titers, opt_regen_vels)
#%% Plot 
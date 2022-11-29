#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:59:00 2022

@author: jayneallen
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
import math
from flexsolve import IQ_interpolation
from matplotlib import pyplot as plt
import pandas as pd

BatchCrystallizer = bst.BatchCrystallizer

ln = math.ln
exp = math.exp

#%% Class

class SuccinicAcidCrystallizer(BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), 
                 target_recovery=0.6,
                 thermo=None,
                 tau=6, N=None, V=None, T=305.15,
                 Nmin=2, Nmax=36, 
                 T_range=(274., 372.5),
                 vessel_material='Stainless steel 316',
                 kW=0.00746):
        
        BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                     tau, N, V, T,
                     Nmin, Nmax, vessel_material,
                     kW)
        self.target_recovery = target_recovery
        self.T_range = T_range
        self.tau = tau
        
    def get_T_from_target_recovery(self, target_recovery):
        # self.target_recovery=target_recovery
        
        in_stream = self.ins[0]
        x_start = in_stream.imol['SuccinicAcid']/sum(in_stream.imol['Water', 'SuccinicAcid'])
        # print(x_start)
        x_end = (1-target_recovery)*x_start
        # T = self.get_T_from_x_end(x_end)
        water_vol = self.ins[0].ivol['Water']
        SA_mass = self.ins[0].imass['SuccinicAcid']
        SA_remaining = (1-target_recovery)*SA_mass
        x_end_gpL = SA_remaining/water_vol
        # print(x_end,x_end_gpL)
        return math.ln(x_end_gpL/29.098)/0.0396 + 273.15
        # return T
    
    def set_effluent_composition_from_recovery(self, recovery):
        in_stream = self.ins[0]
        self.outs[0].imass['s', 'SuccinicAcid'] = recovery*in_stream.imass['SuccinicAcid']
        self.outs[0].imass['l', 'SuccinicAcid'] = (1-recovery)*in_stream.imass['SuccinicAcid']
    
    def get_effective_recovery_from_T(self, T):
        in_stream = self.ins[0]
        out_stream = self.outs[0]
        SA_mass_dissolved_end = 29.098*exp(0.0396*(T-273.15)) * in_stream.F_vol
        SA_mass_total = self.ins[0].imass['SuccinicAcid']
        recovery = 1. - min(1., SA_mass_dissolved_end/SA_mass_total)
        return recovery
    
    def _run(self):
        in_stream, = self.ins
        out_stream, = self.outs
        target_recovery = self.target_recovery
        
        out_stream.copy_like(in_stream)
        # out_stream.sle(T=self.T, solute='SuccinicAcid')
        out_stream.phases=('l', 's')
        self.effective_recovery = target_recovery
        
        Tmin, Tmax = self.T_range
        rec_at_Tmin, rec_at_Tmax = self.get_effective_recovery_from_T(Tmin),\
            self.get_effective_recovery_from_T(Tmax)
        
        if rec_at_Tmin < target_recovery:
            self.T = Tmin
            self.effective_recovery = rec_at_Tmin
        elif rec_at_Tmax > target_recovery:
            self.T = Tmin
            self.effective_recovery = rec_at_Tmax
        else:
            self.T = T = self.get_T_from_target_recovery(target_recovery)
            if T>in_stream.T:
                self.T = T = in_stream.T
                self.effective_recovery = self.get_effective_recovery_from_T(T-273.15)
        
        # self.tau = self.get_t_from_target_recovery(effective_recovery)
        # print(self.ID, effective_recovery, self.effective_recovery)
        self.set_effluent_composition_from_recovery(self.effective_recovery)
        
        if self.effective_recovery>0.:
            out_stream.T =self.T
        else:
            self.T = out_stream.T = in_stream.T
            
#%% Run

SuccinicAcid = tmo.Chemical('SuccinicAcid')
tmo.settings.set_thermo(['Water', 'AceticAcid', SuccinicAcid])

input_stream = tmo.Stream('input_stream')
input_stream.T = 354.91 # K
input_stream.P = 50892 # Pa

input_stream.imass['Water'] = 2.33e+04
input_stream.imass['AceticAcid'] =  6.95e-06
# input_stream.imass['FermMicrobe'] = 0.409
input_stream.imass['SuccinicAcid'] = 437

target_SA_conc_gpL = 437
def SA_conc_obj_fn(SA_mol):
    input_stream.imol['SuccinicAcid'] = SA_mol
    return input_stream.imass['SuccinicAcid']/input_stream.F_vol - target_SA_conc_gpL

IQ_interpolation(SA_conc_obj_fn, 0, 1000)


C1 = SuccinicAcidCrystallizer(ID='C1', ins=input_stream, outs=('output_stream'), N=5)


S1 = bst.FakeSplitter('S1', ins=C1-0, outs=('solids', 'supernatant'))

@S1.add_specification()
def S1_spec():
    S1.outs[0].imol['SuccinicAcid'] = S1.ins[0].imol['s', 'SuccinicAcid'] # currently mock simulation to assume 0 water retained
    S1.outs[1].imol['Water'] = S1.ins[0].imol['s', 'Water'] + S1.ins[0].imol['l', 'Water'] 
    S1.outs[1].imol['SuccinicAcid'] =  S1.ins[0].imol['l', 'SuccinicAcid']


#%% Diagram

flowsheet_sys = bst.main_flowsheet.create_system('flowsheet_sys')
flowsheet_sys.simulate()
flowsheet_sys.diagram(kind='cluster', number=True, format='png')

#%% Get plot data
steps = (100, 100)

heat_utilities = []
power_utilities = []
utility_costs = []

heat_utilities_normalized = []

recoveries = np.linspace(0.5, 0.99, steps[0])
SA_concs = np.linspace(25, 400, steps[1])


for r in recoveries:
    C401.target_recovery = r
    heat_utilities_row = []
    power_utilities_row = []
    utility_costs_row = []
    heat_utilities_normalized_row = []
    for c in SA_concs:
        target_SA_conc_gpL = c
        IQ_interpolation(SA_conc_obj_fn, 0, 1000)
        input_stream.F_vol = 1
        try:
            for i in range(3):
                C401.simulate()
            heat_utilities_row.append(-C1.heat_utilities[0].duty)
            power_utilities_row.append(C1.power_utility.rate)
            utility_costs_row.append(C1.utility_cost)
            heat_utilities_normalized_row.append(-C1.heat_utilities[0].duty/(1e3*C1.outs[0].imass['s', 'SuccinicAcid']))
        except RuntimeError as e:
            if 'no cooling agent' in str(e):
                heat_utilities_row.append(np.inf)
                power_utilities_row.append(np.inf)
                utility_costs_row.append(np.inf)
                heat_utilities_normalized_row.append(np.inf)
            else:
                raise e
    heat_utilities.append(heat_utilities_row)
    power_utilities.append(power_utilities_row)
    utility_costs.append(utility_costs_row)
    heat_utilities_normalized.append(heat_utilities_normalized_row)
    
results_pd = pd.DataFrame(heat_utilities, recoveries, SA_concs)

#%% Plot cooling duty
x = np.array([SA_concs for i in range(len(recoveries))])
y = np.array([recoveries for i in range(len(SA_concs))]).transpose()
ax = plt.subplot()
im = ax.contourf(x, y, heat_utilities, cmap='viridis_r')
plt.ylabel('Target recovery [%]')
plt.xlabel('Initial SA concentration [g/L]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax, label='(-1) * Crystallizer cooling duty [kJ/h]')

plt.figure()

#%% Plot normalized cooling duty
x = np.array([SA_concs for i in range(len(recoveries))])
y = np.array([recoveries for i in range(len(SA_concs))]).transpose()
ax = plt.subplot()
im = ax.contourf(x, y, heat_utilities_normalized, cmap='viridis_r')
plt.ylabel('Target recovery [%]')
plt.xlabel('Initial SA concentration [g/L]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax, label='(-1) * Normalized crystallizer cooling duty [MJ/h/kg-SA-recovered]')

plt.figure()

#%% Plot power utility rate
x = np.array([SA_concs for i in range(len(recoveries))])
y = np.array([recoveries for i in range(len(SA_concs))]).transpose()
ax = plt.subplot()
im = ax.contourf(x, y, power_utilities, cmap='viridis_r')
plt.ylabel('Target recovery [%]')
plt.xlabel('Initial SA concentration [g/L]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax, label='Crystallizer power utility rate [kW]')

plt.figure()

#%% Plot utility cost
x = np.array([SA_concs for i in range(len(recoveries))])
y = np.array([recoveries for i in range(len(SA_concs))]).transpose()
ax = plt.subplot()
im = ax.contourf(x, y, utility_costs, cmap='viridis_r')
plt.ylabel('Target recovery [%]')
plt.xlabel('Initial SA concentration [g/L]')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)

plt.colorbar(im, cax=cax, label='Crystallizer utility cost [$/h]')

plt.figure()


#%% Code for 1-d plots
# #%%

# SA_concs = np.linspace(70, 500, 20)
# Ts = []
# heat_utilities = []
# power_utilities = []
# utility_costs = []
# F_vols = []
# for c in SA_concs:
#     target_SA_conc_gpL = c
#     IQ_interpolation(SA_conc_obj_fn, 0, 1000)
#     input_stream.F_vol = 1
#     C1.simulate()
#     Ts.append(C1.T)
#     heat_utilities.append(C1.heat_utilities[0].duty) # kJ/h
#     power_utilities.append(C1.power_utility.rate) # kW
#     utility_costs.append(C1.utility_cost) # $/h
#     F_vols.append(C1.ins[0].F_vol)

# # Plot
# plt.plot(SA_concs, heat_utilities)

# #%%
# recoveries = np.linspace(0.2, 0.9, 80)

# Ts = []
# heat_utilities = []
# power_utilities = []
# utility_costs = []
# F_vols = []

# for r in recoveries:
#     C1.target_recovery = r
#     C1.simulate()
#     Ts.append(C1.T)
#     heat_utilities.append(C1.heat_utilities[0].duty) # kJ/h
#     power_utilities.append(C1.power_utility.rate) # kW
#     utility_costs.append(C1.utility_cost) # $/h
#     F_vols.append(C1.ins[0].F_vol)

# # Plot
# plt.plot(recoveries, utility_costs)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:12:41 2022

@author: Jayne Allen, Sarang Bhagwat
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
import math
from matplotlib import pyplot as plt

ln = math.log
BatchCrystallizer = bst.BatchCrystallizer

#%% Class 
class SuccinicAcidCrystallizer(BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), 
                 target_recovery=0.6,
                 thermo=None,
                 tau=None, N=None, V=None, T=305.15,
                 Nmin=2, Nmax=36, vessel_material='Carbon steel',
                 kW=0.00746):
        
        BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                     tau, N, V, T,
                     Nmin, Nmax, vessel_material,
                     kW)
        self.target_recovery = target_recovery
    
    def get_T_from_x_end(self, x_end):
        # !!! model here
        # water_vol = self.ins[0].ivol['Water']
        # SA_mass = self.ins[0].imass['SuccinicAcid']
        # SA_remaining = (1-self.target_recovery)*SA_mass
        # x_end_gpL = SA_remaining/water_vol
        # print(x_end,x_end_gpL)
        # return ln(x_end_gpL/29.098)/0.0396 # deg K # mock output for now
        1
        
    def get_T_from_target_recovery(self, target_recovery):
        self.target_recovery=target_recovery
        
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
        return ln(x_end_gpL/29.098)/0.0396 + 273.15
        # return T
    
    def get_t_from_target_recovery(self, target_recovery):
        in_stream = self.ins[0]
        x_start = in_stream.imol['SuccinicAcid']/sum(in_stream.imol['Water', 'SuccinicAcid'])
        x_end = (1-target_recovery)*x_start
        # !!! model here
        return 4 # hours # mock output for now
    
    def set_effluent_composition_from_recovery(self, target_recovery):
        in_stream = self.ins[0]
        self.outs[0].imass['s', 'SuccinicAcid'] = target_recovery*in_stream.imass['SuccinicAcid']
        self.outs[0].imass['l', 'SuccinicAcid'] = (1-target_recovery)*in_stream.imass['SuccinicAcid']
        
    def _run(self):
        in_stream, = self.ins
        out_stream, = self.outs
        target_recovery = self.target_recovery
        self.T = self.get_T_from_target_recovery(target_recovery)
        self.tau = self.get_t_from_target_recovery(target_recovery)
        out_stream.copy_like(in_stream)
        out_stream.sle(T=self.T, solute='SuccinicAcid')
        self.set_effluent_composition_from_recovery(target_recovery)
            
        
#%% Run

SuccinicAcid = tmo.Chemical('SuccinicAcid')
# SuccinicAcid.Cn.l.add_method(153.1)
tmo.settings.set_thermo(['Water', SuccinicAcid])

input_stream = tmo.Stream('input_stream')
input_stream.T = 50 + 273.15

input_stream.imass['SuccinicAcid'] = 250
input_stream.imass['Water'] = 1000

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

N = 50
recoveries = np.linspace(0.55, 0.9, N)
initial_sa_molar_flow = np.linspace(100, 500, N)
Ts = []
heat_utilities = []
power_utilities = []
utility_costs = []

for i in initial_sa_molar_flow:
    for r in recoveries:
        C1.ins[0].imol['SuccinicAcid'] = i
        C1.target_recovery = r
        C1.simulate()
        Ts.append(C1.T)
        heat_utilities.append(C1.heat_utilities[0].duty) # kJ/h
        power_utilities.append(C1.power_utility.rate) # kW
        utility_costs.append(C1.utility_cost) # $/h

heat_utilities = np.array(heat_utilities)
power_utilities = np.array(power_utilities)
utility_costs = np.array(utility_costs)
z = heat_utilities.reshape(N, N)
a = power_utilities.reshape(N, N)
b = utility_costs.reshape(N, N)
x, y = np.meshgrid(recoveries, initial_sa_molar_flow)

#%% Plot

# Contour plots

contour_cooling_utilities = plt.contour(x, y, z)
plt.title('Crys. Cooling Utility Demand vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Cooling Utility Demand [kJ/h]')
cb.ax.tick_params(labelsize=14)
plt.clabel(contour_cooling_utilities, inline=1, fontsize=10)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize = 14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize = 14)
plt.show()

contour_power_utilities = plt.contour(x, y, a)
plt.title('Crys. Power Utility Rate vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Power Utility Rate [kW]')
cb.ax.tick_params(labelsize=14)
plt.clabel(contour_power_utilities, inline=1, fontsize=10)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize = 14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize = 14)
plt.show()

contour_utilities_costs = plt.contour(x, y, b)
plt.title('Crys. Utility Costs vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Utility Costs [$/h]')
cb.ax.tick_params(labelsize=14)
plt.clabel(contour_utilities_costs, inline=1, fontsize=10)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize = 14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize = 14)
plt.show()



######### Colormesh plots

plt.pcolormesh(x,y,z,cmap='jet',shading = 'auto')
plt.title('Crys. Cooling Utility Demand vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Cooling Utility Demand [kJ/h]')
cb.ax.tick_params(labelsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize=14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize=14)

plt.show()


plt.pcolormesh(x,y,a,cmap='jet',shading = 'auto')
plt.title('Crys. Power Utility Rate vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Power Utility Rate [kW]')
cb.ax.tick_params(labelsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize=14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize=14)

plt.show()


plt.pcolormesh(x,y,b,cmap='jet',shading = 'auto')
plt.title('Crys. Utility Costs vs. \n (Target Recovery, Initial SA Rate)', fontsize = 14)
cb = plt.colorbar(label = 'Utility Costs [$/h]')
cb.ax.tick_params(labelsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Target Recovery [mol/mol]', fontsize=14)
plt.ylabel('Initial SA Flow Rate [mol/hr]', fontsize=14)

plt.show()










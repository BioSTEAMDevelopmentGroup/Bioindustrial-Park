# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 16:50:17 2022

@author: sarangbhagwat
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np

BatchCrystallizer = bst.BatchCrystallizer

#%% Class 
class SuccinicAcidCrystallizer(BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), 
                 target_recovery=0.99,
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
        return 273.15 + 5 # deg K # mock output for now
    
    def get_T_from_target_recovery(self, target_recovery):
        in_stream = self.ins[0]
        x_start = in_stream.imol['SuccinicAcid']/sum(in_stream.imol['Water', 'SuccinicAcid'])
        x_end = (1-target_recovery)*x_start
        T = self.get_T_from_x_end(x_end)
        return T
    
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
tmo.settings.set_thermo(['Water', 'SuccinicAcid'])

input_stream = tmo.Stream('input_stream')
input_stream.imass['SuccinicAcid'] = 10
input_stream.imass['Water'] = 90

SAC = SuccinicAcidCrystallizer(ID='SAC', ins=input_stream, outs=('output_stream'), N=5)

SAC.simulate()

SAC.show('cwt100')

print(SAC.results())


    
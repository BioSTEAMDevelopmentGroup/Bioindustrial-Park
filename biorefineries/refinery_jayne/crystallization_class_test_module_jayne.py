#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 16:51:43 2022

@author: jayneallen
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
from biosteam import settings, units

BatchCrystallizer = bst.BatchCrystallizer

#%% Class

class SuccinicAcidCrystallizer(BatchCrystallizer):  #inheritance

    def __init__(self, ID='', ins=None, outs=(),
                 target_recovery=0.99,
                 thermo=None,
                 tau=None, N=None, V=None, T=305.15,
                 Nmin=2, Nmax=36, vessel_material='Carbon steel',
                 kW=0.00746):
        
        BatchCrystallizer.__init__(self, ID, ins, outs, thermo, #why repeat?
                                   tau, N, V, T, Nmin, Nmax, vessel_material,
                                   kW)
        self.target_recovery = target_recovery #why is this line separate?
        
        def get_T_from_x_end(self, x_end):
            # !!! model here
            return 273.15 + 5 #K
        
        def get_T_from_target_recovery(self, target_recovery):
            in_stream = self.ins[0]
            x_start = in_stream.imol['SuccinicAcid']/sum(in_stream.imol['Water', 'SuccinicAcid'])
            x_end = (1-target_recovery)*x_start
            T = self.get_t_from_x_end(x_end)
            return T
        
        def get_t_from_target_recovery(self, target_recovery):
            in_stream= self.ins[0]
            x_start = in_stream.imol['SuccinicAcid']/sum(in_Stream.imol['Water', 'SuccinicAcid'])
            x_end = (1-target_recovery)*x_start
            # !!! model here
            return 4 #hours
        
        def set_effluent_composition_from_recovery(self, target_recovery):
            in_stream = self.ins[0]
            self.outs[0].imass['s', 'SuccinicAcid'] = target_recovery*in_stream.imass['SuccinicAcid']
            self.outs[0].imass['l', 'SuccinicAcid'] = (1-target_recovery)*in_stream.imass['SuccinicAcid']
            
        def _run(self):  #run is last
            in_stream, = self.ins
            out_stream, = self.outs
            target_recovery = self.target_recovery
            self.T = self.get_T_from_target_recovery(target_recovery)
            self.tau = self.get_t_from_target_recovery(target_recovery)
            out_stream.copy_like(in_stream)
            out_stream.sle(T=self.T, solute='SuccinicAcid')
            self.set_effluent_composition_from_recovery(target_recovery)

#%% Streams

settings.set_thermo(['Glucose','Quicklime','Water','H2SO4'])

pure_glucose = bst.Stream(Glucose=1000.)
pure_glucose.show()

quicklime = bst.Stream(Quicklime=10000.,Water=100.)
quicklime.show()

sulfuric_acid = bst.Stream(H2SO4=30.)
sulfuric_acid.show()

elution_water = bst.Stream(Water=40.)
elution_water.show()

#%% Units

# Crystallization
C1 = SuccinicAcidCrystallizer('C1', ins=pure_glucose, outs='solids',
                              target_recovery=0.99, thermo=None,
                              tau=1, N=10,)

#%% Diagram

flowsheet_sys = bst.main_flowsheet.create_system('flowsheet_sys')
flowsheet_sys.simulate()













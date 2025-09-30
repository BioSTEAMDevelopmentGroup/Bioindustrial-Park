# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:48:06 2025

@author: IGB
"""

import biosteam as bst
import thermosteam as tmo
from biosteam.units import Compressor, StirredTankReactor, HXutility, Mixer
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

CEPCI = bst.design_tools.CEPCI_by_year

#%%
@cost
class PowerUnit(bst.Unit):
    _N_ins = 1
    _N_outs = 1
    
    _units = {}
    
    def _init(self, ID='', ins=None, outs=(), thermo=None):
        super()._init()
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
    
    def _cost(self):
        self.add_power_utility(self.power)




class PlasmaReactor(StirredTankReactor):
    
    _N_ins = 2
    _N_outs = 2
    
    auxiliary_unit_names = ('heat_exchanger')
    
    _F_BM_default = {**StirredTankReactor._F_BM_default}
    
    def _init(self, ID='', ins=None, outs=(), thermo=None, *, 
              tau=24,
              T=17+273.15,
              P=101325,
              vessel_material='Stainless steel 316',
              batch=True,
              rigorous_hx=True,
              **args):
        super()._init(tau=tau, T=T, P=P, batch=batch, vessel_material=vessel_material)
    
    def _run(self):
        air, water = self.ins
        vent, effluent = self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        
    
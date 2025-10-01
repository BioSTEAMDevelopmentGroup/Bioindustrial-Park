# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:48:06 2025

@author: IGB
"""

import biosteam as bst
import thermosteam as tmo
from biosteam.units.decorators import cost

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

CEPCI = bst.design_tools.CEPCI_by_year

#%%
class PowerUnit(bst.Unit): # only simulates the power demand and purchase cost
    _N_ins = 1
    _N_outs = 1
    
    def _init(self, ID='', ins=None, outs=(), thermo=None,
              power=None):
        super()._init()
        self.power=power
    
    def _run(self):
        product_in, = self.ins
        product_out, = self.outs
        product_out.copy_like(product_in)
    
    def _cost(self):
        self.add_power_utility(self.power)
        self.purchase_costs['Power unit'] = self.power * 1000 # $1000/kW



class PlasmaReactor(bst.StirredTankReactor):
    
    _N_ins = 2
    _N_outs = 2
    
    T_default = 17+273.15
    P_default = 101325
    
    def _init(self, ID='', ins=None, outs=(), thermo=None, *, 
              tau=24,
              batch=True,
              vessel_material='Stainless steel 316',
              N2_conversion=0.01,
              electricity_consumption_per_mol_N=15, # 15 MJ per mol N fixed
              electricity_to_heat_ratio=0.5, # assumed
              **args):
        super()._init(tau=tau, batch=batch, vessel_material=vessel_material)
        self.reaction = bst.Reaction('2N2 + 5O2 + 2H2O -> 4HNO3', 'N2', X=N2_conversion)
        self.electricity_consumption_per_mol_N = electricity_consumption_per_mol_N
        self.electricity_to_heat_ratio = electricity_to_heat_ratio
    
    def _run(self):
        air, water = self.ins
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.reaction(effluent)
        effluent.T = vent.T = self.T
        effluent.P = vent.P = self.P
        vent.phase = 'g'
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)
    
        self.power = effluent.imol['HNO3'] * 1000 * self.electricity_consumption_per_mol_N / 3.6 # mol * MJ/mol, in kW
        
    def _get_duty(self):
        self.heat_by_electricity = self.effluent.imol['HNO3'] * 1e6 * self.electricity_consumption_per_mol_N * self.electricity_to_heat_ratio
        self.Hnet + self.heat_by_electricity # kJ/hr
        
    
    
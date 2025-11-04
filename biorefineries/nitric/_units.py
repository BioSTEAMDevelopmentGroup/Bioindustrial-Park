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
class PowerUnit(bst.Unit): # only simulates the power
    _N_ins = 1
    _N_outs = 1
    
    def _init(self, ID='', ins=None, outs=(), thermo=None,
              cost_per_power=None, power=None):
        super()._init()
        self.cost_per_power = cost_per_power # $/W
        self.power = power # kW
    
    def _run(self):
        product_in, = self.ins
        product_out, = self.outs
        product_out.copy_like(product_in)
    
    def _design(self):
        # self.add_power_utility(self.power)
        self.power_utility(self.power)
        
    def _cost(self):
        self.purchase_costs['Power unit'] = self.power * 1000 * self.cost_per_power



# class PlasmaReactor(bst.StirredTankReactor):
    
#     _N_ins = 2
#     _N_outs = 2
    
#     T_default = 17+273.15
#     P_default = 101325
    
#     def _init(self, ID='', ins=None, outs=(), thermo=None, *, 
#               tau=None,
#               batch=True,
#               vessel_material='Stainless steel 316',
#               HNO3_scale=None, # kg/day
#               concentration=None, # 4778 NO3 mg/L=4.778g/L=4.778kg/m3
#               electricity_consumption=None, # 15 MJ per mol N fixed
#               electricity_to_heat_ratio=None, # assumed
#               air_scale_ratio=None,
#               **args):
#         super()._init(tau=tau, batch=batch, vessel_material=vessel_material)
#         self.HNO3_scale = HNO3_scale
#         self.concentration = concentration
#         self.electricity_consumption = electricity_consumption
#         self.electricity_to_heat_ratio = electricity_to_heat_ratio
#         self.air_scale_ratio = air_scale_ratio
    
#     def _run(self):
#         effluent = self.outs[1]
        
#         effluent.mix_from(self.ins)
#         effluent.imass['N2'] = 0
#         effluent.imass['O2'] = 0
#         effluent.imass['HNO3'] = self.HNO3_scale / 24
        
#         self.power = effluent.imol['HNO3'] * 1000 * self.electricity_consumption / 3.6 # mol * MJ/mol, in kW
        
#         effluent.T = self.T
#         effluent.P = self.P
    
#     def _get_duty(self):
#         self.heat_by_electricity = -self.effluent.imol['HNO3'] * 1e6 * self.electricity_consumption * self.electricity_to_heat_ratio # kJ
#         return self.heat_by_electricity


class PlasmaReactor(bst.Unit):
    _N_ins = 2
    _N_outs = 2
    
    # _F_BM_default = {'Plasma setup': 2.}
    
    def _init(self, ID='', ins=None, outs=(), thermo=None, HNO3_scale=None, concentration=None,
              electricity_consumption=None,
              cost_per_power=None, heat_by_electricity=None):
        super()._init()
        self.HNO3_scale = HNO3_scale
        self.concentration = concentration
        self.electricity_consumption = electricity_consumption
        self.cost_per_power = cost_per_power
        self.heat_by_electricity = heat_by_electricity
    
    def _run(self):
        effluent = self.outs[1]
        effluent.imass['HNO3'] = self.HNO3_scale / 24
        self.power = effluent.imol['HNO3'] * 1000 * self.electricity_consumption / 3.6 # mol * MJ/mol, in kW
    
    def _design(self):
        q = -self.outs[1].imol['HNO3'] * 1e6 * self.electricity_consumption * self.heat_by_electricity # kJ
        self.add_heat_utility(q, T_in=17+273.15, T_out=10+273.15)
            
        
    def _cost(self):
        self.purchase_costs['Plasma setup'] = self.power * 1000 * self.cost_per_power
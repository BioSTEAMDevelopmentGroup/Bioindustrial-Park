# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 14:17:48 2024

@author: IGB
"""


import numpy as np
import biosteam as bst
import thermosteam as tmo
from biosteam import Stream, Unit, BinaryDistillation
from biosteam.units import HXutility, Mixer, SolidsSeparator, Compressor
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
from thermosteam import MultiStream
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
from qsdsan._sanunit import SanUnit

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

#%% 
# =============================================================================
# H2 production through Liquid Alkaline (LA) electrolysis 
# =============================================================================
# Model H2 production by H2O => H2 + 0.5O2
class Electrolyzer(SanUnit): 
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID="", ins=None, outs=(),
                 stack_lifetime=10,
                 renewable_electricity_price=0.03, # currently possible in U.S. markets with plentiful  wind
                 **kwargs):
        bst.Unit.__init__(self, ID, ins, outs)
        self.stack_lifetime = stack_lifetime
        self.renewable_electricity_price = renewable_electricity_price
        
    def _run(self):
        water, = self.ins
        hydrogen, oxygen = self.outs
        hydrogen.phase = 'g'
        hydrogen.P = 30 * 101325 # Model H2 purification so that output is 99.99% purity and 30 bar
        oxygen.phase = 'g'
        hydrogen.imol['H2'] = water.F_mol
        oxygen.imol['O2'] = water.F_mol / 2
    
    def _cost(self): # based on https://www.osti.gov/biblio/2203367 P17
        self.power_utility.rate = self.outs[0].imass['H2'] * 122000 * 24 / 50000 # base case 50 MTD / 122 MW
        # self.power_utility.cost = self.renewable_electricity_price * self.power_utility.rate
        self.baseline_purchase_costs['Stack'] = self.power_utility.rate * 292 # process equipment, piping, valves, and  instrumentation including temperature, pressure, flow, and level indicators
        self.baseline_purchase_costs['Mechanical BOP'] = self.power_utility.rate * 280
        self.baseline_purchase_costs['Electrical BOP'] = self.power_utility.rate * 155 # for AC-> DC, transformer, power substation
        self.add_OPEX = self._calculate_replacement_cost()
    
    def _calculate_replacement_cost(self):
        stack_replacement_cost  = self.baseline_purchase_costs['Stack'] / self.stack_lifetime / (365*24) # convert from USD/year to USD/hour
        return stack_replacement_cost
    


    

# =============================================================================
# MeOH synthesis process
# =============================================================================

# Reactor is modelled as an adiabatic ideal plug flow reactor (PFR)
class MeOH_SynthesisReactor(bst.units.design_tools.PressureVessel, bst.Unit):
    _N_ins = _N_outs = 1
    
    _units = {**bst.design_tools.PressureVessel._units,
              'Volume': 'ft^3',
              'Catalyst weight': 'kg'}
    
    _F_BM_default = {**bst.design_tools.PressureVessel._F_BM_default,
                     'Catalyst': 1.0} # Cu/ZnO/Al2O3
    
    def __init__(self, ID="", ins=None, outs=(),
                 P=7600000/6894.76,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 length_to_diameter=0.15/0.016, # 10.1016/j.jclepro.2013.06.008
                 WHSV = 0.028 * 3600 / 34.8, # Residence time in kg/hr feed / kg catalyst; 10.1016/j.jclepro.2013.06.008
                 catalyst_density = 1775, # In kg/m^3; 10.1016/j.jclepro.2013.06.008
                 catalyst_price = 13,
                 porosity = 0.5,
                 **kwargs):
            bst.Unit.__init__(self, ID, ins, outs)
            self.P = P
            self.vessel_material = vessel_material # Vessel material
            self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
            self.length_to_diameter = length_to_diameter #: Length to diameter ratio
            self.WHSV = WHSV
            self.catalyst_density = catalyst_density
            self.catalyst_price = 13 # In USD/kg
            self.porosity = 0.5 #  fraction of the bed volume occupied by void space
            self.reactions = ParallelRxn([
                # Reaction definition                 Reactant    Conversion
                Rxn('CO2 + 3H2 -> CH3OH + H2O',       'CO2',      0.210),
                Rxn('CO2 + H2 -> CO + H2O',           'CO2',      0.004),])
                
    
    def _run(self):
        feed, = self.ins
        feed.phase = 'g'
        effluent, = self.outs
        effluent.copy_like(feed)
        self.reactions.adiabatic_reaction(effluent)
        
    def _design(self):
        feed, = self.ins
        effluent, = self.outs
        length_to_diameter = self.length_to_diameter
        P = effluent.get_property('P', 'psi')
        results = self.design_results
        catalyst = feed.F_mass / self.WHSV
        results['Catalyst weight'] = catalyst
        results['Volume'] = volume = catalyst / self.catalyst_density / (1-self.porosity) * 35.3147 # to ft3
        D = cylinder_diameter_from_volume(volume, length_to_diameter)
        L = length_to_diameter * D
        results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
        self.purchase_costs['Catalyst'] = self.catalyst_price * D['Catalyst weight']



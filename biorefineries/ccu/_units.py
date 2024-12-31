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

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

#%% 
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
                     'Catalyst': 1.0}
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 P=7600000,
                 vessel_material='Stainless steel 304',
                 vessel_type='Horizontal',
                 length_to_diameter=3.0,
                 **wargs):
            bst.Unit.__init__(self, ID, ins, outs, thermo)
            self.length_to_diameter = length_to_diameter #: Length to diameter ratio
            self.vessel_material = vessel_material # Vessel material
            self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
            self.catalyst_density = 1059.52 # In kg/m^3
            self.catalyst_price = 13 # In USD/kg
            self.P = P
            self.reactions = ParallelRxn([
                # Reaction definition                 Reactant    Conversion
                Rxn('CO2 + 3H2 -> CH3OH + H2O',       'CO2',      0.210),
                Rxn('CO2 + H2 -> CO + H2O',           'CO2',      0.004),])
                
    
    def _run(self):
        feed, = self.ins
        effluent, = self.outs
        effluent.P = self.P
        effluent.copy_like(feed)
        self.reactions.adiabatic_reaction(effluent)
        
    def _design(self):
        effluent, = self.outs
        length_to_diameter = self.length_to_diameter
        P = effluent.get_property('P', 'psi')
        results = self.design_results
        catalyst = effluent.F_mass / self.WHSV
        results['Catalyst weight'] = catalyst
        results['Volume'] = volume = catalyst / self.catalyst_density * 35.3147 # to ft3
        D = cylinder_diameter_from_volume(volume, length_to_diameter)
        L = length_to_diameter * D
        results.update(self._vessel_design(P, D, L))
        
    def _cost(self):
        D = self.design_results
        self.purchase_costs.update(
            self._vessel_purchase_cost(D['Weight'], D['Diameter'], D['Length'])
        )
        self.purchase_costs['Catalyst'] = self.catalyst_price * D['Catalyst weight']




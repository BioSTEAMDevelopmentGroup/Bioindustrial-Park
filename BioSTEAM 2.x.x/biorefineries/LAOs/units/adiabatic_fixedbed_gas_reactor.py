# -*- coding: utf-8 -*-
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Part of the BioSTEAM project. Under the UIUC open-source license.
# See github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
import biosteam as bst

__all__ = ('AdiabaticFixedbedGasReactor',)

class AdiabaticFixedbedGasReactor(bst.design_tools.PressureVessel, bst.Unit):
    _N_ins = _N_outs = 1
    _N_heat_utilities = 0
    _units = {**bst.design_tools.PressureVessel._units,
              'Volume': 'ft^3',
              'Catalyst weight': 'kg'}
    _BM = {**bst.design_tools.PressureVessel._BM,
           'Catalyst': 1.0}
    
    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 reaction, catalyst_density, catalyst_price, WHSV,
                 P, length_to_diameter=3.0,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical'):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.length_to_diameter = length_to_diameter #: Length to diameter ratio
        self.vessel_material = vessel_material # Vessel material
        self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
        self.reaction = reaction
        self.catalyst_density = catalyst_density # In kg/m^3
        self.catalyst_price = catalyst_price # In USD/kg
        self.WHSV = WHSV #: Residence time in kg/hr feed / kg catalyst
        self.P = P    
    
    def _default_vessel_type(self):
        return 'Vertical'
    
    def _setup(self):
        effluent, = self.outs
        effluent.P = self.P
        effluent.phase = 'g'
    
    def _run(self):
        feed, = self.ins
        assert feed.phase == 'g', 'feed phase must be a gas'
        effluent, = self.outs
        effluent.copy_like(feed)
        self.reaction.adiabatic_reaction(effluent)
        
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

# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import flexsolve as flx
from thermosteam import MultiStream
from biosteam import Unit
from biosteam.units.decorators import cost, design
from biosteam.units.design_tools import size_batch
import thermosteam as tmo
import biosteam as bst

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# %% Constants

_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% Biodiesel production

class BlendingTankWithSkimming(bst.MixTank):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 polar_lipid_fraction=0.1, acetone_recovered=0.9,
                 vessel_type=None, tau=None, V_wf=None, 
                 vessel_material=None, kW_per_m3=1.):
        bst.MixTank.__init__(
            self, ID, ins, outs, thermo,
            vessel_type=vessel_type, tau=tau, V_wf=V_wf, 
            vessel_material=vessel_material, kW_per_m3=kW_per_m3
        )
        self.polar_lipid_fraction = polar_lipid_fraction
        self.acetone_recovered = acetone_recovered
        
    def _run(self):
        feeds = self.ins
        polar_lipids, vegetable_oil = self.outs
        vegetable_oil.mix_from(feeds)
        polar_lipids.copy_thermal_condition(vegetable_oil)
        polar_lipids.imol['Lipid'] = vegetable_oil.imol['Lipid'] * self.polar_lipid_fraction
        polar_lipids.imol['Acetone'] = vegetable_oil.imol['Acetone'] * self.acetone_recovered
        vegetable_oil.mol -= polar_lipids.mol


class GlycerolysisReactor(bst.ContinuousReactor):
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    FFA_fraction_default = 0.9
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 P=101325, T=273.15 + 230, 
                 FFA_fraction=None, FFA_MW=282.46136, 
                 tau=0.5, V_wf=0.8, V_max=355,
                 length_to_diameter=2, kW_per_m3=0.985,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        bst.ContinuousReactor.__init__(self,
            ID, ins, outs, thermo,
            P=P, tau=tau, V_wf=V_wf, V_max=V_max,
            length_to_diameter=length_to_diameter, 
            kW_per_m3=kW_per_m3,
            vessel_material='Stainless steel 316',
            vessel_type='Vertical'
        )
        self.P = P
        self.T = T
        self.FFA_fraction = self.FFA_fraction_default if FFA_fraction is None else FFA_fraction
        self.FFA_MW = FFA_MW
        
    def _run(self):
        feed, *other, N2 = self.ins
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T
        vent.phase = 'g'
        effluent.mix_from(self.ins, energy_balance=False)
        vent.copy_flow(effluent, 'N2', remove=True)
        vent.receive_vent(effluent)
        
    def _design(self):
        bst.ContinuousReactor._design(self)
        # Heat of reaction:
        # Bonds formed:
        # H-OH (467)
        # CO-C (358)
        # Bonds broken:
        # C-OH (358)
        # Very rough estimate...
        Hnet = 467*1000 * self.FFA_fraction * self.ins[0].F_mass / self.FFA_MW
        self.heat_utilities[0](Hnet, self.T)
        
        
        
        

        

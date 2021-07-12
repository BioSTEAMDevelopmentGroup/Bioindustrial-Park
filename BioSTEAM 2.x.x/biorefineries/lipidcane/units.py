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
from thermosteam.reaction import Reaction, ParallelReaction
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
                 vessel_type=None, tau=None, V_wf=None, 
                 vessel_material=None, kW_per_m3=1.):
        bst.MixTank.__init__(
            self, ID, ins, outs, thermo,
            vessel_type=vessel_type, tau=tau, V_wf=V_wf, 
            vessel_material=vessel_material, kW_per_m3=kW_per_m3
        )
        
    def _run(self):
        feeds = self.ins
        polar_lipids, vegetable_oil = self.outs
        vegetable_oil.mix_from(feeds)
        polar_lipids.copy_thermal_condition(vegetable_oil)
        polar_lipids.imol['PL'] = vegetable_oil.imol['PL'] 
        polar_lipids.imass['Acetone'] = vegetable_oil.imass['PL'] * 2
        vegetable_oil.mol -= polar_lipids.mol


class GlycerolysisReactor(bst.ContinuousReactor):
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 P=101325, T=273.15 + 230, 
                 tau=2.0, V_wf=0.8, V_max=355,
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
        
        # End result is near 100% conversion of FFAs with selectivity
        # of 42% MAG, 47% DAG, 11% TAG:
        # [1] Erik Anderson. Superior Process Technologies. 
        # Glycerolysis for Lowering Free Fatty Acid Levels
        
        # [2] Kapil Mamtani, KavehShahbaz Mohammed M.Farid. 
        # Glycerolysis of free fatty acids: A review. 
        # https://doi.org/10.1016/j.rser.2020.110501
        
        self.glycerolysis_baseline = ParallelReaction([
          Reaction('FFA + Glycerol -> MAG + H2O', reactant='FFA',  X=1.00),
          Reaction('DAG + Glycerol -> 2MAG', reactant='DAG',  X=1.00),
          Reaction('TAG + 2Glycerol -> 3MAG', reactant='TAG',  X=1.00),
        ])
        self.glycerolysis = ParallelReaction([
          Reaction('2MAG -> DAG + Glycerol', reactant='MAG',  X=0.49),
          Reaction('3MAG -> TAG + 2Glycerol', reactant='MAG',  X=0.13),
        ])
        
        # Alternative preliminary modeling for legacy purposes:
        # self.glycerolysis = ParallelReaction([
        #   Reaction('FFA + Glycerol -> MAG + H2O', reactant='FFA',  X=0.80),
        #   Reaction('2FFA + Glycerol -> DAG + 2H2O', reactant='FFA',  X=0.15),
        #   Reaction('3FFA + Glycerol -> TAG + 3H2O', reactant='FFA',  X=0.05),
        # ])
        
    def _run(self):
        feed, *other, N2 = self.ins
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T
        vent.phase = 'g'
        effluent.mix_from(self.ins, energy_balance=False)
        self.glycerolysis_baseline.force_reaction(effluent)
        self.glycerolysis(effluent)
        vent.copy_flow(effluent, ('N2', 'Water'), remove=True)
        
    def _design(self):
        bst.ContinuousReactor._design(self)
        self.heat_utilities[0](self.Hnet, self.T)
        
        
        
        

        

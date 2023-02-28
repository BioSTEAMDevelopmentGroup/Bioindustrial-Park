# -*- coding: utf-8 -*-
"""
"""
from thermosteam.reaction import Reaction, ParallelReaction
import thermosteam as tmo
import biosteam as bst
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

__all__ = (
    'GlycerolysisReactor',
)

class GlycerolysisReactor(bst.CSTR):
    _ins_size_is_fixed = False
    _N_ins = 2
    _N_outs = 2
    T_default = 273.15 + 230
    P_default = 101325
    tau_default = 2
    
    def _setup(self):
        super()._setup()
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
        effluent.mol.remove_negatives() # The correct glycerol flow rate is taken care of in a unit specification
        vent.copy_flow(effluent, ('N2', 'Water'), remove=True)
        
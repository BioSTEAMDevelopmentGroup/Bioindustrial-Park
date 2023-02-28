# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Reactors
--------
.. autoclass:: biorefineries.cane.units.SeedTrain
.. autoclass:: biorefineries.cane.units.CoFermentation

Separations
-----------
.. autoclass:: biorefineries.cane.units.OleinCrystallizer
    

"""
import biosteam as bst
from thermosteam import PRxn, Rxn
from biorefineries.cellulosic.units import (
    SeedTrain,
    CoFermentation
)

__all__ = ('SeedTrain', 'CoFermentation', 'OleinCrystallizer')
    
class SeedTrain(SeedTrain):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, reactions=None, saccharification=False):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.saccharification = saccharification
        chemicals = self.chemicals
        self.reactions = reactions or PRxn([
    #   Reaction definition                   Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',       'Glucose',   0.9000, chemicals),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',      'Xylose',    0.8000, chemicals),
    Rxn('Glucose -> Cellmass',                'Glucose',  0.0473, chemicals),
    Rxn('Xylose -> Cellmass',                 'Xylose',  0.0421, chemicals),
        ])
        
    def _setup(self):
        super()._setup()
        self.outs[0].phase = 'g'
        
    def _run(self):
        vent, effluent= self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions.force_reaction(effluent)
        effluent.mol.remove_negatives()
        effluent.T = self.T
        vent.copy_flow(effluent, 'CO2', remove=True)
        

class CoFermentation(CoFermentation):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=36, N=None, V=3785.4118, T=305.15, P=101325,
                 Nmin=2, Nmax=36, cofermentation=None):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau, N, V, T, P, Nmin, Nmax)
        self.P = P
        chemicals = self.chemicals
        self.loss = None
        self.cofermentation = cofermentation or PRxn([
    #   Reaction definition                   Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',       'Glucose',   0.9500, chemicals),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',      'Xylose',    0.8500, chemicals),
    Rxn('Glucose -> Cellmass',                'Glucose',  0.05, chemicals),
    Rxn('Xylose -> Cellmass',                 'Xylose',  0.05, chemicals),
        ])
        
        if 'CSL' in chemicals:
            self.CSL_to_constituents = Rxn(
                'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
            )
            self.CSL_to_constituents.basis = 'mol'
        else:
            self.CSL_to_constituents = None
        if all([i in self.chemicals for i in ('FFA', 'DAG', 'TAG', 'Glycerol')]):
            self.lipid_reaction = self.oil_reaction = PRxn([
                Rxn('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
                Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
            ])
        else:
            self.lipid_reaction = self.oil_reaction = None

            
class OleinCrystallizer(bst.BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T, solid_purity=0.98, melt_purity=0.90,
                 solid_IDs=('TAG', 'FFA', 'PL'), melt_IDs=('AcTAG',),
                 order=None):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self.melt_purity = melt_purity
        self.solid_purity = solid_purity
        self.solid_IDs = solid_IDs
        self.melt_IDs = melt_IDs

    @property
    def Hnet(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        if 's' in feed.phases:
            solid = feed if feed.phase == 's' else feed['s']
            H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, solid.mol) if i.Hfus])
        else:
            H_in = 0.
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out - H_in
        
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        solid_purity = self.solid_purity
        melt_purity = self.melt_purity
        feed = self.ins[0]
        solid_IDs = self.solid_IDs
        melt_IDs = self.melt_IDs
        solid_flows = feed.imass[solid_IDs]
        melt_flows = feed.imass[melt_IDs]
        net_solid_flow = solid_flows.sum()
        net_melt_flow = melt_flows.sum()
        total = net_solid_flow + net_melt_flow
        minimum_melt_purity = net_melt_flow / total
        minimum_solid_purity = net_solid_flow / total
        outlet.empty()
        if solid_purity < minimum_solid_purity:
            outlet.imol['s'] = feed.mol
        elif melt_purity < minimum_melt_purity:
            outlet.imol['l'] = feed.mol
        else: # Lever rule
            solid_purity = (1. - solid_purity)
            melt_fraction = (minimum_melt_purity - solid_purity) / (melt_purity - solid_purity)
            melt = melt_fraction * total
            pure_melt = melt * melt_purity
            impurity_melt = melt - pure_melt
            outlet.imass['l', melt_IDs] = pure_melt * melt_flows / melt_flows.sum()
            outlet.imass['l', solid_IDs] = impurity_melt * solid_flows / solid_flows.sum()
            outlet.imol['s'] = feed.mol - outlet.imol['l']
        outlet.T = self.T
        
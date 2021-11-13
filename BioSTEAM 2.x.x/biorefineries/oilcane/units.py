# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 03:19:53 2021

@author: yrc2
"""
from biosteam import Unit, BatchBioreactor
from thermosteam import PRxn, Rxn
from biorefineries.cornstover.units import (
    SeedTrain,
    CoFermentation
)

__all__ = ('SeedTrain', 'CoFermentation')

class SeedTrain(SeedTrain):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, reactions=None, saccharification=False):
        Unit.__init__(self, ID, ins, outs, thermo)
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
        self.outs[0].phase = 'g'
        
    def _run(self):
        vent, effluent= self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions.force_reaction(effluent)
        effluent.mol[effluent.mol < 0.] = 0.
        effluent.T = self.T
        vent.copy_flow(effluent, 'CO2', remove=True)
        

class CoFermentation(CoFermentation):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=36, N=None, V=3785.4118, T=305.15, P=101325,
                 Nmin=2, Nmax=36, cofermentation=None):
        BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau, N, V, T, P, Nmin, Nmax)
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
    
        self.CSL_to_constituents = Rxn(
            'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
        )
        self.CSL_to_constituents.basis = 'mol'
        
        if all([i in self.chemicals for i in ('FFA', 'DAG', 'TAG', 'Glycerol')]):
            self.lipid_reaction = self.oil_reaction = PRxn([
                Rxn('TAG + 3Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
                Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
            ])
        else:
            self.lipid_reaction = self.oil_reaction = None
            
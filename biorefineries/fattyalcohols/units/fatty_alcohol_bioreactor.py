# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 01:25:56 2020

@author: yoelr
"""
from biosteam import BatchBioreactor

__all__ = ('FattyAlcoholBioreactor',)

class FattyAlcoholBioreactor(BatchBioreactor):
    _N_ins = 1
    _N_outs = 2
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 fermentation_reaction, tau,
                 prefermentaion_reaction=None,
                 postfermentation_reaction=None,
                 N=None, V=None, T, P=101325, Nmin=2, Nmax=32):
        BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau=tau,
                                 N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self.prefermentaion_reaction = prefermentaion_reaction
        self.fermentation_reaction = fermentation_reaction
        self.postfermentation_reaction = postfermentation_reaction
        
    @property
    def efficiency(self):
        return self.fermentation_reaction.X
    @efficiency.setter
    def efficiency(self, efficiency):
        self.fermentation_reaction.X = efficiency

    @property
    def Hnet(self):
        return self.fermentation_reaction.dH * self.ins[0].imol[self.fermentation_reaction.reactant]

    def _run(self):
        vent, effluent = self.outs
        effluent.T = self.T
        effluent.P = self.P
        effluent.copy_flow(self.ins[0])
        if self.prefermentaion_reaction: self.prefermentaion_reaction(effluent)
        self.fermentation_reaction(effluent)
        if self.postfermentation_reaction: self.postfermentation_reaction(effluent)
        vent.copy_flow(effluent, 'CO2', remove=True)
        vent.copy_thermal_condition(effluent)
        
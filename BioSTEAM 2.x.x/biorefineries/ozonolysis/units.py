# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 08:17:38 2021

@author: yrc2
"""

import biosteam as bst
import thermosteam as tmo

class OzonolysisReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    @property
    def effluent(self):
        return self.outs[0]

    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=17, N=None, V=None, T=373.15, P=101325,
                 Nmin=2, Nmax=36):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                   tau = tau , N = N, V = V, T = T, 
                                   P =P ,Nmin =Nmin, Nmax = Nmax)
        c = self.chemicals
        self.reactant_conversion = 0.808
        #self.nonanal_selectivity = 0.094
        self.selectivity_epoxide = 0.071 * c.Oleic_acid.MW / c.Epoxy_stearic_acid.MW 
        #self.oxononanoic_acid_selectivity = 0.029
        #self.nonanoic_acid_selectivity = 0.277
        self.selectivity_azelaic_acid = 0.442 * c.Oleic_acid.MW / c.Azelaic_acid.MW 
        
    def _setup(self):
        # selectivity_azelaic_acid = selectivity_nonanoic_acid = X2 * X1
        # selectivity_nonanal = selectivity_oxonanoic_acid = X1 * (1 - X2)
        # selectivity_epoxide = (1 - X1)
        X1 = 1 - self.selectivity_epoxide
        X2 = self.selectivity_azelaic_acid / X1
        self.reactions = tmo.SeriesReaction([
            # Assumption, every conversion is 0.947 and overall conversion is (0.947^3)
            tmo.Rxn('Oleic_acid + H2O2 -> Epoxy_stearic_acid+ Water ', 'Oleic_acid', X=self.reactant_conversion),
            tmo.Rxn('Epoxy_stearic_acid + H2O2 -> Nonanal + Oxononanoic_acid + H2O', 'Epoxy_stearic_acid', X = X1),
            tmo.Rxn('Nonanal + Oxononanoic_acid + 2H2O2 -> Azelaic_acid + Nonanoic_acid+ 2H2O', 'Nonanal', X = X2)]
        )
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(feed)
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        
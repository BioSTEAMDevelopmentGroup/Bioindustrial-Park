# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from thermosteam import PRxn, Rxn

__all__ = ('OleinCrystallizer', 'Fermentation')

class OleinCrystallizer(bst.BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T, crystal_TAG_purity=0.95, melt_AcTAG_purity=0.90,
                 order=None):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self.melt_AcTAG_purity = melt_AcTAG_purity
        self.crystal_TAG_purity = crystal_TAG_purity

    @property
    def Hnet(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        if 's' in feed.phases:
            H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
        else:
            H_in = 0.
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out - H_in
        
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        crystal_TAG_purity = self.crystal_TAG_purity
        melt_AcTAG_purity = self.melt_AcTAG_purity
        feed = self.ins[0]
        TAG, AcTAG = feed.imass['TAG', 'AcTAG']
        total = TAG + AcTAG
        minimum_melt_purity = AcTAG / total
        minimum_crystal_purity = TAG / total
        outlet.empty()
        if crystal_TAG_purity < minimum_crystal_purity:
            outlet.imol['s'] = feed.mol
        elif melt_AcTAG_purity < minimum_melt_purity:
            outlet.imol['l'] = feed.mol
        else: # Lever rule
            crystal_AcTAG_purity = (1. - crystal_TAG_purity)
            melt_fraction = (minimum_melt_purity - crystal_AcTAG_purity) / (melt_AcTAG_purity - crystal_AcTAG_purity)
            melt = melt_fraction * total
            AcTAG_melt = melt * melt_AcTAG_purity
            TAG_melt = melt - AcTAG
            outlet.imass['l', ('AcTAG', 'TAG')] = [AcTAG_melt, TAG_melt]
            outlet.imol['s'] = feed.mol - outlet.imol['l']
        outlet.T = self.T
        
        
class Fermentation(bst.BatchBioreactor):
    line = 'Fermentation'
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 tau,  N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                 tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax)
        self._load_components()
        chemicals = self.chemicals
        self.hydrolysis_reaction = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.00, chemicals)
        self.fermentation_reaction = PRxn([
        Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156, chemicals),
        Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165, chemicals),
        Rxn('Glucose -> Cells', 'Glucose', 0.10, chemicals, basis='wt').copy(basis='mol'),
            ])
        self.CSL_to_constituents = Rxn(
            'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
        )

    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.CSL_to_constituents(effluent)
        self.hydrolysis_reaction.force_reaction(effluent)
        effluent.mol[effluent.mol < 0.] = 0.
        self.fermentation_reaction(effluent)
        vent.empty()
        vent.receive_vent(effluent)

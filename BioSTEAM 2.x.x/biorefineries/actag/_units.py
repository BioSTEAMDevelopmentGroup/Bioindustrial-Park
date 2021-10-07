# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst

__all__ = ('OleinCrystallizer',)

class OleinCrystallizer(bst.BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T,
                 split={'TriOlein': 0.90, 'AcetylDiOlein': 0.05},
                 order=None):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self._isplit = self.thermo.chemicals.isplit(split, order)

    @property
    def isplit(self):
        """[ChemicalIndexer] Componentwise split of feed to 0th outlet stream."""
        return self._isplit
    @property
    def split(self):
        """[Array] Componentwise split of feed to 0th outlet stream."""
        return self._isplit._data
    @split.setter
    def split(self, values):
        split = self.split
        if split is not values:
            split[:] = values
     
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
        self.ins[0].split_to(outlet['s'], outlet['l'], self.split, energy_balance=False)
        outlet.T = self.T
        
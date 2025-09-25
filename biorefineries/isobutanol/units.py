# -*- coding: utf-8 -*-
"""
Created on Wed Sep 24 13:08:23 2025

@author: saran
"""
import numpy as np
import biosteam as bst
import thermosteam as tmo

cost = bst.decorators.cost
Splitter = bst.Splitter

__all__ = ('SimultaneousSaccharificationFermentationEthanolIsobutanol',
           'SSFEtOHIBO', 'MolecularSieve')

class SimultaneousSaccharificationFermentationEthanolIsobutanol(bst.BatchBioreactor):
    """
    Create a SimultaneousSaccharificationFermentation unit operation that 
    models the simultaneous saccharification and fermentation in the conventional
    dry-grind ethanol process.
    
    Parameters
    ----------
    ins : streams
        Inlet fluids.
    outs : stream
        Outlet fluid.
    yield_: float
        Yield of glucose to ethanol as a fraction of the theoretical yield.
    
    Notes
    -----
    This unit operation doesn't actually model the saccharification process.
    The reactor is modeled by the stoichiometric conversion of glucose to
    ethanol by mol:
        
    .. math:: 
        Glucose -> 2Ethanol + 2CO_2
    
    Yeast is assumed to be produced from any remaining glucose:
        Glucose -> Yeast
    
    A compound with name 'Yeast' must be present. Note that only glucose is 
    taken into account for conversion. Cleaning and unloading time,
    `tau_0`, fraction of working volume, `V_wf`, and number of reactors,
    `N_reactors`, are attributes that can be changed. Cost of a reactor
    is based on the NREL batch fermentation tank cost assuming volumetric
    scaling with a 6/10th exponent [1]_. 
    
    References
    ----------
    .. [1] D. Humbird, R. Davis, L. Tao, C. Kinchin, D. Hsu, and A. Aden
        National. Renewable Energy Laboratory Golden, Colorado. P. Schoen,
        J. Lukas, B. Olthof, M. Worley, D. Sexton, and D. Dudgeon. Harris Group
        Inc. Seattle, Washington and Atlanta, Georgia. Process Design and Economics
        for Biochemical Conversion of Lignocellulosic Biomass to Ethanol Dilute-Acid
        Pretreatment and Enzymatic Hydrolysis of Corn Stover. May 2011. Technical
        Report NREL/TP-5100-47764
    
    
    """
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 tau=60.,  N=None, V=None, T=305.15, P=101325., Nmin=2, Nmax=36,
                 yield_ethanol=0.80, yield_isobutanol=0.15, V_wf=0.83):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
            tau=tau, N=N, V=V, T=T, P=P, Nmin=Nmin, Nmax=Nmax
        )
        self.fermentation_reactions = tmo.ParallelReaction([tmo.Rxn('Glucose -> 2Ethanol + 2CO2',  'Glucose', yield_ethanol),
                                             tmo.Rxn('Glucose -> Isobutanol + 2CO2',  'Glucose', yield_isobutanol),],
                                             )
        self.growth = tmo.Rxn('Glucose -> Yeast',  'Glucose', 1.0)
        self.V_wf = V_wf
    
    def _run(self):
        vent, effluent = self.outs
        effluent.mix_from(self.ins)
        self.fermentation_reactions(effluent)
        self.growth(effluent)
        effluent.imass['Yeast'] += effluent.imass['NH3']
        effluent.imol['NH3'] = 0.
        vent.empty()
        vent.receive_vent(effluent)
        
    @property
    def Hnet(self):
        # X = self.fermentation_reactions.X
        fermentation_reactions = self.fermentation_reactions
        glucose = self.ins[0].imol['Glucose']
        # return self.reaction.dH * glucose + self.growth.dH * (1 - X) + self.H_out - self.H_in
        # return np.sum([i.dH for i in fermentation_reactions]) * glucose + np.sum([self.growth.dH * (1 - i.X) + self.H_out - self.H_in for i in fermentation_reactions])
        return np.sum(fermentation_reactions[0].dH) * glucose + np.sum([self.growth.dH * (1 - i.X) + self.H_out - self.H_in for i in fermentation_reactions]) #!!! reactions[0]

SSFEtOHIBO = SimultaneousSaccharificationFermentationEthanolIsobutanol

#%%
@cost('Flow rate', 'Column', kW=151, BM=1.8,
      cost=2601000, CE=521.9, S=22687, n=0.6)
class MolecularSieve(Splitter):
    """
    Create an isobutanol/water molecular sieve for bioisobutanol plants.
    The molecular sieve is modeled as a component wise separator. Costing
    is based on scaling by the 6/10ths rule from an NREL TEA report [1]_.
    
    Parameters
    ----------
    ins : 
        * [0] Feed (gas)
    outs : 
        * [0] Split stream (gas)
        * [1] Remainder stream (gas)
    split : array_like
            Componentwise split to the 0th output stream

    References
    ----------
    .. [1] Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and
        Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L.
        Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy
        Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof,
        M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle,
        Washington and Atlanta, Georgia)
    
    """
    _units = {'Flow rate': 'kg/hr'}
    def _init(self, split, order=None, P=None, approx_duty=True):
        Splitter._init(self, order=order, split=split)
        self.P = None
        self.approx_duty = approx_duty
        
    def _run(self):
        Splitter._run(self)
        P = self.P
        if P is None: P = self.ins[0].P
        for i in self.outs: i.P = P

    def _design(self):
        self.design_results['Flow rate'] = flow = self._outs[1].F_mass
        if self.approx_duty:
            T = self.outs[0].T
            self.add_heat_utility(1429.65 * flow, T)
            self.add_heat_utility(-55.51 * flow, T)
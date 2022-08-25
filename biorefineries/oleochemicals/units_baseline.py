# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 23:55:59 2022

@author: Lavanya_Kudli
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream
import chemicals_baseline
from biosteam.units.design_tools import size_batch
from biosteam import decorators 
import flexsolve as flx
from math import ceil
import numpy as np
from biosteam import Unit

#Need to add a vaccuum evaportator
class DihydroxylationReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=17, N=None, V=None, T= 273.15, P=None,
                 Nmin=2, Nmax=36):
        
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                   tau = tau , N = N, V = V, T = T, 
                                   P = P ,Nmin = Nmin , Nmax = Nmax)
    def _setup(self):  
                         
            self.reactions = tmo.SeriesReaction([
                tmo.Rxn('Methyl_oleate + H2O2   -> MDHSA ', 'Methyl_oleate', X = 0.9)])
            
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        
class OxidativeCleavageReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 2
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
   #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, '')) 
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=3.5, N=None, V=None, T= 60+273.15, P=101325,
                 Nmin=2, Nmax=36):
        
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                   tau = tau , N = N, V = V, T = T, 
                                   P = P ,Nmin = Nmin , Nmax = Nmax)
    def _setup(self):           
        #oxidative_cleavage_conversion
        X1 = 0.8
        # decarboxylations conversion
        X2 = 0.2
        # X3 = 0.2
        Product_formation = PRxn([Rxn('MDHSA + 1.5 Oxygen  -> Pelargonic_acid + Methyl_azelate','MDHSA', X = X1),
                                  Rxn('MDHSA  ->  Caprylic_acid + Methyl_azelate','MDHSA', X = X2)])                        
        #TODO.xxx check again possible decarboxylation https://doi.org/10.1016/j.renene.2018.01.107
        #TODO.xxx check again Organic Reactions in Strong Alkalis. Part V.l Alkali Fusion of Epoxides and Ethers 

        # Side_reaction = SRxn([Rxn('Pelargonic_acid   -> Caprylic_acid + carbon_dioxide ', 'Pelargonic_acid', X= X2),
                               # Rxn('Methyl_azelate -> suberic_acid + carbon_dioxide', 'Methyl_azelate', X = X3),
                             # ])
        oxidative_cleavage_rxnsys = RxnSys(Product_formation)
        self.reactions = oxidative_cleavage_rxnsys
            
    def _run(self):
        feed = self.ins[0]
        vent, effluent = self.outs    
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        vent.phase = 'g'
        vent.copy_flow(effluent, ('Nitrogen', 'Oxygen'), remove=True)
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P

              
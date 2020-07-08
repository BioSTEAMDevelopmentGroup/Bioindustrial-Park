#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 11:40:48 2020

@author: yalinli_cabbi
"""

'''
TODO:
    NEED CSL AND DAP OR NOT?
    Look for #???

'''



# %% Setup

import numpy as np
import biosteam as bst
from scipy.integrate import odeint
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.designtools import size_batch
from orgacids.EB import EB_simulate

Rxn = bst.reaction.Reaction
ParallelRxn = bst.reaction.ParallelReaction
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3


# %%

# Seed train (5 stages)
@cost('Flow rate', 'Pumps',
      S=43149, CE=522, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=522, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=522, n=0.5, BM=1.5)
class SeedTrain(Unit):
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5
    
    # Number of parallel seed trains
    N_trains = 2
    
    # Cycle time for each batch (hr)
    tau_batch = 24
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)

    def _run(self):
        feed = self.ins[0]
        vent, effluent= self.outs
        effluent.copyflow(feed)
        for i in range(0, 5):
            simulated_flow = EB_simulate(effluent, self.tau_batch)
            effluent.mass[effluent.index('Z_mobilis')] = simulated_flow[0]
            effluent.mass[effluent.index('LacticAcid')] = simulated_flow[1]
            effluent.mass[effluent.index('Glucose')] = simulated_flow[2]
            effluent.mass[effluent.index('Xylose')] = simulated_flow[3]
            # Assume all lost glucose and xylose changed to CO2
            effluent.mass[effluent.index('CO2')] += \
                (simulated_flow[2] + simulated_flow[3]) - \
                (simulated_flow[0] + simulated_flow[1])
        # Assume all DAP and CSL have been consumed
        effluent.mass[effluent.index('DAP')] = 0
        effluent.mass[effluent.index('CSL')] = 0
        effluent.T = self.T
        vent.phase = 'g'
        vent.copyflow(effluent, ('CO2', 'NH3', 'O2'), remove=True)

    def _design(self): 
        maxvol = self.outs[1].volnet*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self._Design
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.massnet for i in self.outs])
        self._heat_utilities[0](self._Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self._Design
        C = self._Cost
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self._power_utility(kW)


# %% Trial

from orgacids.species import species
from orgacids.EB import EB_simulate

bst.Stream.species = species
sugar = bst.Stream('sugar',
                   H2O=1000,
                   Z_mobilis=11,
                   LacticAcid = 53,
                   Glucose=182,
                   Xylose=205
                   )

outs1 = bst.Stream('outs1')
outs2 = bst.Stream('outs2')

R302 = SeedTrain('R302', ins=sugar, outs=('outs1', 'outs2'))
R302.simulate()
R302.show()


# %%

# Saccharification and co-fermentation reactor
@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr, ub=False,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_recirculation_pumps')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900, ub=False,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500, ub=False,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000, ub=False,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr, ub=False,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, ub=False,
      CE=522, n=0.7, BM=2.0, N='N_tanks')
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 3
    _N_outs = 3
    _N_heat_utilities = 2
    
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
    #: Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    #: Saccharification time (hr)
    tau_saccharification = 60
    
    #: Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    #: Number of reactors
    N_reactors = 12
    
    #: Number of continuous saccharification tanks
    N_tanks = 8
    
    #: Number of transfer pumps
    N_transfer_pumps = 5
    
    #: Number of recirculation pumps
    N_recirculation_pumps = 5
    
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    # Split to outs[2]
    saccharified_slurry_split = 0.1

    
    def __init__(self, ID='', ins=None, outs=(), P=101325):
        Unit.__init__(self, ID, ins, outs)

        self.P = P
        
        # Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        # Kept the same as Humbird et al.
        self.saccharification = ParallelRxn([
            #   Reaction definition                   Reactant    Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000),
            Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1.0000)])
        
        self.saccharified_stream = bst.Stream(None)
  
    def _run(self):
        feed, CSL, DAP = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        ss.copyflow(feed)
        self.saccharification(ss.mol)
        # Sidedraw after saccharification but before cofermentation
        sidedraw.mol[:] = ss.mol*self.saccharified_slurry_split
        # Now effluent is after saccharification but before cofermentation
        effluent.mol[:] = ss.mol - sidedraw.mol + CSL.mol + DAP.mol        
        loss = 0
        simulated_flow = EB_simulate(effluent, self.tau_cofermentation)
        effluent.mass[effluent.index('Z_mobilis')] = simulated_flow[0]
        effluent.mass[effluent.index('LacticAcid')] = simulated_flow[1]
        effluent.mass[effluent.index('Glucose')] = simulated_flow[2]
        effluent.mass[effluent.index('Xylose')] = simulated_flow[3]
        loss += (simulated_flow[2] + simulated_flow[3]) - \
            (simulated_flow[0] + simulated_flow[1])
        # Assume all lost glucose and xylose changed to CO2
        effluent.mass[effluent.index('CO2')] = feed.mass[effluent.index('CO2')] + loss
        # Assume all DAP and CSL have been consumed
        effluent.mass[effluent.index('DAP')] = 0
        effluent.mass[effluent.index('CSL')] = 0
        vent.copyflow(effluent, ('CO2', 'NH3', 'O2'), remove=True)
        vent.recieve_vent(effluent)
    
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.volnet
        Design = self._Design
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf/self.N_tanks
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        hu_cooling, hu_fermentation = self._heat_utilities
        hu_cooling(self.saccharified_stream.H_at(T=self.T_fermentation)
                   - self.saccharified_stream.H_at(T=self.T_saccharification),
                   self.T_fermentation)
        ei = effluent.index('LacticAcid')
        LA = (sum([i._mol[ei] for i in self.outs])
                   - sum([i._mol[ei] for i in self.ins]))
        duty = LA*-5568 #??? Where does this come from?
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty


# %% 
        
from orgacids.system import M302, CSL2, DAP2
from orgacids.species import species

bst.Stream.species = species
sugar = bst.Stream('sugar',
                   H2O=1000,
                   Z_mobilis=11,
                   LacticAcid = 53,
                   Glucose=182,
                   Xylose=205
                   )

outs1 = bst.Stream('outs1')
outs2 = bst.Stream('outs2')
outs3 = bst.Stream('outs3')

R301 = SaccharificationAndCoFermentation('R301', ins=(sugar, CSL2, DAP2), 
                                         outs=('outs1', 'outs2', 'outs3')
                                         )
R301.simulate()
R301.show()

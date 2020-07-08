#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:19:18 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

All units are explicitly defined here for transparency and easy reference

@author: yalinli_cabbi
"""

'''
TODO:
   ADD SEPARATION UNITS (HOPEFULLY ONLY USING THE COST DECORATOR)
'''


# %% Setup

from flexsolve import aitken_secant
from biosteam import Unit, MixedStream, Stream
from biosteam.units.decorators import cost
from biosteam.units.designtools import size_batch
import biosteam as bst

Rxn = bst.reaction.Reaction
ParallelRxn = bst.reaction.ParallelReaction
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# TEMPORARY
# Fitted kinetic parameters, g for glucose and x for xylose
n = 3 # [-]
P_max = 100 # [kg/m3]
Y_XS_g = 0.08 # [kg/kg]
Y_XS_x = 0.11 # [kg/kg]
Y_PS_g = 0.33 # [kg/kg]
Y_PS_x = 0.66 # [kg/kg]
mu_max_g = 0.21 # [kg/kg]
mu_max_x = 0.087 # [kg/kg]
K_S_g = 30 # [kg/m3]
K_S_x = 7.39 # [kg/m3]


# %% Feedstock handling

# Feedstock handling system as a whole, capital and operating costs considered in feedstock cost
# CE = 522 is year 2009
@cost(basis='Flow rate', ID='System', units='kg/hr',
      S=94697, ub=False, CE=522, cost=3329690, n=0.6, kW=783, BM=1.7)
class FeedstockHandling(bst.units.Static): pass


# %% Pretreatment

# Sulfuric acid addition tank
# CE = 551 is year 2010
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=551, cost=6210, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=3720, ub=False, CE=522, cost=8000, n=0.8, kW=1, BM=2.3)
class SulfuricAcidAdditionTank(bst.units.Static): pass

# Sulfuric acid in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=136260, ub=False, CE=522, cost=6000, n=0.6, kW=0, BM=1)
class SulfuricAcidMixer(bst.units.Mixer): pass

# Blowdown discharge pump
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=551, cost=25365, n=0.8, kW=93.21, BM=2.3)
class BlowdownDischargePump(bst.units.Static): pass

# Pretreatment flash tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=252891, ub=False, CE=522, cost=511000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=252891, ub=False, CE=522, cost=30000, n=0.8, kW=55.9275, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=252891, ub=False, CE=522, cost=90000, n=0.5, kW=170, BM=1.5)
class PretreatmentFlash(bst.units.Flash): pass

# Oligomer conversion tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=522, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=551, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=522, cost=90000, n=0.5, kW=170, BM=1.5)
class OligomerConversionTank(bst.units.Static): pass

# Ammonia addition tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=410369, ub=False, CE=522, cost=236000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=410369, ub=False, CE=522, cost=21900, n=0.5, kW=7.457, BM=1.5)
class AmmoniaAdditionTank(bst.units.Static): pass

# Ammonia in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=157478, ub=False, CE=522, cost=5000, n=0.5, kW=0, BM=1)
class AmmoniaMixer(bst.units.Mixer): pass

# Waste vapor condenser
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-2, ub=False, CE=522, cost=34000, n=0.7, kW=0, BM=2.2)
class WasteVaporCondenser(bst.units.HXutility): pass

# Steam mixer
class SteamMixer(Unit):
    """
    **ins**
    
        [0] Feed
        
        [1] Steam
    
    **outs**
    
        [0] Mixed
    
    """
    _N_outs = 1
    _N_ins = 2
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), *, P):
        super().__init__(ID, ins, outs)
        self.P = P
    
    @staticmethod
    def _P_at_flow(mol_water, P, mol_array, index_water, mixed, ins):
        mol_array[index_water] = mol_water
        Stream.sum(mixed, ins)
        P_new = mixed.bubble_P()[0]
        return P - P_new
    
    def _run(self):
        ins = self._ins
        steam = ins[1]
        steam_mol = steam.molnet
        mixed = self.outs[0]
        index_water = mixed.index('7732-18-5')
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam.mol, index_water, mixed, ins))
        mixed.P = self.P
        hu = self._heat_utilities[0]
        hu.ID = 'Low pressure steam'
        hu.flow = steam_mol
        hu.cost = steam_mol*bst.HeatUtility.heating_agents['Low pressure steam'].price_kmol
    
    @property
    def installation_cost(self): return 0
    @property
    def purchase_cost(self): return 0
    def _design(self): pass
    def _cost(self): pass
        
# Pretreatment reactor
@cost('Dry flow rate', units='kg/hr', S=83333, CE=522, ub=False,
      cost=19812400, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.Flash._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._mixedstream = MixedStream(None)
        self.reactions = ParallelRxn([
    #            Reaction definition                 Reactant    Conversion
    Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.0990),
    Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.0030),
    Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0030),
    Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0030),
    Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030),
    Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.0030),
    Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030),
    Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.0000),
    Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000),
    Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.0024),
    Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0050),
    Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000),
    Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0024),
    Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050),
    Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000),
    Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0050)])
        vapor, liquid = self.outs
        vapor.phase = 'g'
    
    def _run(self):
        ms = self._mixedstream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copylike(feed)
        self.reactions(liquid.mol) 
        ms.copylike(liquid)
        ms.VLE(T=130+273.15, Q=(liquid.Hf-feed.Hf))
        vapor.mol[:] = ms.vapor_mol
        liquid.mol[:] = ms.liquid_mol
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P


# %% Saccharification and fermentation

# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=40414, ub=False, CE=522, cost=439000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=43149, ub=False, CE=522, cost=8200, n=0.8, kW=7.457, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=40414, ub=False, CE=522, cost=31800, n=0.5, kW=11.1855, BM=1.5)
class SeedHoldTank(bst.units.Static): pass

# Hydrolysate cooler
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-8, ub=False, CE=522, cost=85000, n=0.7, kW=0, BM=2.2)
class HydrolysateCooler(bst.units.HXutility): pass

# Enzyme hydrolysate mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=380000, ub=False, CE=522, cost=109000, n=0.5, kW=74.57, BM=1.7)
class EnzymeHydrolysateMixer(bst.units.Mixer): pass

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
    
    #: Number of parallel seed trains
    N_trains = 2
    
    #: Cycle time for each batch (hr)
    tau_batch = 24
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: wt % media (e.g. corn steep liquor) in each stage 
    # media_loading = 0.50
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def __init__(self, ID='', ins=None, outs=(), iskinetic=False):
        Unit.__init__(self, ID, ins, outs)
        self.iskinetic = iskinetic
        
        self.reactions = ParallelRxn([
    #   Reaction definition                             Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',                      'Glucose',   0.9000),
    Rxn('Glucose + 8.463 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                        'Glucose',   0.0400),
    Rxn('3 Xylose -> 5 LacticAcid',                     'Xylose',    0.8000),
    Rxn('Xylose + 7.052 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                        'Xylose',    0.0400)
                                       ])

    def _run(self):
        feed = self.ins[0]
        vent, effluent= self.outs
        effluent.copyflow(feed)
        self.reactions(effluent.mol)
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
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        self.saccharification = ParallelRxn([
    #   Reaction definition                   Reactant    Conversion
    Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120),
    Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000),
    Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1.0000)])
    
        self.cofermentation = ParallelRxn([
    #   Reaction definition                                         Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',                                  'Glucose',   0.9000),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.0200),
    Rxn('3 Xylose -> 5 LacticAcid',                                 'Xylose',    0.8500),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',    'Xylose',    0.0190)
                                            ])
        # CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.
        # Lactic acid is changed to Glucose to avoid generating extra products
        self.CSL2constituents = Rxn(
        'CSL -> 0.5 H2O + 0.125 Glucose + 0.25 Protein', 'CSL',    1.0000)
    
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
        sidedraw.mol[:] = ss.mol*self.saccharified_slurry_split
        effluent.mol[:] = ss.mol - sidedraw.mol + CSL.mol + DAP.mol
        self.cofermentation(effluent.mol)
        self.CSL2constituents(effluent.mass)
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
        ei = effluent.index('Ethanol')
        ethanol = (sum([i._mol[ei] for i in self.outs])
                   - sum([i._mol[ei] for i in self.ins]))
        duty = ethanol*-5568
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty


# %% Organic acid separation

# Placeholder copying the AmmoniaAdditionTank unit from cornstover biorefinery
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=410369, ub=False, CE=522, cost=236000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=410369, ub=False, CE=522, cost=21900, n=0.5, kW=7.457, BM=1.5)
class SeparationChemicalsAdditionTank(bst.units.Static): pass

# Placeholder copying the BeerTank unit from cornstover biorefinery
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=425878, ub=False, CE=522, cost=636000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=425878, ub=False, CE=522, cost=26800, n=0.8, kW=93.2125, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=425878, ub=False, CE=522, cost=68300, n=0.5, kW=14.914, BM=1.5)
class SeparationSplitter(bst.units.Splitter): pass

# Placeholder copying the heat exchangers of the SaccharificationAndCoFermentation
# unit from cornstover biorefinery
@cost(basis='Duty', ID='Heat exchangers', units='Gcal/hr',
      S=5*_Gcal2kJ, ub=False, CE=522, cost=23900, n=0.7, BM=2.2)
class SeparationHX(bst.units.HXutility): pass


# %% Lignin separation

@cost('Flow rate', 'Flitrate tank agitator', ub=False,
      cost=26e3, CE=551, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump', ub=False,
      cost=13040, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank', ub=False,
      cost=103e3, S=31815, CE=551, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57, ub=False,
      cost= 18173, S=31815, CE=551, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531', ub=False,
      cost=174800, CE=551, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Mafifold flush pump', kW=74.57, ub=False,
      cost=17057, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank', ub=False,
      cost=1520, CE=551, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Lignin wet cake screw',  kW=15*_hp2kW, ub=False,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Lignin wet cake conveyor', kW=10*_hp2kW, ub=False,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter', ub=False,
      cost=3294700, CE=551, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank', ub=False,
      cost=8e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW, ub=False,
      cost=29154, CE=551, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank', ub=False,
      cost=17e3, CE=551, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter', ub=False,
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)', ub=False,
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class PressureFilter(bst.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _design(self):
        self._Design['Flow rate'] = self.outs[0].massnet


# %% Wastewater treatment

# The total cost of wastewater treatment is combined into this placeholder
@cost('Flow rate', 'Wastewater system', units='kg/hr', CE=551, ub=False,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WastewaterSystemCost(bst.Static): pass

class AnaerobicDigestion(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between wastewater and sludge
        
    **ins**
    
        [0] Wastewater
        
        [1] Cool well water
        
    **outs**
    
        [0] Biogas
        
        [1] Wastewater
        
        [2] Sludge
        
        [3] Hot well water
    
    """
    purchase_cost = installation_cost = 0
    _N_ins = 2
    _N_outs = 4
    def __init__(self, ID='', ins=None, outs=(), *, reactions, sludge_split):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
        self.sludge_split = sludge_split
        self.mixed_stream = bst.MixedStream()
    
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link = cool_water
        biogas.T = waste.T = sludge.T = T = 35+273.15
        hot_water.T = feed.T - 5
        cool_water.mol[:] *= (feed.H - feed.H_at(T=T))/(hot_water.H - cool_water.H)
        sludge.copyflow(feed)
        self.reactions(sludge.mol)
        self.mixed_stream.copyflow(sludge)
        self.mixed_stream.VLE(P=101325, Q=0)
        biogas.mol[:] = self.mixed_stream.vapor_mol
        liquid_mol = self.mixed_stream.liquid_mol
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.recieve_vent(waste)     
    
class AerobicDigestion(bst.Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **reactions:** [ReactionSet] Anaerobic digestion reactions.
        
        **sludge_split:** [Array] Split between wastewater and sludge
        
    **ins**
    
        [0] Wastewater
        
        [1] Air
        
        [2] Caustic
        
    **outs**
    
        [0] Vent
        
        [1] Treated wastewater

    """    
    _N_ins = 3
    _N_outs = 2
    purchase_cost = installation_cost = 0
    evaporation = 4/355
    
    def __init__(self, ID='', ins=None, outs=(), *, reactions):
        Unit.__init__(self, ID, ins, outs)
        self.reactions = reactions
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copylike(waste)
        water.mol[:] += air.mol
        water.mol[:] += caustic.mol
        self.reactions(water.mol)
        vent.copyflow(water, ('CO2', 'O2', 'N2'))
        wi = vent.index('Water')
        water_mol = water.mol[wi]
        vent.mol[wi] = water_mol * self.evaporation
        water.mol[:] -= vent.mol


# %% Facilities

# Sulfuric acid storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=551, cost=96000, n=0.7, kW=0, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1981, ub=False, CE=522, cost=7493, n=0.8, kW=0.37285, BM=2.3)
class SulfuricAcidStorage(bst.units.StorageTank): pass

# Ammonia storage tank, no pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1171, ub=False, CE=551, cost=196000, n=0.7, kW=0, BM=2)
class AmmoniaStorageTank(bst.units.StorageTank): pass

# DAP storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=163, ub=False, CE=522, cost=102000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=163, ub=False, CE=522, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=163, ub=False, CE=522, cost=9800, n=0.5, kW=4.10135, BM=1.5)
@cost(basis='Flow rate', ID='Bag', units='kg/hr',
      S=163, ub=False, CE=522, cost=30000, n=0.6, kW=0, BM=1.7)
class DAPStorageTank(bst.units.Static): pass

# CSL storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1393, ub=False, CE=522, cost=70000, n=0.7, kW=0, BM=2.6)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1393, ub=False, CE=551, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=1393, ub=False, CE=522, cost=21200, n=0.5, kW=7.457, BM=1.5)
class CSLStorageTank(bst.units.Static): pass

# Placeholder copying the CSLStorageTank unit from cornstover biorefinery
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1393, ub=False, CE=522, cost=70000, n=0.7, kW=0, BM=2.6)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1393, ub=False, CE=551, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=1393, ub=False, CE=522, cost=21200, n=0.5, kW=7.457, BM=1.5)
class SeparationChemicalsStorageTank(bst.units.Static): pass

# Fire water tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=8343, ub=False, CE=522, cost=803000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=8343, ub=False, CE=522, cost=15000, n=0.8, kW=93.2125, BM=1.7)
class FireWaterTank(bst.units.Static): pass

@cost('Flow rate', units='kg/hr', ub=False,
      S=63, cost=421e3, CE=522, BM=1.8, n=0.6)
class CIPpackage(bst.Facility):
    line = 'Clean in place package'
    _N_ins = 1
    _N_outs = 1
    
    

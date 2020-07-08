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
Considerations for different organic acids
    Fermentation kinetic parameters
    Fermentation conditions
       Temperature
       Headspace, e.g., air vs. O2 vs. anaerobic (currently not considered)
       Nutrient requirements (currently not considered)
       Cell composition (currently not considered)
'''

'''
TODO:
   Separation units
   Check equipment costs
'''


# %% Setup

import biosteam as bst
import thermosteam as tmo
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
from biosteam.units._splitter import run_split_with_mixing
from thermosteam import Stream, MultiStream
from orgacids.EB import EB_simulation

# Used by Sarang's unit
#from scipy.optimize import fsolve
#import math

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
_lb2MT = 0.000453592 # MT is metric tonne
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
_hp2kW = 0.7457
_Gcal2kJ = 4184e3


# %% Feedstock handling

# Feedstock handling system as a whole, capital and operating costs considered in feedstock cost
# CE = 521.9 is year 2009
@cost(basis='Flow rate', ID='System', units='kg/hr',
      S=94697, ub=False, CE=521.9, cost=3329690, n=0.6, kW=783, BM=1.7)
class FeedstockHandling(Unit): pass


# %% Pretreatment

# Sulfuric acid addition tank
# CE = 550.8 is year 2010
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=550.8, cost=6210, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=3720, ub=False, CE=521.9, cost=8000, n=0.8, kW=1, BM=2.3)
class SulfuricAcidAdditionTank(Unit): pass

# Sulfuric acid in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=136260, ub=False, CE=521.9, cost=6000, n=0.6, kW=0, BM=1)
class SulfuricAcidMixer(Unit):
    _N_ins = 2
    _N_outs = 1
    _graphics = bst.units.Mixer._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
    def _run(self):
        acid, water = self.ins
        mixture = self.outs[0]
        # 0.05 is from 1842/36629 from streams 710 and 516 of Humbird et al.
        water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
        mixture.mix_from([water, acid])
        
# Blowdown discharge pump
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=25365, n=0.8, kW=93.21, BM=2.3)
class BlowdownDischargePump(Unit): pass

# Pretreatment flash tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=511000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=30000, n=0.8, kW=55.9275, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=252891, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class PretreatmentFlash(bst.units.Flash): pass

# Oligomer conversion tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class OligomerConversionTank(Unit): pass

# Ammonia in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=157478, ub=False, CE=521.9, cost=5000, n=0.5, kW=0, BM=1)
class AmmoniaMixer(bst.units.Mixer): pass

# Ammonia addition tank
# Assumed all sulfuric acid has been reacted with biomass or neutralized by ammonia
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=236000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=410369, ub=False, CE=521.9, cost=21900, n=0.5, kW=7.457, BM=1.5)
class AmmoniaAdditionTank(Unit): 
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        self.neutralization_rxns = Rxn('2 NH3 + H2SO4 -> NH4SO4 + 2 H2O', 'NH3', 1.)
    
    def _run(self):
        ins = self.ins[0]
        outs = self.outs[0]
        outs.copy_like(ins)
        
        self.neutralization_rxns(outs.mol)
        outs.imol['H2SO4'] = 0
            
# Waste vapor condenser
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-2, ub=False, CE=521.9, cost=34000, n=0.7, kW=0, BM=2.2)
class WasteVaporCondenser(bst.units.HXutility):
    _graphics = bst.units.HXutility._graphics

# Steam mixer
class SteamMixer(Unit):
    """
    **ins**
    
        [0] Feed
        
        [1] Steam
    
    **outs**
    
        [0] Mixed
    
    """
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1
    def __init__(self, ID='', ins=None, outs=(), *, P):
        super().__init__(ID, ins, outs)
        self.P = P
    
    @staticmethod
    def _P_at_flow(mol_water, P, steam, mixed, feed):
        steam.imol['Water'] = mol_water
        mixed.mol[:] = steam.mol + feed.mol
        mixed.H = feed.H + mol_water * 40798 #!!! not sure what 40798 is
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P - P_new
    
    def _run(self):
        feed, steam = self.ins
        steam_mol = steam.F_mol
        mixed = self.outs[0]
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam, mixed, feed))
        mixed.P = self.P
        hu = self.heat_utilities[0]
        hu.agent = hu.get_heating_agent('low_pressure_steam')
        hu.flow = steam_mol
        hu.cost = steam_mol*bst.HeatUtility.get_heating_agent('low_pressure_steam').regeneration_price
    
    @property
    def installation_cost(self): return 0
    def purchase_cost(self): return 0
    def _design(self): pass
    def _cost(self): pass
        
# Pretreatment reactor
@cost('Dry flow rate', units='kg/hr', S=83333, CE=521.9, ub=False,
      cost=19812400, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.units.Flash._graphics
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = MultiStream(None)
        self.pretreatment_rxns = ParallelRxn([
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
    Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0050)
    ])
        vapor, liquid = self.outs
        vapor.phase = 'g'
    
    def _run(self):
        ms = self._multistream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copy_like(feed)
        self.pretreatment_rxns(liquid.mol) 
        ms.copy_like(liquid)
        H = ms.H + ms.Hf - feed.Hf
        ms.vle(T=130+273.15, H=H)
        vapor.mol[:] = ms.imol['g']
        liquid.mol[:] = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P


# %% Saccharification and fermentation

# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=40414, ub=False, CE=521.9, cost=439000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=43149, ub=False, CE=521.9, cost=8200, n=0.8, kW=7.457, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=40414, ub=False, CE=521.9, cost=31800, n=0.5, kW=11.1855, BM=1.5)
class SeedHoldTank(Unit): pass

# Hydrolysate cooler
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      S=-8, ub=False, CE=521.9, cost=85000, n=0.7, kW=0, BM=2.2)
class HydrolysateCooler(bst.units.HXutility): pass

# Enzyme hydrolysate mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      S=380000, ub=False, CE=521.9, cost=109000, n=0.5, kW=74.57, BM=1.7)
class EnzymeHydrolysateMixer(bst.units.Mixer): pass

# Seed train (5 stages)
# Assumes no byproducts in seed train and no need for neutralization
@cost('Flow rate', 'Pumps',
      S=43149, CE=521.9, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=521.9, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=521.9, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=521.9, n=0.5, BM=1.5)
class SeedTrain(Unit):  
    _N_heat_utilities = 1
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    # Number of parallel seed trains
    N_trains = 2

    # Operating temperature (K), updated for the target acid
    T = 30+273.15
    
    # Cycle time for each batch (hr)
    tau_batch = 24
    
    # Whether to use EB model for simulation
    # If not using EB model, substrate conversions can be changed in reactions
    EB = True

    # Intial biomass loading to get things started, won't be needed later
    X_init_g = 0.5 # [g/L] (i.e., [kg/m3])
    X_init_x = 0.5 # [g/L] (i.e., [kg/m3])

    # Kinetic parameters are (n, P_max, Y_XS, Y_PS, mu_max, K_S) 
    kinetic_parameters_g = (3, 61.5, 0.08, 0.27, 0.228, 0)
    kinetic_parameters_x = (3, 61.5, 0.06, 0.66, 0.115, 45)
    
    # Byproduct-to-organic-acid ratios for (EtOH, AceA) (g/g)
    byproduct_ratios_g = (1.04, 0.)
    byproduct_ratios_x = (0., 0.)

    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5

    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    def __init__(self, ID='', ins=None, outs=(), EB=True,
                  X_init_g=0.5, X_init_x=0.5,
                  kinetic_parameters_g=(3, 61.5, 0.08, 0.27, 0.228, 0), 
                  kinetic_parameters_x=(3, 61.5, 0.06, 0.66, 0.115, 45),
                  byproduct_ratios_g=(1.04, 0.), 
                  byproduct_ratios_x=(0., 0.)):
        Unit.__init__(self, ID, ins, outs)

        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.27),
        Rxn('Glucose -> 2 Ethanol + 2 CO2',   'Glucose',   0.55),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> \
                 6 FermMicrobeGlu + 2.4 H2O', 'Glucose',   0.098),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.66),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',  'Xylose',    0),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> \
                 5 FermMicrobeXyl + 2 H2O',   'Xylose',    0.073)
        ])

    def _run(self):
        feed = self.ins[0]
        if feed.imass['FermMicrobeGlu'] == 0: 
            feed.imass['FermMicrobeGlu'] = self.X_init_g * feed.F_vol
        if feed.imass['FermMicrobeXyl'] == 0: 
            feed.imass['FermMicrobeXyl'] = self.X_init_x * feed.F_vol
        vent, effluent= self.outs     
        effluent.copy_flow(feed)
        cofermentation_rxns = self.cofermentation_rxns
        
        # Use EB model to update conversions
        if self.EB == True:
            kinetic_parameters_g = self.kinetic_parameters_g
            kinetic_parameters_x = self.kinetic_parameters_x
            byproduct_ratios_g = self.byproduct_ratios_g
            byproduct_ratios_x = self.byproduct_ratios_x
            ins_concentration = effluent.imass['FermMicrobeGlu', 'FermMicrobeXyl',
                                               'LacticAcid', 'Glucose', 'Xylose']
            mass_in = tuple(ins_concentration/effluent.F_vol)
            mass_in_0 = mass_in

            for i in range(0, 5):
                # Outputs are (X_gf, X_xf, Pf, S_gf, S_xf, P_gf, P_xf) in kg/hr
                # Index         0     1    2    3     4     5     6
                mass_out = EB_simulation(tau=self.tau_batch, feed=mass_in,
                                         kinetic_parameters_g=kinetic_parameters_g,
                                         kinetic_parameters_x=kinetic_parameters_x)
                mass_in = mass_out[0:5]
         
            # Converting from mass basis to molar basis
            cofermentation_rxns.X[0] = (mass_out[5]/(mass_in_0[3]-mass_out[3])) \
                                       / ((2*feed.chemicals.LacticAcid.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[1] = (mass_out[5]*byproduct_ratios_g[0]/(mass_in_0[3]-mass_out[3])) \
                                       / ((2*feed.chemicals.Ethanol.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[2] = (mass_out[5]*byproduct_ratios_g[1]/(mass_in_0[3]-mass_out[3])) \
                                       / ((3*feed.chemicals.AceticAcid.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[3] = ((mass_out[0]-mass_in_0[0])/(mass_in_0[3]-mass_out[3])) \
                                       / ((1*feed.chemicals.FermMicrobeGlu.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[4] = (mass_out[6]/(mass_in_0[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.LacticAcid.MW) / (3*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[5] = (mass_out[6]*byproduct_ratios_x[0]/(mass_in_0[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.Ethanol.MW) / (3*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[6] = (mass_out[6]*byproduct_ratios_x[1]/(mass_in_0[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.AceticAcid.MW) / (2*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[7] = ((mass_out[1]-mass_in_0[1])/(mass_in_0[4]-mass_out[4])) \
                                       / ((1*feed.chemicals.FermMicrobeXyl.MW) / (1*feed.chemicals.Xylose.MW))

        cofermentation_rxns(effluent.mol)
        loss = feed.F_mass - effluent.F_mass
        # Assume all loss is CO2
        effluent.imass['CO2'] += loss
        effluent.T = self.T
        vent.phase = 'g'
        vent.recieve_vent(effluent)

    def _design(self): 
        maxvol = self.outs[1].F_vol*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self.design_results
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.F_mass for i in self.outs])
        self.heat_utilities[0](self.Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self.design_results
        C = self.purchase_costs
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self.power_utility(kW)

# Saccharification and co-fermentation reactor,
# "co-fermentation" means the microbe can simultaneously ferment glucose and xylose
#!!! Did not include heat exchanger (or did?)
@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr, ub=False,
      cost=47200, n=0.8, BM=2.3, CE=521.9, N='N_recirculation_pumps')
@cost('Reactor duty', 'Heat exchangers', CE=521.9, cost=23900, ub=False,
      S=5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=521.9, cost=52500, ub=False,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=521.9, cost=844000, ub=False,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr, ub=False,
      cost=47200/5, CE=521.9, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, ub=False,
      CE=521.9, n=0.7, BM=2.0, N='N_tanks')
class SaccharificationAndCoFermentation(Unit):    
    _N_ins = 4
    _N_outs = 3
    _N_heat_utilities = 2
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    # Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    # Fermentation temperature (K), updated for the target acid
    T_fermentation = 30+273.15
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    # Saccharification time (hr)
    tau_saccharification = 60
    
    # Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    # Unload and clean up time (hr)
    tau_0 = 4
    
    # Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    # Number of reactors
    N_reactors = 12
    
    # Number of continuous saccharification tanks
    N_tanks = 8
    
    # Number of transfer pumps
    N_transfer_pumps = 5
    
    # Number of recirculation pumps
    N_recirculation_pumps = 5
    
    # Split to outs[2]
    saccharified_slurry_split = 0.1
    
    # Whether to use EB model for simulation
    # If not using EB model, substrate conversions can be changed in reactions
    EB = True
    
    # Intial biomass loading to get things started, won't be needed later
    X_init_g = 0.5 # [g/L] (i.e., [kg/m3])
    X_init_x = 0.5 # [g/L] (i.e., [kg/m3])
    
    # Kinetic parameters are (n, P_max, Y_XS, Y_PS, mu_max, K_S) 
    kinetic_parameters_g = (3, 61.5, 0.08, 0.27, 0.228, 0)
    kinetic_parameters_x = (3, 61.5, 0.06, 0.66, 0.115, 45)
    
    # Byproduct-to-organic-acid ratios for (EtOH, AceA) (g/g)
    byproduct_ratios_g = (1.04, 0.)
    byproduct_ratios_x = (0., 0.)
    
    def __init__(self, ID='', ins=None, outs=(), P=101325,
                 EB=True, X_init_g=0.5, X_init_x=0.5,
                 kinetic_parameters_g=(3, 61.5, 0.08, 0.27, 0.228, 0), 
                 kinetic_parameters_x=(3, 61.5, 0.06, 0.66, 0.115, 45),
                 byproduct_ratios_g=(1.04, 0.), 
                 byproduct_ratios_x=(0., 0.)):
        Unit.__init__(self, ID, ins, outs)
    
        self.P = P
        
        # Enzymatic hydrolysis reactions including from downstream batch tank in co-fermentation.
        # Kept the same as Humbird et al.
        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant   Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000),
            Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1.0000)
            ])
        
        self.saccharified_stream = Stream(None)
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.27),
        Rxn('Glucose -> 2 Ethanol + 2 CO2',   'Glucose',   0.55),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> \
                 6 FermMicrobeGlu + 2.4 H2O', 'Glucose',   0.098),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.66),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',  'Xylose',    0),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> \
                 5 FermMicrobeXyl + 2 H2O',   'Xylose',    0.073)
        ])
        
        # # Assume all CSL and DAP are consumed
        # self.CSL_DAP_rxns = ParallelRxn([
        #     #   Reaction definition                              Reactant    Conversion
        #     Rxn('CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL',    1.0000),
        #     # Assume biomass picks up phosphate
        #     Rxn('DAP -> 2 NH3',                                    'DAP',    1.0000)
        #     ])
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
            #   Reaction definition                                             Reactant  Conversion
            Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O', 'LacticAcid',   1.),
            Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O', 'AceticAcid',   1.)
            ])

    def _run(self):
        feed, CSL, DAP, lime = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        ss.copy_flow(feed)
        
        self.saccharification_rxns(ss.mol)
        # Sidedraw after saccharification but before co-fermentation
        sidedraw.mol = ss.mol * self.saccharified_slurry_split
        # Now effluent is after saccharification but before co-fermentation
        effluent.mol = ss.mol - sidedraw.mol + CSL.mol + DAP.mol
        if effluent.imass['FermMicrobeGlu'] == 0: 
            effluent.imass['FermMicrobeGlu'] = self.X_init_g * effluent.F_vol
        if effluent.imass['FermMicrobeXyl'] == 0: 
            effluent.imass['FermMicrobeXyl'] = self.X_init_x * effluent.F_vol
        cofermentation_rxns = self.cofermentation_rxns
        
        # Use EB model to update conversions
        if self.EB == True:
            kinetic_parameters_g = self.kinetic_parameters_g
            kinetic_parameters_x = self.kinetic_parameters_x
            byproduct_ratios_g = self.byproduct_ratios_g
            byproduct_ratios_x = self.byproduct_ratios_x
            ins_concentration = effluent.imass['FermMicrobeGlu', 'FermMicrobeXyl',
                                               'LacticAcid', 'Glucose', 'Xylose']
            mass_in = tuple(ins_concentration/effluent.F_vol)
            # Outputs are (X_gf, X_xf, Pf, S_gf, S_xf, P_gf, P_xf) in kg/hr
            # Index         0     1    2    3     4     5     6
            mass_out = EB_simulation(tau=self.tau_cofermentation, feed=mass_in,
                                     kinetic_parameters_g=kinetic_parameters_g,
                                     kinetic_parameters_x=kinetic_parameters_x)
            
          # Converting from mass basis to molar basis
            cofermentation_rxns.X[0] = (mass_out[5]/(mass_in[3]-mass_out[3])) \
                                       / ((2*feed.chemicals.LacticAcid.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[1] = (mass_out[5]*byproduct_ratios_g[0]/(mass_in[3]-mass_out[3])) \
                                       / ((2*feed.chemicals.Ethanol.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[2] = (mass_out[5]*byproduct_ratios_g[1]/(mass_in[3]-mass_out[3])) \
                                       / ((3*feed.chemicals.AceticAcid.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[3] = ((mass_out[0]-mass_in[0])/(mass_in[3]-mass_out[3])) \
                                       / ((6*feed.chemicals.FermMicrobeGlu.MW) / (1*feed.chemicals.Glucose.MW))
            cofermentation_rxns.X[4] = (mass_out[6]/(mass_in[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.LacticAcid.MW) / (3*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[5] = (mass_out[6]*byproduct_ratios_x[0]/(mass_in[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.Ethanol.MW) / (3*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[6] = (mass_out[6]*byproduct_ratios_x[1]/(mass_in[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.AceticAcid.MW) / (2*feed.chemicals.Xylose.MW))
            cofermentation_rxns.X[7] = ((mass_out[1]-mass_in[1])/(mass_in[4]-mass_out[4])) \
                                       / ((5*feed.chemicals.FermMicrobeXyl.MW) / (1*feed.chemicals.Xylose.MW))
        
        mass_before_fermentation = effluent.F_mass
        cofermentation_rxns(effluent.mol)
        # Assume all lost glucose and xylose changed to CO2
        effluent.imass['CO2'] += (mass_before_fermentation - effluent.F_mass)
        
        # # Changed from effluent.mass to effluent.mol
        # self.CSL_DAP_rxns(effluent.mol)
        
        # Set feed lime mol to match rate of acids production, add 5% extra
        lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
                            + effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1]) \
                            * 1.05
        effluent.mol += lime.mol
        self.neutralization_rxns(effluent.mol)
        vent.recieve_vent(effluent)
    
    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Tank volume'] = v_0*self.tau_tank/self.V_wf/self.N_tanks
        Design['Flow rate'] = v_0/self.N_transfer_pumps
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        hu_cooling, hu_fermentation = self.heat_utilities
        mixture = self.thermo.mixture
        ss = self.saccharified_stream
        mol = ss.mol
        hu_cooling(mixture.H('l', mol, self.T_fermentation, 101325) \
                   - mixture.H('l', mol, self.T_saccharification, 101325),
                   self.T_fermentation
                   )
        diff = sum([i.imol['CalciumLactate'] for i in self.outs]) \
             + sum([i.imol['CalciumAcetate'] for i in self.outs]) \
             - sum([i.imol['CalciumLactate'] for i in self.ins]) \
             - sum([i.imol['CalciumAcetate'] for i in self.ins])
        duty = diff * -5568 #??? Where does this come from?
        hu_fermentation(duty, effluent.T)
        Design['Reactor duty'] = -duty


# %% Organic acid separation

@cost('Flow rate', 'Flitrate tank agitator', ub=False,
      cost=26e3, CE=550.8, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump', ub=False,
      cost=13040, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank', ub=False,
      cost=103e3, S=31815, CE=550.8, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57, ub=False,
      cost= 18173, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531', ub=False,
      cost=174800, CE=550.8, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Manifold flush pump', kW=74.57, ub=False,
      cost=17057, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank', ub=False,
      cost=1520, CE=550.8, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Lignin wet cake screw',  kW=15*_hp2kW, ub=False,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Lignin wet cake conveyor', kW=10*_hp2kW, ub=False,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter', ub=False,
      cost=3294700, CE=550.8, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank', ub=False,
      cost=8e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW, ub=False,
      cost=29154, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank', ub=False,
      cost=17e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter', ub=False,
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)', ub=False,
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
class CellMassFilter(bst.units.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _run(self):
        cell_mass, recycled_water = self.ins
        # No information on recycled_water in flow table in Humbird et al.,
        # therefore 36538 is taken from the scaling basis of equipment T-532
        # 391501 is total flow of stream 508 in Humbird et al.
        recycled_water.imass['Water'] = 36538 * cell_mass.F_mass / 391501
        
        run_split_with_mixing(self)
        retentate, permeate = self.outs
        solids = retentate.F_mass
        mc = self.mositure_content
        retentate.imol['Water'] = water = (solids * mc/(1-mc))/18.01528
        permeate.imol['Water'] -= water
        if permeate.imol['Water'] < water:
            raise ValueError(f'not enough water for {repr(self)}')
        
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass

# Cost copied from OligomerConversionTank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class AcidulationReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.acidulation_rxns = ParallelRxn([
            #   Reaction definition                                    Reactant   Conversion
            Rxn('CalciumLactate + H2SO4 -> 2 LacticAcid + CaSO4', 'CalciumLactate', 0.95),
            Rxn('CalciumAcetate + H2SO4 -> 2 AceticAcid + CaSO4', 'CalciumAcetate', 0.95)
            ])
        
    def _run(self):
        feed, acid = self.ins
        effluent = self.outs[0]
        
        effluent.copy_like(feed)
        # Set feed acid mol to match acidulation needs with 5% extra
        acid.imol['H2SO4'] = (feed.imol['CalciumLactate']/self.acidulation_rxns.X[0] \
                             +feed.imol['CalciumAcetate']/self.acidulation_rxns.X[1]) \
                             * 1.05
        effluent.mix_from([feed, acid])
        self.acidulation_rxns(effluent.mol)
        
#!!! Now costs are the same as the CellMassFilter, but should consider changing to
# S-222, the hydroclone & rotary drum filter as in Aden et al.        
@cost('Flow rate', 'Flitrate tank agitator', ub=False,
      cost=26e3, CE=550.8, kW=7.5*_hp2kW, S=31815, n=0.5, BM=1.5)
@cost('Flow rate', 'Discharge pump', ub=False,
      cost=13040, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Filtrate tank', ub=False,
      cost=103e3, S=31815, CE=550.8, BM=2.0, n=0.7)
@cost('Flow rate', 'Feed pump', kW=74.57, ub=False,
      cost= 18173, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost('Flow rate', 'Stillage tank 531', ub=False,
      cost=174800, CE=550.8, S=31815, n=0.7, BM=2.0)
@cost('Flow rate', 'Manifold flush pump', kW=74.57, ub=False,
      cost=17057, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Recycled water tank', ub=False,
      cost=1520, CE=550.8, S=31815, n=0.7, BM=3.0)
@cost('Flow rate', 'Gypsum wet cake screw',  kW=15*_hp2kW, ub=False,
      cost=2e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Gypsum wet cake conveyor', kW=10*_hp2kW, ub=False,
      cost=7e4, CE=521.9, S=28630, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressure filter', ub=False,
      cost=3294700, CE=550.8, S=31815, n=0.8, BM=1.7)
@cost('Flow rate', 'Pressing air compressor reciever tank', ub=False,
      cost=8e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Cloth wash pump', kW=150*_hp2kW, ub=False,
      cost=29154, CE=550.8, S=31815, n=0.8, BM=2.3)
@cost('Flow rate', 'Dry air compressor reciever tank', ub=False,
      cost=17e3, CE=550.8, S=31815, n=0.7, BM=3.1)
@cost('Flow rate', 'Pressing air pressure filter', ub=False,
      cost=75200, CE=521.9, S=31815, n=0.6, kW=112, BM=1.6)
@cost('Flow rate', 'Dry air pressure filter (2)', ub=False,
      cost=405000, CE=521.9, S=31815, n=0.6, kW=1044, BM=1.6)
#!!! Need to look into Aden et al. and see how the filter was designed there
class GypsumFilter(bst.units.SolidsSeparator):
    _units = {'Flow rate': 'kg/hr'}
    
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass
        
class OrganicAcidsFlash(bst.units.Flash):
    vle_chemicals = ()
    V = None
    P = 101325
    
    def _run(self):
        feed = self.ins[0]
        vap, liq = self.outs
        
        # Vapor-liquid equilibrium
        ms = self._multistream
        ms.empty()
        vle_chemicals = self.vle_chemicals
        ms.imol['l', vle_chemicals] = feed.imol[vle_chemicals]
        ms.vle(P=self.P, V=self.V)

        # Set values
        vap.phase = 'g'
        liq.phase = 'l'
        liq.mol = feed.mol
        vap.imol[vle_chemicals] = ms.imol['g', vle_chemicals]
        liq.imol[vle_chemicals] = ms.imol['l', vle_chemicals]
        vap.T = liq.T = ms.T
        vap.P = liq.P = ms.P
        ms.show(N=100)
        ms.imol['l'] = liq.mol
        ms.show(N=100)
        ms.imol['g'] = vap.mol
        ms.show(N=100)

# Cost copied from OligomerConversionTank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class EsterificationReactor(Unit):
    _N_ins = 2
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.esterification_rxns = ParallelRxn([
            #   Reaction definition                               Reactant  Conversion
            Rxn('LacticAcid + Methanol -> MethylLactate + H2O', 'LacticAcid', 0.98),
            Rxn('AceticAcid + Methanol -> MethylAcetate + H2O', 'AceticAcid', 0.98)
            ])
        
    # def _run(self):
    #     feed, recycled, methanol = self.ins
    #     effluent = self.outs[0]
        
    #     # Add 5% extra
    #     methanol_needed = (feed.imol['LacticAcid']/self.esterification_rxns.X[0] \
    #                         +feed.imol['AceticAcid']/self.esterification_rxns.X[1]) \
    #                       * 1.05
                          
    #     methanol.imol['Methanol'] = max(0, methanol_needed-recycled.imol['Methanol'])
    #     effluent.mix_from([feed, recycled, methanol])
    #     self.esterification_rxns(effluent.mol)
    
    def _run(self):
        feed, methanol = self.ins
        effluent = self.outs[0]
        
        # Add 5% extra
        methanol_needed = (feed.imol['LacticAcid']/self.esterification_rxns.X[0] \
                            +feed.imol['AceticAcid']/self.esterification_rxns.X[1]) \
                          * 1.05
                          
        methanol.imol['Methanol'] = methanol_needed
        effluent.mix_from([feed, methanol])
        self.esterification_rxns(effluent.mol)

# Cost copied from OligomerConversionTank
#!!! Need to add ethanol into consideration?
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=203000, n=0.7, kW=0, BM=2)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=292407, ub=False, CE=550.8, cost=17408, n=0.8, kW=55.9, BM=1.7)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=264116, ub=False, CE=521.9, cost=90000, n=0.5, kW=170, BM=1.5)
class HydrolysisReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=(), outs=()):
        Unit.__init__(self, ID, ins, outs)
        # Conversions are placeholders now
        self.hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                Reactant    Conversion
            Rxn('MethylLactate + H2O -> LacticAcid + Methanol', 'MethylLactate', 0.98),
            Rxn('MethylAcetate + H2O -> AceticAcid + Methanol', 'MethylAcetate', 0.98)
            ])
    
    def _run(self):
        feed, water = self.ins
        effluent = self.outs[0]
        
        # Add 5% extra
        water_needed = (feed.imol['MethylLactate']/self.hydrolysis_rxns.X[0] \
                        +feed.imol['MethylAcetate']/self.hydrolysis_rxns.X[1]) \
                       * 1.05
        water.imol['Water'] = max(0, (water_needed-feed.imol['Water'])) 
        effluent.mix_from([feed, water])
        
        self.hydrolysis_rxns(effluent.mol)


# %% Wastewater treatment

# The total cost of wastewater treatment is combined into this placeholder
@cost('Flow rate', 'Wastewater system', units='kg/hr', CE=550.8, ub=False,
      cost=50280080., n=0.6, BM=1, kW=7139/1.05, S=393100)
class WastewaterSystemCost(Unit): pass

class AnaerobicDigestion(Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **digestion_rxns:** [ReactionSet] Anaerobic digestion reactions.
        
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
    
    def __init__(self, ID='', ins=None, outs=(), *, digestion_rxns, sludge_split):
        Unit.__init__(self, ID, ins, outs)
        self.digestion_rxns = digestion_rxns
        self.sludge_split = sludge_split
        self.multi_stream = MultiStream()
    
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link_with(cool_water, TP=False)
        biogas.T = waste.T = sludge.T = T = 35+273.15
        # Hot water temperature is adjusted to maintain heat balance
        hot_water.T = feed.T - 5
        H_at_35C = feed.thermo.mixture.H(z=feed.mol, phase='l', T=T, P=101325)
        cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        while cool_water.F_mol < 0:
            hot_water.T -= 1
            cool_water.mol[:] *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        sludge.copy_flow(feed)
        self.digestion_rxns(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol[:] = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol[:] = liquid_mol * self.sludge_split
        waste.mol[:] = liquid_mol - sludge.mol
        biogas.recieve_vent(waste)     
    
class AerobicDigestion(Unit):
    """Anaerobic digestion system as modeled by Humbird 2011
    
    **Parameters**
    
        **digestion_rxns:** [ReactionSet] Anaerobic digestion reactions.
        
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
    
    def __init__(self, ID='', ins=None, outs=(), *, digestion_rxns):
        Unit.__init__(self, ID, ins, outs)
        self.digestion_rxns = digestion_rxns
    
    def _run(self):
        waste, air, caustic = self._ins
        vent, water = self.outs
        vent.phase = 'g'
        water.copy_like(waste)
        water.mol[:] += air.mol
        water.mol[:] += caustic.mol
        self.digestion_rxns(water.mol)
        vent.copy_flow(water, ('CO2', 'O2', 'N2'))
        vent.imol['Water'] = water.imol['Water'] * self.evaporation
        water.mol[:] -= vent.mol


# %% Facility units in bst.units

# Sulfuric acid storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1981, ub=False, CE=550.8, cost=96000, n=0.7, kW=0, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1981, ub=False, CE=521.9, cost=7493, n=0.8, kW=0.37285, BM=2.3)
class SulfuricAcidStorageTank(bst.units.StorageTank): pass

# Ammonia storage tank, no pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1171, ub=False, CE=550.8, cost=196000, n=0.7, kW=0, BM=2)
class AmmoniaStorageTank(bst.units.StorageTank): pass

# CSL storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=1393, ub=False, CE=521.9, cost=70000, n=0.7, kW=0, BM=2.6)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=1393, ub=False, CE=550.8, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=1393, ub=False, CE=521.9, cost=21200, n=0.5, kW=7.457, BM=1.5)
class CSLStorageTank(Unit): pass

# DAP storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=102000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=3000, n=0.8, kW=0.37285, BM=3.1)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=9800, n=0.5, kW=4.10135, BM=1.5)
@cost(basis='Flow rate', ID='Bag', units='kg/hr',
      S=163, ub=False, CE=521.9, cost=30000, n=0.6, kW=0, BM=1.7)
class DAPStorageTank(Unit): pass

# For storage of lime used in separation and waste treatment,
# design copied from cornstover biorefinery in Aden et al.
# CE = 386.5 and 389.5 is year 1997 and 1998, respectively
# Base flow from 227 of Aden et al.
@cost(basis='Flow rate', ID='Storage bin', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=136370, n=0.46, kW=0, BM=1.3)
@cost(basis='Flow rate', ID='Dust vent baghouse', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=140707, n=1, kW=0, BM=1.5)
# Power usage scaled based on M-106 of Humbird et al.
@cost(basis='Flow rate', ID='Unloading blower', units='kg/hr',
      S=227, ub=False, CE=389.5, cost=99594, n=0.5, kW=25*_hp2kW*7425/8500, BM=1.4)
# Cost not scaled, thus used n=0
# Power usage scaled based on C-101 of Humbird et al., 
# though C-101 is a belt and here is a rotary type
@cost(basis='Flow rate', ID='Feeder', units='kg/hr',
      S=227, ub=False, CE=386.5, cost=3900, n=0, kW=20*_hp2kW*7000*_lb2MT/160, BM=1.3)
class LimeStorageTank(Unit): pass

# Methanol storage tank, based on Ethanol storage tank from Humbird et al.
# scaled down from 7-d to 14-hr storage
# (similar to sulfuric acid/CSL/DAP storage in Humbird et al.)
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=21808/24/7*14, ub=False, CE=521.9, cost=1340000, n=0.7, kW=0, BM=1.7)
class MethanolStorageTank(Unit): pass

# Fire water tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      S=8343, ub=False, CE=521.9, cost=803000, n=0.7, kW=0, BM=1.8)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      S=8343, ub=False, CE=521.9, cost=15000, n=0.8, kW=93.2125, BM=1.7)
class FireWaterTank(Unit): pass


# %% Sarang's work in progress
'''
class MultiComponentDistillation(Unit, isabstract=True):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 P=101325, *, LHK, rLKD, rHKB, K_HK, K_specs):
        super().__init__(ID, ins, outs)
        self.outs[0].phase = 'g'
        
        self.P = P_in = ins.P
        self.T = T_in = ins.T
        self.LHK = LHK
        self.rLKD = rLKD
        self.rHKB = rHKB
        self.K_HK = K_HK
        self.K_specs = K_specs
        
        # equilibrium_chemicals should be an iterable of Chemical objects
        #??? Why use Water and Ethanol?
        equilibrium_chemicals = tmo.Chemicals(['Ethanol', 'Water'])
        F_l = tmo.equilibrium.LiquidFugacities(equilibrium_chemicals)
        liquid_molar_composition = np.array([0.9, 0.1])
        fugacities = F_l(x=liquid_molar_composition, T=T_in)
        
    def GetRelativeVolatilities(self):
        K_HK = self.K_HK
        K_specs = self.K_specs
        
        rel_vols = []
        for K_i in K_specs:
            rel_vols.append(K_i/K_HK)        
       
        self.rel_vols = rel_vols
        
        return rel_vols
    
    def EstimateMinStages(self):
        LHK = self.LHK
        rLKD = self.rLKD
        rHKB = self.rHKB
        
        alpha_LH = self.rel_vols[self.ins[0].chemicals.index(LHK[0])]
        Nm = math.log10((rLKD/(1-rHKB)) * (rHKB/(1-rLKD))) \
            /math.log10(alpha_LH)
        
        self.alpha_LH = alpha_LH
        self.Nm = Nm
        
        return Nm
    
    def EstimateDBCompositions(self):
        rHKB = self.rHKB
        Nm = self.Nm
        rel_vols = self.rel_vols
        
        A = math.log10((1-rHKB) / rHKB)*1**(-Nm)
        D_by_B_ratios = []
        for rv_i in rel_vols:
            D_by_B_ratios.append(10**(A + Nm*math.log10(rv_i)))        
        
        self.D_by_B_ratios = D_by_B_ratios
        
        DB = []
        mass = self.ins[0].F_mass
        for ind in range(len(D_by_B_ratios)):
            F = self.ins[0].mass[ind] / mass
            d = F / (1 + 1 / (D_by_B_ratios[ind]))
            b = F - d
            DB.append((d,b))
        
        self.DB = DB
        
        return DB
    
    def Underwood_1_single(self, ind, theta):
        return self.rel_vols[ind] * (self.DB[ind][0] + self.DB[ind][1]) \
               / (self.rel_vols[ind] - theta)

    """ Underwood Ean. # 1"""
    #!!! Ean. or Eqn.?
    def Underwood_1(self, theta):
        sum = 0.0000
        
        for ind in range(len(self.DB)):
            sum += self.Underwood_1_single(ind, theta)
            
        return sum

    def Underwood_2_single(self, ind):
        return self.rel_vols[ind] * (self.DB[ind][0]) / (self.rel_vols[ind] - self.theta)

    """Underwood Ean. #2"""
    def Underwood_2(self, Rm):
        sum = 0.0000
        
        for ind in range(len(self.DB)):
            sum += self.Underwood_2_single(ind)
        
        sum -= (Rm + 1)
            
        return sum

    """Uses Underwood Eqns."""
    def EstimateOperatingRefluxRatio(self):
        theta = fsolve(self.Underwood_1, [1, 1, 1])
        self.theta = theta
        
        Rm = fsolve(self.Underwood_2, [1, 1, 1])
        self.Rm = Rm
        
        self.R_heuristic = 1.2
        self.R_op = self.Rm * self.R_heuristic
        
        return self.R_op
    
    """Gilliland Eqn."""
    def Gilliland(self, N):
        N_ratio = ((N-self.Nm) / (N+1))
        fi = (self.R_op-self.Rm) / (self.R_op+1)
        exp = math.exp(((1+54.4*fi)/(11+117.2*fi)) * ((fi-1)/(fi**0.5)))
        
        return N_ratio - 1 + exp
        
    """Uses Gilliland Eqn."""
    def EstimateNumTheoreticalStages(self):
        N = fsolve(self.Gilliland, [4, 4, 4])
        
        return N

    """ WIP: Kirkbride Approximation for Feed Tray Location"""
    def Kirkbride(self): #WIP
        Nf = 1
        
        return Nf

    # """WIP"""
    # def _mass_balance(self):
    #     vap, liq = self.outs
        
    #     # Get all important flow rates (both light and heavy keys and non-keys)
    #     LHK_index = self.LHK_index
    #     LNK_index = self.LNK_index
    #     HNK_index = self.HNK_index
    #     mol = self._mol_in
        
    #     LHK_mol = mol[LHK_index]
    #     LNK_mol = mol[LNK_index]
    #     HNK_mol = mol[HNK_index]
        
    #     # Set light and heavy keys by lever rule
    #     light, heavy = LHK_mol
    #     LHK_mol = light + heavy
    #     zf = light / LHK_mol
    #     split_frac = (zf-self.x_bot) / (self.y_top-self.x_bot)
    #     top_net = LHK_mol * split_frac
        
    #     # Set output streams
    #     vap.mol[LHK_index] = top_net * self._y
    #     liq.mol[LHK_index] = LHK_mol - vap.mol[LHK_index]
    #     vap.mol[LNK_index] = LNK_mol
    #     liq.mol[HNK_index] = HNK_mol
        
    """WIP"""
    def _run(self):
        ins_0 = self.ins[0]
        vap, liq = self.outs
        LHK_0 = self.LHK[0]
        
        rel_vols = self.GetRelativeVolatilities()
        Nm = self.EstimateMinStages()
        DB = self.EstimateDBCompositions()
        R_op = self.EstimateOperatingRefluxRatio()
        
        vap.copy_flow(ins_0)
        liq.copy_flow(ins_0)
        
        FLK = ins_0.mass[ins_0.chemicals.index(LHK_0)]
        F = self.ins[0].F_mass
        
        dLK = DB[vap.chemicals.index(LHK_0)][0]
        bLK = DB[vap.chemicals.index(LHK_0)][1]
        
        D = (FLK-F*bLK) / (dLK-bLK)
        B = F - D
        print(D)
        print(F)
        print(B)
        
        for spec in vap.chemicals.IDs:
            vap.imol[spec] = DB[vap.chemicals.index(spec)][0] * D
            liq.imol[spec] = DB[vap.chemicals.index(spec)][1] * B
'''



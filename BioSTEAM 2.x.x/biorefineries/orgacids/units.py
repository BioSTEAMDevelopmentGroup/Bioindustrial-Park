#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:19:18 2019

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

@author: yalinli_cabbi
"""

'''
TODO:
   Vet separation system design
   Need a Reactor class for Acidulation, Esterification, Hydrolysis
'''


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank
from biosteam.units.decorators import cost
from biosteam.units._tank import compute_number_of_tanks_and_total_purchase_cost
from biosteam.units.design_tools import TankPurchaseCostAlgorithm, vessel_material_factors
from thermosteam import Stream, MultiStream
from orgacids.utils import baseline_feedflow

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction


# %% Feedstock preprocessing

# =============================================================================
# The system as a whole, capital and operating costs already considered in 
# the cost of feedstock cost
# =============================================================================
@cost(basis='Flow rate', ID='System', units='kg/hr',
      kW=511.3205, cost=13329690, S=94697, CE=521.9, n=0.6, BM=1.7)
class FeedstockPreprocessing(Unit):
    _N_ins = 1
    _N_outs = 1
    
    _baseline_feedflow = baseline_feedflow.copy()
    
    # U.S. ton/day, 2000 metric tonne/day as in Humbird et al.
    feedstock_flow_rate = 2205
    
    def _run(self):
        feed_in = self.ins[0]
        feed_out = self.outs[0]
        
        feed_in.mass = self._baseline_feedflow * (self.feedstock_flow_rate/2205)
        feed_out.copy_like(feed_in)


# %% Pretreatment

# =============================================================================
# Sulfuric acid in-line mixer
# =============================================================================
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      cost=6000, S=136260, CE=521.9, n=0.5, BM=1)
class SulfuricAcidMixer(Unit):
    _N_ins = 2
    _N_outs = 1
    _graphics = Mixer._graphics
        
    def _run(self):
        acid, water = self.ins
        mixture = self.outs[0]
        # 0.05 is from 1842/36629 from streams 710 and 516 of Humbird et al.
        water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
        mixture.mix_from([acid, water])

# =============================================================================
# Modified from bst.units.Pump, which won't simulate for 0 flow 
# =============================================================================
class OrganicAcidsPump(Pump):
    def _design(self):
        if self.ins[0].F_mol == 0:
            Design = self.design_results
            Design['Ideal power'] = 0
            Design['Flow rate'] = 0
            Design['Efficiency'] = 0
            Design['Actual power'] = 0
            Design['Pump power'] = 0
            Design['N'] = 0
            Design['Head'] = 0
            Design['Type'] = 'NA'
        else: Pump._design(self)
      
    def _cost(self):
        if self.ins[0].F_mol == 0:
            Cost = self.purchase_costs
            Cost['Pump'] = 0
            Cost['Motor'] = 0
        else: Pump._cost(self)

# =============================================================================
# Adjust pretreatment water loading, 30% from Table 5 on Page 21 of Humbird et al.
# =============================================================================
class PretreatmentMixer(Mixer):
    _N_ins = 3
    _N_outs = 1
    
    solid_loading = 0.3
        
    def _run(self):
        feedstock, acid, water = self.ins
        mixture_out = self.outs[0]
        
        mixture = feedstock.copy()
        mixture.mix_from([feedstock, acid])
        
        total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solid_loading
        water.imass['Water'] = total_mass - mixture.F_mass
        
        mixture_out.mix_from([mixture, water])

# =============================================================================
# Steam mixer
# =============================================================================
class SteamMixer(Unit):
    """
    Parameters
    ----------
    ins :
        [0] Feed        
        [1] Steam
    
    outs : 
        [0] Mixed steam    
        
    """
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), *, P):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        
    @staticmethod
    def _P_at_flow(mol_water, P, steam, mixed, feed):
        steam.imol['Water'] = mol_water
        mixed.mol = steam.mol + feed.mol
        mixed.H = feed.H + steam.H
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P - P_new
    
    def _run(self):
        feed, steam = self.ins
        mixed = self.outs[0]

        steam_mol = steam.F_mol
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam, mixed, feed))
        mixed.P = self.P
    
    @property
    def installation_cost(self): return 0
    
    @property
    def purchase_cost(self): return 0
    
    def _design(self): pass
    def _cost(self): pass

# =============================================================================
# Pretreatment reactor
# =============================================================================
@cost(basis='Dry flow rate', ID='Pretreatment Reactor', units='kg/hr',
      kW=4578, cost=19812400, S=83333, CE=521.9, n=0.6, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = Flash._graphics
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = MultiStream(None)
        vapor, liquid = self.outs
        vapor.phase = 'g'
        
        self.pretreatment_rxns = ParallelRxn([
            #            Reaction definition                 Reactant   Conversion
            # Below from Table 6 on Page 22 of Humbird et al.
            Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.099),
            Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.003),
            Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.003),
            Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1),
            Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9),
            Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.024),
            Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.05),
            Rxn('Acetate -> AceticAcid',                     'Acetate',  1),
            Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.05),
            # Below from Page 106 of Humbird et al.,
            Rxn('Mannan + H2O -> Mannose',                   'Mannan',   0.9),
            Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.024),
            Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.05),
            Rxn('Galactan + H2O -> Galactose',               'Galactan', 0.9),
            Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.024),
            Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.05),
            Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9),
            Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.024),
            Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.05),
            Rxn('Furfural -> Tar',                           'Furfural', 1),
            Rxn('HMF -> Tar',                                'HMF',      1)
            ])
    
    def _run(self):
        ms = self._multistream
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copy_like(feed)
        self.pretreatment_rxns(liquid.mol) 
        ms.copy_like(liquid)
        H = ms.H + ms.Hf - feed.Hf
        ms.vle(T=130+273.15, H=H)
        vapor.mol = ms.imol['g']
        liquid.mol = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P

# =============================================================================
# Modified from bst.units.Flash, add agitator and heat exchanger
# =============================================================================
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
class OrganicAcidsFlash(Flash):
    _units = {**Flash._units, 
              'Flow rate': 'kg/hr'}
    
    def _design(self):
        Flash._design(self)
        self.design_results['Flow rate'] = self.ins[0].F_mass
    
    def _end_decorated_cost_(self):
        rate = self.power_utility.rate
        Flash._cost(self)
        if self.P and self.P < 101325:
            self.power_utility += rate
        else:
            self.power_utility.rate = rate

# =============================================================================
# Ammonia in-line mixer
# =============================================================================
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      # Size basis on the total flow, not just ammonia, 
      # thus assuming difference caused by MWs of NH3 and NH4OH is negligible
      cost=5000, S=157478, CE=521.9, n=0.5, BM=1)
class AmmoniaMixer(Mixer):
    def _run(self):
        ammonia, water = self.ins
        mixture = self.outs[0]
        
        # Make 10% ammonia solution
        water.imass['Water'] = ammonia.imass['AmmoniumHydroxide'] * (1-0.1)/0.1
        mixture.mix_from([ammonia, water])

# =============================================================================
# Modified from bst.units.MixTank, which won't simulate for 0 flow 
# =============================================================================
class OrganicAcidsMixTank(MixTank):
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
        
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.purchase_costs['Tanks'] = 0
        else: MixTank._cost(self)

# =============================================================================
# Add ammonia to neutralize pretreatment slurry
# =============================================================================
class AmmoniaAdditionTank(OrganicAcidsMixTank):
        def _run(self):
            self.neutralization_rxn = Rxn(
                #          Reaction definition      Reactant Conversion
                '2 NH4OH + H2SO4 -> NH4SO4 + 2 H2O', 'H2SO4', 0.95)
            MixTank._run(self)        
            self.neutralization_rxn.adiabatic_reaction(self.outs[0])


# %% Saccharification and fermentation

# =============================================================================
# In-line enzyme hydrolysate mixer
# =============================================================================
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      kW=74.57, cost=109000, S=379938, CE=521.9, n=0.5, BM=1.7)
class EnzymeHydrolysateMixer(Mixer):
    _N_ins = 3
    _N_outs = 1
    
    # 20 mg/g glucan based on Table 18 of Humbird et al.
    enzyme_loading = 20
    solid_loading = 0.2

    def _run(self):
        hydrolysate, enzyme, water = self.ins
        effluent = self.outs[0]
        
        enzyme.imass['Enzyme'] = (self.enzyme_loading/1000) * hydrolysate.imass['Glucan']
        # 13054 and 13836 are from stream 422 in Humbird et al.
        enzyme.imass['Water'] = enzyme.imass['Enzyme'] * (13836-13054)/13836
        mixture = hydrolysate.copy()
        mixture.mix_from([hydrolysate, enzyme])
        
        total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solid_loading
        water.imass['Water'] = total_mass - mixture.F_mass
        
        effluent.mix_from([hydrolysate, enzyme, water])

# =============================================================================
# Saccharification and co-fermentation (both Glucose & Xylose are used in fermentation)
# Not including heat exchanger as saccharificatoin and co-fermentation 
# are at the same temperature now
# =============================================================================
@cost(basis='Saccharification tank size', ID='Prehydrolysis tanks', units='kg',
      # Size basis changed from 421776 as the original design is not sufficient
      cost=3840000, S=283875*24, CE=521.9, n=0.7, BM=2)
@cost(basis='Slurry flow rate', ID='Saccharification transfer pumps', units='kg/hr',
      kW=74.57, cost=47200, S=421776*24, CE=521.9, n=0.8, BM=2.3)
@cost(basis='Fermenter size', ID='Fermentors', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36), CE=521.9, n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Agitators', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      # and total residence time (batch hydrolysis and fermentation)
      kW=22.371, cost=52500, S=(42607+443391+948+116)*(60+36), CE=521.9, n=1, BM=1.5)
@cost(basis='Recirculation flow rate', ID='Recirculation pumps', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=521.9, n=0.8, BM=2.3)
class SaccharificationAndCoFermentation(Unit):    
    _N_ins = 3
    _N_outs = 2
    _units= {'Saccharification tank size': 'kg',
             'Slurry flow rate': 'kg/hr',
             'Fermenter size': 'kg',
             'Recirculation flow rate': 'kg/hr'}             

    # Same T for saccharificatoin and co-fermentation
    T = 50+273.15
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_saccharification = 24
    
    # Co-Fermentation time (hr)
    tau_cofermentation = 120
    
    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.07
    
    CSL_loading = 10 # kg/m3
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        # Based on Table 9 on Page 28 of Humbird et al.
        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.85),
            Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
            ])
        
        self.saccharified_stream = Stream(None)
        
        # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.02),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.07),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.02),
        ])
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
            #   Reaction definition                                             Reactant  Conversion
            Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O', 'LacticAcid',   0.95),
            Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O', 'AceticAcid',   0.95)
            ])

    def _run(self):
        feed, CSL, lime = self.ins
        effluent, sidedraw = self.outs
        ss = self.saccharified_stream
        ss.T = sidedraw.T = effluent.T = self.T
        
        ss.copy_like(feed)
        CSL.imass['CSL'] = feed.F_vol * self.CSL_loading 
        ss.mol += CSL.mol
        self.saccharification_rxns(ss.mol)
        # Sidedraw to SeedTrain
        sidedraw.mol = ss.mol * self.inoculum_ratio
        effluent.mol = ss.mol - sidedraw.mol
        
        self.cofermentation_rxns(effluent.mol)
        effluent.imass['CSL'] = 0
        
        # Set feed lime mol to match rate of acids production, add 5% extra
        lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
                            +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1]) \
                            * 1.05
        effluent.mol += lime.mol
        self.neutralization_rxns(effluent.mol)
    
    def _design(self):
        Design = self.design_results
        total_mass_flow = self.ins[0].F_mass + self.ins[1].F_mass
        Design['Saccharification tank size'] = total_mass_flow * self.tau_saccharification
        Design['Slurry flow rate'] = total_mass_flow
        Design['Fermenter size'] = self.outs[0].F_mass * self.tau_cofermentation
        Design['Recirculation flow rate'] = total_mass_flow

# =============================================================================
# Seed train, 5 stages, 2 trains
# assumes no need for pH control due to low acid concentration
# =============================================================================
@cost(basis='Seed fermenter size', ID='Stage #1 fermenters', units='kg',
      # 44339, 211, and 26 are streams 303, 309, and 310 in Humbird et al.
      cost=75400, S=(44339+211+26)*36, CE=521.9, n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #2 fermenters', units='kg',
      cost=116600, S=(44339+211+26)*36, CE=521.9, n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #3 fermenters', units='kg',
      cost=157600, S=(44339+211+26)*36, CE=521.9, n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #4 fermenters', units='kg',
      cost=352000, S=(44339+211+26)*36, CE=521.9, n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #4 agitators', units='kg',
      kW=11.1855, cost=26000, S=(44339+211+26)*36, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Seed fermenter size', ID='Stage #5 fermenters', units='kg',
      cost=1180000, S=(44339+211+26)*36, CE=521.9, n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #5 agitators', units='kg',
      kW=14.914, cost=43000, S=(44339+211+26)*36, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pumps', units='kg/hr',
      kW=59.656, cost=24300, S=43149, CE=521.9, n=0.8, BM=2.3)
class SeedTrain(Unit):
    _N_ins = 1
    _N_outs = 1
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}

    # Operating temperature (K)
    T = 50+273.15
    
    # Cycle time for each batch (hr), including 12 hr turnaround time 
    tau_batch = 36
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    def __init__(self, ID='', ins=None, outs=(), ferm_ratio=0.9):
        Unit.__init__(self, ID, ins, outs)
        self.ferm_ratio = ferm_ratio

        # FermMicrobe reaction from Table 11 on Page 29 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76*ferm_ratio),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07*ferm_ratio),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.04),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76*ferm_ratio),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.07*ferm_ratio),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.04),
        ])

    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(feed)

        self.cofermentation_rxns(effluent.mol)
        # Assume all CSL is used up
        effluent.imass['CSL'] = 0 
        
        effluent.T = self.T

    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass
        Design['Seed fermenter size'] = self.outs[0].F_mass * self.tau_batch


# %% Organic acid separation

# =============================================================================
# Filter to separate fermentation broth into products liquid and solid
# Scaling values in Humbird et al. was adjusted, as the original design used
# sludge flow rate for scaling, which was much smaller than the feed flow rate
# =============================================================================
@cost(basis='Feed flow rate', ID='Feed tank', units='kg/hr',
      cost=174800, S=31815, CE=550.8, n=0.7, BM=2.0)
@cost(basis='Feed flow rate', ID='Feed pump', units='kg/hr',
      kW=74.57, cost= 18173, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost(basis='Pressing air flow rate', ID='Filter pressing compressor', units='kg/hr',
      kW=111.855, cost=75200, S=808, CE=521.9, n=0.6, BM=1.6)
@cost(basis='Pressing air flow rate', ID='Pressing air compressor reciever', units='kg/hr',
      cost=8000, S=31815, CE=550.8, n=0.7, BM=3.1)
@cost(basis='Drying air flow rate', ID='Filter drying compressor', units='kg/hr',
      kW=1043.98, cost=405000, S=12233, CE=521.9, n=0.6, BM=1.6)
@cost(basis='Drying air flow rate', ID='Dry air compressor reciever', units='kg/hr',
      cost=17000, S=31815, CE=550.8, n=0.7, BM=3.1)
@cost(basis='Feed flow rate', ID='Pressure filter', units='kg/hr',
      cost=3294700, S=31815, CE=550.8, n=0.8, BM=1.7)
@cost(basis='Filtrate flow rate', ID='Filtrate discharge pump', units='kg/hr',
      # Power not specified, based on filtrate tank discharge pump
      kW=55.9275, cost=13040, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost(basis='Filtrate flow rate', ID='Filtrate tank', units='kg/hr',
      cost=103000, S=31815, CE=550.8, n=0.7, BM=2.0)
@cost(basis='Filtrate flow rate', ID='Flitrate tank agitator', units='kg/hr',
      kW=5.59275, cost=26000,  S=337439, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Filtrate flow rate', ID='Filtrate tank discharge pump', units='kg/hr',
      kW=55.9275, cost=13040, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost(basis='Solids flow rate', ID='Cell mass wet cake conveyor', units='kg/hr',
      kW=7.457, cost=70000, S=28630, CE=521.9, n=0.8, BM=1.7)
@cost(basis='Solids flow rate', ID='Cell mass wet cake screw',  units='kg/hr',
      kW=11.1855, cost=20000, S=28630, CE=521.9, n=0.8, BM=1.7)
@cost(basis='Feed flow rate', ID='Recycled water tank', units='kg/hr',
      cost=1520,  S=31815, CE=550.8, n=0.7, BM=3.0)
@cost(basis='Feed flow rate', ID='Manifold flush pump', units='kg/hr',
      kW=74.57, cost=17057, S=31815, CE=550.8, n=0.8, BM=2.3)
@cost(basis='Feed flow rate', ID='Cloth wash pump', units='kg/hr',
      kW=111.855,cost=29154, S=31815, CE=550.8, n=0.8, BM=2.3)
class CellMassFilter(SolidsSeparator):
    _N_ins = 1
    _units= {'Feed flow rate': 'kg/hr',
             'Pressing air flow rate': 'kg/hr',
             'Drying air flow rate': 'kg/hr',
             'Filtrate flow rate': 'kg/hr',
             'Solids flow rate': 'kg/hr'}
            
    def _design(self):
        Design = self.design_results
        Design['Feed flow rate'] = self.ins[0].F_mass
        # 809 is the scailng basis of equipment M-505,
        # 391501 from stream 508 in Humbird et al.
        Design['Pressing air flow rate'] = 809/391501 * self.ins[0].F_mass
        # 12105 and 391501 from streams 559 and 508 in Humbird et al.
        Design['Drying air flow rate'] = 12105/391501 * self.ins[0].F_mass
        Design['Solids flow rate'] = self.outs[0].F_mass
        Design['Filtrate flow rate'] = self.outs[1].F_mass

# =============================================================================
# A reactor class for acidulation, esterification, and hydrolysis
# =============================================================================
class BatchReactor(Unit, isabstract=True):
    '''    
    Create an abstract class for batch reactor unit, purchase cost of the reactor
    is based on volume (0.1-20 m^3 for each tank) calculated by residence time.

    Parameters
    ----------
    ins : stream
        Inlet.
        
    outs : stream
        Outlet.

    tau=1 : float
        Residence time [hr].
        
    V_wf=0.8 : float
        Fraction of working volume over total volume.
        
    kW_per_m3=0.0985: float
        Utility usage of agitator (from MixTank).
        
    vessel_material : str, optional
        Vessel material. Defaults to 'Stainless steel 316'.
        
    References
    ----------
    .. [1] Apostolakou, A. A., Kookos, I. K., Marazioti, C., Angelopoulos, K. C.
        (2009). Techno-economic analysis of a biodiesel production process from
        vegetable oils. Fuel Processing Technology, 90(7–8), 1023–1031.
        https://doi.org/10.1016/j.fuproc.2009.04.017
    '''
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(),
                 tau=1, V_wf=0.8, 
                 kW_per_m3=0.0985, # from MixTank
                 vessel_material='Stainless steel'):
        
        Unit.__init__(self, ID, ins, outs)
        self.tau = tau
        self.V_wf = V_wf
        self.kW_per_m3 = kW_per_m3
        self.vessel_material = vessel_material
        
    purchase_cost_algorithm = TankPurchaseCostAlgorithm(
        lambda V: 15000 * V ** 0.55,
        V_min=0.1, V_max=20, V_units='m^3',
        CE=525.4, material='Stainless steel')
        
    def _design(self):
        Design = self.design_results
        ins_F_vol = sum([i.F_vol for i in self.ins])        
        V = ins_F_vol * self.tau / self.V_wf
        
        Design['Residence time'] = self.tau
        Design['Total volume'] = V
    
    def _cost(self):
        Design = self.design_results
        purchase_costs = self.purchase_costs
        
        if Design['Total volume'] == 0:
            Design['Number of tanks'] = 0
            purchase_costs['Tanks'] = 0
        
        else:
            material_factor = vessel_material_factors[self.vessel_material]
            # The original cost algorith given for stainless steel
            material_factor /= vessel_material_factors['Stainless steel']

            N, Cp = compute_number_of_tanks_and_total_purchase_cost(
                self.design_results['Total volume'],
                self.purchase_cost_algorithm,
                material_factor)
            
            Design['Number of tanks'] = N
            purchase_costs['Tanks'] = Cp
        
        self.power_utility(self.kW_per_m3 * Design['Total volume'])

# =============================================================================
# Reactor to acidify CalciumLactate into lactic acid
# =============================================================================
class AcidulationReactor(BatchReactor):
    _N_ins = 2
    _N_outs = 1
    
    acidulation_rxns = ParallelRxn([
        #   Reaction definition                                    Reactant   Conversion
        Rxn('CalciumLactate + H2SO4 -> 2 LacticAcid + CaSO4', 'CalciumLactate', 1),
        Rxn('CalciumAcetate + H2SO4 -> 2 AceticAcid + CaSO4', 'CalciumAcetate', 1),
        Rxn('2 AmmoniumHydroxide + H2SO4 -> AmmoniumSulfate + 2 H2O', 'AmmoniumHydroxide', 1),
        Rxn('CalciumDihydroxide + H2SO4 -> CaSO4 + 2 H2O', 'CalciumDihydroxide', 1)
        ])
            
    def _run(self):
        feed, acid = self.ins
        effluent = self.outs[0]
        
        # Set feed acid mol to match acidulation needs with 5% extra
        acid.imol['H2SO4'] = (feed.imol['CalciumLactate']/self.acidulation_rxns.X[0] \
                             +feed.imol['CalciumAcetate']/self.acidulation_rxns.X[1] \
                             +0.5*feed.imol['AmmoniumHydroxide']/self.acidulation_rxns.X[2] \
                             +0.5*feed.imol['CalciumDihydroxide']/self.acidulation_rxns.X[3]) \
                             * 1.05
        effluent.mix_from([feed, acid])
        self.acidulation_rxns.adiabatic_reaction(effluent)

# =============================================================================
# Filter to separate gypsum from the acidified fermentation broth
# =============================================================================
@cost(basis='Feed flow rate', ID='Hydroclone & rotary drum filter', units='kg/hr',
      # Size based on stream 239 in Aden et al.,
      # no power as Centrifuge in AerobicDigestion in Humbird et al. has no power
      cost=187567, S=272342, CE=389.5, n=0.39, BM=1.4)
@cost(basis='Filtrate flow rate', ID='Filtered hydrolysate pump', units='kg/hr',
      # Size based on stream 230 in Aden et al.,
      # power based on SaccharificationAndCoFermentation Filtrate Saccharification transfer pumps
      kW=74.57*265125/421776, cost=31862, S=265125, CE=386.5, n=0.79, BM=2.8)
class GypsumFilter(SolidsSeparator):
    _N_ins = 1
    _units = {'Feed flow rate': 'kg/hr',
              'Filtrate flow rate': 'kg/hr'}
    
    def _design(self):
        Design = self.design_results
        Design['Feed flow rate'] = self.ins[0].F_mass
        Design['Filtrate flow rate'] = self.outs[1].F_mass

# # =============================================================================
# # Heat exchanger that uses utility to completely condense vapor into liquid,
# # compare to the originial BioSTEAM HXutility, this one defaults rigours=True,
# # and only have one liquid outs (a MultiStream as in HXutility can causes problem
# # in subsequent units)
# # =============================================================================
# class VaporCondenser(HXutility):
#     def _run(self):
#         self.rigorous = True
#         HXutility._run(self)
#         liquid = Stream(phase='l')
#         liquid.mol = self.outs[0].imol['l']
#         liquid.T = self.outs[0].T
#         self.outs[0] = liquid

# =============================================================================
# Reactor to esterify lactic acid with ethanol for subsequent separation
# =============================================================================
#!!! Change to Reactor class
# Cost copied from Transesterification
@cost('Volume', 'Reactor',
      CE=525.4, cost=15000, n=0.55, kW=1.5, BM=4.3,)
# @cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
#       cost=85000, S=-8, CE=550.8, n=0.7, BM=2.2)
class Esterification(Unit):
    """
    Create an esterification reactor that converts 'LacticAcid', 'AceticAcid', and 'Ethanol'
    to 'EthylLactate', 'EthylAcetate', and 'Water'. Finds the amount of catalyst 'Amberlyst-15'
    required and consumes it.
    
    Parameters
    ----------
    ins : stream sequence
        [0] Organic acids feed
        [1] Ethanol feed (includes catalyst)
        
    outs : stream
        Reactor effluent.
        
    efficiency : float
        Efficiency of conversion (on a 'LacticAcid' basis).
        
    ethanol2LA : float
        Ethanol feed to LacticAcid molar ratio.
        
    T : float
        Operating temperature [K].
    
    """
    _bounds = {'Volume': (0.1, 20)}
    _units = {'Volume': 'm^3'}
    _tau = 1
    _N_ins = 5
    _N_outs = 2
    _N_heat_utilities = 1

    def _more_design_specs(self):
        return (('Residence time', self._tau, 'hr'),
                ('Conversion efficiency', self.X1, ''),
                ('Working volume fraction', 0.8, ''))
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 X1=None, tau=None, X2=None, 
                 ethanol2LA=1.2,
                 cat_load=0.039, 
                 T=363.15, 
                 tau_max=300, 
                 assumeX2equalsX1 = True, 
                 allow_higher_T=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reaction
        

        # self._methanol_composition = chemicals.kwarray(
        #         dict(Methanol=1-catalyst_molfrac,
        #              NaOCH3=catalyst_molfrac))
        self.reactives = \
                ('LacticAcid', 'Ethanol', 'EthylLactate', 'H2O', 'AceticAcid', 'EthylAcetate')
        self.assumeX2equalsX1 = assumeX2equalsX1
        self.X1, self.X2, self._tau = X1, X2, tau
        
        self.T = T #: Operating temperature (K). # 353.40 K
        self.allow_higher_T = allow_higher_T
        self.ethanol2LA = ethanol2LA # 1.134
        
        self.tau_max = tau_max
        
        self.cat_load = cat_load # 2.4%
          # kg
        
        self.K = exp(2.9625 - 515.13/T)
        self.kc = 2.70 * (1e8) * exp(-6011.55/T)
        self.KW = 15.19 * exp(12.01/T)
        self.KEt = 1.22 * exp(359.63/T)
    
    def compute_X1_and_tau(self, time_step=1): # time_step in min
        lle_chemicals = self.chemicals.lle_chemicals
        reactives = self.reactives
        # LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
        #     lle_chemicals.get_index(reactives())
        K = self.K
        KW = self.KW
        KEt = self.KEt
        kc = self.kc
        T = self.T
        cat_load = self.cat_load
        ethanol2LA = self.ethanol2LA
        feed, ethanol1, recycled_LA, supplement_ethanol, ethanol2 = self.ins
        
        tau_max = self.tau_max
        # ethanol = ethanol1.copy()
        # ethanol.mix_from([ethanol,ethanol2])
        ethanol_reqd = self.ethanol_reqd

        # lle_chemicals = self.ins[0].lle_chemicals
        IDs = tuple([i.ID for i in lle_chemicals])
        
        f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
        
        
        temp_flow = feed.copy() 
        # recycle_LA = recycled_LA.copy()
        # temp_flow.mix_from([temp_flow, ])
        
        # temp_flow.imol['Ethanol'] = ethanol_reqd
        # temp_flow.mol += self.ins[1].mol.copy()
        ethanol = self.ethanol
        temp_flow.mix_from([temp_flow, recycled_LA, ethanol, supplement_ethanol])
        # temp_flow.imol['Ethanol'] = ethanol_reqd
        #mcat = cat_load * (self.ins[0].F_mass + ethanol_reqd * self.ins[0].chemicals.LacticAcid.MW)
        mcat = cat_load * temp_flow.F_mass
        self.mcat = mcat
        tau = time_step
        LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
            [IDs.index(ID) for ID in reactives]
            
                
        gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
        gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
            
        curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
        normalized_mol = temp_flow.get_normalized_mol(IDs)
        curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
        activities = gammas_reactives*curr_conc
        # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
        r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
        dX = r*time_step*mcat/(1000) # r is in mol g-1 min-1
        
        new_flows = [1, 1, 1, 1]
        LA_initial = temp_flow.imol['LacticAcid']
        

        while dX/LA_initial>1e-4:
            
            
            if curr_flow[0]<dX or curr_flow[1]<dX:
                dX = min(curr_flow[0], curr_flow[1])
                
            new_flows = [curr_flow[0] - dX, curr_flow[1] - dX, curr_flow[2] + dX, curr_flow[3] + dX]
            temp_flow.set_flow(new_flows, 'kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
            
            if new_flows[0]<=0 or new_flows[1]<=0 or tau>tau_max-time_step: # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
                # dX = min(new_flows)
                break
            
            curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
            normalized_mol = temp_flow.get_normalized_mol(IDs)
            curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
            gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
            gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
            activities = gammas_reactives*curr_conc
            
            # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
            r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
            # print(LA_initial)
            dX = r*time_step*mcat/(1000) # r is in mol g-1 min-1
            
            tau += time_step
            
        X1 = (feed.imol['LacticAcid'] + recycled_LA.imol['LacticAcid']- temp_flow.imol['LacticAcid']) / (feed.imol['LacticAcid']+recycled_LA.imol['LacticAcid'])
        self._tau = tau
        self.X1 = X1
        return X1, tau

    @property
    def tau(self):
        """Residence time (hr)."""
        return self._tau
    # @tau.setter
    # def tau(self, tau):
    #     self._tau = tau
    
    @property
    def efficiency(self):
        """Esterification conversion efficiency."""
        return self.reaction.X[0]


    def _run(self):
        feed, ethanol1, recycled_LA, supplement_ethanol, ethanol2 = self.ins # broth, recycled ethanol, recycled LA, supplementary ethanol
        product, wastewater = self.outs
        ethanol = ethanol1.copy()
        ethanol.mix_from([ethanol1, ethanol2])
        
        
        self.ethanol_reqd = ethanol_reqd = self.ethanol2LA * (feed.imol['LacticAcid'] + recycled_LA.imol['LacticAcid'])
        self.ethanol_difference = ethanol_difference = ethanol_reqd - ethanol.imol['Ethanol']
        
        
        self.added_ethanol = max(0,ethanol_difference)
        
        if self.ethanol_difference<0:
            self.ethanol_ratio = ethanol_ratio = ethanol_reqd/ethanol.imol['Ethanol']
            wastewater.mol = (1-ethanol_ratio)*ethanol.mol
            ethanol.mol = ethanol.mol - wastewater.mol
            supplement_ethanol.empty()
        else:
            wastewater.empty()
            supplement_ethanol.imol['Ethanol'] = self.added_ethanol
        
        self.ethanol = ethanol
        self.wastewater = wastewater
        
        # product.imol['Ethanol'] += self.added_ethanol
        # product.mol[:] = feed.mol + recycled_LA.mol + ethanol.mol + supplement_ethanol.mol
        product.mix_from([feed, recycled_LA, ethanol, supplement_ethanol])
        
        if self.allow_higher_T and self.T<product.T:
            self.T = product.T
            
        if self.X1==None:
            X1, tau = self.compute_X1_and_tau()
        else:
            X1, tau = self.X1, self._tau
        if self.assumeX2equalsX1:
            X2 = X1
            self.X2 = X2
        else:
            X2 = self.X2
            
        self.reaction = ParallelRxn([
            #   Reaction definition                               Reactant  Conversion
            Rxn('LacticAcid + Ethanol -> EthylLactate + H2O', 'LacticAcid', X1),
            Rxn('AceticAcid + Ethanol -> EthylAcetate + H2O', 'AceticAcid', X2)
            ])
        
        self.reaction(product.mol)
        product.T = self.T
        
    def _design(self):
        product = self._outs[0]
        self.design_results['Volume'] = self._tau * product.F_vol / 0.8
        self.heat_utilities[0](self.Hnet, product.T)

# =============================================================================
# Reactor to hydrolyze esters to acids
# =============================================================================
class HydrolysisReactor(BatchReactor):
    _N_ins = 4
    _N_outs = 2
    
    water2EtLA = 12
    
    hydrolysis_rxns = ParallelRxn([
        #              Reaction definition                  Reactant   Conversion
        Rxn('EthylLactate + H2O -> LacticAcid + Ethanol', 'EthylLactate', 0.8),
        Rxn('EthylAcetate + H2O -> AceticAcid + Ethanol', 'EthylAcetate', 0.8)
        ])
    
    def _run(self):
        feed, water, recycle1, recycle2 = self.ins
        effluent, wastewater = self.outs
        water2EtLA = self.water2EtLA
        ester_feed = feed.imol['EthylLactate', 'EthylAcetate'].sum()
        # ester_in = sum([i.imol['EthylLactate', 'EthylAcetate'].sum() for i in self.ins])
        
        recycles = recycle1.copy()
        recycles.mix_from([recycle1, recycle2])
        ester_recycles = recycles.imol['EthylLactate', 'EthylAcetate'].sum()
        
        ester_in = ester_feed + ester_recycles
        water_needed = ester_in * water2EtLA
        
        # water_difference  = water_needed - effluent.imol['Water'] + recycles.imol['Water']
        # We haven't simulate effluent...
        
        water_difference  = water_needed - (feed.imol['Water']+recycles.imol['Water'])
        
        # Influents have more water for reaction, discharge some recycles

        # # Modified calculation to maintain the water2EtLA ratio
        # if water_difference < 0 and recycles.imol['Water'] != 0:
        #     recycle_ratio = (water2EtLA*ester_feed-feed.imol['Water']) / \
        #         (water2EtLA*ester_recycles-recycles.imol['Water'])
        #     wastewater.mol = (1-recycle_ratio) * recycles.mol     
        #     water.empty()

        #!!! Why? To maintain the ratio? Yalin doesn't understand the calculation below
        if water_difference<0 and not recycles.imol['Water'] == 0:
            water_ratio = -(water_difference/recycles.imol['Water'])
            wastewater.mol = water_ratio * recycles.mol            
            recycles.mol -= wastewater.mol
            water.empty()
            
        # Influents don't have enough water, use fresh water
        else:
            water.imol['Water'] = water_difference
            wastewater.empty()
        
        effluent.mix_from([feed, recycles, water])
        # effluent.mol -= wastewater.mol

        wastewater.T = effluent.T        
        
        # effluent.mix_from([effluent, water])
        # effluent.imol['Water'] -= wastewater.imol['Water']
        
        self.hydrolysis_rxns(effluent.mol)


# %% Wastewater treatment

# =============================================================================
# Scrubber to treat waste vapors
# =============================================================================
@cost(basis='Vent flow rate', ID='Scrubber', units='kg/hr',
      cost=215000, S=22608, CE=521.9, n=0.6, BM=2.4)
class VentScrubber(Unit): 
    _N_ins = 2
    _N_outs = 2
    _units= {'Vent flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), *, gas=None):
        Unit.__init__(self, ID, ins, outs)
        if not gas:
            self.gas = sum([(i.ID,) for i in self.ins[0].chemicals 
                            if i.locked_state =='g'], ())
        else:
            self.gas = gas
    
    def _run(self):
        water, vent_in = self.ins
        vent_out, bottom = self.outs

        # Assume outs at Water.Tb (vent_in is predominantly Water)
        vent_out.T = bottom.T = 373.15
        vent_out.copy_flow(vent_in)
        bottom.copy_flow(vent_out, self.gas, remove=True, exclude=True)
        H = vent_in.H - vent_out.H
        H_per_mol = Stream('unit_stream', Water=1, T=373.15).H
        water.imol['Water'] = H / H_per_mol
        bottom.mix_from([water, bottom])
        
    def _design(self):
        self.design_results['Vent flow rate'] = self.ins[1].F_mass

# =============================================================================
# Total cost of wastewater treatment is combined into this placeholder,
# this is because all costs are scaled based on the total flow into WWT
# instead of individual units
# =============================================================================
@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=550.8, n=0.6, BM=1)
class WastewaterSystemCost(Unit): pass

# =============================================================================
# Anaerobic digestion
# =============================================================================
class AnaerobicDigestion(Unit):
    """
    Anaerobic digestion system as modeled by Humbird 2011
    
    Parameters
    ----------  
    ins :    
        [0] Wastewater
        
    outs :   
        [0] Biogas        
        [1] Treated water        
        [2] Sludge
        
    digestion_rxns: [ReactionSet] Anaerobic digestion reactions.  
      
    sludge_split: [Array] Split between wastewater and sludge
    
    """
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 1
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants, split=(), T=35+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.reactants = reactants
        self.split = split
        self.multi_stream = MultiStream(None)
        self.T = T
        chems = self.ins[0].chemicals
        
        # Based on P49 in Humbird et al., 91% of organic components is destroyed,
        # of which 86% is converted to biogas and 5% is converted to sludge,
        # and the biogas is assumed to be 51% CH4 and 49% CO2 on a dry molar basis
        biogas_MW = 0.51*chems.CH4.MW + 0.49*chems.CO2.MW
        f_CH4 = 0.51 * 0.86/0.91/biogas_MW
        f_CO2 = 0.49 * 0.86/0.91/biogas_MW
        f_sludge = 0.05 * 1/0.91/chems.WWTsludge.MW
        
        def anaerobic_rxn(reactant):
            MW = getattr(chems, reactant).MW
            return Rxn(f'{1/MW}{reactant} -> {f_CH4}CH4 + {f_CO2}CO2 + {f_sludge}WWTsludge',
                       reactant, 0.91)
        self.digestion_rxns = ParallelRxn([anaerobic_rxn(i) for i in self.reactants])
                
    def _run(self):
        wastewater = self.ins[0]
        biogas, treated_water, sludge = self.outs
        hu = self.heat_utilities[0]
        T = self.T
        
        # Calculate utility needs to keep digester temperature at 35°C,
        # heat change during reaction is not tracked
        H_at_35C = wastewater.thermo.mixture.H(mol=wastewater.mol, 
                                               phase='l', T=T, P=101325)
        duty = -(wastewater.H - H_at_35C)
        hu(duty, T)

        sludge.copy_flow(wastewater)
        self.digestion_rxns(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, T=T)
        biogas.mol = self.multi_stream.imol['g']
        biogas.phase = 'g'
        liquid_mol = self.multi_stream.imol['l']
        treated_water.mol = liquid_mol * self.split
        sludge.mol = liquid_mol - treated_water.mol
        biogas.receive_vent(treated_water, accumulate=True)
        biogas.T = treated_water.T = sludge.T = T

    
# =============================================================================
# Aerobic digestion
# =============================================================================
class AerobicDigestion(Unit):
    """
    Anaerobic digestion system as modeled by Humbird 2011
    
    Parameters
    ----------
    ins :  
        [0] Wastewater        
        [1] Air
        [2] Caustic, added to neutralize the nitric acid produced by 
            nitrifying bacteria duing nitrification process
        
    outs :    
        [0] Vent
        [1] Treated wastewater
        
    digestion_rxns : [ReactionSet] Anaerobic digestion reactions.   
     
    sludge_split : [Array] Split between wastewater and sludge
    
    """    
    _N_ins = 3
    _N_outs = 2
    # 4350, 4379, 356069, 2252, 2151522, and 109089 are water flows from 
    # streams 622, 630, 611, 632, 621, and 616  in Humbird et al.
    evaporation = 4350/(4379+356069+2252+2151522+109089)
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants, ratio):
        Unit.__init__(self, ID, ins, outs)
        self.reactants = reactants
        self.ratio = ratio
        chems = self.ins[0].chemicals
        
        def growth(reactant):
            f = chems.WWTsludge.MW / getattr(chems, reactant).MW 
            return Rxn(f"{f}{reactant} -> WWTsludge", reactant, 1.)
        
        # Reactions from auto-populated combustion reactions.
        # Based on P49 in Humbird et al, 96% of remaining soluble organic matter 
        # is removed after aerobic digestion, of which 74% is converted to
        # water and CO2 and 22% to cell mass
        combustion_rxns = chems.get_combustion_reactions()
        
        self.digestion_rxns = ParallelRxn([i*0.74 + 0.22*growth(i.reactant)
                                           for i in combustion_rxns
                                           if (i.reactant in reactants)])
        self.digestion_rxns.X[:] = 0.96
        
        #                                      Reaction definition       Reactant Conversion
        self.neutralization_rxn = Rxn('H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O', 'H2SO4', 0.95)
    
    def _run(self):
        influent, air, caustic = self.ins
        vent, effluent = self.outs
        ratio = self.ratio
        vent.phase = 'g'

        # 51061 and 168162 from stream 630 in Humbird et al.
        air.imass['O2'] = 51061 * ratio
        air.imass['N2'] = 168162 * ratio
        # 2252 from stream 632 in Humbird et al
        caustic.imass['NaOH'] = 2252 * ratio
        caustic.imol['NaOH'] += 2 * influent.imol['H2SO4'] / self.neutralization_rxn.X
        caustic.imass['H2O'] = caustic.imass['NaOH']
        effluent.copy_like(influent)
        effluent.mol += air.mol
        effluent.mol += caustic.mol
        self.neutralization_rxn(effluent.mol)
        self.digestion_rxns(effluent.mol)
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        vent.imol['Water'] = effluent.imol['Water'] * self.evaporation
        effluent.imol['Water'] -= vent.imol['Water']
        # Assume NaOH is completely consumed
        effluent.imol['NaOH'] = 0


# %% Facility units in biosteam

# =============================================================================
# Modified from bst.units.StorageTank, which won't simulate for 0 flow 
# =============================================================================
class OrganicAcidsStorageTank(StorageTank):
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.purchase_costs['Tanks'] = 0
        else: StorageTank._cost(self)

# =============================================================================
# For storage of lime used in separation and waste treatment,
# design copied from cornstover biorefinery in Aden et al.
# Base flow from stream 27 of Aden et al.
# =============================================================================
@cost(basis='Flow rate', ID='Storage bin', units='kg/hr',
      cost=136370, S=2395, CE=386.5, n=0.46, BM=1.3)
# Cost not scaled, thus used n=0
# Power usage scaled based on M104 in Humbird et al., 
# (truck dumper hopper for feedstock handling) 
@cost(basis='Flow rate', ID='Feeder', units='kg/hr',
      kW=37.285/3500*140, cost=3900, S=2395, CE=386.5, n=0, BM=1.3)
# Power usage scaled based on M-106 in Humbird et al.
# (dust collection system for feedstock handling) 
@cost(basis='Flow rate', ID='Unloading blower', units='kg/hr',
      kW=18.6425*7425/8500, cost=99594, S=2395, CE=389.5, n=0.5, BM=1.4)
@cost(basis='Flow rate', ID='Dust vent baghouse', units='kg/hr',
      cost=140707, S=2395, CE=386.5, n=1, BM=1.5)
class LimeStorageBin(Unit): pass


# %% Superseded units

'''Blowdown tank'''
# =============================================================================
# # Costs of Tank and Agitator included in the Pump
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=0, S=264116, CE=521.9, n=0.7, BM=2)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=170, cost=0, S=252891, CE=521.9, n=0.5, BM=1.5)
# @cost(basis='Flow rate', ID='Pump3', units='kg/hr',
#       kW=93.2125, cost=25365, S=292407, CE=550.8, n=0.8, BM=2.3)
# class BlowdownTank(Unit): pass
# =============================================================================

'''Oligomer conversion tank'''
# =============================================================================
# # Copied PretreatmentFlash as the original design is not sufficient
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
# @cost(basis='Flow rate', ID='Pump4', units='kg/hr',
#       kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
# class OligomerConversionTank2(Unit): pass
# =============================================================================

'''Waste vapor condenser'''
# =============================================================================
# @cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
#       cost=34000, S=-2, CE=521.9, n=0.7, BM=2.2)
# class WasteVaporCondenser(HXutility):
#     _graphics = HXutility._graphics
# =============================================================================

'''Hydrolysate cooler'''
# =============================================================================
# @cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
#       cost=85000, S=-8, CE=550.8, n=0.7, BM=2.2)
# class HydrolysateCooler(HXutility):
#     _graphics = HXutility._graphics
# =============================================================================

'''Ammonia addition tank'''
# =============================================================================
# # Size basis on the total flow, not just ammonia, 
# # thus assuming size basis difference caused by MWs of NH3 and NH4OH is negligible,
# # pumping is provided by a separate HydrolysatePump
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=236000, S=410369, CE=521.9, n=0.7, BM=2)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=7.457, cost=21900, S=410369, CE=521.9, n=0.5, BM=1.5)
# class AmmoniaAdditionTank2(Unit): 
#     _N_ins = 1
#     _N_outs = 1
#     
#     def __init__(self, ID='', ins=None, outs=()):
#         Unit.__init__(self, ID, ins, outs)
# 
#         #                                      Reaction definition      Reactant Conversion
#         self.neutralization_rxn = Rxn('2 NH4OH + H2SO4 -> NH4SO4 + 2 H2O', 'H2SO4', 0.95)
#     
#     def _run(self):
#         ins = self.ins[0]
#         outs = self.outs[0]
#         outs.copy_like(ins)        
#         self.neutralization_rxn(outs.mol)
# =============================================================================

'''Hydrolysate pump'''
# =============================================================================
# # Spelled as "Hydrolyzate" in Humbird et al.
# @cost(basis='Flow rate', ID='Pump6', units='kg/hr',
#       kW=74.57, cost=22500, S=402194, CE=521.9, n=0.8, BM=2.3)
# class HydrolysatePump(Unit):
#     _units = {'Flow rate': 'kg/hr'}
#     _graphics = Pump._graphics
#     
#     def _design(self):
#         Design = self.design_results
#         Design['Flow rate'] = self.outs[0].F_mass
# =============================================================================

'''Seed hold tank'''
# =============================================================================
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=439000, S=40414, CE=521.9, n=0.7, BM=1.8)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=11.1855, cost=31800, S=40414, CE=521.9, n=0.5, BM=1.5)
# @cost(basis='Flow rate', ID='Pump7', units='kg/hr',
#       kW=7.457, cost=8200, S=43149, CE=521.9, n=0.8, BM=2.3)
# class SeedHoldTank(Unit): pass
# =============================================================================

'''Sulfuric acid storage tank'''
# =============================================================================
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=96000, S=1981, CE=550.8, n=0.7, BM=1.5)
# @cost(basis='Flow rate', ID='Pump11', units='kg/hr',
#       # Size basis changed from 1981 as the original design is not sufficient
#       kW=0.37285, cost=7493, S=1136, CE=550.8, n=0.8, BM=2.3)
# class SulfuricAcidStorageTank(StorageTank): pass
# =============================================================================

'''Ammonia storage tank'''
# =============================================================================
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       # Original size basis for NH3 instead of NH4OH
#       cost=196000, S=1171/17.031*35.046, CE=550.8, n=0.7, BM=2)
# # Design in Humbird et al. has no pump, this one copied from SulfuricAcidStorageTank
# @cost(basis='Flow rate', ID='Pump12', units='kg/hr',
#       # Size basis changed from 1981 as the original design is not sufficient
#       kW=0.37285, cost=7493, S=1136, CE=550.8, n=0.8, BM=2.3)
# class AmmoniaStorageTank(StorageTank): pass
# =============================================================================

'''CSL storage tank'''
# =============================================================================
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=70000, S=1393, CE=521.9, n=0.7, BM=2.6)
# @cost(basis='Flow rate', ID='Agitator', units='kg/hr',
#       kW=7.457, cost=21200, S=1393, CE=521.9, n=0.5, BM=1.5)
# @cost(basis='Flow rate', ID='Pump13', units='kg/hr',
#       kW=0.37285, cost=3000, S=1393, CE=521.9, n=0.8, BM=3.1)
# class CSLstorageTank(Unit): pass
# =============================================================================

'''Ethanol storage tank'''
# =============================================================================
# # Actually not very good to scale down from a big tank/pump
# # to a smaller one, but as the flow can be 0 during simulation, it will cause
# # trouble in biosteam's tank/pump cost calculation (divided by 0)
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=1340000, S=22681, CE=521.9, n=0.7, BM=1.7)
# @cost(basis='Flow rate', ID='Pump14', units='kg/hr',
#       kW=7.457, cost=9200, S=22681, CE=521.9, n=0.8, BM=3.1)
# class EthanolStorageTank(Unit): pass
# =============================================================================

'''Fire water tank'''
# =============================================================================
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=803000, S=8343, CE=521.9, n=0.7, BM=1.7)
# @cost(basis='Flow rate', ID='Pump14', units='kg/hr',
#       kW=93.2125, cost=15000, S=8343, CE=521.9, n=0.8, BM=3.1)
# class FireWaterTank(Unit): pass
# =============================================================================


# %% 

'''Hydrolysis'''    
# =============================================================================
# class Hydrolysis(Unit):
#     """
#     Create an esterification reactor that converts 'LacticAcid', 'AceticAcid', and 'Ethanol'
#     to 'EthylLactate', 'EthylAcetate', and 'Water'. Finds the amount of catalyst 'Amberlyst-15'
#     required and consumes it.
#     
#     Parameters
#     ----------
#     ins : stream sequence
#         * [0] Organic acids feed
#         * [1] Ethanol feed (includes catalyst)
#     outs : stream
#         Reactor effluent.
#     efficiency : float
#         Efficiency of conversion (on a 'LacticAcid' basis).
#     ethanol2LA : float
#         Ethanol feed to LacticAcid molar ratio.
#     T : float
#         Operating temperature [K].
#     
#     """
#     _bounds = {'Volume': (0.1, 20)}
#     _units = {'Volume': 'm^3'}
#     _tau = 1
#     _N_ins = 2
#     _N_outs = 1
#     _N_heat_utilities = 1
# 
#     def _more_design_specs(self):
#         return (('Residence time', self._tau, 'hr'),
#                 ('Conversion efficiency', self.X1, ''),
#                 ('Working volume fraction', 0.8, ''))
# 
#     def __init__(self, ID='', ins=None, outs=(), thermo = None, *, X1 = None, tau = None, X2 = None, cat_load=0.024, T=353.40, assumeX2equalsX1 = True):
#         Unit.__init__(self, ID, ins, outs, thermo)
#         #: [:class:`~thermosteam.ParallelReaction`] Transesterification and catalyst consumption reaction
#         
# 
#         # self._methanol_composition = chemicals.kwarray(
#         #         dict(Methanol=1-catalyst_molfrac,
#         #              NaOCH3=catalyst_molfrac))
#         self.reactives = \
#                 ('LacticAcid', 'Ethanol', 'EthylLactate', 'H2O', 'AceticAcid', 'EthylAcetate')
#         self.assumeX2equalsX1 = assumeX2equalsX1
#         self.X1, self.X2, self._tau = X1, X2, tau
#         
#         self.T = T #: Operating temperature (K). # 353.40 K
#         
#         # self.ethanol2LA = ethanol2LA # 1.134
#         # self.ethanol_reqd = 
#         
#         self.cat_load = cat_load # 2.4%
#           # kg
#         
#         self.K = exp(2.9625 - 515.13/T)
#         self.kc = 2.70 * (1e8) * exp(-6011.55/T)
#         self.KW = 15.19 * exp(12.01/T)
#         self.KEt = 1.22 * exp(359.63/T)
#     
#     def compute_X1_and_tau(self, time_step=450):
#         lle_chemicals = self.chemicals.lle_chemicals
#         reactives = self.reactives
#         # LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#         #     lle_chemicals.get_index(reactives())
#         K = self.K
#         KW = self.KW
#         KEt = self.KEt
#         kc = self.kc
#         T = self.T
#         cat_load = self.cat_load
#         # ethanol2LA = self.ethanol2LA
#         # self.ethanol_reqd = ethanol2LA * self.ins[0].imol['LacticAcid']
#         # ethanol_reqd = self.ethanol_reqd
#         mcat = cat_load * self.ins[0].F_mass
#         
#         
#         # lle_chemicals = self.ins[0].lle_chemicals
#         IDs = tuple([i.ID for i in lle_chemicals])
#         
#         f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
#         
#         
#         temp_flow = self.ins[0].copy()
#         # temp_flow.imol['Ethanol'] = ethanol_reqd
#         tau = time_step
#         LA_index, EtOH_index, EtLA_index, H2O_index, AA_index, EtAA_index = \
#             [IDs.index(ID) for ID in reactives]
#             
# 
#         gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#         gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
#             
#         curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#         normalized_mol = temp_flow.get_normalized_mol(IDs)
#         curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#         activities = gammas_reactives*curr_conc
#         # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#         r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#         dX = r*time_step*mcat/1000
#         
#         new_flows = [1, 1, 1, 1]
#         EtLA_initial = temp_flow.imol['EthylLactate']
#         
# 
#         while dX/EtLA_initial<-1e-5:
#             
# 
#             new_flows = [curr_flow[0] - dX, curr_flow[1] - dX, curr_flow[2] + dX, curr_flow[3] + dX]
#             if new_flows[0]<=0 or new_flows[1]<=0 or tau>7200-time_step: # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
#                 break
#             print(new_flows)
#             temp_flow.set_flow(new_flows, 'kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#             
#             
#             activities = gammas_reactives*curr_conc
#             # r = kc* (activities[EtOH_index]*activities[LA_index] - (activities[EtLA_index]*activities[H2O_index]/K))/(1 + KEt*activities[EtOH_index] + KW*activities[H2O_index])**2
#             r = kc* (activities[1]*activities[0] - (activities[3]*activities[2]/K))/(1 + KEt*activities[3] + KW*activities[2])**2
#             # print(LA_initial)
#             dX = r*time_step*mcat/1000
#             curr_flow = temp_flow.get_flow('kmol/hr', ('LacticAcid', 'Ethanol', 'Water', 'EthylLactate'))
#             normalized_mol = temp_flow.get_normalized_mol(IDs)
#             curr_conc = np.array([normalized_mol[LA_index], normalized_mol[EtOH_index], normalized_mol[H2O_index], normalized_mol[EtLA_index]])
#             gammas = f_gamma(temp_flow.get_normalized_mol(IDs), T)
#             gammas_reactives = np.array([gammas[LA_index], gammas[EtOH_index], gammas[H2O_index], gammas[EtLA_index]])
#             tau += time_step
#             
#         X1 = (self.ins[0].imol['EthylLactate'] - temp_flow.imol['EthylLactate']) / self.ins[0].imol['EthylLactate']
#         self._tau = tau
#         self.X1 = X1
#         return X1, tau
# 
#     @property
#     def tau(self):
#         """Residence time (hr)."""
#         return self._tau
# 
#     
#     # @property
#     def efficiency(self):
#         """Esterification conversion efficiency."""
#         return self.reaction.X[0]
# 
# 
#     def _run(self):
#         if self.X1==None:
#             X1, tau = self.compute_X1_and_tau()
#         else:
#             X1, tau = self.X1, self._tau
#         if self.assumeX2equalsX1:
#             X2 = X1
#             self.X2 = X2
#         else:
#             X2 = self.X2
#             
#         self.reaction = ParallelRxn([
#             #   Reaction definition                               Reactant  Conversion
#             Rxn('EthylLactate + H2O -> LacticAcid + Ethanol', 'EthylLactate', X1),
#             Rxn('EthylAcetate + H2O -> AceticAcid + Ethanol', 'EthylAcetate', X2)
#             ])
#         feed, ethanol = self.ins
#         product, = self.outs
#         # ethanol.imol['Ethanol'] = self.ethanol_reqd
#         
#         product.mol[:] = feed.mol + ethanol.mol
#         self.reaction(product.mol)
#         product.T = self.T
#         
#     def _design(self):
#         effluent = self._outs[0]
#         self.design_results['Volume'] = self._tau * effluent.F_vol / 0.8
#         self.heat_utilities[0](self.Hnet, effluent.T)
# =============================================================================

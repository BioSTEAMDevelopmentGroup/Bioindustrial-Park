#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

@author: sarangbhagwat
"""


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank, LiquidsSplitSettler
from biosteam.units.decorators import cost
from thermosteam import Stream, MultiStream
from HP.process_settings import price
from HP.utils import CEPCI, baseline_feedflow, compute_extra_chemical, adjust_recycle
from HP.chemicals_data import HP_chemicals
tmo.settings.set_thermo(HP_chemicals)
_kg_per_ton = 907.18474
_Gcal_2_kJ = 4.184 * 1e6 # (also MMkcal/hr)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction


# %% 

# =============================================================================
# Feedstock preprocessing
# =============================================================================

# The system as a whole, capital and operating costs already considered in 
# the cost of feedstock cost
@cost(basis='Flow rate', ID='System', units='kg/hr',
      kW=511.3205, cost=13329690, S=94697, CE=CEPCI[2009], n=0.6, BM=1.7)
class FeedstockPreprocessing(Unit):
    # 2205 U.S. ton/day (2000 metric tonne/day) as in Humbird et al.
    _baseline_flow_rate = baseline_feedflow.sum()
    _cached_flow_rate = 2205


# %% 

# =============================================================================
# Pretreatment
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=6210, S=1981, CE=CEPCI[2010],  n=0.7, BM=3)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      cost=8000, S=3720, CE=CEPCI[2009],  n=0.8, BM=2.3)
class SulfuricAcidAdditionTank(Unit):
    _N_ins = 1
    _N_outs = 1
    
    # Bseline is (18+4.1) mg/g dry biomass as in Humbird et al., 93% purity
    acid_loading = 22.1

    def __init__(self, ID='', ins=None, outs=(), *, feedstock_dry_mass):
        Unit.__init__(self, ID, ins, outs)
        self.feedstock_dry_mass = feedstock_dry_mass

    def _run(self):
        self.ins[0].imass['H2SO4'] = self.feedstock_dry_mass * (self.acid_loading/1000)
        self.ins[0].imass['H2O'] = self.ins[0].imass['H2SO4'] / 0.93 * 0.07
        self.outs[0].copy_like(self.ins[0])

# Sulfuric acid in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      cost=6000, S=136260, CE=CEPCI[2009], n=0.5, BM=1)
class SulfuricAcidMixer(Unit):
    _N_ins = 2
    _N_outs = 1
    _graphics = Mixer._graphics
        
    def _run(self):
        acid, water = self.ins
        mixture = self.outs[0]
        # 0.05 is from 1842/36629 from streams 710 and 516 of Humbird et al.
        water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
        mixture.mix_from([water, acid])

# Adjust pretreatment water loading, 30% from Table 5 on Page 21 of Humbird et al.
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

# Steam mixer
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
    def P_at_flow(mol_water, P, steam, mixed, feed):
        steam.imol['Water'] = mol_water
        mixed.mol = steam.mol + feed.mol
        mixed.H = feed.H + steam.H
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P_new-P
    
    def _run(self):
        feed, steam = self.ins
        mixed = self.outs[0]

        steam_mol = steam.F_mol
        steam_mol = max(0, aitken_secant(f=self.P_at_flow,
                                  x0=steam_mol, x1=steam_mol+0.1, 
                                  xtol=1e-4, ytol=1e-4,
                                  args=(self.P, steam, mixed, feed)))
        mixed.P = self.P
    
# Pretreatment reactor
@cost(basis='Dry flow rate', ID='Pretreatment reactor', units='kg/hr',
      kW=5120, cost=19812400, S=83333, CE=CEPCI[2009], n=0.6, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = Flash._graphics
    
    def __init__(self, ID='', ins=None, outs=(), T=130+273.15):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = MultiStream(None)
        self.T = T
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
        ms.vle(T=self.T, H=H)
        vapor.mol = ms.imol['g']
        liquid.mol = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P



# Blowdown tank, costs of Tank and Agitator included in the Pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=0, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=0, S=252891, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=25365, S=292407, CE=CEPCI[2010], n=0.8, BM=2.3)
class BlowdownTank(Unit): pass

# Oligomer conversion tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=203000, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=264116, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=17408, S=292407, CE=CEPCI[2009], n=0.8, BM=2.3)
class OligomerConversionTank(Unit): pass

# Pretreatment flash tank
@cost(basis='Liquid flow', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Liquid flow', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Liquid flow', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=CEPCI[2009], n=0.8, BM=2.3)
class PretreatmentFlash(Flash):
    _units= {'Liquid flow': 'kg/hr'}
    
    def _run(self):
        influent = self.ins[0]
        vapor, liquid = self.outs
        
        ms = tmo.MultiStream('ms')
        ms.copy_like(influent)
        ms.vle(P=101325, H=ms.H)
        
        vapor.mol = ms.imol['g']
        vapor.phase = 'g'
        liquid.mol = ms.imol['l']
        vapor.T = liquid.T = ms.T
        vapor.P = liquid.P = ms.P
        
    def _design(self):
        self.design_results['Liquid flow'] = self.outs[1].F_mass

# Ammonia in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
	  # Size basis on the total flow, not just ammonia, 	
      # thus assuming difference caused by MWs of NH3 and NH4OH is negligible
      cost=5000, S=157478, CE=CEPCI[2009], n=0.5, BM=1)
class AmmoniaMixer(Mixer):
    def _run(self):
        ammonia, water = self.ins
        mixture = self.outs[0]	
        	
        # Make 10% ammonia solution	
        water.imass['Water'] = ammonia.imass['AmmoniumHydroxide'] * (1-0.1)/0.1	
        mixture.mix_from([ammonia, water])

# Ammonia addition tank, size based on the total flow, not just ammonia, 
# thus assuming size difference caused by MWs of NH3 and NH4OH is negligible,
# pumping is provided by a separate HydrolysatePump unit
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=236000, S=410369, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21900, S=410369, CE=CEPCI[2009], n=0.5, BM=1.5)
class AmmoniaAdditionTank(MixTank):
    def _run(self):
        #                                      Reaction definition      Reactant Conversion
        self.neutralization_rxn = Rxn('2 NH4OH + H2SO4 -> NH4SO4 + 2 H2O', 'H2SO4', 1)
        
        MixTank._run(self)   
        self.neutralization_rxn.adiabatic_reaction(self.outs[0])

@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=74.57, cost=22500, S=402194, CE=CEPCI[2009], n=0.8, BM=2.3)
class HydrolysatePump(Unit):
    _units = {'Flow rate': 'kg/hr'}
    _graphics = Pump._graphics
    
    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass


# %% 

# =============================================================================
# Conversion
# =============================================================================

# Enzyme hydrolysate mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      kW=74.57, cost=109000, S=379938, CE=CEPCI[2009], n=0.5, BM=1.7)
class EnzymeHydrolysateMixer(Mixer):	
    _N_ins = 3	
    _N_outs = 1	
    	
    # 20 mg/g glucan based on Table 18 of Humbird et al.	
    enzyme_loading = 20	
    solid_loading = 0.2
    
    def _run(self):	
        hydrolysate, enzyme, water = self.ins	
        effluent = self.outs[0]	
        	
        # Make 10% extra based
        enzyme.imass['Enzyme'] = (self.enzyme_loading/1000*1.1) * hydrolysate.imass['Glucan']
        mixture = hydrolysate.copy()	
        mixture.mix_from([hydrolysate, enzyme])	
        	
        total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solid_loading	
        water.imass['Water'] = max(0, total_mass - mixture.F_mass)
        	
        effluent.mix_from([hydrolysate, enzyme, water])

# # Saccharification and co-fermentation (both Glucose & Xylose are used in fermentation)
# # Not including heat exchanger as saccharificatoin and co-fermentation 
# # are at the same temperature now
# # !!! But inlet is not at T=50+273.15; how is it brought to this temperature?

# @cost(basis='Saccharification tank size', ID='Saccharification tanks', units='kg',
#       cost=3840000, S=421776*24, CE=CEPCI[2009], n=0.7, BM=2)
# @cost(basis='Slurry flow rate', ID='Saccharification transfer pumps', units='kg/hr',
#       kW=74.57, cost=47200, S=421776*24, CE=CEPCI[2009], n=0.8, BM=2.3)
# @cost(basis='Fermenter size', ID='Fermentors', units='kg',
#       cost=10128000, S=(42607+443391+948+116)*(60+36), CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Fermenter size', ID='Agitators', units='kg',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in Humbird et al.)
#       # and total residence time (batch hydrolysis and fermentation)
#       kW=22.371, cost=52500, S=(42607+443391+948+116)*(60+36), CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Recirculation flow rate', ID='Recirculation pumps', units='kg/hr',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in Humbird et al.)
#       kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)

# class SaccharificationAndCoFermentation(Unit):    
#     _N_ins = 2
#     _N_outs = 2
#     _units= {'Saccharification tank size': 'kg',
#              'Slurry flow rate': 'kg/hr',
#              'Fermenter size': 'kg',
#              'Recirculation flow rate': 'kg/hr'}             
    
#     # Residence time of countinuous saccharification tanks (hr)
#     tau_saccharification = 24
    
#     # Co-Fermentation time (hr)
#     tau_cofermentation = 120
    
#     # Equals the split of saccharified slurry to seed train
#     inoculum_ratio = 0.07
    
#     CSL_loading = 10 # kg/m3
    
#     def __init__(self, ID='', ins=None, outs=(), T=50+273.15):
#         Unit.__init__(self, ID, ins, outs)
#         # Same T for saccharificatoin and co-fermentation
#         self.T = T
        
#         # Based on Table 9 on Page 28 of Humbird et al.
#         self.saccharification_rxns = ParallelRxn([
#             #   Reaction definition                   Reactant        Conversion
#             Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
#             Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
#             Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.85),
#             Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
#             ])
        
#         self.saccharified_stream = Stream(None)
        
#         # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
#         self.cofermentation_rxns = ParallelRxn([
#         #      Reaction definition            Reactant    Conversion
#         Rxn('Glucose -> HP',        'Glucose',   .9099),
#         # Rxn('Glucose -> 2,3-Butanediol',        'Glucose',   0.666),
#         Rxn('Glucose -> AceticAcid',               'Glucose',   0.07),
#         Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.02),
#         Rxn('Xylose -> HP',       'Xylose',    0.9099),
#         Rxn('Xylose -> AceticAcid',       'Xylose',    0.07),
#         Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.02),
#         ])
        
#         self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
#         self.xylose_to_HP_rxn = self.cofermentation_rxns[3]
        
#         # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
#         # self.neutralization_rxns = ParallelRxn([
#         # #   Reaction definition                                               Reactant  Conversion
#         # Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'LacticAcid',   1),
#         # Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1),
#         # Rxn('SuccinicAcid + CalciumDihydroxide -> CalciumSuccinate + 2H2O', 'SuccinicAcid', 1)
#         #     ])

#     def _run(self):
#         feed, CSL = self.ins
#         effluent, sidedraw = self.outs
#         ss = self.saccharified_stream
#         ss.T = sidedraw.T = effluent.T = self.T
        
#         ss.copy_like(feed)
#         CSL.imass['CSL'] = feed.F_vol * self.CSL_loading 
#         ss.mol += CSL.mol
#         self.saccharification_rxns(ss.mol)
#         # Sidedraw to SeedTrain
#         sidedraw.mol = ss.mol * self.inoculum_ratio
#         effluent.mol = ss.mol - sidedraw.mol
        
#         self.cofermentation_rxns(effluent.mol)
#         effluent.imass['CSL'] = 0
        
#         # Set feed lime mol to match rate of acids production, add 5% extra
#         # lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
#         #                     +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1] \
#         #                     +effluent.imol['SuccinicAcid']/self.neutralization_rxns.X[2]) \
#         #                     * 1.05
#         # effluent.mol += lime.mol
#         # self.neutralization_rxns(effluent.mol)
        
#         # !!! Temporary alteration for 2,3-HP production
#         # Rxns. can't parse chemical strings with ",", "-"
#         # effluent.imol['2,3-Butanediol'] = effluent.imol['LacticAcid']
#         # effluent.imol['LacticAcid'] = 0
        
#         # effluent.imass['LacticAcid'] = 0
#         # effluent.imass['AceticAcid'] = (8/1000) * effluent.imass['H2O']
#         # effluent.imass['2,3-Butanediol'] = (109.9/1000) * effluent.imass['H2O']
#         # effluent.imass['Glucose'] = 1.5 * effluent.F_vol
#         # effluent.imass['Xylose'] = 1.5 * effluent.F_vol
        
#     def _design(self):
#         Design = self.design_results
#         total_mass_flow = self.ins[0].F_mass + self.ins[1].F_mass
#         Design['Saccharification tank size'] = total_mass_flow * self.tau_saccharification
#         Design['Slurry flow rate'] = total_mass_flow
#         Design['Fermenter size'] = self.outs[0].F_mass * self.tau_cofermentation
#         Design['Recirculation flow rate'] = total_mass_flow



@cost(basis='Saccharification tank size', ID='Saccharification tanks', units='kg',
      cost=3840000, S=421776*24, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Slurry flow rate', ID='Saccharification transfer pumps', units='kg/hr',
      kW=74.57, cost=47200, S=421776*24, CE=CEPCI[2009], n=0.8, BM=2.3)
class Saccharification(Unit):    
    _N_ins = 1
    _N_outs = 2
    _units= {'Saccharification tank size': 'kg',
             'Slurry flow rate': 'kg/hr'}             
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_saccharification = 84 # Humbird
    
    #!!! This needs to be updated
    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.07
    
    # CSL_loading = 10 # kg/m3
    
    def __init__(self, ID='', ins=None, outs=(), T=50+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        
        # Based on Table 9 on Page 28 of Humbird et al.
        #!!! Yalin updated the glucan-to-glucose conversion to match the table
        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.9),
            Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
            ])
        
        # self.saccharified_stream = Stream(None)

    def _run(self):
        feed = self.ins[0]
        ss, sidedraw = self.outs
        # ss = self.saccharified_stream
        
        ss.copy_like(feed)
        
        ss.T = self.T
        # CSL.imass['CSL'] = feed.F_vol * self.CSL_loading 
        # ss.mol += CSL.mol
        self.saccharification_rxns(ss.mol)
        # Sidedraw to SeedTrain
        sidedraw.mol = ss.mol * self.inoculum_ratio
        ss.mol = ss.mol - sidedraw.mol


        
    def _design(self):
        Design = self.design_results
        total_mass_flow = self.ins[0].F_mass 
        Design['Saccharification tank size'] = total_mass_flow * self.tau_saccharification
        Design['Slurry flow rate'] = total_mass_flow
        

@cost(basis='Fermenter size', ID='Fermentors', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36), CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Agitators', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      # and total residence time (batch hydrolysis and fermentation)
      kW=22.371, cost=52500, S=(42607+443391+948+116)*(60+36), CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Recirculation flow rate', ID='Recirculation pumps', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)

class CoFermentation_original(Unit):    
    _N_ins = 3
    _N_outs = 2
    _units= {'Fermenter size': 'kg',
             'Recirculation flow rate': 'kg/hr'}             

    
    # Co-Fermentation time (hr)
    tau_cofermentation = 120
    

    CSL_loading = 10 # kg/m3
    
    def __init__(self, ID='', ins=None, outs=(), T=30+273.15):
        Unit.__init__(self, ID, ins, outs)
        # Same T for saccharificatoin and co-fermentation
        self.T = T
        

        # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2HP + CO2',        'Glucose',   .80),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.03),
        Rxn('Xylose -> 2HP + CO2',       'Xylose',    0.64),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.024),
        ])
        
        self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
        self.xylose_to_HP_rxn = self.cofermentation_rxns[2]
        
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[1]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[3]
        
        # self.cofermentation_rxns[1].X = \
        #     max(0, 1- (.07 + self.cofermentation_rxns[0].X + self.cofermentation_rxns[2].X))
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        # self.neutralization_rxns = ParallelRxn([
        # #   Reaction definition                                               Reactant  Conversion
        # Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'LacticAcid',   1),
        # Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1),
        # Rxn('SuccinicAcid + CalciumDihydroxide -> CalciumSuccinate + 2H2O', 'SuccinicAcid', 1)
        #     ])

    def _run(self):
        feed, sugars, CSL = self.ins
        
        effluent, vapor = self.outs
        effluent.mix_from([feed, sugars])
        # ss = Stream(None)
        # effluent.copy_like(feed)
        effluent.T = vapor.T = self.T
        CSL.imass['CSL'] = feed.F_vol * self.CSL_loading 
        
        self.cofermentation_rxns(effluent.mol)
        vapor.imol['CO2'] = effluent.imol['CO2']
        vapor.phase = 'g'
        
        effluent.imol['CO2'] = 0
        effluent.imass['CSL'] = 0
        
    def _design(self):
        Design = self.design_results
        total_mass_flow = sum([instream.F_mass for instream in self.ins])
        Design['Fermenter size'] = self.outs[0].F_mass * self.tau_cofermentation
        Design['Recirculation flow rate'] = total_mass_flow
        
        
# Seed train, 5 stages, 2 trains
@cost(basis='Seed fermenter size', ID='Stage #1 fermenters', units='kg',
      # 44339, 211, and 26 are streams 303, 309, and 310 in Humbird et al.
      cost=75400, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #2 fermenters', units='kg',
      cost=116600, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #3 fermenters', units='kg',
      cost=157600, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #4 fermenters', units='kg',
      cost=352000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #4 agitators', units='kg',
      kW=11.1855, cost=26000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Seed fermenter size', ID='Stage #5 fermenters', units='kg',
      cost=1180000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #5 agitators', units='kg',
      kW=14.914, cost=43000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pumps', units='kg/hr',
      kW=59.656, cost=24300, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedTrain(Unit):
    _N_ins = 1
    _N_outs = 1
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}
    
    # Cycle time for each batch (hr), including 12 hr turnaround time 
    tau_batch = 36
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    #!!! Should this T be changed to 30°C?
    def __init__(self, ID='', ins=None, outs=(), T=50+273.15, ferm_ratio=0.9):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.ferm_ratio = ferm_ratio

        # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2HP + CO2',        'Glucose',   .53*.8), 
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.03*.8),
        Rxn('Xylose -> 2HP + CO2',       'Xylose',    0.53*.8*.8),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.03*.8*.8),
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

# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=439000, S=40414, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=11.1855, cost=31800, S=40414, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=8200, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedHoldTank(Unit): pass


# %% 

# =============================================================================
# Separation
# =============================================================================

# Filter to separate fermentation broth into products liquid and solid
@cost(basis='Solids flow rate', ID='Feed tank', units='kg/hr',
      cost=174800, S=31815, CE=CEPCI[2010], n=0.7, BM=2.0)
@cost(basis='Solids flow rate', ID='Feed pump', units='kg/hr',
      kW=74.57, cost=18173, S=31815, CE=CEPCI[2010], n=0.8, BM=2.3)
@cost(basis='Pressing air flow rate', ID='Filter pressing compressor', units='kg/hr',
      kW=111.855, cost=75200, S=808, CE=CEPCI[2009], n=0.6, BM=1.6)
@cost(basis='Solids flow rate', ID='Pressing air compressor reciever', units='kg/hr',
      cost=8000, S=31815, CE=CEPCI[2010], n=0.7, BM=3.1)
@cost(basis='Drying air flow rate', ID='Filter drying compressor', units='kg/hr',
      kW=1043.98, cost=405000, S=12233, CE=CEPCI[2009], n=0.6, BM=1.6)
@cost(basis='Solids flow rate', ID='Dry air compressor reciever', units='kg/hr',
      cost=17000, S=31815, CE=CEPCI[2010], n=0.7, BM=3.1)
@cost(basis='Solids flow rate', ID='Pressure filter', units='kg/hr',
      cost=3294700, S=31815, CE=CEPCI[2010], n=0.8, BM=1.7)
@cost(basis='Solids flow rate', ID='Filtrate discharge pump', units='kg/hr',
      # Power not specified, based on filtrate tank discharge pump
      kW=55.9275, cost=13040, S=31815, CE=CEPCI[2010], n=0.8, BM=2.3)
@cost(basis='Solids flow rate', ID='Filtrate tank', units='kg/hr',
      cost=103000, S=31815, CE=CEPCI[2010], n=0.7, BM=2.0)
@cost(basis='Filtrate flow rate', ID='Flitrate tank agitator', units='kg/hr',
      kW=5.59275, cost=26000,  S=337439, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Solids flow rate', ID='Filtrate tank discharge pump', units='kg/hr',
      kW=55.9275, cost=13040, S=31815, CE=CEPCI[2010], n=0.8, BM=2.3)
@cost(basis='Solids flow rate', ID='Cell mass wet cake conveyor', units='kg/hr',
      kW=7.457, cost=70000, S=28630, CE=CEPCI[2009], n=0.8, BM=1.7)
@cost(basis='Solids flow rate', ID='Cell mass wet cake screw',  units='kg/hr',
      kW=11.1855, cost=20000, S=28630, CE=CEPCI[2009], n=0.8, BM=1.7)
@cost(basis='Solids flow rate', ID='Recycled water tank', units='kg/hr',
      cost=1520,  S=31815, CE=CEPCI[2010], n=0.7, BM=3.0)
@cost(basis='Solids flow rate', ID='Manifold flush pump', units='kg/hr',
      kW=74.57, cost=17057, S=31815, CE=CEPCI[2010], n=0.8, BM=2.3)
@cost(basis='Solids flow rate', ID='Cloth wash pump', units='kg/hr',
      kW=111.855,cost=29154, S=31815, CE=CEPCI[2010], n=0.8, BM=2.3)
class CellMassFilter(SolidsSeparator):
    _N_ins = 1
    _units= {'Solids flow rate': 'kg/hr',
             'Pressing air flow rate': 'kg/hr',
             'Drying air flow rate': 'kg/hr',
             'Filtrate flow rate': 'kg/hr'}
            
    def _design(self):
        Design = self.design_results
        # 809 is the scailng basis of equipment M-505,
        # 391501 from stream 508 in Humbird et al.
        Design['Pressing air flow rate'] = 809/391501 * self.ins[0].F_mass
        # 12105 and 391501 from streams 559 and 508 in Humbird et al.
        Design['Drying air flow rate'] = 12105/391501 * self.ins[0].F_mass
        Design['Solids flow rate'] = self.outs[0].F_mass
        Design['Filtrate flow rate'] = self.outs[1].F_mass

from math import pi, ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.exceptions import DesignError
class Reactor(Unit, PressureVessel, isabstract=True):
    '''    
    Create an abstract class for reactor unit, purchase cost of the reactor
    is based on volume calculated by residence time.

    Parameters
    ----------
    ins : stream
        Inlet.        
    outs : stream
        Outlet.
    tau=0.5 : float
        Residence time [hr].        
    V_wf=0.8 : float
        Fraction of working volume over total volume.        
    kW_per_m3=0.985: float
        Power usage of agitator
        (0.985 converted from 5 hp/1000 gal as in [1], for liquid–liquid reaction or extraction).
    wall_thickness_factor=1: float
        A safety factor to scale up the calculated minimum wall thickness.
    vessel_material : str, optional
        Vessel material. Default to 'Stainless steel 316'.
    vessel_type : str, optional
        Vessel type. Can only be 'Horizontal' or 'Vertical'.
        
    References
    ----------
    .. [1] Seider, W. D.; Lewin, D. R.; Seader, J. D.; Widagdo, S.; Gani, R.; 
        Ng, M. K. Cost Accounting and Capital Cost Estimation. In Product 
        and Process Design Principles; Wiley, 2017; pp 470.
    '''
    _N_ins = 2
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3'}
    
    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147 
    
    def __init__(self, ID='', ins=None, outs=(), *, 
                  P=101325, tau=0.5, V_wf=0.8,
                  length_to_diameter=2, kW_per_m3=0.985,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical'):
        
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _design(self):
        Design = self.design_results
        ins_F_vol = self.F_vol_in
        V_total = ins_F_vol * self.tau / self.V_wf
        P = self.P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        wall_thickness_factor = self.wall_thickness_factor
        
        N = ceil(V_total/self._V_max)
        if N == 0:
            V_reactor = 0
            D = 0
            L = 0
        else:
            V_reactor = V_total / N
            D = (4*V_reactor/pi/length_to_diameter)**(1/3)
            D *= 3.28084 # convert from m to ft
            L = D * length_to_diameter

        Design['Residence time'] = self.tau
        Design['Total volume'] = V_total
        Design['Single reactor volume'] = V_reactor
        Design['Number of reactors'] = N
        Design.update(self._vessel_design(P, D, L))
        if wall_thickness_factor == 1: pass
        elif wall_thickness_factor < 1:
            raise DesignError('wall_thickness_factor must be larger than 1')
        else:
              Design['Wall thickness'] *= wall_thickness_factor
              # Weight is proportional to wall thickness in PressureVessel design
              Design['Weight'] = round(Design['Weight']*wall_thickness_factor,2)
            
    def _cost(self):
        Design = self.design_results
        purchase_costs = self.purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in purchase_costs.items():
                purchase_costs[i] = 0
        
        else:
            purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in purchase_costs.items():
                purchase_costs[i] *= Design['Number of reactors']
            
            self.power_utility(self.kW_per_m3 * Design['Total volume'])
    # def _run(self):
    #     PressureVessel._run()
    @property
    def BM(self):
        vessel_type = self.vessel_type
        if not vessel_type:
            raise AttributeError('vessel_type not defined')
        elif vessel_type == 'Vertical':
            return self.BM_vertical
        elif vessel_type == 'Horizontal':
            return self.BM_horizontal 
        else:
            raise RuntimeError("invalid vessel type")

class AcidulationReactor(Reactor):
    _N_ins = 2
    _N_outs = 1
    
    acidulation_rxns = ParallelRxn([
        #   Reaction definition                                           Reactant        Conversion
        Rxn('CalciumLactate + H2SO4 -> 2 HP + CaSO4',         'CalciumLactate',       1),
        Rxn('CalciumDihydroxide + H2SO4 -> CaSO4 + 2 H2O',            'CalciumDihydroxide',   1)
            ])
            
    def _run(self):
        feed, acid = self.ins
        effluent = self.outs[0]
        rxns = self.acidulation_rxns
        chemicals = self.chemicals        
        
        acid_index = chemicals.index('H2SO4')
        reactant_indices = chemicals.indices(rxns.reactants)
        needed_acid = 0
        for i in range(len(reactant_indices)):
            index = reactant_indices[i]
            needed_acid += -(rxns.stoichiometry[i][acid_index])/rxns.X[i] * feed.mol[index]
        
        # Set feed acid mol to match acidulation needs with 5% extra
        acid.imol['H2SO4'] = needed_acid * 1.05
        acid.imass['H2O'] = acid.imass['H2SO4'] / 0.93 * 0.07 # 93% purity
        effluent.mix_from([feed, acid])
        rxns.adiabatic_reaction(effluent)
        
# Filter to separate gypsum from the acidified fermentation broth
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

# class Esterification(Reactor):
#     """
#     Create an esterification reactor that converts organic acids and ethanol
#     to corresponding ethyl esters and water. Finds the amount of catalyst 
#     'Amberlyst-15' required as well as the loss occured by physicochemical attrition.
    
#     Parameters
#     ----------
#     ins : 
#         [0] Main broth
#         [1] Recycled ethanol stream 1
#         [2] Recycled lactic acid stream
#         [3] Supplementary ethanol
#         [4] Recycled ethanol stream 2
    
#     outs : 
#         [0] Main effluent
#         [1] Wastewater stream (discarded recycles)
    
#     ethanol2acids : float
#         Ethanol feed to total acid molar ratio.
#     T : float
#         Operating temperature (K).
#     """
#     _N_ins = 5
#     _N_outs = 2
#     _N_heat_utilities = 1

#     _BM = {**Reactor._BM,
#            'Heat exchangers': 3.17,
#            'Amberlyst-15 catalyst': 1}
    
#     cat_load = 0.039 # wt% of total mass

#     ethanol2acids = 1.5
    
#     # Used in uncertainty analysis to adjust the conversion
#     X_factor = 1
    
#     reactives = ('LacticAcid', 'Ethanol', 'H2O', 'EthylLactate',
#                  'AceticAcid', 'EthylAcetate', 'SuccinicAcid', 'EthylSuccinate')
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
#                  T=351.15, P=101325, tau=None, tau_max=15, 
#                  V_wf=0.8, length_to_diameter=2, kW_per_m3=0.985,
#                  X1=None, X2=None, assumeX2equalsX1=True, allow_higher_T=False,
#                  wall_thickness_factor=1,
#                  vessel_material='Stainless steel 316',
#                  vessel_type='Vertical'):
        
#         Unit.__init__(self, ID, ins, outs)
        
#         self.T = T
#         self.P = P
#         self.tau_max = tau_max
#         self.V_wf = V_wf
#         self.length_to_diameter = length_to_diameter
#         self.kW_per_m3 = kW_per_m3
#         self.X1, self.X2, self._tau = X1, X2, tau
#         self.assumeX2equalsX1 = assumeX2equalsX1
#         self.allow_higher_T = allow_higher_T
#         self.wall_thickness_factor = wall_thickness_factor
#         self.vessel_material = vessel_material
#         self.vessel_type = vessel_type
#         self.heat_exchanger = HXutility(None, None, None, T=T)

#     def compute_coefficients(self, T):
#         K = self.K = exp(2.9625 - 515.13/T)
#         kc = self.kc = 2.70 * (1e7) * exp(-6011.55/T)
#         KW = self.KW = 15.19 * exp(12.01/T)
#         KEt = self.KEt = 1.22 * exp(359.63/T)
#         return K, kc, KW, KEt
    
#     def compute_r(self, flow, reactives, T):
#         lle_chemicals = self.chemicals.lle_chemicals
#         lle_IDs = tuple([i.ID for i in lle_chemicals])
#         f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
#         gammas = f_gamma(flow.get_normalized_mol(lle_IDs), T)
#         gammas_reactives = np.array([gammas[lle_IDs.index(ID)] 
#                                      for ID in reactives])
#         normalized_mol = flow.get_normalized_mol(lle_IDs)
#         curr_conc = np.array([normalized_mol[lle_IDs.index(ID)] 
#                               for ID in reactives])
#         activities = gammas_reactives * curr_conc
#         r_numerator = self.kc * (activities[1]*activities[0]-
#                                  (activities[3]*activities[2]/self.K))
#         r_denominator = (1+self.KEt*activities[3]+self.KW*activities[2])**2
#         r = r_numerator / r_denominator
#         return r

#     def compute_X1_and_tau(self, mixed_stream, time_step):
#         T = self.T
#         cat_load = self.cat_load
#         reactives = self.reactives[0:4]
#         compute_r = self.compute_r
#         time_max = self.tau_max * 60 # tau_max in hr
#         K, kc, KW, KEt = self.compute_coefficients(T)
        
#         temp_flow = mixed_stream.copy()
#         self.mcat = mcat = cat_load * temp_flow.F_mass
#         r = compute_r(temp_flow, reactives, T)
#         dX = r * time_step * mcat / 1000 # r is in mol g-1 min-1
        
#         curr_flow = temp_flow.get_flow('kmol/hr', reactives)
#         new_flows = [1, 1, 1, 1]
#         LA_initial = temp_flow.imol['LacticAcid']
        
#         tau_min = time_step # tau in min
#         while dX/LA_initial>1e-4:            
#             if curr_flow[0]<dX or curr_flow[1]<dX:
#                 dX = min(curr_flow[0], curr_flow[1])
                
#             new_flows = [curr_flow[0]-dX, # LA
#                          curr_flow[1]-dX, # ethanol
#                          curr_flow[2]+dX, # water
#                          curr_flow[3]+dX] # EtLA
            
#             temp_flow.set_flow(new_flows, 'kmol/hr', reactives)
            
#             # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
#             if new_flows[0]<=0 or new_flows[1]<=0 or tau_min>time_max-time_step: 
#                 break

#             r = compute_r(temp_flow, reactives, T)
#             dX = r * time_step * mcat / 1000  # r is in mol g-1 min-1
#             curr_flow = temp_flow.get_flow('kmol/hr', reactives)
#             tau_min += time_step
        
#         LA_in_feeds = mixed_stream.imol['LacticAcid']
#         X1 = (LA_in_feeds-temp_flow.imol['LacticAcid']) / LA_in_feeds
#         tau = tau_min / 60 # convert min to hr
#         return X1, tau

#     @property
#     def tau(self):
#         """Residence time (hr)."""
#         return self._tau

#     def _run(self):
#         # On weight basis, ethanol1 is ~91% ethanol with 8% water,
#         # ethanol2 is ~97% ethanol with ~3% water, so should first recycle ethanol2
#         feed, ethanol1, recycled_LA, supplement_ethanol, ethanol2 = self.ins 
#         effluent, wastewater = self.outs
        
#         acids = ('LacticAcid', 'AceticAcid', 'SuccinicAcid')
#         # Succnic acid is a dicarboxylic acid, needs twice as much ethanol
#         ratios = self.ethanol2acids * np.array([1, 1, 2])
        
#         feeds = feed.copy()
#         feeds.mix_from([feed, recycled_LA])
        
#         # Have enough ethanol in feed and recycle2, discharge some recycle2
#         # and all of recycle1
#         if compute_extra_chemical(feeds, ethanol2, acids, 'Ethanol', ratios) > 0:
#             effluent, ethanol2_discarded = \
#                 adjust_recycle(feeds, ethanol2, acids, 'Ethanol', ratios)
#             wastewater.mix_from([ethanol1, ethanol2_discarded])
#             supplement_ethanol.empty()
        
#         else:
#             # Recycle all of ethanol2 and combine feed and recycle2 as feed2
#             feeds2 = feeds.copy()
#             feeds2.mix_from([feeds, ethanol2])
#             # Have enough ethanol in feed2 and ethanol1
#             if compute_extra_chemical(feeds2, ethanol1, acids, 'Ethanol', ratios) > 0:
#                 effluent, ethanol1_discarded = \
#                     adjust_recycle(feeds2, ethanol1, acids, 'Ethanol', ratios)
#                 wastewater = ethanol1_discarded
#                 supplement_ethanol.empty()
#             # Not have enough ethanol in both recycles, need supplementary ethanol
#             else:
#                 supplement_ethanol.imol['Ethanol'] = \
#                     - compute_extra_chemical(feeds2, ethanol1, acids, 'Ethanol', ratios)
#                 effluent.mix_from(self.ins)
#                 wastewater.empty()

#         if self.allow_higher_T and self.T<effluent.T:
#             self.T = effluent.T
        
#         if self.X1 == None and self.tau == None:
#             X1, tau = self.compute_X1_and_tau(effluent, time_step=1)
#             self.X1 = X1
#             self._tau = tau
#         elif not self.X1 and not self.tau:
#             raise AttributeError('X1 and tau must be both defined, or both as None')
#         else:
#             X1 = self.X1
            
#         if self.assumeX2equalsX1:
#             X2 = self.X2 = X1
#         else:
#             if not self.X2:
#                 raise AttributeError('X2 must be defined if assumeX2equalsX1 is False')
        
#         X1 *= self.X_factor
#         X2 *= self.X_factor

#         self.esterification_rxns = ParallelRxn([
#             #   Reaction definition                                     Reactant  Conversion
#             Rxn('LacticAcid + Ethanol -> EthylLactate + H2O',         'LacticAcid',   X1),
#             Rxn('AceticAcid + Ethanol -> EthylAcetate + H2O',         'AceticAcid',   X2),
#             # Assume succinic acid has the same conversion as acetic acid
#             Rxn('SuccinicAcid + 2 Ethanol -> EthylSuccinate + 2 H2O', 'SuccinicAcid', X2)
#             ])

#         self.esterification_rxns(effluent.mol)
#         effluent.T = self.T
#         self.outs[0].copy_like(effluent)
#         self.outs[1].copy_like(wastewater)
        
#     def _cost(self):
#         super()._cost()
#         hx = self.heat_exchanger
#         N = self.design_results['Number of reactors']
#         single_rx_effluent = self.outs[0].copy()
#         single_rx_effluent.mol[:] /= N
#         hx.simulate_as_auxiliary_exchanger(duty=self.Hnet/N, 
#                                            stream=single_rx_effluent)
#         hu_total = self.heat_utilities[0]
#         hu_single_rx = hx.heat_utilities[0]
#         hu_total.copy_like(hu_single_rx)
#         hu_total.scale(N)
#         self.purchase_costs['Heat exchangers'] = hx.purchase_cost * N
#         self.purchase_costs['Amberlyst-15 catalyst'] = self.mcat * price['Amberlyst15']
        
# class HydrolysisReactor(Reactor):
#     """
#     Create a hydrolysis reactor that hydrolyze organic acid esters into 
#     corresponding acids and ethanol. 
    
#     Parameters
#     ----------
#     ins : 
#         [0] Main broth
#         [1] Supplementary water
#         [2] Recycled water stream 1
#         [3] Recycled water stream 2
    
#     outs : 
#         [0] Main effluent
#         [1] Wastewater stream (discarded recycles)
    
#     water2esters : float
#         Water feed to total ester molar ratio.
#     """
#     _N_ins = 4
#     _N_outs = 2
#     water2esters = 12
    
#     hydrolysis_rxns = ParallelRxn([
#             #   Reaction definition                                       Reactant   Conversion
#             Rxn('EthylLactate + H2O -> LacticAcid + Ethanol',         'EthylLactate',   0.8),
#             Rxn('EthylAcetate + H2O -> AceticAcid + Ethanol',         'EthylAcetate',   0.8),
#             Rxn('EthylSuccinate + 2 H2O -> SuccinicAcid + 2 Ethanol', 'EthylSuccinate', 0.8),
#                 ])
    
#     def _run(self):
#         # On weight basis, recycle2 is near 10% EtLA so will always be recycled,
#         # but recycle1 is >97% water with <1% LA, so will only be used to supply
#         # water needed for the hydrolysis reaction
#         feed, water, recycle1, recycle2 = self.ins
#         effluent, wastewater = self.outs
        
#         esters = ('EthylLactate', 'EthylAcetate', 'EthylSuccinate')
#         # Succnic acid is a dicarboxylic acid, needs twice as much water
#         ratios = self.water2esters * np.array([1, 1, 2])
#         # Have enough water in feed and recycle2, discharge some recycle2
#         # and all of recycle1
#         if compute_extra_chemical(feed, recycle2, esters, 'H2O', ratios) > 0:
#             effluent, recycle2_discarded = \
#                 adjust_recycle(feed, recycle2, esters, 'H2O', ratios)
#             wastewater.mix_from([recycle1, recycle2_discarded])
#             water.empty()        
#         else:
#             # Recycle all of recycle2 and combine feed and recycle2 as feed2
#             feed2 = feed.copy()
#             feed2.mix_from([feed, recycle2])
#             # Have enough water in feed2 and recycle1
#             if compute_extra_chemical(feed2, recycle1, esters, 'H2O', ratios) > 0:
#                 effluent, recycle1_discarded = \
#                     adjust_recycle(feed2, recycle1, esters, 'H2O', ratios)
#                 wastewater = recycle1_discarded
#                 water.empty()
#             # Not have enough water in both recycles, need supplementary water
#             else:
#                 water.imol['H2O'] = \
#                     - compute_extra_chemical(feed2, recycle1, esters, 'H2O', ratios)
#                 effluent.mix_from(self.ins)
#                 wastewater.empty()
        
#         rxns = self.hydrolysis_rxns
#         rxns(effluent.mol)
#         self.outs[0].copy_like(effluent)
#         self.outs[1].copy_like(wastewater)


# %% 

# =============================================================================
# Wastewater treatment
# =============================================================================

# Total cost of wastewater treatment is combined into this placeholder
@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=CEPCI[2010], n=0.6, BM=1)
class WastewaterSystemCost(Unit): pass

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
        	
    digestion_rxns: 
        [ReactionSet] Anaerobic digestion reactions.  	
    sludge_split: 
        [Array] Split between wastewater and sludge	
    	
    """
    auxiliary_unit_names = ('heat_exchanger',)
    _N_ins = 1	
    _N_outs = 3
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants, split=(), T=35+273.15):	
        Unit.__init__(self, ID, ins, outs)	
        self.reactants = reactants	
        self.split = split	
        self.multi_stream = MultiStream(None)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.heat_utilities = hx.heat_utilities
        chems = self.chemicals	
        	
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
        T = self.T	

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
        
    def _design(self):
        wastewater = self.ins[0]
        # Calculate utility needs to keep digester temperature at 35°C,	
        # heat change during reaction is not tracked	
        H_at_35C = wastewater.thermo.mixture.H(mol=wastewater.mol, 	
                                               phase='l', T=self.T, P=101325)	
        duty = -(wastewater.H - H_at_35C)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, wastewater)
        
# class HPLiquidsSplitSettler():
#     def __init__(self, ID='', ins=None, outs=(), *, part_coeffs = [], solvents = (), split = []):
#         self.part_coeffs = part_coeffs
        
#     def run(self):
#         split = self.split
#         solvents = self.solvents
#         for i in len(part_coeffs):
#             k = part_coeffs[i]
            
#             if i not in solvents:
#                 split[i] = part_coeffs[i] * self.ins[0].imass[self.solvents[0]]/self.solvents[1]
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
        
    digestion_rxns : 
        [ReactionSet] Anaerobic digestion reactions
    
    """
    
    _N_ins = 3
    _N_outs = 2
    # 4350, 4379, 356069, 2252, 2151522, and 109089 are water flows from 
    # streams 622, 630, 611, 632, 621, and 616  in Humbird et al.
    evaporation = 4350/(4379+356069+2252+2151522+109089)
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants, ratio=0):
        Unit.__init__(self, ID, ins, outs)
        self.reactants = reactants
        self.ratio = ratio
        chems = self.chemicals
        
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
        
        # Assume NaOH is completely consumed by H2SO4 and digestion products
        effluent.imol['NaOH'] = 0


# %% 

# =============================================================================
# Storage
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2010], n=0.7, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3)
class SulfuricAcidStorageTank(Unit): pass

#!!! Instead of copying the H2SO4, I suggest we use BioSTEAM's storage tank
# Dipotassium hydrogen phosphate storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2010], n=0.7, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3)
class DPHPStorageTank(StorageTank): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      # Original size basis for NH3 instead of NH4OH
      cost=196000, S=1171/17.031*35.046, CE=CEPCI[2010], n=0.7, BM=2)
class AmmoniaStorageTank(Unit): pass

#!!! Humbird's ethanol storage tank is too big
# @cost(basis='Flow rate', ID='Tank', units='kg/hr',
#       cost=1340000, S=22681, CE=CEPCI[2009], n=0.7, BM=1.7)
# @cost(basis='Flow rate', ID='Pump', units='kg/hr',
#       kW=7.457, cost=9200, S=22681, CE=CEPCI[2009], n=0.8, BM=3.1)
# class EthanolStorageTank(StorageTank): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=70000, S=1393, CE=CEPCI[2009], n=0.7, BM=2.6)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21200, S=1393, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=1393, CE=CEPCI[2009], n=0.8, BM=3.1)
class CSLstorageTank(Unit): pass

# For storage of lime used in separation and waste treatment,
# design copied from cornstover biorefinery in Aden et al.
# Base flow from stream 27 of Aden et al.
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

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=803000, S=8343, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=15000, S=8343, CE=CEPCI[2009], n=0.8, BM=3.1)
class FireWaterTank(Unit): pass

class SugarsMixer(Mixer):
    _N_ins = 3
    
# Modified from bst.units.StorageTank, which won't simulate for 0 flow 	
class HPStorageTank(StorageTank):
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.purchase_costs['Tanks'] = 0
        else: StorageTank._cost(self)

# Modified from bst.units.Pump, which won't simulate for 0 flow 
class HPPump(Pump):
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
        

class DPHP_Separation(Unit):
    _N_outs = 2
    def _run(self): pass
    def _cost(self): pass

class DehydrationReactor(Reactor):
    """
    A dehydration reactor.
    """
    _N_ins = 1
    _N_outs = 1
    
    _N_heat_utilities = 1
    _BM = {**Reactor._BM,
            'TiO2 catalyst': 1,
            'Heat exchangers': 3.17}
    mcat_frac = (12/1.5) * (1e3)# kg per m3/h
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('HP -> AA + H2O',         'HP',   0.999) # 99.9% conversion for pure HP stream, assumed theoretical maximum for dilute stream
                ])                                      # Dishisha et al. 2015
    HP_to_AA_rxn = dehydration_rxns[0]
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, T=230+273.15,
                  P=101325, V_wf=0.8, length_to_diameter=2, tau = 1,
                  kW_per_m3=0.0985, # Perry's handbook
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 304',
                  vessel_type='Vertical', X = 0.999):
        Unit.__init__(self, ID, ins, outs)
        
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.heat_exchanger = HXutility(None, None, None, T=T)
        self.tau = tau
        self.X = self.dehydration_rxns[0].X = X
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        effluent.mix_from([feed])
        effluent.T = self.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(effluent.mol)
        
   
    def _cost(self):
        super()._cost()
        hx = self.heat_exchanger
        N = self.design_results['Number of reactors']
        single_rx_effluent = self.ins[0].copy()
        single_rx_effluent.mol[:] /= N
        hx.simulate_as_auxiliary_exchanger(duty=self.Hnet/N, 
                                            stream=single_rx_effluent)
        hu_total = self.heat_utilities[0]
        hu_single_rx = hx.heat_utilities[0]
        hu_total.copy_like(hu_single_rx)
        hu_total.scale(N)
        self.purchase_costs['Heat exchangers'] = hx.purchase_cost * N
        self.purchase_costs['TiO2 catalyst'] =\
            self.mcat_frac * (sum([stream.F_vol for stream in self.outs])) * price['TiO2']
        
        
        
compute_HP_titer = lambda effluent: (effluent.imass['HP'] +
            effluent.imol['CalciumLactate']*2*HP_chemicals.HP.MW)/effluent.F_vol

_316_over_304 = 1.2



@cost(basis='Fermenter size', ID='Fermenter', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Agitator', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [3])
      # and total residence time (batch hydrolysis and fermentation)
      kW=268.452, cost=630000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Recirculation flow rate', ID='Recirculation pump', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [3])
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)  
class CoFermentation(Reactor):
    _N_ins = 3
    _N_outs = 2
    _units= {**Reactor._units,
            'Fermenter size': 'kg',
            'Recirculation flow rate': 'kg/hr',
            'Duty': 'kJ/hr'}
    _BM = {**Reactor._BM,
            'Heat exchangers': 3.17}

    auxiliary_unit_names = ('heat_exchanger',)
    
    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.07
    
    CSL_loading = 10 # g/L (kg/m3)
    
    titer_limit = 50 # in g/L (kg/m3), the maximum titer in collected data
    
    effluent_titer = 0
    
    productivity = 0.89 # in g/L/hr
    
    yield_limit = 0.76 # in g/g-sugar
    
    tau_batch_turnaround = 12 # in hr, the same as the seed train in ref [3]

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, T=50+273.15,
                  P=101325, V_wf=0.8, length_to_diameter=2,
                  kW_per_m3=0.0985, # Perry's handbook
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical',
                  neutralization=False,
                  mode='Continuous', # Batch or Continuous
                  allow_dilution=False,
                  allow_concentration=False):
        
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.neutralization = neutralization
        self.mode = mode        
        self.allow_dilution = allow_dilution
        self.allow_concentration = allow_concentration
        self.mixed_feed = tmo.Stream('mixed_feed')
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.heat_utilities = hx.heat_utilities
        
        # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2HP + CO2 + H2O',        'Glucose',   .53),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.03),
        Rxn('Xylose -> 2HP + CO2 + H2O',       'Xylose',    0.53*0.8),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.03*0.8),
        ])
        
        self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
        self.xylose_to_HP_rxn = self.cofermentation_rxns[2]
        
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[1]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[3]
        
        self._X = self.cofermentation_rxns.X.copy()
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('2 HP + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'HP',   1.),
        # Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1),
        # Rxn('SuccinicAcid + CalciumDihydroxide -> CalciumSuccinate + 2H2O', 'SuccinicAcid', 1)
            ])

    def _run(self):
        
        sugars, feed, CSL, lime = self.ins
        
        effluent, vapor = self.outs
        effluent.mix_from([feed, sugars])
        
        # ss = Stream(None)
        # effluent.copy_like(feed)
        effluent.T = vapor.T = self.T
        CSL.imass['CSL'] = (sugars.F_vol + feed.F_vol) * self.CSL_loading 
        
        self.cofermentation_rxns(effluent.mol)
        vapor.imol['CO2'] = effluent.imol['CO2']
        vapor.phase = 'g'
        
        effluent.imol['CO2'] = 0
        effluent.imass['CSL'] = 0
        
        mixed_feed = self.mixed_feed
        mixed_feed.mix_from([feed, sugars, CSL])
        
        # Need lime to neutralize produced acid
        if self.neutralization:
            self.vessel_material= 'Stainless steel 316'
            # Set feed lime mol to match rate of acids production, add 10% extra
            lime.imol['Lime'] = (effluent.imol['HP']/2/self.neutralization_rxns.X[0]) * 1.1
                                # +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1] \
                                # +effluent.imol['SuccinicAcid']/self.neutralization_rxns.X[2]) \
                                # * 1.1
            effluent.imol['Lime'] = lime.imol['Lime']
            # effluent.mix_from([effluent, lime])
            self.neutralization_rxns.adiabatic_reaction(effluent)
        else:
            self.vessel_material= 'Stainless steel 316'
            lime.empty()
        self.effluent_titer = compute_HP_titer(effluent)
        
    def _design(self):
        mode = self.mode
        Design = self.design_results
        Design.clear()
        self.tau = self.effluent_titer / self.productivity
        _mixture = self._mixture = tmo.Stream(None)
        _mixture.mix_from(self.outs[0:2])
        duty = Design['Duty'] = _mixture.H - self.mixed_feed.H

        if mode == 'Batch':
            tau_cofermentation = self.tau_batch_turnaround + self.tau
            Design['Fermenter size'] = self.outs[0].F_mass * tau_cofermentation
            Design['Recirculation flow rate'] = self.F_mass_in
            self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, _mixture)
        
        elif mode == 'Continuous':
            Reactor._design(self)

        else:
            raise DesignError(f'Fermentation mode must be either Batch or Continuous, not {mode}')

    def _cost(self):
        Design = self.design_results
        purchase_costs = self.purchase_costs
        purchase_costs.clear()
        hx = self.heat_exchanger

        if self.mode == 'Batch':
            Unit._cost()
            self._decorated_cost()
            purchase_costs['Heat exchangers'] = hx.purchase_cost
            # Adjust fermenter cost for acid-resistant scenario
            if not self.neutralization:
                purchase_costs['Fermenter'] *= _316_over_304
                purchase_costs['Agitator'] *= _316_over_304

        elif self.mode == 'Continuous':
            # if not self.neutralization:
            #     self.vessel_material= 'Stainless steel 316'
                
            Reactor._cost(self)
            
            N = Design['Number of reactors']
            single_rx_effluent = self._mixture.copy()
            single_rx_effluent.mol[:] /= N
            hx.simulate_as_auxiliary_exchanger(duty=Design['Duty']/N, 
                                                stream=single_rx_effluent)
            hu_total = self.heat_utilities[0]
            hu_single_rx = hx.heat_utilities[0]
            hu_total.copy_like(hu_single_rx)
            hu_total.scale(N)
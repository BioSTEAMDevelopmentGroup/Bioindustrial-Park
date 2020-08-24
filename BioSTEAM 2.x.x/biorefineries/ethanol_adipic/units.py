#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Fri Jun 26 07:33:04 2020

References:
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040

[4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbons: Dilute-Acid and Enzymatic Deconstruction of Biomass
    to Sugars and Catalytic Conversion of Sugars to Hydrocarbons;
    NREL/TP-5100-62498; National Renewable Energy Lab (NREL), 2015.
    http://www.nrel.gov/docs/fy15osti/62498.pdf
    
[5] Vardon et al., Cis,Cis-Muconic Acid: Separation and Catalysis to Bio-Adipic
    Acid for Nylon-6,6 Polymerization. Green Chem. 2016, 18 (11), 3397–3413.
    https://doi.org/10.1039/C5GC02844B.

@author: yalinli_cabbi
"""


# %% Setup

import thermosteam as tmo
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units import Mixer, Flash, MixTank, HXutility, Pump, SolidsSeparator
from biosteam.units.decorators import cost
from thermosteam import separations
from ethanol_adipic.process_settings import price
from ethanol_adipic.chemicals import total_solids, solubles, insolubles, COD_chemicals
from ethanol_adipic.utils import CEPCI, baseline_feedflow, compute_muconic_titer, \
    compute_COD

_kg_per_ton = 907.18474
_MGD_2_m3hr = (3.78541*1e6/24) / 1e3
_GPM_2_m3hr = (3.78541*60) / 1e3
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
    # 2205 U.S. ton/day (2000 metric tonne/day) as in ref [1]
    _baseline_flow_rate = baseline_feedflow.sum()
    _cached_flow_rate = 2205


# %% 

# =============================================================================
# Pretreatment (acid/acid-base shared ones)
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=6210, S=1981, CE=CEPCI[2010],  n=0.7, BM=3)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      cost=8000, S=3720, CE=CEPCI[2009],  n=0.8, BM=2.3)
class SulfuricAcidAdditionTank(Unit):
    _N_ins = 1
    _N_outs = 1
    
    # Bseline is (18+4.1) mg/g dry biomass as in ref [1], 93% purity
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
        # 0.05 is from 1842/36629 from streams 710 and 516 in ref [1]
        water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
        mixture.mix_from([water, acid])


# Adjust pretreatment water loading, 30% from Table 5 on Page 21 of ref [1]
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

# Mix steams for pretreatment reactor heating
class SteamMixer(Unit):
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
        steam_mol = aitken_secant(f=self.P_at_flow,
                                  x0=steam_mol, x1=steam_mol+0.1, 
                                  xtol=1e-4, ytol=1e-4,
                                  args=(self.P, steam, mixed, feed))
        mixed.P = self.P
    
@cost(basis='Dry flow rate', ID='Pretreatment reactor', units='kg/hr',
      kW=5120, cost=19812400, S=83333, CE=CEPCI[2009], n=0.6, BM=1.5)
class AcidPretreatment(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = Flash._graphics
    
    def __init__(self, ID='', ins=None, outs=(), T=130+273.15):
        Unit.__init__(self, ID, ins, outs)
        self._multistream = tmo.MultiStream(None)
        self.T = T
        vapor, liquid = self.outs
        vapor.phase = 'g'
        
        self.pretreatment_rxns = ParallelRxn([
            #            Reaction definition                 Reactant   Conversion
            # Below from Table 6 on Page 22 in ref [1]
            Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.099),
            Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.003),
            Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.003),
            Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1),
            Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9),
            Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.024),
            Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.05),
            Rxn('Acetate -> AceticAcid',                     'Acetate',  1),
            Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.05),
            # Below from Page 106 of ref [1]
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

# Costs of Tank and Agitator included in the Pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=0, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=0, S=252891, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=25365, S=292407, CE=CEPCI[2010], n=0.8, BM=2.3)
class BlowdownTank(Unit): pass

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
        	
        # Water loading based on stream 275 in ref [1]
        water.imass['Water'] = 150310/1051/(35.046/17.031) * ammonia.imass['NH4OH']
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

@cost(basis='Duty', ID='Condenser', units='kJ/hr',
      # 2 is the duty in MMkca/hr
      cost=34000, S=-2*_Gcal_2_kJ, CE=CEPCI[2009], n=0.7, BM=2.2)
class WasteVaporCondenser(HXutility): pass

# Transport hydrolysate to enzymatic hydrolysis and fermentation
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=74.57, cost=22500, S=402194, CE=CEPCI[2009], n=0.8, BM=2.3)
class HydrolysatePump(Unit):
    _graphics = Pump._graphics
    
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass


# =============================================================================
# Pretreatment (base)
# =============================================================================

# Alkaline extraction
@cost(basis='Duty', ID='Water heater', units='kJ/hr',
      # 8 is duty in Gcal/hr (MMkcal/hr)
      cost=92000, S=8*_Gcal_2_kJ, CE=CEPCI[2010],  n=0.7, BM=2.2)
@cost(basis='Flow rate', ID='Conveyor', units='kg/hr',
      kW=89.484, cost=110000*3, S=277167, CE=CEPCI[2013],  n=0.8, BM=1.7)
@cost(basis='Dry flow', ID='Reactor', units='kg/hr',
      cost=5424000, S=38600, CE=CEPCI[2010],  n=1, BM=1)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=74.57, cost=22500, S=402194, CE=CEPCI[2009],  n=0.8, BM=2.3)
class DeacetylationReactor(Unit):
    _N_ins = 3
    _N_outs = 2
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr',
             'Flow rate': 'kg/hr',
             'Dry flow': 'kg/hr'}
    
    solublization_dict = {'Sucrose': 0.5,
                          'Ash': 0.66,
                          'Lignin': 0.47,
                          'Glucan': 0.02,
                          'Xylan': 0.1,
                          'Arabinan': 0.3}
    
    def __init__(self, ID='', ins=None, outs=(), T=92+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
    
    def _run(self):
        feedstock, caustic, water = self.ins
        liquor, solids = self.outs
        # 70 mg/g dry feedstock
        feedstock_dry_mass = feedstock.F_mass - feedstock.imass['Water']
        caustic.imass['NaOH'] = 70/1000 * feedstock_dry_mass

        liquor.empty()
        solids.mix_from(self.ins)
        liquor.T = solids.T = self.T
        # Inconsistency in black liquor composition in ref [2] between text
        # and streatm tables, assume Protein to be in solids
        liquor.copy_flow(solids, ('Acetate', 'Extractives', 'NaOH'), remove=True)

        for i, j in self.solublization_dict.items():
            liquor.imass[i] = solids.imass[i] * j
            solids.imass[i] -= liquor.imass[i]
        
        # 66% water content
        liquor.imass['Water'] = 0        
        liquor.imass['Water'] = 0.66/(1-0.66) * liquor.F_mass
        
        # 30% total solids content stated in text, but 25% according to stream 301
        # in ref [2], used 25%
        solids.imass['Water'] = 0
        total_solids_mass = solids.imass[total_solids].sum()
        solids.imass['Water'] = total_solids_mass/0.25 - solids.F_mass
        water.imass['Water'] = (liquor.imass['Water']+solids.imass['Water']) \
            - feedstock.imass['Water']        
        mixture = feedstock.copy()
        mixture.mix_from(self.ins)
        self.design_results['Flow rate'] = mixture.F_mass
        
    def _design(self):
        # Use Hnet (as opposed to H_out-H_in when there are reactions/change of materials)
        duty = self.Hnet
        self.heat_utilities[0](duty=duty, T_in=self.ins[0].T)
        self.design_results['Duty'] = duty
        self.design_results['Dry flow'] = \
            self.ins[0].imass[insolubles].sum()
        

# Mill solids to for increased accessible to solids
@cost(basis='Dry flow', ID='Primary disc refiner', units='kg/hr',
      cost=2466700*8, S=62942, CE=CEPCI[2013],  n=0.6, BM=1.5)
@cost(basis='Dry flow', ID='Secondary mill', units='kg/hr',
      cost=578000*11, S=62942, CE=CEPCI[2013],  n=0.6, BM=1.4)
class DiscMill(Unit):
    _units= {'Dry flow': 'kg/hr'}
    
    def _design(self):
        dry_flow_rate = self.outs[0].imass[insolubles].sum()
        self.design_results['Dry flow'] = dry_flow_rate
        # 200 kW per dry tonne
        self.power_utility(200*(dry_flow_rate/1000))

# Transport pretreated liquid for lignin utilization
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=22500, S=204390, CE=CEPCI[2009], n=0.8, BM=2.3)
class BlackLiquorPump(Unit):
    _graphics = Pump._graphics
    
    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass


# %%

# =============================================================================
# Fermentation
# =============================================================================

@cost(basis='Duty', ID='Hydrolysate cooler', units='kJ/hr',
      # 13 is the duty in MMkca/hr
      cost=85000, S=-8*_Gcal_2_kJ, CE=CEPCI[2010], n=0.7, BM=2.2)
class HydrolysateCooler(HXutility): pass

@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      kW=74.57, cost=109000, S=379938, CE=CEPCI[2009], n=0.5, BM=1.7)
class EnzymeHydrolysateMixer(Mixer):
    _N_ins = 3	
    _N_outs = 1	
    
    # Updated from ref [1] values (20 mg/g enzyme loading and 0.2 solid loading)
    # to ref [2]
    enzyme_loading = 10
    solid_loading = 0.25
    
    def _run(self):	
        hydrolysate, enzyme, water = self.ins	
        effluent = self.outs[0]	

        # 10% extra based on Page 23 of ref [2]        
        enzyme.imass['Enzyme'] = (self.enzyme_loading/1000*1.1) * hydrolysate.imass['Glucan']
        mixture = hydrolysate.copy()	
        mixture.mix_from([hydrolysate, enzyme])	
        	
        total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solid_loading	
        water.imass['Water'] = max(0, total_mass-mixture.F_mass)
        	
        effluent.mix_from([hydrolysate, enzyme, water])


@cost(basis='Flow rate', ID='Saccharification tank', units='kg/hr',
      cost=3840000, S=421776, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Saccharification transfer pump', units='kg/hr',
      kW=74.57, cost=47200,  S=421776, CE=CEPCI[2009], n=0.8, BM=2.3)
@cost(basis='Flow rate', ID='Fermentation cooler', units='kg/hr',
      cost=86928, S=421776, CE=CEPCI[2009],  n=1, BM=2.2)
@cost(basis='Flow rate', ID='Fermenter', units='kg/hr',
      cost=10128000, S=421776, CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=268.452, cost=630000, S=421776, CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Flow rate', ID='Recirculation pump', units='kg/hr',
      kW=74.57, cost=47200, S=421776, CE=CEPCI[2009], n=0.8, BM=2.3)
@cost(basis='Duty', ID='Hydrolysate cooler', units='kJ/hr',
      # 13 is the duty in MMkca/hr
      cost=85000, S=-8*_Gcal_2_kJ, CE=CEPCI[2010], n=0.7, BM=2.2)
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 4
    _N_outs = 3
    _N_heat_utilities = 3
    
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
    _units = {'Flow rate': 'kg/hr',
              'Duty': 'kJ/hr'}
    
    # Split to outs[2]
    inoculum_ratio = 0.1
    
    def __init__(self, ID='', ins=None, outs=(), P=101325, C5_saccharification=False):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.C5_saccharification = C5_saccharification
        self.saccharified_stream = tmo.Stream(None)
 
        self.saccharification_rxns_C6 = ParallelRxn([
    #   Reaction definition                   Reactant     Conversion
    Rxn('Glucan -> GlucoseOligomer',          'Glucan',      0.04),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',      0.012),
    Rxn('Glucan + H2O -> Glucose',            'Glucan',      0.9),
    Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1)])
    
        self.saccharification_rxns_C5 = ParallelRxn([
    #   Reaction definition                   Reactant     Conversion
    Rxn('Xylan + H2O -> Xylose',              'Xylan',       0.9),
    Rxn('Arabinan + H2O -> Arabinose',        'Arabinan',    0.85)])
    
        self.loss_rxns = ParallelRxn([
    #   Reaction definition               Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.03),
    Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.03),
    Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.03),
    Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.03),
    Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.03),])
    
        self.fermentation_rxns = ParallelRxn([
    #   Reaction definition                                          Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.95),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.02),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.004),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.006),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.85),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                                    'Xylose',    0.019),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.003),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.046),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.009),
    ])
    
    def _run(self):
        feed, inoculum, CSL, DAP = self.ins
        vent, effluent, sidedraw = self.outs
        vent.P = effluent.P = sidedraw.P = self.P
        ss = self.saccharified_stream
        vent.phase = 'g'
        
        # 0.25 wt% and 0.33 g/L (kg/m3) based on ref [1]
        CSL.imass['CSL'] = feed.imass['CSL'] = 0.0025 * feed.F_mass
        DAP.imass['DAP'] = feed.imass['DAP'] = 0.33 * feed.F_vol
        ss.mix_from(self.ins)
        self.saccharification_rxns_C6(ss.mol)
        if self.C5_saccharification:
            self.saccharification_rxns_C5(ss.mol)
        sidedraw.mol = ss.mol * self.inoculum_ratio
        effluent.mol = ss.mol - sidedraw.mol
        self.loss_rxns(effluent.mol)
        self.fermentation_rxns(effluent.mol)
        vent.receive_vent(effluent)
        
        ss.T = sidedraw.T = self.T_saccharification
        vent.T = effluent.T = self.T_fermentation
    
    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
        effluent = self.outs[1]
        hu_hydrolysate, hu_saccharification, hu_fermentation = self.heat_utilities
        mixture = self.thermo.mixture
        ss = self.saccharified_stream
        ss_in = ss.copy()
        ss_in.mix_from(self.ins)
        self.design_results['Duty'] = ss.H-ss_in.H
        hu_hydrolysate(duty=ss.H-ss_in.H, T_in=ss_in.T)
        
        mol = ss.mol
        duty = (mixture.H('l', mol, self.T_fermentation, 101325.)
                - mixture.H('l', mol, self.T_saccharification, 101325.))
        hu_saccharification(duty, self.T_fermentation)
        ei = effluent.chemicals.index('Ethanol')
        ethanol = (sum([i.mol[ei] for i in self.outs])
                   - sum([i.mol[ei] for i in self.ins]))
        reactor_duty = -5568 * ethanol
        hu_fermentation(reactor_duty, effluent.T)


@cost(basis='Flow rate', ID='Stage #1 reactor', units='kg/hr',
      cost=75400, S=43149, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Stage #2 reactor', units='kg/hr',
      cost=116600, S=43149, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Stage #3 reactor', units='kg/hr',
      cost=157600, S=43149, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Stage #4 reactor', units='kg/hr',
      cost=352000, S=43149, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Stage #4 agitator', units='kg/hr',
      kW=11.1855, cost=26000, S=43149, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Stage #5 reactor', units='kg/hr',
      cost=1180000, S=43149, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Stage #5 agitator', units='kg/hr',
      kW=14.914, cost=43000, S=43149, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Seed transfer pump', units='kg/hr',
      kW=59.656, cost=24800, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedTrain(Unit):
    _N_ins = 3
    _N_outs = 2
    _N_heat_utilities = 1
    
    _units= {'Flow rate': 'kg/hr'}

    #: Operating temperature (K)
    T = 32+273.15
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self.fermentation_rxns = ParallelRxn([
    #   Reaction definition                             Reactant   Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                        'Glucose',   0.04),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.004),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.006),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                        'Xylose',    0.04),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.003),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.046),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.009)])
    
    def _run(self):
        feed, CSL, DAP = self.ins
        vent, effluent = self.outs
        
        # 0.50 wt% and 0.66 g/L (kg/m3) based on ref [1]
        CSL.imass['CSL'] = 0.005 * feed.F_mass
        feed.imass['CSL'] += CSL.imass['CSL']
        DAP.imass['DAP'] = 0.67 * feed.F_vol
        feed.imass['DAP'] += DAP.imass['DAP']
        effluent.copy_flow(feed)

        self.fermentation_rxns(effluent.mol)
        effluent.T = self.T
        vent.phase = 'g'
        vent.copy_flow(effluent, ('CO2', 'NH3', 'O2'), remove=True)

    def _design(self): 
        self.design_results['Flow rate'] = self.ins[0].F_mass
        self.heat_utilities[0](self.Hnet, self.T)

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=439000, S=40414, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=11.1855, cost=31800, S=40414, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=8200, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedHoldTank(Unit): pass



# %%

# =============================================================================
# Ethanol purification
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=636000, S=425878, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=29.828, cost=68300, S=425878, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=26800, S=488719, CE=CEPCI[2009], n=0.8, BM=2.3)
class BeerTank(Unit): pass

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
        # 391501 from stream 508 in ref [1]
        Design['Pressing air flow rate'] = 809/391501 * self.ins[0].F_mass
        # 12105 and 391501 from streams 559 and 508 in ref [1]
        Design['Drying air flow rate'] = 12105/391501 * self.ins[0].F_mass
        Design['Solids flow rate'] = self.outs[0].F_mass
        Design['Filtrate flow rate'] = self.outs[1].F_mass


# %% 

# =============================================================================
# Lignin utilization
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=1317325, S=328984, CE=CEPCI[2011], n=0.7, BM=1.8)
class BlackLiquorStorage(Unit): pass


@cost(basis='Flow rate', ID='Reactor', units='kg/hr',
      cost=16300000, S=323295, CE=CEPCI[2013], n=0.6, BM=1.7)
@cost(basis='Flow rate', ID='Flash/drain tank', units='kg/hr',
      cost=262000, S=323295, CE=CEPCI[2013], n=0.7, BM=2)
class PulpingReactor(Unit):
    _N_ins = 3	
    _N_outs = 2
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr',
             'Flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), T=120+273.15, P=6.32*101325):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        
        self.deconstruction_rxns = ParallelRxn([
            #         Reaction definition        Reactant   Conversion
            Rxn('Glucan + H2O -> Glucose',       'Glucan',      0.48),
            Rxn('Xylan + H2O -> Xylose',         'Xylan',       0.48),
            Rxn('Arabinan + H2O -> Arabinose',   'Arabinan',    0.48),
            Rxn('Lignin -> SolubleLignin',       'Lignin',      0.53)
            ])
        
    def _run(self):
        liquor, residuals, caustic = self.ins
        pulp, solids = self.outs
        
        # Minimum of 2 wt%
        caustic_demand = 0.02 * self.F_mass_in
        caustic_supplied = sum(i.imass['NaOH'] for i in self.ins)
        caustic_needed = caustic_demand - caustic_supplied
        if caustic_needed > 0:
            caustic.imass['NaOH'] = caustic_needed
        else:
            caustic.empty()
        
        mixture = liquor.copy()
        mixture.mix_from(self.ins)
        self.T_in = mixture.T
        self.deconstruction_rxns(mixture)
        
        # removed
        solids.copy_flow(mixture, insolubles, remove=True)
        # Based on stream 713 in ref [2]
        solids.imass['Water'] = 0
        solids.imass['Water'] = solids.F_mass
        mixture.imass['Water'] -= solids.imass['Water']
        
        pulp.copy_flow(mixture)
        pulp.T = solids.T = self.T
        pulp.P = self.P
        
    def _design(self):
        duty = self.Hnet
        self.heat_utilities[0](duty=duty, T_in=self.T_in)
        self.design_results['Flow rate'] = self.F_mass_out


@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=236000, S=410369, CE=CEPCI[2009], n=0.7, BM=2)
class NeutralizationTank(Unit):
    _N_ins = 2
    _N_outs = 1
   
    def __init__(self, ID='', ins=None, outs=(), T=32+273.15):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.neutralization_rxn = Rxn('2 NaOH + H2SO4 -> Na2SO4 + 2 H2O', 
                                      reactant='NaOH', X=1)
    
    def _run(self):
        pulp, acid = self.ins
        mixture = self.outs[0]
        
        acid.imol['H2SO4'] = 0.5*pulp.imol['NaOH'] / self.neutralization_rxn.X
        acid.imass['Water'] = acid.imass['H2SO4'] / 0.93 * 0.07 # 93% purity
        
        mixture.mix_from(self.ins)
        self.neutralization_rxn(mixture)
        self.design_results['Flow rate'] = self.outs[0].F_mass
        
        
# Including seed fermenter, main fermenter, and surge tank
@cost(basis='Flow rate', ID='First seed fermenter', units='kg/hr',
      cost=138000, S=162593, CE=CEPCI[2009], n=1, BM=1.8)
@cost(basis='Flow rate', ID='First seed agitator', units='kg/hr',
      kW=1.677825, cost=10260, S=162593, CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Flow rate', ID='Second seed fermenter', units='kg/hr',
      cost=172500, S=162593, CE=CEPCI[2009], n=1, BM=1.8)
@cost(basis='Flow rate', ID='Second seed agitator', units='kg/hr',
      kW=17.8968, cost=33000, S=162593, CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Flow rate', ID='Bubble seed fermenter', units='kg/hr',
      cost=822300, S=162593, CE=CEPCI[2014], n=1, BM=2.3)
@cost(basis='Flow rate', ID='Seed circulation cooler', units='kg/hr',
      cost=25200, S=162593, CE=CEPCI[2014], n=1, BM=2.2)
@cost(basis='Flow rate', ID='Main fermenter', units='kg/hr',
      cost=28753800, S=162593, CE=CEPCI[2014], n=1, BM=2.3)
@cost(basis='Flow rate', ID='Main fermenter circulation cooler', units='kg/hr',
      cost=817700, S=162593, CE=CEPCI[2014], n=1, BM=2.2)
@cost(basis='Flow rate', ID='Main fermenter circulation pump', units='kg/hr',
      cost=195500, S=162593, CE=CEPCI[2014], n=1, BM=2.3)
@cost(basis='Flow rate', ID='Fermentation air compressor', units='kg/hr',
      cost=1503850, S=162593, CE=CEPCI[2014], n=1, BM=1.6)
@cost(basis='Flow rate', ID='Fermentation air receiver', units='kg/hr',
      cost=119295, S=162593, CE=CEPCI[2014], n=1, BM=2)
@cost(basis='Flow rate', ID='Fermentation surge tank', units='kg/hr',
      cost=32421, S=162593, CE=CEPCI[2011], n=0.6, BM=2.5)
class MuconicFermentation(Unit):
    _N_ins = 7
    _N_outs = 2
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr',
             'Flow rate': 'kg/hr'}
    
    # Assumes 10% inoculum based on the size of the largest seed fermenter (100 m3)
    # and the main fermenter (1000 m3), also consistent with NREL's assumption for 
    # 2,3-BDO fermentation and past assumption on ethanol fermentation
    inoculum_ratio = 0.1
    
    titer_limit = (34.5+68.5) / 2 # in g/L (kg/m3)
    
    effluent_titer = 0

    def __init__(self, ID='', ins=None, outs=(), T=32+273.15, P=1.34*101325,
                 set_titer_limit=False):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        self.set_titer_limit = set_titer_limit
        
        self.seed_fermentation_rxns = ParallelRxn([
    #                           Reaction definition                       Reactant    Conversion
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 P_putidaGrow + 2.4 H2O',    'Glucose',    0.46),
    Rxn('Glucose + 1.94 O2 -> 0.74 MuconicAcid + 1.57 CO2 + 3.78 H2O',    'Glucose',    0.54),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 P_putidaGrow + 2H2O',        'Xylose',     0.46),
    Rxn('Xylose + 1.57 O2 -> 0.62 MuconicAcid + 1.26 CO2 + 3.13 H2O',     'Xylose',     0.54),    
    Rxn('Arabinose + 0.039 CSL + 0.015 DAP -> 5 P_putidaGrow + 2H2O',     'Arabinose',  0.46),
    Rxn('Arabinose + 1.57 O2 -> 0.62 MuconicAcid + 1.26 CO2 + 3.13 H2O',  'Arabinose',  0.54)
            ])
        self._seed_X = self.seed_fermentation_rxns.X.copy()
        
        self.main_fermentation_rxns = ParallelRxn([
    #                           Reaction definition                                 Reactant        Conversion
    Rxn('Glucose + 1.18 O2 + 0.28 NH4OH -> 4.8 P_putida + 1.2 CO2 + 2.26 H2O',      'Glucose',        0.46),
    Rxn('Glucose + 1.94 O2 -> 0.74 MuconicAcid + 1.57 CO2 + 3.78 H2O',              'Glucose',        0.54),
    Rxn('Xylose + 0.98 O2 + 0.23 NH4OH -> 4 P_putida + CO2 + 1.87 H2O',             'Xylose',         0.46),
    Rxn('Xylose + 1.57 O2 -> 0.62 MuconicAcid + 1.26 CO2 + 3.13 H2O',               'Xylose',         0.54),    
    Rxn('Arabinose + 0.98 O2 + 0.23 NH4OH -> 4 P_putida + CO2 + 1.87 H2O',          'Arabinose',      0.46),
    Rxn('Arabinose + 1.57 O2 -> 0.62 MuconicAcid + 1.26 CO2 + 3.13 H2O',            'Arabinose',      0.54),
    Rxn('Sucrose + 2.35 O2 + 0.56 NH4OH -> 9.6 P_putida + 2.4 CO2 + 3.52 H2O',      'Sucrose',        0.46),
    Rxn('Sucrose + 3.8731 O2 -> 1.48 MuconicAcid + 3.13 CO2 + 6.57 H2O',            'Sucrose',        0.54),
    Rxn('Acetate + 0.39 O2 + 0.093 NH4OH -> 1.6 P_putida + 0.4 CO2 + 0.753 H2O',    'Acetate',        1),
    Rxn('Extractives + 0.68 O2 + 0.28 NH4OH -> 4.8 P_putida + 1.2 CO2 + 2.26 H2O',  'Extractives',    0.46),
    Rxn('Extractives + 1.44 O2 -> 0.74 MuconicAcid + 1.57 CO2 + 3.78 H2O',          'Extractives',    0.54),
    # Used a mid point between experimentally proved titer (34.5 g/L in ref [5], ~0% conversion)
    # and the designed titer in ref [2] (68.5/L, ~80% conversion);
    # using 100% conversion in ref [2] would lead to a titer around 77 g/L
    # due to the higher sugar content in ethanol fermentation stillage compared to ref [2],
    # which fermented sugars for 2,3-butanediol
    Rxn('SolubleLignin + 3 O2 -> MuconicAcid + 2 CO2 + H2O',                        'SolubleLignin',  0.4)
    # Rxn('SolubleLignin + 3 O2 -> MuconicAcid + 2 CO2 + H2O',                        'SolubleLignin',  1)
            ])
        self._main_X = self.main_fermentation_rxns.X.copy()

        # Based on chemical usage in ref [2], only need to neutralized to the mono salt        
        self.neutralization_rxn = \
            Rxn('MuconicAcid + NaOH -> MonoSodiumMuconate + 2 H2O',
                reactant='MuconicAcid', X=1)


    def _run(self):
        substrate, water, ammonia, caustic, CSL, DAP, air = self.ins
        vent, broth = self.outs
        
        substrate_seed = substrate.copy()
        substrate_seed.mol = substrate.mol * self.inoculum_ratio
        substrate_main = substrate.copy()
        substrate_main.mol = substrate.mol - substrate_seed.mol
        
        # Assume the same CSL and DAP loading as R301 and R302
        substrate_seed.imass['CSL'] = 0.005 * substrate_seed.F_mass
        substrate_main.imass['CSL'] = 0.0025 * substrate_main.F_mass
        CSL.imass['CSL'] = substrate_seed.imass['CSL'] + substrate_main.imass['CSL']

        substrate_seed.imass['DAP'] = 0.67 * substrate_seed.F_vol
        substrate_main.imass['DAP'] = 0.33 * substrate_main.F_vol
        DAP.imass['DAP'] = substrate_seed.imass['DAP']+substrate_main.imass['DAP']

        # Based on stream 708 in ref [2], however in the text it is stated that 
        # the seed is diluted twofolds
        water.imass['Water'] = substrate_seed.F_mass
        water.imass['Water'] = max(80000, substrate_seed.F_mass)
        
        substrate_seed.mix_from([water, substrate_seed])
        self.seed_fermentation_rxns.force_reaction(substrate_seed)

        substrate_main.imass['P_putida'] = substrate_seed.imass['P_putidaGrow']
        substrate_seed.imass['P_putidaGrow'] = 0
                
        substrate_main.mix_from([substrate_seed, substrate_main])        
        self.main_fermentation_rxns.force_reaction(substrate_main)
        self.neutralization_rxn.force_reaction(substrate_main)
        
        ammonia.imass['NH4OH'] = max(0, -substrate_main.imass['NH4OH'])
        air.imass['O2'] = max(0, -substrate_main.imass['O2'])
        # Air mass ratio based on stream 703 in ref [2]
        air.imass['N2'] = air.imass['O2'] / 0.21 * 0.79
        caustic.imass['NaOH'] = max(0, -substrate_main.imass['NaOH'])
        for i in ('NH4OH', 'O2', 'NaOH'):
            substrate_main.imass[i] = 0        
        
        vent.copy_flow(substrate_main, 'CO2', remove=True)
        vent.imol['N2'] = air.imol['N2']
        broth.copy_flow(substrate_main)
        vent.T = broth.T = self.T
        
        # Avoid getting tiny negatives
        for i in broth.mol.nonzero()[0]:
            if broth.mol[i] < 0:
                broth.mol[i] = min(0, broth.mol[i]+1e-6)
                
        self.effluent_titer = compute_muconic_titer(broth)

    def _design(self):
        mixture = self.ins[0].copy()
        mixture.mix_from(self.ins)
        duty = self.Hnet
        self.heat_utilities[0](duty=duty, T_in=mixture.T)
        self.design_results['Flow rate'] = self.outs[1].F_mass


@cost(basis='Volumetric flow', ID='Membrane separator', units='m3/hr',
      # 1303 in gallon per minute (GPM)
      cost=2048000, S=1303*_GPM_2_m3hr, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Volumetric flow', ID='Ultrafilter', units='m3/hr',
      # Replacement of ultrafilter is 0.0297 $/$ cost
      cost=2048000*0.0297, S=1303*_GPM_2_m3hr, CE=CEPCI[2011], n=0.6, BM=2.5)
class MuconicMembrane(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'Volumetric flow': 'm3/hr'}
    
    def _run(self):
        broth = self.ins[0]
        liquid, solids = self.outs
        
        solids.imass[insolubles] = broth.imass[insolubles]
        # 4.3% loss based on Appendic C in ref [2]
        solids.imass['MonoSodiumMuconate'] = 0.043 * broth.imass['MonoSodiumMuconate']
        # Based on stream 713 in ref [2]
        solids.imass['Water'] = min(broth.imass['Water'], solids.F_mass)
        liquid.mol = broth.mol - solids.mol
        
        # Avoid getting tiny negatives
        for i in liquid.mol.nonzero()[0]:
            if liquid.mol[i] < 0:
                liquid.mol[i] = min(0, liquid.mol[i]+1e-6)
    
    def _design(self):
        self.design_results['Volumetric flow'] = self.F_vol_in


@cost(basis='Volumetric flow in', ID='Filter', units='m3/hr',
      # 1347 in gallon per minute (GPM)
      cost=345234, S=1347*_GPM_2_m3hr, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Volumetric flow crystal', ID='Crystallizer', units='m3/hr',
      # 190 in gallon per minute (GPM)
      cost=7104192, S=190*_GPM_2_m3hr, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Crystal flow', ID='Centrifuge', units='kg/hr',
      cost=327680, S=13403, CE=CEPCI[2011], n=0.6, BM=2.3)
@cost(basis='Dried crystal flow', ID='Dryer', units='kg/hr',
      cost=555008, S=11526, CE=CEPCI[2011], n=0.6, BM=2.6)
class MuconicCrystallizer(Unit):
    _N_ins = 2
    _N_outs = 2
    # Not included in ref [2], but added here for the cooling needs
    _N_heat_utilities = 1
    _units= {'Volumetric flow in': 'm3/hr',
             'Volumetric flow crystal': 'm3/hr',
             'Crystal flow': 'kg/hr',
             'Dried crystal flow': 'kg/hr',
             'Duty': 'kJ/hr'}

    def __init__(self, ID='', ins=None, outs=(), T=15+273.15):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        
        self.reacidification_rxn = \
            Rxn('MonoSodiumMuconate + 0.5 H2SO4 -> MuconicAcid + 0.5 Na2SO4',
                reactant='MonoSodiumMuconate', X=1)

    def _run(self):
        liquid, acid = self.ins
        water, crystal = self.outs
        
        acid.imol['H2SO4'] = 0.5 * liquid.imol['MonoSodiumMuconate']
        acid.imass['Water'] = acid.imass['H2SO4'] / 0.93 * 0.07 # 93% purity
        mixture = liquid.copy()
        mixture.mix_from(self.ins)
        self.T_in = mixture.T
        self.reacidification_rxn(mixture)
        
        # No information on salts and water carried to crystals, assumed none
        wet_crystal = crystal.copy()
        wet_crystal.imass['MuconicAcid'] = 0.988 * mixture.imass['MuconicAcid']
        # Based on equipment sizing in ref [2]
        wet_crystal.imass['Water'] = wet_crystal.imass['MuconicAcid']/10930*(11498-11930)
        self.design_results['Crystal flow'] = wet_crystal.F_mass
        
        crystal.imass['MuconicAcid'] = wet_crystal.imass['MuconicAcid']
        water.mol = mixture.mol - crystal.mol

        crystal.T = water.T = self.T
        crystal.phase = 's'
        
    def _design(self):
        duty = self.H_out - self.H_in
        self.heat_utilities[0](duty=duty, T_in=self.T_in)
        Design = self.design_results
        Design['Volumetric flow in'] = self.F_vol_in
        Design['Volumetric flow crystal'] = self.outs[1].F_vol
        Design['Dried crystal flow'] = self.outs[1].F_mass


@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=1317325, S=328984, CE=CEPCI[2011], n=0.6, BM=2.3)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=59.656, cost=63000, S=328984, CE=CEPCI[2009], n=0.6, BM=2.6)
# Did not simulate the salt as composition is undocumneted in ref [2]
@cost(basis='Salt flow', ID='Centrifuge', units='kg/hr',
      cost=327680, S=13403, CE=CEPCI[2011], n=0.6, BM=2.3)
class MuconicDissolution(Unit):
    _N_ins = 3
    _N_outs = 1
    _units= {'Flow rate': 'kg/hr',
             'Salt flow': 'kg/hr'}
    
    def _run(self):
        crystal, recycled_ethanol, fresh_ethanol = self.ins
        solution = self.outs[0]
        
        fresh_ethanol.imass['Ethanol'] = 4*crystal.F_mass - recycled_ethanol.imass['Ethanol']
        solution.mix_from(self.ins)
        solution.phase = 'l'

    def _design(self):
        self.design_results['Flow rate'] = self.F_mass_in
        # Based on equipment sizing in ref [2]
        self.design_results['Salt flow'] = 149/54079 * self.F_mass_in
        

@cost(basis='Muconic flow', ID='Feed tank', units='kg/hr',
      cost=16721, S=53930, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Muconic flow', ID='Pump', units='kg/hr',
      cost=271926, S=53930, CE=CEPCI[2014], n=0.8, BM=1.4)
@cost(basis='Muconic flow', ID='Feed economizer', units='kg/hr',
      cost=101886, S=53930, CE=CEPCI[2014], n=0.7, BM=2.66)
@cost(basis='Liquid volumetric flow', ID='Reactor', units='m3/hr',
      # 66551 is the liquid volumetric flow rate in L/hr
      cost=6826620, S=66551/1e3, CE=CEPCI[2011], n=0.7, BM=2)
@cost(basis='Liquid volumetric flow', ID='Intercooler #1', units='m3/hr',
      cost=356676, S=66551/1e3, CE=CEPCI[2007], n=0.65, BM=2.21)
@cost(basis='Liquid volumetric flow', ID='Intercooler #2', units='m3/hr',
      cost=484357, S=66551/1e3, CE=CEPCI[2007], n=0.65, BM=2.21)
@cost(basis='H2 flow', ID='H2 makeup compressor', units='kg/hr',
      cost=1661679, S=406, CE=CEPCI[2011], n=0.6, BM=1.09)
@cost(basis='H2 flow', ID='H2 makeup compressor spare', units='kg/hr',
      cost=1661679, S=406, CE=CEPCI[2011], n=0.6, BM=1.08)
@cost(basis='Flow rate', ID='Hot high-pressure separator', units='kg/hr',
      cost=197681, S=54336, CE=CEPCI[2013], n=1, BM=1.5)
@cost(basis='Flow rate', ID='Hot gas cooler', units='kg/hr',
      cost=10980, S=54336, CE=CEPCI[2011], n=0.7, BM=1.66)
@cost(basis='Flow rate', ID='Cold high-pressure separator', units='kg/hr',
      cost=4852, S=54336, CE=CEPCI[2011], n=0.7, BM=2.59)
class MuconicHydrogenation(Unit):
    _N_ins = 2
    _N_outs = 1
    _N_heat_utilities = 1
    _units= {'Muconic flow': 'kg/hr',
             'Liquid volumetric flow': 'm3/hr',
             'H2 flow': 'kg/hr',
             'Flow rate': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), T=78+273.15, P=40*101325):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        self.hydrogenation_rxn = Rxn('MuconicAcid + H2 -> AdipicAcid', 
                                      reactant='MuconicAcid', X=1)

    def _run(self):
        muconic, hydrogen = self.ins
        adipic = self.outs[0]

        # Though H2:muconic == 2.6 on a molar basis, the extra H2 is recycled
        hydrogen.imol['H2'] = muconic.imol['MuconicAcid']
        
        adipic.mix_from(self.ins)
        adipic_cold = adipic.copy()
        adipic.T = self.T
        adipic.P = self.P
        duty = adipic.H - adipic_cold.H
        self.heat_utilities[0](duty=duty, T_in=adipic.T)
        self.hydrogenation_rxn(adipic)
        
    def _design(self):
        Design = self.design_results
        Design['Muconic flow'] = self.ins[0].F_mass
        Design['Liquid volumetric flow'] = self.F_vol_out
        Design['H2 flow'] = self.ins[1].F_mass
        Design['Flow rate'] = self.F_mass_out
    
    # Ru/C catalyst cost not included in capital/variable operating cost in ref [2],
    # but was stated in text/table
    def _cost(self):
        self._decorated_cost()
        self._BM['Ru/C catalyst'] = 1
        # WHSV (feed mass flow/catalyst mass=5/hr)
        self.purchase_costs['Ru/C catalyst'] = self.ins[0].F_mass/5 \
            * price['Hydrogenation catalyst']


@cost(basis='Duty', ID='Feed heater', units='kJ/hr',
      # 13 is the duty in MMkca/hr
      cost=274818, S=13*_Gcal_2_kJ, CE=CEPCI[2011], n=0.6, BM=3)
@cost(basis='Flow rate', ID='Feed tank', units='kg/hr',
      cost=45966, S=290932, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Flow rate', ID='Flash drum', units='kg/hr',
      cost=511000, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
class AdipicEvaporator(Unit):
    _N_ins = 2
    _N_outs = 2
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr',
             'Flow rate': 'kg/hr'}
    
    def _run(self):
        adipic_dilute, adipic_recycled = self.ins
        ethanol, adipic_concentrated = self.outs

        adipic_total = adipic_dilute.copy()
        adipic_total.mix_from(self.ins)
        adipic_concentrated.imass['Ethanol'] = 2.5 * adipic_total.imass['AdipicAcid']
        ethanol.imass['Ethanol'] = adipic_total.imass['Ethanol'] \
            - adipic_concentrated.imass['Ethanol']
        adipic_concentrated.mol = adipic_total.mol - ethanol.mol
            
        ethanol.phase = 'g'
        ethanol.T = adipic_concentrated.T = adipic_total.T

        duty = self.design_results['Duty'] = self.H_out - self.H_in
        self.heat_utilities[0](duty=duty, T_in=adipic_total.T)
        self.design_results['Flow rate'] = self.F_mass_in
        

@cost(basis='Volumetric flow', ID='Crystallizer', units='m3/hr',
      # 190 in gallon per minute (GPM)
      cost=7104192, S=190*_GPM_2_m3hr, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Crystal flow', ID='Centrifuge', units='kg/hr',
      cost=327680, S=13403, CE=CEPCI[2011], n=0.6, BM=2.3)
class AdipicCrystallizer(Unit):
    _N_ins = 1
    _N_outs = 2
    # Not included in ref [2], but added here for the cooling needs
    _N_heat_utilities = 1
    _units= {'Volumetric flow': 'm3/hr',
             'Crystal flow': 'kg/hr'}

    def __init__(self, ID='', ins=None, outs=(), T=15+273.15, P=101325):	
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        
    def _run(self):
        influent = self.ins[0]
        uncrystallized, crystal = self.outs
        
        crystal.imass['AdipicAcid'] = 0.734 * influent.imass['AdipicAcid']
        # Based on stream 710 in ref [2]
        crystal.imass['Ethanol'] = 29/11092 * crystal.imass['AdipicAcid']
        uncrystallized.mol = influent.mol - crystal.mol
        
        crystal.T = uncrystallized.T = self.T
        crystal.phase = 's'
        
    def _design(self):
        Design = self.design_results
        Design['Volumetric flow'] = self.outs[1].F_vol
        Design['Crystal flow'] = self.outs[1].F_mass
        duty = self.H_out - self.H_in
        self.heat_utilities[0](duty=duty, T_in=self.ins[0].T)
        
@cost(basis='Duty', ID='Condenser', units='kJ/hr',
      # -23 is the duty in MMkca/hr
      cost=487000, S=-23*_Gcal_2_kJ, CE=CEPCI[2010], n=0.6, BM=2.8)
class AdipicCondenser(HXutility):
    def _design(self):
        super()._design()
        self.design_results.clear()
        self.design_results['Duty'] = self.Q


# %% 

# =============================================================================
# Wastewater treatment
# =============================================================================

# Cost scaling modified from ref [1] (flow rate basis) to ref [4] (on COD basis)
@cost(basis='COD flow', ID='Anaerobic basin', units='kg-O2/hr',
      kW=2371.326, cost=25800000, S=27211, CE=CEPCI[2012], n=0.6, BM=1.1)
class AnaerobicDigestion(Unit):
    _N_ins = 1	
    _N_outs = 3
    _units= {'COD flow': 'kg-O2/hr'}
    
    auxiliary_unit_names = ('heat_exchanger',)
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants, split=(), T=35+273.15):
        Unit.__init__(self, ID, ins, outs)	
        self.reactants = reactants	
        self.split = split	
        self._multi_stream = tmo.MultiStream(None)	
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.heat_utilities = hx.heat_utilities
        chems = self.chemicals	
        	
        # Based on P49 in ref [1], 91% of organic components is destroyed,	
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

        self.sulfate_rxns = ParallelRxn([
            #   Reaction definition                           Reactant    Conversion
            Rxn('AmmoniumSulfate -> 2 NH3 + H2S + 2 O2',  'AmmoniumSulfate',  1),
            Rxn('H2SO4 -> H2S + 2 O2',                    'H2SO4',            1)
            ])
              	
    def _run(self):	
        wastewater = self.ins[0]	
        biogas, treated_water, sludge = self.outs	
        T = self.T	

        sludge.copy_flow(wastewater)	
        self.digestion_rxns(sludge.mol)
        self.sulfate_rxns(sludge.mol)

        ms = self._multi_stream
        ms.copy_flow(sludge)    
        ms.vle(P=101325, T=T)   
        biogas.mol = ms.imol['g']   
        biogas.phase = 'g'  
        liquid_mol = ms.imol['l']   
        treated_water.mol = liquid_mol * self.split	
        sludge.mol = liquid_mol - treated_water.mol	
        biogas.receive_vent(treated_water, accumulate=True)	
        biogas.T = treated_water.T = sludge.T = T
        
    def _design(self):
        wastewater = self.ins[0]
        self.design_results['COD flow'] = compute_COD(COD_chemicals, self.ins[0])
        # Calculate utility needs to keep digester temperature at 35°C,	
        # heat change during reaction is not tracked	
        H_at_35C = wastewater.thermo.mixture.H(mol=wastewater.mol, 	
                                               phase='l', T=self.T, P=101325)	
        duty = -(wastewater.H - H_at_35C)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, wastewater)


@cost(basis='COD flow', ID='Ammonia addition', units='kg-O2/hr',
      kW=13.4226, cost=195200, S=5600, CE=CEPCI[2012], n=0.6, BM=1.5)
@cost(basis='COD flow', ID='Caustic feed', units='kg-O2/hr',
      kW=4.4742, cost=20000, S=5600, CE=CEPCI[2012], n=0.6, BM=3)
@cost(basis='Volumetric flow', ID='Aerobic basin', units='m3/hr',
      # power usage including polymer addition and feed pumping
      kW=1.4914+134.226,
      # 2.7 in million gallons per day (MGD)
      cost=4804854, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=0.6, BM=2.1)
@cost(basis='COD flow', ID='Blowers', units='kg-O2/hr',
      kW=6711.3, cost=2070000, S=5600, CE=CEPCI[2012], n=0.6, BM=2)
class AerobicDigestion(Unit):    
    _N_ins = 6
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr',
             'Volumetric flow': 'm3/hr'}
    
    # 4350 and 356069 are water flows from streams 622 and 611 in ref [1]
    evaporation = 4350 / 356069
    
    def __init__(self, ID='', ins=None, outs=(), *, reactants,
                 caustic_mass=0, need_ammonia=False):
        Unit.__init__(self, ID, ins, outs)
        self.reactants = reactants
        self.caustic_mass = caustic_mass
        self.need_ammonia = need_ammonia
        
        chemicals = self.chemicals
        def growth(reactant):
            f = chemicals.WWTsludge.MW / getattr(chemicals, reactant).MW 
            return Rxn(f"{f}{reactant} -> WWTsludge", reactant, 1.)
        
        # Reactions from auto-populated combustion reactions.
        # Based on P49 in ref [1], 96% of remaining soluble organic matter 
        # is removed after aerobic digestion, of which 74% is converted to
        # water and CO2 and 22% to cell mass
        combustion_rxns = chemicals.get_combustion_reactions()
        
        self.digestion_rxns = ParallelRxn([i*0.74 + 0.22*growth(i.reactant)
                                            for i in combustion_rxns
                                            if (i.reactant in reactants)])
        self.digestion_rxns.X[:] = 0.96
        
        #                                 Reaction definition       Reactant Conversion
        self.nitrification_rxn = Rxn('NH4OH + 2 O2 -> HNO3 + 2 H2O', 'NH4OH',  1)
        
        self.neutralization_rxns = ParallelRxn([
            #              Reaction definition       Reactant Conversion
            Rxn('H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O',  'H2SO4',    1),
            Rxn('HNO3 + NaOH -> NaNO3 + H2O',        'HNO3',     1)
            ])
    
    def _run(self):
        influent, recycle, caustic, ammonia, polymer, air = self.ins
        vent, effluent = self.outs
        vent.phase = 'g'

        caustic.imass['NaOH'] = self.caustic_mass
        # Ammonia as a nutrient
        if self.need_ammonia is True:
            # Based on Table 33 on Page 73 of ref [2], originally as NH3
            ammonia.imass['NH4OH'] = 36 * 35.046/17.031
        else: ammonia.empty()

        effluent.mix_from(self.ins[0:5])

        self.design_results['Volumetric flow'] = influent.F_vol
        self.nitrification_rxn.force_reaction(effluent.mol)
        self.neutralization_rxns.force_reaction(effluent.mol)
        if effluent.imass['NaOH'] < 0:
            caustic.imass['NaOH'] += -effluent.imol['NaOH']
            effluent.imol['NaOH'] = 0

        # Ratio based on equipment sizing in ref [2]
        ratio = self.design_results['Volumetric flow'] / (2*_MGD_2_m3hr)
        polymer.imass['Polymer'] = 2 * ratio

        # 4693, 54718 and 180206 from stream 601 in ref [2]
        air.imass['O2'] = 54718 * ratio
        air.imass['N2'] = 180206 * ratio
        effluent.mix_from([effluent, air])
        
        self.design_results['COD flow'] = compute_COD(COD_chemicals, effluent)
        self.digestion_rxns(effluent.mol)
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        vent.imol['Water'] = influent.imol['Water'] * self.evaporation
        effluent.imol['Water'] -= vent.imol['Water']
        vent.T = effluent.T

    def _design(self):
        self._decorated_cost()
        if self.need_ammonia == False:
            self.cost_items['Ammonia addition'].cost = 0


@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # power usage including pumps
      kW=63.3845+715.872+14.914,
      # 2.7 in million gallons per day (MGD)
      cost=4898500, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=1, BM=1.6)
@cost(basis='COD flow', ID='Conveyor', units='kg-O2/hr', 
      kW=7.457, cost=7000, S=5600, CE=CEPCI[2012], n=0.6, BM=2.9)
class MembraneBioreactor(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'Volumetric flow': 'm3/hr',
             'COD flow': 'kg-O2/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), *, split):
        Unit.__init__(self, ID, ins, outs)
        self.split = split

    def _run(self):        
        mixture = self.ins[0].copy()
        mixture.mix_from(self.ins)
        separations.split(mixture, *self.outs, self.split)

    def _design(self):
        self.design_results['Volumetric flow'] = self.outs[0].F_vol
        self.design_results['COD flow'] = compute_COD(COD_chemicals, self.ins[0])


@cost(basis='COD flow', ID='Thickeners', units='kg-O2/hr', 
      kW=107.3808, cost=750000, S=5600, CE=CEPCI[2012], n=0.6, BM=1.6)
class BeltThickener(Unit):
    _ins_size_is_fixed = False
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}
    
    def _run(self):
        centrate, solids = self.outs
        
        influent = self.ins[0].copy()
        influent.mix_from(self.ins)
        solids.copy_flow(influent, insolubles)
        # Concentrate sludge to 4% solids
        solids.imass['Water'] = 0.96/0.04 * influent.imass[insolubles].sum()
        if solids.imass['Water'] < influent.imass['Water']:
            ratio = solids.imass['Water'] / influent.imass['Water']
            solids.imass[solubles] = ratio * influent.imass[solubles]
            solids.T = influent.T
            
            centrate.mol = influent.mol - solids.mol
            centrate.T = influent.T
        else:
            centrate.empty()
            solids.copy_like(influent)

    def _design(self):
        self.design_results['COD flow'] = compute_COD(COD_chemicals, self.ins[0])

@cost(basis='COD flow', ID='Centrifuge', units='kg-O2/hr',
      # power usage including feed pumping and centrifuge
      kW=22.371+123.0405, cost=686800, S=5600, CE=CEPCI[2012], n=0.6, BM=2.7)
class SludgeCentrifuge(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}
    
    def _run(self):
        influent = self.ins[0]
        centrate, solids = self.outs

        # Centrifuge captures 95% of the solids at 20% solids
        solids.imass[insolubles] = 0.95 * influent.imass[insolubles]
        solids.imass['Water'] = 0.8/0.2 * (influent.imass[insolubles].sum())
        if solids.imass['Water'] < influent.imass['Water']:
            ratio = solids.imass['Water'] / influent.imass['Water']
            solids.imass[solubles] = ratio * influent.imass[solubles]
            solids.T = influent.T
            
            centrate.mol = influent.mol - solids.mol    
            centrate.T = influent.T
        else:
            centrate.empty()
            solids.copy_like(influent)
    
    def _design(self):
        self.design_results['COD flow'] = compute_COD(COD_chemicals, self.ins[0])

@cost(basis='Volumetric flow', ID='Reactor', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      cost=2450000, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=1, BM=1.8)
@cost(basis='Volumetric flow', ID='Evaporator', units='m3/hr',
      # 2.7 in million gallons per day (MGD)
      kW=1103.636, cost=5000000, S=2.7*_MGD_2_m3hr, CE=CEPCI[2012], n=0.6, BM=1.6)
class ReverseOsmosis(Unit):
    _N_ins = 1
    _N_outs = 2
    _units = {'Volumetric flow': 'm3/hr'}
    
    def _run(self):
        influent = self.ins[0]
        water, brine = self.outs
        
        self.design_results['Volumetric flow'] = self.F_vol_in

        # Based on stream 626 and 627 in ref [1]
        water.imass['Water'] = 376324/(376324+4967) * influent.imass['Water']
        brine.mol = influent.mol - water.mol
        water.T = brine.T = influent.T


@cost(basis='Duty', ID='Feed heater', units='kJ/hr',
      # 13 is the duty in MMkca/hr
      cost=274818, S=13*_Gcal_2_kJ, CE=CEPCI[2011], n=0.6, BM=3)
@cost(basis='Flow rate', ID='Feed tank', units='kg/hr',
      cost=45966, S=290932, CE=CEPCI[2011], n=0.6, BM=2.5)
@cost(basis='Flow rate', ID='Flash drum', units='kg/hr',
      cost=511000, S=264116, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Crude salt flow', ID='Centrifuge', units='kg/hr',
      cost=327680, S=11524, CE=CEPCI[2011], n=0.6, BM=2.3)
@cost(basis='Salt flow', ID='Dryer', units='kg/hr',
      cost=555008, S=11524, CE=CEPCI[2011], n=0.6, BM=2.6)
class SodiumSulfateRecovery(Unit):
    _N_ins = 1
    _N_outs = 3
    _N_heat_utilities = 1
    _units= {'Duty': 'kJ/hr',
             'Flow rate': 'kg/hr',
             'Crude salt flow': 'kg/hr',
             'Salt flow': 'kg/hr'}

    decomposition_rxn = Rxn('NH4OH -> NH3 + H2O', reactant='NH4OH', X=1)

    def _run(self):
        brine = self.ins[0]
        vent, residuals, Na2SO4 = self.outs
        
        influent = brine.copy()
        self.decomposition_rxn(influent.mol)
        Na2SO4.copy_flow(influent, 'Na2SO4')
        residuals.copy_flow(influent, solubles+insolubles)
        vent.mol = influent.mol - residuals.mol
        vent.T = Na2SO4.T = residuals.T = vent.dew_point_at_P(101325).T
        Na2SO4.phase = residuals.phase = 'l'
        vent.phase = 'g'
        Na2SO4.phase = 's'
        self.design_results['Duty'] = self.H_out - self.H_in
    
    def _design(self):
        Design = self.design_results
        self.heat_utilities[0](duty=self.design_results['Duty'], T_in=self.ins[0].T)
        Design['Flow rate'] = self.F_mass_in
        Design['Salt flow'] = (self.F_mass_out-self.outs[0].F_mass)
        # Based on equipment sizing in ref [2]
        Design['Crude salt flow'] = 14871/14163 * Design['Salt flow']




# %% 

# =============================================================================
# Storage
# =============================================================================

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=1340000, S=22681, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=9200, S=22681, CE=CEPCI[2009], n=0.8, BM=3.1)
class EthanolStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=200000, S=473, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=473, CE=CEPCI[2009], n=0.8, BM=3.1)
class DenaturantStorage(Unit): pass

@cost(basis='Flow rate', ID='In-line mixer', units='kg/hr',
      cost=3850, S=23154, CE=CEPCI[2009], n=0.5, BM=1)
class DenaturantMixer(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def _run(self):
        # Based on streams 701 and 515 in ref [1]
        self.ins[1].imass['Denaturant'] = 465/21673 * self.ins[0].imass['Ethanol']
        self.outs[0].mix_from(self.ins)

# Adipic acid and sodium sulfate
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=690900, S=1981, CE=CEPCI[2007], n=0.65, BM=1.85)
class CoproductStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2010], n=0.7, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=7493, S=1981, CE=CEPCI[2010], n=0.8, BM=2.3)
class SulfuricAcidStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      # Original size basis for NH3 instead of NH4OH
      cost=196000, S=1171/17.031*35.046, CE=CEPCI[2010], n=0.7, BM=2)
class AmmoniaStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2011], n=0.7, BM=1.5)
class CausticStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=70000, S=1393, CE=CEPCI[2009], n=0.7, BM=2.6)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21200, S=1393, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=1393, CE=CEPCI[2009], n=0.8, BM=3.1)
class CSLstorage(Unit): pass

@cost(basis='Flow rate', ID='Unloader', units='kg/hr',
      cost=30000, S=163, CE=CEPCI[2009], n=0.6, BM=1.7)
@cost(basis='Flow rate', ID='Make-up tank', units='kg/hr',
      cost=81192, S=142, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=4.10135, cost=9800, S=163, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=163, CE=CEPCI[2009], n=0.8, BM=3.1)
class DAPstorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=803000, S=8343, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=15000, S=8343, CE=CEPCI[2009], n=0.8, BM=3.1)
class FirewaterStorage(Unit): pass




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp
from warnings import warn
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank, LiquidsSplitSettler
from biosteam.units.decorators import cost
from thermosteam import Stream, MultiStream
from biorefineries.make_a_biorefinery.process_settings import price
from biorefineries.make_a_biorefinery.utils import CEPCI, baseline_feedflow, compute_extra_chemical, adjust_recycle
from biorefineries.make_a_biorefinery.chemicals_data import chems
tmo.settings.set_thermo(chems)

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
class FeedstockPreprocessing(Unit): # not used
    # 2205 U.S. ton/day (2000 metric tonne/day) as in Humbird et al.
    _baseline_flow_rate = baseline_feedflow.sum()
    _cached_flow_rate = 2205

   
# %% 

# =============================================================================
# Pretreatment
# =============================================================================
#!!! Replace these fake unit classes & rxns with real pretreatment process unit classes & rxns

@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=CEPCI[2010], n=0.6, BM=1)
class FakePretreatment(Unit):
    def _run(self):
        feedstock, = self.ins
        pretreated_slurry, = self.outs
        pretreated_slurry.copy_like(feedstock)


# %% 

# =============================================================================
# Conversion
# =============================================================================
#!!! Replace these fake unit classes & rxns with real conversion process unit classes & rxns


@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=CEPCI[2010], n=0.6, BM=1)
class FakeSaccharification(Unit):

    _N_ins = 3
    _N_outs = 1
    
    fake_hydrolysis_rxns = ParallelRxn([
        Rxn('Glucan + H2O -> Glucose', 'Glucan', 1.),
        Rxn('Xylan + H2O -> Xylose', 'Xylan', 1.),
        ])
    
    def _run(self):
        pretreated_slurry, water, enzyme = self.ins
        saccharified_slurry, = self.outs
        
        water.imol['Water'] = 5*(pretreated_slurry.imol['Glucan', 'Xylan'].sum())
        
        saccharified_slurry.mix_from([pretreated_slurry, water])
        self.fake_hydrolysis_rxns(saccharified_slurry)
        
        enzyme_mass = 0.001 * saccharified_slurry.F_mass
        enzyme.imass['Enzyme'] = enzyme_mass
        saccharified_slurry.imass['Enzyme'] = enzyme_mass
        
@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=CEPCI[2010], n=0.6, BM=1)
class FakeCoFermentation(Unit):
    
    _N_ins = 1
    _N_outs = 2
    
    fake_cofermentation_rxns = ParallelRxn([
        Rxn('Glucose -> SuccinicAcid', 'Glucose', 0.5),
        Rxn('Xylose -> SuccinicAcid', 'Xylose', 0.5),
        
        Rxn('Glucose -> LacticAcid', 'Glucose', 0.1),
        Rxn('Xylose -> LacticAcid', 'Xylose', 0.15),
        
        Rxn('Glucose -> Ethanol', 'Glucose', 0.1),
        Rxn('Xylose -> Ethanol', 'Xylose', 0.1),
        
        Rxn('Glucose -> FermMicrobe', 'Glucose', 0.05),
        Rxn('Xylose -> FermMicrobe', 'Xylose', 0.05),
        ])
    
    fake_CO2_generation_rxns = ParallelRxn([
        Rxn('Glucose -> 6 CO2 + 6H2O',       'Glucose',   1.-1e-9),
        Rxn('Xylose -> 5 CO2 + 5H2O',        'Xylose',    1.-1e-9),
        ])
    
    def _run(self):
        saccharified_slurry, = self.ins
        fermentation_broth, vented_gases = self.outs
        
        fermentation_broth.copy_like(saccharified_slurry)
        self.fake_cofermentation_rxns(fermentation_broth)
        self.fake_CO2_generation_rxns(fermentation_broth)
        
        CO2_vented = fermentation_broth.imol['CO2']
        vented_gases.imol['CO2'] = CO2_vented
        vented_gases.phase = 'g'
        fermentation_broth.imol['CO2'] = 0.
        
# %% 

# =============================================================================
# Separation
# =============================================================================


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
        self.isplit = isplit = self.thermo.chemicals.isplit(split, None)
        self.split = isplit.data
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
        # biogas.receive_vent(treated_water, accumulate=True)	
        biogas.receive_vent(treated_water)
        biogas.T = treated_water.T = sludge.T = T
        
    def _design(self):
        wastewater = self.ins[0]
        # Calculate utility needs to keep digester temperature at 35°C,	
        # heat change during reaction is not tracked	
        H_at_35C = wastewater.thermo.mixture.H(mol=wastewater.mol, 	
                                               phase='l', T=self.T, P=101325)	
        duty = -(wastewater.H - H_at_35C)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, wastewater)
  
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

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      # Original size basis for NH3 instead of NH4OH
      cost=196000, S=1171/17.031*35.046, CE=CEPCI[2010], n=0.7, BM=2)
class AmmoniaStorageTank(Unit): pass

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



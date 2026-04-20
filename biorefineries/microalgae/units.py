#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat June 23 17:41:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- units

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/API/thermosteam/Chemicals.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP

@author: Xingdong Shi
@version: 0.0.1
"""

import numpy as np
import biosteam as bst
import math
from biosteam import Unit
from biosteam.units.design_tools import pressure_vessel_material_factors
from biosteam.units.decorators import cost
from biosteam.units import HXutility, Mixer, MixTank, Pump, StorageTank, StirredTankReactor
from thermosteam import Stream
import thermosteam as tmo
from ._chemicals import chems
from .utils import baseline_feedflow, CEPCI
from thermosteam import reaction as rxn
from biorefineries.lactic._units import Reactor

def anaerobic_rxn(reactant):
    MW = getattr(chems, reactant).MW
    return rxn.Reaction(f"{1/MW}{reactant} -> {P_ch4}CH4 + {P_co2}CO2 + {P_sludge}WWTsludge", reactant, 0.91)

def find_split(*tuples):

    split = {}
    for ID, flow0, flow1 in tuples:
        total = flow0 + flow1
        split[ID] = flow0 / total if total else 0
    return split

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
ln = math.log
exp = math.exp
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
chems = tmo.settings.chemicals
organics = ['Glucose', 'Xylose', 'Protein', 'Ash', 'AceticAcid', 'Ethanol']
if 'WWTsludge' in organics:
    organics.remove('WWTsludge')
P_sludge = 0.05/0.91/chems.WWTsludge.MW
MW = np.array([chems.CH4.MW, chems.CO2.MW])
CH4_molcomp = 0.60
mass = np.array([CH4_molcomp, 1-CH4_molcomp])*MW
mass /= mass.sum()
mass *= 0.381/(0.91)
P_ch4, P_co2 = mass/MW
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics])

splits = [
    ('Glucose', 3, 42),
    ('Xylose', 7, 85),
    ('Protein', 51, 18),
    ('Ash', 813, 280),
    ('AceticAcid', 0, 5),
    ('Ethanol', 1, 15),
]

tmo.settings.set_thermo(chems)

# Compute caproic acid titer (g/L) from an effluent stream, consistent with
# succinic project's approach of post-calculation based on effluent composition
def compute_caproic_acid_titer(effluent):
    try:
        vol = effluent.F_vol
        return (effluent.imass['CaproicAcid'] / vol) if vol else 0.0
    except Exception:
        return 0.0

# %% 
# =========================
# Microalgae Crushing 
# =========================
@cost(basis='Flow rate', ID='System', units='kg/hr',
      kW=511.3205, cost=13329690, S=94697, CE=CEPCI[2009], n=0.6, BM=1.7)
class FeedstockPreprocessing(Unit):
    _baseline_flow_rate = baseline_feedflow.sum()
    _cached_flow_rate = 2205
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        self.outs[0].copy_like(self.ins[0])
        s = self.outs[0]
        for comp in s.chemicals.IDs:
            if comp != 'Microalgae':
                s.imass[comp] = 0

# %% 
# =============================================================================
# Pretreatment
# =============================================================================
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=6210, S=1981, CE=CEPCI[2010],  n=0.7, BM=3)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      cost=8000, S=3720, CE=CEPCI[2009],  n=0.8, BM=2.3)
class SulfuricAcidAdditionTank(StorageTank):
    _N_ins = 1
    _N_outs = 1
    # acid_volume_per_biomass = 0.8 mL/g
    # acid_density = 1.84  # g/mL
    acid_loading = 1.47  # g H2SO4 / g microalgae

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        microalgae_mass = feed.imass['Microalgae']
        
        pure_H2SO4 = microalgae_mass * self.acid_loading
        
        acid_purity = 0.93
        acid_solution_mass = pure_H2SO4 / acid_purity
        water_mass = acid_solution_mass * (1 - acid_purity)
        out.copy_like(feed)
        out.imass['H2SO4'] = pure_H2SO4
        out.imass['H2O'] = water_mass

# Sulfuric acid in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      cost=6000, S=136260, CE=CEPCI[2009], n=0.5, BM=1)
class SulfuricAcidMixer(Unit):
    _N_ins = 2  
    _N_outs = 1 
    def _run(self):
        microalgae, acid = self.ins
        out = self.outs[0]
        out.mix_from([microalgae, acid])

# Pretreatment reactor
@cost(basis='Dry flow rate', ID='Pretreatment reactor', units='kg/hr',
      kW=5120, cost=19812400, S=83333, CE=CEPCI[2009], n=0.6, BM=1.5)
class AcidPretreatmentReactor(Unit):
    _N_ins = 1  
    _N_outs = 1 
    T = 121 + 273.15  # K
    t = 20 / 60  # h

    microalgae_decomposition = tmo.reaction.Reaction(
        'Microalgae -> 0.5 Protein + 0.159 Lipid + 0.234 Carbohydrate + 0.107 Ash', 
        reactant='Microalgae',
        X=1.0,  
        basis='wt'
    )

    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        self.microalgae_decomposition(out)
        
# pH adjustment
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
	  # Size basis on the total flow, not just ammonia, 	
      # thus assuming difference caused by MWs of NH3 and NH4OH is negligible
      cost=50000, S=157478, CE=CEPCI[2009], n=0.5, BM=1.5)
class NAOHMixer(Mixer):
    def _run(self):
        NAOH, water = self.ins
        mixture = self.outs[0]	
        	
        # Make 10% solution	
        water.imass['Water'] = NAOH.imass['NAOH'] * (1-0.1)/0.1	
        mixture.mix_from([NAOH, water])


@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=236000, S=410369, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21900, S=410369, CE=CEPCI[2009], n=0.5, BM=1.5)
class NAOHAdditionTank(MixTank):
    # Based on experimental data: 0.64kg NaOH consumed per ton of microalgae
    naoh_mass_per_ton_microalgae = 0.64  # kg NaOH per ton microalgae
    
    def _run(self):
        # Calculate required NaOH based on microalgae mass
        if hasattr(self, 'biomass_mass') and self.biomass_mass > 0:
            required_naoh = self.biomass_mass * self.naoh_mass_per_ton_microalgae / 1000  # convert to kg/hr
            # Set NaOH input amount
            if len(self.ins) > 0:
                self.ins[0].imass['NaOH'] = required_naoh
        
        #                                      Reaction definition      Reactant Conversion
        self.neutralization_rxn = Rxn('2 NaOH + H2SO4 -> Na2SO4 + 2 H2O', 'H2SO4', 1)
        
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
    
    # Experimental data: 6kg GlucoAmylase and 266.7kg AlphaAmylase consumed per ton of microalgae
    glucoamylase_per_ton = 6      # kg per ton microalgae
    alphaamylase_per_ton = 266.7  # kg per ton microalgae
    # enzyme_mass_per_ton_microalgae = 172.7  # kg enzyme per ton microalgae (deprecated)
    
    def _run(self):    
        hydrolysate, enzyme, water = self.ins    
        effluent = self.outs[0]    
        
        # Calculate consumption of two enzymes based on microalgae mass
        if hasattr(self, 'biomass_mass') and self.biomass_mass > 0:
            required_glucoamylase = self.biomass_mass * self.glucoamylase_per_ton / 1000  # kg/hr
            required_alphaamylase = self.biomass_mass * self.alphaamylase_per_ton / 1000  # kg/hr
            enzyme.imass['GlucoAmylase'] = required_glucoamylase
            enzyme.imass['AlphaAmylase'] = required_alphaamylase
        else:
            # If no biomass_mass attribute, use default method
            enzyme.imass['GlucoAmylase'] = hydrolysate.F_mass * 0.006  # 0.6%
            enzyme.imass['AlphaAmylase'] = hydrolysate.F_mass * 0.2667 # 26.67%
        
        effluent.mix_from([hydrolysate, enzyme, water])



@cost(basis='Saccharification tank size', ID='Saccharification tanks', units='kg',
      cost=3840000, S=421776*24, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Slurry flow rate', ID='Saccharification transfer pumps', units='kg/hr',
      kW=74.57, cost=47200, S=421776*24, CE=CEPCI[2009], n=0.8, BM=2.3)
class Saccharification(Unit):    
    _N_ins = 1
    _N_outs = 1
    _units= {'Saccharification tank size': 'kg',
             'Slurry flow rate': 'kg/hr'}             
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_saccharification = 3 
        
    def _init(self, ID='', ins=None, outs=(), T=50+273.15):
        
        self.T = T
        
        # Based on Table 9 on Page 28 of Humbird et al.
        #!!! Yalin updated the glucan-to-glucose conversion to match the table
        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Carbohydrate + H2O -> Glucose', 'Carbohydrate', 1)
            ])
        
        # self.saccharified_stream = Stream(None)

    def _run(self):
        feed = self.ins[0]
        ss= self.outs[0]
        ss.copy_like(feed)
        ss.T = self.T
        self.saccharification_rxns(ss.mol)
        neutralization_rxn = Rxn('2 NaOH + H2SO4 -> Na2SO4 + 2 H2O', 'NaOH', 1)
        neutralization_rxn(ss)


        
    def _design(self):
        Design = self.design_results
        total_mass_flow = self.ins[0].F_mass 
        Design['Saccharification tank size'] = total_mass_flow * self.tau_saccharification
        Design['Slurry flow rate'] = total_mass_flow

## MCCA Fermentation
@cost(basis='Fermenter size', ID='Fermenter', units='kg',
      kW=100, cost=5128000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=0.7, BM=1.5)
@cost(basis='Fermenter size', ID='Fermenter agitator', units='kg',
      kW=10,
      cost=630000, 
      S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=0.7, BM=1.5)

class MCCAFermentation(StirredTankReactor):
    _N_ins = 5  
    _N_outs = 4 

    auxiliary_unit_names = ('heat_exchanger',)
    _units= {
        **Reactor._units,
        'Fermenter size': 'kg',
        'Broth flow rate': 'kg/hr',
        }
    
    tau_batch_turnaround = 12 # in hr, seedmentaiton time
    tau_cofermentation = 15 * 24 # in hr, 15 d
    mol_NH4OH_per_acid_pH_control = 0.2988 # from IBRL batch run
    mol_lime_per_acid_pH_control = 0.2988/2 # mol NH4OH / 2
    pH = 5.0

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, T=37+273.15,
                 P=101325, V_wf=0.8, length_to_diameter=0.6,
                 kW_per_m3=0.02 * 0.7457 / 0.075, # ref from succinic projects
                 mixing_intensity=300,
                 #wall_thickness_factor=1,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 neutralization=True,
                 neutralizing_agent='Lime',
                 mode='batch', feed_freq=1,
                 pH_control=True,
                 base_neutralizes_product=True,
                 tau = tau_cofermentation,
                 microalgae_mass_flow: float | None = None,  # kg/hr basis of original algae for yield calc
                 titer: float = 2.003,  # g/L
                 # allow_dilution=False,
                 # allow_concentration=False,
                 # sugars=None
                 ):

        StirredTankReactor.__init__(self, ID, ins, outs)
        
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        #self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.neutralization = neutralization
        self.mode = mode
        self.feed_freq = feed_freq
        self.pH_control=pH_control
        self.base_neutralizes_product = base_neutralizes_product
        self.neutralization = neutralization
        self.neutralizing_agent = neutralizing_agent
        self.tau = tau 
 
        #self.allow_dilution = allow_dilution
        #self.allow_concentration = allow_concentration
        #self.sugars = sugars or tuple(i.ID for i in self.chemicals.sugars)

        # Store optional explicit biomass flow for yield calculations
        self.microalgae_mass_flow = microalgae_mass_flow
        self.titer = titer
        ID = self.ID
        self._mixed_feed = Stream(f'{ID}_mixed_feed')
        self._tot_feed = Stream(f'{ID}_tot_feed')
        # Before reaction, after reaction, with last feed
        self._single_feed0 = Stream(f'{ID}_single_feed0')
        self._single_feed1 = Stream(f'{ID}_single_feed1')
        self._last = Stream(f'{ID}_last')
        self._init = Stream(f'{ID}_init')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)
        
        # the reaction is simplified for simulation
        #self.mcca_rxns = ParallelRxn([
        #Rxn('Protein -> AceticAcid', 'Protein', 0.1),
        #Rxn('Protein -> PropionicAcid', 'Protein', 0.1),
        #Rxn('Protein -> ButyricAcid', 'Protein', 0.1),
        #Rxn('Protein -> CaproicAcid', 'Protein', 0.1),
        #Rxn('Glucose -> Ethanol', 'Glucose', 0.018),
        #Rxn('Glucose -> Butanol', 'Glucose', 0.0022),
        #Rxn('Glucose -> AceticAcid', 'Glucose', 0.10094),
        #Rxn('Glucose -> PropionicAcid', 'Glucose', 0.00977),
        #Rxn('Glucose -> ValericAcid', 'Glucose', 0.012),
        #Rxn('Glucose -> HeptanoicAcid', 'Glucose', 0.004),
        #Rxn('Glucose -> ButyricAcid', 'Glucose', 0.122),
        #Rxn('Glucose -> CaproicAcid', 'Glucose', 0.134),
        #Rxn('Glucose -> CaprylicAcid', 'Glucose', 0.008),
        #])        
        
        # This is a Siplified Neutralization, because we don't want to affect the products yield
        self.lime_neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('H2SO4 + Lime -> CaSO4 + 2 H2O',  'H2SO4',   1)
            ])

        self.lime_pH_control_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('H2SO4 + Lime -> CaSO4 + 2 H2O',  'H2SO4',   1)
            ])
              
 
        # self.mcca_rxns = ParallelRxn([
        #     Rxn('Glucose -> 2 Ethanol + 2 CO2', 'Glucose', 0.7),
        #     Rxn('Glucose -> 2 AceticAcid + 2 CO2 + 2 H2', 'Glucose', 0.2),
        #     Rxn('Glucose -> 2 PropionicAcid + 2 CO2 + 2 H2', 'Glucose', 0.1),
        #     Rxn('5 PropionicAcid + 6 Ethanol -> 5 ValericAcid + AceticAcid + 4 H2O + 2 H2', 'PropionicAcid', 0.7),
        #     Rxn('5 ValericAcid + 6 Ethanol -> 5 HeptanoicAcid + AceticAcid + 4 H2O + 2 H2', 'ValericAcid', 0.7),
        #     Rxn('ButyricAcid + Ethanol -> AceticAcid + Butanol', 'ButyricAcid', 0.1),
        #     Rxn('4 AceticAcid + 6 Ethanol -> 5 ButyricAcid  + 4 H2O + 2 H2', 'AceticAcid', 0.7),
        #     Rxn('5 ButyricAcid + 6 Ethanol -> 5 CaproicAcid + AceticAcid + 4 H2O + 2 H2', 'ButyricAcid', 0.7),
        #     Rxn('5 CaproicAcid + 6 Ethanol -> 5 CaprylicAcid + AceticAcid + 4 H2O + 2 H2', 'CaproicAcid', 0.2),
        # ])

    def _run(self):
        substrate, yeast_seed, FermMicrobe, lime, n2 = self.ins
        broth, gas, fermentation_waste, yeast_recycle = self.outs 
        broth.mix_from(self.ins)
        broth.T = self.T
        broth.P = self.P
        # Determine basis biomass mass (kg/hr)
        microalgae_mass = self.microalgae_mass_flow
        caproic_acid_yield_factor = getattr(self, 'caproic_acid_yield_factor', 1.0)
        
        base_yields = {
            'Ethanol': 0.01,
            'Butanol': 0.004,
            'AceticAcid': 0.1,
            'PropionicAcid': 0.009,
            'ButyricAcid': 0.18,
            'ValericAcid': 0.01,
            'CaproicAcid': 0.27,
            'HeptanoicAcid': 0.006,
            'CaprylicAcid': 0.04
        }
        
        total_base_yield = sum(base_yields.values())
        base_caproic_yield = base_yields['CaproicAcid']
        
        abs_C6 = getattr(self, 'caproic_acid_yield_absolute', None)
        if abs_C6 is not None:
            new_caproic_yield = float(abs_C6)
        else:
            new_caproic_yield = base_caproic_yield * caproic_acid_yield_factor
        delta_caproic = new_caproic_yield - base_caproic_yield
        total_other_base_yield = total_base_yield - base_caproic_yield
        
        target_caproic_mass = microalgae_mass * new_caproic_yield  
        if abs(delta_caproic) > 1e-6 and total_other_base_yield > 1e-6:
            adjustment_factor = max(0.1, (total_other_base_yield - delta_caproic) / total_other_base_yield)
            for compound, base_yield in base_yields.items():
                if compound == 'CaproicAcid':
                    broth.imass[compound] = target_caproic_mass
                else:
                    broth.imass[compound] = microalgae_mass * base_yield * adjustment_factor
        else:
            for compound, base_yield in base_yields.items():
                if compound == 'CaproicAcid':
                    broth.imass[compound] = target_caproic_mass
                else:
                    broth.imass[compound] = microalgae_mass * base_yield
        gas.copy_like(n2)
        gas.imass['H2'] = microalgae_mass * 0.1
        
        # yeast recycle for simulating membrane bioreactor
        yeast_recovery_rate = 0.95
        total_yeast_in_broth = broth.imass['Yeast']
        if total_yeast_in_broth > 0:
            recoverable_yeast = total_yeast_in_broth * yeast_recovery_rate
            yeast_recycle.empty()
            yeast_recycle.imass['Yeast'] = recoverable_yeast
            yeast_recycle.T = self.T
            yeast_recycle.P = self.P
            broth.imass['Yeast'] = total_yeast_in_broth - recoverable_yeast
        else:
            yeast_recycle.empty()
        
        fermentation_waste.empty()
        self.effluent_titer = compute_caproic_acid_titer(broth)
    def _design(self):
        Design = self.design_results
        Design['Fermenter size'] = self.outs[0].F_mass * self.tau
        Design['Broth flow rate'] = self.outs[0].F_mass
        duty = 50 * self.F_mass_in  # Heat duty
        self.add_heat_utility(duty, self.T)

    # Convenience setters for external analyses/Excel loaders
    def set_C6_titer(self, value: float):
        self.titer = float(value)

    def set_C6_yield_factor(self, value: float):
        self.caproic_acid_yield_factor = float(value)

    def set_C6_yield(self, value: float):
        # Absolute C6 yield (fraction of algae mass to CaproicAcid)
        self.caproic_acid_yield_absolute = float(value)

@cost(basis='Fermenter size', ID='Fermenter', units='kg',
      kW=100, cost=5128000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=0.7, BM=1.5)
@cost(basis='Fermenter size', ID='Fermenter agitator', units='kg',
      kW=10,
      cost=630000, 
      S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=0.7, BM=1.5)
class MCCAFermentation_no_yeast(StirredTankReactor):
    _N_ins = 4  
    _N_outs = 3

    auxiliary_unit_names = ('heat_exchanger',)
    _units= {
        **Reactor._units,
        'Fermenter size': 'kg',
        'Broth flow rate': 'kg/hr',
        }
    
    tau_batch_turnaround = 12 # in hr, seedmentaiton time
    tau_cofermentation = 15 * 24 # in hr, 15 d
    mol_NH4OH_per_acid_pH_control = 0.2988 # from IBRL batch run
    mol_lime_per_acid_pH_control = 0.2988/2 # mol NH4OH / 2
    pH = 5.0

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, T=37+273.15,
                 P=101325, V_wf=0.8, length_to_diameter=0.6,
                 kW_per_m3=0.02 * 0.7457 / 0.075, # ref from succinic projects
                 mixing_intensity=300,
                 #wall_thickness_factor=1,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 neutralization=True,
                 neutralizing_agent='Lime',
                 mode='batch', feed_freq=1,
                 pH_control=True,
                 base_neutralizes_product=True,
                 tau = tau_cofermentation,
                 microalgae_mass_flow: float | None = None,  # kg/hr basis of original algae for yield calc
                 titer: float = 1.208,  # g/L
                 # allow_dilution=False,
                 # allow_concentration=False,
                 # sugars=None
                 ):

        StirredTankReactor.__init__(self, ID, ins, outs)
        
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        #self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.neutralization = neutralization
        self.mode = mode
        self.feed_freq = feed_freq
        self.pH_control=pH_control
        self.base_neutralizes_product = base_neutralizes_product
        self.neutralization = neutralization
        self.neutralizing_agent = neutralizing_agent
        self.tau = tau 
 
        #self.allow_dilution = allow_dilution
        #self.allow_concentration = allow_concentration
        #self.sugars = sugars or tuple(i.ID for i in self.chemicals.sugars)

        # Store optional explicit biomass flow for yield calculations
        self.microalgae_mass_flow = microalgae_mass_flow
        self.titer = titer
        
        ID = self.ID
        self._mixed_feed = Stream(f'{ID}_mixed_feed')
        self._tot_feed = Stream(f'{ID}_tot_feed')
        # Before reaction, after reaction, with last feed
        self._single_feed0 = Stream(f'{ID}_single_feed0')
        self._single_feed1 = Stream(f'{ID}_single_feed1')
        self._last = Stream(f'{ID}_last')
        self._init = Stream(f'{ID}_init')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)
        
        # the reaction is simplified for simulation
        #self.mcca_rxns = ParallelRxn([
        #Rxn('Protein -> AceticAcid', 'Protein', 0.1),
        #Rxn('Protein -> PropionicAcid', 'Protein', 0.1),
        #Rxn('Protein -> ButyricAcid', 'Protein', 0.1),
        #Rxn('Protein -> CaproicAcid', 'Protein', 0.1),
        #Rxn('Glucose -> Ethanol', 'Glucose', 0.018),
        #Rxn('Glucose -> Butanol', 'Glucose', 0.0022),
        #Rxn('Glucose -> AceticAcid', 'Glucose', 0.10094),
        #Rxn('Glucose -> PropionicAcid', 'Glucose', 0.00977),
        #Rxn('Glucose -> ValericAcid', 'Glucose', 0.012),
        #Rxn('Glucose -> HeptanoicAcid', 'Glucose', 0.004),
        #Rxn('Glucose -> ButyricAcid', 'Glucose', 0.122),
        #Rxn('Glucose -> CaproicAcid', 'Glucose', 0.134),
        #Rxn('Glucose -> CaprylicAcid', 'Glucose', 0.008),
        #])        
        
        # This is a Siplified Neutralization, because we don't want to affect the products yield
        self.lime_neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('H2SO4 + Lime -> CaSO4 + 2 H2O',  'H2SO4',   1)
            ])

        self.lime_pH_control_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('H2SO4 + Lime -> CaSO4 + 2 H2O',  'H2SO4',   1)
            ])
              
 
        # self.mcca_rxns = ParallelRxn([
        #     Rxn('Glucose -> 2 Ethanol + 2 CO2', 'Glucose', 0.7),
        #     Rxn('Glucose -> 2 AceticAcid + 2 CO2 + 2 H2', 'Glucose', 0.2),
        #     Rxn('Glucose -> 2 PropionicAcid + 2 CO2 + 2 H2', 'Glucose', 0.1),
        #     Rxn('5 PropionicAcid + 6 Ethanol -> 5 ValericAcid + AceticAcid + 4 H2O + 2 H2', 'PropionicAcid', 0.7),
        #     Rxn('5 ValericAcid + 6 Ethanol -> 5 HeptanoicAcid + AceticAcid + 4 H2O + 2 H2', 'ValericAcid', 0.7),
        #     Rxn('ButyricAcid + Ethanol -> AceticAcid + Butanol', 'ButyricAcid', 0.1),
        #     Rxn('4 AceticAcid + 6 Ethanol -> 5 ButyricAcid  + 4 H2O + 2 H2', 'AceticAcid', 0.7),
        #     Rxn('5 ButyricAcid + 6 Ethanol -> 5 CaproicAcid + AceticAcid + 4 H2O + 2 H2', 'ButyricAcid', 0.7),
        #     Rxn('5 CaproicAcid + 6 Ethanol -> 5 CaprylicAcid + AceticAcid + 4 H2O + 2 H2', 'CaproicAcid', 0.2),
        # ])

    def _run(self):
        substrate, FermMicrobe, lime, n2 = self.ins
        broth, gas, fermentation_waste= self.outs 
        broth.mix_from(self.ins)
        broth.T = self.T
        broth.P = self.P
        # Determine basis biomass mass (kg/hr)
        microalgae_mass = self.microalgae_mass_flow
        caproic_acid_yield_factor = getattr(self, 'caproic_acid_yield_factor', 1.0)
        
        base_yields = {
            'Ethanol': 0.01,
            'Butanol': 0.002,
            'AceticAcid': 0.1,
            'PropionicAcid': 0.009,
            'ButyricAcid': 0.086,
            'ValericAcid': 0.01,
            'CaproicAcid': 0.08,
            'HeptanoicAcid': 0.004,
            'CaprylicAcid': 0.005
        }
        
        total_base_yield = sum(base_yields.values())
        base_caproic_yield = base_yields['CaproicAcid']
        
        abs_C6 = getattr(self, 'caproic_acid_yield_absolute', None)
        if abs_C6 is not None:
            new_caproic_yield = float(abs_C6)
        else:
            new_caproic_yield = base_caproic_yield * caproic_acid_yield_factor
        delta_caproic = new_caproic_yield - base_caproic_yield
        total_other_base_yield = total_base_yield - base_caproic_yield
        
        target_caproic_mass = microalgae_mass * new_caproic_yield  
                
        if abs(delta_caproic) > 1e-6 and total_other_base_yield > 1e-6:
            adjustment_factor = max(0.1, (total_other_base_yield - delta_caproic) / total_other_base_yield)
            for compound, base_yield in base_yields.items():
                if compound == 'CaproicAcid':
                    broth.imass[compound] = target_caproic_mass
                else:
                    broth.imass[compound] = microalgae_mass * base_yield * adjustment_factor
        else:
            for compound, base_yield in base_yields.items():
                if compound == 'CaproicAcid':
                    broth.imass[compound] = target_caproic_mass
                else:
                    broth.imass[compound] = microalgae_mass * base_yield

        gas.copy_like(n2)
        gas.imass['H2'] = microalgae_mass * 0.1

        self.effluent_titer = compute_caproic_acid_titer(broth)
    def _design(self):
        Design = self.design_results
        Design['Fermenter size'] = self.outs[0].F_mass * self.tau
        Design['Broth flow rate'] = self.outs[0].F_mass
        duty = 50 * self.F_mass_in  # Heat duty
        self.add_heat_utility(duty, self.T)

    # Convenience setters to align with yeast case
    def set_C6_yield_factor(self, value: float):
        self.caproic_acid_yield_factor = float(value)

    def set_C6_yield(self, value: float):
        self.caproic_acid_yield_absolute = float(value)

# %% 
# =========================
# Solid-Liquid Separation
# =========================
@cost(basis='Flow rate', ID='Filtrate tank', units='kg/hr',
      cost=103000, S=31815, CE=CEPCI[2010], n=0.7, BM=2.0)
@cost(basis='Flow rate', ID='Flitrate tank agitator', units='kg/hr',
      kW=5.59275, cost=26000,  S=337439, CE=CEPCI[2009], n=0.5, BM=1.5)
# ref from Biorefineries.HP
class SolidLiquidSeparation(Unit):
    _N_ins = 1
    _N_outs = 2 
    def _run(self):
        feed = self.ins[0]
        filtrate, residue = self.outs
        residue.empty()
        insolubles = ('Protein', 'Lipid', 'Ash', 'TriOlein', 'AlphaAmylase', 'GlucoAmylase', 'Yeast')
        for comp in feed.chemicals.IDs:
            if comp in insolubles:
                residue.imass[comp] = feed.imass[comp]
                filtrate.imass[comp] = 0
            else:
                filtrate.imass[comp] = feed.imass[comp]
                residue.imass[comp] = 0
        
        filtrate.T = feed.T
        filtrate.P = feed.P
        residue.T = feed.T
        residue.P = feed.P


# %% 
# =========================
# Other Units
# =========================
@cost(basis='Flow rate', ID='NeutralizationTank', units='kg/hr',
	  # Size basis on the total flow, not just ammonia, 	
      # thus assuming difference caused by MWs of NH3 and NH4OH is negligible
      cost=5000, S=157478, CE=CEPCI[2009], n=0.5, BM=1.5)
class NeutralizationTank(Unit):
    _N_ins = 2
    _N_outs = 1
    def _run(self):
        acid_stream, ammonium_hydroxide = self.ins
        self.outs[0].empty()
        self.outs[0].mix_from([acid_stream, ammonium_hydroxide])
        neutralization_rxn = rxn.Rxn('2 NH4OH + H2SO4 -> AmmoniumSulfate + 2 H2O', 'NH4OH', 1)
        neutralization_rxn(self.outs[0])




@cost('Flow rate', 'Anaerobic Digestion Reactor', units= 'kg/hr', kW=100, S=339151,
      cost=8000000, CE=CEPCI[2010], n=0.7, BM=1.8)
class AnaerobicDigestion(bst.Unit):

    _N_ins = 1
    _N_outs = 3
   
    degrade_IDs = [
        'AceticAcid', 'PropionicAcid', 'ButyricAcid', 'ValericAcid',
        'CaproicAcid', 'HeptanoicAcid', 'CaprylicAcid',
        'Ethanol', 'Butanol', 'Glucose', 'Protein',
        'AlphaAmylase', 'GlucoAmylase', 'TriOlein', 'WWTsludge'
    ]
    inorganic_IDs = ['Ash', 'H2SO4', 'Na2SO4', 'AmmoniumSulfate', 'H2O']

    def __init__(self, ID='', ins=None, outs=(), *, tau=15*24, microalgae_mass: float|None=None, **kwargs):
        super().__init__(ID, ins, outs, **kwargs)
        self.tau = tau 
        self.microalgae_mass = microalgae_mass 
        self.T = 37 + 273.15 
        
    V_wf = 0.8
    
    N_reactors = 0
    
    N_transfer_pumps = 1
    
    N_recirculation_pumps = 0

    def _run(self):
        feed = self.ins[0]
        biogas, waste, sludge = self.outs
        if self.microalgae_mass is None:
            raise ValueError('microalgae_mass must be supplied when instantiating AnaerobicDigestion.')
        microalgae_mass = self.microalgae_mass
        H2_mol = microalgae_mass * 0.015
        CH4_mol = microalgae_mass * 0.015
        biogas.empty()
        biogas.imol['H2'] = H2_mol
        biogas.imol['CH4'] = CH4_mol
        biogas.phase = 'g'
        biogas.T = self.T
        waste.empty()
        sludge.empty()
        for ID in self.degrade_IDs:
            amount = feed.imass[ID]
            waste.imass[ID] = amount * 0.05  # 5% into wastewater
            sludge.imass[ID] = amount * 0.25 # 25% into sludge
        for ID in self.inorganic_IDs:
            amount = feed.imass[ID]
            waste.imass[ID] = amount * 0.95
            sludge.imass[ID] = amount * 0.05
        sludge.imass['H2O'] = feed.imass['H2O']*0.3
        waste.imass['H2O'] = feed.imass['H2O']*0.7
        sludge.imass['Microalgae'] = 0
        waste.phase = 'l'
        waste.T = sludge.T = self.T
        
    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = sum([s.F_mass for s in self.outs])  # kg/hr
        duty = 800 * self.microalgae_mass  # Heat duty
        self.add_heat_utility(duty, self.T)







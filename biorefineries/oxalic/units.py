#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>, Wenjun Guo <wenjung2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
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
from warnings import warn
from flexsolve import aitken_secant
from biosteam import Unit
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank, LiquidsSplitSettler, StirredTankReactor, BatchBioreactor, Compressor
from biosteam.units.decorators import cost
from thermosteam import Stream, MultiStream
from biorefineries.oxalic.process_settings import price
from biorefineries.oxalic.utils import CEPCI, baseline_feedflow, compute_extra_chemical, adjust_recycle
from biorefineries.oxalic.chemicals_data import oxalic_chemicals
tmo.settings.set_thermo(oxalic_chemicals)

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

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=6210, S=1981, CE=CEPCI[2010],  n=0.7, BM=3)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      cost=8000, S=3720, CE=CEPCI[2009],  n=0.8, BM=2.3)
class SulfuricAcidAdditionTank(Unit):
    _N_ins = 1
    _N_outs = 1
    
    # Bseline is (18+4.1) mg/g dry biomass as in Humbird et al., 93% purity
    acid_loading = 22.1

    def _init(self, feedstock_dry_mass):
        
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
        # water adjustment currently implemented in H_M201.specification
        
        mixture.mix_from([water, acid])



# Adjust pretreatment water loading, 30% from Table 5 on Page 21 of Humbird et al.
class PretreatmentMixer(Mixer):
    _N_ins = 4
    _N_outs = 1
    
    solid_loading = 0.3
        
    def _run(self):
        feedstock, acid, water, recycled_water = self.ins
        mixture_out = self.outs[0]
        
        mixture = feedstock.copy()
        mixture.mix_from([feedstock, acid])

        # water adjustment implemented in H_M202.specification
        
        mixture_out.mix_from([mixture, water, recycled_water])
 
# Pretreatment reactor
@cost(basis='Dry flow rate', ID='Pretreatment reactor', units='kg/hr',
      kW=5120, cost=19812400, S=83333, CE=CEPCI[2009], n=0.6, BM=1.5)
class PretreatmentReactorSystem(Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = Flash._graphics
    
    def _init(self, T=130+273.15):
        
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
            Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.),
            Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9),
            Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.024),
            Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.05),
            Rxn('Acetate -> AceticAcid',                     'Acetate',  1.),
            Rxn('Lignin -> 10.125 SolubleLignin',            'Lignin',   0.05),
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
            Rxn('Furfural -> Tar',                           'Furfural', 1.),
            Rxn('HMF -> Tar',                                'HMF',      1.)
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
    tau_saccharification = 24. # Humbird
    

    def _init(self, T=50+273.15):
        
        self.T = T
        
        # Based on Table 9 on Page 28 of Humbird et al.
        
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
        ss  = self.outs[0]
        # ss = self.saccharified_stream
        
        ss.copy_like(feed)
        
        ss.T = self.T
        
        self.saccharification_rxns(ss.mol)
        


        
    def _design(self):
        Design = self.design_results
        total_mass_flow = self.ins[0].F_mass 
        Design['Saccharification tank size'] = total_mass_flow * self.tau_saccharification
        Design['Slurry flow rate'] = total_mass_flow
        

        
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
    _N_outs = 2
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}
    
    # Cycle time for each batch (hr), including 12 hr turnaround time 
    tau_batch = 36
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    
    
    
    
    def _init(self, T=30.+273.15, ferm_ratio=0.9):
        
        self.T = T
        self.ferm_ratio = ferm_ratio

        self.regular_succinic_acid_conversion = 0. * ferm_ratio
        
        self.regular_xylitol_conversion = 0.04 * ferm_ratio
        
        self.regular_biomass_conversion = 0.08 * ferm_ratio
        
        # FermMicrobe reaction from Table 14 on Page 31 of Humbird et al.
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 OxalicAcid + 2AceticAcid + O2',        'Glucose',   .49),
        # Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
        Rxn('Glucose -> 2 SuccinicAcid + O2',        'Glucose',   0.),
        Rxn('3 Xylose -> 5 OxalicAcid + 5 AceticAcid + 2.5 O2',       'Xylose',    0.49),
        Rxn('Xylose + -> 1.667 SuccinicAcid + 0.4165 O2',       'Xylose',    1e-5),
        Rxn('Glucose + 0.545 H2O -> 1.091 Xylitol + 0.545 CO2',        'Glucose',   0.040),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',       'Xylose',    0.040),
        Rxn('Glucose -> 6 FermMicrobe + 2.4 H2O',       'Glucose',   0.05-1e-9), 
        Rxn('Xylose -> 5 FermMicrobe + 2 H2O',        'Xylose',    0.05-1e-9), 
        ])
        
        self.CO2_generation_rxns = ParallelRxn([
        Rxn('Glucose -> 6 CO2 + 6H2O',       'Glucose',   1.-1e-9),
        Rxn('Xylose -> 5 CO2 + 5H2O',        'Xylose',    1.-1e-9),
        ])

        
        if 'Sucrose' in self.chemicals:
            self.sucrose_hydrolysis_rxn = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.-1e-9)
            
        
        self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
        self.xylose_to_HP_rxn = self.cofermentation_rxns[2]
        
        self.glucose_to_succinic_acid_rxn = self.cofermentation_rxns[1]
        self.xylose_to_succinic_acid_rxn = self.cofermentation_rxns[3]
        
        self.glucose_to_xylitol_rxn = self.cofermentation_rxns[4]
        self.xylose_to_xylitol_rxn = self.cofermentation_rxns[5]
        
        self.glucose_to_biomass_rxn = self.cofermentation_rxns[6]
        self.xylose_to_biomass_rxn = self.cofermentation_rxns[7]
        self.biomass_generation_rxns = [self.glucose_to_biomass_rxn, self.xylose_to_biomass_rxn]
        
        
    def _run(self):
        feed = self.ins[0]
        effluent, CO2 = self.outs
        effluent.copy_like(feed)
        CO2.phase = 'g'
        effluent.phase='l'
        
        if 'Sucrose' in effluent.chemicals:
            # self.sucrose_hydrolysis_rxn.force_reaction(effluent)
            self.sucrose_hydrolysis_rxn(effluent)
            if effluent.imol['Water'] < 0.: effluent.imol['Water'] = 0.
        
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        
        # Assume all CSL is used up
        effluent.imass['CSL'] = 0 
        
        effluent.T = CO2.T = self.T
        CO2.imass['CO2'] = effluent.imass['CO2']
        effluent.imass['CO2'] = 0

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
              'Reactor volume': 'm3',
              'Single reactor volume': 'm3'}
    
    # For a single reactor, based on diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147 
    
    def _init(self, 
                  P=101325, tau=0.5, V_wf=0.8,
                  length_to_diameter=2, kW_per_m3=0.985,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical'):
        
        
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
        P, D, L = float(P), float(D), float(L)
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
        baseline_purchase_costs = self.baseline_purchase_costs
        
        if Design['Total volume'] == 0:
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] = 0
        
        else:
            baseline_purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in baseline_purchase_costs.items():
                baseline_purchase_costs[i] *= Design['Number of reactors']
            
            self.power_utility(self.kW_per_m3 * Design['Total volume'])
            
            

#%%        
compute_HP_titer = lambda effluent: (effluent.imass['OxalicAcid'] +
            effluent.imol['CalciumOxalate']*oxalic_chemicals.HP.MW)/effluent.F_vol

compute_HP_mass = lambda effluent: effluent.imass['OxalicAcid'] +\
            effluent.imol['CalciumOxalate']*oxalic_chemicals.HP.MW
_316_over_304 = 1.2


#%% Old fermentation reactor code
# @cost(basis='Fermenter size', ID='Fermenter', units='kg',
#       cost=10128000, S=(42607+443391+948+116)*(60+36),
#       CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Fermenter size', ID='Agitator', units='kg',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in ref [3])
#       # and total residence time (batch hydrolysis and fermentation)
#       kW=268.452, cost=630000, S=(42607+443391+948+116)*(60+36),
#       CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Recirculation flow rate', ID='Recirculation pump', units='kg/hr',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in ref [3])
#       kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)  
# class CoFermentation(Reactor):
#     _N_ins = 9
#     _N_outs = 2
#     # _N_heat_utilities = 1
#     _units= {**Reactor._units,
#             'Fermenter size': 'kg',
#             'Recirculation flow rate': 'kg/hr',
#             'Duty': 'kJ/hr'}
#     _F_BM_default = {**Reactor._F_BM_default,
#             # 'Heat exchangers': 3.17,
#             }

#     auxiliary_unit_names = ('heat_exchanger',)
    
#     ferm_ratio = 1.

#     effluent_titer = 0
    
#     productivity = 0.76 # in g/L/hr
    
#     CSL_loading = 10 # g/L (kg/m3)
    
#     magnesium_chloride_loading = 1.
    
#     diammonium_sulfate_loading = 1.
    
#     fraction_of_biomass_C_from_CO2 = 0.
    
#     CO2_safety_factor = 3.
    
#     regular_succinic_acid_conversion = 0.04
    
#     regular_xylitol_conversion = 0.04
    
#     regular_biomass_conversion = 0.08
    
#     def _init(self, thermo=None, *, T=30.+273.15,
#                   P=101325, V_wf=0.8, length_to_diameter=2,
#                   kW_per_m3=0.0985, # Perry's handbook
#                    wall_thickness_factor=1,
#                    # tau=120., #reset by spec.set_productivity
#                   vessel_material='Stainless steel 316',
#                   vessel_type='Vertical',
#                   neutralization=True,
#                   mode='Continuous', # Batch or Continuous
#                   allow_dilution=False,
#                   allow_concentration=False):
        
        
#         Reactor._init(self, P=P, 
#                       # tau=tau,
#                       V_wf=V_wf,
#         length_to_diameter=length_to_diameter, kW_per_m3=kW_per_m3,
#         wall_thickness_factor=wall_thickness_factor,
#         vessel_material=vessel_material,
#         vessel_type=vessel_type)
        
#         self.T = T
#         self.P = P
#         self.V_wf = V_wf
#         self.length_to_diameter = length_to_diameter
#         self.kW_per_m3 = kW_per_m3
#         # self.wall_thickness_factor = wall_thickness_factor
#         self.vessel_material = vessel_material
#         self.vessel_type = vessel_type
#         self.neutralization = neutralization
#         self.mode = mode        
#         self.allow_dilution = allow_dilution
#         self.allow_concentration = allow_concentration
#         self.mixed_feed = tmo.Stream('mixed_feed')

#         # self.heat_exchanger = HXutility(None, None, None, T=T)
        
#         self.cofermentation_rxns = ParallelRxn([
#         #      Reaction definition            Reactant    Conversion
#         Rxn('Glucose -> 2 HP',        'Glucose',   .49),
#         # Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
#         Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.040),
#         Rxn('3 Xylose -> 5 HP',       'Xylose',    0.49),
#         Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.040),
#         Rxn('Glucose + 2 H2O -> 2 Glycerol',        'Glucose',     0.040),
#         Rxn('3Xylose + 5 H2O -> 5 Glycerol',       'Xylose',    0.040),
#         Rxn('Glucose -> 6 FermMicrobe + 2.4 H2O',       'Glucose',   0.05-1e-9), # for consistency with Dunn et al. 2015, which has 5% conversion of algal glycerol to cell mass
#         Rxn('Xylose -> 5 FermMicrobe + 2 H2O',        'Xylose',    0.05-1e-9), # for consistency with Dunn et al. 2015, which has 5% conversion of algal glycerol to cell mass
#         ])
        
        
#         self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
#         self.xylose_to_HP_rxn = self.cofermentation_rxns[2]
        
#         self.glucose_to_succinic_acid_rxn = self.cofermentation_rxns[1]
#         self.xylose_to_succinic_acid_rxn = self.cofermentation_rxns[3]
        
#         self.glucose_to_xylitol_rxn = self.cofermentation_rxns[4]
#         self.xylose_to_xylitol_rxn = self.cofermentation_rxns[5]
        
#         self.CO2_generation_rxns = ParallelRxn([
#         Rxn('Glucose -> 6 CO2 + 6H2O',       'Glucose',   1.-1e-9),
#         Rxn('Xylose -> 5 CO2 + 5H2O',        'Xylose',    1.-1e-9),
#         ])
        
#         self.glucose_to_biomass_rxn = self.cofermentation_rxns[6]
#         self.xylose_to_biomass_rxn = self.cofermentation_rxns[7]
#         self.biomass_generation_rxns = [self.glucose_to_biomass_rxn, self.xylose_to_biomass_rxn]
        
#         self.glucose_to_microbe_rxn = self.biomass_generation_rxns[0]
#         self.xylose_to_microbe_rxn = self.biomass_generation_rxns[1]

#         if 'Sucrose' in self.chemicals:
#             self.sucrose_hydrolysis_rxn = Rxn('Sucrose + Water -> 2Glucose', 'Sucrose', 1.-1e-9)
        
#         self._X = self.cofermentation_rxns.X.copy()
        
#         # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
#         self.neutralization_rxns = ParallelRxn([
#         #   Reaction definition                                               Reactant  Conversion
#         Rxn('2 HP + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'HP',   1.-1e-9),
#         Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1-1e-9),
#             ])
#         self.tau = self.tau_cofermentation = 74 # this value is altered by spec.load_productivity
        
#         self.get_succinic_acid_conc = lambda: self.outs[0].imass['AceticAcid']/self.outs[0].F_vol
#     def _run(self):
        
#         sugars, feed, CSL, magnesium_chloride, diammonium_sulfate, lime, CO2_fresh, CO2_recycled, air_fresh = self.ins
        
#         effluent, vapor = self.outs
#         effluent.mix_from([feed, sugars, CSL])
#         effluent.phase='l'
#         if 'Sucrose' in effluent.chemicals:
#             # self.sucrose_hydrolysis_rxn.force_reaction(effluent)
#             self.sucrose_hydrolysis_rxn(effluent)
#             if effluent.imol['Water'] < 0.: effluent.imol['Water'] = 0.
        
#         effluent.T = vapor.T = self.T
#         CSL.imass['CSL'] = (sugars.F_vol + feed.F_vol) * self.CSL_loading 
        
#         self.cofermentation_rxns(effluent.mol)
#         self.CO2_generation_rxns(effluent.mol)
       
#         self.fresh_CO2_required = fresh_CO2_required = (effluent.imol['FermMicrobe']*\
#                               self.fraction_of_biomass_C_from_CO2*
#                               self.CO2_safety_factor) - CO2_recycled.imol['CO2']
            
#         CO2_fresh.imol['CO2'] = max(0, fresh_CO2_required)
#         effluent.mix_from([effluent, CO2_fresh, CO2_recycled])
        
        
#         vapor.imol['CO2'] = effluent.imol['CO2']
#         vapor.phase = 'g'
        
#         effluent.imol['CO2'] = 0
#         effluent.imass['CSL'] = 0
        
#         mixed_feed = self.mixed_feed
#         mixed_feed.mix_from([feed, sugars, CSL])
        
#         # Need lime to neutralize produced organic acids
#         if self.neutralization:
#             self.vessel_material= 'Stainless steel 316'
#             # Set feed lime mol to match rate of acids production, add 10% extra
#             lime.imol['Lime'] = (effluent.imol['HP']/2/self.neutralization_rxns.X[0]) \
#                                 + effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1] * 1.1
                               
#             effluent.imol['Lime'] = lime.imol['Lime']
#             self.neutralization_rxns(effluent)
            
#         else:
#             self.vessel_material= 'Stainless steel 316'
#             lime.empty()
#         self.effluent_titer = compute_HP_titer(effluent)
        
#         self.abs_mass_bal_diff_CO2_added_to_vent = abs_mass_bal_diff_CO2_added_to_vent = \
#             max(0, self.mol_atom_in('C')-self.mol_atom_out('C'))
        
#         vapor.imol['CO2'] += abs_mass_bal_diff_CO2_added_to_vent
        
#     def _design(self):
#         super()._setup()
#         mode = self.mode
#         Design = self.design_results
#         Design.clear()
       
#         _mixture = self._mixture = tmo.Stream(None)
#         _mixture.mix_from(self.outs[0:2])
#         duty = Design['Duty'] = self.Hnet
#         if mode == 'Batch':
#             raise NotImplementedError('Batch mode is not currently implemented for this biorefinery')

#         elif mode == 'Continuous':
#             Reactor._V_max = 3785.41 
#             Reactor._design(self)
#         else:
#             raise DesignError(f'Fermentation mode must be either Batch or Continuous, not {mode}')

#     def _cost(self):
#         Design = self.design_results
#         # baseline_purchase_costs = self.baseline_purchase_costs
#         baseline_purchase_costs = self.baseline_purchase_costs

#         # hx = self.heat_exchanger

#         if self.mode == 'Batch':
#             raise NotImplementedError('Batch mode is not currently implemented for this biorefinery')
           

#         elif self.mode == 'Continuous':
            
#             Reactor._cost(self)
            
#             # N = Design['Number of reactors']
#             # single_rx_effluent = self._mixture.copy()
#             # hx.simulate_as_auxiliary_exchanger(duty=Design['Duty']/N, 
#             #                                 stream=single_rx_effluent)
#             # self.auxiliary_unit_names = names = tuple([f'heat_exchanger_{i}' for i in range(N)])
#             # for i in names: setattr(self, i, hx)
#             # hu_total = self.heat_utilities[0]
#             # hu_single_rx = hx.heat_utilities[0]
#             # hu_total.copy_like(hu_single_rx)
#             # self.heat_utilities = tuple([self.heat_utilities[0]] * N)
            
#             # hu_total.heat_exchanger = None
            


#     def mol_atom_in(self, atom):
#         return sum([stream.get_atomic_flow(atom) for stream in self.ins])
    
#     def mol_atom_out(self, atom):
#         return sum([stream.get_atomic_flow(atom) for stream in self.outs])

#%%
class BatchCoFermentation(BatchBioreactor):
    # Co-Fermentation time (hr)
    
    _N_ins = 9
    _N_outs = 2
    
    tau_cofermentation = 120 # initial value; updated by spec.load_productivity
    
    productivity = 0.76 # in g/L/hr
    
    CSL_loading = 23 # g/L (kg/m3)
    
    magnesium_chloride_loading = 0. # g/L (kg/m3)
    
    diammonium_sulfate_loading = 5. # g/L (kg/m3)
    
    CSL_and_micronutrient_consumption_fraction = 0.999 # assumed
    
    fraction_of_biomass_C_from_CO2 = 0.5 # % as decimal
    
    CO2_safety_factor = 2.
    
    regular_succinic_acid_conversion = 0. # %theoretical
    
    # regular_xylitol_conversion = 0.085 # %theoretical
    regular_xylitol_conversion = 0. # %theoretical
    
    regular_biomass_conversion = 0.245 # %theoretical
    
    air_m3_per_h_per_m3_reactor = 3.5*60/3 # 3.5 slpm for a 3L bioreactor; Markham et al. 2018 # used when aeration_rate_basis=='fixed rate basis'
    
    DO_saturation_concentration_kg_per_m3 = 7.55e-3 # 7.55 mg/L at 30 degrees C # used when aeration_rate_basis=='DO saturation basis'
    
    DO_saturation_target_level = 0.1 # 10% saturation from Markham et al. 2018 # used when aeration_rate_basis=='DO saturation basis'
    
    air_flow_rate_safety_factor_for_DO_saturation_basis = 2.
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.5 # !!!
    
    # autoselect_N  = True
    
    lime_mol_per_mol_loading = 0.91 # 1.82 mol KOH is for 1 mol oxalic aicd
    
    def _init(self, 
                 #  
                 T=30+273.15,
                 P=101325., 
                 tau=120, # initial value; updated by spec.load_productivity
                 V=3785.,
                 neutralization=False,
                 aeration_rate_basis='DO saturation basis', # 'fixed rate basis' or 'DO saturation basis'
                 ferm_ratio = 1.
                 ):
        BatchBioreactor._init(self,  T=T, P=P, tau=tau, V=V)
        
        self.ferm_ratio = ferm_ratio
        self.aeration_rate_basis = aeration_rate_basis
        self.neutralization = neutralization

        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 OxalicAcid + 2AceticAcid + O2',        'Glucose',   .49),
        # Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
        Rxn('Glucose -> 2 SuccinicAcid + O2',        'Glucose',   0.),
        Rxn('3 Xylose -> 5 OxalicAcid + 5 AceticAcid + 2.5 O2',       'Xylose',    0.49),
        Rxn('Xylose -> 1.667 SuccinicAcid + 0.4165 O2',       'Xylose',    1e-5),
        Rxn('Glucose + 0.545 H2O -> 1.091 Xylitol + 0.545 CO2',        'Glucose',   0.040),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',       'Xylose',    0.040),
        Rxn('Glucose -> 6 FermMicrobe + 2.4 H2O',       'Glucose',   0.05-1e-9), 
        Rxn('Xylose -> 5 FermMicrobe + 2 H2O',        'Xylose',    0.05-1e-9), 
        ])
        
        self.glucose_to_HP_rxn = self.cofermentation_rxns[0]
        self.xylose_to_HP_rxn = self.cofermentation_rxns[2]
        
        self.glucose_to_succinic_acid_rxn = self.cofermentation_rxns[1]
        self.xylose_to_succinic_acid_rxn = self.cofermentation_rxns[3]
        
        self.glucose_to_xylitol_rxn = self.cofermentation_rxns[4]
        self.xylose_to_xylitol_rxn = self.cofermentation_rxns[5]
        
        self.CO2_generation_rxns = ParallelRxn([
        Rxn('Glucose -> 6 CO2 + 6H2O',       'Glucose',   1.-1e-9),
        Rxn('Xylose -> 5 CO2 + 5H2O',        'Xylose',    1.-1e-9),
        ])
        
        self.sucrose_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
        ])
        
        self.glucose_to_biomass_rxn = self.cofermentation_rxns[6]
        self.xylose_to_biomass_rxn = self.cofermentation_rxns[7]
        self.biomass_generation_rxns = [self.glucose_to_biomass_rxn, self.xylose_to_biomass_rxn]
        
        self.glucose_to_microbe_rxn = self.biomass_generation_rxns[0]
        self.xylose_to_microbe_rxn = self.biomass_generation_rxns[1]

        self._X = self.cofermentation_rxns.X.copy()
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('OxalicAcid + CalciumDihydroxide -> CalciumOxalate + 2 H2O',  'OxalicAcid',   1.-1e-9), # !!! update if neutralization scenario is needed
        Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1-1e-9),
            ])

    def _run(self):
        feed, seed, CSL, magnesium_chloride, diammonium_sulfate, lime, CO2_fresh, CO2_recycled, air = self.ins
        vapor, effluent = self.outs
        
        for i in [CSL, magnesium_chloride, diammonium_sulfate, lime, CO2_fresh, air]: i.empty()
        
        vapor.empty()
        effluent.empty()
        
        effluent.mix_from([feed, seed])
        feed_seed_vol = effluent.F_vol
        
        # current_acetate_loading = effluent.imass[self.acetate_ID] / effluent.F_vol
        # required_acetate_spiking = max(0, self.acetate_target_loading - current_acetate_loading)
        # Acetate_spiking.imass[self.acetate_ID] = required_acetate_spiking * effluent.F_vol
        
        # air.P = self.air_pressure
        air.imol['N2'] = 0.79
        air.imol['O2'] = 0.21
        
        if self.aeration_rate_basis == 'DO saturation basis':
            self.air_exit_F_mol_needed = (1./0.21) * (1/32.) * self.DO_saturation_concentration_kg_per_m3 * self.DO_saturation_target_level\
                *(seed.F_vol+feed.F_vol)
            
            air.F_mol = 1e8 # initial value; updated after reactions
        
        elif self.aeration_rate_basis == 'fixed rate basis':
            air.F_vol = self.air_m3_per_h_per_m3_reactor * (seed.F_vol+feed.F_vol) * self.tau
        
        else: raise RuntimeError(f"Unsupported aeration_rate_basis ({self.aeration_rate_basis}); must be 'fixed rate basis' or 'DO saturation basis'.")
        
        CSL.imass['CSL'] = feed_seed_vol * self.CSL_loading 
        magnesium_chloride.imass['MagnesiumChloride'] = feed_seed_vol * self.magnesium_chloride_loading 
        diammonium_sulfate.imass['AmmoniumSulfate'] = feed_seed_vol * self.diammonium_sulfate_loading 
        
        # lime.imol['Lime'] = feed_seed_vol * self.lime_mol_per_L_loading 
        
        effluent.mix_from([effluent, air, CSL, magnesium_chloride, diammonium_sulfate])
        
        effluent.phase = 'l'
        effluent.T = vapor.T = self.T
        
        self.sucrose_hydrolysis_rxns(effluent.mol)
        
        self.fresh_CO2_required = fresh_CO2_required = (
                                      (effluent.imol['FermMicrobe']*
                                      self.fraction_of_biomass_C_from_CO2*
                                      self.CO2_safety_factor) +
                                      (2 * self.glucose_to_succinic_acid_rxn.X * effluent.imol['Glucose'] + 
                                       1.667 * self.xylose_to_succinic_acid_rxn.X * effluent.imol['Xylose'])
                                      - CO2_recycled.imol['CO2']
                                      )
        
        CO2_fresh.imol['CO2'] = max(0, fresh_CO2_required)
        effluent.mix_from([effluent, CO2_fresh, CO2_recycled])
        
        if self.aeration_rate_basis == 'DO saturation basis':
            O2_mol_remaining, N2_mol_remaining = effluent.imol['O2', 'N2']
            O2_mol_excess = O2_mol_remaining - 0.21*self.air_exit_F_mol_needed
            # N2_mol_excess = N2_mol_remaining - 0.79*air_exit_F_mol_needed
            N2_mol_excess = O2_mol_excess * 0.79/0.21
            air_mol_excess = O2_mol_excess + N2_mol_excess
            
            self.air_mol_excess = air_mol_excess
        
            air.F_mol -= air_mol_excess
            effluent.imol['O2'] -= O2_mol_excess
            effluent.imol['N2'] -= N2_mol_excess
        
        aeration_safety_factor = self.air_flow_rate_safety_factor_for_DO_saturation_basis
        effluent.imol['O2'] += (aeration_safety_factor-1) * air.imol['O2']
        effluent.imol['N2'] += (aeration_safety_factor-1) * air.imol['N2']
        air.F_mol += (aeration_safety_factor-1) * air.F_mol
        
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        vapor.imol['CO2', 'O2', 'N2'] = effluent.imol['CO2', 'O2', 'N2']
        vapor.phase = 'g'
        effluent.imol['CO2', 'O2', 'N2'] = 0, 0, 0
        effluent.imass['CSL', 'MagnesiumChloride', 'ZincSulfate'] = \
            (1.-self.CSL_and_micronutrient_consumption_fraction) *\
            effluent.imass['CSL', 'MagnesiumChloride', 'ZincSulfate']
        
        vapor.imol['CO2'] += CSL.get_atomic_flow('C')
        
        # self.effluent_titer = compute_TAL_titer(effluent)
        # self.effluent = effluent
        
        # Need lime to neutralize produced organic acids
        
        if self.neutralization:
            self.vessel_material= 'Stainless steel 316'
            # Set feed lime mol to match rate of acids production, observed ratio was 1.82 moles of KOH per mole of oxalic acid produced
            lime.imol['Lime'] = effluent.imol['OxalicAcid'] * self.lime_mol_per_mol_loading
            # self.neutralization_rxns[0](effluent)
            effluent.imol['OxalicAcid'] -= lime.imol['Lime']
            effluent.imol['CalciumOxalate'] += lime.imol['Lime']
            effluent.imol['H2O'] += 2 * lime.imol['Lime']
            
        else:
            self.vessel_material= 'Stainless steel 316'
            # lime.empty()
        
    @property
    def effluent_titer(self):
        return compute_HP_titer(self.effluent)


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


class AcidulationReactor(Reactor):
    _N_ins = 2
    _N_outs = 1
    
    acid_safety_factor = 1.05 
    def _init(self, *args, **kwargs):
        Reactor._init(self, *args, **kwargs)
        self.acidulation_rxns = ParallelRxn([
            #   Reaction definition                                        Reactant        Conversion
            Rxn('CalciumOxalate + H2SO4 -> OxalicAcid + CaSO4',         'CalciumOxalate',       1.),
            Rxn('CalciumAcetate + H2SO4 -> 2 AceticAcid + CaSO4',         'CalciumAcetate',       1.),
            Rxn('CalciumDihydroxide + H2SO4 -> CaSO4 + 2 H2O',            'CalciumDihydroxide',   1.)
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
        
        # Set feed acid mol to match acidulation needs with some extra
        acid.imol['H2SO4'] = needed_acid * self.acid_safety_factor
        acid.imass['H2O'] = acid.imass['H2SO4'] / 0.93 * 0.07 # 93% purity
        effluent.mix_from([feed, acid])
        # rxns.adiabatic_reaction(effluent)
        rxns[0](effluent)
        
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


class DehydrationReactor(Reactor):
    """
    A dehydration reactor.
    """
    _N_ins = 4
    _N_outs = 2
    
    # _N_heat_utilities = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'TiO2 catalyst': 1,
            # 'Heat exchangers': 3.17,
            }

    mcat_frac = 12/1.5 # kg per kg/h # 1/WHSV # Calculated from Dishisha et al. 2015

    # _equipment_lifetime = {'TiO2 catalyst': 1,}
    def _init(self, thermo=None, *, T=230.+273.15,
                  P=101325., V_wf=0.8, length_to_diameter=2, tau = 1,
                  kW_per_m3=0.0985, # Perry's handbook
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 304',
                  vessel_type='Vertical', X = 0.80):  # 80% based on Dunn et al. 2015
        
        
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        # self.heat_exchanger = HXutility(None, None, None, T=T)
        self.tau = tau
        self.X = X
        
        self.dehydration_reactions = dehydration_reactions = ParallelRxn([
        #   Reaction definition                                       Reactant   Conversion
        Rxn('HP -> AA + H2O',         'HP',   self.X)
            ])     
        
    def _run(self):
        
        HP_to_AA_rxn = self.dehydration_reactions[0]
        
        feed, fresh_catalyst, recycled_HP = self.ins
        effluent, spent_catalyst = self.outs
  
        effluent.mix_from([feed, recycled_HP])
        effluent.T = self.T
        effluent.P = feed.P
        self.dehydration_reactions(effluent.mol)
        effluent.phase = 'l'
        
        fresh_catalyst.imass['TiO2'] = spent_catalyst.imass['TiO2'] =\
            self.mcat_frac * self.ins[0].F_mass/(350.*24.) 
            # assuming a TiO2 lifetime of 1 year
            
    def _cost(self):
        super()._cost()
        # hx = self.heat_exchanger
        # N = self.design_results['Number of reactors']
        # single_rx_effluent = self.outs[0].copy()
        # single_rx_effluent.mol[:] /= N
        
        # hx.simulate_as_auxiliary_exchanger(duty=(self.Hnet)/N, 
        #                                     stream=single_rx_effluent)
        # hu_total = self.heat_utilities[0]
        # hu_single_rx = hx.heat_utilities[0]
        # hu_total.copy_like(hu_single_rx)
        # hu_total.scale(N)
        
        # self.purchase_costs['Heat exchangers'] = hx.purchase_cost * N
        # self.baseline_purchase_costs['Heat exchangers'] = hx.baseline_purchase_cost * N
        # self.installed_costs['Heat exchangers'] = hx.installed_cost * N
        self.baseline_purchase_costs['TiO2 catalyst'] = self.mcat_frac * self.ins[0].F_mass * price['TiO2']
            

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
    
    def _init(self, reactants, split=(), T=35+273.15):	
        	
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
    
    def _init(self, reactants, ratio=0):
        
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
    
# Modified from bst.units.StorageTank, which won't simulate for 0 flow 	
class HPStorageTank(StorageTank):
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.baseline_purchase_costs['Tanks'] = 0
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
            Cost = self.baseline_purchase_costs
            Cost['Pump'] = 0
            Cost['Motor'] = 0
        else: Pump._cost(self)
        

#%% IsothermalCompressor 

# Same as bst.IsothermalCompressor, except that it allows switching with compressible gases 
# in specification while correctly calculating power requirements


class IsothermalCompressor(Compressor, new_graphics=False):
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        self.ideal_power, self.ideal_duty = None, None
        out.copy_like(feed)
        out.P = self.P
        out.T = feed.T
        if self.vle is True: out.vle(T=out.T, P=out.P)
        self.ideal_power, self.ideal_duty = self._calculate_ideal_power_and_duty()

    def _design(self):
        super()._design()
        feed = self.ins[0]
        outlet = self.outs[0]
        (ideal_power, ideal_duty) = self._calculate_ideal_power_and_duty() if not (self.ideal_power and self.ideal_duty) else (self.ideal_power, self.ideal_duty)
        Q = ideal_duty / self.eta
        self.add_heat_utility(unit_duty=Q, T_in=feed.T, T_out=outlet.T)
        self.design_results['Ideal power'] = ideal_power # kW
        self.design_results['Ideal duty'] = ideal_duty # kJ / hr
        self._set_power(ideal_power / self.eta)
        
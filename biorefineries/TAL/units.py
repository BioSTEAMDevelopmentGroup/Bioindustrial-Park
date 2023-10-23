#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

All units are explicitly defined here for transparency and easy reference

@author: sarangbhagwat

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

All units are explicitly defined here for transparency and easy reference

"""


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp, pi, log
from flexsolve import aitken_secant
from biosteam import Unit, BatchCrystallizer
from biosteam.units import Flash, HXutility, Mixer, MixTank, Pump, \
    SolidsSeparator, StorageTank, LiquidsSplitSettler, \
    BatchBioreactor
from biosteam.units.decorators import cost
from biosteam.units.design_tools import CEPCI_by_year as CEPCI
from thermosteam import Stream, MultiStream
from biorefineries.TAL.process_settings import price
# from biorefineries.TAL.utils import CEPCI, baseline_feedflow, compute_extra_chemical, adjust_recycle
from biorefineries.TAL.utils import baseline_feedflow, compute_extra_chemical, adjust_recycle

_kg_per_ton = 907.18474
_Gcal_2_kJ = 4.184 * 1e6 # (also MMkcal/hr)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

compute_TAL_titer = lambda effluent: effluent.imass['TAL'] / effluent.F_vol

compute_TAL_mass = lambda effluent: effluent.imass['TAL']

#%% Reactor
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
                  length_to_diameter=2, 
                  kW_per_m3=0.985,
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
    
    @property
    def F_vol_in(self): # exclude gases
        return sum([i.F_vol for i in self.ins if i.phase=='l' and i.F_mol and not (i.imol['CO2']/i.F_mol==1. 
                                                                       or i.imol['O2']/i.F_mol>0.1
                                                                       or i.imol['H2']/i.F_mol>0.1)])
    
    @property
    def F_vol_out(self): # exclude gases
        return sum([i.F_vol for i in self.outs if i.phase=='l' and i.F_mol and not (i.imol['CO2']/i.F_mol==1. 
                                                                       or i.imol['O2']/i.F_mol>0.1
                                                                       or i.imol['H2']/i.F_mol>0.1)])
    
    
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
        # water.imass['Water'] = acid.imass['SulfuricAcid'] / 0.05
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
        
        # total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solid_loading
        # water.imass['Water'] = total_mass - mixture.F_mass
        # water adjustment currently implemented in H_M202.specification
        
        mixture_out.mix_from([mixture, water, recycled_water])

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
#         Rxn('Glucose -> TAL',        'Glucose',   .9099),
#         # Rxn('Glucose -> 2,3-Butanediol',        'Glucose',   0.666),
#         Rxn('Glucose -> AceticAcid',               'Glucose',   0.07),
#         Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.02),
#         Rxn('Xylose -> TAL',       'Xylose',    0.9099),
#         Rxn('Xylose -> AceticAcid',       'Xylose',    0.07),
#         Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.02),
#         ])
        
#         self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
#         self.xylose_to_TAL_rxn = self.cofermentation_rxns[3]
        
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
        
#         # !!! Temporary alteration for 2,3-TAL production
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
    inoculum_ratio = 0.0 # this ratio not used
    
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
        

# ## CoFermentation
# @cost(basis='Fermenter size', ID='Fermenter', units='kg',
#       cost=10128000, S=(42607+443391+948+116)*(60+36),
#       CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Fermenter size', ID='Fermenter agitator', units='kg',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in ref [1])
#       # and total residence time (batch hydrolysis and fermentation)
#         # kW=268.452, 
#         kW=0, # overwritten; power utility calculated separately based on batch IBRL power per unit volume
#       cost=630000, 
#       S=(42607+443391+948+116)*(60+36),
#       CE=CEPCI[2009], n=1, BM=1.5)
# @cost(basis='Recirculation flow rate', ID='Recirculation pump', units='kg/hr',
#       # Scaling basis based on sum of all streams into fermenter
#       # (304, 306, 311, and 312 in ref [1])
#       kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)
# # Surge tank to hold the fermentation broth
# @cost(basis='Broth flow rate', ID='Surge tank', units='kg/hr',
#       cost=636000, S=425878, CE=CEPCI[2009], n=0.7, BM=1.8)
# @cost(basis='Broth flow rate', ID='Surge tank agitator', units='kg/hr',
#       kW=29.828, cost=68300, S=425878, CE=CEPCI[2009], n=0.5, BM=1.5)
# @cost(basis='Broth flow rate', ID='Surge tank pump', units='kg/hr',
#       kW=93.2125, cost=26800, S=488719, CE=CEPCI[2009], n=0.8, BM=2.3)
# class CoFermentation(Reactor):    
#     _N_ins = 4
#     _N_outs = 2
    
#     _units= {
#         **Reactor._units,
#         'Fermenter size': 'kg',
#         'Recirculation flow rate': 'kg/hr',
#         'Broth flow rate': 'kg/hr',
#         }    
#     auxiliary_unit_names = ('heat_exchanger')
    
#     # _F_BM_default = {
#     #         'Heat exchangers': 3.,}
    
#     # Co-Fermentation time (hr)
#     tau_batch_turnaround = 12 # in hr, the same as the seed train in ref [1]

#     tau_cofermentation = 120 # initial value; updated by spec.load_productivity
    
#     max_batch_reactor_volume = 0.075 # m3 # 75 L - IBRL batch pilot scale
    
#     max_continuous_reactor_volume = 3785.41178 # 1,000,000 gallon from ref [1]
#     CSL_loading = 10 # kg/m3
    
#     acetate_target_loading = 10 # g-AceticAcid-eq / L
    
    
#     def __init__(self, ID='', ins=None, outs=(), T=30+273.15,
#                  P=101325., 
#                  V_wf=0.4, length_to_diameter=0.6,
#                  kW_per_m3=0.02 * 0.7457 / 0.075, 
#                  mixing_intensity=None,
#                  wall_thickness_factor=1,
#                  vessel_material='Stainless steel 316',
#                  vessel_type='Vertical',
#                  mode='batch', feed_freq=1,
#                  mix_in_batch_mode=True,
#                  ):
#         Unit.__init__(self, ID, ins, outs)
#         # Same T for saccharificatoin and co-fermentation
#         self.T = T
#         self.P = P
#         self.V_wf = V_wf
#         self.length_to_diameter = length_to_diameter
#         self.mixing_intensity = mixing_intensity
#         self.kW_per_m3 = kW_per_m3
#         self.wall_thickness_factor = wall_thickness_factor
#         self.vessel_material = vessel_material
#         self.vessel_type = vessel_type
#         self.mode = mode
#         self.feed_freq = feed_freq
#         self.mix_in_batch_mode = mix_in_batch_mode
        
#         self.heat_exchanger = HXutility(None, None, None, T=T) 
        
        
#         self.acetate_hydrolysis_rxns = ParallelRxn([
#         #      Reaction definition            Reactant    Conversion
#         Rxn('SodiumAcetate + H2O -> AceticAcid + NaOH',        'SodiumAcetate',   1.-1e-4), 
#         ])
        
        
#         self.cofermentation_rxns = ParallelRxn([
#         #      Reaction definition            Reactant    Conversion
#         Rxn('Glucose -> 0.6667 TAL + 2 CO2',        'Glucose',   0.19), 
#         Rxn('Glucose -> 0.3 VitaminA',               'Glucose',   1e-8), # retinol
#         Rxn('Glucose -> 0.214286 VitaminD2',               'Glucose',   1e-8), # ergosterol
#         Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.05),
        
#         Rxn('Xylose -> 0.555583 TAL + 1.3334 CO2',       'Xylose',    0.19),
#         Rxn('Xylose -> 0.25 VitaminA',       'Xylose',    1e-8),
#         Rxn('Xylose -> 0.17856 VitaminD2',       'Xylose',    1e-8),
#         Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.05),
        
#         Rxn('AceticAcid -> 0.333333 TAL + H2O',       'AceticAcid',    0.19),
#         Rxn('AceticAcid -> 0.1 VitaminA',       'AceticAcid',    1e-8),
#         Rxn('AceticAcid -> 0.0714 VitaminD2',       'AceticAcid',    1e-8),
#         Rxn('AceticAcid -> 2 FermMicrobe',        'AceticAcid',    0.05),
        
#         ])
        
#         # oxygen currently excluded from balance
#         self.CO2_generation_rxns = ParallelRxn([
#             Rxn('Glucose -> 6CO2 + 6H2O', 'Glucose', 1.),
#             Rxn('Xylose -> 5CO2 + 5H2O', 'Xylose', 1.),
#             Rxn('AceticAcid -> 2CO2 + 2H2O', 'AceticAcid', 1.),
#             ])
        
#         self.sucrose_hydrolysis_rxns = ParallelRxn([
#         #      Reaction definition            Reactant    Conversion
#         Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
#         ])
        
        
#         self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
#         self.xylose_to_TAL_rxn = self.cofermentation_rxns[4]
#         self.acetate_to_TAL_rxn = self.cofermentation_rxns[8]
        
#         self.glucose_to_VitaminA_rxn = self.cofermentation_rxns[1]
#         self.xylose_to_VitaminA_rxn = self.cofermentation_rxns[5]
#         self.acetate_to_VitaminA_rxn = self.cofermentation_rxns[9]
        
#         self.glucose_to_VitaminD2_rxn = self.cofermentation_rxns[2]
#         self.xylose_to_VitaminD2_rxn = self.cofermentation_rxns[6]
#         self.acetate_to_VitaminD2_rxn = self.cofermentation_rxns[10]
        
#         self.glucose_to_microbe_rxn = self.cofermentation_rxns[3]
#         self.xylose_to_microbe_rxn = self.cofermentation_rxns[7]
#         self.acetate_to_microbe_rxn = self.cofermentation_rxns[11]
        
#         self.glucose_to_CO2_rxn = self.CO2_generation_rxns[0]
#         self.xylose_to_CO2_rxn = self.CO2_generation_rxns[1]
#         self.acetate_to_CO2_rxn = self.CO2_generation_rxns[2]

#     def _run(self):
#         feed, sugars, CSL, Acetate_spiking = self.ins
        
#         effluent, vapor = self.outs
#         effluent.mix_from([feed, sugars])
        
#         current_acetate_loading = effluent.imass['AceticAcid'] / effluent.F_vol
#         required_acetate_spiking = max(0, self.acetate_target_loading - current_acetate_loading)
#         Acetate_spiking.imass['SodiumAcetate'] = required_acetate_spiking * effluent.F_vol
        
#         effluent.mix_from([effluent, Acetate_spiking])
        
#         # ss = Stream(None)
#         # effluent.copy_like(feed)
#         effluent.T = vapor.T = self.T
#         CSL.imass['CSL'] = (sugars.F_vol+feed.F_vol) * self.CSL_loading 
        
#         self.acetate_hydrolysis_rxns(effluent.mol)
#         self.sucrose_hydrolysis_rxns(effluent.mol)
#         self.cofermentation_rxns(effluent.mol)
#         self.CO2_generation_rxns(effluent.mol)
        
#         vapor.imol['CO2'] = effluent.imol['CO2']
#         vapor.phase = 'g'
#         effluent.imol['CO2'] = 0
#         effluent.imass['CSL'] = 0
        
#         vapor.imol['CO2'] += CSL.get_atomic_flow('C')
        
#         self.effluent_titer = compute_TAL_titer(effluent)
        
#     def _design(self):
#         mode = self.mode
#         Design = self.design_results
#         Design.clear()

#         if mode == 'batch':
#             self._V_max = self.max_batch_reactor_volume
#             # Reactor._design(self)
#             tau_tot = self.tau_batch_turnaround + self.tau_cofermentation
#             Design['Fermenter size'] = self.outs[0].F_mass * tau_tot / self.V_wf
#             Design['Recirculation flow rate'] = self.F_mass_in
#             Design['Broth flow rate'] = self.outs[0].F_mass
        
#         elif mode=='continuous':
#             self._V_max = self.max_continuous_reactor_volume
#             Reactor._design(self)
#             self.tau_cofermentation = self.tau
#             # Include a backup fermenter for cleaning
#             Design['Number of reactors'] += 1
        
#         duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
#         mixed_feed = tmo.Stream()
#         mixed_feed.mix_from(self.outs)
#         mixed_feed.T=self.ins[0].T
#         # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
#         self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
#                                                             duty=duty,
#                                                             vle=False)
    
#     def _cost(self):
#         mode = self.mode
#         self.power_utility(0)
#         Design = self.design_results
#         purchase_costs = self.baseline_purchase_costs
#         purchase_costs.clear()
        
#         if mode == 'batch':
#             Unit._cost()
#             self._decorated_cost()
#             if self.mix_in_batch_mode: self.power_utility.consumption += self.kW_per_m3*(self.outs[0].F_vol*self.tau)
            
#         elif mode == 'continuous':
#             if not self.neutralization:
#                 self.vessel_material= 'Stainless steel 316'
#             Reactor._cost(self)
#             N_working = Design['Number of reactors'] - 1 # subtract the one back-up
#             # self.power_utility(self.kW_per_m3*Design['Single reactor volume']*N_working)
#             # self.power_utility.consumption += self.kW_per_m3*Design['Single reactor volume']*N_working
        
#             # No power need for the back-up reactor
#         # Reactor._cost(self)
#         # N_working = Design['Number of reactors'] - 1 # subtract the one back-up
        
#         hx = self.heat_exchanger
#         self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
#         self.heat_utilities += hx.heat_utilities
        
#     @property
#     def tau(self):
#         return self.tau_cofermentation
    
#     @tau.setter
#     def tau(self, new_tau):
#         self.tau_cofermentation = new_tau
        
#     @property
#     def mode(self):
#         return self._mode
#     @mode.setter
#     def mode(self, i):
#         if i.lower() in ('fed-batch', 'fedbatch', 'fed batch'):
#             raise ValueError('For fed-batch, set fermentation to "batch" and ' \
#                              'change feed_freq.')
#         elif i.lower() in ('batch', 'continuous'):
#             self._mode = i.lower()
#             self.feed_freq = 1
#         else:
#             raise ValueError(f'Mode can only be "batch" or "continuous", not {i}.')

#     @property
#     def feed_freq(self):
#         return self._feed_freq
#     @feed_freq.setter
#     def feed_freq(self, i):
#         if not (int(i)==i and i >0):
#             raise ValueError('feed_freq can only be positive integers.')
#         elif self.mode == 'continuous' and int(i) != 1:
#             raise ValueError('feed_freq can only be 1 for continuous mode.')
#         else:
#             self._feed_freq = int(i)
    
#     def mol_atom_in(self, atom):
#         return sum([stream.get_atomic_flow(atom) for stream in self.ins])
    
#     def mol_atom_out(self, atom):
#         return sum([stream.get_atomic_flow(atom) for stream in self.outs])
    

class BatchCoFermentation(BatchBioreactor):
    # Co-Fermentation time (hr)
    
    _N_ins = 6
    _N_outs = 2
    
    tau_cofermentation = 120 # initial value; updated by spec.load_productivity
    
    CSL_loading = 32.5 # g/L
    
    acetate_target_loading = 13.667 * (60.05196/82.033789) # g-AceticAcid-eq / L # 13.667 g-sodium acetate /L as in Markham et al. 2018
    
    air_m3_per_h_per_m3_reactor = 3.5*60/3 # 3.5 slpm for a 3L bioreactor; Markham et al. 2018 # used when aeration_rate_basis=='fixed rate basis'
    
    DO_saturation_concentration_kg_per_m3 = 7.8e-3 # 7.8 mg/L at 28 degrees C # used when aeration_rate_basis=='DO saturation basis'
    
    DO_saturation_target_level = 0.5 # 50% saturation from Markham et al. 2018 # used when aeration_rate_basis=='DO saturation basis'
    
    air_flow_rate_safety_factor_for_DO_saturation_basis = 1.
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.5 # !!!
    
    # autoselect_N  = True

    
    def __init__(self, ID='', ins=None, outs=(), 
                 T=28+273.15,
                 P=101325., 
                 tau=120, # initial value; updated by spec.load_productivity
                 V=3785.,
                 acetate_ID='SodiumAcetate',
                 aeration_rate_basis='fixed rate basis', # 'fixed rate basis' or 'DO saturation basis'
                 ):
        BatchBioreactor.__init__(self, ID, ins, outs, T=T, P=P, tau=tau, V=V)
        
        self.aeration_rate_basis = aeration_rate_basis
        
        self.acetate_ID = acetate_ID
        self.acetate_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('SodiumAcetate + H2O -> AceticAcid + NaOH',        'SodiumAcetate',   1.-1e-4), 
        ])
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 0.6667 TAL + 2 CO2',        'Glucose',   0.19), 
        Rxn('Glucose -> 0.3 VitaminA',               'Glucose',   1e-8), # retinol
        Rxn('Glucose + O2 -> CitricAcid + H2O',               'Glucose',   0.08856), # 2H+ excluded # from Markham et al.; 16 g/L citrate from 180 g/L glucose
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.05),
        
        Rxn('Xylose -> 0.555583 TAL + 1.3334 CO2',       'Xylose',    0.19),
        Rxn('Xylose -> 0.25 VitaminA',       'Xylose',    1e-8),
        Rxn('Xylose -> 0.8333 CitricAcid + H2O',               'Xylose',   1e-8),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.05),
        
        Rxn('AceticAcid -> 0.22218 TAL + 2CO2 + 6H2O',       'AceticAcid',    0.19),
        Rxn('AceticAcid -> 0.1 VitaminA',       'AceticAcid',    1e-8),
        Rxn('AceticAcid -> 0.0714 VitaminD2',       'AceticAcid',    1e-8),
        Rxn('AceticAcid -> 2 FermMicrobe',        'AceticAcid',    0.05),
        
        ])
        
        self.CO2_generation_rxns = ParallelRxn([
            Rxn('Glucose + 6O2 -> 6CO2 + 6H2O', 'Glucose', 1.-1e-2),
            Rxn('Xylose + 5O2 -> 5CO2 + 5H2O', 'Xylose', 1.-1e-2),
            Rxn('AceticAcid + 2O2 -> 2CO2 + 2H2O', 'AceticAcid', 1.-1e-2),
            ])
        
        self.sucrose_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
        ])
        
        
        self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
        self.xylose_to_TAL_rxn = self.cofermentation_rxns[4]
        self.acetate_to_TAL_rxn = self.cofermentation_rxns[8]
        
        self.glucose_to_VitaminA_rxn = self.cofermentation_rxns[1]
        self.xylose_to_VitaminA_rxn = self.cofermentation_rxns[5]
        self.acetate_to_VitaminA_rxn = self.cofermentation_rxns[9]
        
        self.glucose_to_CitricAcid_rxn = self.cofermentation_rxns[2]
        self.xylose_to_CitricAcid_rxn = self.cofermentation_rxns[6]
        self.acetate_to_VitaminD2_rxn = self.cofermentation_rxns[10]
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[3]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[7]
        self.acetate_to_microbe_rxn = self.cofermentation_rxns[11]
        
        self.glucose_to_CO2_rxn = self.CO2_generation_rxns[0]
        self.xylose_to_CO2_rxn = self.CO2_generation_rxns[1]
        self.acetate_to_CO2_rxn = self.CO2_generation_rxns[2]

    def _run(self):
        feed, seed, CSL, Acetate_spiking, yeast_extract, air = self.ins
        for i in [CSL, Acetate_spiking]: i.empty()
        
        vapor, effluent = self.outs
        
        effluent.mix_from([feed, seed])
        
        current_acetate_loading = effluent.imass[self.acetate_ID] / effluent.F_vol
        required_acetate_spiking = max(0, self.acetate_target_loading - current_acetate_loading)
        Acetate_spiking.imass[self.acetate_ID] = required_acetate_spiking * effluent.F_vol
        
        # air.P = self.air_pressure
        air.imol['N2'] = 0.79
        air.imol['O2'] = 0.21
        
        if self.aeration_rate_basis == 'DO saturation basis':
            self.air_exit_F_mol_needed = (1./0.21) * (1/32.) * self.air_flow_rate_safety_factor_for_DO_saturation_basis * self.DO_saturation_concentration_kg_per_m3 * self.DO_saturation_target_level\
                *(seed.F_vol+feed.F_vol)
            
            air.F_mol = 1e6 # initial value; updated after reactions
        
        elif self.aeration_rate_basis == 'fixed rate basis':
            air.F_vol = self.air_m3_per_h_per_m3_reactor * (seed.F_vol+feed.F_vol) * self.tau
        
        else: raise RuntimeError(f"Unsupported aeration_rate_basis ({self.aeration_rate_basis}); must be 'fixed rate basis' or 'DO saturation basis'.")
        
        CSL.imass['CSL'] = (seed.F_vol+feed.F_vol) * self.CSL_loading 
        
        effluent.mix_from([effluent, Acetate_spiking, air])
        effluent.T = vapor.T = self.T
        
        self.acetate_hydrolysis_rxns(effluent.mol)
        self.sucrose_hydrolysis_rxns(effluent.mol)
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        
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
        
        vapor.imol['CO2', 'O2', 'N2'] = effluent.imol['CO2', 'O2', 'N2']
        vapor.phase = 'g'
        effluent.imol['CO2', 'O2', 'N2'] = 0, 0, 0
        effluent.imass['CSL'] = 0
        
        vapor.imol['CO2'] += CSL.get_atomic_flow('C')
        
        # self.effluent_titer = compute_TAL_titer(effluent)
        # self.effluent = effluent
        
    @property
    def effluent_titer(self):
        return compute_TAL_titer(self.effluent)


# class AeratedTALCoFermentation(AeratedBioreactor):
    
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
    
    def __init__(self, ID='', ins=None, outs=(), T=30+273.15, ferm_ratio=0.95):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.ferm_ratio = ferm_ratio
        self.heat_exchanger = HXutility(None, None, None, T=T) 
        
        
        self.acetate_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('SodiumAcetate + H2O -> AceticAcid + NaOH',        'SodiumAcetate',   1.-1e-4), 
        ])
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 0.6667 TAL + 2 CO2',        'Glucose',   0.19), 
        Rxn('Glucose -> 0.3 VitaminA',               'Glucose',   1e-8), # retinol
        Rxn('Glucose -> CitricAcid + H2O',               'Glucose',   0.08856), # O2 and 2H+ excluded # from Markham et al.; 16 g/L citrate from 180 g/L glucose
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.05),
        
        Rxn('Xylose -> 0.555583 TAL + 1.3334 CO2',       'Xylose',    0.19),
        Rxn('Xylose -> 0.25 VitaminA',       'Xylose',    1e-8),
        Rxn('Xylose -> 0.8333 CitricAcid + H2O',       'Xylose',    1e-8),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.05),
        
        Rxn('AceticAcid -> 0.22218 TAL + 2CO2 + 6H2O',       'AceticAcid',    0.19),
        Rxn('AceticAcid -> 0.1 VitaminA',       'AceticAcid',    1e-8),
        Rxn('AceticAcid -> 0.0714 VitaminD2',       'AceticAcid',    1e-8),
        Rxn('AceticAcid -> 2 FermMicrobe',        'AceticAcid',    0.05),
        
        ])
        
        # oxygen currently excluded from balance
        self.CO2_generation_rxns = ParallelRxn([
            Rxn('Glucose -> 6CO2 + 6H2O', 'Glucose', 1.),
            Rxn('Xylose -> 5CO2 + 5H2O', 'Xylose', 1.),
            Rxn('AceticAcid -> 2CO2 + 2H2O', 'AceticAcid', 1.),
            ])
        
        self.sucrose_hydrolysis_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Sucrose + H2O -> 2 Glucose',        'Sucrose',   1.-1e-4), 
        ])
        
        
        self.glucose_to_TAL_rxn = self.cofermentation_rxns[0]
        self.xylose_to_TAL_rxn = self.cofermentation_rxns[4]
        self.acetate_to_TAL_rxn = self.cofermentation_rxns[8]
        
        self.glucose_to_VitaminA_rxn = self.cofermentation_rxns[1]
        self.xylose_to_VitaminA_rxn = self.cofermentation_rxns[5]
        self.acetate_to_VitaminA_rxn = self.cofermentation_rxns[9]
        
        self.glucose_to_CitricAcid_rxn = self.cofermentation_rxns[2]
        self.xylose_to_CitricAcid_rxn = self.cofermentation_rxns[6]
        self.acetate_to_CitricAcid_rxn = self.cofermentation_rxns[10]
        
        self.glucose_to_microbe_rxn = self.cofermentation_rxns[3]
        self.xylose_to_microbe_rxn = self.cofermentation_rxns[7]
        self.acetate_to_microbe_rxn = self.cofermentation_rxns[11]
        
        self.glucose_to_CO2_rxn = self.CO2_generation_rxns[0]
        self.xylose_to_CO2_rxn = self.CO2_generation_rxns[1]
        self.acetate_to_CO2_rxn = self.CO2_generation_rxns[2]

    def _run(self):
        feed, = self.ins
        
        effluent, vapor = self.outs
        effluent.mix_from([feed,])
        
        effluent.T = vapor.T = self.T
        
        self.acetate_hydrolysis_rxns(effluent.mol)
        self.sucrose_hydrolysis_rxns(effluent.mol)
        self.cofermentation_rxns(effluent.mol)
        self.CO2_generation_rxns(effluent.mol)
        
        vapor.imol['CO2'] = effluent.imol['CO2']
        vapor.phase = 'g'
        effluent.imol['CO2'] = 0
        effluent.imass['CSL'] = 0
        

    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass
        Design['Seed fermenter size'] = self.outs[0].F_mass * self.tau_batch
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        mixed_feed = tmo.Stream()
        mixed_feed.mix_from(self.outs)
        mixed_feed.T=self.ins[0].T
        # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
                                                            duty=duty,
                                                            vle=False)
        
    def _cost(self):
        Unit._cost()
        self._decorated_cost()
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
    
    def mol_atom_in(self, atom):
        return sum([stream.get_atomic_flow(atom) for stream in self.ins])
    
    def mol_atom_out(self, atom):
        return sum([stream.get_atomic_flow(atom) for stream in self.outs])
    
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

class TALCrystallizer(BatchCrystallizer):
    
    _SA_vol_per_mass = 0.0008252419812169215
    def __init__(self, ID='', ins=None, outs=(), 
                 target_recovery=0.6,
                 thermo=None,
                 tau=6, N=5, 
                 V=None, T=273.15 + 1.,
                 basis='water solubility',
                 get_mol_TAL_dissolved_given_T_and_mol_water=None,
                 fixed_operating_T=273.15 + 1.,
                 Nmin=2, Nmax=36, 
                 T_range=(274., 372.5),
                 vessel_material='Stainless steel 316',
                 output_conc_multiplier=1.,
                 kW=0.00746):
        
        BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                     tau, N, V, T,
                     Nmin, Nmax, vessel_material,
                     kW)
        self.target_recovery = target_recovery
        self.T_range = T_range
        self.tau = tau
        self.basis = basis
        self.get_mol_TAL_dissolved_given_T_and_mol_water = get_mol_TAL_dissolved_given_T_and_mol_water
        self.fixed_operating_T = fixed_operating_T
        self.output_conc_multiplier = output_conc_multiplier
        
    
    def _run(self):
        in_stream, = self.ins
        out_stream, = self.outs
        target_recovery = self.target_recovery
        
        out_stream.copy_like(in_stream)
        # out_stream.sle(T=self.T, solute='SuccinicAcid')
        
        tot_TAL = in_stream.imol['TAL']
        
        if self.fixed_operating_T:
            self.T = self.fixed_operating_T
            
        out_stream.T = self.T
        
        TAL_solubility = self.get_mol_TAL_dissolved_given_T_and_mol_water(out_stream.T, out_stream.imol['Water'])
        out_stream.phases = ('s', 'l')
        TAL_dissolved = min(TAL_solubility, tot_TAL)
        out_stream.imol['l', 'TAL'] = TAL_dissolved
        out_stream.imol['s', 'TAL'] = max(0, tot_TAL - TAL_dissolved)




##




class SorbicAcidCrystallizer(BatchCrystallizer):
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T=250., crystal_SA_purity=0.98, melt_water_purity=0.90,
                 order=None):
        BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                       tau=5, V=1e6, T=T)
        self.melt_AcTAG_purity = melt_AcTAG_purity
        self.crystal_TAG_purity = crystal_TAG_purity

    @property
    def Hnet(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        if 's' in feed.phases:
            H_in = - sum([i.Hfus * j for i,j in zip(self.chemicals, feed['s'].mol) if i.Hfus])
        else:
            H_in = 0.
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out - H_in
        
    def _run(self):
        outlet = self.outs[0]
        outlet.phases = ('s', 'l')
        crystal_TAG_purity = self.crystal_TAG_purity
        melt_AcTAG_purity = self.melt_AcTAG_purity
        feed = self.ins[0]
        TAG, AcTAG = feed.imass['TAG', 'AcTAG'].value
        total = TAG + AcTAG
        minimum_melt_purity = AcTAG / total
        minimum_crystal_purity = TAG / total
        outlet.empty()
        if crystal_TAG_purity < minimum_crystal_purity:
            outlet.imol['s'] = feed.mol
        elif melt_AcTAG_purity < minimum_melt_purity:
            outlet.imol['l'] = feed.mol
        else: # Lever rule
            crystal_AcTAG_purity = (1. - crystal_TAG_purity)
            melt_fraction = (minimum_melt_purity - crystal_AcTAG_purity) / (melt_AcTAG_purity - crystal_AcTAG_purity)
            melt = melt_fraction * total
            AcTAG_melt = melt * melt_AcTAG_purity
            TAG_melt = melt - AcTAG
            outlet.imass['l', ('AcTAG', 'TAG')] = [AcTAG_melt, TAG_melt]
            outlet.imol['s'] = feed.mol - outlet.imol['l']
        outlet.T = self.T
        
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


#TODO:
# 1. Find a source for recovery - the only thing I found for this is in Woods 2007 section 5.3.1 "allow both phases to have > 20% of the diameter and no less than 0.2 m to ensure that the exit phases do not become cross-contaminated"
# 2. Do we need a backup reactor? - need to discuss
# 3. Need electricity? - I think so but I don't know how to account for this? Let's discuss.
# 4. Look into n - mothi's thesis uses n=1, I canot find where I had previously found a value of 0.84
# 5. Look into BM (bare module factor) - Sieder et al 2016 Table 16.11 use 4.16 for  a vertical pressure vessel
# I also had in my notes to check the decanter usable % and I believe that is 90% but it already accounted for in the height sizing of adding an additional 10% to the height.

# @cost(basis='Flow rate', ID='Decanter', units='m3/s',
#        cost=190000, S=0.0012, CE=CEPCI[2008], n=1, BM=2.3)
@cost(basis='Flow rate', ID='Decanter', units='kg/hr',
      cost=588000, S=145930, CE=CEPCI[2013], n=0.5, BM=2)
class Decantation(Unit):
    '''
    A decanter to separate the sorbic acid crystal from the broth.

    Using Stokes’ Law to calculate settling velocity and the desired volumetric flow rate
    of the solution into the decanter.

    Parameters
    ----------
    recovery : float
        Recovery of the crystal.
    V_wf : float
        Working fraction of the decanter.

    References
    ----------
    1. Towler, G., & Sinnott, R. K. (2008). Chemical engineering design: principles, practice and economics of plant and process design.
    2. Peters, M. S., Timmerhaus, K. D., West, R. E., Timmerhaus, K., & West, R. (2003). Plant design and economics for chemical engineers (Vol. 4).
    3. Density: https://pubchem.ncbi.nlm.nih.gov/compound/Sorbic-acid#section=Solubility
    4. Seider, W. D., Lewin, D. R., Seader, J. D., Widagdo, S., Gani, R., Ming Ng, K. (2017) Product and Process Design Principles: Synthesis, Analysis, and Evaluation (4th Edition) Wiley.
    5. M. Viswanathan. (2019) Process generalizations and rules of thumb for scaling up biobased processes.
    6. Woods, D. R. (2007) Rules of Thumb in Engineering Practice, Wiley.
    '''

    _N_ins = 1
    _N_outs = 2

    _units={
        'Settling velocity': 'm/s',
        'Area': 'm2',
        'Diameter': 'm',
        'Length': 'm',
        'Volume': 'm3',
        }

    def __init__(self, ID='', ins=None, outs=(), forced_recovery=None, T=293.15):
        Unit.__init__(self, ID, ins, outs)
        self.forced_recovery = forced_recovery
        self.T = T
        
    def get_SA_solubility(self, T): # mol SA/ (mol SA + mol water)
        return exp(-157.137 + 3505.81/T + 24.1015*log(T)) # Modified Apelblat Equation
    
    def _run(self):
        feed = self.ins[0]
        dil_SA, crystal = self.outs
        dil_SA.copy_like(feed)
        dil_SA.T = crystal.T = self.T
        if self.forced_recovery:
            crystal.imass['SorbicAcid'] = feed.imass['SorbicAcid']*self.forced_recovery
            dil_SA.imass['SorbicAcid'] -= crystal.imass['SorbicAcid']
        else:
            sol_SA = self.get_SA_solubility(self.T)
            still_dissolved_SA = sol_SA * feed.imol['Water'] / (1-sol_SA)
            dil_SA.imass['SorbicAcid'] = still_dissolved_SA
            crystal.imass['SorbicAcid'] = feed.imass['SorbicAcid'] - dil_SA.imass['SorbicAcid']
        

    def _design(self):
        feed = self.ins[0]
        D = self.design_results
        D['Flow rate'] = self.F_mass_in
        # D['Flow rate'] = feed.F_vol/3600. # m3/s
        D['Settling velocity'] = v_settling = \
            (0.00015**2)*9.81*(1204-998.19)/(18*0.001) # m/s
        D['Area'] = A_decanter = feed.F_vol*0.00028/v_settling # m2
        D['Diameter'] = D_decanter = (A_decanter/pi)**(1/2) # m
        D['Length'] = L_decanter = 2.2*D_decanter # m
        D['Volume'] = (pi*D_decanter*L_decanter)/2 # multiply by 2.2 and then divided by 2 makes a 10% extra headspace



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
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol = self.multi_stream.imol['g']	
        biogas.phase = 'g'	
        liquid_mol = self.multi_stream.imol['l']	
        treated_water.mol = liquid_mol * self.split	
        sludge.mol = liquid_mol - treated_water.mol	
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
        
# class TALLiquidsSplitSettler():
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

# # For storage of lime used in separation and waste treatment,
# # design copied from cornstover biorefinery in Aden et al.
# # Base flow from stream 27 of Aden et al.
# @cost(basis='Flow rate', ID='Storage bin', units='kg/hr',
#       cost=136370, S=2395, CE=386.5, n=0.46, BM=1.3)
# # Cost not scaled, thus used n=0
# # Power usage scaled based on M104 in Humbird et al., 
# # (truck dumper hopper for feedstock handling) 
# @cost(basis='Flow rate', ID='Feeder', units='kg/hr',
#       kW=37.285/3500*140, cost=3900, S=2395, CE=386.5, n=0, BM=1.3)
# # Power usage scaled based on M-106 in Humbird et al.
# # (dust collection system for feedstock handling) 
# @cost(basis='Flow rate', ID='Unloading blower', units='kg/hr',
#       kW=18.6425*7425/8500, cost=99594, S=2395, CE=389.5, n=0.5, BM=1.4)
# @cost(basis='Flow rate', ID='Dust vent baghouse', units='kg/hr',
#       cost=140707, S=2395, CE=386.5, n=1, BM=1.5)
# class LimeStorageBin(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=803000, S=8343, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=15000, S=8343, CE=CEPCI[2009], n=0.8, BM=3.1)
class FireWaterTank(Unit): pass

class SugarsMixer(Mixer):
    _N_ins = 3
    
# Modified from bst.units.StorageTank, which won't simulate for 0 flow 	
class TALStorageTank(StorageTank):
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.purchase_costs['Tanks'] = 0
        else: StorageTank._cost(self)

# Modified from bst.units.Pump, which won't simulate for 0 flow 
class TALPump(Pump):
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
        

class TAL_Separation(Unit):
    _N_outs = 2
    def _run(self): pass
    def _cost(self): pass


class Adsorption(PressureVessel):
    '''
    Removal of polar compounds (e.g., amino acids, ions) from the stream. The unit
    is modeled based on description in [1]_. Freundlich isotherm is assumed for TAL
    based on the constants of catechol [2]_, whose structure is similar to TAL.
    
    .. math:: q_e = K_F * C_e^(1/n)
    
    Parameters
    ----------
    adsorbents : str
        Breakthrough time in hr
    breakthrough_time : float
        Breakthrough time in hr
    polars : tuple
        IDs of polar compounds adsorbed to the column
    K_F : float
        Freundlich constant, default at 42.40
    1/n : float
        Freundlich exponent, default at 0.205
    
    References
    ----------
    .. [1] Viswanathan, Process Generalizations and Rules of Thumb for Scaling
        up Biobased Processes. Graduate Theses and Dissertations 2019.
        https://lib.dr.iastate.edu/etd/17115
    .. [2] Kumar et al., Adsorption of Resorcinol and Catechol on Granular
        Activated Carbon: Equilibrium and Kinetics. Carbon 2003, 41 (15), 3015–3025.
        https://doi.org/10.1016/S0008-6223(03)00431-7
    .. [3] Seider, Warren D., et al. (2017). "Cost Accounting and Capital Cost
        Estimation". In Product and Process Design Principles: Synthesis,
        Analysis, and Evaluation (pp. 470, 481). New York: Wiley.
        
    '''
    _N_ins = 2
    _N_outs = 3
    
    #!!! Should choose silica gel
    # 574, 613
    
    # in $/ft3
    adsorbent_cost = {
        'Activated alumina': 72,
        'Activated carbon': 41,
        'Silica gel': 210,
        'Molecular sieves': 85,
        }
    
    def __init__(self, ID='', ins=None, outs=(), *,
                 adsorbent='Activated carbon',
                 P=101325, tau=0.5, V_wf=0.8,
                 length_to_diameter=2, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        self.ins = ins
        self.outs = outs
        self.adsorbent = adsorbent

        self.tau = tau
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.F_M = {}
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
    
    def _cost(self):
        super()._cost()
        Cost = self.purchase_costs
        Design = self.design_results()
        
        # 10 Hp/1,000 gal, or 7.457 kW/3.78541 m3 = 1.97 kW/m3
        hp = 10*Design['Volume']/3.78541
        
        Cost['Agitator'] = 4105*hp**0.57
        
        kW = 1.97*Design['Volume']
        self.power_utility(kW)


class Adsorption_and_Centrifugation(Unit):
    _N_outs = 2
    def _run(self): pass
    def _cost(self): pass


class HydrogenationEstersReactor(Reactor):
    """
    A hydrogenation reactor to produce esters from TAL.
    """
    _N_ins = 5
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default,
            'PdC catalyst': 1}
    mcat_frac = 0.5 # fraction of catalyst by weight in relation to the reactant (TAL)
    mono_di_hydroxy_esters_conversion_ratio = 0.5
    TAL_to_esters_and_DHL_conversion = 0.9
    esters_DHL_conversion_ratio = 0.5
    TAL_to_esters_conversion = esters_DHL_conversion_ratio * TAL_to_esters_and_DHL_conversion
    # catalyst_deactivation_k = 2e-3/10
    catalyst_replacements_per_year = 0. # deprecated variable
    hydrogenation_rxns = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 2H2 + Ethanol -> Ethyl_5_hydroxyhexanoate',         'TAL',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),
            Rxn('TAL + H2 + Ethanol + H2O -> Ethyl_3_5_dihydroxyhexanoate',         'TAL',   TAL_to_esters_conversion*(1.-mono_di_hydroxy_esters_conversion_ratio)),
            Rxn('TAL + 2H2 -> HMTHP',         'TAL',   1e-6), # conversion from Huber group experimental data
            # Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.-1e-5),
            # ])
    # byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 3H2 -> DHL + H2O',         'TAL',   0.5-1e-3), # conversion from Huber group experimental data
            Rxn('TAL + H2 -> HMDHP',         'TAL',   1e-3),  # conversion from Huber group experimental data
            Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.),  # conversion from Huber group experimental data
            Rxn('HMTHP + Ethanol -> Ethyl_5_hydroxyhexanoate',         'HMTHP',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),  # conversion from Huber group experimental data
            Rxn('DHL + Ethanol -> Ethyl_3_5_dihydroxyhexanoate',         'DHL',   TAL_to_esters_conversion*mono_di_hydroxy_esters_conversion_ratio),  # conversion from Huber group experimental data
            ])
    
    TAL_to_esters_rxns = hydrogenation_rxns[:2]
    TAL_to_DHL_rxn = hydrogenation_rxns[3]
    byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + H2O -> Pentenone + CO2',         'TAL',   1.-1e-5),
           ])
    TAL_conversion_rxns = [i for i in hydrogenation_rxns if i.reactant=='TAL']
    
    def __init__(self, ID, ins, outs, 
                 tau = (7./3.) * 2., # from Huber group
                 T=323., # from Huber group
                 vessel_material='Stainless steel 316', **args):
        Reactor.__init__(self, ID, ins, outs, tau=tau, vessel_material=vessel_material)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        
    def _run(self):
        feed, recycle, reagent, fresh_catalyst, recovered_catalyst = self.ins
        effluent, vented_gas, spent_catalyst = self.outs
        
        # effluent.imol['HMDHP'] += 1e-10
        
        H2_mol_needed_for_HMTHP_formation = (sum(feed.imol['TAL', 'HMDHP']) + sum(recycle.imol['TAL', 'HMDHP']))*2*self.hydrogenation_rxns[0].X
        reagent.imol['H2'] = 20*H2_mol_needed_for_HMTHP_formation
        reagent.phase = 'g'
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle, reagent])
        
        # effluent.imol['HMTHP'] -= 1e-10
        
        effluent.phase = 'l'
        
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.hydrogenation_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        extra_H2_mol_to_exclude = effluent.imol['H2']
        
        vented_gas.phase = 'g'
        vented_gas.imol['CO2'] = effluent.imol['CO2']
        # vented_gas.imol['H2'] = effluent.imol['H2']
        
        effluent.imol['CO2'] = 0.
        # effluent.imol['H2'] = 0.
        
        reagent.imol['H2'] -= extra_H2_mol_to_exclude
        effluent.imol['H2'] = 0.
        
        for i in self.outs:
            i.T = self.T
        
        spent_catalyst.imass['Pd'] =\
            self.mcat_frac * self.ins[0].imass['TAL'] * self.tau * self.catalyst_replacements_per_year/ 7884.0
        # spent_catalyst.imass['Pd'] = fresh_catalyst.imass['Pd'] =\
            
            # self.catalyst_deactivation_k * (feed.imol['TAL'] - effluent.imol['TAL'])
            # self.catalyst_deactivation_k * reagent.imol['H2']
            
            # effluent.imol['Ethyl_5_hydroxyhexanoate'])
            # effluent.imol['Ethyl_3_5_dihydroxyhexanoate'])
            # effluent.imol['DHL'])
            
        
            # Alternative 
            # reagent.imol['H2']
        
        # catalyst leaves in the effluent (recovered downstream)
        effluent.imass['Pd'] = self.mcat_frac * self.ins[0].imass['TAL']
        fresh_catalyst.imass['Pd'] = max(0, effluent.imass['Pd'] - recovered_catalyst.imass['Pd'])
        
    def _design(self):
        Reactor._design(self)
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        # mixed_feed = tmo.Stream()
        # mixed_feed.mix_from(self.ins)
        # self.heat_exchanger.simulate_as_auxiliary_exchanger(duty, mixed_feed)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=self.ins, duty=duty)
        
    def _cost(self):
        Reactor._cost(self)
        self.purchase_costs['PdC catalyst'] =\
            self.mcat_frac * self.ins[0].imass['TAL'] * self.tau * price['PdC']
            



class HydrogenationReactor(Reactor):
    """
    A hydrogenation reactor.
    """
    _N_ins = 5
    _N_outs = 2
    auxiliary_unit_names = ('heat_exchanger')
    _F_BM_default = {**Reactor._F_BM_default,
            'Heat exchangers': 3.,
            'Ni-SiO2 catalyst': 1.}
    mcat_frac = 0.5 # fraction of catalyst by weight in relation to the reactant (TAL)
    hydrogenation_rxns = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 2H2 -> HMTHP',         'TAL',   0.968), # conversion from Huber group experimental data
            # Rxn('HMDHP + H2 -> HMTHP',         'HMDHP',   1.-1e-5)
            ])
    byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('TAL + 3H2 -> DHL + H2O',         'TAL',   1.-1e-5), # conversion from Huber group experimental data
            # Rxn('TAL + H2 -> HMDHP',         'TAL',   1-0.2125),  # conversion from Huber group experimental data
            ])
    
    TAL_to_HMTHP_rxn = hydrogenation_rxns[0]
    
    spent_catalyst_replacements_per_year = 1. # number of times the entire catalyst_weight is replaced per year
    
    def __init__(self, ID, ins, outs, 
                 tau = 17., # from Huber group
                 T=100. + 273.15, # from Huber group
                 P=3e6, # 30 bar # from Huber group
                 vessel_material='Stainless steel 316',
                 NiSiO2_catalyst_price=price['Ni-SiO2'],
                 **args):
        Reactor.__init__(self, ID, ins, outs, tau=tau, P=P, vessel_material=vessel_material)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.NiSiO2_catalyst_price = NiSiO2_catalyst_price
        
    def _run(self):
        feed, recycle, reagent, recovered_catalyst, fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        effluent.empty()
        effluent.phase = 'l'
        # effluent.imol['HMDHP'] += 1e-10
        tau = self.tau
        
        self.catalyst_weight = cat_weight = self.mcat_frac * sum(i.imass['TAL'] for i in self.ins) * tau
        
        H2_mol_needed_for_HMTHP_formation = 2 * (sum(feed.imol['TAL', 'HMDHP']) + sum(recycle.imol['TAL', 'HMDHP']))
        # *self.hydrogenation_rxns[0].X
        
        reagent.imol['H2'] = 5*H2_mol_needed_for_HMTHP_formation
        reagent.phase = 'g'
        
        effluent.mix_from([feed, recycle, reagent,])
        
        effluent.phase = 'l'
        
        self.hydrogenation_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        extra_H2_mol_to_exclude = effluent.imol['H2']
        reagent.imol['H2'] -= extra_H2_mol_to_exclude
        effluent.imol['H2'] = 0.
        
        # 
        req_cat_mass_flow = cat_weight/tau
        current_cat_mass_flow = recovered_catalyst.imass['NiSiO2']
        spent_catalyst_mass_flow = self.spent_catalyst_replacements_per_year*cat_weight/self.system.TEA.operating_hours # kg/h
        current_cat_mass_flow-=spent_catalyst_mass_flow
        
        spent_catalyst.phase = 's'
        spent_catalyst.imass['NiSiO2'] = spent_catalyst_mass_flow
        
        fresh_catalyst.phase = 's'
        fresh_catalyst.imass['NiSiO2'] = fresh_cat_mass_flow = max(0, req_cat_mass_flow - current_cat_mass_flow)
        
        effluent.T = self.T
        effluent.P = self.P
        effluent.phases = ('l', 's')
        effluent.imass['s', 'NiSiO2'] = fresh_cat_mass_flow+current_cat_mass_flow
        

        
    def _design(self):
        Reactor._design(self)
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        mixed_feed = tmo.Stream()
        mixed_feed.mix_from(self.outs)
        mixed_feed.T=self.ins[0].T
        # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
                                                            duty=duty,
                                                            vle=False)
    
    def _cost(self):
        super()._cost()
        
        self.purchase_costs['Ni-SiO2 catalyst'] = self.catalyst_weight * self.NiSiO2_catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
    
class DehydrationReactor(Reactor):
    """
    A dehydration reactor.
    """
    _N_ins = 4
    _N_outs = 2
    
    auxiliary_unit_names = ('heat_exchanger')
    
    _F_BM_default = {**Reactor._F_BM_default,
            'Heat exchangers': 3.,
            'Amberlyst-70 catalyst': 1}
    mcat_frac = 0.5 # fraction of catalyst by weight in relation to the reactant (TAL)
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('HMTHP -> PSA',         'HMTHP',   0.842) # conversion from Huber group experimental data
                ])
    byproduct_formation_rxns  = ParallelRxn([
            #   Reaction definition   Reactant   Conversion
            Rxn('HMTHP -> DHL + H2O',         'HMTHP',   1.-1e-5), # conversion from Huber group experimental data
            # Rxn('TAL + H2 -> HMDHP',         'TAL',   1-0.2125),  # conversion from Huber group experimental data
            ])
    HMTHP_to_PSA_rxn = dehydration_rxns[0]
    
    spent_catalyst_replacements_per_year = 5. # number of times the entire catalyst_weight is replaced per year
    
    def __init__(self, ID, ins, outs, 
                 tau = 17.8, # from Huber group
                 T=100. + 273.15, # from Huber group
                 P=3e6, # 30 bar # from Huber group
                 vessel_material='Stainless steel 316',
                 Amberlyst70_catalyst_price=price['Amberlyst-70'],
                 **args):
        Reactor.__init__(self, ID, ins, outs, tau=tau, P=P, vessel_material=vessel_material)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T) 
        self.Amberlyst70_catalyst_price = Amberlyst70_catalyst_price
    
    def _run(self):
        feed, recycle, recovered_catalyst, fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        effluent.empty()
        effluent.phase = 'l'
        
        tau = self.tau
        self.catalyst_weight = cat_weight = self.mcat_frac * sum(i.imass['HMTHP'] for i in self.ins) * tau
        
        effluent.mix_from([feed, recycle])
        
        self.dehydration_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        req_cat_mass_flow = cat_weight/tau
        current_cat_mass_flow = recovered_catalyst.imass['Amberlyst70_']
        spent_catalyst_mass_flow = self.spent_catalyst_replacements_per_year*cat_weight/self.system.TEA.operating_hours # kg/h
        current_cat_mass_flow-=spent_catalyst_mass_flow
        
        spent_catalyst.phase = 's'
        spent_catalyst.imass['Amberlyst70_'] = spent_catalyst_mass_flow
        
        fresh_catalyst.phase = 's'
        fresh_catalyst.imass['Amberlyst70_'] = fresh_cat_mass_flow = max(0, req_cat_mass_flow - current_cat_mass_flow)
        
        effluent.T = self.T
        effluent.P = self.P
        effluent.phases = ('l', 's')
        effluent.imass['s', 'Amberlyst70_'] = fresh_cat_mass_flow+current_cat_mass_flow
        
    def _design(self):
        Reactor._design(self)
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        mixed_feed = tmo.Stream()
        mixed_feed.mix_from(self.outs)
        mixed_feed.T=self.ins[0].T
        # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
                                                            duty=duty,
                                                            vle=False)
    
    def _cost(self):
        super()._cost()
        self.purchase_costs['Amberlyst-70 catalyst'] = self.catalyst_weight * self.Amberlyst70_catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities

            

class RingOpeningHydrolysisReactor(Reactor):
    """
    A ring-opening and hydrolysis reactor.
    """
    _N_ins = 3
    _N_outs = 1
    auxiliary_unit_names = ('heat_exchanger')
    
    _F_BM_default = {**Reactor._F_BM_default,
            'Heat exchangers': 3.,}
    
    ring_opening_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> SA',         'PSA',   0.99) # conversion from Huber group
                ])
    
    PSA_to_SA_rxn = ring_opening_rxns[0]
    hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('SA + KOH -> KSA + H2O',         'SA',   0.99) # conversion from Huber group 
                ])
    SA_to_KSA_rxn = hydrolysis_rxns[0]
    byproduct_formation_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> 0.2PolyPSA',         'PSA',   1.-1e-5) # assumed
                ])
    
    def __init__(self, ID, ins, outs, 
                 tau = 16., # from Huber group
                 T=130. + 273.15, # from Huber group
                 P=3e6, # 30 bar # from Huber group
                 vessel_material='Stainless steel 316', **args):
        Reactor.__init__(self, ID, ins, outs, tau=tau, P=P, vessel_material=vessel_material)
        self.T = T
        self.heat_exchanger = hx = HXutility(None, None, None, T=T)
    
    def _run(self):
        feed, recycle, reagent = self.ins
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.ring_opening_rxns(effluent.mol)
        
        reagent.imol['KOH'] = max(0, effluent.imol['SA'] - effluent.imol['KOH'])
        
        effluent.mix_from([effluent, reagent])
        
        self.hydrolysis_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
        
        effluent.T = self.T
        effluent.P = self.P

    def _design(self):
        Reactor._design(self)
        duty = sum([i.H for i in self.outs]) - sum([i.H for i in self.ins])
        mixed_feed = tmo.Stream()
        mixed_feed.mix_from(self.outs)
        mixed_feed.T=self.ins[0].T
        # mixed_feed.vle(T=mixed_feed.T, P=mixed_feed.P)
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(mixed_feed,), 
                                                            duty=duty,
                                                            vle=False)
    
    def _cost(self):
        super()._cost()
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
            

            
class RingOpeningReactor(Reactor):
    """
    A ring-opening reactor.
    """
    _N_ins = 2
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'RaneyNi catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    ring_opening_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> SA',         'PSA',   1.) # conversion from Chia et al. 2012
                ])
    byproduct_formation_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('PSA -> 0.2PolyPSA',         'PSA',   1.-0.089) # conversion from Chia et al. 2012
                ])
    PSA_to_SA_rxn = ring_opening_rxns[0]
    
    def _run(self):
        feed, recycle = self.ins
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.ring_opening_rxns(effluent.mol)
        self.byproduct_formation_rxns(effluent.mol)
    def _cost(self):
        super()._cost()
        self.purchase_costs['RaneyNi catalyst'] =\
            self.mcat_frac * self.ins[0].imass['HMTHP'] * price['RaneyNi']
            
class DehydrationRingOpeningReactor(Reactor):
    """
    A dehydration reactor.
    """
    _N_ins = 2
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'RaneyNi catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('HMTHP -> SA',         'HMTHP',   0.667) # conversion from Chia et al. 2012
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    
    def _run(self):
        feed, recycle = self.ins
        effluent = self.outs[0]
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(effluent.mol)
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['RaneyNi catalyst'] =\
            self.mcat_frac * self.ins[0].imass['HMTHP'] * price['RaneyNi']


class HydrolysisReactor(Reactor):
    """
    A hydrogenation reactor.
    """
    _N_ins = 3
    _N_outs = 1
    _F_BM_default = {**Reactor._F_BM_default,
            'Amberlyst15 catalyst': 1}
    mcat_frac = 0.01 # fraction of catalyst by weight in relation to the reactant (TAL)
    hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('SA + KOH -> KSA + H2O',         'SA',   1.) # not mentioned in Viswanathan et al. 2020
                ])
    TAL_to_SA_rxn = hydrolysis_rxns[0]
    
    def _run(self):
        feed, recycle, reagent = self.ins
        effluent = self.outs[0]
        reagent.imol['KOH'] = max(0, feed.imol['SA'] + recycle.imol['SA'] - recycle.imol['KOH'])
        
        # effluent = feed.copy()
        effluent.mix_from([feed, recycle, reagent])
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.hydrolysis_rxns(effluent.mol)
        
    def _cost(self):
        super()._cost()
        self.purchase_costs['Amberlyst15 catalyst'] =\
            self.mcat_frac * self.ins[0].imass['SA'] * price['Amberlyst15']
            
class Crystallization(Reactor):
    N_ins = 3
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('KSA + HCl -> SA + KCl',         'KSA',   0.98)
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    def __init__(self, ID='', ins=None, outs=(), tau=1., T=293.15, reagent_fraction=1.):
        Reactor.__init__(self, ID, ins, outs, tau=tau)
        self.reagent_fraction = reagent_fraction
    

    def _run(self):
        feed, reagent, recycle_reagent, recycle_feed = ins = self.ins
        effluent, KCl = outs = self.outs
        reagent.empty()
        reagent.imol['HCl'] = max(0, sum([i.imol['KSA'] for i in ins]) - sum([i.imol['HCl'] for i in ins]))
        # effluent = feed.copy()
        effluent.mix_from([feed, 
                           reagent, 
                           recycle_reagent, 
                           recycle_feed])
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(effluent.mol)
        
        KCl.imol['KCl'] = effluent.imol['KCl']
        KCl.T = effluent.T
        KCl.P = effluent.P
        
        # effluent.imol['H2O'] = effluent.imol['H2O']
        
        effluent.imol['HCl'] = 0.
        effluent.imol['KCl'] = 0
        # effluent.imol['H2O'] = 0
        
    def _cost(self):
        super()._cost()


class HClKOHRecovery(Reactor):
    N_ins = 2
    _N_outs = 2
    _F_BM_default = {**Reactor._F_BM_default}
    dehydration_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('KCl + H2O -> HCl + KOH',         'KCl',   1.)
                ])
    TAL_to_SA_rxn = dehydration_rxns[0]
    
    def _run(self):
        feed, reagent = self.ins
        HCl, KOH = self.outs
        
        reagent.imol['H2O'] = feed.imol['KCl']
        # effluent = feed.copy()
        KOH.mix_from([feed, reagent])
        
        # effluent.T = feed.T
        # effluent.P = feed.P
        
        self.dehydration_rxns(KOH.mol)
        
        HCl.imol['HCl'] = KOH.imol['HCl']
        HCl.T = KOH.T 
        HCl.P = KOH.P
        
        KOH.imol['HCl'] = 0
        
    def _cost(self):
        super()._cost()
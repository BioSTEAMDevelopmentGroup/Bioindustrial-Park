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
'''


# %% Setup

import thermosteam as tmo
from flexsolve import aitken_secant
from biosteam import Unit, HeatUtility
from biosteam.units import Flash, HXutility, Mixer, Pump, SolidsSeparator, StorageTank
from biosteam.units.decorators import cost
from biosteam.units._splitter import run_split_with_mixing
from thermosteam import Stream, MultiStream

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction


# %% Feedstock handling

# Feedstock handling system as a whole, capital and operating costs considered in feedstock cost
@cost(basis='Flow rate', ID='System', units='kg/hr',
      kW=511.3205, cost=13329690, S=94697, CE=521.9, n=0.6, BM=1.7)
class FeedstockHandling(Unit): pass


# %% Pretreatment

# Sulfuric acid in-line mixer
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
        mixture.mix_from([water, acid])

# Sulfuric acid addition tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=6210, S=1981, CE=550.8,  n=0.7, BM=3)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      # Size basis changed from 3720 as the original one is not sufficient
      cost=8000, S=1981, CE=521.9, n=0.8, BM=2.3)
class SulfuricAcidAdditionTank(Unit): pass

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
    
    # Energy (kJ) that can be transfered by 1 kmol of heating agent,
    # estimated conservatively based on the energy of lowest energy-carring 
    # biosteam native heating agent (low_pressure_steam), which is 42759 kJ/kmol
    # lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    # lps_heat_duty_over_mol = lps.H * 0.85 (heat transfer efficiency in Humbird et al.)
    heat_duty_over_mol = 40000
    
    def __init__(self, ID='', ins=None, outs=(), *, P):
        #!!! Why here must use super(). as opposed to Unit.?
        super().__init__(ID, ins, outs)
        self.P = P
        
    @staticmethod
    def _P_at_flow(mol_water, P, steam, mixed, feed, heat_duty_over_mol):
        steam.imol['Water'] = mol_water
        mixed.mol = steam.mol + feed.mol
        mixed.H = feed.H + mol_water * heat_duty_over_mol
        P_new = mixed.chemicals.Water.Psat(mixed.T)
        return P - P_new
    
    def _run(self):
        feed, steam = self.ins
        mixed = self.outs[0]
        mixed.P = self.P
        heat_duty_over_mol = self.heat_duty_over_mol
        hu_heating = self.heat_utilities[0]
        # Humbird et al. used high pressure steam, however as only energy balance 
        # is considered here, the type of steam does not affect simulation results,
        # thus all steams used in the system are set to low_pressure_steam
        hu_heating.load_agent(HeatUtility.get_heating_agent('low_pressure_steam'))

        steam_mol = steam.F_mol
        # Results changes a tiny bit each simulation 
        steam_mol = aitken_secant(self._P_at_flow,
                                  steam_mol, steam_mol+0.1, 
                                  1e-4, 1e-4,
                                  args=(self.P, steam, mixed, feed, heat_duty_over_mol))
        hu_heating.duty = steam_mol * heat_duty_over_mol
        hu_heating.flow = steam_mol
        # Regeneration price is accounted for by CAPEX and OPEX of the OrganicAcidsBT unit
        hu_heating.cost = 0
    
    @property
    def installation_cost(self): return 0
    
    @property
    def purchase_cost(self): return 0
    
    def _design(self): pass
    def _cost(self): pass

# Waste vapor condenser
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      cost=34000, S=-2, CE=521.9, n=0.7, BM=2.2)
class WasteVaporCondenser(HXutility):
    _graphics = HXutility._graphics

# Pretreatment reactor
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

# Blowdown tank, costs of Tank and Agitator included in the Pump
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=0, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=0, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=25365, S=292407, CE=550.8, n=0.8, BM=2.3)
class BlowdownTank(Unit): pass

# Oligomer conversion tank, copied PretreatmentFlash as the original design is not sufficient
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
class OligomerConversionTank(Unit): pass

# Pretreatment flash tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
class PretreatmentFlash(Flash): pass

# Ammonia in-line mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      # Size basis on the total flow, not just ammonia, 
      # thus assuming differences caused by MWs of NH3 and NH4OH
      cost=5000, S=157478, CE=521.9, n=0.5, BM=1)
class AmmoniaMixer(Mixer): pass

# Ammonia addition tank, size basis on the total flow, not just ammonia, 
# thus assuming size basis difference caused by MWs of NH3 and NH4OH is negligible,
# pumping is provided by a separate HydrolysatePump unit
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=236000, S=410369, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21900, S=410369, CE=521.9, n=0.5, BM=1.5)
class AmmoniaAdditionTank(Unit): 
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)

        #                                      Reaction definition      Reactant Conversion
        self.neutralization_rxn = Rxn('2 NH4OH + H2SO4 -> NH4SO4 + 2 H2O', 'H2SO4', 0.95)
    
    def _run(self):
        ins = self.ins[0]
        outs = self.outs[0]
        outs.copy_like(ins)        
        self.neutralization_rxn(outs.mol)

# # Pretreatment water heater, not used as mass and energy balance of 
# # all process water is simulated in OrganicAcidsPWC
# @cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
#       cost=92000, S=8, CE=550.8, n=0.7, BM=2.2)
# class PretreatmentWaterHeater(HXutility):
#     _graphics = HXutility._graphics

# Hydrolysate (spelled as "Hydrolyzate" in Humbird et al.) pump
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=74.57, cost=22500, S=402194, CE=521.9, n=0.8, BM=2.3)
class HydrolysatePump(Unit):
    _units = {'Flow rate': 'kg/hr'}
    _graphics = Pump._graphics
    
    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass


# %% Saccharification and fermentation

# Hydrolysate cooler
@cost(basis='Duty', ID='Heat exchanger', units='Gcal/hr',
      cost=85000, S=-8, CE=550.8, n=0.7, BM=2.2)
class HydrolysateCooler(HXutility):
    _graphics = HXutility._graphics

# Enzyme hydrolysate mixer
@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      kW=74.57, cost=109000, S=379938, CE=521.9, n=0.5, BM=1.7)
class EnzymeHydrolysateMixer(Mixer):
    _N_ins = 2
    _N_outs = 1
    
    # Default loading is 20 mg/g cellulose based on Table 18 of Humbird et al.,
    enzyme_loading = 20/1000

    def _run(self):
        feed, cellulase = self.ins
        effluent = self.outs[0]
        
        cellulase.imass['Enzyme'] = self.enzyme_loading*feed.imass['Glucan']
        # 13054 and 13836 are from stream 422 in Humbird et al.
        cellulase.imass['H2O'] = cellulase.imass['Enzyme'] * (13836-13054)/13836
        cellulase.price *= (cellulase.imass['Enzyme']/cellulase.F_mass)
        effluent.mol = feed.mol + cellulase.mol

# Saccharification and co-fermentation (both Glucose & Xylose are used in fermentation)
# Not including heat exchanger as saccharificatoin and co-fermentation 
# are at the same temperature now
@cost(basis='Prehydrolysis tank size', ID='Prehydrolysis tanks', units='kg',
      # Size basis changed from 421776 as the original design is not sufficient
      cost=3840000, S=283875*24, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Saccharification transfer pumps', units='kg/hr',
      kW=74.57, cost=47200, S=421776*24, CE=521.9, n=0.8, BM=2.3)
@cost(basis='Fermenter size', ID='Fermentors', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36), CE=521.9, n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Agitators', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      # and total residence time (batch hydrolysis and fermentation)
      kW=22.371, cost=52500, S=(42607+443391+948+116)*(60+36), CE=521.9, n=1, BM=1.5)
@cost(basis='Flow rate', ID='Recirculation pumps', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in Humbird et al.)
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=521.9, n=0.8, BM=2.3)
class SaccharificationAndCoFermentation(Unit):    
    _N_ins = 3
    _N_outs = 3
    _units= {'Prehydrolysis tank size': 'kg',
             'Flow rate': 'kg/hr',
             'Fermenter size': 'kg'}             

    # Same T for saccharificatoin and co-fermentation
    T = 50+273.15
    
    # Residence time of countinuous saccharification tanks (hr)
    tau_prehydrolysis = 24
    
    # Co-Fermentation time (hr)
    tau_cofermentation = 120
    
    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.1
    
    CSL_loading = 10 # kg/m3
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        # Based on Table 9 on Page 28 of Humbird et al.
        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.9),
            Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
            ])
        
        self.saccharified_stream = Stream(None)
        
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.069),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.069)
        ])
        
        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
            #   Reaction definition                                             Reactant  Conversion
            Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O', 'LacticAcid',   0.95),
            Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O', 'AceticAcid',   0.95)
            ])

    def _run(self):
        feed, CSL, lime = self.ins
        vent, effluent, sidedraw = self.outs
        vent.phase = 'g'
        ss = self.saccharified_stream
        ss.T = sidedraw.T = vent.T = effluent.T = self.T
        
        ss.copy_like(feed)
        CSL.imass['CSL'] = feed.F_vol * self.CSL_loading 
        ss.mol += CSL.mol
        self.saccharification_rxns(ss.mol)
        # Sidedraw to SeedTrain
        sidedraw.mol = ss.mol * self.inoculum_ratio
        effluent.mol = ss.mol - sidedraw.mol
        
        glucose_xylose_0 = effluent.imass['Glucose'] + effluent.imass['Xylose']
        self.cofermentation_rxns(effluent.mol)
        # Assume all lost Glucose and Xylose changed to CO2
        effluent.imass['CO2'] += glucose_xylose_0 \
                                 - (effluent.imass['Glucose']+effluent.imass['Xylose'])
        effluent.imass['CSL'] = 0
        
        # Set feed lime mol to match rate of acids production, add 5% extra
        lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
                            +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1]) \
                            * 1.05
        effluent.mol += lime.mol
        self.neutralization_rxns(effluent.mol)
        vent.receive_vent(effluent)
    
    def _design(self):
        Design = self.design_results
        Design['Prehydrolysis tank size'] = (self.ins[0].F_mass+self.ins[1].F_mass) * self.tau_prehydrolysis
        Design['Flow rate'] = self.outs[1].F_mass
        Design['Fermenter size'] = (self.outs[0].F_mass+self.outs[1].F_mass) \
                                    * self.tau_cofermentation

# Seed train, 5 stages, 2 trains
# assumes no need for pH control due to low acid concentration
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
    _N_outs = 2
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}

    # Operating temperature (K)
    T = 50+273.15
    
    # Cycle time for each batch (hr), including 12 hr turnaround time 
    tau_batch = 36
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)

        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.68),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.069),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.069),
        ])

    def _run(self):
        feed = self.ins[0]
        vent, effluent= self.outs     
        effluent.copy_like(feed)

        self.cofermentation_rxns(effluent.mol)
        # Assume all lost Glucose and Xylose goes to CO2
        effluent.imass['CO2'] += (feed.imass['Glucose']+feed.imass['Xylose']) \
                                 - (effluent.imass['Glucose']+effluent.imass['Xylose'])
        # Assume all CSL is used up
        effluent.imass['CSL'] = 0 
        
        effluent.T = self.T
        vent.phase = 'g'
        vent.receive_vent(effluent)

    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[1].F_mass
        Design['Seed fermenter size'] = self.outs[1].F_mass * self.tau_batch

# Seed hold tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=439000, S=40414, CE=521.9, n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=11.1855, cost=31800, S=40414, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=8200, S=43149, CE=521.9, n=0.8, BM=2.3)
class SeedHoldTank(Unit): pass


# %% Organic acid separation

@cost(basis='Vent flow rate', ID='Scrubber', units='kg/hr',
      cost=215000, S=22608, CE=521.9, n=0.6, BM=2.4)
@cost(basis='Bottom flow rate', ID='Bottom pump', units='kg/hr',
      # Power based on seed hold transfer pump of slightly larger design
      kW=7.457, cost=6300, S=24527, CE=521.9, n=0.8, BM=2.3)
class VentScrubber(Unit): 
    _N_ins = 2
    _N_outs = 2
    _units= {'Vent flow rate': 'kg/hr',
             'Bottom flow rate': 'kg/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), *, gas):
        Unit.__init__(self, ID, ins, outs)
        self.gas = gas
    
    def _run(self):
        water, vent_entry = self.ins
        vent_exit, bottoms = self.outs
        
        vent_exit.copy_like(vent_entry)
        bottoms.copy_flow(vent_exit, self.gas, remove=True, exclude=True)
        bottoms.mol += water.mol
        
    def _design(self):
        Design = self.design_results
        Design['Vent flow rate'] = self.ins[1].F_mass
        Design['Bottom flow rate'] = self.outs[1].F_mass

# Scaling values in Humbird et al. was adjusted, as the original design used
# sludge flow rate for scaling, which was much smaller than the feed flow rate
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
    
    def _run(self):      
        run_split_with_mixing(self)
        retentate, permeate = self.outs
        solids = retentate.F_mass
        mc = self.mositure_content
        retentate.imol['Water'] = water = (solids * mc/(1-mc))/18.01528
        permeate.imol['Water'] -= water
        if permeate.imol['Water'] < water:
            raise ValueError(f'not enough water for {repr(self)}')
            
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

# Cost copied from PretreatmentFlash
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
class AcidulationReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
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

# Cost copied from PretreatmentFlash
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
class EsterificationReactor(Unit):
    _N_ins = 3
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        self.esterification_rxns = ParallelRxn([
            #   Reaction definition                               Reactant  Conversion
            Rxn('LacticAcid + Methanol -> MethylLactate + H2O', 'LacticAcid', 0.95),
            Rxn('AceticAcid + Methanol -> MethylAcetate + H2O', 'AceticAcid', 0.95)
            ])

    def _run(self):
        feed, recycled_methanol, supplement_methanol = self.ins
        effluent = self.outs[0]
        
        # Add 5% extra
        methanol_needed = (feed.imol['LacticAcid']/self.esterification_rxns.X[0] \
                           +feed.imol['AceticAcid']/self.esterification_rxns.X[1]) \
                          * 1.05 - (feed.imol['Methanol']+recycled_methanol.imol['Methanol'])
        supplement_methanol.imol['Methanol'] = max(0, methanol_needed)
        effluent.mix_from([feed, recycled_methanol, supplement_methanol])
        self.esterification_rxns(effluent.mol)

# Cost copied from PretreatmentFlash
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=511000, S=264116, CE=521.9, n=0.7, BM=2)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=170, cost=90000, S=252891, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=55.9275, cost=30000, S=204390, CE=521.9, n=0.8, BM=2.3)
class HydrolysisReactor(Unit):
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=(), outs=()):
        Unit.__init__(self, ID, ins, outs)
        
        self.hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                 Reactant   Conversion
            Rxn('MethylLactate + H2O -> LacticAcid + Methanol', 'MethylLactate', 0.95),
            Rxn('MethylAcetate + H2O -> AceticAcid + Methanol', 'MethylAcetate', 0.95)
            ])
    
    def _run(self):
        feed, water = self.ins
        feed.phase = 'l'
        effluent = self.outs[0]
        
        # Add 5% extra
        water_needed = (feed.imol['MethylLactate']/self.hydrolysis_rxns.X[0] \
                        +feed.imol['MethylAcetate']/self.hydrolysis_rxns.X[1]) \
                       * 1.05
        water.imol['Water'] = max(0, (water_needed-feed.imol['Water'])) 
        effluent.mix_from([feed, water])
        
        self.hydrolysis_rxns(effluent.mol)


# %% Wastewater treatment

# Total cost of wastewater treatment is combined into this placeholder
@cost(basis='Flow rate', ID='Wastewater system', units='kg/hr', 
      kW=7018.90125, S=393100, cost=50280080, CE=550.8, n=0.6, BM=1)
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
        self.multi_stream = MultiStream(None)
    
    def _run(self):
        feed, cool_water = self.ins
        biogas, waste, sludge, hot_water = self.outs
        biogas.phase = 'g'
        hot_water.link_with(cool_water, TP=False)
        biogas.T = waste.T = sludge.T = hot_water.T = 35+273.15
        H_at_35C = feed.thermo.mixture.H(z=feed.mol, phase='l', T=35+273.15, P=101325)
        # Water flow is adjusted to maintain heat balance
        cool_water.mol *= (feed.H - H_at_35C)/(hot_water.H - cool_water.H)
        sludge.copy_flow(feed)
        self.digestion_rxns(sludge.mol)
        self.multi_stream.copy_flow(sludge)
        self.multi_stream.vle(P=101325, H=self.multi_stream.H)
        biogas.mol = self.multi_stream.imol['g']
        liquid_mol = self.multi_stream.imol['l']
        sludge.mol = liquid_mol * self.sludge_split
        waste.mol = liquid_mol - sludge.mol
        biogas.receive_vent(waste, accumulate=True)     
    
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
    # 4350, 4379, 356069, 2252, 2151522, and 109089 are water flows from 
    # streams 622, 630, 611, 632, 621, and 616  in Humbird et al.
    evaporation = 4350/(4379+356069+2252+2151522+109089)
    
    def __init__(self, ID='', ins=None, outs=(), *, digestion_rxns, ratio):
        Unit.__init__(self, ID, ins, outs)
        self.digestion_rxns = digestion_rxns
        self.ratio = ratio
        
        #                                      Reaction definition       Reactant Conversion
        self.neutralization_rxn = Rxn('H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O', 'H2SO4', 0.95)
    
    def _run(self):
        influent, air, caustic = self.ins
        vent, effluent = self.outs
        ratio = self.ratio
        vent.phase = 'g'

        # 51061 and 168162 from stream 630 in Humbird et al.
        air.imol['O2'] = 51061 * ratio
        air.imol['N2'] = 168162 * ratio
        # 2252 from stream 632 in Humbird et al
        caustic.imol['NaOH'] = (2252*ratio) + (2*influent.imol['H2SO4']/self.neutralization_rxn.X)
        caustic.imol['H2O'] = caustic.imol['NaOH']
        effluent.copy_like(influent)
        effluent.mol += air.mol
        effluent.mol += caustic.mol
        self.neutralization_rxn(effluent.mol)
        self.digestion_rxns(effluent.mol)
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'))
        vent.imol['Water'] = effluent.imol['Water'] * self.evaporation
        effluent.imol['Water'] -= vent.imol['Water']


# %% Facility units in biosteam

# Sulfuric acid storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=550.8, n=0.7, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      # Size basis changed from 1981 as the original design is not sufficient
      kW=0.37285, cost=7493, S=1136, CE=550.8, n=0.8, BM=2.3)
class SulfuricAcidStorageTank(StorageTank): pass

# Ammonia storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      # Original size basis for NH3 instead of NH4OH
      cost=196000, S=1171/17.031*35.046, CE=550.8, n=0.7, BM=2)
# Design in Humbird et al. has no pump, this one copied from SulfuricAcidStorageTank
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      # Size basis changed from 1981 as the original design is not sufficient
      kW=0.37285, cost=7493, S=1136, CE=550.8, n=0.8, BM=2.3)
class AmmoniaStorageTank(StorageTank): pass

# CSL storage tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=70000, S=1393, CE=521.9, n=0.7, BM=2.6)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21200, S=1393, CE=521.9, n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=1393, CE=521.9, n=0.8, BM=3.1)
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
class LimeStorageTank(Unit): pass

# Fire water tank
@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=803000, S=8343, CE=521.9, n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=15000, S=8343, CE=521.9, n=0.8, BM=3.1)
class FireWaterTank(Unit): pass


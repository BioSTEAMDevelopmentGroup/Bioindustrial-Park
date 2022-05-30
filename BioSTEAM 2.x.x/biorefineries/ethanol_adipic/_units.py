#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234
[3] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbons: Dilute-Acid and Enzymatic Deconstruction of Biomass
    to Sugars and Catalytic Conversion of Sugars to Hydrocarbons;
    NREL/TP-5100-62498; National Renewable Energy Lab (NREL), 2015.
    http://www.nrel.gov/docs/fy15osti/62498.pdf

'''


# %% Setup

import thermosteam as tmo
from biosteam import Unit
from biosteam.units import HXutility, Pump
from biosteam.units.decorators import cost
from biorefineries.ethanol_adipic._settings import price
from biorefineries.ethanol_adipic._chemicals import total_solids, solubles, insolubles
from biorefineries.ethanol_adipic._utils import CEPCI, compute_muconic_titer

# These units are the same as those from the lactic acid biorefinery
from biorefineries.lactic._units import (
    SulfuricAcidAdditionTank,
    SulfuricAcidMixer,
    PretreatmentMixer,
    AcidPretreatment,
    BlowdownTank,
    OligomerConversionTank,
    PretreatmentFlash,
    AmmoniaMixer,
    AmmoniaAdditionTank,
    WasteVaporCondenser,
    HydrolysatePump,
    EnzymeHydrolysateMixer,
    CellMassFilter,
    AnaerobicDigestion,
    AerobicDigestion,
    MembraneBioreactor,
    BeltThickener,
    SludgeCentrifuge,
    ReverseOsmosis,
    SulfuricAcidStorage,
    AmmoniaStorage,
    CausticStorage,
    CSLstorage,
    FirewaterStorage,
    )

_GPM_2_m3hr = (3.78541*60) / 1e3
_Gcal_2_kJ = 4.184 * 1e6 # (also MMkcal/hr)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction


# %%

# =============================================================================
# Pretreatment (base-specific ones)
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
        self.heat_utilities[0](unit_duty=duty, T_in=self.ins[0].T)
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
# Conversion
# =============================================================================

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
      cost=85000, S=-8*_Gcal_2_kJ, CE=CEPCI[2010], n=0.7, BM=2.2)
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 4
    _N_outs = 3
    _N_heat_utilities = 2

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
        hu_hydrolysate, hu_saccharification = self.heat_utilities
        mixture = self.thermo.mixture
        ss = self.saccharified_stream

        ss_in = ss.copy()
        ss_in.mix_from(self.ins)
        # To cool or heat hydrolysate to the fermentation temperature
        hu_hydrolysate(unit_duty=ss.H-ss_in.H, T_in=ss_in.T)

        mol = ss.mol
        duty = (mixture.H('l', mol, self.T_fermentation, 101325.)
                - mixture.H('l', mol, self.T_saccharification, 101325.))
        self.design_results['Duty'] = duty
        # To cool the hydrolysate down to fermentation temperature
        hu_saccharification(duty, self.T_fermentation)


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
        self.heat_utilities[0](unit_duty=duty, T_in=self.T_in)
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

    target_titer = 100
    # target_titer = (34.5+68.5) / 2 # in g/L (kg/m3)

    effluent_titer = 0

    def __init__(self, ID='', ins=None, outs=(), T=32+273.15, P=1.34*101325):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
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

    # Lower conversions of extratives and lignin than ref [3] to better reflect
    # current state of technology, resulting in a titer of ~31.5 g/L for the
    # composition in ref [3]. Titers from real hydrolysate remains low (< 15 g/L)
    # based on refs [5] and [6]
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
    Rxn('Extractives + 0.68 O2 + 0.28 NH4OH -> 4.8 P_putida + 1.2 CO2 + 2.26 H2O',  'Extractives',    0),
    Rxn('Extractives + 1.44 O2 -> 0.74 MuconicAcid + 1.57 CO2 + 3.78 H2O',          'Extractives',    0),
    Rxn('SolubleLignin + 3 O2 -> MuconicAcid + 2 CO2 + H2O',                        'SolubleLignin',  0.5)
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
        self.heat_utilities[0](unit_duty=duty, T_in=mixture.T)
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
        self.heat_utilities[0](unit_duty=duty, T_in=self.T_in)
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
        self.heat_utilities[0](unit_duty=duty, T_in=adipic.T)
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
        self._F_BM_default['Ru/C catalyst'] = 1
        # WHSV (feed mass flow/catalyst mass=5/hr)
        self.baseline_purchase_costs['Ru/C catalyst'] = self.ins[0].F_mass/5 \
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
        self.heat_utilities[0](unit_duty=duty, T_in=adipic_total.T)
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
        self.heat_utilities[0](unit_duty=duty, T_in=self.ins[0].T)

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

    def __init__(self, ID='', ins=None, outs=()):
        Unit.__init__(self, ID, ins, outs)
        self.decomposition_rxn = Rxn('NH4OH -> NH3 + H2O', reactant='NH4OH', X=1)

    def _run(self):
        brine = self.ins[0]
        vent, residuals, Na2SO4 = self.outs

        influent = brine.copy()
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
        self.heat_utilities[0](unit_duty=self.design_results['Duty'], T_in=self.ins[0].T)
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

@cost(basis='Flow rate', ID='Unloader', units='kg/hr',
      cost=30000, S=163, CE=CEPCI[2009], n=0.6, BM=1.7)
@cost(basis='Flow rate', ID='Make-up tank', units='kg/hr',
      cost=81192, S=142, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=4.10135, cost=9800, S=163, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=163, CE=CEPCI[2009], n=0.8, BM=3.1)
class DAPstorage(Unit): pass
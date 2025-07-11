#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
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
[3] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.
[4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbons: Dilute-Acid and Enzymatic Deconstruction of Biomass
    to Sugars and Catalytic Conversion of Sugars to Hydrocarbons;
    NREL/TP-5100-62498; National Renewable Energy Lab (NREL), 2015.
    http://www.nrel.gov/docs/fy15osti/62498.pdf

'''


# %% Setup

import numpy as np
import thermosteam as tmo
from math import exp, pi, ceil
from biosteam import Stream, Unit, main_flowsheet
from biosteam.exceptions import DesignError
from biosteam.units import HXutility, Mixer, SolidsSeparator, StorageTank
from biosteam.units.design_tools import PressureVessel, pressure_vessel_material_factors as factors
from biosteam.units.decorators import cost
from . import (
    CEPCI, get_baseline_feedflow, compute_lactic_titer,
    compute_extra_chemical, adjust_recycle, compute_COD,
    )

_MGD_2_m3hr = 157.7255 # auom('gallon').conversion_factor('m3')*1e6/24
_GPM_2_m3hr = 0.2271 # auom('gallon').conversion_factor('m3')*60
_Gcal_2_kJ = 4184000 # auom('kcal').conversion_factor('kJ')*1e6 , also MMkcal/hr
_316_over_304 = factors['Stainless steel 316'] / factors['Stainless steel 304']
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# Feedstock flow rate in dry U.S. ton per day, 907.1847 is auom('ton').conversion_factor('kg')
def get_flow_tpd(flowsheet=None):
    flowsheet = flowsheet or main_flowsheet
    U101 = flowsheet.unit.U101
    feedstock = U101.ins[0]
    return (feedstock.F_mass-feedstock.imass['H2O'])*24/907.1847*(1-U101.divert_ratio)


# %%

# =============================================================================
# Feedstock preprocessing
# =============================================================================

# Capital and operating costs already considered in
# the cost of feedstock cost
class FeedstockPreprocessing(Unit):
    _N_outs = 2
    # 2205 U.S. ton/day (2000 metric tonne/day) as in ref [1]
    _cached_flow_rate = 2205

    def __init__(self, ID='', ins=None, outs=(), thermo=None, divert_ratio=0):
        Unit.__init__(self, ID, ins, outs, thermo)
        self._baseline_flow_rate = get_baseline_feedflow(self.chemicals).sum()
        self.divert_ratio = divert_ratio

    def _run(self):
        total = self.ins[0]
        processed, burned = self.outs
        burned.mass = self.divert_ratio * total.mass
        processed.mass = total.mass - burned.mass

    def _design(self):
        self.design_results['Flow rate'] = self.outs[0].F_mass


# %%

# =============================================================================
# Conversion
# =============================================================================

@cost(basis='Flow rate', ID='Mixer', units='kg/hr',
      kW=74.57, cost=109000, S=379938, CE=CEPCI[2009], n=0.5, BM=1.7)
class EnzymeHydrolysateMixer(Mixer):
    _N_ins = 3
    _N_outs = 1
    _graphics = Mixer._graphics
    auxiliary_unit_names = ('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(),
                 enzyme_loading=20, solids_loading=0.2, T=None):
        Unit.__init__(self, ID, ins, outs)
        self.enzyme_loading = enzyme_loading
        self.solids_loading = solids_loading
        self.T = T
        ID = self.ID
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)


    def _run(self):
        hydrolysate, enzyme, water = self.ins
        effluent = self.outs[0]

        # 10% extra based on Page 23 of ref [2]
        enzyme.imass['Enzyme'] = (self.enzyme_loading/1000*1.1) * hydrolysate.imass['Glucan']
        mixture = hydrolysate.copy()
        mixture.mix_from([hydrolysate, enzyme])

        total_mass = (mixture.F_mass-mixture.imass['Water'])/self.solids_loading
        water.imass['Water'] = max(0, total_mass-mixture.F_mass)

        effluent.mix_from([hydrolysate, enzyme, water])
        if self.T: effluent.T = self.T


    def _design(self):
        self.design_results['Flow rate'] = self.F_mass_in
        hx = self.heat_exchanger
        hx.ins[0].mix_from(self.ins)
        hx.outs[0].copy_like(self.outs[0])
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)


# Saccharification and co-fermentation (both glucose & xylose are used in fermentation)
# Not including heat exchanger as saccharification and co-fermentation
# are at the same temperature now
@cost(basis='Saccharification tank size', ID='Saccharification tank', units='kg',
      cost=3840000, S=421776*24, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Slurry flow rate', ID='Saccharification transfer pump', units='kg/hr',
      kW=74.57, cost=47200, S=421776, CE=CEPCI[2009], n=0.8, BM=2.3)
@cost(basis='Fermenter size', ID='Fermenter', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Fermenter agitator', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [1])
      # and total residence time (batch hydrolysis and fermentation)
      kW=268.452, cost=630000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Recirculation flow rate', ID='Recirculation pump', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [1])
      # Pumps already in SS316
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)
# Surge tank to hold the fermentation broth
@cost(basis='Broth flow rate', ID='Surge tank', units='kg/hr',
      cost=636000, S=425878, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Broth flow rate', ID='Surge tank agitator', units='kg/hr',
      kW=29.828, cost=68300, S=425878, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Broth flow rate', ID='Surge tank pump', units='kg/hr',
      kW=93.2125, cost=26800, S=488719, CE=CEPCI[2009], n=0.8, BM=2.3)
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 5
    _N_outs = 2
    _units= {
        'Saccharification tank size': 'kg',
        'Slurry flow rate': 'kg/hr',
        'Fermenter size': 'kg',
        'Recirculation flow rate': 'kg/hr',
        'Broth flow rate': 'kg/hr',
        }

    auxiliary_unit_names = ('heat_exchanger',)

    tau_saccharification = 24 # in hr

    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.07

    CSL_loading = 10 # g/L (kg/m3)

    target_titer = 130 # in g/L (kg/m3), the maximum titer in collected data

    effluent_titer = 0

    productivity = 0.89 # in g/L/hr

    target_yield = 0.76 # in g/g-sugar

    tau_cofermentation = 0 # will be calculated based on titer and productivity

    tau_turnaround = 12 # in hr, the same as the seed train in ref [1]

    def __init__(self, ID='', ins=None, outs=(), T=50+273.15,
                 neutralization=True, allow_dilution=False):
        Unit.__init__(self, ID, ins, outs)
        # Same T for saccharification and co-fermentation
        self.T = T
        ID = self.ID
        self.saccharified_stream = Stream(f'{ID}_ss')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)

        self.neutralization = neutralization
        self.allow_dilution = allow_dilution

        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.9),
            Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
            ])

        # FermMicrobe reaction from ref [1]
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.02),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.07),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.02)
        ])
        self._X = self.cofermentation_rxns.X.copy()

        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'LacticAcid',   1),
        Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1),
        Rxn('SuccinicAcid + CalciumDihydroxide -> CalciumSuccinate + 2H2O', 'SuccinicAcid', 1)
            ])

    def _run(self):
        feed, inoculum, CSL, lime, water = self.ins
        effluent, sidedraw = self.outs
        ss = self.saccharified_stream

        CSL.imass['CSL'] = feed.F_vol * self.CSL_loading
        if not self.allow_dilution:
            water.empty()
        ss.mix_from((feed, inoculum, CSL, water))
        ss.T = sidedraw.T = effluent.T = self.T

        self.saccharification_rxns(ss.mol)
        # Sidedraw to SeedTrain
        sidedraw.mol = ss.mol * self.inoculum_ratio
        effluent.mol = ss.mol - sidedraw.mol
        # effluent_copy = effluent.copy()

        self.cofermentation_rxns(effluent.mol)
        # Assume all CSL is consumed
        effluent.imass['CSL'] = 0

        # Need lime to neutralize produced acid
        if self.neutralization:
            # Set feed lime mol to match rate of acids production, add 10% extra
            lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
                                +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1] \
                                +effluent.imol['SuccinicAcid']/self.neutralization_rxns.X[2]) \
                                * 1.1
            effluent.mix_from((effluent, lime))
            self.neutralization_rxns.adiabatic_reaction(effluent)
        else:
            lime.empty()
        self.effluent_titer = compute_lactic_titer(effluent)

    def _design(self):
        Design = self.design_results
        total_mass_flow = self.ins[0].F_mass + self.ins[1].F_mass
        Design['Saccharification tank size'] = \
            total_mass_flow * self.tau_saccharification
        Design['Slurry flow rate'] = total_mass_flow
        tau = self.tau_cofermentation = self.effluent_titer/self.productivity
        Design['Fermenter size'] = self.outs[0].F_mass * (tau+self.tau_turnaround)
        Design['Recirculation flow rate'] = total_mass_flow # internal circulation
        Design['Broth flow rate'] = self.outs[0].F_mass

        hx = self.heat_exchanger
        hx.ins[0].mix_from(self.ins[0:3])
        hx.outs[0].copy_like(hx.ins[0])
        hx.outs[0].T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)


    def _cost(self):
        super()._cost()
        self._decorated_cost()
        # Adjust fermenter cost for acid-resistant scenario
        if not self.neutralization:
            self.baseline_purchase_costs['Fermenter'] *= _316_over_304
            self.baseline_purchase_costs['Agitator'] *= _316_over_304

    @property
    def lactic_yield(self):
        X = self.cofermentation_rxns.X
        if not X[0] == X[3]:
            raise ValueError('Glucose and xylose has different yields.')
        return X[0]


# Saccharification and co-fermentation (both glucose & xylose are used in fermentation)
# Not including heat exchanger as saccharification and co-fermentation
# are at the same temperature now
@cost(basis='Saccharification tank size', ID='Saccharification tank', units='kg',
      cost=3840000, S=421776*24, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Slurry flow rate', ID='Saccharification transfer pump', units='kg/hr',
      kW=74.57, cost=47200, S=421776, CE=CEPCI[2009], n=0.8, BM=2.3)
class Saccharification(Unit):
    _N_ins = 1
    _N_outs = 1
    _units= {'Saccharification tank size': 'kg',
             'Slurry flow rate': 'kg/hr'}

    # Extend saccharification time based on ref [1] as saccharification and
    # co-fermentation are separated now
    tau_saccharification = 84 # in hr

    def __init__(self, ID='', ins=None, outs=(), T=50+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.T = T

        self.saccharification_rxns = ParallelRxn([
            #   Reaction definition                   Reactant        Conversion
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',         0.04),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',         0.012),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',         0.9),
            Rxn('Cellobiose + H2O -> 2 Glucose',      'Cellobiose',     1)
            ])

    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(feed)
        self.saccharification_rxns(effluent.mol)

    def _design(self):
        Design = self.design_results
        Design['Saccharification tank size'] = self.F_mass_in * self.tau_saccharification
        Design['Slurry flow rate'] = self.F_mass_in


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
    tau : float
        Residence time [hr].
    V_wf : float
        Fraction of working volume over total volume.
    mixing_intensity: float
        Mechanical mixing intensity, [/s].
    kW_per_m3: float
        Power usage of agitator
        (converted from 0.5 hp/1000 gal as in [1]).
        If mixing_intensity is provided, this will be calculated based on
        the mixing_intensity and viscosity of the influent mixture as in [2]_
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
    .. [2] Shoener et al. Energy Positive Domestic Wastewater Treatment:
        The Roles of Anaerobic and Phototrophic Technologies.
        Environ. Sci.: Processes Impacts 2014, 16 (6), 1204–1222.
        `<https://doi.org/10.1039/C3EM00711A>`_.

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
    _Vmax = pi/4*(20**2)*40/35.3147

    def __init__(self, ID='', ins=None, outs=(), *,
                 P=101325, tau=0.5, V_wf=0.8,
                 length_to_diameter=2, mixing_intensity=None, kW_per_m3=0.0985,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):

        Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
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

        N = ceil(V_total/self._Vmax)
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
             Design['Weight'] = round(Design['Weight']*wall_thickness_factor, 2)

    def _cost(self):
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs

        if Design['Total volume'] == 0:
            for i, j in purchase_costs.items():
                purchase_costs[i] = 0
        else:
            purchase_costs.update(self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']))
            for i, j in purchase_costs.items():
                purchase_costs[i] *= Design['Number of reactors']

            self.power_utility(self.kW_per_m3*Design['Total volume'])


    @property
    def kW_per_m3(self):
        G = self.mixing_intensity
        if G is None:
            return self._kW_per_m3
        else:
            if not hasattr(self, '_mixed_in'):
                mixed_in = self._mixed_in = Stream(f'{self.ID}_mixed_in')
            else: mixed_in = self._mixed_in
            mixed_in.mix_from(self.ins)
            kW_per_m3 = mixed_in.mu*(G**2)/1e3
            return kW_per_m3
    @kW_per_m3.setter
    def kW_per_m3(self, i):
        if self.mixing_intensity and i is not None:
            raise AttributeError('`mixing_intensity` is provided, kw_per_m3 will be calculated.')
        else:
            self._kW_per_m3 = i


@cost(basis='Fermenter size', ID='Fermenter', units='kg',
      cost=10128000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Fermenter size', ID='Fermenter agitator', units='kg',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [1])
      # and total residence time (batch hydrolysis and fermentation)
      kW=268.452, cost=630000, S=(42607+443391+948+116)*(60+36),
      CE=CEPCI[2009], n=1, BM=1.5)
@cost(basis='Recirculation flow rate', ID='Recirculation pump', units='kg/hr',
      # Scaling basis based on sum of all streams into fermenter
      # (304, 306, 311, and 312 in ref [1])
      kW=74.57, cost=47200, S=(42607+443391+948+116), CE=CEPCI[2009], n=0.8, BM=2.3)
# Surge tank to hold the fermentation broth
@cost(basis='Broth flow rate', ID='Surge tank', units='kg/hr',
      cost=636000, S=425878, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Broth flow rate', ID='Surge tank agitator', units='kg/hr',
      kW=29.828, cost=68300, S=425878, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Broth flow rate', ID='Surge tank pump', units='kg/hr',
      kW=93.2125, cost=26800, S=488719, CE=CEPCI[2009], n=0.8, BM=2.3)
class CoFermentation(Reactor):
    _N_ins = 6
    _N_outs = 2
    _units= {
        **Reactor._units,
        'Fermenter size': 'kg',
        'Recirculation flow rate': 'kg/hr',
        'Broth flow rate': 'kg/hr',
        }

    auxiliary_unit_names = ('heat_exchanger',)

    # Equals the split of saccharified slurry to seed train
    inoculum_ratio = 0.07

    CSL_loading = 10 # g/L (kg/m3)

    target_titer = 130 # in g/L (kg/m3), the maximum titer in collected data

    effluent_titer = 0

    productivity = 0.89 # in g/L/hr

    target_yield = 0.76 # in g/g-sugar

    tau_batch_turnaround = 12 # in hr, the same as the seed train in ref [1]

    tau_cofermentation = 0 # in hr

    max_sugar = 0 # g/L

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, T=50+273.15,
                 P=101325, V_wf=0.8, length_to_diameter=0.6,
                 kW_per_m3=None, mixing_intensity=300,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 neutralization=True,
                 mode='batch', feed_freq=1,
                 allow_dilution=False,
                 allow_concentration=False,
                 sugars=None):

        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.neutralization = neutralization
        self.mode = mode
        self.feed_freq = feed_freq
        self.allow_dilution = allow_dilution
        self.allow_concentration = allow_concentration
        self.sugars = sugars or tuple(i.ID for i in self.chemicals.sugars)
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

        # FermMicrobe reaction from ref [1]
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.02),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.07),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.02)
        ])
        self._X = self.cofermentation_rxns.X.copy()

        # Neutralization of lactic acid and acetic acid by lime (Ca(OH)2)
        self.neutralization_rxns = ParallelRxn([
        #   Reaction definition                                               Reactant  Conversion
        Rxn('2 LacticAcid + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'LacticAcid',   1),
        Rxn('2 AceticAcid + CalciumDihydroxide -> CalciumAcetate + 2 H2O',  'AceticAcid',   1),
        Rxn('SuccinicAcid + CalciumDihydroxide -> CalciumSuccinate + 2H2O', 'SuccinicAcid', 1)
            ])

    def _run(self):
        # add_feed for fed-batch
        init_feed, inoculum, CSL, lime, water, add_feed = self.ins
        effluent, sidedraw = self.outs
        mixed_feed = self._mixed_feed
        ferm_rxns = self.cofermentation_rxns
        feed_freq = self.feed_freq

        tot_feed = self._tot_feed
        tot_feed.mix_from((init_feed, add_feed))
        inoculum_r = self.inoculum_ratio - 1e-6 # this is to prevent tiny rounding errors
        CSL.imass['CSL'] = tot_feed.F_vol * self.CSL_loading
        if not self.allow_dilution:
            water.empty()

        mixed_feed.mix_from((tot_feed, inoculum, CSL, water))
        # For seed preparation
        sidedraw.mol = mixed_feed.mol * inoculum_r

        if not 1/feed_freq >= inoculum_r:
            raise ValueError('Not enough initial feed for seed inoculum.')

        effluent.copy_like(mixed_feed)
        effluent.mol = mixed_feed.mol - sidedraw.mol
        self.influent_titer = titer0 = compute_lactic_titer(effluent)
        sugars = self.sugars
        if feed_freq == 1:
            self.max_sugar = effluent.imass[sugars].sum()/effluent.F_vol
            ferm_rxns(effluent.mol)
            self._init.empty()
        else:
            init = self._init
            init.mix_from((*self.ins[:3], self.ins[4]))
            init.mol -= sidedraw.mol
            ferm_rxns(init.mol)

            single_feed0 = self._single_feed0
            single_feed1 = self._single_feed1
            last = self._last
            single_feed0.copy_like(add_feed)
            if single_feed0.F_mass != 0:
                single_feed0.F_mass /= (feed_freq-1)
            single_feed1.copy_like(single_feed0)

            ferm_rxns(single_feed1.mol)
            last.mix_from((init, *[single_feed1.copy()]*max(0, (feed_freq-2)),
                           single_feed0))
            self.max_sugar = last.imass[sugars].sum()/last.F_vol
            effluent.mix_from((init, *[single_feed1.copy()]*(feed_freq-1)))

        # Assume all CSL is consumed
        effluent.imass['CSL'] = 0

        # Need lime to neutralize produced acid
        if self.neutralization:
            self.vessel_material= 'Stainless steel 304'
            # Set feed lime mol to match rate of acids production, add 10% extra
            lime.imol['Lime'] = (effluent.imol['LacticAcid']/2/self.neutralization_rxns.X[0] \
                                +effluent.imol['AceticAcid']/2/self.neutralization_rxns.X[1] \
                                +effluent.imol['SuccinicAcid']/self.neutralization_rxns.X[2]) \
                                * 1.1
            effluent.mix_from((effluent, lime))
            self.neutralization_rxns.adiabatic_reaction(effluent)
        else:
            self.vessel_material= 'Stainless steel 316'
            lime.empty()

        self.effluent_titer = titer1 = compute_lactic_titer(effluent)
        self.tau_cofermentation = (titer1-titer0)/self.productivity

        mixed_feed.T = sidedraw.T = effluent.T = self.T


    def _design(self):
        mode = self.mode
        Design = self.design_results
        Design.clear()

        if mode == 'batch':
            tau_tot = self.tau_batch_turnaround + self.tau_cofermentation
            Design['Fermenter size'] = self.outs[0].F_mass * tau_tot
            Design['Recirculation flow rate'] = self.F_mass_in
            Design['Broth flow rate'] = self.outs[0].F_mass
        else:
            self._Vmax = 3785.41178 # 1,000,000 gallon from ref [1]
            Reactor._design(self)
            self.tau_cofermentation = self.tau
            # Include a backup fermenter for cleaning
            Design['Number of reactors'] += 1
            
        hx = self.heat_exchanger
        hx.ins[0].mix_from((*self.ins[0:3], *self.ins[4:]))
        hx.outs[0].copy_like(hx.ins[0])
        hx.outs[0].T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)


    def _cost(self):
        Design = self.design_results
        purchase_costs = self.baseline_purchase_costs
        purchase_costs.clear()
        
        if self.mode == 'batch':
            Unit._cost()
            self._decorated_cost()
            # Adjust fermenter cost for acid-resistant scenario
            if not self.neutralization:
                purchase_costs['Fermenter'] *= _316_over_304
                purchase_costs['Agitator'] *= _316_over_304
        else:
            if not self.neutralization:
                self.vessel_material= 'Stainless steel 316'
            Reactor._cost(self)
            N_working = Design['Number of reactors'] - 1 # subtract the one back-up
            # No power need for the back-up reactor
            self.power_utility(self.kW_per_m3*Design['Single reactor volume']*N_working)

    @property
    def tau(self):
        return self.tau_cofermentation

    @property
    def mode(self):
        return self._mode
    @mode.setter
    def mode(self, i):
        if i.lower() in ('fed-batch', 'fedbatch', 'fed batch'):
            raise ValueError('For fed-batch, set fermentation to "batch" and ' \
                             'change feed_freq.')
        elif i.lower() in ('batch', 'continuous'):
            self._mode = i.lower()
            self.feed_freq = 1
        else:
            raise ValueError(f'Mode can only be "batch" or "continuous", not {i}.')

    @property
    def feed_freq(self):
        return self._feed_freq
    @feed_freq.setter
    def feed_freq(self, i):
        if not (int(i)==i and i >0):
            raise ValueError('feed_freq can only be positive integers.')
        elif self.mode == 'continuous' and int(i) != 1:
            raise ValueError('feed_freq can only be 1 for continuous mode.')
        else:
            self._feed_freq = int(i)

    @property
    def lactic_yield(self):
        X = self.cofermentation_rxns.X
        if not X[0] == X[3]:
            raise ValueError('Glucose and xylose has different yields.')
        return X[0]


# Seed train, 5 stages, 2 trains
@cost(basis='Seed fermenter size', ID='Stage #1 fermenter', units='kg',
      # 44339, 211, and 26 are streams 303, 309, and 310 in ref [1]
      cost=75400*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #2 fermenter', units='kg',
      cost=116600*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #3 fermenter', units='kg',
      cost=157600*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Seed fermenter size', ID='Stage #4 fermenter', units='kg',
      cost=352000*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #4 agitator', units='kg',
      kW=11.1855*_316_over_304, cost=26000, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Seed fermenter size', ID='Stage #5 fermenter', units='kg',
      cost=1180000*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.7, BM=2)
@cost(basis='Seed fermenter size', ID='Stage #5 agitator', units='kg',
      kW=14.914, cost=43000*_316_over_304, S=(44339+211+26)*36, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Seed train pump', units='kg/hr',
      kW=59.656, cost=24300, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedTrain(Unit):
    _N_ins = 1
    _N_outs = 1
    _units= {'Seed fermenter size': 'kg',
             'Flow rate': 'kg/hr'}

    tau_turnaround = 12 # in hr, the same as ref [1]

    effluent_titer = 0

    productivity = 0.89*0.95 # in g/L/hr

    # Yield as a ratio of the yield in the main fermenter
    ferm_ratio = 0.95

    def __init__(self, ID='', ins=None, outs=(), T=50+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.T = T

        # FermMicrobe reaction from ref [1]
        self.cofermentation_rxns = ParallelRxn([
        #      Reaction definition            Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',        'Glucose',   0.76*0.95),
        Rxn('Glucose -> 3 AceticAcid',        'Glucose',   0.07*0.95),
        Rxn('Glucose -> 6 FermMicrobe',       'Glucose',   0.04),
        Rxn('3 Xylose -> 5 LacticAcid',       'Xylose',    0.76*0.95),
        Rxn('2 Xylose -> 5 AceticAcid',       'Xylose',    0.07*0.95),
        Rxn('Xylose -> 5 FermMicrobe',        'Xylose',    0.04)
        ])
        self._X = self.cofermentation_rxns.X.copy()

    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]
        effluent.copy_like(feed)

        self.cofermentation_rxns(effluent.mol)
        self.effluent_titer = compute_lactic_titer(effluent)

        # Assume all CSL is consumed
        effluent.imass['CSL'] = 0
        effluent.T = self.T

    def _design(self):
        Design = self.design_results
        Design['Flow rate'] = self.outs[0].F_mass
        tau_total = self.tau_turnaround + self.effluent_titer/self.productivity
        Design['Seed fermenter size'] = self.outs[0].F_mass * tau_total


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
        # 809 is the scaling basis of equipment M-505,
        # 391501 from stream 508 in ref [1]
        Design['Pressing air flow rate'] = 809/391501 * self.ins[0].F_mass
        # 12105 and 391501 from streams 559 and 508 in ref [1]
        Design['Drying air flow rate'] = 12105/391501 * self.ins[0].F_mass
        Design['Solids flow rate'] = self.outs[0].F_mass
        Design['Filtrate flow rate'] = self.outs[1].F_mass

class AcidulationReactor(Reactor):
    _N_ins = 2
    _N_outs = 1

    def _setup(self):
        super()._setup()
        self.acidulation_rxns = ParallelRxn([
        #   Reaction definition                                           Reactant        Conversion
        Rxn('CalciumLactate + H2SO4 -> 2 LacticAcid + CaSO4',         'CalciumLactate',       1),
        Rxn('CalciumAcetate + H2SO4 -> 2 AceticAcid + CaSO4',         'CalciumAcetate',       1),
        Rxn('CalciumSuccinate + H2SO4 -> SuccinicAcid + CaSO4',       'CalciumSuccinate',     1),
        Rxn('2 AmmoniumHydroxide + H2SO4 -> AmmoniumSulfate + 2 H2O', 'AmmoniumHydroxide',    1),
        Rxn('CalciumDihydroxide + H2SO4 -> CaSO4 + 2 H2O',            'CalciumDihydroxide',   1)
            ])

    bypass = False

    def _run(self):
        if self.bypass:
            self.ins[1].empty()
            self.outs[0].copy_like(self.ins[0])
        else:
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

            # Set feed acid mol to match acidulation needs with 10% extra
            acid.imol['H2SO4'] = needed_acid * 1.1
            acid.imass['H2O'] = acid.imass['H2SO4'] / 0.93 * 0.07 # 93% purity
            effluent.mix_from([feed, acid])
            rxns.adiabatic_reaction(effluent)

    def _design(self):
        if self.bypass: self.design_results.clear()
        else: super()._design()

    def _cost(self):
        if self.bypass: self.baseline_purchase_costs.clear()
        else: super()._cost()


# Filter to separate gypsum from the acidified fermentation broth
@cost(basis='Feed flow rate', ID='Hydrocyclone & rotary drum filter', units='kg/hr',
      # Size based on stream 239 in ref [3]
      cost=187567, S=272342, CE=CEPCI[1998], n=0.39, BM=1.4)
@cost(basis='Filtrate flow rate', ID='Filtered hydrolysate pump', units='kg/hr',
      # Size based on stream 230 in ref [3], power based on
      # SaccharificationAndCoFermentation Filtrate Saccharification transfer pumps
      kW=74.57*265125/421776, cost=31862, S=265125, CE=CEPCI[1997], n=0.79, BM=2.8)
class GypsumFilter(SolidsSeparator):
    _N_ins = 1
    _units = {'Feed flow rate': 'kg/hr',
              'Filtrate flow rate': 'kg/hr'}
    bypass = False

    def _run(self):
        if self.bypass:
            self.outs[0].empty()
            self.outs[1].copy_like(self.ins[0])
        else: super()._run()

    def _design(self):
        if self.bypass: self.design_results.clear()
        else:
            self.design_results['Feed flow rate'] = self.ins[0].F_mass
            self.design_results['Filtrate flow rate'] = self.outs[1].F_mass

    def _cost(self):
        if self.bypass: self.baseline_purchase_costs.clear()
        else: self._decorated_cost()


class Esterification(Reactor):
    """
    Create an esterification reactor that converts organic acids and ethanol
    to corresponding ethyl esters and water. Finds the amount of catalyst
    'Amberlyst-15' required as well as the loss occurred by physicochemical attrition.

    Parameters
    ----------
    ins :
        [0] Main broth
        [1] Recycled ethanol stream 1
        [2] Recycled lactic acid stream
        [3] Supplementary ethanol
        [4] Recycled ethanol stream 2

    outs :
        [0] Main effluent
        [1] Wastewater stream (discarded recycles)

    ethanol2acids : float
        Ethanol feed to total acid molar ratio.
    T : float
        Operating temperature, [K].
    catalyst_price : float
        Unit price of the Amberlyst-15 catalyst, [$/kg].
    """
    _N_ins = 5
    _N_outs = 2

    _F_BM_default = {
        **Reactor._F_BM_default,
        'Amberlyst-15 catalyst': 1,
        }

    cat_load = 0.039 # wt% of total mass

    ethanol2acids = 1.5

    # Used in uncertainty analysis to adjust the conversion
    X_factor = 1

    reactives = ('LacticAcid', 'Ethanol', 'H2O', 'EthylLactate',
                 'AceticAcid', 'EthylAcetate', 'SuccinicAcid', 'EthylSuccinate')

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=351.15, P=101325, tau=None, tau_max=15,
                 V_wf=0.8, length_to_diameter=2, mixing_intensity=None, kW_per_m3=1.97,
                 X1=None, X2=None, assumeX2equalsX1=True, allow_higher_T=False,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 catalyst_price=138.03, # price['Amberlyst15']
                 ):

        Unit.__init__(self, ID, ins, outs)

        self.T = T
        self.P = P
        self.tau_max = tau_max
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.X1, self.X2, self._tau = X1, X2, tau
        self.assumeX2equalsX1 = assumeX2equalsX1
        self.allow_higher_T = allow_higher_T
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.catalyst_price = catalyst_price
        ID = self.ID
        self._mixed = Stream(f'{ID}_mixed')
        self._tmp_flow = Stream(f'{ID}_tmp_flow')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)

    def compute_coefficients(self, T):
        K = self.K = exp(2.9625 - 515.13/T)
        kc = self.kc = 2.70 * (1e7) * exp(-6011.55/T)
        KW = self.KW = 15.19 * exp(12.01/T)
        KEt = self.KEt = 1.22 * exp(359.63/T)
        return K, kc, KW, KEt

    def compute_r(self, flow, reactives, T):
        lle_chemicals = self.chemicals.lle_chemicals
        lle_IDs = tuple([i.ID for i in lle_chemicals])
        f_gamma = tmo.equilibrium.DortmundActivityCoefficients(lle_chemicals)
        gammas = f_gamma(flow.get_normalized_mol(lle_IDs), T)
        gammas_reactives = np.array([gammas[lle_IDs.index(ID)]
                                     for ID in reactives])
        normalized_mol = flow.get_normalized_mol(lle_IDs)
        curr_conc = np.array([normalized_mol[lle_IDs.index(ID)]
                              for ID in reactives])
        activities = gammas_reactives * curr_conc
        r_numerator = self.kc * (activities[1]*activities[0]-
                                 (activities[3]*activities[2]/self.K))
        r_denominator = (1+self.KEt*activities[3]+self.KW*activities[2])**2
        r = r_numerator / r_denominator
        return r

    def compute_X1_and_tau(self, mixed_stream, time_step):
        T = self.T
        cat_load = self.cat_load
        reactives = self.reactives[0:4]
        compute_r = self.compute_r
        time_max = self.tau_max * 60 # tau_max in hr
        K, kc, KW, KEt = self.compute_coefficients(T)

        tmp_flow = self._tmp_flow
        tmp_flow.copy_like(mixed_stream)
        self.mcat = mcat = cat_load * tmp_flow.F_mass
        r = compute_r(tmp_flow, reactives, T)
        dX = r * time_step * mcat / 1000 # r is in mol g-1 min-1

        curr_flow = tmp_flow.get_flow('kmol/hr', reactives)
        new_flows = [1, 1, 1, 1]
        LA_initial = tmp_flow.imol['LacticAcid']

        tau_min = time_step # tau in min
        while dX/LA_initial>1e-4:
            if curr_flow[0]<dX or curr_flow[1]<dX:
                dX = min(curr_flow[0], curr_flow[1])

            new_flows = [curr_flow[0]-dX, # LA
                         curr_flow[1]-dX, # ethanol
                         curr_flow[2]+dX, # water
                         curr_flow[3]+dX] # EtLA

            tmp_flow.set_flow(new_flows, 'kmol/hr', reactives)

            # Zhao et al. 2008 reported 96% conversion of NH4LA -> BuLA in 6h
            if new_flows[0]<=0 or new_flows[1]<=0 or tau_min>time_max-time_step: break

            r = compute_r(tmp_flow, reactives, T)
            dX = r * time_step * mcat / 1000  # r is in mol g-1 min-1
            curr_flow = tmp_flow.get_flow('kmol/hr', reactives)
            tau_min += time_step

        LA_in_feeds = mixed_stream.imol['LacticAcid']
        X1 = (LA_in_feeds-tmp_flow.imol['LacticAcid']) / LA_in_feeds
        tau = tau_min / 60 # convert min to hr
        return X1, tau

    @property
    def tau(self):
        """Residence time [hr]."""
        return self._tau

    def _run(self):
        # On weight basis, ethanol1 is ~91% ethanol with 8% water,
        # ethanol2 is ~97% ethanol with ~3% water, so should first recycle ethanol2
        feed, ethanol1, recycled_LA, supplement_ethanol, ethanol2 = self.ins
        effluent, wastewater = self.outs

        acids = ('LacticAcid', 'AceticAcid', 'SuccinicAcid')
        # Succinic acid is a dicarboxylic acid, needs twice as much ethanol
        ratios = self.ethanol2acids * np.array([1, 1, 2])

        mixed = self._mixed
        mixed.mix_from([feed, recycled_LA])

        # Have enough ethanol in feed and ethanol2, discharge some ethanol2
        # and all of ethanol1
        if compute_extra_chemical(mixed, ethanol2, acids, 'Ethanol', ratios) > 0:
            effluent, ethanol2_discarded = \
                adjust_recycle(mixed, ethanol2, acids, 'Ethanol', ratios)
            wastewater.mix_from([ethanol1, ethanol2_discarded])
            supplement_ethanol.empty()
        else:
            # Recycle all of ethanol2 and some ethanol1
            mixed.mix_from([mixed, ethanol2])
            # Have enough ethanol in feed2 and ethanol1
            if compute_extra_chemical(mixed, ethanol1, acids, 'Ethanol', ratios) > 0:
                effluent, ethanol1_discarded = \
                    adjust_recycle(mixed, ethanol1, acids, 'Ethanol', ratios)
                wastewater = ethanol1_discarded
                supplement_ethanol.empty()
            # Not have enough ethanol in both ethanol2 and ethanol1,
            # need supplementary ethanol
            else:
                supplement_ethanol.imol['Ethanol'] = \
                    - compute_extra_chemical(mixed, ethanol1, acids, 'Ethanol', ratios)
                effluent.mix_from(self.ins)
                wastewater.empty()

        if self.allow_higher_T and self.T<effluent.T:
            self.T = effluent.T

        if self.X1 == None and self.tau == None:
            X1, tau = self.compute_X1_and_tau(effluent, time_step=1)
            self.X1 = X1
            self._tau = tau
        elif not self.X1 and not self.tau:
            raise AttributeError('X1 and tau must be both defined, or both as None')
        else:
            X1 = self.X1

        if self.assumeX2equalsX1:
            X2 = self.X2 = X1
        else:
            if not self.X2:
                raise AttributeError('X2 must be defined if assumeX2equalsX1 is False')

        X1 *= self.X_factor
        X2 *= self.X_factor

        self.esterification_rxns = ParallelRxn([
            #   Reaction definition                                     Reactant  Conversion
            Rxn('LacticAcid + Ethanol -> EthylLactate + H2O',         'LacticAcid',   X1),
            Rxn('AceticAcid + Ethanol -> EthylAcetate + H2O',         'AceticAcid',   X2),
            # Assume succinic acid has the same conversion as acetic acid
            Rxn('SuccinicAcid + 2 Ethanol -> EthylSuccinate + 2 H2O', 'SuccinicAcid', X2)
            ])

        self.esterification_rxns(effluent.mol)
        effluent.T = self.T
        self.outs[0].copy_like(effluent)
        self.outs[1].copy_like(wastewater)

    def _cost(self):
        super()._cost()
        
        hx = self.heat_exchanger
        # It does not matter if a fraction of the stream is discarded, it does not affect the energy balance.
        hx.ins[0].mix_from(self.ins)
        hx.outs[0].mix_from(self.outs)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)
        
        self.baseline_purchase_costs['Amberlyst-15 catalyst'] = self.mcat * self.catalyst_price

class HydrolysisReactor(Reactor):
    """
    Create a hydrolysis reactor that hydrolyze organic acid esters into
    corresponding acids and ethanol.

    Parameters
    ----------
    ins :
        [0] Main broth
        [1] Supplementary water
        [2] Recycled water stream 1
        [3] Recycled water stream 2

    outs :
        [0] Main effluent
        [1] Wastewater stream (discarded recycles)

    water2esters : float
        Water feed to total ester molar ratio.
    """
    _N_ins = 4
    _N_outs = 2
    water2esters = 12

    # Avoid having to write the _init again
    def _setup(self):
        super()._setup()
        self.hydrolysis_rxns = ParallelRxn([
            #   Reaction definition                                       Reactant   Conversion
            Rxn('EthylLactate + H2O -> LacticAcid + Ethanol',         'EthylLactate',   0.8),
            Rxn('EthylAcetate + H2O -> AceticAcid + Ethanol',         'EthylAcetate',   0.8),
            Rxn('EthylSuccinate + 2 H2O -> SuccinicAcid + 2 Ethanol', 'EthylSuccinate', 0.8),
            ])
        if not hasattr(self, '_mixed'): self._mixed = Stream(f'{self.ID}_mixed')


    def _run(self):
        # On weight basis, recycle2 is near 10% EtLA so will always be recycled,
        # but recycle1 is >97% water with <1% LA, so will only be used to supply
        # water needed for the hydrolysis reaction
        feed, water, recycle1, recycle2 = self.ins
        effluent, wastewater = self.outs

        esters = ('EthylLactate', 'EthylAcetate', 'EthylSuccinate')
        # Succnic acid is a dicarboxylic acid, needs twice as much water
        ratios = self.water2esters * np.array([1, 1, 2])
        # Have enough water in feed and recycle2, discharge some recycle2
        # and all of recycle1
        if compute_extra_chemical(feed, recycle2, esters, 'H2O', ratios) > 0:
            effluent, recycle2_discarded = \
                adjust_recycle(feed, recycle2, esters, 'H2O', ratios)
            wastewater.mix_from([recycle1, recycle2_discarded])
            water.empty()
            self._mixed.empty()
        else:
            # Recycle all of recycle2 and combine feed and recycle2 as feed2
            mixed = self._mixed
            mixed.mix_from([feed, recycle2])
            # Have enough water in feed2 and recycle1
            if compute_extra_chemical(mixed, recycle1, esters, 'H2O', ratios) > 0:
                effluent, recycle1_discarded = \
                    adjust_recycle(mixed, recycle1, esters, 'H2O', ratios)
                wastewater = recycle1_discarded
                water.empty()
            # Not have enough water in both recycles, need supplementary water
            else:
                water.imol['H2O'] = \
                    - compute_extra_chemical(mixed, recycle1, esters, 'H2O', ratios)
                effluent.mix_from(self.ins)
                wastewater.empty()

        rxns = self.hydrolysis_rxns
        rxns(effluent.mol)
        self.outs[0].copy_like(effluent)
        self.outs[1].copy_like(wastewater)


# %%

# =============================================================================
# Wastewater treatment
# =============================================================================

@cost(basis='COD flow', ID='Anaerobic basin', units='kg-O2/hr',
      kW=2371.326, cost=25800000, S=27211, CE=CEPCI[2012], n=0.6, BM=1.1)
class AnaerobicDigestion(Unit):
    _N_ins = 1
    _N_outs = 3
    _units= {'COD flow': 'kg-O2/hr'}

    auxiliary_unit_names = ('heat_exchanger',)

    def __init__(self, ID='', ins=None, outs=(), *, reactants, COD_chemicals,
                 split=(), T=35+273.15):
        Unit.__init__(self, ID, ins, outs)
        self.reactants = reactants if isinstance(reactants[0], str) else [i.ID for i in reactants]
        self.COD_chemicals = COD_chemicals
        self.split = split
        self.T = T
        ID = self.ID
        self._multi_stream = tmo.MultiStream(f'{ID}_ms')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        # Add '.' in ID for auxiliary units
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=T)

        # Based on P49 in ref [1], 91% of organic components is destroyed,
        # of which 86% is converted to biogas and 5% is converted to sludge,
        # and the biogas is assumed to be 51% CH4 and 49% CO2 on a dry molar basis
        chems = self.chemicals
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
        biogas.receive_vent(treated_water)
        sludge.mol = liquid_mol * (1 - self.split)
        biogas.T = treated_water.T = sludge.T = T

    def _design(self):
        self.design_results['COD flow'] = compute_COD(self.COD_chemicals, self.ins[0])
        
        hx = self.heat_exchanger
        hx.ins[0].mix_from(self.ins)
        hx.outs[0].copy_like(hx.ins[0])
        hx.outs[0].T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)


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

    def __init__(self, ID='', ins=None, outs=(), *, reactants, COD_chemicals,
                 caustic_mass=None, need_ammonia=False):
        Unit.__init__(self, ID, ins, outs)
        # In the case that the actual chemicals are provided
        self.reactants = reactants if isinstance(reactants[0], str) else [i.ID for i in reactants]
        self.COD_chemicals = COD_chemicals
        self.caustic_mass = caustic_mass or get_flow_tpd()/2205 * 2252
        self.need_ammonia = need_ammonia

        chems = self.chemicals
        def growth(reactant):
            f = chems.WWTsludge.MW / getattr(chems, reactant).MW
            return Rxn(f"{f}{reactant} -> WWTsludge", reactant, 1.)

        # Reactions from auto-populated combustion reactions.
        # Based on P49 in ref [1], 96% of remaining soluble organic matter
        # is removed after aerobic digestion, of which 74% is converted to
        # water and CO2 and 22% to cell mass
        combustion_rxns = chems.get_combustion_reactions()
        self.digestion_rxns = ParallelRxn([i*0.74 + 0.22*growth(i.reactant)
                                           for i in combustion_rxns
                                           if (i.reactant in self.reactants)])
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

        self.design_results['COD flow'] = compute_COD(self.COD_chemicals, effluent)
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

    def __init__(self, ID='', ins=None, outs=(), *, COD_chemicals, split):
        Unit.__init__(self, ID, ins, outs)
        self.COD_chemicals = COD_chemicals
        self.split = split

    def _run(self):
        mixture = self.ins[0].copy()
        mixture.mix_from(self.ins)
        mixture.split_to(*self.outs, self.split)

    def _design(self):
        self.design_results['Volumetric flow'] = self.outs[0].F_vol
        self.design_results['COD flow'] = compute_COD(self.COD_chemicals, self.ins[0])


@cost(basis='COD flow', ID='Thickeners', units='kg-O2/hr',
      kW=107.3808, cost=750000, S=5600, CE=CEPCI[2012], n=0.6, BM=1.6)
class BeltThickener(Unit):
    _ins_size_is_fixed = False
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}

    def __init__(self, ID='', ins=None, outs=(),
                 solubles=None, insolubles=None, *, COD_chemicals):
        Unit.__init__(self, ID, ins, outs)
        self.solubles = solubles or tuple(i.ID for i in self.chemicals.solubles)
        self.insolubles = insolubles or tuple(i.ID for i in self.chemicals.insolubles)
        self.COD_chemicals = COD_chemicals

    def _run(self):
        centrate, solids = self.outs
        solubles = self.solubles
        insolubles = self.insolubles

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
        self.design_results['COD flow'] = compute_COD(self.COD_chemicals, self.ins[0])


@cost(basis='COD flow', ID='Centrifuge', units='kg-O2/hr',
      # power usage including feed pumping and centrifuge
      kW=22.371+123.0405, cost=686800, S=5600, CE=CEPCI[2012], n=0.6, BM=2.7)
class SludgeCentrifuge(Unit):
    _N_ins = 1
    _N_outs = 2
    _units= {'COD flow': 'kg-O2/hr'}

    __init__ = BeltThickener.__init__

    def _run(self):
        influent = self.ins[0]
        centrate, solids = self.outs
        solubles = self.solubles
        insolubles = self.insolubles

        # Centrifuge captures 95% of the solids at 20% solids
        solids.imass[insolubles] = 0.95 * influent.imass[insolubles]
        solids.imass['Water'] = 0.8/0.2 * (influent.imass[insolubles].sum())
        if solids.imass['Water'] < influent.imass['Water']: # have enough water
            # Assume solubles have the same split as water
            ratio = solids.imass['Water'] / influent.imass['Water']
            solids.imass[solubles] = ratio * influent.imass[solubles]
            solids.T = influent.T
            centrate.mol = influent.mol - solids.mol
            centrate.T = influent.T
        else:
            centrate.empty()
            solids.copy_like(influent)

    _design = BeltThickener._design


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


# %%

# =============================================================================
# Storage
# =============================================================================

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
      cost=70000, S=1393, CE=CEPCI[2009], n=0.7, BM=2.6)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=7.457, cost=21200, S=1393, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=0.37285, cost=3000, S=1393, CE=CEPCI[2009], n=0.8, BM=3.1)
class CSLstorage(Unit): pass

# For storage of lime used in separation and waste treatment as in ref [3]
@cost(basis='Flow rate', ID='Storage bin', units='kg/hr',
      cost=136370, S=2395, CE=CEPCI[1997], n=0.46, BM=1.3)
# Cost not scaled, thus used n=0
# Power usage scaled based on M104 (truck dumper hopper) in ref [1]
@cost(basis='Flow rate', ID='Feeder', units='kg/hr',
      kW=37.285/3500*140, cost=3900, S=2395, CE=CEPCI[1997], n=0, BM=1.3)
# Power usage scaled based on M-106 (dust collection system) in ref [1]
@cost(basis='Flow rate', ID='Unloading blower', units='kg/hr',
      kW=18.6425*7425/8500, cost=99594, S=2395, CE=CEPCI[1998], n=0.5, BM=1.4)
@cost(basis='Flow rate', ID='Dust vent baghouse', units='kg/hr',
      cost=140707, S=2395, CE=CEPCI[1997], n=1, BM=1.5)
class LimeStorage(Unit): pass

@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=96000, S=1981, CE=CEPCI[2011], n=0.7, BM=1.5)
class CausticStorage(Unit): pass

# Modified from bst.units.StorageTank, which won't simulate for empty flow
class SpecialStorage(StorageTank):
    def _cost(self):
        if self.ins[0].F_mol == 0:
            self.design_results['Number of tanks'] = 0
            self.baseline_purchase_costs['Tanks'] = 0
        else: StorageTank._cost(self)
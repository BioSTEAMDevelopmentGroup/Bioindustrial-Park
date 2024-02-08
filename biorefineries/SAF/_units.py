#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 08:16:49 2023

@author: wenjun
"""

#%%
import numpy as np
import biosteam as bst
import thermosteam as tmo
from biosteam import Stream, Unit, BinaryDistillation
from biosteam.units import HXutility, Mixer, SolidsSeparator, Compressor
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
from thermosteam import MultiStream
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume

_lb2kg = 0.453592
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124 # auom('gallon').conversion_factor('m3')*60
_Gcal2kJ = 4184000 # auom('kcal').conversion_factor('kJ')*1e6 , also MMkcal/hr
# _316_over_304 = factors['Stainless steel 316'] / factors['Stainless steel 304']

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

CEPCI = bst.design_tools.CEPCI_by_year
#%% 

# Pretreatment
@cost('Dry flow rate', 'Pretreatment reactor system', units='kg/hr',
      S=83333, CE=522, cost=19812400 * 0.993, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactor(bst.units.design_tools.PressureVessel, Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.Flash._graphics
    _units = {'Residence time': 'hr',
              'Reactor volume': 'm3'}
    
    def _init(self, T=130+273.15, tau=0.166, V_wf=0.8, length_to_diameter=2, 
                 vessel_material='Stainless steel 316', 
                 vessel_type='Horizontal',
                 reactions=None,
                 run_vle=True):
        self._load_components()
        vapor, liquid = self.outs
        vapor.phase = 'g'
        self.T = T
        chemicals = self.chemicals
        if reactions is None:
            self.reactions = ParallelRxn([
        #            Reaction definition                 Reactant    Conversion
        Rxn('Glucan + H2O -> Glucose',                   'Glucan',   0.0990, chemicals),
        Rxn('Glucan + H2O -> GlucoseOligomer',           'Glucan',   0.0030, chemicals),
        Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0030, chemicals),
        # Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0240, chemicals),
        # Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030, chemicals),
        # Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.0030, chemicals),
        # Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030, chemicals),
        Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.0000, chemicals),
        Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000, chemicals),
        Rxn('Xylan + H2O -> XyloseOligomer',             'Xylan',    0.0240, chemicals),
        Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0500, chemicals),
        Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000, chemicals),
        Rxn('Arabinan + H2O -> ArabinoseOligomer',       'Arabinan', 0.0240, chemicals),
        Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050, chemicals),
        Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000, chemicals),
        Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0500, chemicals)
            ])
            # self.glucan_to_glucose = self.reactions[0]
            # self.xylan_to_xylose = self.reactions[8]
            # self.glucose_to_byproducts = self.reactions[1:3]
            # self.xylose_to_byproducts = self.reactions[9:12]
        else:
            self.reactions = reactions
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.run_vle = run_vle
        
    def _load_components(self):
        thermo = self.thermo
        self._multistream = MultiStream(None, thermo=thermo)
    
    def _run(self):
        feed = self.ins[0]
        vapor, liquid = self.outs
        liquid.copy_like(feed)
        self.reactions.adiabatic_reaction(liquid) 
        if self.T:
            if self.run_vle:
                ms = self._multistream
                ms.copy_like(liquid)
                ms.vle(T=self.T, H=ms.H)
                vapor.mol[:] = ms.imol['g']
                liquid.mol[:] = ms.imol['l']
                vapor.T = liquid.T = ms.T
                vapor.P = liquid.P = ms.P
            else:
                liquid.T = self.T

    def _design(self):
        Design = self.design_results
        ins_F_vol = self.F_vol_in
        V_reactor = ins_F_vol * self.tau / self.V_wf
        P = self.outs[0].P * 0.000145038 # Pa to psi
        length_to_diameter = self.length_to_diameter
        D = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter)
        D *= 3.28084 # convert from m to ft
        L = D * length_to_diameter
        Design['Residence time'] = self.tau
        Design['Reactor volume'] = V_reactor
        Design.update(self._vessel_design(float(P), float(D), float(L)))
        self._decorated_design()
            
    def _cost(self):
        Design = self.design_results
        self.baseline_purchase_costs.update(
            self._vessel_purchase_cost(
                Design['Weight'], Design['Diameter'], Design['Length']
            )
        )
        self._decorated_cost()


# Fermentation

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





@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_reactors')
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 5
    _N_outs = 3
    

    # #: Saccharification temperature (K)
    # T_saccharification = 48+273.15

    # #: Fermentation temperature (K)
    # T_fermentation = 32+273.15
    
    #: Saccharification and Co-Fermentation time (hr)
    tau = 72
    
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    #: Number of reactors
    N_reactors = 12
    
    _units = {'Flow rate': 'm3/hr',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    auxiliary_unit_names = ('heat_exchanger',)

    # Split to outs[2]
    inoculum_ratio = 0.1

    def __init__(self, ID='', ins=None, outs=(), T=34+273.15, P=101325):
        Unit.__init__(self, ID, ins, outs)
        self.T = T
        self.P = P
        ID = self.ID
        self.saccharified_stream = Stream(f'{ID}_ss')
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out,T=T)

        self.saccharification_rxns = ParallelRxn([
    #   Reaction definition                   Reactant     Conversion
    Rxn('Glucan -> GlucoseOligomer',          'Glucan',      0.04),
    Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',      0.012),
    Rxn('Glucan + H2O -> Glucose',            'Glucan',      0.9),
    Rxn('Cellobiose + H2O -> Glucose',        'Cellobiose',  1),
    Rxn('Xylan + H2O -> Xylose',              'Xylan',       0.9),
    Rxn('Arabinan + H2O -> Arabinose',        'Arabinan',    0.85)])

        self.loss_rxns = ParallelRxn([
    #   Reaction definition               Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.03),
    Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.03),
    Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.03),
    Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.03),
    Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.03),])

        self.cofermentation_rxns = ParallelRxn([
    #   Reaction definition                                          Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.95),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.02),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.004),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.006),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.85),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',    'Xylose',    0.019),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.003),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.046),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.009),
    ])

    def _run(self):
        feed, inoculum, CSL, DAP, concentrated_juice = self.ins
        vent, effluent, sidedraw = self.outs
        ss = self.saccharified_stream
        
        # 0.25 wt% and 0.33 g/L (kg/m3) based on ref [1]
        CSL.imass['CSL'] = 0.0025 * (feed.F_mass + concentrated_juice.F_mass)
        DAP.imass['DAP'] = 0.33 * (feed.F_vol + concentrated_juice.F_vol)
        
        ss.mix_from((feed, inoculum, CSL, DAP, concentrated_juice))
        ss.T = sidedraw.T = effluent.T = self.T
        
        vent.P = effluent.P = sidedraw.P = self.P
        vent.phase = 'g'
        
        
        self.saccharification_rxns(ss.mol)
        # Sidedraw to seedtrain
        sidedraw.mol = ss.mol * self.inoculum_ratio
        effluent.mol = ss.mol - sidedraw.mol
        
        
        self.loss_rxns(effluent.mol)
        self.cofermentation_rxns(effluent.mol)
        
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Flow rate'] = v_0 / self.N_reactors
        Design.update(size_batch(v_0, self.tau, self.tau_0, self.N_reactors, self.V_wf))

        
        hx = self.heat_exchanger
        hx.ins[0].mix_from(self.ins[0:4])
        hx.outs[0].mix_from(self.outs[0:2])
        hx.outs[0].T = self.T
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins,outs=hx.outs)



@cost('Flow rate', 'Pumps',
      S=43149, CE=522, cost=24800, n=0.8, kW=40, BM=2.3)
@cost('Stage #1 reactor volume', 'Stage #1 reactors',
      cost=37700, S=20*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #2 reactor volume', 'Stage #2 reactors',
      cost=58300, S=200*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #3 reactor volume', 'Stage #3 reactors',
      cost=78800, S=2e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 reactors',
      cost=176e3, S=20e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #4 reactor volume', 'Stage #4 agitators',
      cost=26e3/2, S=20e3*_gal2m3, kW=7.5, CE=522, n=0.5, BM=1.5)
@cost('Stage #5 reactor volume', 'Stage #5 reactors',
      cost=590e3, S=200e3*_gal2m3, CE=522, n=0.7, BM=1.8)
@cost('Stage #5 reactor volume', 'Stage #5 agitators',
      cost=43e3/2, S=200e3*_gal2m3, kW=10, CE=522, n=0.5, BM=1.5)
class SeedTrain(Unit):
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    
    _units= {'Flow rate': 'kg/hr',
             'Stage #1 reactor volume': 'm3',
             'Stage #2 reactor volume': 'm3',
             'Stage #3 reactor volume': 'm3',
             'Stage #4 reactor volume': 'm3',
             'Stage #5 reactor volume': 'm3'}
    
    @property
    def N_stages(self): 
        """Number of stages in series."""
        return 5
    
    #: Number of parallel seed trains
    N_trains = 2
    
    #: Cycle time for each batch (hr)
    tau_batch = 24
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
   
    def _init(self,T = 32+273.15,fermentation_rxns=None):
        self.T = T
        self.fermentation_rxns = ParallelRxn([
        #   Reaction definition                             Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9000),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                            'Glucose',   0.0400),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8000),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                            'Xylose',    0.0400),
        Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460),
        Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090)
            ])
    _setup = Unit._setup
    
    def _run(self):
        feed, CSL, DAP = self.ins
        vent, effluent= self.outs
        
        # 0.50 wt% and 0.66 g/L (kg/m3) from ref [1]
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
        maxvol = self.outs[1].F_vol*self.tau_turnover
        vol = maxvol*10**-self.N_stages
        Design = self.design_results
        for i in range(1, self.N_stages+1):
            Design[f'Stage #{i} reactor volume'] = vol
            vol *= 10 
        Design['Flow rate'] = sum([i.F_mass for i in self.outs])
        self.add_heat_utility(self.Hnet, self.T)

    def _cost(self):
        N = self.N_trains
        D = self.design_results
        C = self.baseline_purchase_costs
        kW = 0
        for i, x in self.cost_items.items():
            S = D[x._basis]
            q = S/x.S
            C[i] = N*bst.CE/x.CE*x.cost*q**x.n
            kW += N*x.kW*q
        self.power_utility(kW)





@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=439000, S=40414, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=11.1855, cost=31800, S=40414, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=7.457, cost=8200, S=43149, CE=CEPCI[2009], n=0.8, BM=2.3)
class SeedHoldTank(Unit): pass





@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=636000, S=425878, CE=CEPCI[2009], n=0.7, BM=1.8)
@cost(basis='Flow rate', ID='Agitator', units='kg/hr',
      kW=29.828, cost=68300, S=425878, CE=CEPCI[2009], n=0.5, BM=1.5)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=26800, S=488719, CE=CEPCI[2009], n=0.8, BM=2.3)
class BeerTank(Unit): pass





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
      # Power not specified, from filtrate tank discharge pump
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

#%%

from math import pi, ceil
from biosteam.units.design_tools import PressureVessel
from biosteam.exceptions import DesignError
class Reactor(Unit, PressureVessel, isabstract=True):
    '''    
    Create an abstract class for reactor unit, purchase cost of the reactor
    is from volume calculated by residence time.

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
    
    # For a single reactor, from diameter and length from PressureVessel._bounds,
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
        Design['Reactor volume'] = V_reactor
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
        purchase_costs = self.baseline_purchase_costs
        
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



    






#%% Upgrading

### Dehydration
class AdiabaticFixedbedDehydrationReactor(Reactor):
    _N_ins=2
    _N_outs=2
    
    auxiliary_unit_names=('heat exchanger')
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3'}
    # from Rules of thumb in engineering practice,
    # for a single fixed bed of solid catalyst reactor, 
    # volume of reactor 1–10 000 L
    _V_max = 10
    
    _F_BM_default={**Reactor._F_BM_default,
                   'Heat exchanger':3.17,
                   'SynDol catalyst':1.0} 
    # from Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation : 
    # main components of Al2O3-MgO/SiO2
    
    def __init__(self, ID, ins, outs,
                 T=350+273.15, # from Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
                 P=1.4*101325, # from Bioethylene Production from Ethanol: A Review and Techno-economical Evaluation
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=0.43, # from 	Catalytic dehydration of bioethanol to ethylene
                 tau=3.14/3600, # from Process Design for the Production of Ethylene from Ethanol
                 catalyst_price=20.48, # Price ['Syndol price']
                 **wargs): 
        Reactor.__init__(self, ID, ins, outs)
        self.T=T
        self.P=P
        self.tau=tau
        self.length_to_diameter = length_to_diameter #: Length to diameter ratio
        self.vessel_material = vessel_material # Vessel material
        self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
        self.WHSV = WHSV # Residence time in kg/hr feed / kg catalyst
        self.catalyst_price = catalyst_price # In USD/kg
        self.heat_exchanger = hx = HXutility(None, None,None, T=T)
    
    def _setup(self):
        super()._setup()
        self.overall_C2H5OH_conversion=overall_C2H5OH_conversion=0.995
        self.C2H5OH_to_C2H4_selectivity=C2H5OH_to_C2H4_selectivity=0.988
        self.C2H5OH_to_CH3OCH3_selectivity=C2H5OH_to_CH3OCH3_selectivity=0.00052
        self.C2H5OH_to_CH3CHO_selectivity=C2H5OH_to_CH3CHO_selectivity=0.002
        self.C2H5OH_to_C2H6_selectivity=C2H5OH_to_C2H6_selectivityy=0.0027
        self.C2H5OH_to_C3H6_selectivity=C2H5OH_to_C3H6_selectivity=0.0006
        self.C2H5OH_to_C4H8_selectivity=C2H5OH_to_C4H8_selectivity=0.005
        self.C2H5OH_to_CO_selectivity=C2H5OH_to_CO_selectivity=0.00007
        self.C2H5OH_to_CO2_selectivity=C2H5OH_to_CO2_selectivity=0.0011
        
        self.dehydration_rxns=ParallelRxn([
        #   Reaction definition                    Reactant    Conversion
        Rxn('C2H5OH -> C2H4 + H2O',               'C2H5OH',   C2H5OH_to_C2H4_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH -> CH3OCH3 + H2O',          'C2H5OH',   C2H5OH_to_CH3OCH3_selectivity*overall_C2H5OH_conversion),
        Rxn('C2H5OH -> CH3CHO +H2',               'C2H5OH',   C2H5OH_to_CH3CHO_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH + H2 -> 2 C2H6 +2 H2O',     'C2H5OH',   C2H5OH_to_C2H6_selectivityy*overall_C2H5OH_conversion),
        Rxn('3 C2H5OH -> C3H6 + 3 H2O',           'C2H5OH',   C2H5OH_to_C3H6_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH -> C4H8 + 2 H2O',           'C2H5OH',   C2H5OH_to_C4H8_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH -> CO + CH4 + H2',          'C2H5OH',   C2H5OH_to_CO_selectivity*overall_C2H5OH_conversion),
        Rxn('C2H5OH +H2O -> CO2 + CH4 + H2',      'C2H5OH',   C2H5OH_to_CO2_selectivity*overall_C2H5OH_conversion),
            ])
        
    def _run(self):
        feed, fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        
        effluent.copy_like(feed)
        self.dehydration_rxns(effluent.mol)
        effluent.T=self.T
        effluent.P=self.P
        
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.mix_from(self.outs)
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=False)
        
    def _cost(self):
        super()._cost()
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming SynDol lifetime of 12 months, enough for one year
        self.purchase_costs['SynDol catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities = hx.heat_utilities
       
    




class CausticTower(bst.Unit):
    _N_ins=2
    _N_outs=1
    
    _outs_size_is_fixed=False
    
    _units={'Feed volumetric flow':'m3/hr',
            'Caustic flow':'kg/hr'}

    def __init__(self, ID='',ins=None,outs=(),P=27*101325):
        Unit.__init__(self, ID,ins,outs)
        self.P=P
                                     #   Reaction definition          Reactant    Conversion
        self.neutralization_rxn=Rxn('2 NaOH + CO2 -> Na2CO3 + H2O',   'CO2',     1)

    def _run(self):
        influent,caustic=self.ins
        effluent=self.outs[0]
        
        caustic.imol['NaOH'] = 2*influent.imol['CO2']
        caustic.imass['H2O'] = caustic.imass['NaOH']
        
        effluent.copy_like(self.ins[0])
        effluent.P=self.P
        effluent.phase='g'
        effluent.imol['CO2']=0
    
    def _design(self):
        Design=self.design_results
        Design['Feed volumetric flow']=self.F_vol_in
        Design['Caustic flow']=self.ins[1].F_mass
        


            
        
   
    
    
    
    
    
    
### Oligomerization 

# First and second oligomerization are from 'Aviation fuel production pathways from lignocellulosic biomass via alcohol intermediates – A technical analysis' and 'Systems and processes for conversion of ethylene feedstocks to hydrocarbon fuels'

class EthyleneDimerizationReactor(Reactor): 
    _N_ins=2
    _N_outs=2
    
    ins_size_is_fixed = False 
    
    auxiliary_unit_names = ('heat_exchanger')
    
    _units={**PressureVessel._units,
            'Residence time':'hr',
            'Total Volume':'m^3',
            'Reactor volume':'m^3',
            'Catalyst weight':'kg'}
    
    # For a single reactor, from diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147 
    
    _F_BM_default={**Reactor._F_BM_default,
                   'Heat exchanger':3.17,
                   'Ni-loaded aluminosilicate':1} # a catalyst from nickel on an amorphous aluminum silicate support (SiO2-Al2O3)
     
    def __init__(self, ID, ins, outs,
                 T=85+273.15,
                 P=101325, # reactor pressure = 22 bar, maintained by N2
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=5, 
                 tau=1/5,
                 catalyst_price=12.3, # Price['Ni-loaded aluminosilicate catalyst']
                 **kwargs): 
        Reactor.__init__(self, ID, ins,outs)
        self.T=T
        self.P=P
        self.length_to_diameter=length_to_diameter
        self.vessel_material=vessel_material
        self.vessel_type=vessel_type
        self.WHSV=WHSV
        self.tau=tau
        self.catalyst_price=catalyst_price
        self.heat_exchanger=hx=HXutility(None,None,None,T=T)
       
    def _setup(self):
        super()._setup()
        self.C2H4_to_C4H8_conversion=C2H4_to_C4H8_conversion=0.988
        self.C2H4_to_C6H12_conversion=C2H4_to_C6H12_conversion=0.00052
        self.C2H4_to_C8H16_conversion=C2H4_to_C8H16_conversion=0.002
        self.C2H4_to_C10H20_conversion=C2H4_to_C10H20_conversion=0.003
        self.oligomerization_rxns=ParallelRxn([
            #   Reaction definition               Reactant       Conversion
            Rxn('2 C2H4 -> C4H8',                 'C2H4',     C2H4_to_C4H8_conversion),
            Rxn('3 C2H4 -> C6H12',                'C2H4',     C2H4_to_C6H12_conversion),
            Rxn('4 C2H4 -> C8H16',                'C2H4',     C2H4_to_C8H16_conversion),
            Rxn('5 C2H4 -> C10H20',               'C2H4',     C2H4_to_C10H20_conversion),
                ])
        self.C2H4_to_C4H8_rxn=self.oligomerization_rxns[0]
        self.C2H4_to_C6H12_rxn=self.oligomerization_rxns[1]
        self.C2H4_to_C8H16_rxn=self.oligomerization_rxns[2]
        self.C2H4_to_C10H20_rxn=self.oligomerization_rxns[3]
        
    def _run(self):
        feed,fresh_catalyst=self.ins
        effluent,spent_catalyst=self.outs
        
        effluent.copy_like(feed)
        self.oligomerization_rxns(effluent.mol)
        effluent.T=self.T
        effluent.P=self.P
        
        spent_catalyst.phase='s'
            
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.mix_from(self.outs)
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=True)
        
    def _cost(self):
        super()._cost()
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming catalyst lifetime of 12 months, enough for one year
        self.purchase_costs['1st catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities = hx.heat_utilities
    
       










# Intermediate to C12~ 
class OligomerizationReactor(Reactor):
    _N_ins=3
    _N_outs=2
    
    auxiliary_unit_names = ('heat_exchanger')
    
    _units={**PressureVessel._units,
            'Residence time':'hr',
            'Total Volume':'m^3',
            'Reactor volume':'m^3',
            'Catalyst weight':'kg'}
    
    # For a single reactor, from diameter and length from PressureVessel._bounds,
    # converted from ft3 to m3
    _V_max = pi/4*(20**2)*40/35.3147 
    
    _F_BM_default={**Reactor._F_BM_default,
                   'Heat exchanger':3.17,
                   'SiO2-Al2O3 catalyst':1} # Acid catalyst 
    def __init__(self,ID,ins=None,outs=(),thermo=None,*,
                 T=225+273.15,
                 P=101325,
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=3,
                 tau=1/3,
                 catalyst_price=12.3,
                 ):
        Reactor.__init__(self,ID,ins,outs)
        self.T=T
        self.P=P
        self.length_to_diameter=length_to_diameter
        self.vessel_material=vessel_material
        self.vessel_type=vessel_type
        self.WHSV=WHSV
        self.tau=tau
        self.catalyst_price=catalyst_price
        self.heat_exchanger=hx=HXutility(None, None,None,T=T)
        
    def _run(self):
        feed, recycle, fresh_catalyst=self.ins
        effluent, spent_catalyst=self.outs
        # effluent. = ['C4H8','C5H10','C6H12','C7H14','C8H16','C9H18','C10H20','C11H22',
        # 'C12H24','C13H26','C14H28','C16H32','C18H36']
        
        effluent.imass['C4H8']=0.097*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C5H10']=0.016*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C6H12']=0.169*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C7H14']=0.002*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C8H16']=0.21*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C9H18']=0.022*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C10H20']=0.175*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C11H22']=0.012*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C12H24']=0.123*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C13H26']=0.011*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C14H28']=0.079*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C16H32']=0.052*(self.ins[0].F_mass+self.ins[1].F_mass)
        effluent.imass['C18H36']=0.013*(self.ins[0].F_mass+self.ins[1].F_mass)
        
        effluent.T=self.T
        effluent.phase='l'
        spent_catalyst.phase='s'
            
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=True)
        
    def _cost(self):
        super()._cost()
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming catalyst lifetime of 12 months, enough for one year
        self.purchase_costs['2nd catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
         
        
        
        
    
    
    


### Hydrogenation 

# Two catalysts commonly used: 
# 1. Cobalt molybdenum (Como) catalyst. Reaction parameters are from 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks'
# 2. Pd/Al2O3 catalyst. Reaction parameters are from 'comparative techno-economic analysis and process design for indirect liquefaction pathways to distillate-range fuels via biomass-derived oxygenated intermediates upgrading'

# 1 has lower price of catalyst but higher temperature
# 2 has higher price of catalyst but longer lifetime and lower temperature
# Twp ways are considered here as a variable

class HydrogenationReactor(Reactor): # 1
    _N_ins = 3
    _N_outs = 2
    
    auxiliary_unit_names=('heat_exchanger')


    _F_BM_default = {**Reactor._F_BM_default,
            'Heat exchangers': 3.17,
            'Como catalyst': 1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=45.+273.15,
                 P=3.6*10E6, 
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=3, 
                 tau=1/3,
                 catalyst_price=19/_lb2kg, # lifetime 3 yrs 
                 **kargs):
        Reactor.__init__(self, ID,ins,outs)
        self.T = T
        self.P = P
        self.length_to_diameter = length_to_diameter
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.WHSV=WHSV
        self.tau=tau
        self.catalyst_price=catalyst_price
        self.heat_exchanger=hx=bst.HXutility(None,None,None,T=T)
    
    def _setup(self):
        super()._setup()
        self.hydrogenation_rxns = ParallelRxn([
            #   Reaction definition           Reactant   Conversion
            Rxn('C6H12 + H2 -> C6H14',         'C6H12',    1.00),
            Rxn('C7H14 + H2 -> C7H16',         'C7H14',    1.00),
            Rxn('C8H16 + H2 -> C8H18',         'C8H16',    1.00),
            Rxn('C9H18 + H2 -> C9H20',         'C9H18',    1.00),
            Rxn('C10H20 + H2 -> C10H22',       'C10H20',   1.00),
            Rxn('C11H22 + H2 -> C11H24',       'C11H22',   1.00),
            Rxn('C12H24 + H2 -> C12H26',       'C12H24',   1.00),
            Rxn('C13H26 + H2 -> C13H28',       'C13H26',   1.00),
            Rxn('C14H28 + H2 -> C14H30',       'C14H28',   1.00),
            Rxn('C16H32 + H2 -> C16H34',       'C16H32',   1.00),
            Rxn('C18H36 + H2 -> C18H38',       'C18H36',   1.00),
                ])                                      
        
    def _run(self):
        influent, h2 ,fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        
        h2.imol['H2'] = self.ins[0].F_mol
        effluent.mix_from([influent,h2])
        effluent.T=self.T
        self.hydrogenation_rxns(effluent.mol)
        
        fresh_catalyst.phase='s'
        
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=True)
        
    def _cost(self):
        super()._cost()
        
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming catalyst lifetime of 12 months
        self.purchase_costs['como_catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities



class HydrogenationReactor2(Reactor): # 2
    _N_ins = 3
    _N_outs = 2
    
    auxiliary_unit_names=('heat_exchanger')


    _F_BM_default = {**Reactor._F_BM_default,
            'Heat exchangers': 3.17,
            # 'Compressor':2.15,
            'Pd/Al2O3 catalyst': 1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 T=45.+273.15,
                 P=34.5*101325, 
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=1, 
                 tau=1,
                 catalyst_price=55.20/_lb2kg, # lifetime 3 yrs 
                 **kargs):
        Reactor.__init__(self, ID,ins,outs)
        self.T = T
        self.P = P
        self.length_to_diameter = length_to_diameter
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.WHSV=WHSV
        self.tau=tau
        self.catalyst_price=catalyst_price
        self.heat_exchanger=hx=bst.HXutility(None,None,None,T=T)
    
    def _setup(self):
        super()._setup()
        self.hydrogenation_rxns = ParallelRxn([
            #   Reaction definition           Reactant   Conversion
            Rxn('C6H12 + H2 -> C6H14',         'C6H12',    1.00),
            Rxn('C7H14 + H2 -> C7H16',         'C7H14',    1.00),
            Rxn('C8H16 + H2 -> C8H18',         'C8H16',    1.00),
            Rxn('C9H18 + H2 -> C9H20',         'C9H18',    1.00),
            Rxn('C10H20 + H2 -> C10H22',       'C10H20',   1.00),
            Rxn('C11H22 + H2 -> C11H24',       'C11H22',   1.00),
            Rxn('C12H24 + H2 -> C12H26',       'C12H24',   1.00),
            Rxn('C13H26 + H2 -> C13H28',       'C13H26',   1.00),
            Rxn('C14H28 + H2 -> C14H30',       'C14H28',   1.00),
            Rxn('C16H32 + H2 -> C16H34',       'C16H32',   1.00),
            Rxn('C18H36 + H2 -> C18H38',       'C18H36',   1.00),
                ])                                      
        
    def _run(self):
        influent, h2 ,fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        
        effluent.mix_from([influent,h2])
        
        # effluent.imol['C6H14']=self.ins[0].imol['C6H12']
        # effluent.imol['C7H16']=self.ins[0].imol['C7H14']
        # effluent.imol['C8H18']=self.ins[0].imol['C8H16']
        # effluent.imol['C9H20']=self.ins[0].imol['C9H18']
        # effluent.imol['C10H22']=self.ins[0].imol['C10H20']
        # effluent.imol['C11H24']=self.ins[0].imol['C11H22']
        # effluent.imol['C12H26']=self.ins[0].imol['C12H24']
        # effluent.imol['C13H28']=self.ins[0].imol['C13H26']
        # effluent.imol['C14H30']=self.ins[0].imol['C14H28']
        # effluent.imol['C16H34']=self.ins[0].imol['C16H32']
        # effluent.imol['C18H38']=self.ins[0].imol['C18H36']
       
        self.hydrogenation_rxns(effluent.mol)
        effluent.T=self.T
        
        h2.imol['H2'] = self.ins[0].F_mol
        
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=True)
        
    def _cost(self):
        super()._cost()
        
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming catalyst lifetime of 36 months
        self.purchase_costs['Pd_Al2O3_catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities += hx.heat_utilities
        








@cost(basis='Flow rate', ID='Tank', units='kg/hr',
      cost=803000, S=8343, CE=CEPCI[2009], n=0.7, BM=1.7)
@cost(basis='Flow rate', ID='Pump', units='kg/hr',
      kW=93.2125, cost=15000, S=8343, CE=CEPCI[2009], n=0.8, BM=3.1)
class FireWaterTank(Unit): pass





# Units in configuration 2 : Guerbet coupling of ethanol

class GuerbetCouplingReactor(Reactor):
# From "Advanced fuels from ethanol – a superstructure optimization approach."
    _N_ins=2
    _N_outs=2
    
    auxiliary_unit_names=('heat exchanger')
    
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Total volume': 'm3',
              'Reactor volume': 'm3'}
    # from Rules of thumb in engineering practice,
    # for a single fixed bed of solid catalyst reactor, 
    # volume of reactor 1–10 000 L
    _V_max = 10
    
    _F_BM_default={**Reactor._F_BM_default,
                   'Heat exchanger':3.17,
                   'Ca-HAP catalyst':1.0} 
    # from 
    
    def __init__(self, ID, ins, outs,
                 T=350+273.15, 
                 P=1.*101325, 
                 length_to_diameter=3,
                 vessel_material='Stainless steel 304',
                 vessel_type='Vertical',
                 WHSV=0.05, 
                 tau=1.78/3600, # From "Review of catalytic systems and thermodynamics for the Guerbet condensation reaction and challenges for biomass valorization."
                 catalyst_price=8, 
                 **wargs): 
        Reactor.__init__(self, ID, ins, outs)
        self.T=T
        self.P=P
        self.tau=tau
        self.length_to_diameter = length_to_diameter #: Length to diameter ratio
        self.vessel_material = vessel_material # Vessel material
        self.vessel_type = vessel_type # 'Horizontal' or 'Vertical'
        self.WHSV = WHSV # Residence time in kg/hr feed / kg catalyst
        self.catalyst_price = catalyst_price # In USD/kg
        self.heat_exchanger = hx = HXutility(None, None,None, T=T)
    
    def _setup(self):
        super()._setup()
        self.overall_C2H5OH_conversion=overall_C2H5OH_conversion=0.75
        self.C2H5OH_to_C2H4_selectivity=C2H5OH_to_C2H4_selectivity=0.030/0.75
        self.C2H5OH_to_C4H8_selectivity=C2H5OH_to_C4H8_selectivity=0.034/0.75
        self.C2H5OH_to_C6H12_selectivity=C2H5OH_to_C6H12_selectivity=0.016/0.75
        self.C2H5OH_to_C4H6_selectivity=C2H5OH_to_C4H6_selectivity=0.047/0.75
        self.C2H5OH_to_C6H10_selectivity=C2H5OH_to_C6H10_selectivity=0.035/0.75
        self.C2H5OH_to_C4H9OH_selectivity=C2H5OH_to_C4H9OH_selectivity=0.031/0.75
        self.C2H5OH_to_C5H11OH_selectivity=C2H5OH_to_C5H11OH_selectivity=0.267/0.75
        self.C2H5OH_to_C6H13OH_selectivity=C2H5OH_to_C6H13OH_selectivity=0.112/0.75
        self.C2H5OH_to_C7H15OH_selectivity=C2H5OH_to_C7H15OH_selectivity=0.096/0.75
        self.C2H5OH_to_C8H17OH_selectivity=C2H5OH_to_C8H17OH_selectivity=0.142/0.75
        
        self.Guerbet_rxn=ParallelRxn([
        #   Reaction definition                    Reactant    Conversion
        Rxn('C2H5OH -> C2H4 + H2O',               'C2H5OH',   C2H5OH_to_C2H4_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH -> C4H8 + 2 H2O',           'C2H5OH',   C2H5OH_to_C4H8_selectivity*overall_C2H5OH_conversion),
        Rxn('3 C2H5OH -> C6H12 + 3H2O',           'C2H5OH',   C2H5OH_to_C6H12_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH + H2 -> 2 C2H6 +2 H2O',     'C2H5OH',   C2H5OH_to_C4H6_selectivity*overall_C2H5OH_conversion),
        Rxn('3 C2H5OH -> C3H6 + 3 H2O',           'C2H5OH',   C2H5OH_to_C6H10_selectivity*overall_C2H5OH_conversion),
        Rxn('',           'C2H5OH',   C2H5OH_to_C4H9OH_selectivity*overall_C2H5OH_conversion),
        Rxn('2 C2H5OH -> CO + CH4 + H2',          'C2H5OH',   C2H5OH_to_C5H11OH_selectivity*overall_C2H5OH_conversion),
        Rxn('C2H5OH +H2O -> CO2 + CH4 + H2',      'C2H5OH',   C2H5OH_to_C6H13OH_selectivity*overall_C2H5OH_conversion),
            ])
        
    def _run(self):
        feed, fresh_catalyst = self.ins
        effluent, spent_catalyst = self.outs
        
        effluent.copy_like(feed)
        self.dehydration_rxns(effluent.mol)
        effluent.T=self.T
        effluent.P=self.P
        
    def _design(self):
        Reactor._design(self)
        duty=sum([i.H for i in self.outs])-sum([i.H for i in self.ins])
        feed=tmo.Stream()
        feed.mix_from(self.outs)
        feed.T=self.ins[0].T
        self.heat_exchanger.simulate_as_auxiliary_exchanger(ins=(feed,),
                                                            duty=duty,
                                                            vle=False)
        
    def _cost(self):
        super()._cost()
        self.catalyst_weight = self.ins[0].F_mass/self.WHSV 
        # assuming SynDol lifetime of 12 months, enough for one year
        self.purchase_costs['SynDol catalyst']=self.catalyst_weight * self.catalyst_price
        
        hx = self.heat_exchanger
        self.baseline_purchase_costs['Heat exchangers'] = hx.purchase_cost
        self.heat_utilities = hx.heat_utilities
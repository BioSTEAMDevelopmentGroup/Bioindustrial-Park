# -*- coding: utf-8 -*-
"""
.. contents:: :local:

Reactors
--------
.. autoclass:: biorefineries.cellulosic.units.PretreatmentReactorSystem
.. autoclass:: biorefineries.cellulosic.units.SeedTrain
.. autoclass:: biorefineries.cellulosic.units.ContinuousPresaccharification
.. autoclass:: biorefineries.cellulosic.units.Saccharification
.. autoclass:: biorefineries.cellulosic.units.CoFermentation
.. autoclass:: biorefineries.cellulosic.units.SaccharificationAndCoFermentation
.. autoclass:: biorefineries.cellulosic.units.SimultaneousSaccharificationAndCoFermentation

Separations
-----------
.. autoclass:: biorefineries.cellulosic.units.ReverseOsmosis
.. autoclass:: biorefineries.cellulosic.units.Nanofilter
.. autoclass:: biorefineries.cellulosic.units.HydrolyzateSolidLiquidSeparator
.. autoclass:: biorefineries.cellulosic.units.PretreatmentFlash

Feedstock handling
------------------
.. autoclass:: biorefineries.cellulosic.units.FeedStockHandling

Tanks
-----
.. autoclass:: biorefineries.cellulosic.units.AmmoniaStorageTank
.. autoclass:: biorefineries.cellulosic.units.AmmoniaReacidificationTank
.. autoclass:: biorefineries.cellulosic.units.BeerTank
.. autoclass:: biorefineries.cellulosic.units.CSLStorageTank
.. autoclass:: biorefineries.cellulosic.units.DAPStorageTank
.. autoclass:: biorefineries.cellulosic.units.AmmoniaAdditionTank
.. autoclass:: biorefineries.cellulosic.units.AmmoniaMixer
.. autoclass:: biorefineries.cellulosic.units.EnzymeHydrolysateMixer
.. autoclass:: biorefineries.cellulosic.units.SulfuricAcidMixer
.. autoclass:: biorefineries.cellulosic.units.SeedHoldTank
.. autoclass:: biorefineries.cellulosic.units.SulfuricAcidStorageTank
.. autoclass:: biorefineries.cellulosic.units.SulfuricAcidTank
.. autoclass:: biorefineries.cellulosic.units.OligomerConversionTank

Heat exchange
-------------
.. autoclass:: biorefineries.cellulosic.units.HydrolysateHeatExchanger
.. autoclass:: biorefineries.cellulosic.units.PretreatmentWasteHeater
.. autoclass:: biorefineries.cellulosic.units.WasteVaporCondenser

Pumps
-----
.. autoclass:: biorefineries.cellulosic.units.HydrolysatePump
.. autoclass:: biorefineries.cellulosic.units.BlowdownDischargePump
.. autoclass:: biorefineries.cellulosic.units.HydrolyzatePump
.. autoclass:: biorefineries.cellulosic.units.ReacidifiedHydrolyzatePump

"""
from biosteam.units.design_tools.geometry import cylinder_diameter_from_volume
from thermosteam import MultiStream
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
import thermosteam as tmo
import biosteam as bst
import numpy as np

__all__ = (
    'PretreatmentReactorSystem',
    'SeedTrain',
    'ContinuousPresaccharification',
    'Saccharification',
    'CoFermentation',
    'SaccharificationAndCoFermentation',
    'SimultaneousSaccharificationAndCoFermentation',
    'ReverseOsmosis',
    'Nanofilter',
    'HydrolysatePump',
    'AmmoniaStorageTank',
    'AmmoniaReacidificationTank',
    'BeerTank',
    'BlowdownDischargePump',
    'CSLStorageTank',
    'DAPStorageTank',
    'FeedStockHandling',
    'PretreatmentFlash',
    'HydrolysateHeatExchanger',
    'PretreatmentWasteHeater',
    'WasteVaporCondenser',
    'HydrolyzatePump',
    'HydrolyzateSolidLiquidSeparator',
    'AmmoniaAdditionTank',
    'AmmoniaMixer',
    'EnzymeHydrolysateMixer',
    'SulfuricAcidMixer',
    'OligomerConversionTank',
    'ReacidifiedHydrolyzatePump',
    'SeedHoldTank',
    'SulfuricAcidStorageTank',
    'SulfuricAcidTank',
)

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# %% Constants

_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% Pretreatment

@cost('Dry flow rate', 'Pretreatment reactor system', units='kg/hr',
      S=83333, CE=522, cost=19812400 * 0.993, n=0.6, kW=4578, BM=1.5)
class PretreatmentReactorSystem(bst.units.design_tools.PressureVessel, Unit):
    _N_ins = 1
    _N_outs = 2
    _graphics = bst.Flash._graphics
    _units = {'Residence time': 'hr',
              'Reactor volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), T=130+273.15, thermo=None, 
                 tau=0.166, V_wf=0.8, length_to_diameter=2, 
                 vessel_material='Stainless steel 316', 
                 vessel_type='Horizontal',
                 reactions=None,
                 run_vle=True):
        Unit.__init__(self, ID, ins, outs, thermo)
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
        Rxn('Galactan + H2O -> GalactoseOligomer',       'Galactan', 0.0240, chemicals),
        Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030, chemicals),
        Rxn('Mannan + H2O -> MannoseOligomer',           'Mannan',   0.0030, chemicals),
        Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030, chemicals),
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
            self.glucan_to_glucose = self.reactions[0]
            self.xylan_to_xylose = self.reactions[8]
            self.glucose_to_byproducts = self.reactions[1:3]
            self.xylose_to_byproducts = self.reactions[9:12]
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
    _N_ins = 1
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
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: wt % media (e.g. corn steep liquor) in each stage 
    # media_loading = 0.50
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 reactions=None, saccharification=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        chemicals = self.chemicals
        if reactions is None:
            self.reactions = ParallelRxn([
        #   Reaction definition                             Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9000, chemicals),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                            'Glucose',   0.0400, chemicals),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040, chemicals),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060, chemicals),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8000, chemicals),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                            'Xylose',    0.0400, chemicals),
        Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030, chemicals),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460, chemicals),
        Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090, chemicals)
            ])
            self.glucose_to_ethanol = self.reactions[0]
            self.xylose_to_ethanol = self.reactions[4]
            self.glucose_to_byproducts = self.reactions[1:4]
            self.xylose_to_byproducts = self.reactions[5:]
        else:
            self.reactions = reactions
        if callable(saccharification):
            self.saccharification = saccharification
        elif saccharification:
            self.saccharification = ParallelRxn([
                Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400, chemicals),
                Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120, chemicals),
                Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000, chemicals),
                Rxn('Cellobiose + H2O -> 2Glucose',       'Cellobiose',  1.0000, chemicals)]
            )
        else:
            self.saccharification = None
    
    _setup = Unit._setup
    
    def _run(self):
        vent, effluent= self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        if self.saccharification:
            self.saccharification(effluent)
        self.reactions.force_reaction(effluent)
        effluent.empty_negative_flows()
        effluent.T = self.T
        vent.phase = 'g'
        vent.copy_flow(effluent, ('CO2', 'O2'), remove=True)

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
        

# %% Saccharification and fermentation (consolidated bioprocess)

@cost('Flow rate', 'Transfer pumps', kW=58, S=352*_gpm2m3hr,
      cost=47200/5, CE=522, n=0.8, BM=2.3, N='N_transfer_pumps')
@cost('Tank volume', 'Tanks', cost=3840e3/8, S=250e3*_gal2m3, 
      CE=522, n=0.7, BM=2.0, N='N_tanks')
class ContinuousPresaccharification(Unit):
    _N_ins = 1
    _N_outs = 1
    
    #: Residence time of countinuous saccharification tanks (hr)
    tau_tank = 24
    
    #: Number of continuous saccharification tanks
    N_tanks = 8
    
    #: Number of transfer pumps
    N_transfer_pumps = 5
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    _units = {'Flow rate': 'm3/hr',
              'Tank volume': 'm3'}
    
    def __init__(self, ID='', ins=None, outs=(), P=101325, reactions=None):
        Unit.__init__(self, ID, ins, outs)
        self.P = P
        self.reactions = reactions 
        
    def _run(self):
        inlet = self.ins[0]
        outlet = self.outs[0]
        outlet.copy_flow(inlet)
        outlet.T = inlet.T
        outlet.P = self.P
        if self.reactions: self.reactions.adiabatic_reaction(outlet)
        
    def _design(self):
        inlet = self.ins[0]
        v_0 = inlet.F_vol
        Design = self.design_results
        Design['Tank volume'] = v_0 * self.tau_tank / self.V_wf / self.N_tanks
        Design['Flow rate'] = v_0 / self.N_transfer_pumps


class Saccharification(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
        
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    _units = {'Flow rate': 'm3/hr',
              'Reactor volume': 'm3',
              'Reactor duty': 'kJ/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 tau=72, N=None, V=3785.4118, T=48+273.15, P=101325,
                 Nmin=2, Nmax=36, reactions=None):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau, N, V, T, P, Nmin, Nmax)
        chemicals = self.chemicals
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from 
        #: downstream batch tank in co-fermentation.
        self.reactions = reactions or ParallelRxn([
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400, chemicals),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120, chemicals),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000, chemicals),
            Rxn('Cellobiose + H2O -> 2Glucose',       'Cellobiose',  1.0000, chemicals)]
        )
        
    _setup = Unit._setup
        
    @property
    def saccharification(self):
        return self.reactions
        
    @property
    def vent(self):
        return None
    @property
    def effluent(self):
        return self.outs[0]
        
    def _run(self):
        feed, = self.ins
        effluent, = self.outs
        effluent.copy_flow(feed)
        effluent.T = self.T
        effluent.P = self.P
        self.reactions(effluent)

class CoFermentation(bst.BatchBioreactor):
    _N_ins = 1
    _ins_size_is_fixed = False
        
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 tau=36, N=None, V=3785.4118, T=305.15, P=101325,
                 Nmin=2, Nmax=36, cofermentation=None, loss=None):
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo, tau, N, V, T, P, Nmin, Nmax)
        self.P = P
        chemicals = self.chemicals
        self.loss = loss or ParallelRxn([
    #   Reaction definition               Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.0300, chemicals),
    Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.0300, chemicals),
    Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.0300, chemicals),
    Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.0300, chemicals),
    Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.0300, chemicals),
        ])
        if cofermentation is None:
            self.cofermentation = ParallelRxn([
        #   Reaction definition                                          Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.9500, chemicals),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.0200, chemicals),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.0040, chemicals),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.0060, chemicals),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.8500, chemicals),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                                        'Xylose',    0.0190, chemicals),
        Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.0030, chemicals),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.0460, chemicals),
        Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.0090, chemicals),
            ])
            self.glucose_to_ethanol = self.cofermentation[0]
            self.xylose_to_ethanol = self.cofermentation[4]
            self.glucose_to_byproducts = self.cofermentation[1:4]
            self.xylose_to_byproducts = self.cofermentation[5:]
        else:
            self.cofermentation = cofermentation
        
        if 'CSL' in self.chemicals:
            self.CSL_to_constituents = Rxn(
                'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
            )
            self.CSL_to_constituents.basis = 'mol'
        else:
            self.CSL_to_constituents = None
        
        if all([i in self.chemicals for i in ('FFA', 'DAG', 'TAG', 'Glycerol')]):
            self.oil_reaction = self.lipid_reaction = ParallelRxn([
                Rxn('TAG + 3 Water -> 3FFA + Glycerol', 'TAG', 0.23, chemicals),
                Rxn('TAG + Water -> FFA + DAG', 'TAG', 0.02, chemicals)
            ])
        else:
            self.oil_reaction = self.lipid_reaction = None
        
    def _run(self):
        feeds = self.ins
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T
        vent.phase = 'g'
        effluent.mix_from(feeds, energy_balance=False)
        if self.loss: self.loss(effluent)
        if self.lipid_reaction: 
            self.lipid_reaction.force_reaction(effluent)
            effluent.empty_negative_flows()
        self.cofermentation(effluent)
        if self.CSL_to_constituents: self.CSL_to_constituents(effluent)
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)
        

@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=-5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Batch duty', 'Fermentor batch cooler', CE=522, cost=86928,
      S=-5*_Gcal2kJ, n=0.7, BM=1.8) # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
class SaccharificationAndCoFermentation(Unit):
    _N_ins = 1
    _N_outs = 3
    _ins_size_is_fixed = False
        
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
    #: Saccharification time (hr)
    tau_saccharification = 60
    
    #: Co-Fermentation time (hr)
    tau_cofermentation = 36
    
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    #: Working volume fraction (filled tank to total tank volume)
    V_wf = 0.9
    
    #: Number of reactors
    N_reactors = 12
    
    _units = {'Flow rate': 'm3/hr',
              'Reactor volume': 'm3',
              'Batch duty': 'kJ/hr',
              'Reactor duty': 'kJ/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, P=101325, saccharification_split=0.1,
                 saccharification=None, loss=None, cofermentation=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P
        self.saccharification_split = saccharification_split
        chemicals = self.chemicals
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from 
        #: downstream batch tank in co-fermentation.
        self.saccharification = saccharification or ParallelRxn([
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400, chemicals),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120, chemicals),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000, chemicals),
            Rxn('Cellobiose + H2O -> 2Glucose',       'Cellobiose',  1.0000, chemicals)]
        )
        self.loss = loss or ParallelRxn([
    #   Reaction definition               Reactant    Conversion
    Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.0300, chemicals),
    Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.0300, chemicals),
    Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.0300, chemicals),
    Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.0300, chemicals),
    Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.0300, chemicals),])
        self.cofermentation = cofermentation or ParallelRxn([
    #   Reaction definition                                          Reactant    Conversion
    Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.9500, chemicals),
    Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.0200, chemicals),
    Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.0040, chemicals),
    Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.0060, chemicals),
    Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.8500, chemicals),
    Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                                    'Xylose',    0.0190, chemicals),
    Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.0030, chemicals),
    Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.0460, chemicals),
    Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.0090, chemicals),
    ])
        self.CSL_to_constituents = Rxn(
            'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, chemicals, basis='wt',
        )
        self.CSL_to_constituents.basis = 'mol'
        
    def _run(self):
        feed, *other = self.ins
        vent, effluent, sidedraw = self.outs
        effluent.copy_like(feed)
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        self.saccharification(effluent)
        self._batch_duty = effluent.Hnet - feed.Hnet
        sidedraw.copy_like(effluent)
        sidedraw.mol *= self.saccharification_split
        effluent.mol -= sidedraw.mol
        effluent.mix_from([effluent, *other], energy_balance=False)
        self.loss(effluent)
        self.cofermentation(effluent)
        self.CSL_to_constituents(effluent)
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Flow rate'] = v_0 / self.N_reactors
        tau = self.tau_saccharification + self.tau_cofermentation
        Design.update(size_batch(v_0, tau, self.tau_0, self.N_reactors, self.V_wf))
        Design['Batch duty'] = batch_duty = self._batch_duty
        Design['Reactor duty'] = reactor_duty = self.Hnet - batch_duty
        self.add_heat_utility(reactor_duty + batch_duty, effluent.T)
   


@cost('Flow rate', 'Recirculation pumps', kW=30, S=340*_gpm2m3hr,
      cost=47200, n=0.8, BM=2.3, CE=522, N='N_reactors')
@cost('Reactor duty', 'Heat exchangers', CE=522, cost=23900,
      S=-5*_Gcal2kJ, n=0.7, BM=2.2, N='N_reactors') # Based on a similar heat exchanger
@cost('Reactor volume', 'Agitators', CE=522, cost=52500,
      S=1e6*_gal2m3, n=0.5, kW=90, BM=1.5, N='N_reactors')
@cost('Reactor volume', 'Reactors', CE=522, cost=844000,
      S=1e6*_gal2m3, n=0.5, BM=1.5, N='N_reactors')
class SimultaneousSaccharificationAndCoFermentation(Unit):
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
        
    #: Saccharification temperature (K)
    T_saccharification = 48+273.15
    
    #: Fermentation temperature (K)
    T_fermentation = 32+273.15
    
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
              'Batch duty': 'kJ/hr',
              'Reactor duty': 'kJ/hr'}
    
    def __init__(self, ID='', ins=None, outs=(), P=101325, thermo=None,
                 saccharification=None, loss=None, cofermentation=None):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P = P
        chemicals = self.chemicals
        #: [ParallelReaction] Enzymatic hydrolysis reactions including from 
        #: downstream batch tank in co-fermentation.
        self.saccharification = saccharification or ParallelRxn([
            Rxn('Glucan -> GlucoseOligomer',          'Glucan',   0.0400, chemicals),
            Rxn('Glucan + 0.5 H2O -> 0.5 Cellobiose', 'Glucan',   0.0120, chemicals),
            Rxn('Glucan + H2O -> Glucose',            'Glucan',   0.9000, chemicals),
            Rxn('Cellobiose + H2O -> 2Glucose',       'Cellobiose',  1.0000, chemicals)
        ])
        self.loss = loss or ParallelRxn([
        #   Reaction definition               Reactant    Conversion
        Rxn('Glucose -> 2 LacticAcid',       'Glucose',   0.0300, chemicals),
        Rxn('3 Xylose -> 5 LacticAcid',      'Xylose',    0.0300, chemicals),
        Rxn('3 Arabinose -> 5 LacticAcid',   'Arabinose', 0.0300, chemicals),
        Rxn('Galactose -> 2 LacticAcid',     'Galactose', 0.0300, chemicals),
        Rxn('Mannose -> 2 LacticAcid',       'Mannose',   0.0300, chemicals),])
        self.cofermentation = cofermentation or ParallelRxn([
        #   Reaction definition                                          Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                             'Glucose',   0.9500, chemicals),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O', 'Glucose',   0.0200, chemicals),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',                       'Glucose',   0.0040, chemicals),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',                   'Glucose',   0.0060, chemicals),
        Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                            'Xylose',    0.8500, chemicals),
        Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
                                                                        'Xylose',    0.0190, chemicals),
        Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',                  'Xylose',    0.0030, chemicals),
        Rxn('Xylose + H2O -> Xylitol + 0.5 O2',                         'Xylose',    0.0460, chemicals),
        Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',              'Xylose',    0.0090, chemicals),
        ])
        self.CSL_to_constituents = Rxn(
            'CSL -> 0.5 H2O + 0.25 LacticAcid + 0.25 Protein', 'CSL', 1.0000, basis='wt',
        )
        self.CSL_to_constituents.basis = 'mol'
        
    def _run(self):
        vent, effluent = self.outs
        vent.P = effluent.P = self.P
        vent.T = effluent.T = self.T_fermentation
        vent.phase = 'g'
        effluent.mix_from(self.ins, energy_balance=False)
        self.saccharification(effluent)
        self.loss(effluent)
        self.cofermentation(effluent)
        self.CSL_to_constituents(effluent)
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)

    def _design(self):
        effluent = self.outs[1]
        v_0 = effluent.F_vol
        Design = self.design_results
        Design['Flow rate'] = v_0 / self.N_reactors
        Design.update(size_batch(v_0, self.tau, self.tau_0, self.N_reactors, self.V_wf))
        Design['Reactor duty'] = reactor_duty = self.Hnet
        self.add_heat_utility(reactor_duty, effluent.T)


# %% Pretreatment separations 

# Membrane separation processes. Perry's Chemical Engineer's Handbook 7th Edition. 
@cost('Flow rate', 'Nanofiltration', kW=18000, S=0.25, units='m3/s',
      cost=17300, n=0.6, BM=1., CE=bst.units.design_tools.CEPCI_by_year[1996])
@cost('Annual flow rate', 'Membranes and maintenance', S=1., units='m3/yr',
      cost=0.13, n=1., BM=1., CE=bst.units.design_tools.CEPCI_by_year[1996],
      annual=True)
class ReverseOsmosis(bst.SolidsSeparator):
    pass

Europe_investment_site_factor = 1.2
euro_to_dollar = 1.04
installation_cost = 115000 # euro
volume_treated = 10 * 24 * 30 * 18 # m3
membrane_area = 27.5 # m2
membrane_cost = 95 # euro / m2
cleaning_cost = 50 # euro / m2
maintenance_cost = installation_cost * 0.02
yearly_operating_hours = 8000
operating_cost = membrane_area * (membrane_cost + cleaning_cost) + maintenance_cost * yearly_operating_hours / (18 * 30 * 24)
operating_cost_per_volume_treated = operating_cost / volume_treated
US_operating_cost = euro_to_dollar * operating_cost_per_volume_treated / Europe_investment_site_factor
US_installation_cost = euro_to_dollar * installation_cost / Europe_investment_site_factor
CEPCI2013 = bst.units.design_tools.CEPCI_by_year[2013]
electricity_cost = 2e3 # euro / yr
electricity_price = 0.04 # euro / kWh
electricity_demand = electricity_cost / yearly_operating_hours / electricity_price

@cost('Flow rate', 'Nanofiltration', kW=electricity_demand, S=10, units='m3/hr',
      cost=US_installation_cost, n=0.6, BM=1., CE=CEPCI2013)
@cost('Annual flow rate', 'Membranes and maintenance', S=1., units='m3/yr',
      cost=US_operating_cost, n=1., BM=1., CE=CEPCI2013,
      annual=True)
class Nanofilter(bst.Unit):
    _N_ins = 1
    _N_outs = 2
    volume_reduction = 0.80
    lignin_retention = 0.60
    NaOH_retention = 0.10
    sugar_retention = 0.95
    sugars_IDs = ('Glucose', 'Xylose', 'Arabinose', 'Galactose', 'Mannose',
              'GlucoseOligomer', 'GalactoseOligomer', 'MannoseOligomer',
              'XyloseOligomer', 'ArabinoseOligomer', 'HMF', 'Furfural', 'Xylobiose')
    lignin_IDs = ('SolubleLignin',)
    retentate_only_IDs = ('Lignin', 'Glucan', 'Xylan', 'Arabinan', 'Galactan',
                          'Mannan', 'Solids', 'Ash')
    
    def _run(self):
        feed, = self.ins
        permeate, retentate = self.outs
        permeate.empty()
        retentate.empty()
        lignin_retention = self.lignin_retention
        sugar_retention = self.sugar_retention
        NaOH_retention = self.NaOH_retention
        volume_reduction = self.volume_reduction
        water = feed.imass['Water']
        permeate.imass['Water'] = volume_reduction * water
        retentate.imass['Water'] = feed.imass['Water'] - permeate.imass['Water']
        retentate.copy_flow(feed, self.retentate_only_IDs)
        not_permeate_unique = (*self.retentate_only_IDs, *self.lignin_IDs, *self.sugars_IDs, 'NaOH', 'Water')
        not_permeate_unique = set([self.chemicals[i].ID for i in not_permeate_unique])
        permeate_unique = tuple([i.ID for i in self.chemicals if i.ID not in not_permeate_unique])
        permeate.copy_flow(feed, permeate_unique)
        F_mass = feed.F_mass
        lignin_comps = feed.imass[self.lignin_IDs]
        sugar_comps = feed.imass[self.sugars_IDs]
        NaOH = feed.imass['NaOH']
        lignin = lignin_comps.sum()
        sugar = sugar_comps.sum()
        lignin_def = lignin_comps / lignin
        sugar_def = sugar_comps / sugar
        NaOH = feed.imass['NaOH']
        feed_mass = np.array([lignin, sugar, NaOH])
        z_mass = feed_mass / F_mass
        K = np.array([lignin_retention, sugar_retention, NaOH_retention])
        phi_guess = 0.5
        Fa = permeate.F_mass
        Fb = retentate.F_mass
        assert abs((Fa + Fb + lignin + sugar + NaOH) - F_mass) < 1e-6
        permeate_mass = feed_mass / water * (1 - K) * permeate.imass['Water']
        
        # phi = tmo.separations.compute_phase_fraction(z_mass, K, phi_guess, Fa/F_mass, Fb/F_mass)
        # x_retentate = z_mass / (phi * K + (1. - phi))
        # retentate_mass = x_retentate * (1. - phi) * F_mass
        
        permeate.imass[self.lignin_IDs] = permeate_mass[0] * lignin_def
        permeate.imass[self.sugars_IDs] = permeate_mass[1] * sugar_def
        permeate.imass['NaOH'] = permeate_mass[2]
        retentate.mol = feed.mol - permeate.mol
        assert (permeate.mol >= 0.).all()
        retentate.T = permeate.T = feed.T


# %% Simple unit operations

@cost('Flow rate', 'Pump', units='kg/hr',
      S=402194, CE=522, cost=22500, n=0.8, kW=74.57, BM=2.3)
class HydrolysatePump(Unit): pass

@cost('Flow rate', 'Tank', S=1171, units='kg/hr',
      CE=522, cost=196000, n=0.7, BM=2)
class AmmoniaStorageTank(Unit): pass

@cost('Flow rate', 'Agitator', S=410846, units='kg/hr',
      CE=522, cost=None, n=0.5, BM=1.5, kW=14.914)
@cost('Flow rate', 'Tank', S=410369, units='kg/hr',
      CE=522, cost=None, n=0.7, BM=2)
class AmmoniaReacidificationTank(Unit): pass

@cost('Flow rate', 'Tank', S=425878, units='kg/hr',
      CE=522, cost=636000, n=0.7, BM=1.8)
@cost('Flow rate', 'Pump', S=425878, units='kg/hr',
      CE=522, cost=26800, n=0.8, BM=2.3, kW=93.2125)
@cost('Flow rate', 'Agitator', S=425878, units='kg/hr',
      CE=522, cost=68300, n=0.5, BM=1.5, kW=14.914)
class BeerTank(Unit): pass

@cost('Flow rate', 'Pump', S=292407, units='kg/hr',
      CE=551, cost=25365, n=0.8, BM=2.3, kW=93.2125)
class BlowdownDischargePump(Unit): pass

@cost('Flow rate', 'Tank', S=1393, units='kg/hr',
      CE=522, cost=70000, n=0.7, BM=2.6)
@cost('Flow rate', 'Pump', S=1393, units='kg/hr',
      CE=522, cost=3000, n=0.8, BM=3.1, kW=0.37285)
@cost('Flow rate', 'Agitator', S=1393, units='kg/hr',
      CE=522, cost=21200, n=0.5, BM=1.5, kW=7.457)
class CSLStorageTank(Unit): pass

@cost('Flow rate', 'Tank', S=163, units='kg/hr',
      CE=522, cost=102000, n=0.7, BM=1.8)
@cost('Flow rate', 'Pump', S=163, units='kg/hr',
      CE=522, cost=3000, n=0.8, BM=3.1, kW=0.37735)
@cost('Flow rate', 'Agitator', S=163, units='kg/hr',
      CE=522, cost=9800, n=0.5, BM=1.5, kW=4.10135)
@cost('Flow rate', 'Bag unloader', S=163, units='kg/hr',
      CE=522, cost=30000, n=0.6, BM=1.7)
class DAPStorageTank(Unit): pass

@cost('Flow rate', 'System', S=94697, units='kg/hr',
      CE=522, cost=13329690, n=0.6, BM=1.7, kW=783)
class FeedStockHandling(Unit): pass

@cost('Flow rate', 'Tank', S=252891, units='kg/hr',
      CE=522, cost=511000, n=0.7, BM=2)
@cost('Flow rate', 'Pump', S=252891, units='kg/hr',
      CE=522, cost=30000, n=0.8, BM=1.7, kW=55.9)
@cost('Flow rate', 'Agitator', S=252891, units='kg/hr',
      CE=522, cost=90000, n=0.5, BM=1.5, kW=170)
class PretreatmentFlash(bst.Flash): pass

@cost('Duty', 'Heat exchanger', S=8, units='Gcal/hr',
      CE=522, cost=85000, n=0.7, BM=2.2)
class HydrolysateHeatExchanger(bst.HXutility): pass

@cost('Duty', 'Heat Exchanger', S=8, units='Gcal/hr',
      CE=551, cost=92000, n=0.7, BM=2.2)
class PretreatmentWasteHeater(bst.HXutility): pass

@cost('Duty', 'Heat Exchanger', S=2, units='Gcal/hr',
      CE=522, cost=34000, n=0.7, BM=2.2)
class WasteVaporCondenser(bst.HXutility): pass

@cost('Flow rate', 'Pump', S=402194, units='kg/hr',
      CE=522, cost=22500, n=0.8, BM=2.3, kW=74.57)
class HydrolyzatePump(Unit): pass

@cost('Flow rate', 'Separator', S=39000, units='kg/hr',
      CE=522, cost=35000000, n=0.7, BM=1.7)
class HydrolyzateSolidLiquidSeparator(Unit): pass

@cost('Flow rate', 'Tank', S=410369, units='kg/hr',
      CE=522, cost=236000, n=0.7, BM=2)
@cost('Flow rate', 'Agitator', S=410369, units='kg/hr',
      CE=522, cost=21900, n=0.5, BM=1.5, kW=7.457)
class AmmoniaAdditionTank(bst.MixTank): pass

@cost('Flow rate', 'Mixer', S=157478, units='kg/hr',
      CE=522, cost=5000, n=0.5, BM=1)
class AmmoniaMixer(bst.Mixer): pass

@cost('Flow rate', 'Mixer', S=380000, units='kg/hr',
      CE=522, cost=109000, n=0.5, BM=1.7, kW=74.57)
class EnzymeHydrolysateMixer(bst.Mixer): pass

@cost('Flow rate', 'Mixer', S=136260, units='kg/hr',
      CE=522, cost=6000, n=0.5, BM=1)
class SulfuricAcidMixer(bst.Mixer): pass

@cost('Flow rate', 'Tank', S=264116, units='kg/hr',
      CE=522, cost=203000, n=0.7, BM=2)
@cost('Flow rate', 'Pump', S=264116, units='kg/hr',
      CE=551, cost=17480, n=0.8, BM=1.7, kW=55.9)
@cost('Flow rate', 'Agitator', S=264116, units='kg/hr',
      CE=522, cost=90000, n=0.5, BM=1.5, kW=170)
class OligomerConversionTank(Unit): pass

@cost('Flow rate', 'Pump', S=402195, units='kg/hr',
      CE=522, cost=None, n=0.8, BM=2.3, kW=44.742)
class ReacidifiedHydrolyzatePump(Unit): pass

@cost('Flow rate', 'Pump', S=43149, units='kg/hr',
      CE=522, cost=8200, n=0.8, BM=2.3, kW=7.457)
@cost('Flow rate', 'Tank', S=40414, units='kg/hr',
      CE=522, cost=439000, n=0.7, BM=1.8)
@cost('Flow rate', 'Agitator', S=40414, units='kg/hr',
      CE=522, cost=31800, n=0.5, BM=1.5, kW=11.3205)
class SeedHoldTank(Unit): pass

@cost('Flow rate', 'Tank', S=1981, units='kg/hr',
      CE=522, cost=96000, n=0.7, BM=1.5)
@cost('Flow rate', 'Pump', S=1981, units='kg/hr',
      CE=522, cost=7493, n=0.8, BM=2.3, kW=0.5)
class SulfuricAcidStorageTank(Unit): pass

@cost('Flow rate', 'Tank', S=1981, units='kg/hr',
      CE=551, cost=6210, n=0.7, BM=2)
@cost('Flow rate', 'Pump', S=3720, units='kg/hr',
      CE=522, cost=8000, n=0.8, BM=2.3, kW=1)
class SulfuricAcidTank(Unit): pass

# -*- coding: utf-8 -*-
"""
Created on 2025-04-18 15:20:45

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import biosteam as bst
import thermosteam as tmo

from biosteam.units.decorators import cost, copy_algorithm
from biosteam.units.design_tools import CEPCI_by_year, cylinder_diameter_from_volume, cylinder_area
from biosteam import tank_factory
from thermosteam import MultiStream
from biosteam import Unit
from biosteam.units.decorators import cost
from biosteam.units.design_tools import size_batch
import thermosteam as tmo
import biosteam as bst
import numpy as np

import numpy as np
import flexsolve as flx

__all__ = (
    ##### Upstream #####    
    # Feedstock Preparation
    #'PretreatmentReactorSystem',

    ## glucose as beginning

    # Bioreactor Fermentaion
    'SeedTrain',
    'CoFermentation',
    #'PSA',

    ##### Downstream #####
    # Cell Harvesting

    # Purification
    #'Centrifuge', 'Concentration','Filtration_Agarose',
    #'Sedimentation','UltraFiltration_Desalt','SprayDring',

    #'ScrewPress','CellDisruption',
)
Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# %% Constants
_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% 
##### UpStream #####

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
    tau_batch = 16
    
    @property
    def tau_turnover(self):
        """Turnover time (hr) calculated by batch time divided by number of trains."""
        return self.tau_batch/self.N_trains
    
    #: Operating temperature (K)
    T = 32+273.15
    
    # #: Diammonium phosphate loading in g/L of fermentation broth
    # DAP = 0.67 
    
    def _init(self, reactions=None, saccharification=None, T=None):
        chemicals = self.chemicals

        if reactions is None:
            self.reactions = ParallelRxn([
        #   Reaction definition                             Reactant    Conversion
        Rxn('Glucose -> 2 Ethanol + 2 CO2',                 'Glucose',   0.9000, chemicals),
        Rxn('Glucose + 0.047 CSL + 0.018 DAP -> 6 Z_mobilis + 2.4 H2O',
                                                            'Glucose',   0.0400, chemicals),
        Rxn('Glucose + 2 H2O -> 2 Glycerol + O2',           'Glucose',   0.0040, chemicals),
        Rxn('Glucose + 2 CO2 -> 2 SuccinicAcid + O2',       'Glucose',   0.0060, chemicals),
        # Rxn('3 Xylose -> 5 Ethanol + 5 CO2',                'Xylose',    0.8000, chemicals),
        # Rxn('Xylose + 0.039 CSL + 0.015 DAP -> 5 Z_mobilis + 2 H2O',
        #                                                     'Xylose',    0.0400, chemicals),
        # Rxn('3 Xylose + 5 H2O -> 5 Glycerol + 2.5 O2',      'Xylose',    0.0030, chemicals),
        # Rxn('Xylose + H2O -> Xylitol + 0.5 O2',             'Xylose',    0.0460, chemicals),
        # Rxn('3 Xylose + 5 CO2 -> 5 SuccinicAcid + 2.5 O2',  'Xylose',    0.0090, chemicals)
            ])
            self.glucose_to_ethanol = self.reactions[0]
            # self.xylose_to_ethanol = self.reactions[4]
            self.glucose_to_byproducts = self.reactions[1:4]
            # self.xylose_to_byproducts = self.reactions[5:]
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
        if T: self.T = T
    
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

class CoFermentation(bst.BatchBioreactor):
    _N_ins = 1
    _ins_size_is_fixed = False
        
    #: Unload and clean up time (hr)
    tau_0 = 4
    
    def _init(self, tau=36, N=None, V=3785.4118, T=305.15, P=101325,
                Nmin=2, Nmax=36, cofermentation=None, loss=None):
        super()._init(tau, N, V, T, P, Nmin, Nmax)
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

# %%
class PSA(bst.Flash):
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


# class CellDisruption(bst.Homogenizer): pass


# # %%

# #@cost('Flow rate', units='kg/hr', CE=CEPCI_by_year[2010], cost=100000, S=100000, n=0.6, kW=100)
# #@copy_algorithm(bst.SolidLiquidsSplitCentrifuge, run=False)       
# class Centrifuge(bst.SpliSolidLiquidsSplitCentrifugetter): pass

# # %%
# class Evaporator(bst.MultiEffectEvaporator): pass

# # %%

# class Concentration(bst.PressureFilter): pass
        
# # class Filtration_Agarose(bst.PressureFilter): pass

# # %%

# class Sedimentation(bst.Clarifier): pass

# class UltraFiltration_Desalt(bst.PressureFilter): pass
# # %%
# class SprayDring(bst.SprayDryer): pass

# # %%

# class ScrewPress(bst.ScrewPress): pass
        
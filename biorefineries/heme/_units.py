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

from biorefineries import cellulosic

import numpy as np
import flexsolve as flx

__all__ = (
    ##### Upstream #####    
    # Feedstock Preparation
    ## glucose as beginning
    
    # Bioreactor Fermentaion
    #'SeedTrain',
    'Fermentation',
    #'PSA',

    ##### Downstream #####
    # Cell Harvesting


    # Purification
    #'Centrifuge', 'Concentration','Filtration_Agarose',
    #'Sedimentation','UltraFiltration_Desalt','SprayDring',

    #'ScrewPress','CellDisruption',
)


# %% 
##### UpStream #####

class SeedTrain(bst.SeedTrain):
    _N_ins = 2
    _N_outs = 1
    _graphics = bst.SeedTrain._graphics

    V_max_default = 500
    
    def _init(self, reactions=None, saccharification=False):
        self.saccharification = saccharification
        chemicals = self.chemicals
        self.reactions = reactions or tmo.PRxn([
    #   Reaction definition                   Reactant    Conversion
    tmo.Rxn('Glucose -> 2 Ethanol + 2 CO2',       'Glucose',   0.90, chemicals, correct_mass_balance=True),
    tmo.Rxn('Glucose -> Cellmass',                'Glucose',  0.05, chemicals, correct_mass_balance=True),
        ])

    def _setup(self):
        super()._setup()
        self.outs[0].phase = 'g'

    def _run(self):
        vent, effluent= self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions.force_reaction(effluent)
        effluent.mol.remove_negatives()
        effluent.T = self.T
        vent.empty()
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)

class Fermentation(bst.AeratedBioreactor):
    V_max_default = 500


class cFermentation(bst.AeratedBioreactor):
    V_max_default = 500
    def _init(
            self, fermentation_reaction, cell_growth_reaction, 
            dT_hx_loop=8,
            Q_O2_consumption=-460240, 
            # [kJ/kmol] equivalent to 110 kcal / mol 
            # as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
            batch=True,
            **kwargs,
        ):
        bst.AeratedBioreactor._init(self, reactions=None, batch=batch, dT_hx_loop=dT_hx_loop, 
                                    Q_O2_consumption=Q_O2_consumption,
                                    optimize_power=True, **kwargs)
        chemicals = self.chemicals
    
        self.fermentation_reaction = fermentation_reaction
        self.cell_growth_reaction = cell_growth_reaction

        self.LegH_reaction = bst.SRxn([
            bst.PRxn([
                bst.Rxn('Glucose + FeSO4 + NH3 + O2 -> Heme_b + H2O + H2SO4', 'TAG', 0.90, chemicals,correct_mass_balance=True),
                bst.Rxn('Glucose +  (NH4)2SO4 + O2 -> Protein + H2O + H2SO4', 'TAG', 0.90, chemicals,correct_mass_balance=True)]),
            bst.Rxn('Heme_b + Protein -> Leghemoglobin', 'TAG', 0.95, chemicals,correct_mass_balance=True)
        ])
    
    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        assert not effluent.imol['CO2', 'O2', 'N2'].any()
    
    def _run_reactions(self, effluent):
        self.LegH_reaction.force_reaction(effluent)
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.fermentation_reaction.force_reaction(effluent)
        self.cell_growth_reaction.force_reaction(effluent) 

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
        
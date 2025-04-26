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
    'SeedTrain',
    'Fermentation',

    'Flash',

    'Centrifuge', 'Concentration','Filtration_Agarose',
    'Sedimentation','UltraFiltration_Desalt','SprayDring',

    'ScrewPress','CellDisruption',
)


# %% 
##### UpStream #####

class SeedTrain(bst.SeedTrain):
    _N_ins = 2
    _N_outs = 1
    _graphics = bst.SeedTrain._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.SeedTrain.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.SeedTrain._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.SeedTrain._cost(self)
        

# class SeedTrain(cellulosic.SeedTrain): pass


# class SeedTrain(Unit):

class Fermentation(bst.AeratedBioReactor): pass

# %%
class Flash(bst.Flash): 
    N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 1
    _graphics = bst.Flash._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 P=101325, Q=0, Q_in=False,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.Flash.__init__(self, ID, ins, outs, thermo, P=P, Q=Q, Q_in=Q_in)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.Flash._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.Flash._cost(self)


class CellDisruption(bst.Homogenizer): pass


# %%

@cost('Flow rate', units='kg/hr', CE=CEPCI_by_year[2010], cost=100000, S=100000, n=0.6, kW=100)
@copy_algorithm(bst.SolidLiquidsSplitCentrifuge, run=False)
class Centrifuge(bst.Splitter):
    _graphics = bst.SolidLiquidsSplitCentrifuge._graphics
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 split=None,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.Splitter.__init__(self, ID, ins, outs, thermo, split=split)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.Splitter._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.Splitter._cost(self)
        
# class Centrifuge(bst.SpliSolidLiquidsSplitCentrifugetter): pass

# %%
class Evaporator(bst.MultiEffectEvaporator):
    _N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 0
    _graphics = bst.MultiEffectEvaporator._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 P = (101325, 73580, 50892, 32777),
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.MultiEffectEvaporator.__init__(self, ID, ins, outs, thermo, P=P)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.MultiEffectEvaporator._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.MultiEffectEvaporator._cost(self)

# %%

class Concentration(bst.PressureFilter):
    _N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 0
    _graphics = bst.PressureFilter._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.PressureFilter.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.PressureFilter._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.PressureFilter._cost(self)
        
# class Filtration_Agarose(bst.PressureFilter): pass

# %%

class Sedimentation(bst.Clarifier): pass

class UltraFiltration_Desalt(bst.PressureFilter): pass
# %%
class SprayDring(bst.SprayDryer): 
    N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 1
    _graphics = bst.SprayDryer._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.SprayDryer.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.SprayDryer._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.SprayDryer._cost(self)

# %%

class ScrewPress(bst.ScrewPress):     
    N_ins = 1
    _N_outs = 2
    _N_heat_utilities = 0
    _graphics = bst.ScrewPress._graphics
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_material='Carbon steel', vessel_type='Vertical',
                 wall_thickness_factor=1):
        bst.ScrewPress.__init__(self, ID, ins, outs, thermo)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.wall_thickness_factor = wall_thickness_factor
    
    def _design(self):
        bst.ScrewPress._design(self)
        self.design_results['Wall thickness'] *= self.wall_thickness_factor
        self.design_results['Weight'] *= self.wall_thickness_factor
    
    def _cost(self):
        bst.ScrewPress._cost(self)
        
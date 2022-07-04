# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
import thermosteam as tmo
from . import units

__all__ = ('FattyalcoholProcessSpecifications',)

ReactorSpecification = bst.process_tools.ReactorSpecification

# %% Constant specifications

fermentation_products = ('Hexanol','Octanol', 'Decanol', 
                         'Dodecanol', 'Tetradecanol')

substrates = ('Glucose',)

# %% Process specification

class FattyalcoholProcessSpecifications:
    """
    Create an FattyalcoholProcessSpecifications object that handles overall
    process and technology specifications of the fattyalcohol production system.
    
    """
    __slots__ = ('operating_days',
                 'plant_size',
                 'titer',
                 'productivity',
                 'yield_',
                 'reaction_name',
                 'unit',
    )
    
    products = fermentation_products
    substrates = substrates
    
    def __init__(self, system, tea,
                 plant_size=150000,
                 operating_days=350.4, 
                 titer=80,
                 productivity=1,
                 yield_=0.5):
        self.plant_size = plant_size
        self.operating_days = operating_days
        self.titer = titer
        self.productivity = productivity
        self.yield_ = yield_
        self.reaction_name = 'fermentation_reaction'
        self.unit = tmo.utils.get_instance(system.units, units.FattyAlcoholBioreactor)
        
    def load_specifications(self, plant_size=None, operating_days=None,
                            titer=None, productivity=None, yield_=None):
        self.load_yield(yield_ or self.yield_)
        self.load_operating_days(operating_days or self.operating_days)
        self.load_productivity(productivity or self.productivity)
        self.load_titer(titer or self.titer)
        self.load_plant_size(plant_size or self.plant_size)
        
    def load_plant_size(self, plant_size):
        production = plant_size * 907.18474 / 24. / self.operating_days
        F_mass_substrates = production / self.products_over_substrates
        feed = self.feed
        z_mass_substrates = feed.get_mass_composition(self.substrates).sum()
        feed.F_mass = F_mass_substrates / z_mass_substrates
        self.plant_size = plant_size
    
    def load_operating_days(self, operating_days):
        self.operating_days = self.tea.operating_days = operating_days
    
    load_titer = ReactorSpecification.load_titer
    load_yield = ReactorSpecification.load_yield
    load_productivity = ReactorSpecification.load_productivity
    _calculate_titer = ReactorSpecification._calculate_titer
    _titer_objective_function = ReactorSpecification._titer_objective_function
    evaluate_across_productivity = ReactorSpecification.evaluate_across_productivity
    evaluate_across_productivity = ReactorSpecification.evaluate_across_productivity
    evaluate_across_TRY = ReactorSpecification.evaluate_across_TRY
    feed = ReactorSpecification.feed
    vent = ReactorSpecification.vent
    effluent = ReactorSpecification.effluent
    products_over_substrates = ReactorSpecification.products_over_substrates
    reaction = ReactorSpecification.reaction
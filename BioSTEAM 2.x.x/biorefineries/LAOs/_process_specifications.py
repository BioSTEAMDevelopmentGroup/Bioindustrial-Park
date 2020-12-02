# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
import biosteam as bst
from thermosteam.utils import get_instance
from . import units as units

__all__ = ('LAOsProcessSpecifications',
           'fermentation_products',
           'load_process_settings')

# %% Constant specifications

fermentation_products = ('Hexanol','Octanol', 'Decanol', 
                         'Dodecanol', 'Tetradecanol')

# Not used; all steam produced on-site
# prices = [0.2378, 0.2756, 0.3171]
# Ts = [412.189, 454.77, 508.991]
# coef = np.polyfit(Ts, prices, 1)
# calculate_steam_price_at_T = np.poly1d(coef)

# %% Overall process specifications

def load_process_settings():
    bst.process_tools.default_utilities()
    
    # CEPCI
    bst.CE = 567.5 # 2017
    
    # High pressure steam conditions is a process specification that
    # depends on dehydration reactor temperature
    high_pressure_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    high_pressure_steam.heat_transfer_efficiency = 0.95
    high_pressure_steam.regeneration_price = 0.
    
    # About 10 degC higher than boiling point of tridecane
    medium_pressure_steam = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
    medium_pressure_steam.heat_transfer_efficiency = 0.95
    medium_pressure_steam.T = 510
    medium_pressure_steam.P = medium_pressure_steam.chemicals.Water.Psat(510)
    medium_pressure_steam.regeneration_price = 0.
    # medium_pressure_steam.regeneration_price = calculate_steam_price_at_T(510)
    
    # Low pressure steam, T is default
    low_pressure_steam = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    low_pressure_steam.T = 412.189
    low_pressure_steam.P = 344738.0
    low_pressure_steam.regeneration_price = 0.
    # low_pressure_steam.regeneration_price = 0.2378
    # low_pressure_steam.heat_transfer_efficiency = 0.95
    
    cooling_water = bst.HeatUtility.get_cooling_agent('cooling_water')
    chilled_water = bst.HeatUtility.get_cooling_agent('chilled_water')
    cooling_water.regeneration_price = cooling_water.heat_transfer_price = 0.
    chilled_water.regeneration_price = chilled_water.heat_transfer_price = 0.

class LAOsProcessSpecifications:
    """
    Create an LAOsProcessSpecifications object that handles overall
    process and technology specifications of the LAOs production system.
    
    """
    __slots__ = ('system',
                 'dehydration_reactor',
                 'TEA',
                 'fermentation_specification',
                 'dehydration_reactor_temperature',
                 'dehydration_reactor_conversion',
                 'LAOs_over_alcohol_arr',
                 'operating_days',
                 'plant_size',
    )
    
    def __init__(self, system, TEA,
                 plant_size=150000,
                 operating_days=350.4, 
                 fermentation_titer=40,
                 fermentation_productivity=1.0,
                 fermentation_yield=0.90,
                 dehydration_reactor_temperature=623.15,
                 dehydration_reactor_mass_fraction=0.2,
                 dehydration_reactor_conversion=0.9999,
        ):
        #: [System] Should simulate the complete process.
        self.system = system
        
        #: [TEA] Should include settings for the complete process.
        self.TEA = TEA
        
        fermentation_reactor = get_instance(system.units, units.FattyAlcoholBioreactor)
        
        #: [AdiabaticFixedbedGasReactor] Dehydration reactor.
        self.dehydration_reactor = get_instance(system.units, units.AdiabaticFixedbedGasReactor)
        
        #: [ReactorSpecification] Specifications for fementor.
        self.fermentation_specification = bst.process_tools.ReactorSpecification(
            reactor = fermentation_reactor, 
            reaction_name = 'fermentation_reaction',
            substrates = ('Glucose',),
            products = fermentation_products,
            titer = fermentation_titer, 
            productivity = fermentation_productivity,
            yield_ = fermentation_yield, 
            production = 100., # Temporary
        )
        
        MW_alcohols = np.array([i.MW for i in fermentation_reactor.chemicals[fermentation_products]]) 
        #: [1d array] Stoichiometric coefficients by weight of fatty alcohol to LAOs.
        self.LAOs_over_alcohol_arr = (MW_alcohols - 18.01528) / MW_alcohols
        
        self.dehydration_reactor_temperature = dehydration_reactor_temperature 
        self.operating_days = operating_days 
        self.plant_size = plant_size 
        
        #: [float] Mass fraction of alcohols entering dehydration reaction.
        #: Setting this parameter does nothing directly, but it is used
        #: elsewhere in the process simulation.
        self.dehydration_reactor_mass_fraction = dehydration_reactor_mass_fraction
        
        #: [float]Conversion of alcohols in dehydration reaction.
        self.dehydration_reactor_conversion = dehydration_reactor_conversion
    
    def run_specifications(self, *args, **kwargs):
        """Load specifications and simulate system."""
        self.load_specifications(*args, **kwargs)
        self.system.simulate()
      
    def load_specifications(self,
                           operating_days=None,
                           dehydration_reactor_temperature=None,
                           fermentation_titer=None,
                           fermentation_productivity=None,
                           fermentation_yield=None,
                           plant_size=None,
                           dehydration_reactor_mass_fraction=None,
                           dehydration_reactor_conversion=None):
        """Load all specifications without simulating the system."""
        self.load_operating_days(operating_days or self.operating_days)
        self.load_dehydration_reactor_temperature(dehydration_reactor_temperature or self.dehydration_reactor_temperature)
        fspec = self.fermentation_specification
        fspec.load_yield(fermentation_yield or self.fermentation_yield)
        fspec.load_titer(fermentation_titer or self.fermentation_titer)
        fspec.load_productivity(fermentation_productivity or self.fermentation_productivity)
        self.load_dehydration_reactor_conversion(dehydration_reactor_conversion or self.dehydration_reactor_conversion)
        self.load_dehydration_reactor_mass_fraction(dehydration_reactor_mass_fraction or self.dehydration_reactor_mass_fraction)
        self.load_plant_size(plant_size or self.plant_size)
        
    @property
    def LAOs_over_alcohol(self):
        """[float] LAOs per fatty alcohol ratio in g LAOs / g fatty alcohol."""
        fermentation_specificaiton = self.fermentation_specification
        reaction = fermentation_specificaiton.reaction
        if reaction.basis != 'wt': reaction = reaction.copy('wt')
        products = reaction.istoichiometry[fermentation_specificaiton.products]
        products /= products.sum()
        return (products * self.LAOs_over_alcohol_arr).sum()
    
    @property
    def fermentation_yield(self):
        return self.fermentation_specification.yield_
    @fermentation_yield.setter
    def fermentation_yield(self, yield_):
        self.fermentation_specification.yield_ = yield_
    
    @property
    def fermentation_titer(self):
        return self.fermentation_specification.titer
    @fermentation_titer.setter
    def fermentation_titer(self, titer):
        self.fermentation_specification.titer = titer
    
    @property
    def fermentation_productivity(self):
        return self.fermentation_specification.productivity
    @fermentation_productivity.setter
    def fermentation_productivity(self, productivity):
        self.fermentation_specification.productivity = productivity
    
    @property
    def dehydration_reactor_mass_fraction(self):
        """
        [float] Mass fraction of alcohols entering dehydration reaction.
        Setting this parameter does nothing directly, but it is used
        elsewhere in the process simulation.
        """
        return self.dehydration_reactor.dehydration_reactor_mass_fraction
    @dehydration_reactor_mass_fraction.setter
    def dehydration_reactor_mass_fraction(self, dehydration_reactor_mass_fraction):
        self.dehydration_reactor.dehydration_reactor_mass_fraction = dehydration_reactor_mass_fraction
    
    def load_operating_days(self, operating_days):
        """
        Load operating days specification.
        
        Parameters
        ----------
        operating_days : float
            Number of operating days per year for TEA.
            
        """
        self.operating_days = self.TEA.operating_days = operating_days
    
    def load_plant_size(self, plant_size):
        """
        Load plant size specification.
        
        plant_size : float 
            Plant size in ton product / yr of LAOs.
        
        Notes
        -----
        The net flow rate of feed to bioreactor is adjusted to achieve this
        specification.
        
        """
        self.plant_size = plant_size
        self.fermentation_specification.load_production(
            plant_size / self.LAOs_over_alcohol * 907.18474 / 24. / self.operating_days
        )
    
    def load_dehydration_reactor_conversion(self, dehydration_reactor_conversion):
        self.dehydration_reactor.reaction.X[:] = self.dehydration_reactor_conversion = dehydration_reactor_conversion
    
    def load_dehydration_reactor_temperature(self, dehydration_reactor_temperature):
        """
        Load dehydration reactor temperature specification.
        
        Parameters
        ----------
        dehydration_reactor_temperature : float
            Alcohol dehydration fixed bed reactor temperature in Kelvin.
            
        Notes
        -----
        The temperature of high pressure steam and its regeneration price
        are adjusted so that its 10 K higher that reactor operating temperature.
        
        """
        self.dehydration_reactor.T = T_reactor = dehydration_reactor_temperature
        
        # Higher pressure steam is required to heatup reactor feed
        high_pressure_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
        high_pressure_steam.T = T_steam = T_reactor + 10
        # high_pressure_steam.regeneration_price =  calculate_steam_price_at_T(T_reactor)
        high_pressure_steam.P = high_pressure_steam.chemicals.Water.Psat(T_steam)
        
    def load_dehydration_reactor_mass_fraction(self, dehydration_reactor_mass_fraction):
        """
        Load dehydration reactor mass fraction specification.

        Parameters
        ----------
        dehydration_reactor_mass_fraction : float
            Mass fraction of alcohols entering dehydration reaction.
        
        Notes
        -----
        Setting this parameter does nothing directly, but it is used
        elsewhere in the process simulation.

        """
        self.dehydration_reactor_mass_fraction = dehydration_reactor_mass_fraction
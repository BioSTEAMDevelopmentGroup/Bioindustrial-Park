# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import flexsolve as flx
import biosteam as bst
from math import sqrt, log, exp
from biosteam.utils import MockStream
import thermosteam as tmo
from biorefineries.cane.composition import (
    set_lipid_fraction as set_oil_fraction,
    set_line_composition_parameters,
    set_composition_by_line,
)
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
import numpy as np
from .feature_mockups import (
    all_parameter_mockups,
    all_metric_mockups, 
)
from .process_settings import load_process_settings
from .chemicals import create_cellulosic_oilcane_chemicals as create_chemicals
from biorefineries.cellulosic import PretreatmentReactorSystem as PRS
from .systems import (
    create_sugarcane_to_ethanol_system,
    create_sugarcane_to_sugar_and_ethanol_system,
    create_oilcane_to_biodiesel_and_ethanol_1g,
    create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    create_sugarcane_to_ethanol_combined_1_and_2g,  
    create_oilcane_to_crude_oil_and_ethanol_1g,
    create_oilcane_to_crude_oil_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    create_oilcane_to_biodiesel_1g,
    create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation,
    create_oilcane_to_biodiesel_and_actag_1g,
    create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation,
)
from .data import microbial_oil_baseline as perf
from .parse_configuration import (
    parse_configuration,
    format_configuration,
)
from biorefineries.tea import (
    create_cellulosic_ethanol_tea as create_tea
)
from .composition import (
    CaneCompositionSpecification,
)
from .oil_extraction import (
    OilExtractionSpecification,
    MockExtractionSpecification,
)
from .data.lca_characterization_factors import (
    GWP_characterization_factors,
    set_GWPCF,
    GWP,
)

__all__ = (
    'Biorefinery',
    'YRCP2023',
)

PRS_cost_item = PRS.cost_items['Pretreatment reactor system']

energycane_dct = dict(
    Water=60,
    Sucrose=9.0,
    Glucose=0.1,
    Fiber=29.3,
    Ash=2.8,
)

kg_per_ton = 907.18474
kg_per_MT = 1000
L_per_gal = 3.7854
biodiesel_kg_per_gal = 3.3111
biodiesel_gal_per_kg = 1. / biodiesel_kg_per_gal
ethanol_kg_per_gal = 2.98668849
ethanol_gal_per_kg = 1. / ethanol_kg_per_gal
biodiesel_L_per_kg = biodiesel_gal_per_kg * L_per_gal
biodiesel_kg_per_L = 1. / biodiesel_L_per_kg
ethanol_L_per_kg = ethanol_gal_per_kg * L_per_gal
ethanol_kg_per_L = 1. / ethanol_L_per_kg

screwpress_microbial_oil_recovery = frozenset([7, 9])
cellulosic_configurations = frozenset([-2, 2, 4, 6, 8, 9])
biodiesel_configurations = frozenset([1, 2, 5, 6, 7, 8, 9])
ethanol_configurations = frozenset([-2, -1, 1, 2, 3, 4])
actag_configurations = frozenset([10, 11])
conventional_ethanol_configurations = ethanol_configurations.difference(cellulosic_configurations)
cellulosic_ethanol_configurations = cellulosic_configurations.intersection(ethanol_configurations)
ethanol_biodiesel_configurations = ethanol_configurations.intersection(biodiesel_configurations)
cellulosic_ethanol_biodiesel_configurations = ethanol_biodiesel_configurations.intersection(cellulosic_ethanol_configurations)
conventional_ethanol_biodiesel_configurations = ethanol_biodiesel_configurations.difference(cellulosic_ethanol_configurations)
oil_configurations = biodiesel_configurations.difference(ethanol_configurations)
cellulosic_oil_configurations = oil_configurations.intersection(cellulosic_configurations)

system_factories = {
    -1: create_sugarcane_to_ethanol_system,
    -2: create_sugarcane_to_ethanol_combined_1_and_2g,
    -3: create_sugarcane_to_sugar_and_ethanol_system, # Not yet optimized for sugar production, not recommended to use
    1: create_oilcane_to_biodiesel_and_ethanol_1g,
    2: create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    3: create_oilcane_to_crude_oil_and_ethanol_1g,
    4: create_oilcane_to_crude_oil_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    5: create_oilcane_to_biodiesel_1g,
    6: create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation,
    7: create_oilcane_to_biodiesel_1g,
    8: create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation,
    9: create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation,
    1j: create_oilcane_to_biodiesel_and_actag_1g,
    2j: create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation,
}
system_factory_options = {
    -1: dict(use_area_convention=True, pellet_bagasse=False, dry_bagasse=True),
    5: dict(fed_batch=True),
    6: dict(fed_batch=True),
    7: dict(fed_batch=False),
    8: dict(fed_batch=False),
    9: dict(fed_batch=False, oil_extraction='screwpress'),
}
area_names = {
    -1: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'CH&P', 'Utilities', 
         'HXN', 'Storage', 'Wastewater treatment'],
    -2: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.',
         'Wastewater treatment', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    -3: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'CH&P', 'Utilities', 
         'HXN', 'Storage', 'Wastewater treatment'],
    1: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'Oil ext.', 'Biod. prod.', 
        'CH&P', 'Utilities', 'HXN', 'Storage', 'Wastewater treatment'],
    2: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.', 'Wastewater treatment',
        'Oil ext.', 'CH&P',  'Biod. prod.', 'Utilities', 'HXN', 'Storage'],
    3: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'Oil ext.',  'CH&P', 
        'Utilities', 'HXN', 'Storage', 'Wastewater treatment'],
    4: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.', 
        'Wastewater treatment', 'Oil ext.', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    5: ['Juicing', 'Oil production', 'Biodiesel production',
        'Utilities', 'HXN', 'CH&P', 'Wastewater treatment', 'Storage'],
    6: ['Juicing', 'Pretreatment', 'Oil production', 'Biodiesel production',
        'Utilities', 'HXN', 'CH&P', 'Wastewater treatment', 'Storage'],
    # 9: ['Feedstock handling', 'Juicing', 'TAG prod.', 'Oil ext.', 'AcTAG sep.', 
    #     'Biod. prod.', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    # 10: ['Feedstock handling', 'Juicing', 'Pretreatment', 'TAG prod.', 'AcTAG sep.', 
    #      'Wastewater treatment', 'Oil ext.', 'CH&P',  'Biod. prod.', 'Utilities', 'HXN', 'Storage'],
}

area_names[7] = area_names[5]
area_names[9] = area_names[8] = area_names[6]

price_distributions_2022 = None
price_distributions_2023 = None
def get_price_distributions_module(year):
    if year == 2022:
        global price_distributions_2022
        if price_distributions_2022 is None: 
            from .data import price_distributions_2022
        return price_distributions_2022
    elif year == 2023:
        global price_distributions_2023
        if price_distributions_2023 is None: 
            from .data import price_distributions_2023
        return price_distributions_2023
    else:
        raise ValueError(f'invalid year {year}; year must be either 2022 or 2023')

def rename_storage_units(units, storage_area):
    bst.rename_units([i for i in units if bst.is_storage_unit(i)], storage_area)

def linear_fit(x0, x1, y0, y1):
    m = (y1 - y0) / (x1 - x0)
    b = y0 - m * x0
    return m, b

def linear_val(x, mb):
    m, b = mb
    return x * m + b

def exponential_fit(x0, x1, y0, y1): # A x ^ n
    x0 = log(x0)
    x1 = log(x1)
    y0 = log(y0)
    y1 = log(y1)
    n, logA = linear_fit(x0, x1, y0, y1)
    return n, exp(logA)

def exponential_val(x, nA): # A x ^ n
    n, A = nA    
    return A * x ** n

def YRCP2023():
    Biorefinery.default_conversion_performance_distribution = 'longterm'
    Biorefinery.default_prices_correleted_to_crude_oil = True
    Biorefinery.default_oil_content_range = [1.8, 5.4]
    Biorefinery.default_income_tax_range = [21, 21 + 1e-6] # Work around for constant value
    Biorefinery.default_baseline_oil_content = 1.8
    Biorefinery.default_year = 2023
    Biorefinery.default_WWT = 'high-rate'
    Biorefinery.default_dry_biomass_yield_distribution = shape.Uniform(
        0.19 * Biorefinery.baseline_dry_biomass_yield, Biorefinery.baseline_dry_biomass_yield
    )
    Biorefinery._derivative_disabled = True

class Biorefinery:
    cache = {}
    baseline_dry_biomass_yield = 25.62 # dry MT / ha / y
    baseline_available_land = 1600000 * 0.3 / baseline_dry_biomass_yield # ha
    set_feedstock_line = set_line_composition_parameters
    set_composition_by_line = set_composition_by_line
    default_prices_correleted_to_crude_oil = False
    default_conversion_performance_distribution = 'longterm'
    default_year = 2022
    default_WWT = None
    default_oil_content_range = [5, 15]
    default_income_tax_range = [21, 28] # Davis et al. 2018; https://www.nrel.gov/docs/fy19osti/71949.pdf
    default_baseline_oil_content = 10
    baseline_moisture_content = 0.70
    baseline_feedstock_CF = GWP_characterization_factors['sugarcane'] # [kg CO2e / kg sugarcane] with transportation
    baseline_feedstock_price = 0.035 # [USD / kg sugarcane] with transportation
    baseline_transportation_cost = 0.11 * baseline_feedstock_price # Energy cane usage for cellulosic ethanol: estimation of feedstock costs; https://ageconsearch.umn.edu/record/46837/files/Energy%20Cane%20Feedstock%20Estimation_SAEA1.pdf
    baseline_transportation_CF = 0.029 * baseline_feedstock_CF # Estimated from GREET 2020
    default_dry_biomass_yield_distribution = shape.Uniform(
        baseline_dry_biomass_yield, baseline_dry_biomass_yield + 1e-6 # Work around for constant value
    )
    
    cases = {
        'constant biomass yield',
    } # Special cases for Monte Carlo
    
    def constant_biomass_yield(self):
        self.set_dry_biomass_yield.distribution = shape.Uniform(
            self.baseline_dry_biomass_yield, 
            self.baseline_dry_biomass_yield + 1e-6 # Work around for constant value
        )
    
    def update_feedstock(self):
        feedstock = self.feedstock
        feedstock.price = self.feedstock_price
        feedstock.set_CF(GWP, self.feedstock_CF)
        self.sys.rescale(feedstock, self.feedstock_flow / feedstock.F_mass)
    
    @property
    def feedstock_price(self): # [USD / kg] with transportation
        return (
            (self.baseline_feedstock_price 
             - self.baseline_transportation_cost) * self.baseline_dry_biomass_yield / self.dry_biomass_yield
            + self.transportation_cost
        )
    
    @property
    def feedstock_CF(self): # [kg CO2e / kg] with transportation
        return (
            (self.baseline_feedstock_CF
             - self.baseline_transportation_CF) * self.baseline_dry_biomass_yield / self.dry_biomass_yield
            + self.transportation_CF
        )
    
    @property
    def transportation_cost(self): # [USD / kg] 
        return (
            self.baseline_transportation_cost
            * sqrt(self.available_land / self.baseline_available_land)
        )
    
    @property
    def transportation_CF(self): # [kg CO2e / kg] 
        return (
            self.baseline_transportation_CF
            * sqrt(self.available_land / self.baseline_available_land)
        )
    
    @property
    def chemicals(self):
        try:
            chemicals = self._chemicals
        except AttributeError:
            self._chemicals = chemicals = create_chemicals()
        return chemicals
    
    _derivative_disabled = False
    @classmethod
    def disable_derivative(cls, disable=True):
        cls._derivative_disabled = disable
    @classmethod
    def enable_derivative(cls, enable=True):
        cls._derivative_disabled = not enable
    
    @property
    def feedstock_flow(self): # kg / hr
        annual_crushing_capacity = self.annual_crushing_capacity
        feedstock = self.feedstock
        moisture_content = feedstock.get_mass_fraction('Water')
        dry_content = 1 - moisture_content
        return 1000 * annual_crushing_capacity / self.tea.operating_hours / dry_content
    
    @property
    def annual_crushing_capacity(self): # dry MT / y
        return self.available_land * self.dry_biomass_yield
    
    @property
    def baseline_annual_crushing_capacity(self): # dry MT / y
        return self.baseline_available_land * self.baseline_dry_biomass_yield
    
    @property
    def wet_biomass_yield(self):
        dry_biomass_yield = self.dry_biomass_yield
        moisture_content = self.feedstock.get_mass_fraction('Water')
        dry_fraction = 1 - moisture_content
        return dry_biomass_yield / dry_fraction
        
    def use_maximum_theoretical_microbial_oil_performance(self):
        max_lipid_yield_glucose = perf.max_lipid_yield_glucose
        max_lipid_yield_xylose = perf.max_lipid_yield_xylose
        self.set_glucose_to_microbial_oil_yield.setter(max_lipid_yield_glucose * 100)
        self.set_xylose_to_microbial_oil_yield.setter(max_lipid_yield_xylose * 100)
        self.set_fermentation_microbial_oil_productivity.setter(1.902) # Similar to sugarcane ethanol
        self.set_fermentation_microbial_oil_titer.setter(137) # Similar to sugarcane ethanol
    
    def use_upper_microbial_oil_performance(self):
        for i in (self.set_glucose_to_microbial_oil_yield,
                  self.set_xylose_to_microbial_oil_yield,
                  self.set_fermentation_microbial_oil_productivity,
                  self.set_fermentation_microbial_oil_titer):
            i.setter(i.bounds[1])
    
    def __new__(cls,
            name, chemicals=None,
            avoid_natural_gas=True, 
            conversion_performance_distribution=None,
            year=None, cache=cache, 
            prices_correleted_to_crude_oil=None,
            WWT_kwargs=None,
            oil_content_range=None,
            remove_biodiesel_production=None,
            update_feedstock_price=None,
            simulate=True,
        ):
        if update_feedstock_price is None: update_feedstock_price = True
        if year is None: year = cls.default_year
        if conversion_performance_distribution is None: 
            conversion_performance_distribution = cls.default_conversion_performance_distribution
        else:
            conversion_performance_distribution = conversion_performance_distribution.replace(' ', '').replace('-', '').lower()
        if prices_correleted_to_crude_oil is None:
            prices_correleted_to_crude_oil = cls.default_prices_correleted_to_crude_oil
        number, agile, feedstock_line, case = configuration = parse_configuration(name)
        if WWT_kwargs is None:
            if cls.default_WWT is None:
                WWT_key = None
            else:
                WWT_kwargs = dict(kind=cls.default_WWT)
                WWT_key = cls.default_WWT
        else:
            WWT_key = WWT_kwargs.get('kind', None)
        key = (number, agile, feedstock_line,
               WWT_key, conversion_performance_distribution,
               year, case)
        if cache and key in cache: 
            return cache[key]
        else:
            self = super().__new__(cls)
        self.configuration = configuration
        flowsheet_name = format_configuration(configuration, latex=False)
        flowsheet = bst.Flowsheet(flowsheet_name)
        main_flowsheet.set_flowsheet(flowsheet)
        if chemicals: self._chemicals = chemicals
        else: chemicals = self.chemicals
        bst.settings.set_thermo(chemicals)
        load_process_settings()
        s = flowsheet.stream
        u = flowsheet.unit
        get_stream = lambda ID: s[ID] if ID in flowsheet.stream else MockStream(ID) # Retrieve mock-stream in case stream does not exist in configuration
        
        ## System
        if number in system_factories:
            self.cane_sys = cane_sys = system_factories[number](operating_hours=24 * 200, WWT_kwargs=WWT_kwargs,
                                                                **system_factory_options.get(number, {}), autorename=True)
            if remove_biodiesel_production:
                unit = flowsheet(bst.BlendingTankWithSkimming)
                microbial_oil = unit.ins[0]
                microbial_oil.ID = 'microbial_oil'
                microbial_oil.price = 0.661
                remove_units = {unit, *unit.get_downstream_units()}
                self.cane_sys = cane_sys = bst.System.from_units(
                    'cane_sys',
                    units=set(cane_sys.units).difference(remove_units),
                )
                s.catalyst.empty()
                s.HCl.empty()
                s.methanol.empty()
        else:
            raise NotImplementedError(number)
        cane_sys.set_tolerance(
            rmol=1e-5, mol=1e-2, subsystems=True,
        )
        
        if number == 1:
            cane_sys.prioritize_unit(u.T508)
        elif number == 2:
            cane_sys.prioritize_unit(u.T808)
        
        if number < 0:
            s.sugarcane.register_alias('feedstock')
        else:
            s.oilcane.register_alias('feedstock')
        feedstock = s.feedstock
        
        chemicals.define_group(
            name='Fiber',
            IDs=['Cellulose', 'Hemicellulose', 'Lignin'],
            composition=[0.4704 , 0.2775, 0.2520],
            wt=True, # Composition is given as weight
        )
        chemicals.define_group(
            name='Sugar',
            IDs=['Sucrose', 'Glucose'],
            composition=[13.7, 1.21],
            wt=True,
        )
        
        BT = flowsheet(bst.BoilerTurbogenerator)
        if number not in cellulosic_configurations: BT.boiler_efficiency = 0.89
        
        ## Unit groups
        areas = area_names[number]
        rename_storage_units(cane_sys.units, (areas.index('Storage') + 1) * 100)
        unit_groups = UnitGroup.group_by_area(cane_sys.units)
        for i, j in zip(unit_groups, areas): i.name = j
        for i in unit_groups: i.autofill_metrics(shorthand=True)
        
        for HXN_group in unit_groups:
            if HXN_group.name == 'HXN':
                HXN_group.filter_savings = False # Allow negative values in heat utilities
                HXN = HXN_group.units[0]
                try:
                    assert isinstance(HXN, bst.HeatExchangerNetwork)
                except: breakpoint()
        HXN.raise_energy_balance_error = False
        HXN.vle_quenched_streams = False
        
        ## Split bagasse to boiler meet energy demand
        if avoid_natural_gas:
            for splitter in flowsheet.unit:
                if getattr(splitter, 'isbagasse_splitter', False): break
            else:
                avoid_natural_gas = False
        if avoid_natural_gas:
            cane_sys.add_specification(self.update_feedstock)
            self.bagasse_splitter = splitter
            minimum_fraction_burned = 0
            maximum_fraction_burned = 0.7
            @cane_sys.add_bounded_numerical_specification(
                x0=minimum_fraction_burned, x1=maximum_fraction_burned, 
                xtol=1e-4, ytol=100, args=(splitter,)
            )
            def adjust_bagasse_to_boiler(fraction_burned, splitter):
                # Returns energy consumption at given fraction processed (not sent to boiler).
                splitter.split[:] = 1 - fraction_burned
                cane_sys.simulate()
                excess = BT._excess_electricity_without_natural_gas
                if fraction_burned == minimum_fraction_burned and excess > 0:
                    splitter.neglect_natural_gas_streams = False # No need to neglect
                    return 0 # No need to burn bagasse
                elif fraction_burned == maximum_fraction_burned and excess < 0: 
                    splitter.neglect_natural_gas_streams = False # Cannot be neglected
                    return 0 # Cannot satisfy energy demand even at 30% sent to boiler (or minimum fraction processed)
                else:
                    splitter.neglect_natural_gas_streams = True
                    return excess
            @cane_sys.add_specification(args=(splitter,))
            def assume_negligible_natural_gas_streams(splitter):
                if splitter.neglect_natural_gas_streams:
                    for i in natural_gas_streams: i.empty()
        else:
            cane_sys.add_specification(self.update_feedstock, simulate=True)
        
        if abs(number) in cellulosic_configurations:
            prs = flowsheet(bst.PretreatmentReactorSystem)
            saccharification = flowsheet(bst.Saccharification)
            seed_train = flowsheet(bst.SeedTrain)
            try:
                fermentor = flowsheet(bst.CoFermentation)
            except:
                fermentor = flowsheet(bst.AeratedBioreactor)
            self.pretreatment_rxnsys = tmo.ReactionSystem(
                prs.reactions, saccharification.saccharification
            )
            self.fermentation_rxnsys = tmo.ReactionSystem(
                seed_train.reactions, fermentor.cofermentation
            )
            self.cellulosic_rxnsys = tmo.ReactionSystem(
                prs.reactions, saccharification.saccharification,
                seed_train.reactions, fermentor.cofermentation
            )
            saccharification.saccharification.X[0] = 0.0 # Baseline
            prs.reactions.X[10] = 0.0 # baseline
        else:
            try:
                fermentor = flowsheet(bst.Fermentation)
            except:
                fermentor = flowsheet(bst.AeratedBioreactor)
        self.fermentor = fermentor
        
        def set_glucose_yield(glucose_yield):
            if number in cellulosic_configurations:
                glucose_yield *= 0.01
                X1 = prs.reactions.X[0]
                X1_side = prs.reactions.X[1:3].sum()
                X2_side = saccharification.saccharification.X[:2].sum()
                saccharification.saccharification.X[2] = X2 = (glucose_yield - X1) / (1 - X1_side)
                X_excess = (X2_side + X2) - 1
                if X_excess > 0: breakpoint()
                
        def set_xylose_yield(xylose_yield):
            if number in cellulosic_configurations:
                xylose_yield *= 0.01
                X1_side = prs.reactions.X[9:11].sum()
                prs.reactions.X[8] = X1 = xylose_yield
                X_excess = (X1_side + X1) - 1
                if X_excess > 0.: breakpoint()
        
        if agile: 
            self.oilsorghum = oilsorghum = bst.Stream(
                ID='oilsorghum', phase='l', T=298.15, P=101325, 
                Water=2.333e+05, Glucose=3703, Sucrose=4.196e+04, Ash=2000, 
                Cellulose=2.227e+04, Hemicellulose=1.314e+04, Lignin=1.193e+04, 
                Solids=5000, units='kg/hr'
            )
            
            sys = bst.AgileSystem()
            @sys.operation_parameter(mode_dependent=True)
            def set_oil_content(oil_content, mode):
                F_mass = feedstock.F_mass
                feedstock.copy_flow(mode.feedstock)
                feedstock.F_mass = F_mass
                if number > 0: 
                    set_oil_fraction(oil_content, feedstock,
                                     FFA_fraction=mode.FFA_content,
                                     z_mass_carbs_baseline=mode.z_mass_carbs_baseline,
                                     PL_fraction=mode.PL_content)
            sys.operation_parameter(set_glucose_yield)
            sys.operation_parameter(set_xylose_yield)
            
            self.cane_mode = cane_mode = sys.operation_mode(cane_sys,
                operating_hours=180*24, oil_content=0.10, feedstock=feedstock.copy(),
                z_mass_carbs_baseline=0.1491, glucose_yield=85, xylose_yield=65, 
                FFA_content=0.10, PL_content=0.10
            )
            self.sorghum_mode = sorghum_mode = sys.operation_mode(cane_sys, 
                operating_hours=60*24, oil_content=0.07, glucose_yield=79, xylose_yield=86,
                feedstock=oilsorghum,
                z_mass_carbs_baseline=0.1371, FFA_content=0.10, PL_content=0.10,
            )
            tea = create_tea(sys)
            tea.operating_days = 260 
            tea.IRR = 0.10
        else:
            sys = cane_sys
            tea = create_tea(sys)
            tea.operating_days = 200
            tea.IRR = 0.10
        
        ## Specifications for analysis
        self.composition_specification = composition_specification = CaneCompositionSpecification(feedstock)
        if number < 0:
            oil_extraction_specification = MockExtractionSpecification()
        else:
            crushing_mill = flowsheet(bst.CrushingMill) # Separates bagasse
            _, screen = flowsheet(bst.VibratingScreen)
            juice = screen.outs[0]
            microbial_oil_centrifuge = None
            if number in cellulosic_configurations:
                pressure_filter = flowsheet(bst.PressureFilter) # Separates lignin
            else:
                pressure_filter  = None
            if 4 < number < 10:
                cellmass_centrifuge = flowsheet(bst.SolidLiquidsSplitCentrifuge) # Separates cell mass, oil, and water
                if isinstance(cellmass_centrifuge, list):
                    cellmass_centrifuge = sorted([i for i in cellmass_centrifuge if type(i) is bst.SolidsCentrifuge], key=lambda x: x.ID)[0]
                cellmass_centrifuge.strict_moisture_content = False
            else:
                cellmass_centrifuge = None
            if number in screwpress_microbial_oil_recovery:
                screw_press = flowsheet(bst.ScrewPress) # Separates oil from cell mass
            else:
                screw_press = None
            if number == 8:
                for microbial_oil_centrifuge in flowsheet.unit:
                    if microbial_oil_centrifuge.line == 'Microbial oil centrifuge': break
            else:
                microbial_oil_centrifuge = None
            
            oil_extraction_specification = OilExtractionSpecification(
                sys, crushing_mill, pressure_filter, cellmass_centrifuge, 
                microbial_oil_centrifuge, juice, screw_press, fermentor,
            )
        
        ## Account for cellulosic vs advanced RINs
        if number in cellulosic_ethanol_configurations:
            # Note that GREET cellulosic ethanol from corn stover results in a 
            # GWP of 0.41 kg CO2e / L-ethanol. So the cellulosic ethanol from
            # bagasse (~0.34, 0.32 for configurations S2 and O2) can certainly
            # apply as an EPA cellulosic ethanol pathway.
            RIN_splitter = bst.Splitter('RIN_splitter',
                ins=s.ethanol,
                outs=['advanced_ethanol', 'cellulosic_ethanol'],
                split=0.5
            )
            @RIN_splitter.add_specification(run=True)
            def adjust_ethanol_split():
                # outs: stream sequence
                # [0] Advanced biofuel ethanol
                # [1] Cellulosic biofuel ethanol
                glucose_yield = fermentor.cofermentation[0].product_yield('Ethanol', 'wt')
                xylose_yield = fermentor.cofermentation[1].product_yield('Ethanol', 'wt')
                juice_sugar = s.juice.imass['Glucose', 'Sucrose'].sum() * glucose_yield
                hydrolysate_sugar = (
                    s.slurry.imass['Glucose'] * glucose_yield
                    + s.slurry.imass['Xylose', 'Arabinose'].sum() * xylose_yield
                )
                RIN_splitter.split[:] = juice_sugar / (juice_sugar + hydrolysate_sugar)
            
            cane_sys.update_configuration([*cane_sys.units, RIN_splitter])
            assert RIN_splitter in cane_sys.units
        elif number in conventional_ethanol_configurations:
            s.ethanol.ID = 'advanced_ethanol'
            s.advanced_ethanol.register_alias('ethanol')
        elif number in conventional_ethanol_configurations:
            s.ethanol.ID = 'advanced_ethanol'
            s.advanced_ethanol.register_alias('ethanol')
        if number in biodiesel_configurations:
            if (number not in cellulosic_ethanol_configurations
                and number in cellulosic_configurations):
                # Note that GREET cellulosic ethanol from corn stover results in a 
                # GWP of 0.41 kg CO2e / L-ethanol. So the cellulosic ethanol from
                # bagasse (~0.34, 0.32 for configurations S2 and O2) can certainly
                # apply as an EPA cellulosic ethanol pathway.
                RIN_splitter = bst.Splitter('RIN_splitter',
                    ins=s.biodiesel,
                    outs=['biomass_based_diesel', 'cellulosic_based_diesel'],
                    split=0.5
                )
                @RIN_splitter.add_specification(run=True)
                def adjust_biodiesel_split():
                    # outs: stream sequence
                    # [0] Biomass based diesel
                    # [1] Cellulosic based diesel
                    glucose_yield = fermentor.cofermentation[0].product_yield('TriOlein', 'wt')
                    xylose_yield = fermentor.cofermentation[1].product_yield('TriOlein', 'wt')
                    juice_sugar = s.juice.imass['Glucose', 'Sucrose'].sum() * glucose_yield
                    hydrolysate_sugar = (
                        s.slurry.imass['Glucose'] * glucose_yield
                        + s.slurry.imass['Xylose', 'Arabinose'].sum() * xylose_yield
                    )
                    RIN_splitter.split[:] = juice_sugar / (juice_sugar + hydrolysate_sugar)
                
                cane_sys.update_configuration([*cane_sys.units, RIN_splitter])
                assert RIN_splitter in cane_sys.units
            else:
                # A biodiesel stream should already exist
                s.biodiesel.register_alias('biomass_based_diesel')
        
        ## Additional modifications for speed
        for i in sys.units:
            if isinstance(i, bst.MultiEffectEvaporator): i.flash = False
        HXN.acceptable_energy_balance_error = 0.02
        HXN.replace_unit_heat_utilities = False
        HXN.force_ideal_thermo = True
        HXN.cache_network = True
        HXN.avoid_recycle = True
        try: AD = flowsheet(bst.AerobicDigestion)
        except: pass
        else:
            WWTsys = cane_sys.find_system(AD)
            WWTsys.set_tolerance(mol=10, method='fixed-point')
        
        ## LCA
        # Set non-negligible characterization factors
        stream_names = ('FGD_lime', 'cellulase', 'urea', 'caustic', 'cellulosic_ethanol',
                        'catalyst', 'methanol', 'HCl', 'NaOH', 'crude_glycerol', 'pure_glycerine',
                        'denaturant', 'ethanol', 'advanced_ethanol', 'cellulosic_ethanol',
                        'biomass_based_diesel', 'cellulosic_based_diesel', 'natural_gas',
                        'biodiesel', 'lime', 'H3PO4')
        (FGD_lime, cellulase, urea, caustic, cellulosic_ethanol,
         catalyst, methanol, HCl, NaOH, crude_glycerol, pure_glycerine,
         denaturant, ethanol, advanced_ethanol, cellulosic_ethanol,
         biomass_based_diesel, cellulosic_based_diesel, natural_gas, 
         biodiesel, lime, H3PO4) = [get_stream(i) for i in stream_names]
        
        def get_renamed_stream(search, ID):
            if search in flowsheet.stream:
                stream = s[search]  
                stream.ID = ID
            else:
                stream = MockStream(ID)
            return stream
        
        NaOCl = get_renamed_stream('naocl_R602', 'NaOCl')
        citric_acid = get_renamed_stream('citric_R602', 'citric_acid')
        bisulfite = get_renamed_stream('bisulfite_R602', 'bisulfite')
        polishing_filter_air = get_renamed_stream('air_R603', 'polishing_filter_air')
        polishing_filter_vent = get_renamed_stream('vent_R603', 'polishing_filter_vent')
        set_GWPCF(NaOCl, 'NaOCl')
        set_GWPCF(citric_acid, 'citric acid')
        set_GWPCF(bisulfite, 'bisulfite')
        set_GWPCF(feedstock, 'sugarcane')
        set_GWPCF(H3PO4, 'H3PO4')
        set_GWPCF(lime, 'lime', dilution=0.046) # Diluted with water
        set_GWPCF(denaturant, 'gasoline')
        set_GWPCF(FGD_lime, 'lime', dilution=0.451)
        set_GWPCF(cellulase, 'cellulase', dilution=0.05) 
        set_GWPCF(urea, 'urea')
        set_GWPCF(caustic, 'NaOH', 0.5)
        set_GWPCF(catalyst, 'methanol catalyst mixture')
        set_GWPCF(methanol, 'methanol')
        set_GWPCF(HCl, 'HCl')
        set_GWPCF(NaOH, 'NaOH')
        set_GWPCF(pure_glycerine, 'pure-glycerol')
        set_GWPCF(crude_glycerol, 'crude-glycerol', dilution=0.80)
        set_GWPCF(biodiesel, 'biodiesel displacement')
        bst.PowerUtility.set_CF(GWP, GWP_characterization_factors['Electricity'])
        self.natural_gas_streams = natural_gas_streams = [natural_gas]
        for stream in natural_gas_streams: set_GWPCF(stream, 'CH4')
        
        ## Add BioSTEAM objects to module for easy access
        self.price_distribution_module = dist = get_price_distributions_module(year)
        
        ## Model
        model = bst.Model(sys, exception_hook='raise', retry_evaluation=False)
        parameter = model.parameter
        metric = model.metric
        
        # DO NOT DELETE:
        # natural_gas.phase = 'g'
        # natural_gas.set_property('T', 60, 'degF')
        # natural_gas.set_property('P', 14.73, 'psi')
        # original_value = natural_gas.imol['CH4']
        # natural_gas.imass['CH4'] = 1 
        # V_ng = natural_gas.get_total_flow('m3/hr')
        # natural_gas.imol['CH4'] = original_value
        V_ng = 1.473318463076884 # Natural gas volume at 60 F and 14.73 psi [m3 / kg]
        
        def potential_shortterm_gain_in_performance(lb, ub, *args, **kwargs):
            distribution = shape.Uniform(lb, lb + 0.15 * (ub - lb))
            return parameter(*args, distribution=distribution, bounds=(lb, ub), baseline=lb, **kwargs)
        
        def potential_longterm_gain_in_performance(lb, ub, *args, **kwargs):
            distribution = shape.Uniform(lb, ub)
            return parameter(*args, distribution=distribution, bounds=(lb, ub), baseline=lb, **kwargs)
        
        def uniform(lb, ub, *args, **kwargs):
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), **kwargs)
        
        def default(baseline, *args, **kwargs):
            lb = 0.75*baseline
            ub = 1.25*baseline
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub),
                             baseline=baseline, **kwargs)
        
        def default_gwp(baseline, *args, **kwargs):
            lb = 0.90*baseline
            ub = 1.10*baseline
            return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), 
                             baseline=baseline, **kwargs)
        
        def triangular(lb, mid, ub, *args, **kwargs):
            return parameter(*args, distribution=shape.Triangle(lb, mid, ub), bounds=(lb, ub), **kwargs)
        
        if oil_content_range is None:
            oil_content_range = self.default_oil_content_range
        
        if conversion_performance_distribution == 'shortterm':
            performance = potential_shortterm_gain_in_performance
        elif conversion_performance_distribution == 'longterm':
            performance = potential_longterm_gain_in_performance
        else:
            raise ValueError("`conversion_performance_distribution` must be either 'longterm' or 'shortterm")
            
        @performance(60, 90, units='%', kind='coupled')
        def set_juicing_oil_recovery(juicing_oil_recovery):
            oil_extraction_specification.load_juicing_oil_recovery(juicing_oil_recovery / 100.)
        
        @performance(70 if number in screwpress_microbial_oil_recovery else 50,
                     90, units='%', kind='coupled')
        def set_microbial_oil_recovery(microbial_oil_recovery):
            oil_extraction_specification.load_microbial_oil_recovery(microbial_oil_recovery / 100.)
        
        @performance(70.0, 90, units='%', kind='coupled')
        def set_bagasse_oil_recovery(bagasse_oil_recovery):
            oil_extraction_specification.load_bagasse_oil_recovery(bagasse_oil_recovery / 100.)
    
        # Baseline from Huang's 2016 paper, but distribution more in line with Florida sugarcane harvesting (3-5 months)
        @uniform(4 * 30, 6 * 30, units='day/y', baseline=180, kind='coupled')
        def set_cane_operating_days(cane_operating_days):
            if agile:
                cane_mode.operating_hours = cane_operating_days * 24
            else:
                tea.operating_days = cane_operating_days
    
        # From Ed Cahoon and Huang 2017
        @uniform(30, 60, units='day/y', baseline=45, kind='coupled')
        def set_sorghum_operating_days(sorghum_operating_days):
            if agile: sorghum_mode.operating_hours = sorghum_operating_days * 24
        
        @default(self.baseline_available_land, element='feedstock', units='ha', kind='coupled')
        def set_available_land(available_land):
            self.available_land = available_land
            
        @parameter(
            element='feedstock', 
            baseline=self.baseline_dry_biomass_yield,
            distribution=self.default_dry_biomass_yield_distribution,
            units='DMT/ha/y', 
            kind='coupled',
        )
        def set_dry_biomass_yield(dry_biomass_yield):
            if number < 0:
                self.dry_biomass_yield = self.baseline_dry_biomass_yield 
            else:
                self.dry_biomass_yield = dry_biomass_yield
    
        @parameter(distribution=dist.copd, element='crude oil', 
                   baseline=dist.mcop, units='USD/L')
        def set_crude_oil_price(price):
            self.crude_oil_price = price
        
        @default(self.baseline_feedstock_price, element='feedstock', units='USD/kg')
        def set_baseline_feedstock_price(price):
            self.baseline_feedstock_price = price
        
        if prices_correleted_to_crude_oil:
            predict = lambda name, scalar: dist.models[name].predict(np.array([[scalar]]))[0]
                                                                     
            @parameter(distribution=dist.residual_distributions['Cellulosic ethanol'],
                       element='Cellulosic ethanol', baseline=0., units='USD/L')
            def set_cellulosic_ethanol_price(price): 
                cellulosic_ethanol.price = predict('Cellulosic ethanol', self.crude_oil_price + price) * ethanol_L_per_kg
                
            @parameter(distribution=dist.residual_distributions['Advanced ethanol'],
                       element='Advanced ethanol', baseline=0., units='USD/L')
            def set_advanced_ethanol_price(price): 
                advanced_ethanol.price =  predict('Advanced ethanol', self.crude_oil_price + price) * ethanol_L_per_kg
                
            # USDA ERS historical price data
            @parameter(distribution=dist.residual_distributions['Biomass based diesel'], 
                       element='Biomass based diesel', units='USD/L', baseline=0.)
            def set_biomass_based_diesel_price(price):
                biomass_based_diesel.price =  predict('Biomass based diesel', self.crude_oil_price + price) * biodiesel_L_per_kg
        
            @parameter(distribution=dist.residual_distributions['Cellulosic based diesel'],
                       element='Cellulosic based diesel', units='USD/L', baseline=0.)
            def set_cellulosic_based_diesel_price(price):
                cellulosic_based_diesel.price =  predict('Cellulosic based diesel', self.crude_oil_price + price) * biodiesel_L_per_kg
        
            # https://www.eia.gov/energyexplained/natural-gas/prices.php
            @parameter(distribution=dist.residual_distributions['Natural gas'],
                       element='Natural gas', units='USD/m3', baseline=0.)
            def set_natural_gas_price(price): 
                BT.natural_gas_price =  predict('Natural gas', self.crude_oil_price + price) * V_ng
        
            @parameter(distribution=dist.residual_distributions['Electricity'],
                       element='Electricity', units='USD/kWh', baseline=0.)
            def set_electricity_price(price): 
                bst.PowerUtility.price = predict('Electricity', self.crude_oil_price + price)
                
        else:
            # USDA ERS historical price data with EPA RIN prices
            @parameter(distribution=dist.cepd, element='Cellulosic ethanol', 
                       baseline=dist.mcep, units='USD/L')
            def set_cellulosic_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
                cellulosic_ethanol.price = price * ethanol_L_per_kg
                
            @parameter(distribution=dist.aepd, element='Advanced ethanol', 
                       baseline=dist.maep, units='USD/L')
            def set_advanced_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
                advanced_ethanol.price = price * ethanol_L_per_kg
                
            # USDA ERS historical price data
            @parameter(distribution=dist.bpd, element='Biomass based diesel', units='USD/L', baseline=dist.mbp)
            def set_biomass_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
                biomass_based_diesel.price = price * biodiesel_L_per_kg
        
            @parameter(distribution=dist.cbpd, element='Cellulosic based diesel', units='USD/L', baseline=dist.mcbp)
            def set_cellulosic_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
                cellulosic_based_diesel.price = price * biodiesel_L_per_kg
        
            # https://www.eia.gov/energyexplained/natural-gas/prices.php
            @parameter(distribution=dist.natural_gas_price_distribution, element='Natural gas', units='USD/m3',
                       baseline=4.73 * 35.3146667/1e3)
            def set_natural_gas_price(price): 
                BT.natural_gas_price = price * V_ng
        
            @parameter(distribution=dist.electricity_price_distribution, units='USD/kWh',
                       element='electricity', baseline=dist.mean_electricity_price)
            def set_electricity_price(price): 
                bst.PowerUtility.price = price
        
        # 10% is suggested for waste reducing, but 15% is suggested for investment
        @uniform(10., 15., units='%', baseline=10, element='')
        def set_IRR(IRR):
            tea.IRR = IRR / 100.
        
        @uniform(0.10, 0.22, units='USD/kg', element=crude_glycerol.ID)
        def set_crude_glycerol_price(price):
            crude_glycerol.price = price
    
        @default(0.65, units='USD/kg', element=pure_glycerine.ID)
        def set_pure_glycerol_price(price):
            pure_glycerine.price = price
        
        @default(72, units='h', element='Saccharification')
        def set_saccharification_reaction_time(reaction_time):
            if abs(number) in cellulosic_configurations: saccharification.tau = reaction_time
        
        @default(0.212, units='USD/kg', element=cellulase.ID)
        def set_cellulase_price(price):
            cellulase.price = price
    
        if abs(number) in cellulosic_configurations:
            cellulase_mixer, = [i for i in flowsheet.unit if hasattr(i, 'enzyme_loading')]
       
        @default(0.02, units='wt. % cellulose', element='cellulase', kind='coupled')
        def set_cellulase_loading(cellulase_loading):
            if abs(number) in cellulosic_configurations: cellulase_mixer.enzyme_loading = cellulase_loading
        
        @default(PRS_cost_item.cost / 1e6, units='million USD', element='Pretreatment reactor system')
        def set_reactor_base_cost(base_cost):
            PRS_cost_item.cost = base_cost * 1e6
        
        @performance(85, 97.5, units='%', element='Pretreatment and saccharification', kind='coupled')
        def set_cane_glucose_yield(cane_glucose_yield):
            if agile:
                cane_mode.glucose_yield = cane_glucose_yield
            elif abs(number) in cellulosic_configurations:
                set_glucose_yield(cane_glucose_yield)
        
        @uniform(79, 97.5, units='%', element='Pretreatment and saccharification',
                 baseline=79, kind='coupled')
        def set_sorghum_glucose_yield(sorghum_glucose_yield):
            if not agile: return
            sorghum_mode.glucose_yield = sorghum_glucose_yield
            
        @performance(65, 97.5, units='%', element='Pretreatment and saccharification', kind='coupled')
        def set_cane_xylose_yield(cane_xylose_yield):
            if agile:
                cane_mode.xylose_yield = cane_xylose_yield
            elif abs(number) in cellulosic_configurations:
                set_xylose_yield(cane_xylose_yield)
        
        @performance(86, 97.5, units='%', element='Pretreatment and saccharification', kind='coupled')
        def set_sorghum_xylose_yield(sorghum_xylose_yield):
            if not agile: return
            sorghum_mode.xylose_yield = sorghum_xylose_yield
        
        @performance(90, 95, units='%', element='Cofermenation', kind='coupled')
        def set_glucose_to_ethanol_yield(glucose_to_ethanol_yield):
            if number in cellulosic_ethanol_configurations:
                glucose_to_ethanol_yield *= 0.01
                # fermentor.cofermentation[2].X = 0.004 # Baseline
                # fermentor.cofermentation[3].X = 0.006 # Baseline
                # fermentor.loss[0].X = 0.03 # Baseline
                split = u.S401.split.mean()
                X1 = split * seed_train.reactions.X[0]
                X2 = split * seed_train.reactions.X[2]
                X3 = (glucose_to_ethanol_yield - X1) / (1 - X1 - X2)
                X_excess = X3 * 1.0526 - 1
                if X_excess > 0.: breakpoint()
                fermentor.cofermentation.X[0] = X3
                fermentor.cofermentation.X[2] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass
        
        @performance(50, 95, units='%', element='Cofermenation', kind='coupled')
        def set_xylose_to_ethanol_yield(xylose_to_ethanol_yield):
            if number in cellulosic_ethanol_configurations:
                xylose_to_ethanol_yield *= 0.01
                split = u.S401.split.mean()
                X1 = split * seed_train.reactions.X[1]
                X2 = split * seed_train.reactions.X[3]
                X3 = (xylose_to_ethanol_yield - X1) / (1 - X1 - X2)
                X_excess = X3 * 1.0526 - 1
                if X_excess > 0.: breakpoint()
                fermentor.cofermentation.X[1] = X3
                fermentor.cofermentation.X[3] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass
    
        @performance(68.5, 137, units='g/L', element='Cofermentation', kind='coupled')
        def set_cofermentation_ethanol_titer(ethanol_titer):
            if number in cellulosic_ethanol_configurations: fermentor.titer = ethanol_titer
    
        @performance(0.951, 1.902, units='g/L/h', element='Cofermentation')
        def set_cofermentation_ethanol_productivity(ethanol_productivity):
            if number in cellulosic_ethanol_configurations: fermentor.productivity = ethanol_productivity
    
        fed_batch = system_factory_options.get(number, {}).get('fed_batch')
        if fed_batch:
            lipid_yield = perf.fed_batch_lipid_yield_mean
            titer = perf.fed_batch_titer_mean
            productivity = perf.fed_batch_productivity_mean
        else:
            lipid_yield = perf.batch_lipid_yield_mean
            titer = perf.batch_titer_mean
            productivity = perf.batch_productivity_mean
        
        hydrolysate_lipid_yield = perf.hydrolysate_lipid_yield
        hydrolysate_titer = perf.hydrolysate_titer
        hydrolysate_productivity = perf.hydrolysate_productivity
            
        @performance(hydrolysate_lipid_yield * 100, lipid_yield * 100, units='%', element='Cofermenation', kind='coupled')
        def set_glucose_to_microbial_oil_yield(glucose_to_microbial_oil_yield):
            if number in cellulosic_oil_configurations:
                fermentation_reaction = fermentor.cofermentation[0]
                cell_growth_reaction = fermentor.cofermentation[2]
            elif number in oil_configurations:
                fermentation_reaction = fermentor.fermentation_reaction
                cell_growth_reaction = fermentor.cell_growth_reaction
            else:
                return
            glucose_to_microbial_oil_yield *= 0.01
            fermentation_reaction.product_yield(
                product='TAG', product_yield=glucose_to_microbial_oil_yield, basis='wt'
            )
            # Almost all the rest goes towards cell mass and CO2
            if number in cellulosic_oil_configurations:
                # Reaction in parallel
                cell_growth_reaction.X = 0.999 - fermentation_reaction.X
            else:
                # Reaction in series
                cell_growth_reaction.X = 0.999 
            
        @performance(hydrolysate_lipid_yield * 100, lipid_yield * 100, units='%', element='Cofermenation', kind='coupled')
        def set_xylose_to_microbial_oil_yield(xylose_to_microbial_oil_yield):
            if number in cellulosic_oil_configurations:
                xylose_to_microbial_oil_yield *= 0.01
                fermentation_reaction = fermentor.cofermentation[1]
                cell_growth_reaction = fermentor.cofermentation[3]
                fermentation_reaction.product_yield(
                    product='TAG', product_yield=xylose_to_microbial_oil_yield, basis='wt'
                )
                cell_growth_reaction.X = 0.999 - fermentation_reaction.X # Almost all the rest goes towards cell mass and CO2
    
        @performance(hydrolysate_titer, titer, units='g/L', element='Fermentation', kind='coupled')
        def set_fermentation_microbial_oil_titer(microbial_oil_titer):
            if number in oil_configurations: fermentor.titer = microbial_oil_titer
    
        @performance(hydrolysate_productivity, productivity, units='g/L/h', element='Fermentation', kind='coupled')
        def set_fermentation_microbial_oil_productivity(microbial_oil_productivity):
            if number in oil_configurations: fermentor.productivity = microbial_oil_productivity
    
        @default(10, element='oilcane', units='% oil', kind='coupled')
        def set_cane_PL_content(cane_PL_content):
            if agile: cane_mode.PL_content = cane_PL_content / 100.
            else: composition_specification.PL = cane_PL_content / 100.
        
        @default(10, element='oilsorghum', units='% oil', kind='coupled')
        def set_sorghum_PL_content(sorghum_PL_content):
            if agile: sorghum_mode.PL_content = sorghum_PL_content / 100.
        
        @default(10, element='oilcane', units='% oil', kind='coupled')
        def set_cane_FFA_content(cane_FFA_content):
            if agile: cane_mode.FFA_content = cane_FFA_content / 100.
            else: composition_specification.FFA = cane_FFA_content / 100.
        
        @default(10, element='oilsorghum', units='% oil', kind='coupled')
        def set_sorghum_FFA_content(sorghum_FFA_content):
            if agile: sorghum_mode.FFA_content = sorghum_FFA_content / 100. 
    
        @uniform(*oil_content_range, element='oilcane', units='dry wt. %', kind='coupled')
        def set_cane_oil_content(cane_oil_content):
            if number < 0: return
            if agile:
                cane_mode.oil_content = cane_oil_content / 100.
            else:
                composition_specification.load_oil_content(
                    cane_oil_content / 100.
                )
    
        @uniform(-3., 0., element='oilsorghum', units='dry wt. %', kind='coupled',
                 baseline=0.)
        def set_relative_sorghum_oil_content(relative_sorghum_oil_content):
            if number < 0: return
            if agile:
                sorghum_mode.oil_content = cane_mode.oil_content + relative_sorghum_oil_content / 100.
    
        @default(23, units='% oil', kind='coupled', name='TAG to FFA conversion')
        def set_TAG_to_FFA_conversion(TAG_to_FFA_conversion):
            if number in conventional_ethanol_biodiesel_configurations:
                u.R301.oil_reaction.X[0] = TAG_to_FFA_conversion / 100.
            elif number in cellulosic_ethanol_biodiesel_configurations:
                u.R401.oil_reaction.X[0] = TAG_to_FFA_conversion / 100.
        
        @default_gwp(feedstock.characterization_factors[GWP], name='GWP', 
                 element='sugarcane', units='kg*CO2e/kg')
        def set_feedstock_GWP(value):
            self.baseline_feedstock_CF = value
        
        @default_gwp(methanol.characterization_factors[GWP], name='GWP', 
                     element=methanol.ID, units='kg*CO2e/kg')
        def set_methanol_GWP(value):
            methanol.characterization_factors[GWP] = value
        
        # @default(crude_glycerol.characterization_factors[GWP], name='GWP', 
        #          element=crude_glycerol, units='kg*CO2e/kg')
        # def set_crude_glycerol_GWP(value):
        #     crude_glycerol.characterization_factors[GWP] = value
        
        @default_gwp(pure_glycerine.characterization_factors[GWP], name='GWP', 
                     element=pure_glycerine.ID, units='kg*CO2e/kg')
        def set_pure_glycerine_GWP(value):
            pure_glycerine.characterization_factors[GWP] = value
        
        @default_gwp(cellulase.characterization_factors[GWP], name='GWP', 
                     element=cellulase.ID, units='kg*CO2e/kg')
        def set_cellulase_GWP(value):
            cellulase.characterization_factors[GWP] = value * 0.02
        
        @default_gwp(natural_gas.characterization_factors[GWP], name='GWP', 
                     element=natural_gas.ID, units='kg*CO2e/kg')
        def set_natural_gas_GWP(value):
            for ng in natural_gas_streams:
                ng.characterization_factors[GWP] = value
        
        income_tax_lb, income_tax_ub = self.default_income_tax_range
        tea.income_tax = income_tax_lb / 100
        @uniform(income_tax_lb, income_tax_ub, baseline=income_tax_lb, name='Income tax', units='%')
        def set_income_tax(income_tax):
            tea.income_tax = income_tax * 0.01
                    
        if agile:
            feedstock_flow = lambda: sys.flow_rates[feedstock] / kg_per_MT # MT / y
            biodiesel_flow = lambda: (sys.flow_rates.get(cellulosic_based_diesel, 0.) + sys.flow_rates.get(biomass_based_diesel, 0.)) * biodiesel_L_per_kg # L / y
            ethanol_flow = lambda: (sys.flow_rates.get(cellulosic_ethanol, 0.) + sys.flow_rates.get(advanced_ethanol, 0.)) * ethanol_L_per_kg # L / y
            natural_gas_flow = lambda: sum([sys.flow_rates[i] for i in natural_gas_streams]) * V_ng # m3 / y
            crude_glycerol_flow = lambda: sys.flow_rates.get(crude_glycerol, 0.) # kg / y
            
            @sys.operation_metric(annualize=True)
            def direct_nonbiogenic_emissions(mode):
                return sum([i.F_mol for i in natural_gas_streams]) * chemicals.CO2.MW
            
        else:
            feedstock_flow = lambda: sys.operating_hours * feedstock.F_mass / kg_per_MT # MT / y
            biodiesel_flow = lambda: sys.operating_hours * biodiesel.F_mass * biodiesel_L_per_kg # L / y
            ethanol_flow = lambda: sys.operating_hours * ethanol.F_mass * ethanol_L_per_kg # L / y
            crude_glycerol_flow = lambda: sys.operating_hours * crude_glycerol.F_mass # kg / y
            natural_gas_flow = lambda: sum([i.F_mass for i in natural_gas_streams]) * sys.operating_hours * V_ng # m3 / y
            direct_nonbiogenic_emissions = lambda: sum([i.F_mol for i in natural_gas_streams]) * chemicals.CO2.MW * sys.operating_hours
        electricity = lambda: sys.operating_hours * sys.power_utility.rate
        
        def get_GWP_mean_ethanol_price():
            if number in cellulosic_ethanol_configurations:
                f = float(RIN_splitter.split.mean())
                return f * dist.mean_advanced_ethanol_price + (1 - f) * dist.mean_cellulosic_ethanol_price
            else:
                return dist.mean_advanced_ethanol_price
        
        def get_GWP_mean_biodiesel_price():
            if number in cellulosic_configurations:
                f = float(RIN_splitter.split.mean())
                return f * dist.mean_biomass_based_diesel_price + (1 - f) * dist.mean_cellulosic_based_diesel_price
            else:
                return dist.mean_advanced_ethanol_price
        
        sys.define_process_impact(
            key=GWP,
            name='Direct non-biogenic emissions',
            basis='kg',
            inventory=direct_nonbiogenic_emissions,
            CF=1.,
        )
        
        self.flows = {
            'feedstock': feedstock_flow,
            'biodiesel': biodiesel_flow,
            'ethanol': ethanol_flow,
            'natural_gas': natural_gas_flow,
            'crude_glycerol': crude_glycerol_flow,
            'electricity': electricity
        }
        @metric(units='USD/MT')
        def MFPP():
            price = tea.solve_price(feedstock)
            return kg_per_MT * price
        
        @metric(units='USD/L')
        def MESP():
            ethanol_streams = [advanced_ethanol, cellulosic_ethanol]
            ethanol_streams = [i for i in ethanol_streams if not i.isempty()]
            if not ethanol_streams: return 0.
            price = tea.solve_price(ethanol_streams)
            return price / ethanol_L_per_kg
        
        @metric(units='USD/L')
        def MBSP():
            biodiesel_streams = [biomass_based_diesel, cellulosic_based_diesel]
            biodiesel_streams = [i for i in biodiesel_streams if not i.isempty()]
            if not biodiesel_streams: return 0.
            price = tea.solve_price(biodiesel_streams)
            return price / biodiesel_L_per_kg
        
        @metric(units='MT/y')
        def feedstock_consumption():
            return feedstock_flow()
        
        @metric(units='L/MT')
        def biodiesel_production():
            return biodiesel_flow() / feedstock_consumption.get()
        
        @metric(units='L/ha')
        def biodiesel_yield():
            return biodiesel_production.get() * self.wet_biomass_yield
        
        @metric(units='L/MT')
        def ethanol_production():
            return ethanol_flow() / feedstock_consumption.get()
        
        @metric(units='kWh/MT')
        def electricity_production():
            value = - electricity() / feedstock_consumption.get()
            if value < 0.: value = 0.
            return value
        
        @metric(units='GGE/MT')
        def net_energy_production():
            GGE_biodiesel_annual = biodiesel_production.get() / 0.9536 / L_per_gal
            GGE_ethanol_annual = ethanol_production.get() / 1.5 / L_per_gal
            GEE_electricity_production = max(-electricity() * 3600 / 114000, 0.) / feedstock_consumption.get()
            GEE_crude_glycerol = crude_glycerol_flow() * 0.1059 / feedstock_consumption.get()
            return (GGE_biodiesel_annual + GGE_ethanol_annual + GEE_electricity_production + GEE_crude_glycerol)
        
        @metric(units='m3/MT')
        def natural_gas_consumption():
            value = natural_gas_flow() / feedstock_consumption.get()
            return value
        
        @metric(units='10^6*USD')
        def TCI():
            return tea.TCI / 1e6 # 10^6*$
        
        @metric(units='%')
        def heat_exchanger_network_error():
            return HXN.energy_balance_percent_error if HXN else 0.
    
        @metric(name='GWP', element='Economic allocation', units='kg*CO2e / USD')
        def GWP_economic(): # Cradle to gate
            GWP_material = sys.get_total_feeds_impact(GWP) # kg CO2 eq. / y
            GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / y
            sales = (
                biodiesel_flow() * get_GWP_mean_biodiesel_price()
                + ethanol_flow() * get_GWP_mean_ethanol_price()
                + crude_glycerol_flow() * dist.mean_glycerol_price
                + max(-electricity(), 0) * dist.mean_electricity_price
            )
            return (GWP_material + GWP_emissions) / sales
    
        @metric(name='Ethanol GWP', element='Economic allocation', units='kg*CO2e / L')
        def GWP_ethanol(): # Cradle to gate
            return GWP_economic.get() * get_GWP_mean_ethanol_price()
        
        @metric(name='Biodiesel GWP', element='Economic allocation', units='kg*CO2e / L')
        def GWP_biodiesel(): # Cradle to gate
            if number > 0:
                return GWP_economic.get() * get_GWP_mean_biodiesel_price()
            else:
                return 0.
        
        @metric(name='Crude glycerol GWP', element='Economic allocation', units='kg*CO2e / kg')
        def GWP_crude_glycerol(): # Cradle to gate
            if number in biodiesel_configurations:
                return GWP_economic.get() * dist.mean_glycerol_price
            else:
                return 0.
        
        @metric(name='Electricity GWP', element='Economic allocation', units='kg*CO2e / MWh')
        def GWP_electricity(): # Cradle to gate
            if electricity_production.get():
                return GWP_economic.get() * dist.mean_electricity_price * 1000.
    
        @metric(name='Ethanol GWP', element='Displacement allocation', units='kg*CO2e / L')
        def GWP_ethanol_displacement(): # Cradle to gate
            if number in ethanol_configurations:
                GWP_material = sys.get_total_feeds_impact(GWP)
                GWP_electricity_production = GWP_characterization_factors['Electricity'] * electricity_production.get() * feedstock_consumption.get()
                GWP_coproducts = sys.get_total_products_impact(GWP)
                GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / y
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production - GWP_coproducts # kg CO2 eq. / y
                return GWP_total / (ethanol_production.get() * feedstock_consumption.get())
            else:
                return 0.
        
        @metric(name='Biodiesel GWP', element='Displacement allocation', units='kg*CO2e / L')
        def GWP_biodiesel_displacement(): # Cradle to gate
            if number in biodiesel_configurations:
                GWP_material = sys.get_total_feeds_impact(GWP)
                GWP_electricity_production = GWP_characterization_factors['Electricity'] * electricity_production.get() * feedstock_consumption.get()
                GWP_coproducts = sys.get_total_products_impact(GWP)
                GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / y
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production - GWP_coproducts # kg CO2 eq. / y
                return GWP_total / (biodiesel_production.get() * feedstock_consumption.get())
            else:
                return 0.
        
        # DO NOT DELETE:
        # import thermosteam as tmo
        # glycerol = tmo.Chemical('Glycerol')
        # ethanol = tmo.Chemical('Ethanol')
        # glycerol_GGE = 0.80 * (glycerol.LHV / glycerol.MW) / 121300 # 0.1059 GGE / kg crude-glycerol
        
        @metric(name='Biofuel GWP', element='Energy allocation', units='kg*CO2e / GGE')
        def GWP_biofuel_allocation(): # Cradle to gate
            GWP_material = sys.get_total_feeds_impact(GWP)
            GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / y
            GWP_total = GWP_material + GWP_emissions # kg CO2 eq. / y
            GGE_biodiesel_annual = (biodiesel_production.get() * feedstock_consumption.get()) / 0.9536 / L_per_gal
            GGE_ethanol_annual = (ethanol_production.get() * feedstock_consumption.get()) / 1.5 / L_per_gal
            GEE_electricity_production = max(-electricity() * 3600 / 114000, 0.) 
            GEE_crude_glycerol = crude_glycerol_flow() * 0.1059
            return GWP_total / (GGE_biodiesel_annual + GGE_ethanol_annual + GEE_electricity_production + GEE_crude_glycerol)
        
        @metric(name='Ethanol GWP', element='Energy allocation', units='kg*CO2e / L')
        def GWP_ethanol_allocation(): # Cradle to gate
            return GWP_biofuel_allocation.get() / 1.5 / L_per_gal
        
        @metric(name='Biodiesel GWP', element='Energy allocation', units='kg*CO2e / L')
        def GWP_biodiesel_allocation(): # Cradle to gate
            if number > 0:
                return GWP_biofuel_allocation.get() / 0.9536 / L_per_gal
            else:
                return 0.
        
        @metric(name='Crude-glycerol GWP', element='Energy allocation', units='kg*CO2e / kg')
        def GWP_crude_glycerol_allocation(): # Cradle to gate
            if number > 0:
                return GWP_biofuel_allocation.get() * 0.1059
            else:
                return 0.
    
        @metric(units='USD/MT')
        def MFPP_derivative():
            if number < 0:
                return 0.
            if self._derivative_disabled:
                return np.nan
            if agile:
                cane_mode.oil_content += 0.01
                sorghum_mode.oil_content += 0.01
            else:
                composition_specification.load_oil_content(composition_specification.oil + 0.01)
            sys.simulate()
            return MFPP.difference()
        
        @metric(units='L/MT')
        def biodiesel_production_derivative():
            if number < 0: return 0.
            if self._derivative_disabled: return np.nan
            return biodiesel_production.difference()
        
        @metric(units='L/MT')
        def ethanol_production_derivative():
            if number < 0: return 0.
            if self._derivative_disabled: return np.nan
            return ethanol_production.difference()
        
        @metric(units='kWh/MT')
        def electricity_production_derivative():
            if number < 0: return 0.
            if self._derivative_disabled: return np.nan
            return electricity_production.difference()
        
        @metric(units='cf/MT')
        def natural_gas_consumption_derivative():
            if number < 0: return 0.
            if self._derivative_disabled: return np.nan
            return natural_gas_consumption.difference()
        
        @metric(units='10^6*USD')
        def TCI_derivative():
            if number < 0: return 0.
            if self._derivative_disabled: return np.nan
            return TCI.difference()
        
        @metric(name='GWP derivative', element='Economic allocation', units='kg*CO2e / USD')
        def GWP_economic_derivative(): # Cradle to gate
            if number < 0: return 0.
            if self._derivative_disabled: return 0.
            return GWP_economic.difference()
    
        @metric(name='GWP derivative', element='Ethanol', units='kg*CO2e / L')
        def GWP_ethanol_derivative(): # Cradle to gate
            return GWP_economic_derivative.get() * get_GWP_mean_ethanol_price()
        
        @metric(name='GWP derivative', element='Biodiesel', units='kg*CO2e / L')
        def GWP_biodiesel_derivative(): # Cradle to gate
            if number > 0:
                return GWP_economic_derivative.get() * get_GWP_mean_biodiesel_price()
            else:
                return 0.
        
        @metric(name='GWP derivative', element='Crude glycerol', units='kg*CO2e / kg')
        def GWP_crude_glycerol_derivative(): # Cradle to gate
            if number > 0:
                return GWP_economic_derivative.get() * dist.mean_glycerol_price
            else:
                return 0.
        
        @metric(name='GWP derivative', element='Electricity', units='kg*CO2e / MWh')
        def GWP_electricity_derivative(): # Cradle to gate
            if electricity_production.get():
                return GWP_economic_derivative.get() * dist.mean_electricity_price * 1000.
            else:
                return 0.
    
        @metric(units='%')
        def ROI():
            return 100. * tea.ROI
        
        # def competitive_biomass_yield_objective(biomass_yield, target):
        #     self.dry_biomass_yield = biomass_yield
        #     F_mass = self.feedstock.F_mass
        #     self.update_feedstock()
        #     specifications = self.sys.specifications
        #     self.sys.specifications = []
        #     try:
        #         self.sys.simulate()
        #     finally:
        #         self.sys.specifications = specifications
        #     self.feedstock.F_mass = F_mass
        #     return 100. * tea.ROI - target
        
        # @metric(name='Competitive biomass yield', element='Feedstock', units='dry MT/ha')
        # def competitive_biomass_yield():
        #     # assert self.ROI_target is None
        #     if self.ROI_target is None: return np.nan
        #     if composition_specification.oil == 0: return self.baseline_dry_biomass_yield
        #     f = competitive_biomass_yield_objective
        #     x0 = 0.1 * self.baseline_dry_biomass_yield
        #     x1 = 3 * self.baseline_dry_biomass_yield
        #     args = (self.ROI_target,)
        #     if (y0:=f(x0, *args)) > 0. or (y1:=f(x1, *args)) < 0.:
        #         return np.nan
        #     else:
        #         return flx.IQ_interpolation(
        #             f, x0, x1, y0, y1, x=competitive_biomass_yield.get(evaluate=False),
        #             args=args, 
        #         )
        
        def NE_and_TCI_at_biomass_yield(dry_biomass_yield):
            self.dry_biomass_yield = dry_biomass_yield
            self.update_feedstock()
            sys.simulate()
            return tea.net_earnings, tea.TCI
        
        def competitive_biomass_yield_objective(biomass_yield, target, mb_NE, An_TCI):
            return 100. * linear_val(biomass_yield, mb_NE) / exponential_val(max(biomass_yield, 1), An_TCI) - target
        
        @metric(name='Competitive biomass yield', element='Feedstock', units='dry MT/ha')
        def competitive_biomass_yield():
            # assert self.ROI_target is None
            if self.ROI_target is None: return np.nan
            if composition_specification.oil == 0: return self.baseline_dry_biomass_yield
            x0 = self.dry_biomass_yield
            assert x0 < 100, "dry biomass yield over 100 dry MT / ha"
            x1 = 2 * x0
            NE0, TCI0 = tea.net_earnings, tea.TCI
            NE1, TCI1 = NE_and_TCI_at_biomass_yield(x1)
            self.dry_biomass_yield = x0 # Reset dry biomass yield
            mb_NE = linear_fit(x0, x1, NE0, NE1)
            An_TCI = exponential_fit(x0, x1, TCI0, TCI1)
            
            # DO NOT DELETE: This is to validate the functional approximation made here
            # from matplotlib import pyplot as plt
            # N = 10
            # xs = np.linspace(x0, x1, N)
            # NEs, TCIs = zip(*[NE_and_TCI_at_biomass_yield(x) for x in xs])
            # NEs = np.array(NEs)
            # TCIs = np.array(TCIs)
            # plt.figure()
            # plt.plot(xs, NEs, label='NE actual')
            # NEs_model = np.array([linear_val(x, mb_NE) for x in xs])
            # plt.scatter(xs, NEs_model, label='NE model')
            
            # plt.figure()
            # plt.plot(xs, TCIs, label='TCI actual')
            # TCIs_model = np.array([exponential_val(x, An_TCI) for x in xs])
            # plt.scatter(xs, TCIs_model, label='TCI model')
            # dNEs = (NEs - NEs_model)
            
            # TCIs_log = np.log(TCIs)
            # TCIs_model_log = np.log(TCIs_model)
            # dTCIs_log = TCIs_log - TCIs_model_log
            # TCI_log_SSR = (dTCIs_log * dTCIs_log).sum()
            # TCI_ave_log = np.mean(TCIs_log)
            # dTCIs_log_mean = (TCIs_log - TCI_ave_log)
            # TCI_log_SST = (dTCIs_log_mean * dTCIs_log_mean).sum()
            # R2_TCI_log = 1 - TCI_log_SSR / TCI_log_SST
            # print(R2_TCI_log) # 0.9995
            # assert R2_TCI_log > 0.999
            
            # NE_ave = np.mean(NEs)
            # dNEs = NEs - NEs_model
            # dNEs_mean = NEs - NE_ave
            # NE_SSR = (dNEs * dNEs).sum()
            # NE_SST = (dNEs_mean * dNEs_mean).sum()
            # R2_NE = 1 - NE_SSR / NE_SST
            
            # print(R2_NE) # 0.9999
            # assert R2_NE > 0.999
            
            f = competitive_biomass_yield_objective
            args = (self.ROI_target, mb_NE, An_TCI)
            value = flx.aitken_secant(
                f, x0, x1, args=args, 
            )
            return value
        
        @metric(name='Energy competitive biomass yield', element='Feedstock', units='dry MT/ha')
        def energy_competitive_biomass_yield():
            # assert self.net_energy_target is None
            if self.net_energy_target is None: return np.nan
            NEP_target = self.net_energy_target
            if composition_specification.oil == 0: return self.baseline_dry_biomass_yield
            NEP = net_energy_production.cache
            return self.baseline_dry_biomass_yield * NEP / NEP_target  
            
        @metric(units='%', name='Breakeven IRR')
        def IRR():
            finance_interest = tea.finance_interest
            if finance_interest: 
                IRR_original = tea.IRR
                tea.IRR = finance_interest
                financing = tea.NPV > 0 # Financing does not help otherwise
                IRR = tea.solve_IRR(financing=financing)
                tea.IRR = IRR_original
            else:
                IRR = tea.solve_IRR()
            return 100. * IRR
        
        # def competitive_microbial_oil_yield_objective(microbial_oil_yield, target):
        #     self.update_dry_biomass_yield(self.baseline_dry_biomass_yield)
        #     set_glucose_to_microbial_oil_yield.setter(microbial_oil_yield)
        #     set_xylose_to_microbial_oil_yield.setter(microbial_oil_yield)
        #     self.sys.simulate()
        #     return 100 * tea.ROI - target
        
        # @metric(name='Competitive microbial oil yield', element='Feedstock', units='wt. %')
        # def competitive_microbial_oil_yield():
        #     if self.ROI_target is None or self.microbial_oil_analysis_disactivated: return np.nan
        #     f = competitive_microbial_oil_yield_objective
        #     x0 = 20
        #     x1 = 95
        #     args = (self.ROI_target,)
        #     if (y0:=f(x0, *args)) > 0. or (y1:=f(x1, *args)) < 0.:
        #         return np.nan
        #     else:
        #         flx.IQ_interpolation(
        #             f, x0, x1, y0, y1, xtol=1e-2, ytol=1e-3, args=args,
        #         )
        #         if number in cellulosic_oil_configurations:
        #             return 100. * fermentor.cofermentation[0].product_yield(product='TriOlein', basis='wt')
        #         elif number in oil_configurations:
        #             return 100. * fermentor.fermentation_reaction[0].product_yield(product='TriOlein', basis='wt')
        
        # @metric(name='Energy competitive microbial oil yield', element='Feedstock', units='wt. %')
        # def energy_competitive_microbial_yield():
        #     if self.net_energy_target is None or self.microbial_oil_analysis_disactivated: return np.nan
        #     f = competitive_microbial_oil_yield_objective
        #     x0 = 20
        #     x1 = 95
        #     if (y0:=f(x0)) > 0. or (y1:=f(x1)) < 0.:
        #         return np.nan
        #     else:
        #         flx.IQ_interpolation(
        #             f, x0, x1, y0, y1, xtol=1e-2, ytol=1e-3, args=(self.net_energy_target,)
        #         )
        #         if number in cellulosic_oil_configurations:
        #             return 100. * fermentor.cofermentation[0].product_yield(product='TriOlein', basis='wt')
        #         elif number in oil_configurations:
        #             return 100. * fermentor.fermentation_reaction[0].product_yield(product='TriOlein', basis='wt')
        
        # Single point evaluation for detailed design results
        def set_baseline(p, x=None):
            if x is None:
                x = p.baseline
            else:
                p.baseline = x
            p.setter(x)
        
        if number in cellulosic_ethanol_configurations:
            set_baseline(set_sorghum_glucose_yield, 79)
            set_baseline(set_sorghum_xylose_yield, 86)
            set_baseline(set_cane_glucose_yield, 91.0)
            set_baseline(set_cane_xylose_yield, 97.5)
            set_baseline(set_glucose_to_ethanol_yield, 90)
            set_baseline(set_xylose_to_ethanol_yield, 42)
        set_baseline(set_cane_oil_content, self.default_baseline_oil_content)
        set_baseline(set_bagasse_oil_recovery, 70)
        set_baseline(set_juicing_oil_recovery, 60)
        set_baseline(set_microbial_oil_recovery)
        set_baseline(set_crude_oil_price) 
        set_baseline(set_advanced_ethanol_price) 
        set_baseline(set_cellulosic_ethanol_price) 
        set_baseline(set_biomass_based_diesel_price)
        set_baseline(set_cellulosic_based_diesel_price)
        set_baseline(set_crude_glycerol_price, dist.mean_glycerol_price)
        set_baseline(set_natural_gas_price)
        set_baseline(set_electricity_price)
        if number in cellulosic_ethanol_configurations:
            get_stream('ethanol').price = 0.789
        if number > 0:
            set_baseline(set_cane_PL_content, 10)
            set_baseline(set_cane_FFA_content, 10)
        
        for i in model._parameters: setattr(self, i.setter.__name__, i)
        for i in model._metrics: setattr(self, i.getter.__name__, i)
        self.sys = sys
        self.tea = tea
        self.model = model
        self.flowsheet = flowsheet
        self.unit_groups = unit_groups
        self.configuration = configuration
        self.composition_specification = composition_specification
        self.oil_extraction_specification = oil_extraction_specification
        self.available_land = self.baseline_available_land # ha
        self.dry_biomass_yield = self.baseline_dry_biomass_yield # MT / ha
        self.ROI_target = None
        self.net_energy_target = None
        self.microbial_oil_analysis_disactivated = True
        self.__dict__.update(flowsheet.to_dict())
        if feedstock_line: self.set_feedstock_line(feedstock_line)
        if cache is not None: cache[key] = self
        
        if case is not None:
            getattr(self, case)()
        
        # Avoid erros in Monte Carlo of microbial oil production with huge cell
        # mass production
        if number in cellulosic_configurations:
            flowsheet('SludgeCentrifuge').strict_moisture_content = False
            flowsheet(bst.PressureFilter).isplit[
                'Yeast', 'Glucan', 'Xylan', 'Arabinan',
                'Galactan', 'Lignin', 'Solids',
            ] = 1.
        
        try:
            PolishingFilter = flowsheet(bst.wastewater.high_rate.PolishingFilter)
        except:
            pass
        else:
            def adjust_system_convergence(system):
                system.converge_method = 'fixed-point'
                system.maxiter = 500
                system.molar_tolerance = 5
            
            PolishingFilter.recycle_system_hook = adjust_system_convergence
        
        ## Simulation
        if simulate:
            sys.simulate()
            if update_feedstock_price:
                feedstock.price = tea.solve_price(feedstock)
        
        ## Tests
        if feedstock_line is None:
            try:
                assert len(all_parameter_mockups) == len(model.parameters)
                assert len(all_metric_mockups) == len(model.metrics)
                for mockup, real in zip(all_parameter_mockups, model.parameters):
                    assert mockup.index == real.index
                for mockup, real in zip(all_metric_mockups, model.metrics):
                    assert mockup.index == real.index
            except:
                breakpoint()
        return self
    
    def oil_recovery(self):
        f = self.flowsheet
        tank = f(bst.BlendingTankWithSkimming)
        recovered_oil = tank.ins[0].F_mass
        microbial_oil = 0
        for reactor in f(bst.StirredTankReactor):
            microbial_oil += sum([i.imass['Oil'] for i in reactor.outs]) - sum([i.imass['Oil'] for i in reactor.ins])
        biomass_oil = self.feedstock.imass['Oil']
        total_oil = microbial_oil + biomass_oil
        return {
            '% Recovery': 100 * recovered_oil / total_oil,
            'Microbial fraction': microbial_oil / total_oil,
            'Biomass fraction': biomass_oil / total_oil,
            'Oil recovery [kg/s]': recovered_oil / 3600,
            'Total oil [kg/s]': total_oil / 3600,
        }
        

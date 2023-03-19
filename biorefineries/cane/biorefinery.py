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
from biosteam.utils import MockStream
import thermosteam as tmo
from biorefineries.cane.composition import (
    set_lipid_fraction as set_oil_fraction,
    set_line_composition_parameters,
)
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
import numpy as np
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

cellulosic_configurations = frozenset([-2, 2, 4, 6, 8])
biodiesel_configurations = frozenset([1, 2, 5, 6, 7, 8])
ethanol_configurations = frozenset([-2, -1, 1, 2, 3, 4])
actag_configurations = frozenset([9, 10])
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
    9: create_oilcane_to_biodiesel_and_actag_1g,
    10: create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation,
}
system_factory_options = {
    -1: dict(use_area_convention=True, pellet_bagasse=False, dry_bagasse=True),
    5: dict(fed_batch=True),
    6: dict(fed_batch=True),
    7: dict(fed_batch=False),
    8: dict(fed_batch=False),
}
area_names = {
    -1: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'CH&P', 'Utilities', 
         'HXN', 'Storage'],
    -2: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.',
         'Wastewater treatment', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    -3: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'CH&P', 'Utilities', 
         'HXN', 'Storage'],
    1: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'Oil ext.', 'Biod. prod.', 
        'CH&P', 'Utilities', 'HXN', 'Storage'],
    2: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.', 'Wastewater treatment',
        'Oil ext.', 'CH&P',  'Biod. prod.', 'Utilities', 'HXN', 'Storage'],
    3: ['Feedstock handling', 'Juicing', 'EtOH prod.', 'Oil ext.',  'CH&P', 
        'Utilities', 'HXN', 'Storage'],
    4: ['Feedstock handling', 'Juicing', 'Pretreatment', 'EtOH prod.', 
        'Wastewater treatment', 'Oil ext.', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    5: ['Feedstock handling', 'Juicing', 'Oil prod. & ext.', 'Biod. prod.', 
        'CH&P', 'Utilities', 'HXN', 'Storage'],
    6: ['Feedstock handling',  'Juicing', 'Pretreatment', 'Oil prod. & ext.',
        'Wastewater treatment', 'Biod. prod.', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    # 9: ['Feedstock handling', 'Juicing', 'TAG prod.', 'Oil ext.', 'AcTAG sep.', 
    #     'Biod. prod.', 'CH&P', 'Utilities', 'HXN', 'Storage'],
    # 10: ['Feedstock handling', 'Juicing', 'Pretreatment', 'TAG prod.', 'AcTAG sep.', 
    #      'Wastewater treatment', 'Oil ext.', 'CH&P',  'Biod. prod.', 'Utilities', 'HXN', 'Storage'],
}
area_names[7] = area_names[5]
area_names[8] = area_names[6]

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

def YRCP2023():
    Biorefinery.default_conversion_performance_distribution = 'shortterm'
    Biorefinery.default_prices_correleted_to_crude_oil = True
    Biorefinery.default_year = 2023

class Biorefinery:
    cache = {}
    dry_biomass_yield = None # MT / ha
    baseline_dry_biomass_yield = 25.62 # MT / ha 
    set_feedstock_line = set_line_composition_parameters
    default_prices_correleted_to_crude_oil = False
    default_conversion_performance_distribution = 'longterm'
    default_year = 2022
    baseline_moisture_content = 0.70
    
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
    def wet_biomass_yield(self):
        dry_biomass_yield = self.dry_biomass_yield
        if dry_biomass_yield is None: return None
        moisture_content = self.feedstock.get_mass_fraction('Water')
        dry_fraction = 1 - moisture_content
        return dry_biomass_yield / dry_fraction
    
    def update_dry_biomass_yield(self, dry_biomass_yield):
        self.dry_biomass_yield = dry_biomass_yield
        # Set preliminary feedstock price assuming 10% is due to transportation
        # and 90% is based on productivity (a function of height)
        feedstock = self.feedstock
        baseline_dry_biomass_yield = self.baseline_dry_biomass_yield
        baseline_moisture_content = self.baseline_moisture_content
        f_yield = baseline_dry_biomass_yield / dry_biomass_yield
        moisture_content = feedstock.get_mass_fraction('Water')
        dry_content = 1. - moisture_content
        baseline_dry_content = 1. - baseline_moisture_content
        f_dry = dry_content / baseline_dry_content
        f = f_yield * f_dry
        feedstock.price = (0.9 * f + 0.10) * 0.035
        # Make the same assumption for the cane feedstock (subject to change).
        feedstock.set_CF(
            GWP, GWP_characterization_factors['sugarcane'] * (0.9 * f + 0.10)
        )
        
    def __new__(cls, name, chemicals=None, reduce_chemicals=False, 
                 avoid_natural_gas=True, conversion_performance_distribution=None,
                 year=None, cache=cache, feedstock_line=None,
                 prices_correleted_to_crude_oil=None):
        if year is None: year = cls.default_year
        if conversion_performance_distribution is None: 
            conversion_performance_distribution = cls.default_conversion_performance_distribution
        else:
            conversion_performance_distribution = conversion_performance_distribution.replace(' ', '').replace('-', '').lower()
        if prices_correleted_to_crude_oil is None:
            prices_correleted_to_crude_oil = cls.default_prices_correleted_to_crude_oil
        number, agile, energycane = configuration = parse_configuration(name)
        key = (number, agile, energycane, conversion_performance_distribution, year)
        if cache and key in cache: 
            return cache[key]
        else:
            self = super().__new__(cls)
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
            self.cane_sys = cane_sys = system_factories[number](operating_hours=24 * 200, **system_factory_options.get(number, {}))
        else:
            raise NotImplementedError(number)
        cane_sys.set_tolerance(rmol=1e-5, mol=1e-2, subsystems=True, subfactor=1.5)
        
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
            # Default composition as equimolar
        )
        
        if energycane:
            feedstock.reset_flow(**energycane_dct, total_flow=feedstock.F_mass, units='kg/hr')
            feedstock.ID = 'energycane'
        
        BT = flowsheet(bst.BoilerTurbogenerator)
        if number not in cellulosic_configurations: BT.boiler_efficiency = 0.89
        
        ## Unit groups
        areas = area_names[number]
        rename_storage_units(cane_sys.units, len(areas) * 100)
        unit_groups = UnitGroup.group_by_area(cane_sys.units)
        for i, j in zip(unit_groups, areas): i.name = j
        for i in unit_groups: i.autofill_metrics(shorthand=True)
        
        for HXN_group in unit_groups:
            if HXN_group.name == 'HXN':
                HXN_group.filter_savings = False # Allow negative values in heat utilities
                HXN = HXN_group.units[0]
                assert isinstance(HXN, bst.HeatExchangerNetwork)
        HXN.raise_energy_balance_error = True
        HXN.vle_quenched_streams = False
        
        ## Split bagasse to boiler meet energy demand
        if avoid_natural_gas:
            for splitter in flowsheet.unit:
                if getattr(splitter, 'isbagasse_splitter', False):
                    self.bagasse_splitter = splitter
                    minimum_fraction_burned = 0
                    maximum_fraction_burned = 0.7
                    self.recycle_data = recycle_data = {}
                    @cane_sys.add_bounded_numerical_specification(
                        x0=minimum_fraction_burned, x1=maximum_fraction_burned, 
                        xtol=1e-4, ytol=100, args=(splitter,)
                    )
                    def adjust_bagasse_to_boiler(fraction_burned, splitter):
                        # Returns energy consumption at given fraction processed (not sent to boiler).
                        splitter.split[:] = 1 - fraction_burned
                        operation_mode = getattr(sys, 'active_operation_mode', None)
                        if fraction_burned in (minimum_fraction_burned, maximum_fraction_burned):
                            key = (operation_mode, fraction_burned)
                        else:
                            key = (operation_mode, 'last')
                        if key in recycle_data: 
                            material_data = recycle_data[key]
                        else:
                            recycle_data[key] = material_data = cane_sys.get_material_data()
                        try:
                            cane_sys.simulate(material_data=material_data, update_material_data=True)
                        except:
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
            
        if abs(number) in cellulosic_configurations:
            prs = flowsheet(bst.PretreatmentReactorSystem)
            saccharification = flowsheet(bst.Saccharification)
            seed_train = flowsheet(bst.SeedTrain)
            fermentor = flowsheet(bst.CoFermentation)
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
            
            if energycane:
                self.cane_mode = cane_mode = sys.operation_mode(cane_sys,
                    operating_hours=180*24, oil_content=0.02, feedstock=feedstock.copy(),
                    z_mass_carbs_baseline=0.091,
                    z_mass_solids_baseline=0., 
                    z_mass_ash_baseline=0.028,
                    z_mass_water_baseline=0.60,
                    glucose_yield=85, xylose_yield=65, 
                    FFA_content=0.10, PL_content=0.10
                )
            else:
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
        
        tea.income_tax = 0.21 # Davis et al. 2018; https://www.nrel.gov/docs/fy19osti/71949.pdf
        
        ## Specifications for analysis
        self.composition_specification = composition_specification = CaneCompositionSpecification(feedstock)
        if number < 0:
            oil_extraction_specification = MockExtractionSpecification()
        else:
            crushing_mill = flowsheet(bst.CrushingMill) # Separates bagasse
            if number in cellulosic_configurations:
                pressure_filter = flowsheet(bst.PressureFilter) # Separates lignin
            else:
                pressure_filter  = None
            cellmass_centrifuge = flowsheet(bst.SolidsCentrifuge) # Separates cell mass
            microbial_oil_recovery = 0.9 # Baseline
            oil_extraction_specification = OilExtractionSpecification(
                sys, crushing_mill, pressure_filter, cellmass_centrifuge, microbial_oil_recovery, 
            )
        
        ## Account for cellulosic vs advanced RINs
        if number in cellulosic_ethanol_configurations:
            # Note that GREET cellulosic ethanol from corn stover results in a 
            # GWP of 0.41 kg CO2-eq / L-ethanol. So the cellulosic ethanol from
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
                # GWP of 0.41 kg CO2-eq / L-ethanol. So the cellulosic ethanol from
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
        HXN.replace_unit_heat_utilities = True
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
        dist = get_price_distributions_module(year)
        
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
        
        if conversion_performance_distribution == 'shortterm':
            performance = potential_shortterm_gain_in_performance
            oil_content_range = [2, 5]
        elif conversion_performance_distribution == 'longterm':
            performance = potential_longterm_gain_in_performance
            oil_content_range = [5, 15]
        else:
            raise ValueError("`conversion_performance_distribution` must be either 'longterm' or 'shortterm")
            
        @performance(60, 90, units='%', kind='coupled')
        def set_crushing_mill_oil_recovery(oil_recovery):
            oil_extraction_specification.load_crushing_mill_oil_recovery(oil_recovery / 100.)
        
        @performance(90.0, 95, units='%', kind='coupled')
        def set_microbial_oil_recovery(microbial_oil_recovery):
            oil_extraction_specification.load_microbial_oil_recovery(microbial_oil_recovery / 100.)
        
        @performance(70.0, 90, units='%', kind='coupled')
        def set_bagasse_oil_recovery(bagasse_oil_recovery):
            oil_extraction_specification.load_bagasse_oil_recovery(bagasse_oil_recovery / 100.)
    
        # Baseline from Huang's 2016 paper, but distribution more in line with Florida sugarcane harvesting (3-5 months)
        @uniform(4 * 30, 6 * 30, units='day/yr', baseline=180)
        def set_cane_operating_days(cane_operating_days):
            if agile:
                cane_mode.operating_hours = cane_operating_days * 24
            else:
                tea.operating_days = cane_operating_days
    
        # From Ed Cahoon and Huang 2017
        @uniform(30, 60, units='day/yr', baseline=45)
        def set_sorghum_operating_days(sorghum_operating_days):
            if agile: sorghum_mode.operating_hours = sorghum_operating_days * 24
        
        @default(1600000, units='MT/yr', kind='isolated')
        def set_annual_crushing_capacity(annual_crushing_capacity):
            sys.rescale(feedstock, kg_per_MT * annual_crushing_capacity / tea.operating_hours / feedstock.F_mass)
    
        @parameter(distribution=dist.copd, element='Crude oil', baseline=dist.mcop, units='USD/barrel')
        def set_crude_oil_price(price):
            self.crude_oil_price = price
        
        if prices_correleted_to_crude_oil:
            @parameter(distribution=dist.cepd_offset, element=cellulosic_ethanol, baseline=0., units='USD/L')
            def set_cellulosic_ethanol_price(price): 
                cellulosic_ethanol.price = (dist.f_cep * self.crude_oil_price + price) * ethanol_L_per_kg
                
            @parameter(distribution=dist.aepd_offset, element=advanced_ethanol, baseline=0., units='USD/L')
            def set_advanced_ethanol_price(price): 
                advanced_ethanol.price =  (dist.f_aep * self.crude_oil_price + price) * ethanol_L_per_kg
                
            # USDA ERS historical price data
            @parameter(distribution=dist.bpd_offset, element=biomass_based_diesel, units='USD/L', baseline=0.)
            def set_biomass_based_diesel_price(price):
                biomass_based_diesel.price =  (dist.f_bp * self.crude_oil_price + price) * biodiesel_L_per_kg
        
            @parameter(distribution=dist.cbpd_offset, element=cellulosic_based_diesel, units='USD/L', baseline=0.)
            def set_cellulosic_based_diesel_price(price):
                cellulosic_based_diesel.price =  (dist.f_cbp * self.crude_oil_price + price) * biodiesel_L_per_kg
        
            # https://www.eia.gov/energyexplained/natural-gas/prices.php
            @parameter(distribution=dist.ngpd_offset, element=s.natural_gas, units='USD/m3', baseline=0.)
            def set_natural_gas_price(price): 
                BT.natural_gas_price =  (dist.f_ngp * self.crude_oil_price + price) * V_ng
        
            @parameter(distribution=dist.elecpd_offset, units='USD/kWhr', baseline=0.)
            def set_electricity_price(electricity_price): 
                bst.PowerUtility.price = dist.f_elecp * self.crude_oil_price + electricity_price
                
        else:
            # USDA ERS historical price data with EPA RIN prices
            @parameter(distribution=dist.cepd, element=cellulosic_ethanol, 
                       baseline=dist.mcep, units='USD/L')
            def set_cellulosic_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
                cellulosic_ethanol.price = price * ethanol_L_per_kg
                
            @parameter(distribution=dist.aepd, element=advanced_ethanol, 
                       baseline=dist.maep, units='USD/L')
            def set_advanced_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
                advanced_ethanol.price = price * ethanol_L_per_kg
                
            # USDA ERS historical price data
            @parameter(distribution=dist.bpd, element=biomass_based_diesel, units='USD/L', baseline=dist.mbp)
            def set_biomass_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
                biomass_based_diesel.price = price * biodiesel_L_per_kg
        
            @parameter(distribution=dist.cbpd, element=cellulosic_based_diesel, units='USD/L', baseline=dist.mcbp)
            def set_cellulosic_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
                cellulosic_based_diesel.price = price * biodiesel_L_per_kg
        
            # https://www.eia.gov/energyexplained/natural-gas/prices.php
            @parameter(distribution=dist.natural_gas_price_distribution, element=s.natural_gas, units='USD/m3',
                       baseline=4.73 * 35.3146667/1e3)
            def set_natural_gas_price(price): 
                BT.natural_gas_price = price * V_ng
        
            @parameter(distribution=dist.electricity_price_distribution, units='USD/kWhr', baseline=dist.mean_electricity_price)
            def set_electricity_price(electricity_price): 
                bst.PowerUtility.price = electricity_price
        
        # 10% is suggested for waste reducing, but 15% is suggested for investment
        @uniform(10., 15., units='%', baseline=10)
        def set_IRR(IRR):
            tea.IRR = IRR / 100.
        
        @uniform(0.10, 0.22, units='USD/kg', element=crude_glycerol)
        def set_crude_glycerol_price(price):
            crude_glycerol.price = price
    
        @default(0.65, units='USD/kg', element=pure_glycerine)
        def set_pure_glycerol_price(price):
            pure_glycerine.price = price
        
        @default(72, units='hr', element='Saccharification')
        def set_saccharification_reaction_time(reaction_time):
            if abs(number) in cellulosic_configurations: saccharification.tau = reaction_time
        
        @default(0.212, units='USD/kg', element=cellulase)
        def set_cellulase_price(price):
            cellulase.price = price
    
        if abs(number) in cellulosic_configurations:
            cellulase_mixer, = [i for i in flowsheet.unit if hasattr(i, 'enzyme_loading')]
       
        @default(0.02, units='wt. % cellulose', element='cellulase', kind='coupled')
        def set_cellulase_loading(cellulase_loading):
            if abs(number) in cellulosic_configurations: cellulase_mixer.enzyme_loading = cellulase_loading
        
        @default(PRS_cost_item.cost, units='million USD', element='Pretreatment reactor system')
        def set_reactor_base_cost(base_cost):
            PRS_cost_item.cost = base_cost
        
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
    
        @performance(0.951, 1.902, units='g/L', element='Cofermentation')
        def set_cofermentation_ethanol_productivity(ethanol_productivity):
            if number in cellulosic_ethanol_configurations: fermentor.productivity = ethanol_productivity
    
        fed_batch = system_factory_options.get(number, {}).get('fed_batch')
        @performance(60 if fed_batch else 49.5, 95, units='%', element='Cofermenation', kind='coupled')
        def set_glucose_to_microbial_oil_yield(glucose_to_microbial_oil_yield):
            glucose_to_microbial_oil_yield *= 0.01
            cell_growth = min(0.99 - glucose_to_microbial_oil_yield, 0.4 * glucose_to_microbial_oil_yield) # Almost all the rest goes towards cell mass
            if number in cellulosic_oil_configurations:
                seed_train.reactions.X[0] = fermentor.cofermentation.X[0] = glucose_to_microbial_oil_yield 
                seed_train.reactions.X[2] = fermentor.cofermentation.X[2] = cell_growth
            elif number in oil_configurations:
                fermentor.fermentation_reaction.X[0] = glucose_to_microbial_oil_yield
                fermentor.cell_growth_reaction.X[0] = cell_growth
        
        @performance(60 if fed_batch else 49.5, 95, units='%', element='Cofermenation', kind='coupled')
        def set_xylose_to_microbial_oil_yield(xylose_to_microbial_oil_yield):
            if number in cellulosic_oil_configurations:
                xylose_to_microbial_oil_yield *= 0.01
                cell_growth = min(0.99 - xylose_to_microbial_oil_yield, 0.4 * xylose_to_microbial_oil_yield)
                seed_train.reactions.X[1] = fermentor.cofermentation.X[1] = xylose_to_microbial_oil_yield 
                seed_train.reactions.X[3] = fermentor.cofermentation.X[3] = cell_growth # Almost all the rest goes towards cell mass
    
        @performance(89.4 if fed_batch else 27.4, 137, units='g/L', element='Cofermentation', kind='coupled')
        def set_cofermentation_microbial_oil_titer(microbial_oil_titer):
            if number in cellulosic_oil_configurations: fermentor.titer = microbial_oil_titer
    
        @performance(0.61 if fed_batch else 0.31, 1.902, units='g/L', element='Cofermentation')
        def set_cofermentation_microbial_oil_productivity(microbial_oil_productivity):
            if number in cellulosic_oil_configurations: fermentor.productivity = microbial_oil_productivity
    
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
            elif energycane:
                composition_specification.load_oil_content(
                    cane_oil_content / 100.,
                    z_mass_carbs_baseline=0.091,
                    z_mass_solids_baseline=0., 
                    z_mass_ash_baseline=0.028,
                    z_mass_water_baseline=0.60,
                )
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
                 element=feedstock, units='kg*CO2-eq/kg')
        def set_feedstock_GWP(value):
            if number > 0 and self.dry_biomass_yield is None:
                feedstock.characterization_factors[GWP] = value
        
        @default_gwp(methanol.characterization_factors[GWP], name='GWP', 
                     element=methanol, units='kg*CO2-eq/kg')
        def set_methanol_GWP(value):
            methanol.characterization_factors[GWP] = value
        
        # @default(crude_glycerol.characterization_factors[GWP], name='GWP', 
        #          element=crude_glycerol, units='kg*CO2-eq/kg')
        # def set_crude_glycerol_GWP(value):
        #     crude_glycerol.characterization_factors[GWP] = value
        
        @default_gwp(pure_glycerine.characterization_factors[GWP], name='GWP', 
                     element=pure_glycerine, units='kg*CO2-eq/kg')
        def set_pure_glycerine_GWP(value):
            pure_glycerine.characterization_factors[GWP] = value
        
        @default_gwp(cellulase.characterization_factors[GWP], name='GWP', 
                     element=cellulase, units='kg*CO2-eq/kg')
        def set_cellulase_GWP(value):
            cellulase.characterization_factors[GWP] = value * 0.02
        
        @default_gwp(natural_gas.characterization_factors[GWP], name='GWP', 
                     element=natural_gas, units='kg*CO2-eq/kg')
        def set_natural_gas_GWP(value):
            for ng in natural_gas_streams:
                ng.characterization_factors[GWP] = value
        
        @uniform(21., 28., baseline=21., name='Income tax', units='%')
        def set_income_tax(income_tax):
            tea.income_tax = income_tax * 0.01
        
        if agile:
            feedstock_flow = lambda: sys.flow_rates[feedstock] / kg_per_MT # MT / yr
            biodiesel_flow = lambda: (sys.flow_rates.get(cellulosic_based_diesel, 0.) + sys.flow_rates.get(biomass_based_diesel, 0.)) * biodiesel_L_per_kg # L / yr
            ethanol_flow = lambda: (sys.flow_rates.get(cellulosic_ethanol, 0.) + sys.flow_rates.get(advanced_ethanol, 0.)) * ethanol_L_per_kg # L / yr
            natural_gas_flow = lambda: sum([sys.flow_rates[i] for i in natural_gas_streams]) * V_ng # m3 / yr
            crude_glycerol_flow = lambda: sys.flow_rates.get(crude_glycerol, 0.) # kg / yr
            
            @sys.operation_metric(annualize=True)
            def direct_nonbiogenic_emissions(mode):
                return sum([i.F_mol for i in natural_gas_streams]) * chemicals.CO2.MW
            
        else:
            feedstock_flow = lambda: sys.operating_hours * feedstock.F_mass / kg_per_MT # MT / yr
            biodiesel_flow = lambda: sys.operating_hours * biodiesel.F_mass * biodiesel_L_per_kg # L / yr
            ethanol_flow = lambda: sys.operating_hours * ethanol.F_mass * ethanol_L_per_kg # L / yr
            crude_glycerol_flow = lambda: sys.operating_hours * crude_glycerol.F_mass # kg / yr
            natural_gas_flow = lambda: sum([i.F_mass for i in natural_gas_streams]) * sys.operating_hours * V_ng # m3 / yr
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
        
        @metric(units='MT/yr')
        def feedstock_consumption():
            return feedstock_flow()
        
        @metric(units='L/MT')
        def biodiesel_production():
            return biodiesel_flow() / feedstock_consumption.get()
        
        @metric(units='L/ha')
        def biodiesel_yield():
            if self.dry_biomass_yield is None: 
                return None
            else:
                return biodiesel_production.get() * self.wet_biomass_yield
        
        @metric(units='L/MT')
        def ethanol_production():
            return ethanol_flow() / feedstock_consumption.get()
        
        @metric(units='kWhr/MT')
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
    
        @metric(name='GWP', element='Economic allocation', units='kg*CO2*eq / USD')
        def GWP_economic(): # Cradle to gate
            GWP_material = sys.get_total_feeds_impact(GWP) # kg CO2 eq. / yr
            GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
            sales = (
                biodiesel_flow() * get_GWP_mean_biodiesel_price()
                + ethanol_flow() * get_GWP_mean_ethanol_price()
                + crude_glycerol_flow() * dist.mean_glycerol_price
                + max(-electricity(), 0) * dist.mean_electricity_price
            )
            return (GWP_material + GWP_emissions) / sales
    
        @metric(name='Ethanol GWP', element='Economic allocation', units='kg*CO2*eq / L')
        def GWP_ethanol(): # Cradle to gate
            return GWP_economic.get() * get_GWP_mean_ethanol_price()
        
        @metric(name='Biodiesel GWP', element='Economic allocation', units='kg*CO2*eq / L')
        def GWP_biodiesel(): # Cradle to gate
            if number > 0:
                return GWP_economic.get() * get_GWP_mean_biodiesel_price()
            else:
                return 0.
        
        @metric(name='Crude glycerol GWP', element='Economic allocation', units='kg*CO2*eq / kg')
        def GWP_crude_glycerol(): # Cradle to gate
            if number in biodiesel_configurations:
                return GWP_economic.get() * dist.mean_glycerol_price
            else:
                return 0.
        
        @metric(name='Electricity GWP', element='Economic allocation', units='kg*CO2*eq / MWhr')
        def GWP_electricity(): # Cradle to gate
            if electricity_production.get():
                return GWP_economic.get() * dist.mean_electricity_price * 1000.
    
        @metric(name='Ethanol GWP', element='Displacement allocation', units='kg*CO2*eq / L')
        def GWP_ethanol_displacement(): # Cradle to gate
            if number in ethanol_configurations:
                GWP_material = sys.get_total_feeds_impact(GWP)
                GWP_electricity_production = GWP_characterization_factors['Electricity'] * electricity_production.get() * feedstock_consumption.get()
                GWP_coproducts = sys.get_total_products_impact(GWP)
                GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production - GWP_coproducts # kg CO2 eq. / yr
                return GWP_total / (ethanol_production.get() * feedstock_consumption.get())
            else:
                return 0.
        
        @metric(name='Biodiesel GWP', element='Displacement allocation', units='kg*CO2*eq / L')
        def GWP_biodiesel_displacement(): # Cradle to gate
            if number in biodiesel_configurations:
                GWP_material = sys.get_total_feeds_impact(GWP)
                GWP_electricity_production = GWP_characterization_factors['Electricity'] * electricity_production.get() * feedstock_consumption.get()
                GWP_coproducts = sys.get_total_products_impact(GWP)
                GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
                GWP_total = GWP_material + GWP_emissions - GWP_electricity_production - GWP_coproducts # kg CO2 eq. / yr
                return GWP_total / (biodiesel_production.get() * feedstock_consumption.get())
            else:
                return 0.
        
        # DO NOT DELETE:
        # import thermosteam as tmo
        # glycerol = tmo.Chemical('Glycerol')
        # ethanol = tmo.Chemical('Ethanol')
        # glycerol_GGE = 0.80 * (glycerol.LHV / glycerol.MW) / 121300 # 0.1059 GGE / kg crude-glycerol
        
        @metric(name='Biofuel GWP', element='Energy allocation', units='kg*CO2*eq / GGE')
        def GWP_biofuel_allocation(): # Cradle to gate
            GWP_material = sys.get_total_feeds_impact(GWP)
            GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
            GWP_total = GWP_material + GWP_emissions # kg CO2 eq. / yr
            GGE_biodiesel_annual = (biodiesel_production.get() * feedstock_consumption.get()) / 0.9536 / L_per_gal
            GGE_ethanol_annual = (ethanol_production.get() * feedstock_consumption.get()) / 1.5 / L_per_gal
            GEE_electricity_production = max(-electricity() * 3600 / 114000, 0.) 
            GEE_crude_glycerol = crude_glycerol_flow() * 0.1059
            return GWP_total / (GGE_biodiesel_annual + GGE_ethanol_annual + GEE_electricity_production + GEE_crude_glycerol)
        
        @metric(name='Ethanol GWP', element='Energy allocation', units='kg*CO2*eq / L')
        def GWP_ethanol_allocation(): # Cradle to gate
            return GWP_biofuel_allocation.get() / 1.5 / L_per_gal
        
        @metric(name='Biodiesel GWP', element='Energy allocation', units='kg*CO2*eq / L')
        def GWP_biodiesel_allocation(): # Cradle to gate
            if number > 0:
                return GWP_biofuel_allocation.get() / 0.9536 / L_per_gal
            else:
                return 0.
        
        @metric(name='Crude-glycerol GWP', element='Energy allocation', units='kg*CO2*eq / kg')
        def GWP_crude_glycerol_allocation(): # Cradle to gate
            if number > 0:
                return GWP_biofuel_allocation.get() * 0.1059
            else:
                return 0.
    
        @metric(units='%')
        def IRR():
            return 100. * tea.solve_IRR()
    
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
        
        @metric(units='kWhr/MT')
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
        
        @metric(name='GWP derivative', element='Economic allocation', units='kg*CO2*eq / USD')
        def GWP_economic_derivative(): # Cradle to gate
            if number < 0: return 0.
            if self._derivative_disabled: return 0.
            return GWP_economic.difference()
    
        @metric(name='Ethanol GWP derivative', element='Ethanol', units='kg*CO2*eq / L')
        def GWP_ethanol_derivative(): # Cradle to gate
            return GWP_economic_derivative.get() * get_GWP_mean_ethanol_price()
        
        @metric(name='Biodiesel GWP derivative', element='Biodiesel', units='kg*CO2*eq / L')
        def GWP_biodiesel_derivative(): # Cradle to gate
            if number > 0:
                return GWP_economic_derivative.get() * get_GWP_mean_biodiesel_price()
            else:
                return 0.
        
        @metric(name='Crude glycerol GWP derivative', element='Crude glycerol', units='kg*CO2*eq / kg')
        def GWP_crude_glycerol_derivative(): # Cradle to gate
            if number > 0:
                return GWP_economic_derivative.get() * dist.mean_glycerol_price
            else:
                return 0.
        
        @metric(name='Electricity GWP derivative', element='Electricity', units='kg*CO2*eq / MWhr')
        def GWP_electricity_derivative(): # Cradle to gate
            if electricity_production.get():
                return GWP_economic_derivative.get() * dist.mean_electricity_price * 1000.
            else:
                return 0.
    
        @metric(units='%')
        def ROI():
            if self.dry_biomass_yield is None: return np.nan
            return 100. * tea.ROI
        
        def competitive_biomass_yield_objective(biomass_yield, target):
            self.update_dry_biomass_yield(biomass_yield)
            return 100. * tea.ROI - target
        
        @metric(name='Competitive biomass yield', element='Feedstock', units='dry MT/ha')
        def competitive_biomass_yield():
            if self.ROI_target is None: return np.nan
            if composition_specification.oil == 0: return self.baseline_dry_biomass_yield
            f = competitive_biomass_yield_objective
            x0 = 0.1 * self.baseline_dry_biomass_yield
            x1 = 3 * self.baseline_dry_biomass_yield
            args = (self.ROI_target,)
            if (y0:=f(x0, *args)) > 0. or (y1:=f(x1, *args)) < 0.:
                return np.nan
            else:
                return flx.IQ_interpolation(
                    f, x0, x1, y0, y1, args=args
                )
        
        @metric(name='Energy competitive biomass yield', element='Feedstock', units='dry MT/ha')
        def energy_competitive_biomass_yield():
            if self.net_energy_target is None: return np.nan
            NEP_target = self.net_energy_target
            if composition_specification.oil == 0: return self.baseline_dry_biomass_yield
            NEP = net_energy_production.cache
            return self.baseline_dry_biomass_yield * NEP / NEP_target  
            
        def competitive_microbial_oil_yield_objective(microbial_oil_yield, target):
            self.update_dry_biomass_yield(self.baseline_dry_biomass_yield)
            set_glucose_to_microbial_oil_yield.setter(microbial_oil_yield)
            set_xylose_to_microbial_oil_yield.setter(microbial_oil_yield)
            self.sys.simulate()
            return 100 * tea.ROI - target
        
        @metric(name='Competitive microbial oil yield', element='Feedstock', units='wt. %')
        def competitive_microbial_oil_yield():
            if self.ROI_target is None or self.microbial_oil_analysis_disactivated: return np.nan
            f = competitive_microbial_oil_yield_objective
            x0 = 20
            x1 = 95
            args = (self.ROI_target,)
            if (y0:=f(x0, *args)) > 0. or (y1:=f(x1, *args)) < 0.:
                return np.nan
            else:
                flx.IQ_interpolation(
                    f, x0, x1, y0, y1, xtol=1e-2, ytol=1e-3, args=args,
                )
                if number in cellulosic_oil_configurations:
                    return 100. * fermentor.cofermentation[0].product_yield(product='TriOlein', basis='wt')
                elif number in oil_configurations:
                    return 100. * fermentor.fermentation_reaction[0].product_yield(product='TriOlein', basis='wt')
        
        @metric(name='Energy competitive microbial oil yield', element='Feedstock', units='wt. %')
        def energy_competitive_microbial_yield():
            if self.net_energy_target is None or self.microbial_oil_analysis_disactivated: return np.nan
            f = competitive_microbial_oil_yield_objective
            x0 = 20
            x1 = 95
            if (y0:=f(x0)) > 0. or (y1:=f(x1)) < 0.:
                return np.nan
            else:
                flx.IQ_interpolation(
                    f, x0, x1, y0, y1, xtol=1e-2, ytol=1e-3, args=(self.net_energy_target,)
                )
                if number in cellulosic_oil_configurations:
                    return 100. * fermentor.cofermentation[0].product_yield(product='TriOlein', basis='wt')
                elif number in oil_configurations:
                    return 100. * fermentor.fermentation_reaction[0].product_yield(product='TriOlein', basis='wt')
        
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
        if energycane:
            set_baseline(set_cane_oil_content, 2)
        else:
            set_baseline(set_cane_oil_content, 10)
        set_baseline(set_bagasse_oil_recovery, 70)
        set_baseline(set_crushing_mill_oil_recovery, 60)
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
        self.ROI_target = None
        self.net_energy_target = None
        self.microbial_oil_analysis_disactivated = True
        self.__dict__.update(flowsheet.to_dict())
        if feedstock_line: self.set_feedstock_line(feedstock_line)
        if cache is not None: cache[key] = self
        
        # Avoid erros in Monte Carlo of microbial oil production with huge cell
        # mass production
        if number in cellulosic_configurations:
            flowsheet('SludgeCentrifuge').strict_moisture_content = False
        
        ## Simulation
        sys.simulate()
        feedstock.price = tea.solve_price(feedstock)
        if reduce_chemicals: 
            cane_sys.reduce_chemicals()
            cane_sys._load_stream_links()
            sys.simulate()
        return self

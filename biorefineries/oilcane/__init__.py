# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import (
    units,
    _process_settings,
    _chemicals,
    systems,
    _tea,
    utils,
    _distributions,
    _evaluation,
    _uncertainty_plots,
    _contour_plots,
    _lca_characterization_factors,
    _load_data,
    _parse_configuration,
)
from .units import *
from ._process_settings import *
from ._chemicals import *
from ._contour_plots import *
from .systems import *
from ._tea import *
from .utils import *
from ._distributions import *
from ._evaluation import *
from ._uncertainty_plots import *
from ._contour_plots import *
from ._lca_characterization_factors import *
from ._load_data import *
from ._parse_configuration import *
from ._feature_mockups import *

__all__ = (
    *units.__all__,
    *_process_settings.__all__,
    *_chemicals.__all__,
    *systems.__all__,
    *_tea.__all__,
    *utils.__all__,
    *_distributions.__all__,
    *_evaluation.__all__,
    *_uncertainty_plots.__all__,
    *_contour_plots.__all__,
    *_lca_characterization_factors.__all__,
    *_load_data.__all__,
    *_parse_configuration.__all__,
    'sys',
    'tea', 
    'flowsheet',
)

import biosteam as bst
from biosteam.utils import MockStream
import thermosteam as tmo
from biorefineries.sugarcane import create_sugarcane_to_ethanol_system, create_sugarcane_to_sugar_and_ethanol_system
from biorefineries.cane.composition import (
    set_lipid_fraction as set_oil_fraction,
    get_composition,
)
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
import numpy as np
from ._process_settings import load_process_settings
from ._chemicals import create_chemicals
from .systems import (
    create_oilcane_to_biodiesel_and_ethanol_1g,
    create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    create_sugarcane_to_ethanol_combined_1_and_2g,  
    create_oilcane_to_crude_oil_and_ethanol_1g,
    create_oilcane_to_crude_oil_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    create_oilcane_to_biodiesel_1g,
    create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation
)
from ._parse_configuration import (
    parse_configuration,
    format_configuration,
)
from ._tea import (
    create_tea,
)
from .utils import (
    OilExtractionSpecification,
    CaneCompositionSpecification,
    MockExtractionSpecification,
)
from ._distributions import (
    ethanol_no_RIN_price_distribution,
    natural_gas_price_distribution,
    biomass_based_diesel_price_distribution,
    mean_biomass_based_diesel_price,
    cellulosic_based_diesel_price_distribution,
    mean_cellulosic_based_diesel_price,
    cellulosic_ethanol_price_distribution,
    advanced_ethanol_price_distribution,
    mean_ethanol_no_RIN_price,
    mean_advanced_ethanol_price,
    mean_cellulosic_ethanol_price,
    mean_glycerol_price,
    mean_natural_gas_price,
    mean_electricity_price,
)
from ._lca_characterization_factors import (
    GWP_characterization_factors,
    set_GWPCF,
    GWP,
)
from ._tables import (
    save_detailed_expenditure_tables, 
    save_detailed_life_cycle_tables
)
from biorefineries import cornstover as cs

energycane_dct = dict(
    Water=60,
    Sucrose=9.0,
    Glucose=0.1,
    Fiber=29.3,
    Ash=2.8,
)

_system_loaded = False
_chemicals_loaded = False

PRS = cs.PretreatmentReactorSystem
PRS_cost_item = PRS.cost_items['Pretreatment reactor system']
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
dry_biomass_yield = None # MT / hc
baseline_dry_biomass_yield = 39 # MT / hc 

cellulosic_configurations = frozenset([-2, 2, 4, 6, 8])
biodiesel_configurations = frozenset([1, 2, 5, 6, 7, 8])
ethanol_configurations = frozenset([-2, -1, 1, 2, 3, 4])
actag_configurations = frozenset([9, 10])
conventional_ethanol_configurations = ethanol_configurations.difference(cellulosic_configurations)
cellulosic_ethanol_configurations = cellulosic_configurations.intersection(ethanol_configurations)
ethanol_biodiesel_configurations = ethanol_configurations.intersection(biodiesel_configurations)
cellulosic_ethanol_biodiesel_configurations = ethanol_biodiesel_configurations.intersection(cellulosic_ethanol_configurations)
conventional_ethanol_biodiesel_configurations = ethanol_biodiesel_configurations.difference(cellulosic_ethanol_configurations)
cellulosic_oil_configurations = cellulosic_configurations.intersection(biodiesel_configurations).difference(ethanol_configurations)
oil_configurations = biodiesel_configurations.difference(cellulosic_configurations).difference(ethanol_configurations)

configuration_names = (
    'S1', 'O1', 'S2', 'O2', 'S1*', 'O1*', 'S2*', 'O2*', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', # 'O9', 'O10',
    'S2+', 'O2+', 'S2*+', 'O2*+', 'O4+', 'O6+', 'O8+', # 'O9'
)
comparison_names = ( 
    'O1 - S1', 'O2 - S2',  'O2 - O1', 'O1* - O1',  'O2* - O2', 'O5 - O7',
    'O6 - O8', 'O5 - O7', 'O6 - O5', 'O8 - O7',
)
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
    # 9: create_oilcane_to_biodiesel_and_actag_1g,
    # 10: create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation,
}
system_factory_options = {
    -1: dict(use_area_convention=True, pellet_bagasse=False),
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

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def disable_derivative(disable=True):
    global _derivative_disabled
    _derivative_disabled = disable
    for dct in cache.values(): dct['_derivative_disabled'] = disable
    
def enable_derivative(enable=True):
    global _derivative_disabled
    _derivative_disabled = disable = not enable
    for dct in cache.values(): dct['_derivative_disabled'] = disable

def rename_storage_units(units, storage_area):
    bst.rename_units([i for i in units if bst.is_storage_unit(i)], storage_area)
    
_derivative_disabled = False
cache = {}

def load(name, cache=cache, reduce_chemicals=False, RIN=True,
         avoid_natural_gas=True):
    dct = globals()
    number, agile, energycane = dct['configuration'] = configuration = parse_configuration(name)
    key = (number, agile, energycane, RIN)
    if cache is not None and key in cache:
        dct.update(cache[key])
        return
    global cane_sys, sys, tea, specs, flowsheet, _system_loaded
    global oil_extraction_specification, model, unit_groups
    global HXN, BT
    if not _chemicals_loaded: load_chemicals()
    flowsheet_name = format_configuration(configuration, latex=False)
    flowsheet = bst.Flowsheet(flowsheet_name)
    main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    s = flowsheet.stream
    u = flowsheet.unit
    get_stream = lambda ID: s[ID] if ID in flowsheet.stream else MockStream(ID) # Retrieve mock-stream in case stream does not exist in configuration
    
    ## System
    if number in system_factories:
        dct['cane_sys'] = cane_sys = system_factories[number](operating_hours=24 * 200, **system_factory_options.get(number, {}))
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
                dct['bagasse_splitter'] = splitter
                minimum_fraction_burned = 0
                maximum_fraction_burned = 0.7
                recycle_data = {}
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
                    cane_sys.simulate(material_data=material_data, update_material_data=True)
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
        prs = flowsheet(cs.units.PretreatmentReactorSystem)
        saccharification = flowsheet(cs.units.Saccharification)
        seed_train = flowsheet(cs.units.SeedTrain)
        fermentor = flowsheet(cs.units.CoFermentation)
        dct['pretreatment_rxnsys'] = tmo.ReactionSystem(
            prs.reactions, saccharification.saccharification
        )
        dct['fermentation_rxnsys'] = tmo.ReactionSystem(
            seed_train.reactions, fermentor.cofermentation
        )
        dct['cellulosic_rxnsys'] = tmo.ReactionSystem(
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
        dct['oilsorghum'] = oilsorghum = bst.Stream(
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
            dct['cane_mode'] = cane_mode = sys.operation_mode(cane_sys,
                operating_hours=180*24, oil_content=0.02, feedstock=feedstock.copy(),
                z_mass_carbs_baseline=0.091,
                z_mass_solids_baseline=0., 
                z_mass_ash_baseline=0.028,
                z_mass_water_baseline=0.60,
                glucose_yield=85, xylose_yield=65, 
                FFA_content=0.10, PL_content=0.10
            )
        else:
            dct['cane_mode'] = cane_mode = sys.operation_mode(cane_sys,
                operating_hours=180*24, oil_content=0.10, feedstock=feedstock.copy(),
                z_mass_carbs_baseline=0.1491, glucose_yield=85, xylose_yield=65, 
                FFA_content=0.10, PL_content=0.10
            )
        dct['sorghum_mode'] = sorghum_mode = sys.operation_mode(cane_sys, 
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
    dct['composition_specification'] = composition_specification = CaneCompositionSpecification(feedstock)
    if number < 0:
        isplit_a = None
        isplit_b = None
        oil_extraction_specification = MockExtractionSpecification()
    else:
        isplit_b = isplit_a = None
        for i in cane_sys.cost_units:
            if getattr(i, 'tag', None) == 'oil extraction':
                isplit_a = i.isplit
                break
        
        for i in cane_sys.cost_units:
            if getattr(i, 'tag', None) == 'bagasse oil extraction':
                isplit_b = i.isplit
                break
        
        oil_extraction_specification = OilExtractionSpecification(
            sys, isplit_a, isplit_b, 
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
            glucose_yield = fermentor.cofermentation.X[0]
            xylose_yield = fermentor.cofermentation.X[1]
            juice_sugar = s.juice.imass['Glucose', 'Sucrose'].sum() * glucose_yield
            hydrolysate_sugar = (
                s.slurry.imass['Glucose'] * glucose_yield
                + s.slurry.imass['Xylose', 'Arabinose'].sum() * xylose_yield
            )
            RIN_splitter.split[:] = juice_sugar / (juice_sugar + hydrolysate_sugar)
        
        cane_sys.update_configuration([*sys.units, RIN_splitter])
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
                glucose_yield = fermentor.cofermentation.X[0]
                xylose_yield = fermentor.cofermentation.X[1]
                juice_sugar = s.juice.imass['Glucose', 'Sucrose'].sum() * glucose_yield
                hydrolysate_sugar = (
                    s.slurry.imass['Glucose'] * glucose_yield
                    + s.slurry.imass['Xylose', 'Arabinose'].sum() * xylose_yield
                )
                RIN_splitter.split[:] = juice_sugar / (juice_sugar + hydrolysate_sugar)
            
            cane_sys.update_configuration([*sys.units, RIN_splitter])
            assert RIN_splitter in cane_sys.units
        else:
            # A biodiesel stream should already exist
            s.biodiesel.register_alias('biomass_based_diesel')
    
    if RIN: 
        cbpd = cellulosic_based_diesel_price_distribution 
        mcbp = mean_cellulosic_based_diesel_price
        bpd = biomass_based_diesel_price_distribution 
        mbp = mean_biomass_based_diesel_price
        cepd = cellulosic_ethanol_price_distribution 
        mcep = mean_cellulosic_ethanol_price
        aepd = advanced_ethanol_price_distribution 
        maep = mean_advanced_ethanol_price
    else:
        raise NotImplementedError('biodiesel price distribution without RINs has not been implemented in `_distributions.py` module')
        cepd = aepd = ethanol_no_RIN_price_distribution 
        mcep = maep = mean_ethanol_no_RIN_price
    
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
    dct['natural_gas_streams'] = natural_gas_streams = [natural_gas]
    for stream in natural_gas_streams: set_GWPCF(stream, 'CH4')
    
    ## Add BioSTEAM objects to module for easy access
    dct.update(flowsheet.to_dict())
    
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
    
    @uniform(60, 90, units='%', kind='coupled')
    def set_crushing_mill_oil_recovery(oil_recovery):
        oil_extraction_specification.load_crushing_mill_oil_recovery(oil_recovery / 100.)
    
    @uniform(70.0, 90, units='%', kind='coupled')
    def set_saccharification_oil_recovery(saccharification_oil_recovery):
        oil_extraction_specification.load_saccharification_oil_recovery(saccharification_oil_recovery / 100.)

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

    # USDA ERS historical price data with EPA RIN prices
    @parameter(distribution=cepd, element=cellulosic_ethanol, 
               baseline=mcep, units='USD/L')
    def set_cellulosic_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
        cellulosic_ethanol.price = price * ethanol_L_per_kg
        
    @parameter(distribution=aepd, element=advanced_ethanol, 
               baseline=maep, units='USD/L')
    def set_advanced_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
        advanced_ethanol.price = price * ethanol_L_per_kg
        
    # USDA ERS historical price data
    @parameter(distribution=bpd, element=biomass_based_diesel, units='USD/L', baseline=mbp)
    def set_biomass_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        biomass_based_diesel.price = price * biodiesel_L_per_kg

    @parameter(distribution=cbpd, element=cellulosic_based_diesel, units='USD/L', baseline=mcbp)
    def set_cellulosic_based_diesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        cellulosic_based_diesel.price = price * biodiesel_L_per_kg

    # https://www.eia.gov/energyexplained/natural-gas/prices.php
    @parameter(distribution=natural_gas_price_distribution, element=s.natural_gas, units='USD/m3',
               baseline=4.73 * 35.3146667/1e3)
    def set_natural_gas_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        BT.natural_gas_price = price * V_ng

    # https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
    # Data from historical prices, 2010-2020
    @triangular(0.0583, 0.065, 0.069, units='USD/kWhr',
                baseline=0.0637)
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
    
    @uniform(85, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=85, kind='coupled')
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
        
    @uniform(65, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=65, kind='coupled')
    def set_cane_xylose_yield(cane_xylose_yield):
        if agile:
            cane_mode.xylose_yield = cane_xylose_yield
        elif abs(number) in cellulosic_configurations:
            set_xylose_yield(cane_xylose_yield)
    
    @uniform(86, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=86, kind='coupled')
    def set_sorghum_xylose_yield(sorghum_xylose_yield):
        if not agile: return
        sorghum_mode.xylose_yield = sorghum_xylose_yield
    
    @uniform(90, 95, units='%', element='Cofermenation',
             baseline=90, kind='coupled')
    def set_glucose_to_ethanol_yield(glucose_to_ethanol_yield):
        if number in cellulosic_ethanol_configurations:
            glucose_to_ethanol_yield *= 0.01
            # fermentor.cofermentation[2].X = 0.004 # Baseline
            # fermentor.cofermentation[3].X = 0.006 # Baseline
            # fermentor.loss[0].X = 0.03 # Baseline
            split = np.mean(u.S401.split)
            X1 = split * seed_train.reactions.X[0]
            X2 = split * seed_train.reactions.X[2]
            X3 = (glucose_to_ethanol_yield - X1) / (1 - X1 - X2)
            split = np.mean(u.S401.split)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: breakpoint()
            fermentor.cofermentation.X[0] = X3
            fermentor.cofermentation.X[2] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass
    
    @uniform(50, 95, units='%', element='Cofermenation',
             baseline=50, kind='coupled')
    def set_xylose_to_ethanol_yield(xylose_to_ethanol_yield):
        if number in cellulosic_ethanol_configurations:
            xylose_to_ethanol_yield *= 0.01
            split = np.mean(u.S401.split)
            X1 = split * seed_train.reactions.X[1]
            X2 = split * seed_train.reactions.X[3]
            X3 = (xylose_to_ethanol_yield - X1) / (1 - X1 - X2)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: breakpoint()
            fermentor.cofermentation.X[1] = X3
            fermentor.cofermentation.X[3] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass

    @uniform(68.5, 137, units='g/L', element='Cofermentation',
             baseline=68.5, kind='coupled')
    def set_cofermentation_titer(titer):
        if number in cellulosic_ethanol_configurations: fermentor.titer = titer

    @uniform(0.951, 1.902, units='g/L', element='Cofermentation',
             baseline=0.951)
    def set_cofermentation_productivity(productivity):
        if number in cellulosic_ethanol_configurations: fermentor.productivity = productivity

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

    @uniform(5., 15., element='oilcane', units='dry wt. %', kind='coupled')
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
        if number > 0:
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
        biodiesel_flow = lambda: sys.flow_rates.get(biodiesel, 0.) * biodiesel_L_per_kg # L / yr
        ethanol_flow = lambda: (sys.flow_rates.get(cellulosic_ethanol, 0.) + sys.flow_rates.get(advanced_ethanol, 0.))* ethanol_L_per_kg # L / yr
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
            return f * mean_advanced_ethanol_price + (1 - f) * mean_cellulosic_ethanol_price
        else:
            return mean_advanced_ethanol_price
    
    def get_GWP_mean_biodiesel_price():
        if number in cellulosic_configurations:
            f = float(RIN_splitter.split.mean())
            return f * mean_biomass_based_diesel_price + (1 - f) * mean_cellulosic_based_diesel_price
        else:
            return mean_advanced_ethanol_price
    
    sys.define_process_impact(
        key=GWP,
        name='Direct non-biogenic emissions',
        basis='kg',
        inventory=direct_nonbiogenic_emissions,
        CF=1.,
    )
    
    dct['flows'] = {
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
    
    @metric(units='L/hc')
    def biodiesel_yield():
        if dry_biomass_yield is None: 
            return None
        else:
            return biodiesel_flow() / feedstock_consumption.get() * dry_biomass_yield
    
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
            + crude_glycerol_flow() * mean_glycerol_price
            + max(-electricity(), 0) * mean_electricity_price
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
            return GWP_economic.get() * mean_glycerol_price
        else:
            return 0.
    
    @metric(name='Electricity GWP', element='Economic allocation', units='kg*CO2*eq / MWhr')
    def GWP_electricity(): # Cradle to gate
        if electricity_production.get():
            return GWP_economic.get() * mean_electricity_price * 1000.

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
        if dry_biomass_yield is None: return None
        # Set prelimiary feedstock price assuming 10% is due to transportation
        # and 90% is based on productivity (a function of height)
        feedstock.price = (0.9 * baseline_dry_biomass_yield / dry_biomass_yield + 0.10) * 0.035
        return 100. * tea.solve_IRR()

    @metric(units='USD/MT')
    def MFPP_derivative():
        if number < 0:
            return 0.
        if _derivative_disabled:
            return np.nan
        if agile:
            cane_mode.oil_content += 0.01
            sorghum_mode.oil_content += 0.01
        else:
            composition_specification.load_oil_content(composition_specification.oil + 0.01)
        sys.simulate()
        # value = (kg_per_MT * tea.solve_price(feedstock) - MFPP.cache)
        # feedstock.price = tea.solve_price(feedstock)
        # print('AFTER')
        # print('MFPP', kg_per_MT * tea.solve_price(feedstock))
        # print('VOC', tea.VOC / 1e3)
        # print('TCI', tea.TCI / 1e6)
        # print('sales', tea.sales / 1e3)
        # print('NPV', tea.NPV)
        return MFPP.difference()
    
    @metric(units='L/MT')
    def biodiesel_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return biodiesel_production.difference()
    
    @metric(units='L/MT')
    def ethanol_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return ethanol_production.difference()
    
    @metric(units='kWhr/MT')
    def electricity_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return electricity_production.difference()
    
    @metric(units='cf/MT')
    def natural_gas_consumption_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        # print('natural gas production derivative', value)
        return natural_gas_consumption.difference()
    
    @metric(units='10^6*USD')
    def TCI_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return TCI.difference()
    
    @metric(name='GWP derivative', element='Economic allocation', units='kg*CO2*eq / USD')
    def GWP_economic_derivative(): # Cradle to gate
        if number < 0: return 0.
        if _derivative_disabled: return 0.
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
            return GWP_economic_derivative.get() * mean_glycerol_price
        else:
            return 0.
    
    @metric(name='Electricity GWP derivative', element='Electricity', units='kg*CO2*eq / MWhr')
    def GWP_electricity_derivative(): # Cradle to gate
        if electricity_production.get():
            return GWP_economic_derivative.get() * mean_electricity_price * 1000.
        else:
            return 0.
    
    # @metric
    # def MASP():
    #     return tea.solve_price(s.acTAG) if number in actag_configurations else 0.
    
    # @metric(units='MMGGE/yr')
    # def productivity():
    #     GGE = (ethanol_flow() / 1.5
    #         + biodiesel_flow() / 0.9536
    #         - electricity() * 3600 / 131760
    #         - natural_gas_flow() / 126.67)
    #     return GGE / 1e6
    
    # Single point evaluation for detailed design results
    def set_baseline(p, x):
        p.setter(x)
        p.baseline = x
    
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
    set_baseline(set_saccharification_oil_recovery, 70)
    set_baseline(set_crushing_mill_oil_recovery, 60)
    set_baseline(set_advanced_ethanol_price, maep) 
    set_baseline(set_cellulosic_ethanol_price, mcep) 
    set_baseline(set_biomass_based_diesel_price, mbp)
    set_baseline(set_cellulosic_based_diesel_price, mcbp)
    set_baseline(set_crude_glycerol_price, mean_glycerol_price)
    set_baseline(set_natural_gas_price, mean_natural_gas_price)
    set_baseline(set_electricity_price, mean_electricity_price)
    if number in cellulosic_ethanol_configurations:
        get_stream('ethanol').price = 0.789
    if number > 0:
        set_baseline(set_cane_PL_content, 10)
        set_baseline(set_cane_FFA_content, 10)
    # set_fermentation_solids_loading(20) # Same as Humbird
    # set_feedstock_oil_content(10) # Consistent with SI of Huang's 2016 paper
    # set_ethanol_price(2.356) # Consistent with Huang's 2016 paper
    # set_biodiesel_price(4.569) # Consistent with Huang's 2016 paper
    # set_natural_gas_price(4.198) # Consistent with Humbird's 2012 paper
    # set_electricity_price(0.0572) # Consistent with Humbird's 2012 paper
    # set_operating_days(200) # Consistent with Huang's 2016 paper
    
    for i in sys.units:
        if isinstance(i, bst.MultiEffectEvaporator): i.flash = False
    
    for i in model._parameters:
        dct[i.setter.__name__] = i
    for i in model._metrics:
        dct[i.getter.__name__] = i
    if cache is not None: cache[key] = dct.copy()
    
    ## Simulation
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
        # sys.track_recycle(WWTsys.recycle)
    sys.simulate()
    feedstock.price = tea.solve_price(feedstock)
    if reduce_chemicals: 
        cane_sys.reduce_chemicals()
        cane_sys._load_stream_links()
        HXN.simulate()

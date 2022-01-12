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
    _system,
    _tea,
    _oil_extraction_specification,
    _distributions,
    _evaluation,
    _box_plots,
    _contour_plots,
    _lca_characterization_factors,
    _load_data,
    _parse_configuration,
)
from .units import *
from ._process_settings import *
from ._chemicals import *
from ._contour_plots import *
from ._system import *
from ._tea import *
from ._oil_extraction_specification import *
from ._distributions import *
from ._evaluation import *
from ._box_plots import *
from ._contour_plots import *
from ._lca_characterization_factors import *
from ._load_data import *
from ._parse_configuration import *
from ._variable_mockups import *

__all__ = (
    units.__all__,
    *_process_settings.__all__,
    *_chemicals.__all__,
    *_system.__all__,
    *_tea.__all__,
    *_oil_extraction_specification.__all__,
    *_distributions.__all__,
    *_evaluation.__all__,
    *_box_plots.__all__,
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
from biorefineries.sugarcane import create_sugarcane_to_ethanol_system
from biorefineries.lipidcane import (
    set_lipid_fraction as set_oil_fraction, 
)
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
import numpy as np
from ._process_settings import load_process_settings
from ._chemicals import create_chemicals
from ._system import (
    create_oilcane_to_biodiesel_and_ethanol_1g,
    create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation,
    create_sugarcane_to_ethanol_combined_1_and_2g,    
)
from ._parse_configuration import (
    parse,
    format_configuration,
)
from ._tea import (
    create_tea,
)
from ._oil_extraction_specification import (
    OilExtractionSpecification,
    MockExtractionSpecification,
)
from ._distributions import (
    ethanol_price_distribution,
    biodiesel_minus_ethanol_price_distribution,
    natural_gas_price_distribution,
    mean_biodiesel_price,
    mean_ethanol_price,
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

_system_loaded = False
_chemicals_loaded = False

PRS = cs.PretreatmentReactorSystem
PRS_cost_item = PRS.cost_items['Pretreatment reactor system']
kg_per_ton = 907.18474


configuration_names = (
    'S1', 'O1', 'S2', 'O2', 'S1*', 'O1*', 'S2*', 'O2*',
)
comparison_names = (
    # 'I - ∅', 
    'O1 - S1', 
    'O2 - S2', 
    'O2 - O1', 
    'O1* - O1', 
    'O2* - O2',  
)

other_comparison_names = (
    'O1* - S1*', 'O2* - S2*', 
)

across_oil_content_names = (
    'O1', 'O2', 
)

across_oil_content_agile_names = (
    'O1*', 'O2*', 
)

across_oil_content_comparison_names = (
    'O1 - S1', 'O2 - S2', 'O2 - O1', 
)

across_oil_content_agile_direct_comparison_names = (
    'O1* - O1', 'O2* - O2', 
)

across_oil_content_agile_comparison_names = (
    'O1* - S1*', 'O2* - S2*', 'O2* - O1*', 
)

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def disable_derivative(disable=True):
    global _derivative_disabled
    _derivative_disabled = disable
    
def enable_derivative(enable=True):
    global _derivative_disabled
    _derivative_disabled = not enable
    
_derivative_disabled = False

def load(name, cache={}, reduce_chemicals=True, 
         enhanced_cellulosic_performance=False,
         enhanced_biodiesel_production=False):
    dct = globals()
    number, agile = dct['configuration'] = configuration = parse(name)
    key = (number, agile, enhanced_cellulosic_performance)
    if key in cache:
        dct.update(cache[key])
        return
    global oilcane_sys, sys, tea, specs, flowsheet, _system_loaded
    global oil_extraction_specification, model, unit_groups
    global HXN, BT
    if not _chemicals_loaded: load_chemicals()
    flowsheet_name = format_configuration(configuration)
    if enhanced_cellulosic_performance:
        flowsheet_name += '_enhanced_fermentation'
    flowsheet = bst.Flowsheet(flowsheet_name)
    main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    s = flowsheet.stream
    u = flowsheet.unit
    operating_hours = 24 * 200
    
    ## System
    
    area_names = None
    def rename_storage_units(storage):
        bst.rename_units([i for i in oilcane_sys.units if bst.is_storage_unit(i)], storage)
    
    if number == -1:
        isplit_efficiency_is_reversed = None
        # starting_chemicals = create_starting_chemicals()
        # bst.settings.set_thermo(starting_chemicals)
        oilcane_sys = create_sugarcane_to_ethanol_system(
            operating_hours=operating_hours,
            use_area_convention=True,
            pellet_bagasse=True,
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'EtOH prod.', 
            'CH&P',
            'Utilities',
            'HXN',
            'Storage',
        ]
        rename_storage_units(700)
    elif number == -2:
        isplit_efficiency_is_reversed = None
        oilcane_sys = create_sugarcane_to_ethanol_combined_1_and_2g(
            operating_hours=operating_hours,
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Pretreatment',
            'EtOH prod.',
            'Wastewater treatment',
            'CH&P', 
            'Utilities',
            'HXN',
            'Storage',
        ]
        rename_storage_units(900)
    elif number == 1:
        isplit_efficiency_is_reversed = False
        oilcane_sys = create_oilcane_to_biodiesel_and_ethanol_1g(
            operating_hours=operating_hours,
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'EtOH prod.', 
            'Oil ext.',
            'Biod. prod.', 
            'CH&P',
            'Utilities',
            'HXN',
            'Storage',
        ]
        rename_storage_units(1000)
    elif number == 2:
        isplit_efficiency_is_reversed = True
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Pretreatment',
            'EtOH prod.',
            'Wastewater treatment',
            'Oil ext.',
            'CH&P', 
            'Biod. prod.',
            'Utilities',
            'HXN',
            'Storage',
        ]
        oilcane_sys = create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1100)
    else:
        raise NotImplementedError(number)
    oilcane_sys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    dct.update(flowsheet.to_dict())
    
    def get_stream(ID):
        if ID in flowsheet.stream:
            return s[ID]
        else:
            return MockStream(ID)
    
    crude_glycerol = get_stream('crude_glycerol')
    pure_glycerine = get_stream('pure_glycerine')
    if number == 1:
        oilcane_sys.prioritize_unit(u.T608)
    elif number == 2:
        oilcane_sys.prioritize_unit(u.T808)
    if number < 0:
        dct['feedstock'] = feedstock = s.sugarcane
        dct['oilcane_sys'] = oilcane_sys
    else:
        dct['feedstock'] = feedstock = s.oilcane
    unit_groups = UnitGroup.group_by_area(oilcane_sys.units)
    if area_names:
        for i, j in zip(unit_groups, area_names): i.name = j
    for i in unit_groups: i.autofill_metrics(shorthand=True)
    
    for BT in oilcane_sys.units:
        if isinstance(BT, bst.BoilerTurbogenerator): break

    HXN = None
    for HXN_group in unit_groups:
        if HXN_group.name == 'HXN':
            HXN_group.filter_savings = False
            HXN = HXN_group.units[0]
            assert isinstance(HXN, bst.HeatExchangerNetwork)
    unit_groups[-1].metrics[-1].getter = lambda: 0.    
    
    
    if abs(number) == 2:
        prs, = flowsheet(cs.units.PretreatmentReactorSystem)
        saccharification, = flowsheet(cs.units.Saccharification)
        seed_train, = flowsheet(cs.units.SeedTrain)
        fermentor, = flowsheet(cs.units.CoFermentation)
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
        if abs(number) == 2:
            glucose_yield *= 0.01
            X1 = prs.reactions.X[0]
            X1_side = prs.reactions.X[1:3].sum()
            X2_side = saccharification.saccharification.X[:2].sum()
            saccharification.saccharification.X[2] = X2 = (glucose_yield - X1) / (1 - X1_side)
            X_excess = (X2_side + X2) - 1
            if X_excess > 0: breakpoint()
            
    def set_xylose_yield(xylose_yield):
        if abs(number) == 2:
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
            if number > 0: 
                set_oil_fraction(oil_content, feedstock,
                                 FFA_fraction=mode.FFA_content,
                                 z_mass_carbs_baseline=mode.z_mass_carbs_baseline,
                                 PL_fraction=mode.PL_content)
            else:
                F_mass = feedstock.F_mass
                feedstock.copy_flow(mode.feedstock)
                feedstock.F_mass = F_mass
        sys.operation_parameter(set_glucose_yield)
        sys.operation_parameter(set_xylose_yield)
        
        dct['cane_mode'] = cane_mode = sys.operation_mode(oilcane_sys,
            operating_hours=200*24, oil_content=0.05, feedstock=feedstock.copy(),
            z_mass_carbs_baseline=0.1491, glucose_yield=85, xylose_yield=65, 
            FFA_content=0.10, PL_content=0.10
        )
        dct['sorghum_mode'] = sorghum_mode = sys.operation_mode(oilcane_sys, 
            operating_hours=60*24, oil_content=0.05, glucose_yield=79, xylose_yield=86,
            feedstock=oilsorghum,
            z_mass_carbs_baseline=0.1371, FFA_content=0.10, PL_content=0.10,
        )
        tea = create_tea(sys)
        tea.operating_days = 260 
        tea.IRR = 0.10
    else:
        sys = oilcane_sys
        tea = create_tea(sys)
        tea.operating_days = 200
        tea.IRR = 0.10
        
        # def electricity_GWP():
        #     power_utility = bst.PowerUtility.sum([i.power_utility for i in sys.cost_units])
        #     return power_utility.get_impact(production_key=GWP) * sys.operating_hours
        
    tea.income_tax = 0.21 # Davis et al. 2018; https://www.nrel.gov/docs/fy19osti/71949.pdf
    
    ## Specification for analysis
    if number < 0:
        dct['biodiesel'] = bst.Stream('biodiesel')
        isplit_a = None
        isplit_b = None
        oil_extraction_specification = MockExtractionSpecification()
    else:
        for i in oilcane_sys.cost_units:
            if getattr(i, 'tag', None) == 'oil extraction efficiency':
                isplit_a = i.isplit
                break
        
        for i in oilcane_sys.cost_units:
            if getattr(i, 'tag', None) == 'bagasse oil retention':
                isplit_b = i.isplit
                break
        
        oil_extraction_specification = OilExtractionSpecification(
            sys, [feedstock], isplit_a, isplit_b, isplit_efficiency_is_reversed
        )
    
    ## LCA
    
    # Set non-negligible characterization factors
    if abs(number) != 2:
        for i in ('FGD_lime', 'cellulase', 'DAP', 'CSL', 'caustic'): MockStream(i)
    if number < 0:
        for i in ('catalyst', 'methanol', 'HCl', 'NaOH', 'crude_glycerol', 'pure_glycerine'): MockStream(i)
    if abs(number) != 1:
        MockStream('dryer_natural_gas')
        
    set_GWPCF(feedstock, 'sugarcane')
    set_GWPCF(s.H3PO4, 'H3PO4')
    set_GWPCF(s.lime, 'lime', dilution=0.046) # Diluted with water
    set_GWPCF(s.denaturant, 'gasoline')
    set_GWPCF(s.FGD_lime, 'lime', dilution=0.451)
    set_GWPCF(s.cellulase, 'cellulase', dilution=0.02) 
    set_GWPCF(s.DAP, 'DAP')
    set_GWPCF(s.CSL, 'CSL')
    set_GWPCF(s.caustic, 'NaOH', 0.5)
    set_GWPCF(s.catalyst, 'NaOH', 0.5)
    set_GWPCF(s.catalyst, 'methanol catalyst mixture')
    set_GWPCF(s.methanol, 'methanol')
    set_GWPCF(s.HCl, 'HCl')
    set_GWPCF(s.NaOH, 'NaOH')
    set_GWPCF(s.pure_glycerine, 'pure-glycerol')
    set_GWPCF(s.dryer_natural_gas, 'CH4')
    set_GWPCF(s.crude_glycerol, 'crude-glycerol', dilution=0.80)
    set_GWPCF(s.biodiesel, 'biodiesel displacement')
    bst.PowerUtility.set_CF(GWP, GWP_characterization_factors['Electricity'])
    dct['natural_gas_streams'] = natural_gas_streams = [s.natural_gas]
    if abs(number) == 1: natural_gas_streams.append(s.dryer_natural_gas)
    for stream in natural_gas_streams:
        set_GWPCF(stream, 'CH4')
    
    ## Model
    model = bst.Model(sys, exception_hook='raise', retry_evaluation=False)
    parameter = model.parameter
    metric = model.metric
    
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
    
    @uniform(40, 70, units='%', kind='coupled')
    def set_bagasse_oil_retention(oil_retention):
        oil_extraction_specification.load_oil_retention(oil_retention / 100.)
    
    def oil_extraction_efficiency_hook(x):
        if number < 0:
            return x
        elif number == 1:
            return 50.0 + x
        elif number == 2:
            return 70.0 + x
    
    @uniform(0., 20, units='%', kind='coupled', 
             hook=oil_extraction_efficiency_hook)
    def set_bagasse_oil_extraction_efficiency(bagasse_oil_extraction_efficiency):
        oil_extraction_specification.load_efficiency(bagasse_oil_extraction_efficiency / 100.)

    capacity = feedstock.F_mass / kg_per_ton
    @default(capacity, units='ton/hr', kind='coupled')
    def set_plant_capacity(capacity):
        feedstock.F_mass = F_mass = kg_per_ton * capacity
        if agile: oilsorghum.F_mass = F_mass

    # USDA ERS historical price data
    @parameter(distribution=ethanol_price_distribution, element=s.ethanol, 
               baseline=1.90, units='USD/gal')
    def set_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
        s.ethanol.price = price / 2.98668849
        
    # USDA ERS historical price data
    @parameter(distribution=biodiesel_minus_ethanol_price_distribution, element=s.biodiesel, units='USD/gal',
               baseline=2.47, hook=lambda x: s.ethanol.price * 2.98668849 + x)
    def set_biodiesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        s.biodiesel.price = price / 3.3111

    # https://www.eia.gov/energyexplained/natural-gas/prices.php
    @parameter(distribution=natural_gas_price_distribution, element=s.natural_gas, units='USD/cf',
               baseline=4.73)
    def set_natural_gas_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        BT.natural_gas_price = 51.92624700383502 * price / 1000. 

    # https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
    # Data from historical prices, 2010-2020
    @triangular(0.0583, 0.065, 0.069, units='USD/kWhr',
                baseline=0.0637)
    def set_electricity_price(electricity_price): 
        bst.PowerUtility.price = electricity_price
        
    # From Huang's 2016 paper
    @uniform(6 * 30, 7 * 30, units='day/yr', baseline=200)
    def set_operating_days(operating_days):
        if agile:
            cane_mode.operating_hours = operating_days * 24
        else:
            tea.operating_days = operating_days
    
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
        if abs(number) == 2: saccharification.tau = reaction_time
    
    @default(0.212, units='USD/kg', element='cellulase')
    def set_cellulase_price(price):
        if abs(number) == 2: s.cellulase.price = price
    
    @default(0.02, units='wt. % cellulose', element='cellulase')
    def set_cellulase_loading(cellulase_loading):
        if abs(number) == 2: u.M302.cellulase_loading = cellulase_loading
    
    @default(PRS_cost_item.cost, units='million USD', element='Pretreatment reactor system')
    def set_reactor_base_cost(base_cost):
        PRS_cost_item.cost = base_cost
    
    @uniform(85, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=85)
    def set_cane_glucose_yield(cane_glucose_yield):
        if agile:
            cane_mode.glucose_yield = cane_glucose_yield
        elif abs(number) == 2:
            set_glucose_yield(cane_glucose_yield)
    
    @uniform(79, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=79)
    def set_sorghum_glucose_yield(sorghum_glucose_yield):
        if not agile: return
        sorghum_mode.glucose_yield = sorghum_glucose_yield
        
    @uniform(65, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=65)
    def set_cane_xylose_yield(cane_xylose_yield):
        if agile:
            cane_mode.xylose_yield = cane_xylose_yield
        elif abs(number) == 2:
            set_xylose_yield(cane_xylose_yield)
    
    @uniform(86, 97.5, units='%', element='Pretreatment and saccharification',
             baseline=86)
    def set_sorghum_xylose_yield(sorghum_xylose_yield):
        if not agile: return
        sorghum_mode.xylose_yield = sorghum_xylose_yield
    
    @uniform(90, 95, units='%', element='Cofermenation',
             baseline=90)
    def set_glucose_to_ethanol_yield(glucose_to_ethanol_yield):
        if abs(number) == 2:
            glucose_to_ethanol_yield *= 0.01
            # fermentor.cofermentation[2].X = 0.004 # Baseline
            # fermentor.cofermentation[3].X = 0.006 # Baseline
            # fermentor.loss[0].X = 0.03 # Baseline
            split = np.mean(u.S403.split)
            X1 = split * seed_train.reactions.X[0]
            X2 = split * seed_train.reactions.X[2]
            X3 = (glucose_to_ethanol_yield - X1) / (1 - X1 - X2)
            split = np.mean(u.S403.split)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: breakpoint()
            fermentor.cofermentation.X[0] = X3
            fermentor.cofermentation.X[2] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass
    
    @uniform(50, 95, units='%', element='Cofermenation',
             baseline=50)
    def set_xylose_to_ethanol_yield(xylose_to_ethanol_yield):
        if abs(number) == 2:
            # fermentor.cofermentation[6].X = 0.004 # Baseline
            # fermentor.cofermentation[7].X = 0.046 # Baseline
            # fermentor.cofermentation[8].X = 0.009 # Baseline
            # fermentor.loss[1].X = 0.03 # Baseline
            xylose_to_ethanol_yield *= 0.01
            split = np.mean(u.S403.split)
            X1 = split * seed_train.reactions.X[1]
            X2 = split * seed_train.reactions.X[3]
            X3 = (xylose_to_ethanol_yield - X1) / (1 - X1 - X2)
            X_excess = X3 * 1.0526 - 1
            if X_excess > 0.: breakpoint()
            fermentor.cofermentation.X[1] = X3
            fermentor.cofermentation.X[3] = X3 * 0.0526 # 95% towards ethanol, the other 5% goes towards cell mass

    @uniform(68.5, 137, units='g/L', element='Cofermentation',
             baseline=68.5)
    def set_cofermentation_titer(titer):
        if abs(number) == 2: fermentor.titer = titer

    @uniform(0.951, 1.902, units='g/L', element='Cofermentation',
             baseline=0.951)
    def set_cofermentation_productivity(productivity):
        if abs(number) == 2: fermentor.productivity = productivity

    @default(10, element='oilcane', units='% oil', kind='coupled')
    def set_cane_PL_content(cane_PL_content):
        if agile: cane_mode.PL_content = cane_PL_content / 100.
        else: oil_extraction_specification.PL_content = cane_PL_content / 100.
    
    @default(10, element='oilsorghum', units='% oil', kind='coupled')
    def set_sorghum_PL_content(sorghum_PL_content):
        if agile: sorghum_mode.PL_content = sorghum_PL_content / 100.
    
    @default(10, element='oilcane', units='% oil', kind='coupled')
    def set_cane_FFA_content(cane_FFA_content):
        if agile: cane_mode.FFA_content = cane_FFA_content / 100.
        else: oil_extraction_specification.FFA_content = cane_FFA_content / 100.
    
    @default(10, element='oilsorghum', units='% oil', kind='coupled')
    def set_sorghum_FFA_content(sorghum_FFA_content):
        if agile: sorghum_mode.FFA_content = sorghum_FFA_content / 100. 

    @uniform(5., 15., element='oilcane', units='dry wt. %', kind='coupled')
    def set_cane_oil_content(cane_oil_content):
        if agile:
            cane_mode.oil_content = cane_oil_content / 100.
        else:
            oil_extraction_specification.load_oil_content(cane_oil_content / 100.)

    @uniform(-3., 0., element='oilsorghum', units='dry wt. %', kind='coupled',
             baseline=0.)
    def set_relative_sorghum_oil_content(relative_sorghum_oil_content):
        if agile:
            sorghum_mode.oil_content = cane_mode.oil_content + relative_sorghum_oil_content / 100.

    @default(23, units='% oil', kind='coupled', name='TAG to FFA conversion')
    def set_TAG_to_FFA_conversion(TAG_to_FFA_conversion):
        if number == 1:
            u.R301.oil_reaction.X[0] = TAG_to_FFA_conversion / 100.
        elif number == 2:
            u.R401.oil_reaction.X[0] = TAG_to_FFA_conversion / 100.
    
    @default_gwp(feedstock.characterization_factors[GWP], name='GWP', 
             element=feedstock, units='kg*CO2-eq/kg')
    def set_feedstock_GWP(value):
        if number > 0:
            feedstock.characterization_factors[GWP] = value
    
    @default_gwp(s.methanol.characterization_factors[GWP], name='GWP', 
                 element=s.methanol, units='kg*CO2-eq/kg')
    def set_methanol_GWP(value):
        s.methanol.characterization_factors[GWP] = value
    
    # @default(crude_glycerol.characterization_factors[GWP], name='GWP', 
    #          element=crude_glycerol, units='kg*CO2-eq/kg')
    # def set_crude_glycerol_GWP(value):
    #     crude_glycerol.characterization_factors[GWP] = value
    
    @default_gwp(s.pure_glycerine.characterization_factors[GWP], name='GWP', 
                 element=pure_glycerine, units='kg*CO2-eq/kg')
    def set_pure_glycerine_GWP(value):
        s.pure_glycerine.characterization_factors[GWP] = value
    
    @default_gwp(s.cellulase.characterization_factors[GWP], name='GWP', 
                 element=s.cellulase, units='kg*CO2-eq/kg')
    def set_cellulase_GWP(value):
        s.cellulase.characterization_factors[GWP] = value * 0.02
    
    @default_gwp(s.natural_gas.characterization_factors[GWP], name='GWP', 
                 element=s.natural_gas, units='kg*CO2-eq/kg')
    def set_natural_gas_GWP(value):
        for ng in natural_gas_streams:
            ng.characterization_factors[GWP] = value
    
    s.natural_gas.phase = 'g'
    s.natural_gas.set_property('T', 60, 'degF')
    s.natural_gas.set_property('P', 14.73, 'psi')
    original_value = s.natural_gas.imol['CH4']
    s.natural_gas.imass['CH4'] = 1 
    V_ng = s.natural_gas.get_total_flow('ft3/hr')
    s.natural_gas.imol['CH4'] = original_value
    
    if agile:
        feedstock_flow = lambda: sys.flow_rates[feedstock] / kg_per_ton # ton/yr
        biodiesel_flow = lambda: sys.flow_rates.get(s.biodiesel, 0.) / 3.3111 # gal/yr
        ethanol_flow = lambda: sys.flow_rates[s.ethanol] / 2.98668849 # gal/yr
        natural_gas_flow = lambda: sum([sys.flow_rates[i] for i in natural_gas_streams]) * V_ng # cf/yr
        crude_glycerol_flow = lambda: sys.flow_rates.get(s.crude_glycerol, 0.) # kg/yr
        
        @sys.operation_metric(annualize=True)
        def direct_nonbiogenic_emissions(mode):
            return sum([i.F_mol for i in natural_gas_streams]) * chemicals.CO2.MW
        
    else:
        feedstock_flow = lambda: sys.operating_hours * feedstock.F_mass / kg_per_ton # ton/yr
        biodiesel_flow = lambda: sys.operating_hours * s.biodiesel.F_mass / 3.3111 # gal/yr
        ethanol_flow = lambda: sys.operating_hours * s.ethanol.F_mass / 2.98668849 # gal/yr
        crude_glycerol_flow = lambda: sys.operating_hours * s.crude_glycerol.F_mass # kg/yr
        natural_gas_flow = lambda: sum([i.F_mass for i in natural_gas_streams]) * sys.operating_hours * V_ng # cf/yr
        direct_nonbiogenic_emissions = lambda: sum([i.F_mol for i in natural_gas_streams]) * chemicals.CO2.MW * sys.operating_hours
    electricity = lambda: sys.operating_hours * sys.power_utility.rate
    
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
    @metric(units='USD/ton')
    def MFPP():
        price = tea.solve_price(feedstock)
        return kg_per_ton * price
    
    @metric(units='ton/yr')
    def feedstock_consumption():
        return feedstock_flow()
    
    @metric(units='Gal/ton')
    def biodiesel_production():
        return biodiesel_flow() / feedstock_consumption.get()
    
    @metric(units='Gal/ton')
    def ethanol_production():
        return ethanol_flow() / feedstock_consumption.get()
    
    @metric(units='kWhr/ton')
    def electricity_production():
        value = - electricity() / feedstock_consumption.get()
        if value < 0.: value = 0.
        return value
    
    @metric(units='cf/ton')
    def natural_gas_consumption():
        value = natural_gas_flow() / feedstock_consumption.get()
        return value
    
    @metric(units='10^6*USD')
    def TCI():
        return tea.TCI / 1e6 # 10^6*$
    
    @metric(units='%')
    def heat_exchanger_network_error():
        return HXN.energy_balance_percent_error if HXN else 0.    

    def GWP_displacement(): # Cradle to gate
        GWP_total = sys.get_total_feeds_impact(GWP) + min(electricity(), 0) * GWP_characterization_factors['Electricity'] # kg CO2 eq. / yr
        sales = (
            biodiesel_flow() * mean_biodiesel_price
            + ethanol_flow() * mean_ethanol_price
            + crude_glycerol_flow() * mean_glycerol_price
        )
        GWP_per_USD = GWP_total / sales
        return {
            'Ethanol': GWP_per_USD * mean_ethanol_price,
            'Biodiesel': GWP_per_USD * mean_biodiesel_price,
            'Crude glycerol': GWP_per_USD * mean_glycerol_price,
        }
    dct['GWP_displacement'] = GWP_displacement

    @metric(name='GWP', element='Economic allocation', units='kg*CO2*eq / USD')
    def GWP_economic(): # Cradle to gate
        GWP_material = sys.get_total_feeds_impact(GWP) # kg CO2 eq. / yr
        GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
        sales = (
            biodiesel_flow() * mean_biodiesel_price
            + ethanol_flow() * mean_ethanol_price
            + crude_glycerol_flow() * mean_glycerol_price
            + max(-electricity(), 0) * mean_electricity_price
        )
        return (GWP_material + GWP_emissions) / sales

    @metric(name='Ethanol GWP', element='Economic allocation', units='kg*CO2*eq / gal')
    def GWP_ethanol(): # Cradle to gate
        return GWP_economic.get() * mean_ethanol_price
    
    @metric(name='Biodiesel GWP', element='Economic allocation', units='kg*CO2*eq / gal')
    def GWP_biodiesel(): # Cradle to gate
        if number > 0:
            return GWP_economic.get() * mean_biodiesel_price
        else:
            return 0.
    
    @metric(name='Crude glycerol GWP', element='Economic allocation', units='kg*CO2*eq / kg')
    def GWP_crude_glycerol(): # Cradle to gate
        if number > 0:
            return GWP_economic.get() * mean_glycerol_price
        else:
            return 0.
    
    @metric(name='Electricity GWP', element='Economic allocation', units='kg*CO2*eq / MWhr')
    def GWP_electricity(): # Cradle to gate
        if abs(number) == 1:
            return GWP_economic.get() * mean_electricity_price * 1000.
        else:
            return 0.

    @metric(name='Ethanol GWP', element='Displacement allocation', units='kg*CO2*eq / gal')
    def GWP_ethanol_displacement(): # Cradle to gate
        GWP_material = sys.get_total_feeds_impact(GWP)
        GWP_electricity_production = GWP_characterization_factors['Electricity'] * electricity_production.get() * feedstock_consumption.get()
        GWP_coproducts = sys.get_total_products_impact(GWP)
        GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
        GWP_total = GWP_material + GWP_emissions - GWP_electricity_production - GWP_coproducts # kg CO2 eq. / yr
        return GWP_total / (ethanol_production.get() * feedstock_consumption.get())
    
    # import thermosteam as tmo
    # glycerol = tmo.Chemical('Glycerol')
    # ethanol = tmo.Chemical('Ethanol')
    # glycerol_GGE = 0.80 * (glycerol.LHV / glycerol.MW) / 121300 # 0.1059 GGE / kg crude-glycerol
    
    @metric(name='Biofuel GWP', element='Energy allocation', units='kg*CO2*eq / GGE')
    def GWP_biofuel_allocation(): # Cradle to gate
        GWP_material = sys.get_total_feeds_impact(GWP)
        GWP_emissions = sys.get_process_impact(GWP) # kg CO2 eq. / yr
        GWP_total = GWP_material + GWP_emissions # kg CO2 eq. / yr
        GGE_biodiesel_annual = (biodiesel_production.get() * feedstock_consumption.get()) / 0.9536
        GGE_ethanol_annual = (ethanol_production.get() * feedstock_consumption.get()) / 1.5
        GEE_electricity_production = max(-electricity() * 3600 / 114000, 0.) 
        GEE_crude_glycerol = crude_glycerol_flow() * 0.1059
        return GWP_total / (GGE_biodiesel_annual + GGE_ethanol_annual + GEE_electricity_production + GEE_crude_glycerol)
    
    @metric(name='Ethanol GWP', element='Energy allocation', units='kg*CO2*eq / gal')
    def GWP_ethanol_allocation(): # Cradle to gate
        return GWP_biofuel_allocation.get() / 1.5
    
    @metric(name='Biodiesel GWP', element='Energy allocation', units='kg*CO2*eq / gal')
    def GWP_biodiesel_allocation(): # Cradle to gate
        if number > 0:
            return GWP_biofuel_allocation.get() / 0.9536
        else:
            return 0.
    
    @metric(name='Crude-glycerol GWP', element='Energy allocation', units='kg*CO2*eq / kg')
    def GWP_crude_glycerol_allocation(): # Cradle to gate
        if number > 0:
            return GWP_biofuel_allocation.get() * 0.1059
        else:
            return 0.

    @metric(units='USD/ton')
    def MFPP_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        if agile:
            cane_mode.oil_content += 0.01
            sorghum_mode.oil_content += 0.01
        else:
            oil_extraction_specification.load_oil_content(oil_extraction_specification.oil_content + 0.01)
        sys.simulate()  
        # value = (kg_per_ton * tea.solve_price(feedstock) - MFPP.cache)
        # feedstock.price = tea.solve_price(feedstock)
        # print('AFTER')
        # print('MFPP', kg_per_ton * tea.solve_price(feedstock))
        # print('VOC', tea.VOC / 1e3)
        # print('TCI', tea.TCI / 1e6)
        # print('sales', tea.sales / 1e3)
        # print('NPV', tea.NPV)
        return MFPP.difference()
    
    @metric(units='Gal/ton')
    def biodiesel_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return biodiesel_production.difference()
    
    @metric(units='Gal/ton')
    def ethanol_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return ethanol_production.difference()
    
    @metric(units='kWhr/ton')
    def electricity_production_derivative():
        if number < 0: return 0.
        if _derivative_disabled: return np.nan
        return electricity_production.difference()
    
    @metric(units='cf/ton')
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
        if _derivative_disabled: return np.nan
        return GWP_economic.difference()

    @metric(name='Ethanol GWP derivative', element='Ethanol', units='kg*CO2*eq / gal')
    def GWP_ethanol_derivative(): # Cradle to gate
        return GWP_economic_derivative.get() * mean_ethanol_price
    
    @metric(name='Biodiesel GWP derivative', element='Biodiesel', units='kg*CO2*eq / gal')
    def GWP_biodiesel_derivative(): # Cradle to gate
        if number > 0:
            return GWP_economic_derivative.get() * mean_biodiesel_price
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
        if abs(number) == 1:
            return GWP_economic_derivative.get() * mean_electricity_price * 1000.
        else:
            return 0.
    
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
    
    if abs(number) == 2:
        if enhanced_cellulosic_performance:
            set_baseline(set_sorghum_glucose_yield, 97.5)
            set_baseline(set_sorghum_xylose_yield, 97.5)
            set_baseline(set_cane_glucose_yield, 97.5)
            set_baseline(set_cane_xylose_yield, 97.5)
            set_baseline(set_glucose_to_ethanol_yield, 95)
            set_baseline(set_xylose_to_ethanol_yield, 95)
            set_baseline(set_cofermentation_titer, 120.)
            set_baseline(set_cofermentation_productivity, 2.0)
        else:
            set_baseline(set_sorghum_glucose_yield, 79)
            set_baseline(set_sorghum_xylose_yield, 86)
            set_baseline(set_cane_glucose_yield, 91.0)
            set_baseline(set_cane_xylose_yield, 97.5)
            set_baseline(set_glucose_to_ethanol_yield, 90)
            set_baseline(set_xylose_to_ethanol_yield, 42)
    if number == 1 and enhanced_biodiesel_production:
        set_baseline(set_cane_oil_content, 15)
        set_baseline(set_bagasse_oil_extraction_efficiency, oil_extraction_efficiency_hook(20))
        set_baseline(set_bagasse_oil_retention, 40)
    else:
        set_baseline(set_cane_oil_content, 5)
        set_baseline(set_bagasse_oil_extraction_efficiency, oil_extraction_efficiency_hook(0.))
        set_baseline(set_bagasse_oil_retention, 70)
    set_baseline(set_ethanol_price, mean_ethanol_price) 
    set_baseline(set_crude_glycerol_price, mean_glycerol_price)
    set_baseline(set_biodiesel_price, mean_biodiesel_price)
    set_baseline(set_natural_gas_price, mean_natural_gas_price)
    set_baseline(set_electricity_price, mean_electricity_price)
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
    
    for i in model._parameters:
        dct[i.setter.__name__] = i
    for i in model._metrics:
        dct[i.getter.__name__] = i
    cache[key] = dct.copy()
    
    ## Simulation
    try: 
        sys.simulate()
    except Exception as e:
        raise e
    if reduce_chemicals:
        oilcane_sys.reduce_chemicals()
    oilcane_sys._load_stream_links()
    HXN.force_ideal_thermo = True
    HXN.cache_network = True
    HXN.simulate()

# DO NOT DELETE: For removing ylabel and yticklabels and combine plots
# import biorefineries.oilcane as oc
# import matplotlib.pyplot as plt
# oc.plot_configuration_breakdown('O2')
# ax, *_ = plt.gcf().get_axes()
# yticks = ax.get_yticks()
# plt.yticks(yticks, ['']*len(yticks))
# plt.ylabel('')
# plt.show()

# DO NOT DELETE: For better TEA tickmarks
# import biorefineries.oilcane as oc
# import numpy as np
# oc.plot_monte_carlo(derivative=True, comparison=False,
#     tickmarks=np.array([
#         [-3.0, -1.5, 0, 1.5, 3.0, 4.5],
#         [-6, -3, 0, 3, 6, 9],
#         [-2.25, -1.5, -0.75, 0, 0.75, 1.5],
#         [-10, 0, 10, 20, 30, 40],
#         [-300, -225, -150, -75, 0, 75]
#     ]),
#     labels=['Conventional', 'Cellulosic']
# )

# DO NOT DELETE: For better LCA tickmarks
# import biorefineries.oilcane as oc
# import numpy as np
# oc.plot_monte_carlo(kind='LCA', comparison=True, agile=False,
#     tickmarks=np.array([
#         [-2, -1, 0, 1, 2, 3, 4, 5],
#         [-2, 0, 2, 4, 6, 8, 10],
#         [-1, 0, 1, 2, 3, 4, 5],
#         [-150, -75, -50, -25, 0., 75, 150, 225, 300],
#     ]),
# )

# DO NOT DELETE: For SI Monte Carlo
# import biorefineries.oilcane as oc
# import numpy as np
# oc.plot_monte_carlo(
#     absolute=True, comparison=False,
#     labels=['Sugarcane\nconventional', 'Oilcane\nconventional',
#             'Sugarcane\ncellulosic', 'Oilcane\ncellulosic',
#             'Sugarcane\nconventional\nagile', 'Oilcane\nconventional\nagile',
#             'Sugarcane\ncellulosic\nagile', 'Oilcane\ncellulosic\nagile'],
# )

# # DO NOT DELETE: For GWP tables
# from biorefineries import oilcane as oc
# import biosteam as bst
# def get_sys(ID):
#     oc.load(ID)
#     return oc.sys
# def get_ethanol(ID):
#     oc.load(ID)
#     return oc.ethanol
# IDs = ('S1', 'S2', 'O1', 'O2')
# systems = [get_sys(i) for i in IDs]
# items = [get_ethanol(i) for i in IDs]
# bst.settings.define_impact_indicator('GWP', 'kg*CO2e')
# bst.report.lca_table_displacement_allocation(systems, 'GWP', items, 'ethanol', system_names=IDs)

# # Calculate xylose conversion based on net conversion of sugars
# import biosteam as bst
# import biorefineries.oilcane as oc
# oc.load('L2')
# feed = bst.Stream.sum(oc.R401.ins)
# total_glucose, total_xylose = feed.imass['Glucose', 'Xylose'] + feed.imass['Glucan', 'Xylan'] * 1.11
# total_sugars = total_glucose + total_xylose
# glucose = 0.91 * total_glucose
# xylose = total_xylose
# sugar_conversion = 0.83 * 0.95 
# xylose_conversion = (sugar_conversion * total_sugars - 0.95 * glucose) / xylose
# print(xylose_conversion)

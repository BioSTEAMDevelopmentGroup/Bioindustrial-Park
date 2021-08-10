# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from colorpalette import Palette
from math import floor, ceil
import os
import biosteam as bst
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from biorefineries.sugarcane import create_sugarcane_to_ethanol_system
from biorefineries.lipidcane import (
    create_chemicals as create_starting_chemicals,
)
from thermosteam.units_of_measure import format_units
from biosteam.utils import CABBI_colors, colors
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
from warnings import warn
import pandas as pd
import numpy as np
from typing import NamedTuple
from ..lipidcane import set_lipid_fraction, get_lipid_fraction
from .. import PY37
from . import (_process_settings,
               _chemicals,
               _system,
               _tea,
               _lipid_extraction_specification,
               _distributions,
)
from ._contours import *
from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *
from ._lipid_extraction_specification import *
from ._distributions import *
from biorefineries.sugarcane import create_tea as create_conventional_ethanol_tea
import numpy as np

__all__ = [*_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           *_lipid_extraction_specification.__all__,
           *_distributions.__all__,
           'lipidcane_sys',
           'lipidcane_tea', 
           'flowsheet',
]
_system_loaded = False
_chemicals_loaded = False

derivative_disabled = True
kg_per_ton = 907.18474
lipid_content = np.linspace(0.05, 0.15, 5)

area_colors = {
    'Feedstock handling': CABBI_colors.teal, 
    'Juicing': CABBI_colors.green_dirty,
    'EtOH prod.': CABBI_colors.blue,
    'Oil ext.': CABBI_colors.brown,
    'Biod. prod.': CABBI_colors.orange,
    'Pretreatment': CABBI_colors.green,
    'Wastewater treatment': colors.purple,
    'CH&P': CABBI_colors.yellow,
    'Utilities': colors.red,
    'Storage': CABBI_colors.grey,
    'HXN': colors.orange,
}

area_hatches = {
    'Feedstock handling': 'x', 
    'Juicing': '-',
    'EtOH prod.': '/',
    'Oil ext.': '\\',
    'Biod. prod.': '/|',
    'Pretreatment': '//',
    'Wastewater treatment': r'\\',
    'CH&P': '',
    'Utilities': '\\|',
    'Storage': '',
    'HXN': '+',
}

for i in area_colors: area_colors[i] = area_colors[i].tint(20)
palette = Palette(**area_colors)


configuration_names = (
    'S1', 'L1', 'S2', 'L2', 'S1*', 'L1*', 'S2*', 'L2*',
)
comparison_names = (
    # 'I - ∅', 
    'L1 - S1', 
    'L2 - S2', 
    'L2 - L1', 
    'L1* - L1', 
    'L2* - L2',  
)

other_comparison_names = (
    'L1* - S1*', 'L2* - S2*', 
)

across_lipid_content_names = (
    'L1', 'L2', 
)

across_lipid_content_agile_names = (
    'L1*', 'L2*', 
)

across_lipid_content_comparison_names = (
    'L1 - S1', 'L2 - S2', 'L2 - L1', 
)

across_lipid_content_agile_direct_comparison_names = (
    'L1* - L1', 'L2* - L2', 
)

across_lipid_content_agile_comparison_names = (
    'L1* - S1*', 'L2* - S2*', 'L2* - L1*', 
)

(set_cane_lipid_content, set_relative_sorghum_lipid_content, set_bagasse_lipid_retention, 
 set_bagasse_lipid_extraction_efficiency, 
 set_plant_capacity, set_ethanol_price,
 set_biodiesel_price, set_natural_gas_price, 
 set_electricity_price, set_operating_days, 
 set_IRR) = all_parameter_mockups = (
    bst.MockVariable('Cane lipid content', 'dry wt. %', 'Stream-sugarcane'),
    bst.MockVariable('Relative sorghum lipid content', 'dry wt. %', 'Stream-sugarcane'),
    bst.MockVariable('Lipid retention', '%', 'Stream-sugarcane'),
    bst.MockVariable('Bagasse lipid extraction efficiency', '%', 'Stream-sugarcane'),
    bst.MockVariable('Capacity', 'ton/hr', 'Stream-sugarcane'),
    bst.MockVariable('Price', 'USD/gal', 'Stream-ethanol'),
    bst.MockVariable('Price', 'USD/gal', 'Stream-biodiesel'),
    bst.MockVariable('Price', 'USD/cf', 'Stream-natural gas'),
    bst.MockVariable('Electricity price', 'USD/kWh', 'biorefinery'),
    bst.MockVariable('Operating days', 'day/yr', 'biorefinery'),
    bst.MockVariable('IRR', '%', 'biorefinery'),
)
(MFPP, biodiesel_production, ethanol_production, 
 electricity_production, natural_gas_consumption, 
 TCI, feedstock_consumption, heat_exchanger_network_error,
 MFPP_derivative, biodiesel_production_derivative, ethanol_production_derivative, 
 electricity_production_derivative, natural_gas_consumption_derivative, 
 TCI_derivative) = all_metric_mockups = (
    bst.MockVariable('MFPP', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('Biodiesel production', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Ethanol production', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Electricity production', 'kWhr/ton', 'Biorefinery'),
    bst.MockVariable('Natural gas consumption', 'cf/ton', 'Biorefinery'),
    bst.MockVariable('TCI', '10^6*USD', 'Biorefinery'),
    bst.MockVariable('Feedstock consumption', 'ton/yr', 'Biorefinery'),
    bst.MockVariable('Heat exchanger network error', '%', 'Biorefinery'),
    bst.MockVariable('MFPP derivative', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('Biodiesel production derivative', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Ethanol production derivative', 'Gal/ton', 'Biorefinery'),
    bst.MockVariable('Electricity production derivative', 'kWhr/ton', 'Biorefinery'),
    bst.MockVariable('Natural gas consumption derivative', 'cf/ton', 'Biorefinery'),
    bst.MockVariable('TCI derivative', '10^6*USD', 'Biorefinery'),
)

metric_mockups = (
    MFPP, 
    TCI,
    ethanol_production,
    biodiesel_production,
    electricity_production,
    natural_gas_consumption
)

derivative_metric_mockups = (
    MFPP_derivative, 
    TCI_derivative,
    ethanol_production_derivative,
    biodiesel_production_derivative,
    electricity_production_derivative,
    natural_gas_consumption_derivative,
)

def asconfiguration(x):
    number, agile = x
    return Configuration(int(number), bool(agile))

def ascomparison(x):
    a, b = x
    return ConfigurationComparison(asconfiguration(a), asconfiguration(b))

class Configuration(NamedTuple):
    number: int
    agile: bool = False
   
class ConfigurationComparison(NamedTuple):
    a: Configuration
    b: Configuration

def name_to_configuration(name):
    name = name.replace(' ', '')
    return Configuration((-1 if name.startswith('S') else 1) * int(name[1]), '*' in name)

def parse(x):
    try:
        return asconfiguration(x)
        return ascomparison(x)
    except:
        pass
    if isinstance(x, int):
        return Configuration(x)
    elif isinstance(x, str):
        x = x.upper()
        left, *right = x.split('-')
        if right:
            if len(right) == 1:
                right = right[0]
                return ConfigurationComparison(
                    name_to_configuration(left),
                    name_to_configuration(right)
                )
            else:
                raise RuntimeError('cannot parse multiple subtractions')
        else:
            factor = -1 if x.startswith('S') else 1
            return Configuration(factor * int(x[1]), '*' in x)
    raise ValueError(f'could not parse {x}')
    

def format_name(name):
    key = parse(name)
    if isinstance(key, Configuration):
        return format_configuration(key)
    elif isinstance(key, ConfigurationComparison):
        return format_comparison(key)
    else:
        raise Exception('unknown error')

def format_configuration(configuration):
    number, agile = configuration
    if number == -2:
        name = 'S2'
    elif number == -1:
        name = 'S1'
    elif number == 0:
        name = '∅'
    elif number == 1:
        name = 'L1'
    elif number == 2:
        name = 'L2'
    else:
        raise ValueError(f'invalid configuration {configuration}')
    name = r'$\mathtt{' + name + '}$'
    if agile: name += '*'
    return name

def format_comparison(comparison):
    return ' $-$ '.join([format_configuration(i) for i in comparison])

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def sorghum_feedstock(sugarcane, ID):
    return bst.Stream(
        ID=ID, phase='l', T=298.15, P=101325, 
        Water=2.333e+05, Glucose=3703, Sucrose=4.196e+04, Ash=2000, 
        Cellulose=2.227e+04, Hemicellulose=1.314e+04, Lignin=1.193e+04, 
        Solids=5000, units='kg/hr'
    )

def load(name, cache={}):
    dct = globals()
    number, agile = dct['configuration'] = parse(name)
    if (number, agile) in cache:
        dct.update(cache[number, agile])
        return
    global lipidcane_sys, lipidcane_tea, specs, flowsheet, _system_loaded
    global lipid_extraction_specification, model, unit_groups
    global HXN, BT
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('lipidcane2g')
    main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    u = flowsheet.unit
    s = flowsheet.stream
    operating_hours = 24 * 200
    
    ## System
    
    area_names = None
    def rename_storage_units(storage):
        bst.rename_units([i for i in lipidcane_sys.units if bst.is_storage_unit(i)], storage)
    
    if number == -1:
        isplit_efficiency_is_reversed = None
        # starting_chemicals = create_starting_chemicals()
        # bst.settings.set_thermo(starting_chemicals)
        lipidcane_sys = create_sugarcane_to_ethanol_system(
            operating_hours=operating_hours,
            use_area_convention=True,
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
        lipidcane_sys = create_sugarcane_to_ethanol_combined_1_and_2g(
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
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_1g(
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
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1100)
    else:
        raise NotImplementedError(number)
    lipidcane_sys.set_tolerance(rmol=1e-3)
    lipidcane_sys.N_runs = 3
    dct.update(flowsheet.to_dict())
    if number == 1:
        lipidcane_sys.prioritize_unit(T608)
    elif number == 2:
        lipidcane_sys.prioritize_unit(T808)
    if number < 0:
        dct['lipidcane'] = sugarcane
        dct['lipidcane_sys'] = lipidcane_sys
    unit_groups = UnitGroup.group_by_area(lipidcane_sys.units)
    if area_names:
        for i, j in zip(unit_groups, area_names): i.name = j
    for i in unit_groups: i.autofill_metrics(shorthand=True)
    
    for BT in lipidcane_sys.units:
        if isinstance(BT, bst.BoilerTurbogenerator): break

    HXN = None
    for HXN_group in unit_groups:
        if HXN_group.name == 'HXN':
            HXN_group.filter_savings = False
            HXN = HXN_group.units[0]
            assert isinstance(HXN, bst.HeatExchangerNetwork)
    
    unit_groups[-1].metrics[-1].getter = lambda: 0.    
    
    if agile: 
        dct['lipidsorghum'] = lipidsorghum = sorghum_feedstock(lipidcane, ID='lipidsorghum')
        original_lipidcane = lipidcane.copy()
        
        class AgileLipidcaneSystem(bst.AgileSystem):
            
            def set_parameters(self, switch2sorghum):
                F_mass_original = lipidcane.F_mass
                if switch2sorghum:
                    lipidcane.copy_like(lipidsorghum)
                    if number >= 0: set_lipid_fraction(self.sorghum_lipid_content, lipidcane)
                else:
                    lipidcane.copy_like(original_lipidcane)
                    if number >= 0: set_lipid_fraction(self.cane_lipid_content, lipidcane)
                lipidcane.F_mass = F_mass_original
                
        lipidcane_tea = create_agile_tea(lipidcane_sys.units)      
        lipidcane_sys = AgileLipidcaneSystem(lipidcane_sys, [0, 1], [200 * 24, 60 * 24],
                                             tea=lipidcane_tea)
        lipidcane_sys.cane_lipid_content = 0.1
        lipidcane_sys.sorghum_lipid_content = 0.07
    else:
        lipidcane_tea = create_tea(lipidcane_sys)
        lipidcane_sys.operating_hours = 24 * 200
    lipidcane_tea.income_tax = 0.21 # Davis et al. 2018; https://www.nrel.gov/docs/fy19osti/71949.pdf
    
    ## Specification for analysis
    if number < 0:
        dct['biodiesel'] = bst.Stream('biodiesel')
        isplit_a = None
        isplit_b = None
        lipid_extraction_specification = MockExtractionSpecification()
    else:
        for i in lipidcane_sys.units:
            if getattr(i, 'tag', None) == 'lipid extraction efficiency':
                isplit_a = i.isplit
                break
        
        for i in lipidcane_sys.units:
            if getattr(i, 'tag', None) == 'bagasse lipid retention':
                isplit_b = i.isplit
                break
        
        lipid_extraction_specification = LipidExtractionSpecification(
            lipidcane_sys, lipidcane, isplit_a, isplit_b, isplit_efficiency_is_reversed
        )
    
    ## Model
    
    model = bst.Model(lipidcane_sys, exception_hook='raise', retry_evaluation=False)
    parameter = model.parameter
    metric = model.metric
    
    def uniform(lb, ub, *args, **kwargs):
        return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), **kwargs)
    
    def triangular(lb, mid, ub, *args, **kwargs):
        return parameter(*args, distribution=shape.Triangle(lb, mid, ub), bounds=(lb, ub), **kwargs)
    
    # Currently at ~5%, but total lipid content is past 10%
    
    @uniform(5., 15., element=lipidcane, units='dry wt. %', kind='coupled')
    def set_cane_lipid_content(cane_lipid_content):
        if agile:
            lipidcane_sys.cane_lipid_content = cane_lipid_content / 100.
        else:
            lipid_extraction_specification.load_lipid_content(cane_lipid_content / 100.)
    
    @uniform(-3., 0., element=lipidcane, units='dry wt. %', kind='coupled')
    def set_relative_sorghum_lipid_content(relative_sorghum_lipid_content):
        if agile:
            lipidcane_sys.sorghum_lipid_content = lipidcane_sys.cane_lipid_content + relative_sorghum_lipid_content / 100.
    
    @uniform(65, 75., element=lipidcane, units='%', kind='coupled')
    def set_bagasse_lipid_retention(lipid_retention):
        lipid_extraction_specification.load_lipid_retention(lipid_retention / 100.)
    
    def lipid_extraction_efficiency_hook(x):
        if number < 0:
            return x
        elif number == 1:
            return 50.0 + x
        elif number == 2:
            return 70.0 + x
    
    @uniform(0., 20, units='%', element=lipidcane, kind='coupled', 
             hook=lipid_extraction_efficiency_hook)
    def set_bagasse_lipid_extraction_efficiency(bagasse_lipid_extraction_efficiency):
        lipid_extraction_specification.load_efficiency(bagasse_lipid_extraction_efficiency / 100.)

    capacity = lipidcane.F_mass / kg_per_ton
    @uniform(0.9 * capacity, 1.1 * capacity, units='ton/hr',
             element=lipidcane, kind='coupled')
    def set_plant_capacity(capacity):
        lipidcane.F_mass = kg_per_ton * capacity

    # USDA ERS historical price data
    @parameter(distribution=ethanol_price_distribution, element=ethanol, units='USD/gal')
    def set_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
        ethanol.price = price / 2.98668849

    # USDA ERS historical price data
    @parameter(distribution=biodiesel_minus_ethanol_price_distribution, element=biodiesel, units='USD/gal',
               hook=lambda x: ethanol.price * 2.98668849 + x)
    def set_biodiesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        biodiesel.price = price / 3.3111

    # https://www.eia.gov/energyexplained/natural-gas/prices.php
    @parameter(distribution=natural_gas_price_distribution, element=natural_gas, units='USD/cf')
    def set_natural_gas_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        BT.natural_gas_price = 51.92624700383502 * price / 1000. 

    # https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
    # Data from historical prices, 2010-2020
    @triangular(0.0583, 0.065, 0.069, units='USD/kWh')
    def set_electricity_price(electricity_price): 
        bst.PowerUtility.price = electricity_price
        
    # From Huang's 2016 paper
    @uniform(6 * 30, 7 * 30, units='day/yr')
    def set_operating_days(operating_days):
        if agile:
            lipidcane_sys.time_steps[0] = operating_days * 24
        else:
            lipidcane_tea.operating_days = operating_days
    
    # 10% is suggested for waste reducing, but 15% is suggested for investment
    @uniform(10., 15., units='%')
    def set_IRR(IRR):
        lipidcane_tea.IRR = IRR / 100.

    # # Initial solids loading in Humbird is 20% (includes both dissolved and insoluble solids)
    # @uniform(20, 25, units='%', element='Fermentation')
    # def set_fermentation_solids_loading(solids_loading):
    #     if number == 2: M402.solids_loading = solids_loading / 100.
    
    natural_gas.phase = 'g'
    natural_gas.set_property('T', 60, 'degF')
    natural_gas.set_property('P', 14.73, 'psi')
    original_value = natural_gas.imol['CH4']
    natural_gas.imass['CH4'] = 1 
    V_ng = natural_gas.get_total_flow('ft3/hr')
    natural_gas.imol['CH4'] = original_value
    
    if agile:
        feedstock_flow = lambda: lipidcane_sys.flow_rates[lipidcane] / kg_per_ton # ton/yr
        biodiesel_flow = lambda: lipidcane_sys.flow_rates.get(biodiesel, 0.) / 3.3111 # gal/yr
        ethanol_flow = lambda: lipidcane_sys.flow_rates[ethanol] / 2.98668849 # gal/yr
        natural_gas_flow = lambda: lipidcane_sys.flow_rates[natural_gas] * V_ng # cf/yr
        if number <= 1:
            electricity = lambda: lipidcane_sys.utility_cost / bst.PowerUtility.price
        elif number == 2:
            electricity = lambda: 0.
    else:
        feedstock_flow = lambda: lipidcane_sys.operating_hours * lipidcane.F_mass / kg_per_ton # ton/yr
        biodiesel_flow = lambda: lipidcane_sys.operating_hours * biodiesel.F_mass / 3.3111 # gal/yr
        ethanol_flow = lambda: lipidcane_sys.operating_hours * ethanol.F_mass / 2.98668849 # gal/yr
        natural_gas_flow = lambda: lipidcane_sys.operating_hours * natural_gas.F_mass * V_ng # cf/yr
        if number <= 1:
            electricity = lambda: lipidcane_sys.operating_hours * sum([i.rate for i in lipidcane_sys.power_utilities])
        elif number == 2:
            electricity = lambda: 0.
    dct['flows'] = flows = {'feedstock': feedstock_flow,
                            'biodiesel': biodiesel_flow,
                            'ethanol': ethanol_flow,
                            'natural_gas': natural_gas_flow,
                            'electricity': electricity}
    @metric(units='USD/ton')
    def MFPP():
        lipidcane.price = price = lipidcane_tea.solve_price(lipidcane)
        return kg_per_ton * price
    
    @metric(units='Gal/ton')
    def biodiesel_production():
        return biodiesel_flow() / feedstock_flow()
    
    @metric(units='Gal/ton')
    def ethanol_production():
        return ethanol_flow() / feedstock_flow()
    
    @metric(units='kWhr/ton')
    def electricity_production():
        value = - electricity() / feedstock_flow() 
        if value < 0.: value = 0.
        return value
    
    @metric(units='cf/ton')
    def natural_gas_consumption():
        return natural_gas_flow() / feedstock_flow()
    
    @metric(units='10^6*USD')
    def TCI():
        return lipidcane_tea.TCI / 1e6 # 10^6*$
    
    @metric(units='ton/yr')
    def feedstock_consumption():
        return feedstock_flow()
    
    @metric(units='%')
    def heat_exchanger_network_error():
        return HXN.energy_balance_percent_error if HXN else 0.    
    
    # @metric(units='MMGGE/yr')
    # def productivity():
    #     GGE = (ethanol_flow() / 1.5
    #        + biodiesel_flow() / 0.9536
    #        - electricity() * 3600 / 131760
    #        - natural_gas_flow() / 126.67)
    #     return GGE / 1e6

    @metric(units='USD/ton')
    def MFPP_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        if agile:
            lipidcane_sys.cane_lipid_content += 0.01
            lipidcane_sys.sorghum_lipid_content += 0.01
        else:
            lipid_extraction_specification.load_lipid_content(lipid_extraction_specification.lipid_content + 0.01)
        lipidcane_sys.simulate()  
        # value = (kg_per_ton * lipidcane_tea.solve_price(lipidcane) - MFPP.cache) / 1.
        # lipidcane.price = lipidcane_tea.solve_price(lipidcane)
        # print('AFTER')
        # print('MFPP', kg_per_ton * lipidcane_tea.solve_price(lipidcane))
        # print('VOC', lipidcane_tea.VOC / 1e3)
        # print('TCI', lipidcane_tea.TCI / 1e6)
        # print('sales', lipidcane_tea.sales / 1e3)
        # print('NPV', lipidcane_tea.NPV)
        return (kg_per_ton * lipidcane_tea.solve_price(lipidcane) - MFPP.cache) / 1.
    
    @metric(units='Gal/ton')
    def biodiesel_production_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        value = (biodiesel_flow() / feedstock_flow() - biodiesel_production.cache) / 1.
        # print('biodiesel production derivative', value)
        return value
    
    @metric(units='Gal/ton')
    def ethanol_production_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        value = (ethanol_flow() / feedstock_flow() - ethanol_production.cache) / 1.
        # print('ethanol production derivative', value)
        return value
    
    @metric(units='kWhr/ton')
    def electricity_production_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        value = (- electricity() / feedstock_flow() - electricity_production.cache) / 1.
        # print('electricity production derivative', value)
        return value
    
    @metric(units='cf/ton')
    def natural_gas_consumption_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        value =(natural_gas_flow() / feedstock_flow() - natural_gas_consumption.cache) / 1.
        # print('natural gas production derivative', value)
        return value
    
    @metric(units='10^6*USD')
    def TCI_derivative():
        if number < 0: return 0.
        if derivative_disabled: return None
        value = (lipidcane_tea.TCI / 1e6  - TCI.cache) / 1. # 10^6*$
        # print('TCI production derivative', value)
        return value
        
    
    
    # # Single point evaluation for detailed design results
    lipid_extraction_specification.load_lipid_retention(0.70)
    lipid_extraction_specification.load_lipid_content(0.05)
    set_bagasse_lipid_extraction_efficiency(0.)
    set_ethanol_price.setter(1.898) 
    set_biodiesel_price.setter(4.363)
    set_natural_gas_price.setter(4.3)
    set_electricity_price.setter(0.0641)
    if number == 2: M402.solids_loading = 0.20
    # set_fermentation_solids_loading(20) # Same as Humbird
    # set_feedstock_lipid_content(10) # Consistent with SI of Huang's 2016 paper
    # set_ethanol_price(2.356) # Consistent with Huang's 2016 paper
    # set_biodiesel_price(4.569) # Consistent with Huang's 2016 paper
    # set_natural_gas_price(4.198) # Consistent with Humbird's 2012 paper
    # set_electricity_price(0.0572) # Consistent with Humbird's 2012 paper
    set_operating_days(200) # Consistent with Huang's 2016 paper
    
    for i in model.parameters:
        dct[i.setter.__name__] = i
    for i in model.metrics:
        dct[i.getter.__name__] = i
    cache[number, agile] = dct.copy()
    
    ## Simulation
    try: 
        lipidcane_sys.simulate()
    except Exception as e:
        raise e
    finally:
        lipidcane_tea.IRR = 0.10

def evaluate_configurations_across_extraction_efficiency_and_lipid_content(
        efficiency, lipid_content, lipid_retention, agile, configurations,
    ):
    A = len(agile)
    C = len(configurations)
    M = len(all_metric_mockups)
    data = np.zeros([A, C, M])
    for ia in range(A):
        for ic in range(C):    
            load([int(configurations[ic]), agile[ia]])
            lipid_extraction_specification.load_specifications(
                efficiency=efficiency, 
                lipid_content=lipid_content, 
                lipid_retention=lipid_retention
            )
            if agile[ia]:
                lipidcane_sys.cane_lipid_content = lipidcane_sys.sorghum_lipid_content = lipid_content
            lipidcane_sys.simulate()
            data[ia, ic, :] = [j() for j in model.metrics]
    return data

N_metrics = len(all_metric_mockups)
evaluate_configurations_across_extraction_efficiency_and_lipid_content = np.vectorize(
    evaluate_configurations_across_extraction_efficiency_and_lipid_content, 
    excluded=['lipid_retention', 'agile', 'configurations'],
    signature=f'(),(),(),(a),(c)->(a,c,{N_metrics})'
)

def evaluate_configurations_across_sorghum_and_cane_lipid_content(
        sorghum_lipid_content, cane_lipid_content, configurations, relative,
    ):
    C = len(configurations)
    M = len(all_metric_mockups)
    data = np.zeros([C, M])
    for ic in range(C):
        load([int(configurations[ic]), True])
        lipidcane_sys.cane_lipid_content = cane_lipid_content
        if relative:
            lipidcane_sys.sorghum_lipid_content = cane_lipid_content + sorghum_lipid_content
        else:
            lipidcane_sys.sorghum_lipid_content = sorghum_lipid_content
        lipidcane_sys.simulate()
        data[ic, :] = [j() for j in model.metrics]
    return data

evaluate_configurations_across_sorghum_and_cane_lipid_content = np.vectorize(
    evaluate_configurations_across_sorghum_and_cane_lipid_content, 
    excluded=['configurations', 'relative'],
    signature=f'(),(),(c),()->(c,{N_metrics})'
)              

def evaluate_MFPP_uncertainty_across_ethanol_and_biodiesel_prices(name, ethanol_price, biodiesel_price):
    table = get_monte_carlo(name)
    lipidcane_price = table[MFPP.index].to_numpy()[:, np.newaxis] # USD/ton
    biodiesel_flow = table[biodiesel_production.index].to_numpy()[:, np.newaxis] * 1e6 # gal/yr
    ethanol_price_baseline = table[set_ethanol_price.index].to_numpy()[:, np.newaxis]
    biodiesel_price_baseline = table[set_biodiesel_price.index].to_numpy()[:, np.newaxis]
    ethanol_flow = table[ethanol_production.index].to_numpy()[:, np.newaxis] * 1e6 # gal/yr
    feedstock_flow = table[feedstock_consumption.index].to_numpy()[:, np.newaxis] # ton/yr
    baseline_price = (
        lipidcane_price
        - (ethanol_price_baseline * ethanol_flow + biodiesel_price_baseline * biodiesel_flow) / feedstock_flow
    )
    return (
        baseline_price 
        + (ethanol_price[np.newaxis, :] * ethanol_flow + biodiesel_price[np.newaxis, :] * biodiesel_flow) / feedstock_flow
    )

def evaluate_MFPP_benefit_uncertainty_across_ethanol_and_biodiesel_prices(name, ethanol_price, biodiesel_price, baseline=None):
    if baseline is None:
        configuration = parse(name)
        number, agile = configuration
        assert number > 0
        baseline = Configuration(-number, agile)
    MFPP_baseline = evaluate_MFPP_uncertainty_across_ethanol_and_biodiesel_prices(baseline, ethanol_price, biodiesel_price)
    MFPP = evaluate_MFPP_uncertainty_across_ethanol_and_biodiesel_prices(name, ethanol_price, biodiesel_price)
    return MFPP - MFPP_baseline

def evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, configuration=None):
    if configuration is not None: load(configuration)
    feedstock_flow = flows['feedstock']()
    biodiesel_flow = flows['biodiesel']()
    ethanol_flow = flows['ethanol']()
    baseline_price = (
        lipidcane_tea.solve_price(lipidcane) * kg_per_ton
        - (ethanol.price * ethanol_flow  * 2.98668849 + biodiesel.price * 3.3111 * biodiesel_flow) / feedstock_flow
    )
    return (
        baseline_price 
        + (ethanol_price * ethanol_flow + biodiesel_price * biodiesel_flow) / feedstock_flow
    )

def evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, baseline=None, configuration=None):
    if configuration is None: configuration = globals()['configuration']
    if baseline is None:
        number, agile = configuration
        assert number > 0
        baseline = Configuration(-number, agile)
    MFPP_baseline = evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, baseline)
    MFPP = evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, configuration)
    return MFPP - MFPP_baseline

def spearman_file(name):
    number, agile = parse(name)
    folder = os.path.dirname(__file__)
    filename = f'lipidcane_spearman_{number}'
    if agile: filename += '_agile'
    filename += '.xlsx'
    return os.path.join(folder, filename)

def monte_carlo_file(name, across_lipid_content=False):
    number, agile = parse(name)
    folder = os.path.dirname(__file__)
    filename = f'lipidcane_monte_carlo_{number}'
    if agile: filename += '_agile'
    if across_lipid_content: filename += '_across_lipid_content'
    filename += '.xlsx'
    return os.path.join(folder, filename)

def run_uncertainty_and_sensitivity(name, N, rule='L',
                                    across_lipid_content=False, 
                                    sample_cache={}):
    np.random.seed(1)
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.utils.DesignWarning)
    load(name)
    key = (N, rule)
    if key in sample_cache:
        samples = sample_cache[key]
    else:
        sample_cache[key] = samples = model.sample(N, rule)
    model.load_samples(samples)
    file = monte_carlo_file(name, across_lipid_content)
    if across_lipid_content:
        if parse(name).number < 0:
            model.evaluate(notify=int(N/10))
            model.evaluate_across_coordinate(
                name='Lipid content',
                f_coordinate=lambda x: None,
                coordinate=lipid_content,
                notify=int(N/10), 
                notify_coordinate=True,
                xlfile=file,
                re_evaluate=False,
            )
        else:
            def f(x):
                lipid_extraction_specification.locked_lipid_content = False
                lipid_extraction_specification.load_lipid_content(x)
                lipid_extraction_specification.locked_lipid_content = True
            model.evaluate_across_coordinate(
                name='Lipid content',
                f_coordinate=f,
                coordinate=lipid_content,
                notify=int(N/10), 
                notify_coordinate=True,
                xlfile=file,
            )
    else:
        model.evaluate(notify=int(N/10))
        model.table.to_excel(file)
        rho, p = model.spearman_r()
        file = spearman_file(name)
        rho.to_excel(file)
    
def run_all(N, across_lipid_content=False, rule='L', configurations=None):
    if configurations is None: configurations = configuration_names
    for name in configurations:
        print(f"Running {name}:")
        run_uncertainty_and_sensitivity(
            name, N, rule, across_lipid_content,
        )

def get_monte_carlo_across_lipid_content(name, metric, derivative=False):
    key = parse(name)
    if isinstance(key, Configuration):
        df = pd.read_excel(
            monte_carlo_file(key, True),
            sheet_name=metric if isinstance(metric, str) else metric.short_description,
            index_col=0
        )
    elif isinstance(key, ConfigurationComparison):
        df = (
            get_monte_carlo_across_lipid_content(key.a, metric)
            - get_monte_carlo_across_lipid_content(key.b, metric)
        )
    else:
        raise Exception('unknown error')
    if derivative: 
        arr = np.diff(df.values) / np.diff(df.columns.values) / 100.
    else:
        arr = df.values
    return arr
        

def get_monte_carlo(name):
    key = parse(name)
    if isinstance(key, Configuration):
        file = monte_carlo_file(key)
        return pd.read_excel(file, header=[0, 1], index_col=[0])
    elif isinstance(key, ConfigurationComparison):
        index = [i.index for i in metric_mockups + derivative_metric_mockups]
        df_a = get_monte_carlo(key.a)[index]
        df_b = get_monte_carlo(key.b)[index]
        row_a = df_a.shape[0]
        row_b = df_b.shape[0]
        if row_a != row_b:
            length = min(row_a, row_b)
            return df_a.iloc[:length] - df_b.iloc[:length]
        else:
            return df_a - df_b
    else:
        raise Exception('unknown error')

def plot_monte_carlo_across_coordinate(coordinate, data, color_wheel):
    if isinstance(data, list):
        return [plot_monte_carlo_across_coordinate(coordinate, i, color_wheel) for i in data]
    else:
        color = color_wheel.next()
        return bst.plots.plot_montecarlo_across_coordinate(
            coordinate, data,
            light_color=color.tint(50).RGBn,
            dark_color=color.shade(50).RGBn,
        )

def plot_monte_carlo_across_lipid_content(kind=0, derivative=False):
    MFPP, TCI, *production, electricity_production, natural_gas_consumption = metric_mockups
    rows = [MFPP, TCI, production]
    if kind == 0:
        columns = across_lipid_content_names
    elif kind == 1:
        columns = across_lipid_content_agile_names
    elif kind == 2:
        columns = across_lipid_content_comparison_names
    elif kind == 3:
        columns = across_lipid_content_agile_comparison_names
    elif kind == 4:
        columns = across_lipid_content_agile_direct_comparison_names
    else:
        raise NotImplementedError(str(kind))
    if derivative:
        x = 100 * (lipid_content[:-1] + np.diff(lipid_content) / 2.)
        ylabels = [
            f"MFPP der. [{format_units('USD/ton')}]",
            f"TCI der. [{format_units('10^6*USD')}]",
            f"Production der. [{format_units('gal/ton')}]"
        ]
    else:
        x = 100 * lipid_content
        ylabels = [
            f"MFPP$\backprime$ [{format_units('USD/ton')}]",
            f"TCI [{format_units('10^6*USD')}]",
            f"Production [{format_units('gal/ton')}]"
        ]
    N_cols = len(columns)
    N_rows = len(rows)
    fig, axes = plt.subplots(ncols=N_cols, nrows=N_rows)
    data = np.zeros([N_rows, N_cols], dtype=object)
    
    def get_data(metric, name):
        if isinstance(metric, bst.Variable):
            return get_monte_carlo_across_lipid_content(name, metric, derivative)
        else:
            return [get_data(i, name) for i in metric]
    
    data = np.array([[get_data(i, j) for j in columns] for i in rows])
    tickmarks = [None] * N_rows
    get_max = lambda x: max([i.max() for i in x]) if isinstance(x, list) else x.max()
    get_min = lambda x: min([i.min() for i in x]) if isinstance(x, list) else x.min()
    N_ticks = 5
    for r in range(N_rows):
        lb = min(min([get_min(i) for i in data[r, :]]), 0)
        ub = max([get_max(i) for i in data[r, :]])
        diff = 0.1 * (ub - lb)
        ub += diff
        if derivative:
            lb = floor(lb)
            ub = ceil(ub)
            step = (ub - lb) / (N_ticks - 1)
            tickmarks[r] = [0, 1] if step == 0 else [int(lb + step * i) for i in range(N_ticks)]
        else:
            if rows[r] is MFPP:
                if kind == 0 or kind == 1:
                    tickmarks[r] = [-20, 0, 20, 40, 60]
                elif kind == 2:
                    tickmarks[r] = [-20, -10, 0, 10, 20]
                elif kind == 3:
                    tickmarks[r] = [-10, 0, 10, 20, 30]
                elif kind == 4:
                    tickmarks[r] = [-5, 0, 5, 10, 15]
                continue
            lb = floor(lb / 15) * 15
            ub = ceil(ub / 15) * 15
            step = (ub - lb) / (N_ticks - 1)
            tickmarks[r] = [0, 1] if step == 0 else [int(lb + step * i) for i in range(N_ticks)]
    color_wheel = CABBI_colors.wheel()
    for j in range(N_cols):
        color_wheel.restart()
        for i in range(N_rows):
            arr = data[i, j]
            ax = axes[i, j]
            plt.sca(ax)
            percentiles = plot_monte_carlo_across_coordinate(x, arr, color_wheel)
            if i == 0: ax.set_title(format_name(columns[j]))
            xticklabels = i == N_rows - 1
            yticklabels = j == 0
            if xticklabels: plt.xlabel('Lipid content [dry wt. %]')
            if yticklabels: plt.ylabel(ylabels[i])
            bst.plots.style_axis(ax,  
                                 xticks = [5, 10, 15],
                                 yticks = tickmarks[i],
                                 xticklabels= xticklabels, 
                                 yticklabels= yticklabels,
                                 ytick0=False)
    for i in range(N_cols): fig.align_ylabels(axes[:, i])
    plt.subplots_adjust(hspace=0.1, wspace=0.1)

def monte_carlo_box_plot(data, positions, light_color, dark_color):
    return plt.boxplot(x=data, positions=positions, patch_artist=True,
                     widths=0.8, whis=[5, 95],
                     boxprops={'facecolor':light_color,
                               'edgecolor':dark_color},
                     medianprops={'color':dark_color,
                                  'linewidth':1.5},
                     flierprops = {'marker':'D',
                                   'markerfacecolor': light_color,
                                   'markeredgecolor': dark_color,
                                   'markersize':6})

def monte_carlo_results():
    results = {}
    for name in configuration_names + comparison_names + other_comparison_names:
        df = get_monte_carlo(name)
        results[name] = dct = {}
        for metric in metric_mockups + derivative_metric_mockups:
            index = metric.index
            data = df[index].values
            q05, q25, q50, q75, q95 = percentiles = np.percentile(data, [5,25,50,75,95], axis=0)
            dct[index[1]] = subdct = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
    return results

def plot_monte_carlo(derivative=False, configuration_names=configuration_names, comparison_names=comparison_names):
    if derivative:
        # columns = configurations = ['L1', 'L1*', 'L1* - L1']
        configuration_names = ['L1', 'L2', 'L1*', 'L2*']
        comparison_names = ['L2 - L1', 'L2* - L2', 'L1* - L1']
        MFPP, TCI, *production, electricity_production, natural_gas_consumption = derivative_metric_mockups
    else:
        MFPP, TCI, *production, electricity_production, natural_gas_consumption = metric_mockups
    columns = configurations = configuration_names + comparison_names
    N_cols = len(columns)
    rows = metrics = [
        MFPP, 
        TCI, 
        production,
        electricity_production,
        natural_gas_consumption,
    ]
    N_rows = len(rows)
    fig, axes = plt.subplots(ncols=1, nrows=N_rows)
    axes = axes.flatten()
    xtext = [format_name(i).replace(' ', '') for i in configurations]
    N_marks = len(xtext)
    xticks = tuple(range(N_marks))
    color_wheel = CABBI_colors.wheel()
    ylabels = [
        f"MFPP\n[{format_units('USD/ton')}]",
        f"TCI\n[{format_units('10^6*USD')}]",
        f"Production\n[{format_units('Gal/ton')}]",
        f"Elec. prod.\n[{format_units('kWhr/ton')}]",
        f"NG cons.\n[{format_units('cf/ton')}]"
    ]
    if derivative:
        ylabels = [i.replace("\n[", " der. \n[") for i in ylabels]
    def get_data(metric, name):
        if isinstance(metric, bst.Variable):
            df = get_monte_carlo(name)
            values = df[metric.index].values
            return values
        else:
            return [get_data(i, name) for i in metric]
    
    def plot(arr, position):
        if isinstance(arr, list):
            return [plot(i, position) for i in arr]
        else:
            color = color_wheel.next()
            light_color = color.RGBn
            dark_color = color.shade(60).RGBn
            return monte_carlo_box_plot(
                    arr, (position,),
                    light_color=light_color,
                    dark_color=dark_color,
            )
    
    data = np.zeros([N_rows, N_cols], dtype=object)
    data = np.array([[get_data(i, j) for j in columns] for i in rows])
    step_min = 1 if derivative else 30
    tickmarks = [bst.plots.rounded_tickmarks_from_data(i, step_min=step_min, N_ticks=5, lb_max=0) for i in data]
    color_wheel = CABBI_colors.wheel()

    x0 = len(configuration_names) - 0.5
    xf = len(columns) - 0.5
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        bst.plots.plot_vertical_line(x0)
        ax.axvspan(x0, xf, color=colors.purple_tint.tint(60).RGBn)
        plt.xlim(-0.5, xf)

    for j in range(N_cols):
        color_wheel.restart()
        for i in range(N_rows):
            ax = axes[i]
            plt.sca(ax)
            plot(data[i, j], j)
            plt.ylabel(ylabels[i])
    
    for i in range(N_rows):
        ax = axes[i]
        plt.sca(ax)
        yticks = tickmarks[i]
        plt.ylim([yticks[0], yticks[1]])
        if yticks[0] < 0.:
            bst.plots.plot_horizontal_line(0, color=CABBI_colors.black.RGBn, linestyle='--')
        bst.plots.style_axis(ax,  
            xticks = xticks,
            yticks = yticks,
            xticklabels= xtext, 
            ytickf=False,
        )
    
    fig.align_ylabels(axes)
    plt.subplots_adjust(hspace=0)
    plt.sca(axes[1])
    # legend = plt.legend(
    #     handles=[
    #         mpatches.Patch(facecolor=color_wheel[0].RGBn, 
    #                        edgecolor=CABBI_colors.black.RGBn,
    #                        label='Lipid-cane only'),
    #         mpatches.Patch(facecolor=color_wheel[1].RGBn, 
    #                        edgecolor=CABBI_colors.black.RGBn,
    #                        label='Lipid-cane & lipid-sorghum'),
    #     ], 
    #     bbox_to_anchor=(0, 1, 1, 0), 
    #     loc="lower right", 
    #     # mode="expand", 
    #     # ncol=2
    # )
    # legend.get_frame().set_linewidth(0.0)

def plot_spearman_MFPP(top=None):
    MFPPs = []
    configuration_names = ['L1', 'L1*', 'L2', 'L2*']
    for name in configuration_names:
        file = spearman_file(name)
        try: 
            df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
        except: 
            warning = RuntimeWarning(f"file '{file}' not found")
            warn(warning)
            continue
        MFPPs.append(df['Biorefinery', 'MFPP [USD/ton]'])
    stream_price = format_units('USD/gal')
    ng_price = format_units('USD/cf')
    electricity_price = format_units('USD/kWhr')
    operating_days = format_units('day/yr')
    capacity = format_units('ton/hr')
    index = [
         'Cane lipid content [5 $-$ 15 dry wt. %]',
         'Relative sorghum lipid content [-3 $-$ 0 dry wt. %]',
         'Bagasse lipid retention [65 $-$ 75 %]',
         '$^a$Lipid extraction efficiency [baseline + 0 $-$ 25 %]',
        f'Plant capacity [330 $-$ 404 {capacity}]',
        f'Ethanol price [1.02, 1.80, 2.87 {stream_price}]',
        f'Biodiesel price relative to ethanol [0.31, 2.98, 4.11 {stream_price}]',
        f'Natural gas price [3.71, 4.73, 6.18 {ng_price}]',
        f'Electricity price [0.0583, 0.065, 0.069 {electricity_price}]',
        f'Operating days [180 $-$ 210 {operating_days}]',
         'IRR [10 $-$ 15 %]',
       # '$^a$Fermentation solids loading [20% $-$ 25%]',
    ]
    color_wheel = CABBI_colors.wheel()
    fig, ax = bst.plots.plot_spearman_2d(MFPPs, top=top, index=index, 
                                         color_wheel=color_wheel,
                                         name='MFPP')
    plt.legend(
        handles=[
            mpatches.Patch(
                color=color_wheel[i].RGBn, 
                label=format_name(configuration_names[i])
            )
            for i in range(len(configuration_names))
        ], 
        loc='upper left'
    )
    return fig, ax

def plot_configuration_breakdown(name, across_coordinate=False, **kwargs):
    load(name)
    if across_coordinate:
        return bst.plots.plot_unit_groups_across_coordinate(
            set_feedstock_lipid_content,
            [5, 10, 15],
            'Feedstock lipid content [dry wt. %]',
            unit_groups,
            colors=[area_colors[i.name].RGBn for i in unit_groups],
            hatches=[area_hatches[i.name] for i in unit_groups],
            **kwargs,
        )
    else:
        return bst.plots.plot_unit_groups(
            unit_groups,
            colors=[area_colors[i.name].RGBn for i in unit_groups],
            hatches=[area_hatches[i.name] for i in unit_groups],
            fraction=True,
            **kwargs,
        )
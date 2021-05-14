# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from chaospy import distributions as shape
from warnings import warn
import pandas as pd
import biosteam as bst
import numpy as np
from .. import PY37
from . import (utils,
               _process_settings,
               _chemicals,
               _system,
               _tea,
               _lipid_extraction_specification,
)

__all__ = [*utils.__all__,
           *_process_settings.__all__,
           *_chemicals.__all__,
           *_system.__all__,
           *_tea.__all__,
           *_lipid_extraction_specification.__all__,
           'lipidcane_sys',
           'lipidcane_tea', 
           'flowsheet',
]

from .utils import *
from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *
from ._lipid_extraction_specification import *
from biorefineries.sugarcane import create_tea as create_conventional_ethanol_tea

_system_loaded = False
_chemicals_loaded = False

name_keys = {
    'divided 1 and 2g front end oil separation': 5,
    'divided 1 and 2g hydrolyzate oil separation': 4,
    'combined 1 and 2g post fermentation oil separation': 3,
    'divided 1 and 2g post fermentation oil separation': 2,
    '1g': 1,
    'divided 1 and 2g bagasse expression': 0,
}

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def load(name, agile=False, cache={}):
    if name in cache:
        globals().update(cache[name])
        return
    import biosteam as bst
    from biosteam import main_flowsheet as F, UnitGroup
    global lipidcane_sys, lipidcane_tea, specs, flowsheet, _system_loaded
    global lipid_extraction_specification, model, unit_groups
    global MFPP, TCI, installed_equipment_cost, productivity
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('lipidcane2g')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    dct = globals()
    u = flowsheet.unit
    s = flowsheet.stream
    operating_hours = 24 * 200
    if name in name_keys: name = name_keys[name]
    
    ## System
    
    area_names = None
    def rename_storage_units(storage):
        bst.rename_units([i for i in lipidcane_sys.units if bst.is_storage_unit(i)], storage)
    
    if name == 5:
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,    
        )
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Biod. prod.',
            'Conv. ferm.', 
            'Pretreatment',
            'Cofementation',
            'Ethanol sep.', 
            'Wastewater treatment', 
            'Boiler turbogenerator',
            'Utilities',
            'Storage'
        ]
        rename_storage_units()
    elif name == 4:
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Biod. prod.', 
            'Conv. ferm.', 
            'Pretreatment', 
            'Cofementation',
            'Ethanol sep.', 
            'Wastewater treatment', 
            'Boiler turbogenerator',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_hydrolyzate_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units()
    elif name == 3:
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'Oil ext.',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1200)
    elif name == 2:
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Conv. ferm.',
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'Oil ext.',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1300)
    elif name == 1:
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
            'Storage',
        ]
        rename_storage_units(900)
    elif name == 0:
        area_names = [
            'Feedstock handling', 
            'Juicing', 
            'Conv. ferm.',
            'Oil ext.',
            'Pretreatment',
            'Cofementation',
            'Wastewater treatment',
            'CH&P',
            'Ethanol sep.', 
            'Biod. prod.',
            'Utilities',
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_1_and_2g_bagasse_expression(
            operating_hours=operating_hours,
        )
        rename_storage_units(1300)
    
    else:
        raise NotImplementedError(name)
    unit_groups = UnitGroup.group_by_area(lipidcane_sys.units)
    if area_names:
        for i, j in zip(unit_groups, area_names): i.name = j
    
    ## Not currently used; agile for next study
    
    if name in (1, 2, 3, 4) if isinstance(name, int) else '2g' in name:
        if agile:
            lipidcane_tea = create_agile_tea(lipidcane_sys.units)
            lipidcane_sys.simulate()
            scenario_a = sugarcane_tea.create_scenario(lipidcane_sys)
            cornstover_sys = trim_to_cornstover_hot_water_cellulosic_ethanol(
                lipidcane_sys,
                operating_hours=24 * 100,
            )
            cornstover_sys.simulate()
            scenario_b = lipidcane_tea.create_scenario(cornstover_sys)
            sugarcane_tea.compile_scenarios([scenario_a, scenario_b])
            dct.update(flowsheet.to_dict())
            lipidcane_tea.IRR = 0.10
            s.ethanol.price = lipidcane_tea.solve_price(s.ethanol)
            return
    
    ## TEA
    lipidcane_tea = create_tea(lipidcane_sys)
    lipidcane_tea.operating_days = 200
    
    ## Simulation
    try: 
        lipidcane_sys.simulate()
    except Exception as e:
        raise e
    else:
        lipidcane_tea.IRR = 0.10
    finally:
        dct.update(flowsheet.to_dict())
    
    ## Specification for analysis
    
    for i in lipidcane_sys.units:
        if getattr(i, 'tag', None) == 'lipid extraction efficiency':
            isplit_a = i.isplit
            break
    
    for i in lipidcane_sys.units:
        if getattr(i, 'tag', None) == 'bagasse lipid retention':
            isplit_b = i.isplit
            break
    
    lipid_extraction_specification = LipidExtractionSpecification(
            lipidcane_sys, lipidcane, isplit_a, isplit_b,
        )
    
    ## Model
        
    system = lipidcane_sys
    model = bst.Model(system)
    parameter = model.parameter
    ethanol_price = ethanol.price * 2.98668849 # USD/gal
    biodiesel_price = biodiesel.price * 3.3036 # USD/gal
    
    def uniform(lb, ub, *args, **kwargs):
        return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), **kwargs)
    
    uniform(0.05, 0.10, lipid_extraction_specification.load_lipid_content,
            element=lipidcane, units='dry wt. %')
    if name in (0, 1):
        lb = 0.45
        ub = 0.55
    elif name in (2, 3):
        lb = 0.75
        ub = 0.85
    else:
        lb = 0.2
        ub = 1.0
    uniform(lb, ub, lipid_extraction_specification.load_efficiency,
            element=lipidcane, units='%')
    uniform(0.6, 0.99, lipid_extraction_specification.load_lipid_retention,
            element=lipidcane, units='%')

    @uniform(0.9 * ethanol_price, 1.1 * ethanol_price, element=ethanol, units='USD/gal')
    def set_ethanol_price(price):
        ethanol.price = price / 2.98668849

    @uniform(0.9 * biodiesel_price, 1.1 * biodiesel_price, element=biodiesel, units='USD/gal')
    def set_biodiesel_price(price):
        biodiesel.price = price / 3.3111

    @uniform(0.0572, 0.07014, units='USD/kWh')
    def set_electricity_price(electricity_price):
        bst.PowerUtility.price = electricity_price
        
    @uniform(6 * 30, 7 * 30, units='day/yr')
    def set_operating_days(operating_days):
        lipidcane_tea.operating_days = operating_days
    
    @uniform(10., 15., units='%')
    def set_IRR(IRR):
        lipidcane_tea.IRR = IRR / 100.

    metric = model.metric
    kg_per_ton = 907.18474
    
    @metric(units='USD/ton')
    def MFPP():
        return kg_per_ton * lipidcane_tea.solve_price(lipidcane)

    @metric(units='10^6*USD')
    def TCI():
        return lipidcane_tea.TCI / 1e6 # 10^6*$

    @metric(units='10^6*USD')
    def installed_equipment_cost():
        return lipidcane_tea.installed_equipment_cost / 1e6 # 10^6*$
    
    @metric(units='GGE/ton')
    def productivity():
        feedstock = lipidcane.get_total_flow('ton/hr')
        GGE = (ethanol.F_mass * 2.98668849 / 1.5
           + biodiesel.get_total_flow('gal/hr') / 0.9536
           - sum([i.rate for i in lipidcane_sys.power_utilities]) * 3600 / feedstock / 131760
           - natural_gas.get_total_flow('ft3/hr') / 126.67)
        return GGE / feedstock
    
    cache[name] = dct.copy()

def evaluate_across_configurations(
        efficiency, lipid_content, lipid_retention, metrics, configurations,
    ):
    M = len(metrics)
    P = len(configurations)
    data = np.zeros([M, P])
    for i in range(P):
        load(int(configurations[i]))
        lipid_extraction_specification.load_specifications(
            efficiency=efficiency, 
            lipid_content=lipid_content, 
            lipid_retention=lipid_retention
        )
        lipid_extraction_specification.system.simulate()
        data[:, i] = [j() for j in metrics]
    return data

evaluate_across_configurations = np.vectorize(
    evaluate_across_configurations, 
    excluded=['metrics', 'configurations'],
    signature='(),(),(),(m),(p)->(m,p)'
)       

def run_uncertainty_and_sensitivity(name, N, rule='L'):
    import os
    folder = os.path.dirname(__file__)
    file = os.path.join(folder, f'lipidcane_monte_carlo_{name}.xlsx')
    load(name)
    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()
    model.table.to_excel(file)
    rho, p = model.spearman_r()
    file = os.path.join(folder, f'lipidcane_spearman_{name}.xlsx')
    rho.to_excel(file)

def plot_spearman_MFPP():
    import os
    import matplotlib.pyplot as plt
    from thermosteam.units_of_measure import format_units
    import matplotlib.patches as mpatches
    from biosteam.utils import CABBI_wheel
    folder = os.path.dirname(__file__)
    MFPPs = []
    for name in (1, 2, 3):
        file = os.path.join(folder, f'lipidcane_spearman_{name}.xlsx')
        try: 
            df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
        except: 
            warning = RuntimeWarning(f"file '{file}' not found")
            warn(warning)
            continue
        MFPPs.append(df['Biorefinery', 'MFPP [USD/ton]'])
    stream_price = format_units('USD/gal')
    electricity_price = format_units('USD/kWhr')
    operating_days = format_units('day/yr')
    index = [
        f'Feedstock lipid content [5 $-$ 10 dry wt. %]',
        f'*Lipid extraction efficiency [45 $-$ 55; 75 $-$ 80 %]',
        f'Bagasse lipid retention [60 $-$ 99 %]',
        f'Ethanol price [2.12 $-$ 2.59 {stream_price}]',
        f'Biodiesel price [4.10 $-$ 5.01 {stream_price}]',
        f'Electricity price [0.0572 $-$ 0.07 {electricity_price}]',
        f'Operating days [180 $-$ 210 {operating_days}]',
         'IRR [10 $-$ 15 %]',
    ]
    CABBI_wheel = tuple(CABBI_wheel)
    fig, ax = bst.plots.plot_spearman_2d(MFPPs, index=index, color_wheel=CABBI_wheel,
                                         name='MFPP')
    plt.legend(
        handles=[
            mpatches.Patch(color=CABBI_wheel[0].RGBn, label='Configuration I'),
            mpatches.Patch(color=CABBI_wheel[1].RGBn, label='Configuration II'),
            mpatches.Patch(color=CABBI_wheel[2].RGBn, label='Configuration III'),
        ], 
        # loc='upper left'
    )
    return fig, ax

def ethanol_price():
    return lipidcane_tea.solve_price(flowsheet.stream.ethanol) * 2.98668849
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biosteam import main_flowsheet, UnitGroup
from chaospy import distributions as shape
from warnings import warn
import pandas as pd
import numpy as np
from ..lipidcane import set_lipid_fraction, get_lipid_fraction
from .. import PY37
from . import (_process_settings,
               _chemicals,
               _system,
               _tea,
               _lipid_extraction_specification,
               _distributions,
)

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

from ._process_settings import *
from ._chemicals import *
from ._system import *
from ._tea import *
from ._lipid_extraction_specification import *
from ._distributions import *
from biorefineries.sugarcane import create_tea as create_conventional_ethanol_tea

_system_loaded = False
_chemicals_loaded = False

name_keys = {
    'combined 1 and 2g post fermentation oil separation': 2,
    '1g': 1,
}

def load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def load(name, agile=False, cache={}):
    if (name, agile) in cache:
        globals().update(cache[name, agile])
        return
    global lipidcane_sys, lipidcane_tea, specs, flowsheet, _system_loaded
    global lipid_extraction_specification, model, unit_groups
    global MFPP, TCI, installed_equipment_cost, productivity
    if not _chemicals_loaded: load_chemicals()
    flowsheet = bst.Flowsheet('lipidcane2g')
    main_flowsheet.set_flowsheet(flowsheet)
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
    
    if name == 1:
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
    elif name == 2:
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
            'Storage'
        ]
        lipidcane_sys = create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(
            operating_hours=operating_hours,
        )
        rename_storage_units(1000)
    else:
        raise NotImplementedError(name)
    unit_groups = UnitGroup.group_by_area(lipidcane_sys.units)
    if area_names:
        for i, j in zip(unit_groups, area_names): i.name = j
    
    for BT in lipidcane_sys.units:
        if isinstance(BT, bst.BoilerTurbogenerator): break
    
    ## TEA
    lipidcane_sys.set_tolerance(rmol=1e-3)
    units = lipidcane_sys.units
    dct.update(flowsheet.to_dict())
    if agile: 
        dct['lipidsorghum'] = lipidsorghum = lipidsorghum_feedstock(lipidcane)
        original_lipidcane = lipidcane.copy()
        
        class AgileLipidcaneSystem(bst.AgileSystem):
            
            def set_parameters(self, switch2sorghum):
                F_mass_original = lipidcane.F_mass
                lipid_content = get_lipid_fraction(lipidcane)
                if switch2sorghum:
                    lipidcane.copy_like(lipidsorghum)
                else:
                    lipidcane.copy_like(original_lipidcane)
                set_lipid_fraction(lipid_content, lipidcane)
                lipidcane.F_mass = F_mass_original
                np.testing.assert_allclose(get_lipid_fraction(lipidcane), lipid_content)
                
        lipidcane_tea = create_agile_tea(lipidcane_sys.units)        
        lipidcane_sys = AgileLipidcaneSystem(lipidcane_sys, [0, 1], [200 * 24, 60 * 24],
                                             tea=lipidcane_tea)
    else:
        lipidcane_tea = create_tea(lipidcane_sys)
        lipidcane_sys.operating_hours = 24 * 200
    feedstocks = [lipidcane]
    
    ## Specification for analysis
    
    for i in units:
        if getattr(i, 'tag', None) == 'lipid extraction efficiency':
            isplit_a = i.isplit
            break
    
    for i in units:
        if getattr(i, 'tag', None) == 'bagasse lipid retention':
            isplit_b = i.isplit
            break
    
    lipid_extraction_specification = LipidExtractionSpecification(
            lipidcane_sys, feedstocks, isplit_a, isplit_b,
        )
    
    ## Model
    
    kg_per_ton = 907.18474
    model = bst.Model(lipidcane_sys, exception_hook='raise')
    parameter = model.parameter
    metric = model.metric
    
    def uniform(lb, ub, *args, **kwargs):
        return parameter(*args, distribution=shape.Uniform(lb, ub), bounds=(lb, ub), **kwargs)
    
    capacity = lipidcane.F_mass / kg_per_ton
    @uniform(0.5 * capacity, 2. * capacity, units='ton/hr', element=lipidcane)
    def set_plant_capacity(capacity):
        lipidcane.F_mass = kg_per_ton * capacity
    
    # Currently at ~5%, but total lipid content is past 10%
    
    set_feedstock_lipid_content = uniform(0.05, 0.10, lipid_extraction_specification.load_lipid_content,
            element=lipidcane, units='dry wt. %')
    if name == 1:
        lb = 0.475
        ub = 0.525
    elif name == 2:
        lb = 0.75
        ub = 0.80
    set_lipid_extraction_efficiency = uniform(lb, ub, lipid_extraction_specification.load_efficiency,
            element=lipidcane, units='%')
    set_bagasse_lipid_retention = uniform(0.80, 0.90, lipid_extraction_specification.load_lipid_retention,
            element=lipidcane, units='%')

    # USDA ERS historical price data
    @parameter(distribution=ethanol_price_distribution, element=ethanol, units='USD/gal')
    def set_ethanol_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to Nov 2020
        ethanol.price = price / 2.98668849

    # USDA ERS historical price data
    @parameter(distribution=biodiesel_price_distribution, element=biodiesel, units='USD/gal')
    def set_biodiesel_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        biodiesel.price = price / 3.3111

    # https://www.eia.gov/energyexplained/natural-gas/prices.php
    @parameter(distribution=natural_gas_price_distribution, element=natural_gas, units='USD/cf')
    def set_natural_gas_price(price): # Triangular distribution fitted over the past 10 years Sep 2009 to March 2021
        BT.natural_gas_price = 51.92624700383502 * price / 1000. 

    # From Emma's literature search; EIA electricity prices vary much too widely across states.
    # Best use prices similar to previous TEAs for a more accurate comparison.
    @uniform(0.0572, 0.07014, units='USD/kWh')
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

    # Initial solids loading in Humbird is 20% (includes both dissolved and insoluble solids)
    @uniform(20, 25, units='%', element='Fermentation')
    def set_fermentation_solids_loading(solids_loading):
        if name == 2: M402.solids_loading = solids_loading / 100.
    
    natural_gas.phase = 'g'
    natural_gas.set_property('T', 60, 'degF')
    natural_gas.set_property('P', 14.73, 'psi')
    original_value = natural_gas.imol['CH4']
    natural_gas.imass['CH4'] = 1 
    V_ng = natural_gas.get_total_flow('ft3/hr')
    natural_gas.imol['CH4'] = original_value
    
    if agile:
        feedstock_flow = lambda: lipidcane_sys.flow_rates[lipidcane] / kg_per_ton # ton/yr
        biodiesel_flow = lambda: lipidcane_sys.flow_rates[biodiesel] / 3.3111 # gal/yr
        ethanol_flow = lambda: lipidcane_sys.flow_rates[ethanol] / 2.98668849 # gal/yr
        natural_gas_flow = lambda: lipidcane_sys.flow_rates[natural_gas] * V_ng # cf/yr
        if name == 1:
            electricity = lambda: lipidcane_sys.utility_cost / bst.PowerUtility.price
        elif name == 2:
            electricity = lambda: 0.
    else:
        feedstock_flow = lambda: lipidcane_sys.operating_hours * lipidcane.F_mass / kg_per_ton
        biodiesel_flow = lambda: lipidcane_sys.operating_hours * biodiesel.F_mass / 3.3111 # gal/yr
        ethanol_flow = lambda: lipidcane_sys.operating_hours * ethanol.F_mass / 2.98668849 # gal/yr
        natural_gas_flow = lambda: lipidcane_sys.operating_hours * natural_gas.F_mass * V_ng # cf/yr
        if name == 1:
            electricity = lambda: lipidcane_sys.operating_hours * sum([i.rate for i in lipidcane_sys.power_utilities])
        elif name == 2:
            electricity = lambda: 0.
    
    @metric(units='USD/ton')
    def MFPP():
        return kg_per_ton * lipidcane_tea.solve_price(lipidcane)
    
    @metric(units='gal/ton')
    def biodiesel_production():
        return biodiesel_flow() / feedstock_flow()
    
    @metric(units='gal/ton')
    def ethanol_production():
        return ethanol_flow() / feedstock_flow()
    
    @metric(units='kWhr/ton')
    def electricity_production():
        return - electricity() / feedstock_flow()
    
    @metric(units='cf/ton')
    def natural_gas_consumption():
        return natural_gas_flow() / feedstock_flow()
    
    @metric(units='GGE/ton')
    def productivity():
        GGE = (ethanol_flow() / 1.5
           + biodiesel_flow() / 0.9536
           - electricity() * 3600 / 131760
           - natural_gas_flow() / 126.67)
        return GGE / feedstock_flow()

    @metric(units='10^6*USD')
    def TCI():
        return lipidcane_tea.TCI / 1e6 # 10^6*$
    
    # # Single point evaluation for detailed design results
    set_feedstock_lipid_content(0.10) # Consistent with SI of Huang's 2016 paper
    if name == 1: 
        efficiency = 0.5
    elif name == 2: 
        efficiency = 0.775
    set_lipid_extraction_efficiency(efficiency) # Middle point
    set_bagasse_lipid_retention(0.9)  # Middle point
    set_fermentation_solids_loading(20) # Same as Humbird
    set_ethanol_price(2.356) # Consistent with Huang's 2016 paper
    set_biodiesel_price(4.569) # Consistent with Huang's 2016 paper
    set_natural_gas_price(4.198) # Consistent with Humbird's 2012 paper
    set_electricity_price(0.0572) # Consistent with Humbird's 2012 paper
    set_operating_days(200) # Consistent with Huang's 2016 paper
    
    for i in model.parameters + model.metrics:
        dct[i.name] = i
    cache[name, agile] = dct.copy()
    
    ## Simulation
    try: 
        lipidcane_sys.simulate()
    except Exception as e:
        raise e
    finally:
        lipidcane_tea.IRR = 0.10

def evaluate_across_configurations(
        efficiency, lipid_content, lipid_retention, agile, configurations, metrics,
    ):
    A = len(agile)
    C = len(configurations)
    M = len(metrics)
    data = np.zeros([A, C, M])
    for ia in range(A):
        for ic in range(C):    
            load(int(configurations[ic]), agile[ia])
            lipid_extraction_specification.load_specifications(
                efficiency=efficiency, 
                lipid_content=lipid_content, 
                lipid_retention=lipid_retention
            )
            lipid_extraction_specification.system.simulate()
            data[ia, ic, :] = [j() for j in metrics]
    return data

evaluate_across_configurations = np.vectorize(
    evaluate_across_configurations, 
    excluded=['agile', 'configurations', 'metrics'],
    signature='(),(),(),(a),(c),(m)->(a,c,m)'
)       

def spearman_file(name, agile):
    import os
    folder = os.path.dirname(__file__)
    filename = f'lipidcane_spearman_{name}'
    if agile: filename += '_agile'
    filename += '.xlsx'
    return os.path.join(folder, filename)

def monte_carlo_file(name, agile):
    import os
    folder = os.path.dirname(__file__)
    filename = f'lipidcane_monte_carlo_{name}'
    if agile: filename += '_agile'
    filename += '.xlsx'
    return os.path.join(folder, filename)

def run_uncertainty_and_sensitivity(name, agile, N, rule='L'):
    import os
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.utils.DesignWarning)
    load(name, agile)
    samples = model.sample(N, rule)
    model.load_samples(samples)
    model.evaluate()
    file = monte_carlo_file(name, agile)
    model.table.to_excel(file)
    rho, p = model.spearman_r()
    file = spearman_file(name, agile)
    rho.to_excel(file)
    
def run_all(N, rule='L'):
    for name, agile in ((1, False), (2, False), (1, True), (2, True)):
        run_uncertainty_and_sensitivity(name, agile, N, rule)

def plot_spearman_MFPP():
    import matplotlib.pyplot as plt
    from thermosteam.units_of_measure import format_units
    import matplotlib.patches as mpatches
    from biosteam.utils import CABBI_wheel
    MFPPs = []
    for name, agile in ((1, False), (2, False), (1, True), (2, True)):
        file = spearman_file(name, agile)
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
        f'Plant capacity [330 $-$ 404 {capacity}]',
          'Feedstock lipid content [5 $-$ 10 dry wt. %]',
          '$^b$Lipid extraction efficiency [47.5 $-$ 52.5; 75 $-$ 80 %]',
          'Bagasse lipid retention [80 $-$ 90 %]',
        f'Ethanol price [0.89, 2.17, 3.18 {stream_price}]',
        f'Biodiesel price [2.55, 3.62, 5.74 {stream_price}]',
        f'Natural gas price [3.71, 4.73, 6.18 {ng_price}]',
        f'Electricity price [0.0572 $-$ 0.07 {electricity_price}]',
        f'Operating days [180 $-$ 210 {operating_days}]',
          'IRR [10 $-$ 15 %]',
          '$^a$Fermentation solids loading [20% $-$ 25%]',
    ]
    color_wheel = tuple(CABBI_wheel)
    fig, ax = bst.plots.plot_spearman_2d(MFPPs, index=index, 
                                         color_wheel=color_wheel,
                                         name='MFPP')
    plt.legend(
        handles=[
            mpatches.Patch(color=color_wheel[0].RGBn, label='Configuration I'),
            mpatches.Patch(color=color_wheel[1].RGBn, label='Configuration II'),
            mpatches.Patch(color=color_wheel[2].RGBn, label='Agile Configuration I'),
            mpatches.Patch(color=color_wheel[3].RGBn, label='Agile Configuration II'),
        ], 
        # loc='upper left'
    )
    return fig, ax

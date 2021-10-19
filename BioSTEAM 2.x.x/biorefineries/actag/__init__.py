# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
import numpy as np
import os
from biorefineries.cornstover import create_tea
from ._chemicals import create_conventional_chemicals, create_cellulosic_chemicals
from ._system import (
    create_cellulosic_acTAG_system,
    create_conventional_acTAG_system,
)
from thermosteam.units_of_measure import format_units
from ._units import *
from ._contours import *
from biosteam.plots import (
    plot_contour_2d,
    MetricBar,
    plot_scatter_points,
    plot_contour_single_metric,
)

# TODO: find composition for minimum product specification (24 mm2/s viscosity)

MPSP, TCI, heating_duty, cooling_duty = all_metric_mockups = (
    bst.MockVariable('MPSP', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('TCI', '10^6*USD', 'Biorefinery'),
    bst.MockVariable('Heating duty', 'GJ/yr', 'Biorefinery'),
    bst.MockVariable('Cooling duty', 'GJ/yr', 'Biorefinery'),
)
N_metrics = len(all_metric_mockups)
kg_per_ton = 907.18474
_conventional_chemicals_loaded = False
_cellulosic_chemicals_loaded = False

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 607.5 # 2019
    bst.PowerUtility.price = 0.065
    HeatUtility = bst.HeatUtility
    steam_utility = HeatUtility.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 0.9
    steam_utility.regeneration_price = 0.30626
    steam_utility.T = 529.2
    steam_utility.P = 44e5
    HeatUtility.get_agent('cooling_water').regeneration_price = 0
    HeatUtility.get_agent('chilled_water').heat_transfer_price = 0

def load_conventional_chemicals():
    global conventional_chemicals, _chemicals_loaded
    conventional_chemicals = create_conventional_chemicals()
    _conventional_chemicals_loaded = True

def load_cellulosic_chemicals():
    global cellulosic_chemicals, _chemicals_loaded
    cellulosic_chemicals = create_cellulosic_chemicals()
    _cellulosic_chemicals_loaded = True

def load(configuration, cache={}):
    global tea, sys, sys_no_dry_fractionation, fermentation, model, set_selectivity
    configuration = int(configuration)
    dct = globals()
    if configuration in cache:
        dct.update(cache[configuration])
        return
    if configuration == 1:
        if not _conventional_chemicals_loaded: 
            load_conventional_chemicals()
            chemicals = conventional_chemicals
        flowsheet = bst.Flowsheet('conventional_acTAG')
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(chemicals)
        sys = create_conventional_acTAG_system(f'actag_sys_{configuration}')
        tea = create_tea(sys)
        operating_hours = 200 * 24
    elif configuration == 2:
        if not _cellulosic_chemicals_loaded: 
            load_cellulosic_chemicals()
            chemicals = cellulosic_chemicals
        flowsheet = bst.Flowsheet('cellulosic_acTAG')
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(chemicals)
        sys = create_cellulosic_acTAG_system(f'actag_sys_{configuration}')
        operating_hours = 330 * 24
    else:
        raise ValueError(f"invalid configuration '{configuration}'; only 1 and 2 are valid")
    u = flowsheet.unit
    s = flowsheet.stream
    dct.update(flowsheet.to_dict())
    lipid = u.F401.outs[1]
    downstream_units = lipid.sink.get_downstream_units()
    downstream_units.add(lipid.sink)
    u.F401.outs[1] = s.acTAG
    sys_no_dry_fractionation = bst.System.from_units(
        f'actag_sys_{configuration}_no_dry_fractionation', 
        [i for i in sys.units if i not in downstream_units], 
    )
    sys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    sys_no_dry_fractionation.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    tea = create_tea(sys)
    tea_no_dry_fractionation = create_tea(sys_no_dry_fractionation)
    tea.operating_hours = tea_no_dry_fractionation.operating_hours = operating_hours
    load_process_settings()
    fermentation = u.R301
    fermentation.product_yield = 0.06
    fermentation.selectivity = 0.5
    ## Model
    model = bst.Model(sys, exception_hook='raise', retry_evaluation=False)
    parameter = model.parameter
    metric = model.metric
    
    @parameter
    def set_selectivity(selectivity):
        fermentation.selectivity = selectivity / 100.
    
    @metric(units='USD/ton')
    def MPSP(): 
        if fermentation.selectivity == 1.:
            return tea_no_dry_fractionation.solve_price(s.acTAG) * kg_per_ton
        else:
            return tea.solve_price(s.acTAG) * kg_per_ton
    
    @metric(units='10^6*USD')
    def TCI(): 
        if fermentation.selectivity == 1.:
            return tea_no_dry_fractionation.TCI / 1e6 # 10^6*$
        else:
            return tea.TCI / 1e6 # 10^6*$
    
    @metric(units='GJ/yr')
    def heating_duty(): 
        if fermentation.selectivity == 1.:
            return sys_no_dry_fractionation.get_heating_duty() / 1e6
        else:
            return sys.get_heating_duty() / 1e6
    
    @metric(units='GJ/yr')
    def cooling_duty(): 
        if fermentation.selectivity == 1.:
            return sys_no_dry_fractionation.get_cooling_duty() / 1e6
        else:
            return sys.get_cooling_duty() / 1e6
    
    for i in model._parameters:
        dct[i.setter.__name__] = i
    for i in model._metrics:
        dct[i.getter.__name__] = i
    cache[configuration] = dct.copy()
    
def evaluate_across_yield_titer_selectivity_and_productivity(product_yield, titer, selectivity, productivities, configuration):
    load(configuration)
    metrics = model.metrics
    P = len(productivities)
    M = len(model.metrics)
    data = np.zeros([P, M])
    fermentation.titer = titer
    fermentation.product_yield = product_yield / 100.
    fermentation.selectivity = selectivity / 100.
    try:
        if selectivity == 100.:
            sys_no_dry_fractionation.empty_recycles()
            sys_no_dry_fractionation.simulate()
        else:
            sys.empty_recycles()
            sys.simulate()
    except RuntimeError as e:
        if str(e) == 'infeasible to evaporate any more water':
            for i in range(P):
                for j in range(M):
                    data[i, j] = np.nan
        else:
            raise e
        
    for i in range(P):
        fermentation.tau = titer / productivities[i]
        fermentation._summary()
        for j in range(M):
            data[i, j] = metrics[j]()
    return data

evaluate_across_yield_titer_selectivity_and_productivity = np.vectorize(
    evaluate_across_yield_titer_selectivity_and_productivity, 
    excluded=['productivities', 'configuration'],
    signature=f'(),(),(),(p),()->(p,{N_metrics})'
)

def fermentation_data(configuration, load):
    # Generate contour data
    x = np.linspace(40, 90, 20)
    y = np.linspace(5, 80, 20)
    z = np.array([50, 75, 100])
    w = np.array([0.033, 1.00])
    X, Y, Z = np.meshgrid(x, y, z)
    folder = os.path.dirname(__file__)
    file = f'fermentation_data_{configuration}.npy'
    file = os.path.join(folder, file)
    if load:
        data = np.load(file)
    else:
        data = evaluate_across_yield_titer_selectivity_and_productivity(
            X, Y, Z, w, configuration
        )
    # M, N = data.shape[-2:]
    # for i in range(1, M):
    #     for j in range(1, N):
    #         if 
    np.save(file, data)
    return X, Y, z, w, data
    

    
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

MPSP, TCI = all_metric_mockups = (
    bst.MockVariable('MPSP', 'USD/ton', 'Biorefinery'),
    bst.MockVariable('TCI', '10^6*USD', 'Biorefinery'),
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
    global tea, sys, fermentation, model, set_selectivity
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
        sys = create_conventional_acTAG_system()
        tea = create_tea(sys)
        tea.operating_hours = 200 * 24
    elif configuration == 2:
        if not _cellulosic_chemicals_loaded: 
            load_cellulosic_chemicals()
            chemicals = cellulosic_chemicals
        flowsheet = bst.Flowsheet('cellulosic_acTAG')
        bst.main_flowsheet.set_flowsheet(flowsheet)
        bst.settings.set_thermo(chemicals)
        sys = create_cellulosic_acTAG_system()
        tea = create_tea(sys)
        tea.operating_hours = 330 * 24
    else:
        raise ValueError(f"invalid configuration '{configuration}'; only 1 and 2 are valid")
    u = flowsheet.unit
    s = flowsheet.stream
    sys.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
    dct.update(flowsheet.to_dict())
    load_process_settings()
    fermentation = u.R301
    
    ## Model
    model = bst.Model(sys, exception_hook='raise', retry_evaluation=False)
    parameter = model.parameter
    metric = model.metric
    
    @parameter
    def set_selectivity(selectivity):
        if configuration == 1:
            X = fermentation.fermentation_reaction.X
        elif configuration == 2:
            X = fermentation.cofermentation.X
        X[0] = fermentation.product_yield * selectivity
        X[1] = fermentation.product_yield * (1 - selectivity)
    
    @metric(units='USD/ton')
    def MPSP(): return tea.solve_price(s.acTAG) * kg_per_ton
    
    @metric(units='10^6*USD')
    def TCI(): return tea.TCI / 1e6 # 10^6*$
    
    for i in model._parameters:
        dct[i.setter.__name__] = i
    for i in model._metrics:
        dct[i.getter.__name__] = i
    cache[configuration] = dct.copy()
    
def evaluate_across_titer_yield_selectivity_and_productivity(titer, product_yield, selectivity, productivities, configuration):
    load(configuration)
    metrics = model.metrics
    P = len(productivities)
    M = N_metrics
    data = np.zeros([P, M])
    fermentation.titer = titer
    fermentation.product_yield = product_yield / 100.
    set_selectivity(selectivity / 100.)
    sys.simulate()
    for i in range(P):
        fermentation.tau = titer / productivities[i]
        fermentation._summary()
        for j in range(M):
            data[i, j] = metrics[j]()
    return data

evaluate_across_titer_yield_selectivity_and_productivity = np.vectorize(
    evaluate_across_titer_yield_selectivity_and_productivity, 
    excluded=['productivities', 'configuration'],
    signature=f'(),(),(),(p),()->(p,{N_metrics})'
)

def fermentation_data(configuration, load):
    # Generate contour data
    x = np.linspace(5, 40, 15)
    y = np.linspace(40, 90, 15)
    z = np.array([50, 70, 90])
    w = np.array([0.01, 0.10, 1.00])
    X, Y, Z = np.meshgrid(x, y, z)
    folder = os.path.dirname(__file__)
    file = f'fermentation_data_{configuration}.npy'
    file = os.path.join(folder, file)
    if load:
        data = np.load(file)
    else:
        data = evaluate_across_titer_yield_selectivity_and_productivity(
            X, Y, Z, w, configuration
        )
    np.save(file, data)
    return X, Y, z, w, data
    
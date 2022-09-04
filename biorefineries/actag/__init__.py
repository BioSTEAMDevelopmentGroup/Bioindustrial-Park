# -*- coding: utf-8 -*-
"""
"""
from biorefineries import oilcane as oc
import biosteam as bst
import numpy as np
import pandas as pd
import os
from biorefineries.cornstover import create_tea
from ._chemicals import *
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

# %% Baseline and targets
# L_per_gal = 3.7854
# biodiesel_L_per_kg = 1.143245447132373
# mean_biodiesel_price_5yr = np.mean([3.13, 3.23, 3.31, 3.33, 3.3, 3.34, 
#                                     3.23, 3.45, 3.16, 3.04, 3.02, 3.13, 
#                                     3.11, 3.02, 3.01, 3.04, 3.09, 3.09,
#                                     3.06, 3.12, 3.24, 3.08, 2.98, 2.86, 
#                                     2.87, 2.91, 2.99, 3.06, 3.13, 3.21, 
#                                     3.26, 3.35, 3.11, 2.97, 2.77, 2.74, 
#                                     2.87, 2.92, 3.12, 3.3, 2.88, 3.51,
#                                     3.61, 4.08, 4.91, 5, 5.99, 5.91, 5.71,
#                                     5.49, 4.88, 4.83, 6.03, 6.5, 6.99, 7.53, 
#                                     7.35]) /  L_per_gal * biodiesel_L_per_kg

baseline = {
    1: {'Yield': 34.,
        'Titer': 12.,
        'Selectivity': 50,
        'Productivity': 0.0625},
    2: {'Yield': 34., # TODO: Update with actual experimental values
        'Titer': 12.,
        'Selectivity': 50.,
        'Productivity': 0.0625}
}
target = {
    1: {'Yield': 55.,
        'Titer': 30.,
        'Selectivity': 75.,
        'Productivity': 0.7},
    2: {'Yield': 55.,
        'Titer': 30.,
        'Selectivity': 75.,
        'Productivity': 0.7},
}
for i in baseline, target:
    i[3] = i[1]
    i[4] = i[2]

# %% Mock variables and settings

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

# %% System loading

def load_process_settings():
    bst.process_tools.default_utilities()
    bst.CE = 607.5 # 2019
    bst.PowerUtility.price = 0.065
    HeatUtility = bst.HeatUtility
    steam_utility = HeatUtility.get_agent('low_pressure_steam')
    steam_utility.heat_transfer_efficiency = 1.0
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

def load(configuration, simulate=None, cache={}):
    global tea, sys, sys_no_dry_fractionation, model, set_selectivity
    configuration = int(configuration)
    key = (configuration, simulate)
    dct = globals()
    if key in cache:
        dct.update(cache[key])
        if configuration in (3, 4):
            flowsheet = oc.flowsheet
            fermentation = flowsheet('Fermentation')
    else:
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
        elif configuration == 3:
            oc.load('O9', cache={})
            sys = oc.sys
            sys.N_runs = 2
            flowsheet = oc.flowsheet
            fermentation = flowsheet('Fermentation')
            # oc.biodiesel.price = mean_biodiesel_price_5yr
        elif configuration == 4:
            oc.load('O10', cache={})
            sys = oc.sys
            sys.N_runs = 2
            flowsheet = oc.flowsheet
            fermentation = flowsheet('CoFermentation')
            # oc.biodiesel.price = mean_biodiesel_price_5yr
        else:
            raise ValueError(f"invalid configuration '{configuration}'; only 1 and 2 are valid")
        u = flowsheet.unit
        s = flowsheet.stream
        dct.update(flowsheet.to_dict())
        if configuration in (1, 2):
            lipid = u.F401.outs[1]
            downstream_units = lipid.sink.get_downstream_units()
            downstream_units.add(lipid.sink)
            u.F401.outs[1] = s.acTAG
            sys_no_dry_fractionation = bst.System.from_units(
                f'actag_sys_{configuration}_no_dry_fractionation', 
                [i for i in sys.units if i not in downstream_units], 
            )
            sys_no_dry_fractionation.set_tolerance(rmol=1e-5, mol=1e-3, subsystems=True)
            tea = create_tea(sys)
            tea_no_dry_fractionation = create_tea(sys_no_dry_fractionation)
            tea.operating_hours = tea_no_dry_fractionation.operating_hours = operating_hours
            tea.income_tax = 0.21
            tea_no_dry_fractionation.income_tax = 0.21
            load_process_settings()
            fermentation = u.R301
        else:
            tea = oc.tea
        sys.set_tolerance(rmol=1e-6, mol=1e-3, subsystems=True)
        dct['fermentation'] = fermentation 
        fermentation.product_yield = 0.34
        fermentation.selectivity = 0.5
        ## Model
        model = bst.Model(sys, exception_hook='raise', retry_evaluation=False)
        parameter = model.parameter
        metric = model.metric
        
        @parameter
        def set_selectivity(selectivity):
            fermentation.selectivity = selectivity / 100.
        
        @metric(units='USD/MT')
        def MPSP(): 
            if fermentation.selectivity == 1. and configuration in (1, 2):
                return tea_no_dry_fractionation.solve_price(s.acTAG) * 1000
            else:
                return tea.solve_price(s.acTAG) * 1000
        
        @metric(units='10^6*USD')
        def TCI(): 
            if fermentation.selectivity == 1. and configuration in (1, 2):
                return tea_no_dry_fractionation.TCI / 1e6 # 10^6*$
            else:
                return tea.TCI / 1e6 # 10^6*$
        
        for i in model._parameters:
            dct[i.setter.__name__] = i
        for i in model._metrics:
            dct[i.getter.__name__] = i
        cache[key] = dct.copy()
        
    if simulate == 'baseline':
        dct = baseline[configuration]
    elif simulate == 'target':
        dct = target[configuration]
    elif simulate is None:
        return
    else:
        raise ValueError(f"simulate must be 'target' or 'baseline'; not {simulate}")
    
    fermentation.selectivity = selectivity = dct['Selectivity'] / 100.
    fermentation.productivity = dct['Productivity']
    fermentation.titer = dct['Titer']
    fermentation.product_yield = dct['Yield'] / 100.
    if selectivity == 1. and configuration in (1, 2):
        sys_no_dry_fractionation.simulate()
    else:
        sys.simulate()
    
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
        if selectivity == 100. and configuration in (1, 2):
            # sys_no_dry_fractionation.empty_recycles()
            sys_no_dry_fractionation.simulate()
        else:
            # sys.empty_recycles()
            sys.simulate()
    except RuntimeError as e:
        if str(e) == 'infeasible to evaporate any more water':
            for i in range(P):
                for j in range(M):
                    data[i, j] = np.nan
        else:
            raise e
        
    for i in range(P):
        fermentation.productivity = productivities[i]
        fermentation.tau = titer / productivities[i]
        fermentation._reevaluate()
        for j in range(M):
            data[i, j] = metrics[j]()
    return data

evaluate_across_yield_titer_selectivity_and_productivity = np.vectorize(
    evaluate_across_yield_titer_selectivity_and_productivity, 
    excluded=['productivities', 'configuration'],
    signature=f'(),(),(),(p),()->(p,{N_metrics})'
)
N = 15
def fermentation_data(configuration, load):
    # Generate contour data
    x = np.linspace(30, 90, N)
    y = np.linspace(5, 60, N)
    z = np.array([50, 75, 100])
    w = np.array([0.0625, 0.7])
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

def save_contour_plot_data(configuration):
    X, Y, selectivities, productivities, data = fermentation_data(configuration, load)
    yields = np.linspace(30, 90, N)
    titers = np.linspace(5, 60, N)
    MPSP_data = data[:, :, :, :, 0]
    folder = os.path.dirname(__file__)
    if configuration == 1:
        configuration_name = 'sugarcane'
    elif configuration == 2:
        configuration_name = 'switchgrass'
    for i, selectivity in enumerate(selectivities):
        for j, productivity in enumerate(productivities):
            file = f'{configuration_name}_biorefinery_selectivity_{selectivity}_productivity_{productivity}.xlsx'
            file = os.path.join(folder, file)
            pd.DataFrame(
                data=MPSP_data[:, :, i, j],
                index=pd.MultiIndex.from_tuples([(i,) for i in yields], names=['Yield [% theoretical]']),
                columns=pd.MultiIndex.from_tuples([(i,) for i in titers], names=['Titer [g / L]']),
            ).to_excel(file, startrow=1)
    

def save_tables():
    import biorefineries.actag as actag
    import biorefineries.cornstover as cs
    import biosteam as bst
    actag.load(3, 'baseline')
    tea1_baseline = actag.tea
    sys1_baseline = actag.sys
    actag.load(3, 'target')
    tea1_target = actag.tea
    sys1_target = actag.sys
    actag.load(4, 'baseline')
    tea2_baseline = actag.tea
    sys2_baseline = actag.sys
    actag.load(4, 'target')
    tea2_target = actag.tea
    sys2_target = actag.sys
    folder = os.path.dirname(__file__)
    file = 'tables.xlsx'
    file = os.path.join(folder, file)
    writer = pd.ExcelWriter(file)
    bst.report.voc_table(
        [sys1_baseline, sys2_baseline, sys1_target, sys2_target], 
        ['TAG', 'acTAG'],
        ['Baseline conventional', 'Baseline cellulosic', 
         'Target conventional', 'Target cellulosic']).to_excel(writer, 'VOC')
    cs.capex_table(
        [tea1_baseline, tea2_baseline, tea1_target, tea2_target], 
        ['Baseline conventional', 'Baseline cellulosic', 
         'Target conventional', 'Target cellulosic']).to_excel(writer, 'CAPEX')
    cs.foc_table(
        [tea1_baseline, tea2_baseline, tea1_target, tea2_target], 
        ['Baseline conventional', 'Baseline cellulosic', 
         'Target conventional', 'Target cellulosic']).to_excel(writer, 'FOC')
    writer.save()
    
def save_plots(load=True, single_bar=False):
    import matplotlib.pyplot as plt
    if single_bar:
        set_plot_style()
        plot_yield_titer_selectivity_productivity_contours([3, 4], load=load)
        folder = os.path.dirname(__file__)
        file = 'contour_plots.svg'
        file = os.path.join(folder, file)
        plt.savefig(file, transparent=True)
    else:
        set_plot_style(half=True)
        folder = os.path.dirname(__file__)
        fig0, axes0 = plot_yield_titer_selectivity_productivity_contours(3, load=load)
        plt.subplots_adjust(left=0.17, right=0.86)
        # plt.tight_layout()
        # fig0.set_constrained_layout_pads(hspace=0.02, wspace=0.02)
        file = 'contour_plots_left.svg'
        file = os.path.join(folder, file)
        plt.savefig(file, transparent=True)
        fig1, axes1 = plot_yield_titer_selectivity_productivity_contours(4, load=load)
        plt.subplots_adjust(left=0.17, right=0.86)
        # plt.tight_layout()
        # fig1.set_constrained_layout_pads(hspace=0.02, wspace=0.02)
        file = 'contour_plots_right.svg'
        file = os.path.join(folder, file)
        plt.show()
        plt.savefig(file, transparent=True)
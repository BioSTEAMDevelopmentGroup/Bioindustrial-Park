# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 00:04:53 2021

@author: yrc2
"""
import numpy as np
import biosteam as bst
from warnings import warn
from biorefineries import oilcane as oc
from ._feature_mockups import (
    all_metric_mockups, 
)
from ._load_data import (
    monte_carlo_file,
    autoload_file_name,
    spearman_file,
)

__all__ = (
    'evaluate_configurations_across_recovery_and_oil_content',
    'evaluate_configurations_across_sorghum_and_cane_oil_content',
    'evaluate_metrics_across_composition',
    'run_uncertainty_and_sensitivity',
    'save_pickled_results',
    'run_all',
)

def no_derivative(f):
    def f_derivative_disabled(*args, **kwargs):
        oc.disable_derivative()
        try:
            return f(*args, **kwargs)
        finally:
            oc.enable_derivative()
    return f_derivative_disabled

def evaluate_configurations_across_recovery_and_oil_content(
        recovery, oil_content, configurations, 
    ):
    A, B = configurations.shape
    data = np.zeros([A, B, N_metrics])
    for index, configuration in np.ndenumerate(configurations):
        oc.load(configuration)
        if oc.configuration.agile:
            oc.cane_mode.oil_content = oc.sorghum_mode.oil_content = oil_content
            oc.oil_extraction_specification.load_crushing_mill_oil_recovery(recovery)
        else:
            oc.oil_extraction_specification.load_specifications(
                crushing_mill_oil_recovery=recovery, 
                oil_content=oil_content, 
            )
        try:
            oc.sys.simulate()
        except:
            oc.sys.empty_recycles()
            oc.sys.simulate()
        data[index] = [j() for j in oc.model.metrics]
    return data

N_metrics = len(all_metric_mockups)
evaluate_configurations_across_recovery_and_oil_content = no_derivative(
    np.vectorize(
        evaluate_configurations_across_recovery_and_oil_content,
        signature=f'(),(),(a, b)->(a, b, {N_metrics})'
    )
)
def evaluate_configurations_across_sorghum_and_cane_oil_content(
        sorghum_oil_content, cane_oil_content, configurations, relative,
    ):
    C = len(configurations)
    M = len(all_metric_mockups)
    data = np.zeros([C, M])
    for ic in range(C):
        oc.load([int(configurations[ic]), True])
        oc.cane_mode.oil_content = cane_oil_content
        if relative:
            oc.sorghum_mode.oil_content = cane_oil_content + sorghum_oil_content
        else:
            oc.sorghum_mode.oil_content = sorghum_oil_content
        oc.sys.simulate()
        data[ic, :] = [j() for j in oc.model.metrics]
    return data

evaluate_configurations_across_sorghum_and_cane_oil_content = no_derivative(
    np.vectorize(
        evaluate_configurations_across_sorghum_and_cane_oil_content, 
        excluded=['configurations', 'relative'],
        signature=f'(),(),(c),()->(c,{N_metrics})'
    )
)

def evaluate_metrics_at_composition(oil, fiber, water, configuration):
    oc.load(configuration)
    cs = oc.composition_specification
    try:
        oc.load_composition(oc.feedstock, oil, fiber, water, cs.FFA, cs.PL)
    except ValueError:
        return np.array([np.nan for i in oc.model.metrics])
    oc.sys.simulate()
    # print('oil:', oil, ', fiber:', fiber, ', water:', water)
    # print(np.array([oc.model.metrics[0](), oc.model.metrics[2]()]))
    return np.array([i() for i in oc.model.metrics])

evaluate_metrics_across_composition = no_derivative(
    np.vectorize(
        evaluate_metrics_at_composition, 
        excluded=['configuration'],
        signature=f'(),(),(),()->({N_metrics})'
    )
)

def save_pickled_results(N, configurations=None, rule='L', optimize=True):
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    if configurations is None: configurations = oc.configuration_names
    # liter_per_gal = 3.7854
    for name in configurations:
        oc.load(name)
        np.random.seed(1)
        samples = oc.model.sample(N, rule)
        oc.model.load_samples(samples, optimize=optimize, ss=False)
        file = monte_carlo_file(name, False)
        oc.model.load_pickled_results(
            file=autoload_file_name(name),
            safe=False
        )
        # The energy allocated GWP was on a per gallon basis in the last Monte Carlo simulation
        # TODO: Rerun simulations and to remove this commented code block
        # oc.model.table[GWP_ethanol_allocation.index] /= liter_per_gal
        # oc.model.table[GWP_biodiesel_allocation.index] /= liter_per_gal
        oc.model.table.to_excel(file)
        oc.model.table = oc.model.table.dropna(how='all', axis=1)
        for i in oc.model.metrics:
            if i.index not in oc.model.table: oc.model._metrics.remove(i)
        oc.model.table = oc.model.table.dropna(how='any', axis=0)
        rho, p = oc.model.spearman_r()
        file = spearman_file(name)
        rho.to_excel(file)

def run_uncertainty_and_sensitivity(name, N, rule='L',
                                    across_lines=False, 
                                    sample_cache={},
                                    autosave=True,
                                    autoload=True,
                                    optimize=True):
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    oc.load(name, cache=None)
    key = (N, rule)
    if key in sample_cache:
        samples = sample_cache[key]
    else:
        np.random.seed(1)
        sample_cache[key] = samples = oc.model.sample(N, rule)
    oc.model.load_samples(samples, optimize=optimize)
    file = monte_carlo_file(name, across_lines)
    N_notify = min(int(N/10), 20)
    autosave = N_notify if autosave else False
    if across_lines:
        df = oc.get_composition_data()
        def set_line(line):
            np.random.seed(1)
            oc.set_line_composition_parameters(line)
            oc.line = line
            oc.model.load_samples(oc.model.sample(N, rule=rule), optimize=optimize)
        
        @no_derivative
        def evaluate(**kwargs):
            autoload_file = autoload_file_name(oc.line)
            oc.model.evaluate(
                autosave=autosave, 
                autoload=autoload,
                file=autoload_file,
                **kwargs,
            )
            
        oc.model.evaluate_across_coordinate(
            name='Line',
            notify=int(N/10),
            f_coordinate=set_line,
            f_evaluate=evaluate,
            coordinate=df.index,
            notify_coordinate=True,
            xlfile=file,
        )
    else:
        autoload_file = autoload_file_name(name)
        success = False
        for i in range(3):
            try:
                if name not in ('O1', 'O2'): oc.disable_derivative()
                oc.model.evaluate(
                    notify=N,
                    autosave=autosave,
                    autoload=autoload,
                    file=autoload_file
                )
            except Exception as e:
                raise e from None
                warn('failed evaluation; restarting without cache')
            else:
                success = True
                oc.enable_derivative()
                break
        if not success:
            raise RuntimeError('evaluation failed')
        oc.model.table.to_excel(file)
        oc.model.table = oc.model.table.dropna(how='all', axis=1)
        for i in oc.model.metrics:
            if i.index not in oc.model.table: oc.model._metrics.remove(i)
        oc.model.table = oc.model.table.dropna(how='any', axis=0)
        rho, p = oc.model.spearman_r(filter='omit nan')
        file = spearman_file(name)
        rho.to_excel(file)

run = run_uncertainty_and_sensitivity
    
def run_all(N, across_lines=False, rule='L', configurations=None,
            filter=None,**kwargs):
    if configurations is None: configurations = oc.configuration_names
    for name in configurations:
        if filter and not filter(name): continue
        print(f"Running {name}:")
        run_uncertainty_and_sensitivity(
            name, N, rule, across_lines, **kwargs
        )
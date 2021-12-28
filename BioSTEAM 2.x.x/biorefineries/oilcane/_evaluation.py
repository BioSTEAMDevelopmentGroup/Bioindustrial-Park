# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 00:04:53 2021

@author: yrc2
"""
import numpy as np
import biosteam as bst
from biorefineries import oilcane as oc
from ._variable_mockups import (
    all_metric_mockups
)
from ._parse_configuration import (
    parse,
    Configuration,
)
from ._load_data import (
    get_monte_carlo,
    monte_carlo_file,
    autoload_file_name,
    spearman_file,
)

__all__ = (
    'evaluate_configurations_across_extraction_efficiency_and_oil_content',
    'evaluate_configurations_across_sorghum_and_cane_oil_content',
    'evaluate_MFPP_uncertainty_across_ethanol_and_biodiesel_prices',
    'evaluate_MFPP_benefit_uncertainty_across_ethanol_and_biodiesel_prices',
    'evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices',
    'evaluate_MFPP_across_ethanol_and_biodiesel_prices',
    'evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices',
    'run_uncertainty_and_sensitivity',
    'run_all',
)

def no_derivative(f):
    def f_derivative_disabled(*args, **kwargs):
        oc.disable_derivative()
        try:
            return f(*args, **kwargs)
        finally:
            oc.enable_derivative()
    return f

def evaluate_configurations_across_extraction_efficiency_and_oil_content(
        efficiency, oil_content, oil_retention, agile, configurations,
    ):
    A = len(agile)
    C = len(configurations)
    M = len(all_metric_mockups)
    data = np.zeros([A, C, M])
    for ia in range(A):
        for ic in range(C):    
            oc.load([int(configurations[ic]), agile[ia]])
            if agile[ia]:
                oc.cane_mode.oil_content = oc.sorghum_mode.oil_content = oil_content
                oc.oil_extraction_specification.load_efficiency(efficiency)
                oc.oil_extraction_specification.load_oil_retention(oil_retention)
            else:
                oc.oil_extraction_specification.load_specifications(
                    efficiency=efficiency, 
                    oil_content=oil_content, 
                    oil_retention=oil_retention
                )
            oc.sys.simulate()
            data[ia, ic, :] = [j() for j in oc.model.metrics]
    return data

N_metrics = len(all_metric_mockups)
evaluate_configurations_across_extraction_efficiency_and_oil_content = no_derivative(
    np.vectorize(
        evaluate_configurations_across_extraction_efficiency_and_oil_content, 
        excluded=['oil_retention', 'agile', 'configurations'],
        signature=f'(),(),(),(a),(c)->(a,c,{N_metrics})'
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

def evaluate_MFPP_uncertainty_across_ethanol_and_biodiesel_prices(name, ethanol_price, biodiesel_price):
    table = get_monte_carlo(name)
    oilcane_price = table[oc.MFPP.index].to_numpy()[:, np.newaxis] # USD/ton
    biodiesel_flow = table[oc.biodiesel_production.index].to_numpy()[:, np.newaxis] * 1e6 # gal/yr
    ethanol_price_baseline = table[oc.set_ethanol_price.index].to_numpy()[:, np.newaxis]
    biodiesel_price_baseline = table[oc.set_biodiesel_price.index].to_numpy()[:, np.newaxis]
    ethanol_flow = table[oc.ethanol_production.index].to_numpy()[:, np.newaxis] * 1e6 # gal/yr
    feedstock_flow = table[oc.feedstock_consumption.index].to_numpy()[:, np.newaxis] # ton/yr
    baseline_price = (
        oilcane_price
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
@no_derivative
def evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, configuration=None):
    if configuration is not None: oc.load(configuration)
    feedstock_flow = oc.flows['feedstock']()
    biodiesel_flow = oc.flows['biodiesel']()
    ethanol_flow = oc.flows['ethanol']()
    baseline_price = (
        oc.tea.solve_price(oc.oilcane) * oc.kg_per_ton
        - (oc.ethanol.price * ethanol_flow  * 2.98668849 + oc.biodiesel.price * 3.3111 * biodiesel_flow) / feedstock_flow
    )
    return (
        baseline_price 
        + (ethanol_price * ethanol_flow + biodiesel_price * biodiesel_flow) / feedstock_flow
    )

def evaluate_MFPP_benefit_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, baseline=None, configuration=None):
    if configuration is None: configuration = oc.configuration
    if baseline is None:
        number, agile = configuration
        assert number > 0
        baseline = Configuration(-number, agile)
    MFPP_baseline = evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, baseline)
    MFPP = evaluate_MFPP_across_ethanol_and_biodiesel_prices(ethanol_price, biodiesel_price, configuration)
    return MFPP - MFPP_baseline

def run_uncertainty_and_sensitivity(name, N, rule='L',
                                    across_oil_content=False, 
                                    sample_cache={},
                                    autosave=True,
                                    autoload=True):
    np.random.seed(1)
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.utils.DesignWarning)
    oc.load(name)
    key = (N, rule)
    if key in sample_cache:
        samples = sample_cache[key]
    else:
        sample_cache[key] = samples = oc.model.sample(N, rule)
    oc.model.load_samples(samples)
    file = monte_carlo_file(name, across_oil_content)
    if across_oil_content:
        if parse(name).number < 0:
            oc.model.evaluate(notify=int(N/10))
            oc.model.evaluate_across_coordinate(
                name='Oil content',
                f_coordinate=lambda x: None,
                coordinate=oc.oil_content,
                notify=int(N/10), 
                notify_coordinate=True,
                xlfile=file,
                re_evaluate=False,
            )
        else:
            def f(x):
                oc.oil_extraction_specification.locked_oil_content = False
                oc.oil_extraction_specification.load_oil_content(x)
                oc.oil_extraction_specification.locked_oil_content = True
            oc.model.evaluate_across_coordinate(
                name='Oil content',
                f_coordinate=f,
                coordinate=oc.oil_content,
                notify=int(N/10), 
                notify_coordinate=True,
                xlfile=file,
            )
    else:
        N = min(int(N/10), 50)
        oc.model.evaluate(notify=N,
                       autosave=N if autosave else N,
                       autoload=autoload,
                       file=autoload_file_name(name))
        oc.model.table.to_excel(file)
        rho, p = oc.model.spearman_r()
        file = spearman_file(name)
        rho.to_excel(file)

run = run_uncertainty_and_sensitivity
    
def run_all(N, across_oil_content=False, rule='L', configurations=None, **kwargs):
    if configurations is None: configurations = oc.configuration_names
    for name in configurations:
        print(f"Running {name}:")
        run_uncertainty_and_sensitivity(
            name, N, rule, across_oil_content, **kwargs
        )
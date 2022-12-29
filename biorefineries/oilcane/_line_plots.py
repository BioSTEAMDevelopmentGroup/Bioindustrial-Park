# -*- coding: utf-8 -*-
"""
"""
import biorefineries.oilcane as oc
import biosteam as bst
import numpy as np
import itertools
from functools import cache
import pandas as pd
import os
from thermosteam.utils import style_axis
from warnings import filterwarnings
from thermosteam.utils import GG_colors, colors
from matplotlib import pyplot as plt

results_folder = os.path.join(os.path.dirname(__file__), 'results')
indicator_colors = (
    colors.neutral, GG_colors.purple, GG_colors.orange, GG_colors.green,
    GG_colors.red, GG_colors.blue, GG_colors.yellow,
)
indicator_sizes = (3, 3, 3,
                   2, 2, 2,
                   1, 1, 1,)
indicators = (
    'IRR', # %
    'GWP', # kg-CO2-eq / L
    'BY', # Biodiesel yield [MT / hc]
)
indicator_configurations = {
    'IRR-GWP-BY': [1, 1, 1],
    'GWP-BY': [0, 1, 1],
    'IRR-BY': [1, 0, 1],
    'IRR-GWP': [1, 1, 0],
    'IRR': [1, 0, 0],
    'GWP': [0, 1, 0],
    'BY': [0, 0, 1],
}
indicator_configurations = {i: np.array(j).reshape([3, 1]) for i, j in indicator_configurations.items()}
lines = ('WT', 'Energycane', '1566')
maximize_indicator_guide = {
    indicators[0]: True,
    indicators[1]: False,
    indicators[2]: True,
}

def create_scenarios(ranges, n=5):
    linspaces = [np.linspace(*i, n) for i in ranges]
    scenarios = itertools.product(*linspaces)
    return np.array(list(scenarios))

@cache
def create_model(configuration=None):
    if configuration is None:
        configuration = 'O6'
    if configuration == 'mockup':
        model = bst.Model(None)
    else:
        oc.load(configuration)
        model = bst.Model(oc.sys)
    parameter = model.parameter
    metric = model.metric
    
    @parameter(bounds=(50, 100), units='% sugarcane', baseline=60)
    def set_oilcane_biomass_yield(oilcane_biomass_yield):
        if oc.feedstock.imass['Lipid'] / (oc.feedstock.F_mass - oc.feedstock.imass['Water']) > 0.001:
            oc.dry_biomass_yield = oc.baseline_dry_biomass_yield * oilcane_biomass_yield * 0.01
    
    @parameter(bounds=(50, 90), units='%', baseline=60)
    def set_fermentation_yield(fermentation_yield):
        fermentation_yield /= 100
        for i in oc.fermentation_rxnsys.reactions:
            i.X[:2] = fermentation_yield
            i.X[2:] = 0.99 - fermentation_yield
    
    parameter(oc.set_cane_glucose_yield, name='Glucose yield')
    parameter(oc.set_cane_xylose_yield, name='Xylose yield')
    parameter(oc.set_saccharification_oil_recovery, name='Oil recovery')

    metric(oc.IRR)
    metric(oc.GWP_biodiesel_displacement)
    metric(oc.biodiesel_yield)
    return model

def set_cane_line(line):
    df = oc.get_composition_data()
    try:
        line = df.loc[line]
    except:
        line = df.loc[int(line)]
    oil = line['Stem Oil (dw)']
    moisture = line['Water (wt)']
    fiber = line['Fiber (dw)']
    biomass = line['Dry biomass yield (MT/hc)']
    oc.composition_specification.oil = oil['Mean']
    oc.composition_specification.fiber = fiber['Mean']
    oc.composition_specification.moisture = moisture['Mean']
    oc.composition_specification.load_composition()
    oc.dry_biomass_yield = biomass['Mean']
    GWP_sugarcane = oc.GWP_characterization_factors['sugarcane']
    oc.feedstock.set_CF(
        oc.GWP, GWP_sugarcane * (0.9 * oc.baseline_dry_biomass_yield / oc.dry_biomass_yield + 0.10)
    )

indicator_results = None   
def evaluate_winners():
    global indicator_results
    filterwarnings('ignore')
    model = create_model()
    scenarios = create_scenarios([i.bounds for i in model.parameters])
    model.load_samples(scenarios)
    if indicator_results is None:
        indices = [i.index for i in model.metrics]
        indicator_results = {i: [] for i in indicators}
        for line in lines:
            set_cane_line(line)
            model.evaluate()
            for indicator, index in zip(indicators, indices):
                indicator_results[indicator].append(model.table[index].copy())
    indicator_winners = np.vstack([
        (np.argmax if maximize_indicator_guide[i] else np.argmin)(
            np.vstack(indicator_results[i]), axis=0, keepdims=True,
        )
        for i in indicators
    ])
    
    winner_results = {i: {j: [] for j in indicator_configurations} for i in lines}
    for n, line in enumerate(lines):
        score_board = indicator_winners == n
        for indicator, configuration in indicator_configurations.items():
            winner_results[line][indicator] = scenarios[np.where((score_board == configuration).all(axis=0))]
        
    for line in lines:
        line_results = winner_results[line]
        xlfile = os.path.join(results_folder, f"{line}_winner.xlsx")
        with pd.ExcelWriter(xlfile) as writer:
            for indicator in indicator_configurations:
                pd.DataFrame(
                    line_results[indicator], columns=[i.name for i in model.parameters]
                ).to_excel(
                    writer, 
                    sheet_name=indicator
                )
        
def plot_winner_lines():
    model = create_model()
    fig, axes = plt.subplots(len(lines))
    parameters = model.parameters
    x = list(range(len(parameters)))
    
    def normalize_parameter_values(values):
        lbs = np.array([i.bounds[0] for i in parameters])
        ubs = np.array([i.bounds[1] for i in parameters])
        return (values - lbs) / (ubs - lbs)
    
    for ax, line in zip(axes, lines):
        plt.sca(ax)
        xlfile = os.path.join(results_folder, f"{line}_winner.xlsx")
        for indicator, color, size in zip(indicator_configurations, indicator_colors, indicator_sizes):
            df = pd.read_excel(xlfile, index_col=0, sheet_name=indicator)
            values = normalize_parameter_values(df.values)
            for y in values: 
                plt.plot(x, y, c=color.RGBn, lw=size)
        plt.ylim([0, 1])
        style_axis(
            ax, xticks=x, yticks=[0, 1],
            xticklabels=[f"{i.bounds[0]}\n{i.name}\n[{i.units}]" for i in parameters],
            yticklabels=['Min', 'Max'],
        )
        plt.sca(ax._cached_ytwin)
        plt.xticks(x, [f"{i.bounds[1]}" for i in parameters])
    plt.subplots_adjust(right=0.92, hspace=1, top=0.9, bottom=0.10)
            
                
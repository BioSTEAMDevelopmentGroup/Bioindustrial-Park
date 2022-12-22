# -*- coding: utf-8 -*-
"""
"""
import biorefineries.oilcane as oc
import biosteam as bst
import thermosteam as tmo
import numpy as np
import itertools
from functools import cache
import pandas as pd
import os
from warnings import filterwarnings
from thermosteam.utils import GG_colors
from matplotlib import pyplot as plt

results_folder = os.path.join(os.path.dirname(__file__), 'results')
indicator_colors = (GG_colors.red, GG_colors.blue, GG_colors.yellow)
indicators = ('IRR [%]', 'GWP [kg-CO2-eq / L]', 'biodiesel yield [MT / hc]')
lines = ('WT', 'Energycane', '1566')
maximize_indicator_guide = {
    indicators[0]: True,
    indicators[1]: False,
    indicators[2]: True,
}

def create_scenarios(ranges, n=3):
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
    
    @parameter(bounds=(50, 100), units='% Sugarcane', baseline=60)
    def set_oilcane_biomass_yield(biomass_yield):
        oc.dry_biomass_yield = oc.baseline_dry_biomass_yield * biomass_yield * 0.01
    
    @parameter(bounds=(50, 90), units='%', baseline=60)
    def set_fermentation_yield(fermentation_yield):
        fermentation_yield /= 100
        for i in oc.fermentation_rxnsys.reactions:
            i.X[:2] = fermentation_yield
            i.X[2:] = 0.99 - fermentation_yield
    
    parameter(oc.set_cane_glucose_yield)
    parameter(oc.set_cane_xylose_yield)
    parameter(oc.set_saccharification_oil_recovery)

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
    
def evaluate_winners():
    filterwarnings('ignore')
    model = create_model()
    scenarios = create_scenarios([i.bounds for i in model.parameters])
    model.load_samples(scenarios)
    indicator_results = {i: [] for i in indicators}
    indices = [i.index for i in model.metrics]
    for line in lines:
        set_cane_line(line)
        model.evaluate()
        for indicator, index in zip(indicators, indices):
            indicator_results[indicator].append(model.table[index].copy())
    winner_results = {i: {j: [] for j in indicators} for i in lines}
    for indicator in indicators:
        f = np.argmax if maximize_indicator_guide[indicator] else np.argmin
        winners = f(np.hstack(indicator_results[indicator]), axis=0)
        for scenario, winner in zip(scenarios, winners):
            line = lines[winner]
            winner_results[line][indicator].append(scenario)
    for line in lines:
        line_results = winner_results[line]
        xlfile = os.path.join(results_folder, f"{line}_winner.xlsx")
        with pd.ExcelWriter(xlfile) as writer:
            for indicator in indicators:
                pd.DataFrame(
                    line_results[indicator], columns=[i.name for i in model.parameters]
                ).to_excel(
                    writer, 
                    sheet_name=indicator.split('[')[0].strip(' ')
                )
        
def plot_winner_lines():
    model = bst.create_model('mockup')
    axes = plt.subplots(len(lines))
    parameters = model.parameters
    x = list(range(len(parameters)))
    
    def normalize_parameter_values(values):
        lbs = np.array([i.bounds[0] for i in parameters])
        ubs = np.array([i.bounds[1] for i in parameters])
        return (values - lbs) / (ubs - lbs)
    
    for ax, line in zip(axes, lines):
        plt.sca(ax)
        xlfile = os.path.join(results_folder, f"{line}_winner.xlsx")
        for color, indicator in zip(indicator_colors, indicators):
            df = pd.read_excel(xlfile, sheet_name=pd.read_excel())
            values = normalize_parameter_values(df.values())
            for y in values: 
                plt.plot(x, y, c=color.RGBn, alpha=0.5)
        
        tmo.plots.style_axis(
            ax, xticks=x, yticks=[0, 1],
            xticklabels=[f"{i.min} {i.units}\n{i.name}" for i in parameters],
            yticklabels=['Min', 'Max'],
        )
        plt.sca(ax._cached_ytwin)
        plt.xticks(x, [f"{i.max} {i.units}" for i in parameters])
            
                
                    
if __name__ == '__main__':
    evaluate_winners()
        
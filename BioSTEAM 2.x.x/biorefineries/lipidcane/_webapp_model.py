# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import matplotlib.pyplot as plt
import biosteam as bst
from biorefineries import lipidcane as lc
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import triang
from biosteam.process_tools import UnitGroup
from biosteam import plots
import pandas as pd
import numpy as np

np.random.seed(1234)

# %% Create metrics 

ethanol_density = 2.9866 # kg/gallon
tea = lc.lipidcane_tea
ugroup = UnitGroup('Biorefinery', tea.units)
get_steam_demand = lambda: sum([i.flow for i in lc.BT.steam_utilities]) * 18.01528 * tea._operating_hours / 1000
get_electricity_consumption = ugroup.get_electricity_consumption
get_electricity_production = ugroup.get_electricity_production
get_excess_electricity = lambda: get_electricity_production() - get_electricity_consumption()
get_ethanol_production = lambda: lc.ethanol.F_mass * tea._operating_hours / 907.185
get_biodiesel_production = lambda: lc.biodiesel.F_mass * tea._operating_hours / 907.185
get_MESP = lambda: tea.solve_price(lc.ethanol) * ethanol_density
get_MFPP = lambda: tea.solve_price(lc.lipidcane)
get_IRR = lambda: tea.solve_IRR() * 100
get_FCI = lambda: tea.FCI / 10**6

metrics = (Metric('IRR', get_IRR, '%'),
           Metric('MFPP', get_MFPP, 'USD/kg'),
           Metric('MESP', get_MESP, 'USD/gal'), 
           Metric('FCI', get_FCI, 'million USD'),
           Metric('Biodiesel production', get_biodiesel_production, 'ton/yr'),
           Metric('Ethanol production', get_ethanol_production, 'ton/yr'),
           Metric('Steam demand', get_steam_demand, 'MT/yr'),
           Metric('Consumed electricity', get_electricity_consumption, 'MW'),
           Metric('Excess electricity', get_excess_electricity, 'MW'))

# Used to index metrics by name
metrics_by_name = {i.name: i for i in metrics}

# %% Create model and populate parameters

#: [float] Plant size process specification in ton / yr
plant_size_ = lc.lipidcane.F_mass * lc.lipidcane_tea.operating_days * 24 / 907.185
constant_lipid_content = False
def load_specifications():
    # Plant size
    lc.lipidcane.F_mass = plant_size_ / (lc.lipidcane_tea.operating_days * 24 / 907.185)

model = Model(lc.lipidcane_sys, metrics, load_specifications)
param = model.parameter
lipid_content = 100. * lc.utils.get_lipid_fraction()

@param(units='dry wt. %', 
       distribution=triang(lipid_content),
       baseline=lipid_content)
def set_lipic_content(lipid_content):
    if constant_lipid_content: return
    lc.utils.set_lipid_fraction(lipid_content / 100.)

@param(units='ton/yr',
       distribution=triang(plant_size_),
       baseline=plant_size_)
def set_plant_size(plant_size):
    global plant_size_
    plant_size_ = plant_size

@param(units='days/yr', 
       distribution=triang(tea.operating_days),
       baseline=tea.operating_days)
def set_operating_days(operating_days):
    tea.operating_days = operating_days
    
@param(element=lc.ethanol, units='USD/kg',
       distribution=triang(lc.ethanol.price),
       baseline=lc.ethanol.price)
def set_ethanol_price(price):
    lc.ethanol.price = price
    
@param(element=lc.lipidcane, units='USD/kg',
       distribution=triang(lc.lipidcane.price),
       baseline=lc.lipidcane.price)
def set_lipidcane_price(price):
    lc.lipidcane.price = price

@param(element=bst.PowerUtility, units='USD/kWhr',
       distribution=triang(bst.PowerUtility.price),
       baseline=bst.PowerUtility.price)
def set_electricity_price(price):
    bst.PowerUtility.price = price

@param(units='%',
       distribution=triang(100. * tea.IRR),
       baseline=100. * tea.IRR)
def set_IRR(IRR):
    tea.IRR = IRR / 100.

model.get_parameters()[-1].name = 'IRR'

# %% Interface

def evaluate_sample(lipid_content, plant_size, operating_days, ethanol_price,
                    lipidcane_price, electricity_price, IRR):
    """
    Return a pandas.Series object with metric values.
    
    Examples
    --------
    >>> # baseline_sample = model.get_baseline_sample()
    >>> # baseline_sample
    array([1.012e+01, 1.764e+06, 2.000e+02, 7.890e-01, 3.455e-02, 7.820e-02, 1.804e+01])
    >>> # metric_values = evaluate_sample(*baseline_sample)
    >>> # metric_values
    IRR [%]                                     0.18
    MFPP [USD/kg]                               0.0346
    MESP [USD/gal]                              2.36
    FCI [million USD]                           222
    Biodiesel production [ton/yr]               4.62e+04
    Ethanol production [ton/yr]                 6.14e+04
    Steam demand [MT/yr]                        5.67e+05
    Consumed electricity [MW]                   4.63e+04
    Excess electricity [MW]                     3.98e+05
    dtype: float64
    >>> metric_values['Consumed electricity [MW]']
    46295.52888978263
    
    """
    return model([lipid_content, plant_size, operating_days, ethanol_price,
                  lipidcane_price, electricity_price, IRR])['Biorefinery']

def plot_montecarlo(metric, N=100):
    """
    Create a matplotlib boxplot and return plot objects.
    
    Examples
    --------
    >>> # plot_montecarlo('Excess electricity', 100) -> creates plot
    {'whiskers': [<matplotlib.lines.Line2D at 0x27ac8469b88>,
      <matplotlib.lines.Line2D at 0x27ab70a9c08>],
     'caps': [<matplotlib.lines.Line2D at 0x27ac852ed88>,
      <matplotlib.lines.Line2D at 0x27ac8551f48>],
     'boxes': [<matplotlib.patches.PathPatch at 0x27ac84692c8>],
     'medians': [<matplotlib.lines.Line2D at 0x27ac855b2c8>],
     'fliers': [<matplotlib.lines.Line2D at 0x27ac8551fc8>],
     'means': []}
    
    """
    if isinstance(metric, str): metric = metrics_by_name[metric]
    samples = model.sample(N, 'L')
    model.load_samples(samples)
    model.evaluate()
    df = model.table[[metric.index]]
    bx = plots.plot_montecarlo(df, transpose=True)
    plt.xticks([], [])
    plt.ylabel(metric.name_with_units)
    return bx

def plot_montecarlo_across_lipid_content(metric,
                                         N_coordinate_samples = 50,
                                         N_coordinates = 10,
                                         rule = 'L'):
    """
    Create a montecarlo plot across lipid content and return figure.
    
    Examples
    --------
    >>> # plot_montecarlo_across_lipid_content('IRR') -> creates plot
    
    """
    global constant_lipid_content
    if isinstance(metric, str): metric = metrics_by_name[metric]
    
    ### Evaluate ###
    
    # Monte Carlo across lipid fraction    
    coordinate = np.linspace(0.11, 0.01, N_coordinates)
    samples = model.sample(N_coordinate_samples, rule)
    model.load_samples(samples)
    
    constant_lipid_content = True
    model.evaluate_across_coordinate('Lipid fraction',
          lc.utils.set_lipid_fraction, coordinate,
          xlfile='Monte Carlo across lipid fraction.xlsx')
    constant_lipid_content = False    
    
    # TODO: Add model for sugarcane
    # # Sugar cane Monte Carlo    
    # samples = model_sc.sample(N_coordinate_samples, rule)
    # model_sc.load_samples(samples)
    # model_sc.evaluate()
    # model_sc.table.to_excel('Monte Carlo sugarcane.xlsx')

    ### Plot ###

    fig = plt.figure()
    
    # Plot lipid-cane    
    readxl = lambda sheet: pd.read_excel('Monte Carlo across lipid fraction.xlsx',
                                         sheet_name=sheet, index_col=0)
    df_metric = readxl(metric.name)
    lipid_content = np.array(df_metric.columns) * 100
    plt.ylabel(metric.name_with_units)
    plots.plot_montecarlo_across_coordinate(lipid_content, df_metric)
    plt.xlim([1, 11])
    plt.xlabel("Lipid content [dry wt. %]")
    
    # # Plot sugarcane
    # data_sc = pd.read_excel('Monte Carlo sugarcane.xlsx', header=[0,1])
    # get_metric = lambda name: np.asarray(data_sc['Biorefinery', name]).flatten()
    # metric_values = get_metric(metric.name)
    # plots.plot_montecarlo(metric_values)
    # plots.plot_vertical_line(1)
    # plt.xlim(-1, 10.5)
    
    return fig

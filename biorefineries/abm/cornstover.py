# -*- coding: utf-8 -*-
"""
This module exposes the `ABM_TEA_model`, a biosteam.Model object
for evaluating multiple scenarios for the agent based model project.

Examples
--------
>>> import numpy as np
>>> from biorefineries.cornstover.abm import ABM_TEA_model
>>> cornstover_price = [0.05159, 0.06] # Corn stover price [USD/kg]
>>> miscanthus_price = [0.08, 0.09] # Miscanthus price [USD/kg]
>>> ethanol_price = [0.80, 0.85] # Ethanol price [USD/kg]
>>> operating_days = [350.4, 350.4] # Operating days [day/yr]
>>> IRR = [0.1, 0.15] # Internal rate of return
>>> start_year = [2020, 2020]
>>> end_year = [2050, 2050]
>>> plant_capacity = [876072883, 876072883] # Plant capacity [kg/yr]
>>> x_cornstover = [0.1, 0.9] # Corn stover fraction [by wt.]
>>> samples = np.array([
...     cornstover_price, 
...     miscanthus_price,
...     ethanol_price,
...     operating_days,
...     IRR,
...     start_year,
...     end_year,
...     plant_capacity,
...     x_cornstover,
... ]).T
>>> ABM_TEA_model.load_samples(samples)
>>> ABM_TEA_model.evaluate()
>>> table = ABM_TEA_model.table # A pandas data frame
>>> MESP_index = ABM_TEA_model.metrics[0].index
>>> ABM_TEA_model.table[MESP_index] 
0    0.75
1   0.927
Name: (Biorefinery, MESP [USD/kr]), dtype: float64
>>> parameters = ABM_TEA_model.get_parameters()
>>> cornstover_price_index = parameters[0].index
>>> ABM_TEA_model.table[cornstover_price_index] 
0   0.0516
1     0.06
Name: (Corn stover, Price [USD/kg]), dtype: float64

"""
from biorefineries import cornstover as cs
import biosteam as bst
import os
import pandas as pd

__all__ = ('ABM_TEA_model',)

cs.load()
bst.CE = 607.5

# %% Composition utilities

# From NREL/TP-5100-47764
cornstover_dry_composition = cs.chemicals.kwarray(
    dict(Glucan=0.319,
         Xylan=0.189,
         Galactan=0.015, 
         Arabinan=0.028, 
         Mannan=0.0030, 
         Lignin=0.133,
         Acetate=0.022,
         Protein=0.037,
         Extract=0.143, 
         Ash=0.07, 
         Sucrose=0.036)
)
cornstover_dry_composition /= cornstover_dry_composition.sum()
cellulosic_moisture_content = cs.chemicals.kwarray(
    dict(Water=0.20)
)
cornstover_composition = 0.8 * cornstover_dry_composition + cellulosic_moisture_content

# Miscanthus compostion from DOI 10.1007/s12257-016-0485-x
# 5.3% was from other unnamed compounds, here we assume that that 5.3% of 
# compounds present are in the same ratios as in cornstover.
miscanthus_dry_composition = cs.chemicals.kwarray(
    dict(Glucan=0.434,
         Xylan=0.256,
         Galactan=0.015, 
         Arabinan=0.017,  
         Lignin=0.202)
)
others_composition = cs.chemicals.kwarray(
    dict(Acetate=0.022,
         Protein=0.037,
         Extract=0.143, 
         Ash=0.07)
)
others_composition /= others_composition.sum()
miscanthus_dry_composition += (1. - miscanthus_dry_composition.sum()) * others_composition
miscanthus_composition = 0.8 * miscanthus_dry_composition + cellulosic_moisture_content

def set_mixed_cornstover_miscanthus_feedstock(x_cornstover):
    x_miscanthus = 1 - x_cornstover
    composition = (x_cornstover * cornstover_composition
                   + x_miscanthus * miscanthus_composition)
    cs.cornstover.mass = cs.cornstover.F_mass * composition

# %% Functional ABM TEA model for legacy purposes

def ABM_TEA_function(
        cornstover_fraction=1.0,
        operating_days=350.4,
        plant_capacity=876072883.4242561, 
        price_cornstover=0.05159, 
        price_miscanthus=0.08, 
        price_ethanol=0.80,
        IRR=0.10,
        duration=(2007, 2037),
    ):
    """
    Return a dictionary of biorefinery metrics for the production of cellulosic
    ethanol from mixed feedstocks.

    Parameters
    ----------
    cornstover_fraction : float
        Fractino of cornstover in feedstock.
    operating_days : float
        Number of operating days per year.
    plant_capacity : float
        Plat capacity in kg/yr of feedstock.
    price_cornstover : float
        Price of cornstover in USD/kg.
    price_miscanthus : float
        Price of miscanthus in USD/kg.
    price_ethanol : float
        Price of ethanol in USD/kg.
    IRR : float
        Internal rate of return as a fraction (not percent!).
    duration : tuple(int, int)
        Years of operation.

    Returns
    -------
    metrics: dict
        Includes MESP [USD/kg], MFPP [USD/kg], IRR [-], NPV [USD], 
        TCI [USD], FOC [USD/yr], VOC [USD/yr], Electricity consumption [MWhr/yr], 
        Electricity production [MWhr/yr], and Production [kg/yr].
    
    """
    x_cornstover = cornstover_fraction
    if not 0. <= x_cornstover <= 1.:
        raise ValueError('cornstover fraction must be between 0 to 1; {x_cornstover} given')
    set_mixed_cornstover_miscanthus_feedstock(x_cornstover)
    cs.cornstover.price = (price_cornstover * x_cornstover 
                           + price_miscanthus * (1 - x_cornstover))
    cs.ethanol.price = price_ethanol
    hours = operating_days * 24 
    cs.cornstover.F_mass = plant_capacity / hours
    cs.cornstover_tea.operating_days = operating_days
    cs.cornstover_tea.duration = duration
    cs.cornstover_sys.simulate()
    cs.cornstover_tea.IRR = IRR
    unit_group = cs.AllAreas
    return {
        'MESP': cs.cornstover_tea.solve_price(cs.ethanol),
        'MFPP': cs.cornstover_tea.solve_price(cs.cornstover),
        'IRR': cs.cornstover_tea.solve_IRR(),
        'NPV': cs.cornstover_tea.NPV,
        'TCI': cs.cornstover_tea.TCI,
        'VOC': cs.cornstover_tea.VOC,
        'FOC': cs.cornstover_tea.FOC,
        'Electricity consumption [MWhr/yr]': hours * unit_group.get_electricity_consumption(), 
        'Electricity production [MWhr/yr]': hours * unit_group.get_electricity_production(),
        'Production': cs.ethanol.F_mass * hours,
    }

# %% ABM Model object

metrics = [
    bst.Metric('MESP', lambda:cs.cornstover_tea.solve_price(cs.ethanol), 'USD/kr'),    
    bst.Metric('MFPP', lambda:cs.cornstover_tea.solve_price(cs.cornstover), 'USD/kr'),
    bst.Metric('IRR', cs.cornstover_tea.solve_IRR),
    bst.Metric('NPV', lambda:cs.cornstover_tea.NPV, 'USD'),
    bst.Metric('TCI', lambda:cs.cornstover_tea.TCI, 'USD'),
    bst.Metric('VOC', lambda:cs.cornstover_tea.VOC, 'USD/yr'),
    bst.Metric('FOC', lambda:cs.cornstover_tea.FOC, 'USD/yr'),
    bst.Metric('Electricity consumption', lambda:operating_hours * cs.AllAreas.get_electricity_consumption(), 'MWhr/yr'),
    bst.Metric('Electricity production', lambda:operating_hours * cs.AllAreas.get_electricity_production(), 'MWhr/yr'),
    bst.Metric('Production', lambda:operating_hours * cs.ethanol.F_mass, 'kg/yr')
]

ABM_TEA_model = bst.Model(cs.cornstover_sys, metrics)

@ABM_TEA_model.parameter(element='Corn stover', units='USD/kg')
def set_cornstover_price(price):
    global price_cornstover
    price_cornstover = price

@ABM_TEA_model.parameter(element='Miscanthus', units='USD/kg')
def set_miscanthus_price(price):
    global price_miscanthus
    price_miscanthus = price

@ABM_TEA_model.parameter(element='Ethanol', units='USD/kg')
def set_ethanol_price(price):
    cs.ethanol.price = price

operating_hours = cs.cornstover_tea.operating_days * 24
@ABM_TEA_model.parameter(units='day/yr')
def set_operating_days(operating_days):
    global operating_hours
    cs.cornstover_tea.operating_days = operating_days
    operating_hours = operating_days * 24

@ABM_TEA_model.parameter
def set_IRR(IRR):
    cs.cornstover_tea.IRR = IRR

@ABM_TEA_model.parameter
def set_start_year(start_year):
    cs.cornstover_tea.duration = (start_year, cs.cornstover_tea.duration[1])
    
@ABM_TEA_model.parameter
def set_end_year(end_year):
    cs.cornstover_tea.duration = (cs.cornstover_tea.duration[0], end_year)

@ABM_TEA_model.parameter(units='kg/yr')
def set_plant_capacity(plant_capacity):
    cs.cornstover.F_mass = plant_capacity / operating_hours

ABM_TEA_model.parameter(set_mixed_cornstover_miscanthus_feedstock,
                        name='Corn stover fraction', units='by wt.')


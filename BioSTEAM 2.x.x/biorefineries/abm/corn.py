# -*- coding: utf-8 -*-
"""
This module exposes the `ABM_TEA_model`, a biosteam.Model object
for evaluating multiple scenarios for the agent based model project.

Examples
--------
>>> import numpy as np
>>> from biorefineries.corn.abm import ABM_TEA_model
>>> corn_price = [0.0847, 0.09] # Corn price [USD/kg]
>>> DDGS_price = [0.05159, 0.0968] # DDGS price [USD/kg]
>>> corn_oil_price = [0.5, 0.6] # Corn oil price [USD/kg]
>>> ethanol_price = [0.485, 0.5] # Ethanol price [USD/kg]
>>> operating_days = [330., 330.] # Operating days [day/yr]
>>> IRR = [0.1, 0.15] # Internal rate of return
>>> start_year = [2020, 2020]
>>> end_year = [2040, 2040]
>>> plant_capacity = [876072883, 876072883] # Plant capacity [kg/yr]
>>> samples = np.array([
...     corn_price, 
...     DDGS_price,
...     corn_oil_price,
...     ethanol_price,
...     operating_days,
...     IRR,
...     start_year,
...     end_year,
...     plant_capacity,
... ]).T
>>> ABM_TEA_model.load_samples(samples)
>>> ABM_TEA_model.evaluate()
>>> table = ABM_TEA_model.table # A pandas data frame
>>> MESP_index = ABM_TEA_model.metrics[0].index
>>> ABM_TEA_model.table[MESP_index] 
0   0.385
1   0.386
Name: (Biorefinery, MESP [USD/kr]), dtype: float64
>>> parameters = ABM_TEA_model.get_parameters()
>>> cornstover_price_index = parameters[0].index
>>> ABM_TEA_model.table[cornstover_price_index] 
0   0.0847
1     0.09
Name: (Corn, Price [USD/kg]), dtype: float64

"""
from biorefineries import corn as cn
import biosteam as bst

__all__ = ('ABM_TEA_model',)

cn.load()
bst.CE = 607.5

# %% Functional ABM TEA model for legacy purposes

def ABM_TEA_function(
        operating_days=330,
        plant_capacity=876072883.4242561, 
        price_corn=0.13227735731092652, 
        price_DDGS=0.12026, 
        price_corn_oil=0.56,
        price_ethanol=0.48547915353569393,
        IRR=0.15,
        duration=(2007, 2027),
    ):
    """
    Return a dictionary of biorefinery metrics for the production of cellulosic
    ethanol from mixed feedstocks.

    Parameters
    ----------
    operating_days : float
        Number of operating days per year.
    plant_capacity : float
        Plat capacity in kg/yr of feedstock.
    price_corn : float
        Price of corn in USD/kg.
    price_DDGS : float
        Price of DDGS in USD/kg.
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
    cn.corn.price = price_corn
    cn.DDGS.price = price_DDGS
    cn.crude_oil.price = price_corn_oil
    cn.ethanol.price = price_ethanol
    hours = operating_days * 24 
    cn.corn.F_mass = plant_capacity / hours
    cn.corn_tea.operating_days = operating_days
    cn.corn_tea.duration = duration
    cn.corn_sys.simulate()
    cn.corn_tea.IRR = IRR
    unit_group = cn.all_areas
    return {
        'MESP': cn.corn_tea.solve_price(cn.ethanol),
        'MFPP': cn.corn_tea.solve_price(cn.corn),
        'IRR': cn.corn_tea.solve_IRR(),
        'NPV': cn.corn_tea.NPV,
        'TCI': cn.corn_tea.TCI,
        'VOC': cn.corn_tea.VOC,
        'FOC': cn.corn_tea.FOC,
        'Electricity consumption [MWhr/yr]': hours * unit_group.get_electricity_consumption(), 
        'Electricity production [MWhr/yr]': hours * unit_group.get_electricity_production(),
        'Production': hours * cn.ethanol.F_mass,
    }

# %% ABM Model object

metrics = [
    bst.Metric('MESP', lambda:cn.corn_tea.solve_price(cn.ethanol), 'USD/kr'),    
    bst.Metric('MFPP', lambda:cn.corn_tea.solve_price(cn.corn), 'USD/kr'),
    bst.Metric('IRR', cn.corn_tea.solve_IRR),
    bst.Metric('NPV', lambda:cn.corn_tea.NPV, 'USD'),
    bst.Metric('TCI', lambda:cn.corn_tea.TCI, 'USD'),
    bst.Metric('VOC', lambda:cn.corn_tea.VOC, 'USD/yr'),
    bst.Metric('FOC', lambda:cn.corn_tea.FOC, 'USD/yr'),
    bst.Metric('Electricity consumption', lambda:operating_hours * cn.all_areas.get_electricity_consumption(), 'MWhr/yr'),
    bst.Metric('Electricity production', lambda:operating_hours * cn.all_areas.get_electricity_production(), 'MWhr/yr'),
    bst.Metric('Production', lambda:operating_hours * cn.ethanol.F_mass, 'kg/yr')
]

ABM_TEA_model = bst.Model(cn.corn_sys, metrics)

@ABM_TEA_model.parameter(element='Corn', units='USD/kg')
def set_corn_price(price):
    cn.corn.price = price

@ABM_TEA_model.parameter(element='DDGS', units='USD/kg')
def set_DDGS_price(price):
    cn.DDGS.price = price
    
@ABM_TEA_model.parameter(element='Corn oil', units='USD/kg')
def set_corn_oil_price(price):
    cn.crude_oil.price = price

@ABM_TEA_model.parameter(element='Ethanol', units='USD/kg')
def set_ethanol_price(price):
    cn.ethanol.price = price

operating_hours = cn.corn_tea.operating_days * 24
@ABM_TEA_model.parameter(units='day/yr')
def set_operating_days(operating_days):
    global operating_hours
    cn.corn_tea.operating_days = operating_days
    operating_hours = operating_days * 24

@ABM_TEA_model.parameter
def set_IRR(IRR):
    cn.corn_tea.IRR = IRR

@ABM_TEA_model.parameter
def set_start_year(start_year):
    cn.corn_tea.duration = (start_year, cn.corn_tea.duration[1])
    
@ABM_TEA_model.parameter
def set_end_year(end_year):
    cn.corn_tea.duration = (cn.corn_tea.duration[0], end_year)

@ABM_TEA_model.parameter(units='kg/yr')
def set_plant_capacity(plant_capacity):
    cn.corn.F_mass = plant_capacity / operating_hours

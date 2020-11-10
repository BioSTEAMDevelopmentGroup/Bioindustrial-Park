# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:48:17 2020

@author: yrc2
"""
from biorefineries import corn as cn

__all__ = ('ABM_TEA_model',)


def ABM_TEA_model(
        operating_days=350.4,
        plant_capacity=876072883.4242561, 
        price_corn=0.08476585075177462, 
        price_DDGS=0.09687821462905594, 
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
        'Production': cn.ethanol.F_mass * operating_days * 24.,
    }
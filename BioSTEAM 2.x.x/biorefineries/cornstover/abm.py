# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 18:48:17 2020

@author: yrc2
"""
from biorefineries import cornstover as cs

__all__ = ('ABM_TEA_model',)

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

def ABM_TEA_model(
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
        'Production': cs.ethanol.F_mass * operating_days * 24.,
    }
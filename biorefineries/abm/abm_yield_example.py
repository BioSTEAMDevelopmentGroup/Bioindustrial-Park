# -*- coding: utf-8 -*-
"""
Example module for calculating ethanol yield

"""
from biorefineries.abm.cornstover import ABM_TEA_model as lignocellulose_model
from biorefineries.abm.corn import ABM_TEA_model as corn_model

lignocellulosic_capacity = 876072883.4242561
corn_capacity = 876072883.4242561
specific_ethanol_volume = 1.267 # L/kg
kg_per_ton = 907.185
def compute_ethanol_yield(plant_capacity, production):
    return (
        specific_ethanol_volume
        * kg_per_ton
        * production
        / plant_capacity
    )

# Yield of Corn
results = corn_model(plant_capacity=corn_capacity)
ethanol_yield = compute_ethanol_yield(corn_capacity, results['Production'])
print('Corn ethanol yield:', round(ethanol_yield), 'L/ton')

# Yield of Corn Stover
results = lignocellulose_model(plant_capacity=lignocellulosic_capacity, cornstover_fraction=1)
ethanol_yield = compute_ethanol_yield(lignocellulosic_capacity, results['Production'])
print('Corn stover ethanol yield:', round(ethanol_yield), 'L/ton')

# Yield of Minscanthus
results = lignocellulose_model(plant_capacity=lignocellulosic_capacity, cornstover_fraction=0)
ethanol_yield = compute_ethanol_yield(lignocellulosic_capacity, results['Production'])
print('Miscanthus ethanol yield:', round(ethanol_yield), 'L/ton')
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 22:48:34 2021

@author: yrc2
"""
from matplotlib import pyplot as plt
import numpy as np
from chaospy import distributions as shape

__all__ = ('ethanol_price_distribution', 
           'biodiesel_price_distribution',
           'natural_gas_price_distribution')

ethanol_prices = [ # USD / gal
    1.97, 2.02, 2.23, 2.11, 2.17, 2.54, 3.18, 3.01, 2.81, 2.76, 3.13, 2.86,
    2.87, 2.89, 3.06, 2.97, 2.77, 2.72, 2.92, 2.84, 2.46, 1.64, 1.97, 2.28,
    1.76, 1.18, 1.59, 1.72, 1.67, 1.72, 1.78, 1.72, 1.96, 1.99, 2.18, 2.24, 
    2.1, 1.59, 1.59, 1.98, 1.93, 1.81, 0.94, 1.34, 0.89
]

biodiesel_prices = [ # USD / gal
    3.02, 3.12, 3.38, 3.38, 3.27, 3.3, 3.32, 3.33, 3.21, 3.18, 3.18, 3.3, 
    3.35, 3.66, 3.9, 4.25, 4.49, 4.6, 4.82, 5.07, 5.23, 5.4, 5.45, 5.44, 5.74, 
    5.58, 5.34, 4.92, 4.67, 4.67, 4.79, 4.97, 4.76, 4.4, 4.39, 4.47, 4.53,
    4.49, 4.12, 3.94, 4.37, 4.51, 4.73, 4.75, 4.84, 4.85, 5.01, 4.9, 4.84,
    4.56, 3.96, 3.9, 3.39, 3.37, 3.67, 3.75, 3.76, 3.7, 3.66, 3.48, 3.36,
    3.14, 3.15, 3.03, 2.83, 2.88, 2.9, 2.89, 3.08, 3.18, 2.96, 2.55, 2.6,
    2.67, 2.65, 2.84, 2.79, 2.92, 3.17, 3.21, 3.28, 3.31, 3.11, 3.23, 3.35,
    3.37, 3.49, 3.64, 3.52, 3.45, 3.14, 3.11, 3.18, 3.13, 3.23, 3.31, 3.33,
    3.3, 3.34, 3.23, 3.45, 3.16, 3.04, 3.02, 3.13, 3.11, 3.02, 3.01, 3.04,
    3.09, 3.09, 3.06, 3.12, 3.24, 3.08, 2.98, 2.86, 2.87, 2.91, 2.99, 3.06, 
    3.13, 3.21, 3.26, 3.33, 3.11, 2.97, 2.77, 2.74, 2.87, 2.9, 3.07, 3.02,
    3.02, 3.47, 3.61, 3.85, 4.92, 5.0
]

natural_gas_prices = [ # City gate [USD / cf] 2010 to 2019
    6.18, 5.63, 4.73, 4.88, 5.71, 4.26, 3.71, 4.16, 4.23, 3.81,
]

electricity_prices = []
def triangular_distribution(x):
    a, b, c = fit_triangular_distribution(x)
    return shape.Triangle(a, c, b)

def fit_triangular_distribution(x):
    a = min(x)
    b = max(x)
    c = np.mean(x)
    return a, b, c

def plot_triangular_distribution(a, b, c):
    return plt.plot([a, c, b], [0, 2 / (b - a), 0])

def plot_histogram(x, *args, bins=10, density=True, **kwargs):
    return plt.hist(x, *args, **kwargs)

ethanol_price_distribution = triangular_distribution(ethanol_prices)
biodiesel_price_distribution = triangular_distribution(biodiesel_prices)
natural_gas_price_distribution = triangular_distribution(natural_gas_prices)
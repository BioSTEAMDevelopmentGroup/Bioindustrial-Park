# -*- coding: utf-8 -*-
"""
Created on Fri May 14 22:48:34 2021

@author: yrc2
"""
from matplotlib import pyplot as plt
import numpy as np
from chaospy import distributions as shape
from math import sqrt

__all__ = (
    'ethanol_price_distribution', 
    'biodiesel_price_distribution',
    'biodiesel_minus_ethanol_price_distribution',
    'natural_gas_price_distribution',
    'mean_glycerol_price',
    'mean_ethanol_price',
    'mean_biodiesel_price',
    'mean_natural_gas_price',
    'mean_electricity_price',
    'mean_soymeal_price',
)

ethanol_prices = [ # USD / gal Dec 2009 - Nov 2020, by quarter
    2.08, 1.66, 1.71, 2.33, 2.37, 2.67, 2.87, 2.83, 2.30, 2.26, 2.47,
    2.48, 2.39, 2.58, 2.57, 2.31, 2.19, 2.76, 2.28, 2.16, 1.83, 1.63, 1.64,
    1.61, 1.47, 1.52, 1.59, 1.58, 1.57, 1.60, 1.62, 1.24, 1.02, 1.30, 1.37, 
    1.20, 1.09, 1.09, 1.32, 1.35
]

biodiesel_prices = [ # USD / gal Dec 2009 - Nov 2020, by month
    3.38, 3.27, 3.3, 3.32, 3.33, 3.21, 3.18, 3.18, 3.3, 
    3.35, 3.66, 3.9, 4.25, 4.49, 4.6, 4.82, 5.07, 5.23, 5.4, 5.45, 5.44, 5.74, 
    5.58, 5.34, 4.92, 4.67, 4.67, 4.79, 4.97, 4.76, 4.4, 4.39, 4.47, 4.53,
    4.49, 4.12, 3.94, 4.37, 4.51, 4.73, 4.75, 4.84, 4.85, 5.01, 4.9, 4.84,
    4.56, 3.96, 3.9, 3.39, 3.37, 3.67, 3.75, 3.76, 3.7, 3.66, 3.48, 3.36,
    3.14, 3.15, 3.03, 2.83, 2.88, 2.9, 2.89, 3.08, 3.18, 2.96, 2.55, 2.6,
    2.67, 2.65, 2.84, 2.79, 2.92, 3.17, 3.21, 3.28, 3.31, 3.11, 3.23, 3.35,
    3.37, 3.49, 3.64, 3.52, 3.45, 3.14, 3.11, 3.18, 3.13, 3.23, 3.31, 3.33,
    3.3, 3.34, 3.23, 3.45, 3.16, 3.04, 3.02, 3.13, 3.11, 3.02, 3.01, 3.04,
    3.09, 3.09, 3.06, 3.12, 3.24, 3.08, 2.98, 2.86, 2.87, 2.91, 2.99, 3.06, 
    3.13, 3.21, 
]


N_quarters = int(len(biodiesel_prices)/3)
biodiesel_prices_quarter = [ # By quarter
    np.mean(biodiesel_prices[3*i: 3*i+3]) for i in range(N_quarters)
]

biodiesel_minus_ethanol_prices = [i - j for i, j in zip(biodiesel_prices, ethanol_prices)]

natural_gas_prices = [ # City gate [USD / cf] 2010 to 2019
    6.18, 5.63, 4.73, 4.88, 5.71, 4.26, 3.71, 4.16, 4.23, 3.81,
]

electricity_prices = []
def triangular_distribution(x):
    a, b, c = fit_triangular_distribution(x)
    return shape.Triangle(a, c, b)

def triangular_distribution_median(a, b, c):
    if c >= (a + b) / 2: 
        return a + sqrt((b - a) * (c - a) * 0.5)
    else:
        return b - sqrt((b - a) * (b - c) * 0.5)

def fit_triangular_distribution(x, median=False):
    a = min(x)
    b = max(x)
    if median:
        median = np.median(x)
        if median >= (a + b) / 2: 
            dummy = a - median
            c = a + 2 * dummy * dummy / (b - a)
        else:
            dummy = b - median
            c = b - 2 * dummy * dummy / (b - a)
        assert triangular_distribution_median(a, b, c) == median 
    else:
        c = 3 * np.mean(x) - (a + b)
        if c < a: c = a + 1e-6
    return a, b, c

def plot_triangular_distribution(a, b, c):
    return plt.plot([a, c, b], [0, 2 / (b - a), 0])

def plot_histogram(x, *args, bins=10, density=True, **kwargs):
    return plt.hist(x, *args, **kwargs)


# https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
# Data from historical prices, 2010-2020
electricity_price_distribution = shape.Triangle(0.0583, 0.065, 0.069)
ethanol_price_distribution = triangular_distribution(ethanol_prices)
biodiesel_price_distribution = triangular_distribution(biodiesel_prices)
biodiesel_minus_ethanol_price_distribution = triangular_distribution(biodiesel_minus_ethanol_prices)
natural_gas_price_distribution = triangular_distribution(natural_gas_prices)

mean_glycerol_price = (0.10 + 0.22) * 0.5 
mean_ethanol_price = np.mean(ethanol_prices)
mean_biodiesel_price = np.mean(biodiesel_prices)
mean_natural_gas_price = np.mean(natural_gas_prices)
mean_electricity_price = sum([0.0583, 0.065, 0.069]) / 3.
mean_soymeal_price = 0.33 # 10 yr average; https://markets.businessinsider.com/commodities/soybean-meal-price?op=1
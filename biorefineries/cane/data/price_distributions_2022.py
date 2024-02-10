# -*- coding: utf-8 -*-
"""
Created on Fri May 14 22:48:34 2021

@author: yrc2
"""
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import os
from chaospy import distributions as shape
from math import sqrt

__all__ = (
    'ethanol_no_RIN_price_distribution', 
    'biodiesel_no_RIN_price_distribution',
    'advanced_ethanol_price_distribution',
    'cellulosic_ethanol_price_distribution',
    'biomass_based_diesel_price_distribution',
    'cellulosic_based_diesel_price_distribution',
    'natural_gas_price_distribution',
    'mean_glycerol_price',
    'mean_ethanol_no_RIN_price',
    'mean_advanced_ethanol_price',
    'mean_cellulosic_ethanol_price',
    'mean_biomass_based_diesel_price',
    'mean_cellulosic_based_diesel_price',
    'mean_natural_gas_price',
    'mean_electricity_price',
    'mean_soymeal_price',
)

liter_per_gal = 3.7854
RINfile = os.path.join(os.path.dirname(__file__), 'RIN_prices_2022.csv')
df_RIN = pd.read_csv(RINfile)
df_D3 = df_RIN[df_RIN['Fuel (D Code)'] == 'D3'] # Cellulosic ethanol
df_D4 = df_RIN[df_RIN['Fuel (D Code)'] == 'D4'] # Biodiesel
df_D5 = df_RIN[df_RIN['Fuel (D Code)'] == 'D5'] # Advanced ethanol
df_D6 = df_RIN[df_RIN['Fuel (D Code)'] == 'D6'] # Conventional ethanol

def summarize_RIN_credit_by_year(df, month_offset):
    data = {}
    for date, price in zip(df['Transfer Date by Week'], df['RIN Price']):
        month, day, year = date.split('/')
        year = int(year)
        if year in data:
            data_by_year = data[year]
        else:
            data_by_year = data[year] = []
        if int(month) - month_offset < 0: year -= 1
        data_by_year.append(float(price.strip('$')))
    
    for year, data_by_year in data.items():
        data[year] = np.mean(data_by_year)
    return data

def get_RIN_credit_by_year_June2017_to_May2022(df):
    offset = 6
    dct = summarize_RIN_credit_by_year(df, offset)
    years = (2017, 2018, 2019, 2020, 2021)
    return np.array([dct[y] for y in years])

# June 2017 to May 2022 (start of year is offset by 6 months)
RIN_D3_prices = get_RIN_credit_by_year_June2017_to_May2022(df_D3) / liter_per_gal
RIN_D4_prices = get_RIN_credit_by_year_June2017_to_May2022(df_D4) / liter_per_gal
RIN_D5_prices = get_RIN_credit_by_year_June2017_to_May2022(df_D5) / liter_per_gal
RIN_D6_prices = get_RIN_credit_by_year_June2017_to_May2022(df_D6) / liter_per_gal


# %% Raw data

ethanol_no_RIN_prices = np.array([ # USD / gal May 2017 - Dec 2020, by month
    1.6, 1.62, 1.64, 1.55, 1.19, 0.98, 0.95, 1.0, 1.11, 1.25, 1.31,
    1.34, 1.4, 1.41, 1.31, 1.2, 1.21, 1.19, 1.08, 1.09, 1.12, 1.18, 1.17,
    1.21, 1.36, 1.31, 1.3, 1.22, 1.39, 1.44, 1.32, 1.15, 1.03, 0.84, 0.59,
    0.75, 0.85, 0.91, 0.87, 0.95, 0.94, 0.78, 0.6, 
]) / liter_per_gal

ethanol_with_D6RIN_prices = np.array([ # USD / gal Jan 2021  - May 2022, by month
    1.55, 1.73, 1.87, 2.15, 2.62, 2.53, 2.37, 2.35, 2.52, 2.62, 3.42, 3.11, 
    2.22, 2.17, 2.49, 2.64, 2.82
]) / liter_per_gal

assert len(ethanol_no_RIN_prices) + len(ethanol_with_D6RIN_prices) == 5 * 12

# Includes RIN D4
biomass_based_diesel_prices = np.array([ # USD / gal June 2017 - June 2022, by month
    3.13, 3.23, 3.31, 3.33, 3.3, 3.34, 3.23, 3.45, 3.16, 3.04, 3.02, 3.13,
    3.11, 3.02, 3.01, 3.04, 3.09, 3.09, 3.06, 3.12, 3.24, 3.08, 2.98, 2.86,
    2.87, 2.91, 2.99, 3.06, 3.13, 3.21, 3.26, 3.35, 3.11, 2.97, 2.77, 2.74, 
    2.87, 2.92, 3.12, 3.3, 2.88, 3.51, 3.61, 4.08, 4.5, 4.91, 5.0, 5.5, 5.99, 
    5.95, 5.91, 5.71, 5.6, 5.49, 4.88, 4.83, 6.03, 6.5, 6.99, 7.53, 7.35
]) / liter_per_gal

def by_month_to_year(arr):
    N_years = int(len(arr)/12)
    return np.array([ 
        np.mean(arr[12*i: 12*i+12]) for i in range(N_years)
    ])

biomass_based_diesel_prices = by_month_to_year(biomass_based_diesel_prices)
biodiesel_no_RIN_prices = biomass_based_diesel_prices - RIN_D4_prices * 1.5
cellulosic_based_diesel_prices = biodiesel_no_RIN_prices + RIN_D3_prices * 1.5
ethanol_no_RIN_prices = by_month_to_year(ethanol_no_RIN_prices)
ethanol_with_D6RIN_prices = by_month_to_year(ethanol_with_D6RIN_prices)
ethanol_no_RIN_prices = np.array([
    *ethanol_no_RIN_prices, 
    *(ethanol_with_D6RIN_prices - RIN_D6_prices[-2:]),
])
assert len(ethanol_no_RIN_prices) == 5

natural_gas_prices = np.array([ # City gate [USD / mcf] 2017 to 2021, mcf = 1,000 cf
    4.16, 4.23, 3.81, 3.43, 6.11
]) * 35.3146667/1e3 # To USD / m3

crude_oil_prices = np.array([ # USD / barrel West Texas Intermediate 2018 - 2022; https://www.macrotrends.net/1369/crude-oil-price-history-chart
    65.23, 56.99, 39.68, 68.17, 94.53
])

# %% Shorten data to 2017 - 2022 because of data limitations for cellulosic ethanol RINs 

# Account for RIN credits
cellulosic_ethanol_prices = ethanol_no_RIN_prices + RIN_D3_prices
advanced_ethanol_prices = ethanol_no_RIN_prices + RIN_D5_prices

# plt.plot(cellulosic_ethanol_prices)
# plt.plot(advanced_ethanol_prices)
# plt.plot(biomass_based_diesel_prices)

def triangular_distribution(x, median=False):
    a, b, c = fit_triangular_distribution(x, median)
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
    return plt.plot([a, b, c], [0, 2 / (b - a), 1])

def plot_histogram(x, *args, bins=10, density=True, **kwargs):
    return plt.hist(x, *args, **kwargs)

# Price distributions
electricity_price_distribution = shape.Triangle(0.0583, 0.065, 0.069) # Data from historical prices, 2010-2020; https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
RIN_D3_price_distribution = triangular_distribution(RIN_D3_prices)
RIN_D4_price_distribution = triangular_distribution(RIN_D4_prices)
RIN_D5_price_distribution = triangular_distribution(RIN_D5_prices)
ethanol_no_RIN_price_distribution = triangular_distribution(ethanol_no_RIN_prices)
advanced_ethanol_price_distribution = triangular_distribution(advanced_ethanol_prices)
cellulosic_ethanol_price_distribution = triangular_distribution(cellulosic_ethanol_prices)
biodiesel_no_RIN_price_distribution = triangular_distribution(biodiesel_no_RIN_prices)
biomass_based_diesel_price_distribution = triangular_distribution(biomass_based_diesel_prices)
cellulosic_based_diesel_price_distribution = triangular_distribution(cellulosic_based_diesel_prices)
natural_gas_price_distribution = triangular_distribution(natural_gas_prices)
crude_oil_price_distribution = triangular_distribution(crude_oil_prices)

# https://www.eia.gov/outlooks/aeo/pdf/00%20AEO2021%20Chart%20Library.pdf
# Data from historical prices, 2010-2020
electricity_price_distribution = shape.Triangle(0.0583, 0.065, 0.069)

# Mean values
mean_glycerol_price = (0.10 + 0.22) * 0.5 
mean_RIN_D3_price = np.mean(RIN_D3_prices)
mean_RIN_D4_price = np.mean(RIN_D4_prices)
mean_RIN_D5_price = np.mean(RIN_D5_prices)
mean_ethanol_no_RIN_price = np.mean(ethanol_no_RIN_prices)
mean_biodiesel_no_RIN_price = np.mean(biodiesel_no_RIN_prices)
mean_cellulosic_ethanol_price = np.mean(cellulosic_ethanol_prices)
mean_advanced_ethanol_price = np.mean(advanced_ethanol_prices)
mean_biomass_based_diesel_price = np.mean(biomass_based_diesel_prices)
mean_cellulosic_based_diesel_price = np.mean(cellulosic_based_diesel_prices)
mean_natural_gas_price = np.mean(natural_gas_prices)
mean_electricity_price = sum([0.0583, 0.065, 0.069]) / 3.
mean_soymeal_price = 0.33 # 10 yr average; https://markets.businessinsider.com/commodities/soybean-meal-price?op=1
mean_crude_oil_price = np.mean(crude_oil_prices)
mcop = mean_crude_oil_price

# Short hand distributions names
cbpd = cellulosic_based_diesel_price_distribution 
mcbp = mean_cellulosic_based_diesel_price
bpd = biomass_based_diesel_price_distribution 
mbp = mean_biomass_based_diesel_price
cepd = cellulosic_ethanol_price_distribution 
mcep = mean_cellulosic_ethanol_price
aepd = advanced_ethanol_price_distribution 
maep = mean_advanced_ethanol_price
copd = crude_oil_price_distribution
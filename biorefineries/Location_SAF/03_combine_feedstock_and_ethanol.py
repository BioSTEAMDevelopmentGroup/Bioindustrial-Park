#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 11:44:30 2024

@author: bianco3

This file is used to combine uncertainty results from the feedstock transport model with the 
BioSTEAM ethanol refinery model simulation uncertainty results, to get ethanol MSP and CI 
for all candidate locations, with uncertainty.


"""

#%% Import packages
import pandas as pd
import geopandas as gpd
import numpy as np
#from scipy.spatial.distance import cdist
#import libpysal as ps
#import matplotlib.pyplot as plt
#from pyproj import Transformer
#import math
#import matplotlib as mpl
import os
#import ast
#from scipy.optimize import curve_fit

#%% Load results from previos model files

# load existing arrays for chosen feedstock

feedstock = 'switchgrass'

if feedstock == 'switchgrass':
    name = ''
    
elif feedstock == 'miscanthus':
    name = '_mis'

jet_producers_states_supplied = np.load(f'jet_producers_states_supplied{name}.npy')
feedstock_delivered_price = np.load(f'Feedstock_price{name}.npy') # in $/ wet kg crop (miscanthus)
ethanol_unit_transp_cost_each_jet = np.load(f'Ethanol_transport_cost{name}.npy') # in $/Mg ethanol
sent_ethanol_no_uncertainty = np.load(f'sent_ethanol_no_uncertainty{name}.npy') # in Mg ethanol/yr

possible_locations_with_state = gpd.read_file(f"possible_locations_with_state{name}.shp") # per feedstock


#%% Function to calculate stats from uncertatinty results

def calculate_stats_uncertainty(array):
    array_data = array
    # Calculate the required statistics for each location (across the 1000 samples)
    mean_values = np.mean(array_data, axis=0)
    variance_values = np.var(array_data, axis=0)
    percentile_5th = np.percentile(array_data, 5, axis=0)
    percentile_95th = np.percentile(array_data, 95, axis=0)

    # Combine the stats into a single dictionary for easy access
    stats = {
        'mean': mean_values,
        'variance': variance_values,
        '5th_percentile': percentile_5th,
        '95th_percentile': percentile_95th
    }
    return stats

#%% Combine with ethanol price equations


# Get folder path of this file
folder = os.path.dirname(__file__)
# Join the folder path with the folder name where Ethanol results are located
input_data_folder = os.path.join(folder, 'Ethanol_results')

# read dataframes of a0 and b0 coefficients by state
file_name = f'Ethanol_price_for_python{name}.xlsx'
file_path = os.path.join(input_data_folder, file_name)
sheet_name = 'a0'
a0 = pd.read_excel(file_path, sheet_name=sheet_name)

sheet_name = 'b0'
b0 = pd.read_excel(file_path, sheet_name=sheet_name)

#%% Calculate ethanol MSP for 1000 locations and 1000 samples
ethanol_price_y_list = []
for i in range(1000):
    y_poss_locations_i = []
    for j in range(1000):
        State_name = possible_locations_with_state.iloc[j].NAME
        #print(State_name)
        a0_j = a0[State_name][i]
        b0_j = b0[State_name][i]
        eth_price_y_j = feedstock_delivered_price[i][j] * a0_j + b0_j
        y_poss_locations_i.append(eth_price_y_j)
        
    ethanol_price_y_list.append(np.array(y_poss_locations_i)) 

ethanol_price_y_array = np.array(ethanol_price_y_list)

stats_ethanol_price = calculate_stats_uncertainty(ethanol_price_y_array)

df = pd.DataFrame(stats_ethanol_price)

# Sort the DataFrame by the 'mean' column
df_sorted = df.sort_values(by='mean').reset_index()

#%% Conversion units

liters_per_gal = 3.785 # L/gal
liters_per_m3 = 1000 # L/m3
density_of_ethanol = 747.58 # from biosteam in kg/m3
kg_per_ton = 1000 # kg/Mg (Metric tonnes)
_conversion_USDperMg_to_gal_ethanol = 1/kg_per_ton * density_of_ethanol * 1/liters_per_m3 * liters_per_gal

#%% Weighted ethanol price for all jet producers

# Compute the weighted prices
weighted_prices = ethanol_unit_transp_cost_each_jet * sent_ethanol_no_uncertainty

# Sum the weighted prices and weights
sum_weighted_prices = np.sum(weighted_prices, axis=2)
sum_weights = np.sum(sent_ethanol_no_uncertainty, axis=1)

# Compute the weighted average prices
# Avoid division by zero by using np.where to handle cases where sum_weights is zero
weighted_average_cost_Mg = np.where(sum_weights != 0, sum_weighted_prices / sum_weights, 0)

# weighted_average_cost_Mg is in S/Mg ethanol, I need to transform it to $/ kg ethanol (input units to BioSTEAM next model)
weighted_average_cost_kg = weighted_average_cost_Mg / 1000 # in $/ kg ethanol

weighted_average_cost_gal_array = weighted_average_cost_Mg * _conversion_USDperMg_to_gal_ethanol # in $/gal

ethanol_delivered_price = np.array(ethanol_price_y_list) + weighted_average_cost_gal_array

ethanol_unit_transport_per_gal = ethanol_unit_transp_cost_each_jet * _conversion_USDperMg_to_gal_ethanol

ethanol_delivered_price_each_jet = ethanol_unit_transport_per_gal + np.expand_dims(np.array(ethanol_price_y_list), axis=-1) # in $/gal ethanol

np.save(f'ethanol_delivered_price_each_jet{name}.npy', ethanol_delivered_price_each_jet) # save 

#%% Import files to combine feedstock CI with ethanol CI with uncertainty

# load arrays
feedstock_delivered_GHG = np.load(f'Feedstock_GHG{name}.npy') # in kgCO2eq/ wet kg crop 
ethanol_unit_transp_GHG_each_jet = np.load(f'Ethanol_transport_CI{name}.npy') # in kgCO2eq/Mg ethanol

#%% Combine with ethanol CI equations


# read dataframes of a0 and b0 coefficients by state
file_name = f'Ethanol_GWP_for_python{name}.xlsx'
file_path = os.path.join(input_data_folder, file_name)
sheet_name = 'a0'
a0 = pd.read_excel(file_path, sheet_name=sheet_name)

sheet_name = 'b0'
b0 = pd.read_excel(file_path, sheet_name=sheet_name)

def calculate_stats_uncertainty_array(array):
    array_data = array
    # Calculate the required statistics for each location (across the 1000 samples)
    mean_values = np.mean(array_data, axis=0)
    variance_values = np.var(array_data, axis=0)
    percentile_5th = np.percentile(array_data, 5, axis=0)
    percentile_95th = np.percentile(array_data, 95, axis=0)

    # Combine the stats into a single dictionary for easy access
    stats = {
        'mean': mean_values,
        'variance': variance_values,
        '5th_percentile': percentile_5th,
        '95th_percentile': percentile_95th
    }
    return stats

# Initialize an empty array with the shape (1000, 1000)
ethanol_GHG_y_array = np.zeros((1000, 1000))

state_names = possible_locations_with_state['NAME'].values  # Extract state names as a NumPy array

for i in range(1000):  # Loop over each sample
    # Extract the a0 and b0 values for the i-th sample (row) across all states (columns)
    a0_i = a0.loc[i, state_names].values  # Select the i-th row and state name columns for a0
    b0_i = b0.loc[i, state_names].values  # Same for b0
    
    # Perform the vectorized calculation for all possible locations
    eth_GHG_y_i = feedstock_delivered_GHG[i, :] * a0_i + b0_i
    
    # Assign the result to the i-th row of the final array
    ethanol_GHG_y_array[i, :] = eth_GHG_y_i

# Now ethanol_GHG_y_array is a NumPy array with the shape (1000, 1000)
# ethanol_GHG is in gCO2e/MJ ethanol

# calculate stats for ethanol GHG
stats_ethanol_GHG = calculate_stats_uncertainty_array(ethanol_GHG_y_array)

df = pd.DataFrame(stats_ethanol_GHG)

# Sort the DataFrame by the 'mean' column
df_sorted = df.sort_values(by='mean').reset_index()

#%% Compute the weighted CI for all jet producers

# 
weighted_CI = ethanol_unit_transp_GHG_each_jet * sent_ethanol_no_uncertainty

# Sum the weighted prices and weights
sum_weighted_CI = np.sum(weighted_CI, axis=2)
sum_weights = np.sum(sent_ethanol_no_uncertainty, axis=1)

# Compute the weighted average prices
# Avoid division by zero by using np.where to handle cases where sum_weights is zero
weighted_average_GHG_Mg = np.where(sum_weights != 0, sum_weighted_CI / sum_weights, 0)

# weighted_average_GHG_Mg is in kgCO2e/Mg ethanol, I need to transform it to kgCO2eq/ kg ethanol for BioSTEAM next model input
weighted_average_GHG_kg = weighted_average_GHG_Mg / 1000 # in kgCO2eq/ kg ethanol

# ethanol_GHG_y_array is in gCO2 eq/MJ ethanol, I have to transform it to kgCO2e/kg ethanol
ethanol_GHG_kgCO2_per_kg = ethanol_GHG_y_array / 1000 * 29.7 # 29.7 is the energy density of ethanol. Source: https://energyeducation.ca/encyclopedia/Ethanol#cite_note-5
# kgCO2 eq/ kg ethanol
np.save(f'ethanol_GHG_kgCO2_per_kg{name}.npy', ethanol_GHG_kgCO2_per_kg)

ethanol_delivered_GHG = ethanol_GHG_kgCO2_per_kg + weighted_average_GHG_kg # in kgCO2 eq/ kg ethanol

#%% Save weighted average results

np.save(f'ethanol_delivered_price{name}.npy', ethanol_delivered_price) 
np.save(f'ethanol_delivered_CI{name}.npy', ethanol_delivered_GHG) 

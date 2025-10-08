#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 14:24:15 2024

@author: bianco3

# This model creates the final results for SAF price and CI, and Pareto Front

Choose feedstock in line 25
"""

#%% Import packages
import numpy as np
import pandas as pd
#from scipy.optimize import curve_fit
import time
#import matplotlib.pyplot as plt
#import geopandas as gpd
import os

#%% Load necessary arrays

#choose feedstock
feedstock = 'switchgrass' # or switchgrass

if feedstock == 'switchgrass':
    name = ''
elif feedstock == 'miscanthus':
    name = '_mis'
    
jet_producers_states_supplied = np.load(f'jet_producers_states_supplied{name}.npy')
ethanol_delivered_price_each_jet = np.load(f'ethanol_delivered_price_each_jet{name}.npy') # $/gal
sent_ethanol_no_uncertainty = np.load(f'sent_ethanol_no_uncertainty{name}.npy') # in Mg ethanol/yr

#%% Conversion units

# $/gal ethanol to $/kg

liters_per_gal = 3.785 # L/gal
liters_per_m3 = 1000 # L/m3
density_of_ethanol = 747.58 # biosteam # from internet # kg/m3
kg_per_ton = 1000 # kg/Mg (Metric tonnes)
_conversion_USDperGal_to_kg_ethanol = 1/liters_per_gal * liters_per_m3 * 1/density_of_ethanol

#%% Definition of function to load states data

def load_all_state_data(file_path, price_type, sizes, start_col=48):
    data_dict = {}
    
    # Loop through each size and load the corresponding Excel sheet
    for size in sizes:
        sheet_name = f'Price{price_type}_{size}'
        try:
            # Load the specific sheet into a DataFrame
            df = pd.read_excel(file_path, sheet_name=sheet_name, header=0)
            
            # Convert column titles to lowercase
            df.columns = df.columns.str.lower()
            
            # Select only the columns starting from start_col
            df_selected = df.iloc[:, start_col:]
            
            # Store the selected DataFrame in the dictionary with size as the key
            data_dict[size] = df_selected
        except Exception as e:
            print(f"Error loading sheet {sheet_name}: {e}")
    
    return data_dict

#%%

ethanol_delivered_price_each_jet_perkg = ethanol_delivered_price_each_jet/ liters_per_gal * liters_per_m3/density_of_ethanol
# in $/kg

#%%

# Load the entire Excel file into a dictionary of DataFrames
# Get folder path of this file
folder = os.path.dirname(__file__)
# Join the folder path with the folder name where Ethanol results are located
input_data_folder = os.path.join(folder, 'SAF_results')

# read dataframes of a0 and b0 coefficients by state
file_name = 'Jet fuel price for python.xlsx'
sizes = [5, 14, 31, 32, 44, 45, 62, 71, 80]  # Refinery sizes in MM Gal/yr
file_path = os.path.join(input_data_folder, file_name)

# Load all state data once
y_values_price_0 = load_all_state_data(file_path, price_type=0, sizes=sizes)
y_values_price_25 = load_all_state_data(file_path, price_type=25, sizes=sizes)

#%%

# Define number of candidate locations
num_points = 1000

# Load all necessary sheets into memory
sheet_names = {
    226367.: ('a1_80', 'b1_80'),
    210909.: ('a1_71', 'b1_71'),
    184545.: ('a1_62', 'b1_62'),
    134094.: ('a1_45', 'b1_45'),
    131818.: ('a1_44', 'b1_44'),
    94549.: ('a1_32', 'b1_32'),
    92273.: ('a1_31', 'b1_31'),
    41822.: ('a1_14', 'b1_14'),
    15458.: ('a1_5', 'b1_5'),
}

# Dictionary to store loaded DataFrames
data_frames = {}

for key, (a1_sheet, b1_sheet) in sheet_names.items():
    data_frames[a1_sheet] = pd.read_excel(file_path, sheet_name=a1_sheet)
    data_frames[b1_sheet] = pd.read_excel(file_path, sheet_name=b1_sheet)

# Initialize the jet_prices array
jet_prices = np.zeros((1000, 1000, 2))

# Dictionary to store previously calculated a1 and b1 values
calculated_values = {}

# Start processing
for loc in range(num_points):  # num_points = 1000
    for jet in range(2):  # Two refineries
        start_time = time.time() # Start timing
        state_name = jet_producers_states_supplied[loc][jet].lower()
        #size = np.round(sent_ethanol_array[0], 0)[loc][jet]
        size = np.round(sent_ethanol_no_uncertainty,0)[loc][jet]

        if size == 0.:
            a1 = np.zeros(1000)
            b1 = np.zeros(1000)

        elif size in calculated_values and state_name in calculated_values[size]: # verificar
            # Retrieve cached values
            a1, b1 = calculated_values[size][state_name]
            
        elif size in sheet_names:
            # Access the preloaded data from the dictionary
            a1_sheet = data_frames[sheet_names[size][0]]
            b1_sheet = data_frames[sheet_names[size][1]]
            
            # Ensure the state name matches case-insensitively
            matching_columns_a1 = [col for col in a1_sheet.columns if col.lower() == state_name]
            matching_columns_b1 = [col for col in b1_sheet.columns if col.lower() == state_name]

            if matching_columns_a1 and matching_columns_b1:
                a1 = np.array(a1_sheet[matching_columns_a1[0]])
                b1 = np.array(b1_sheet[matching_columns_b1[0]])
            else:
                raise ValueError(f"State '{state_name}' not found in sheets for size {size}.")
        

        # Generate jet prices for this location and refinery
        x = ethanol_delivered_price_each_jet_perkg[:, loc, jet]  # Shape (1000, )
        jet_prices[:, loc, jet] = a1 * x + b1

        elapsed_time = time.time() - start_time  # End timing
        print(f"Elapsed time for location {loc}, jet {jet}: {elapsed_time:.2f} seconds")

# save results
np.save(f'jet_prices_separated{name}.npy', jet_prices)

#%% 
# compute SAF weighted average price for each jet producer
# Compute the weighted prices
weighted_prices = jet_prices * sent_ethanol_no_uncertainty

# Sum the weighted prices and weights
sum_weighted_prices = np.sum(weighted_prices, axis=2)
sum_weights = np.sum(sent_ethanol_no_uncertainty, axis=1)

# Compute the weighted average prices
# Avoid division by zero by using np.where to handle cases where sum_weights is zero
jet_prices_average = np.where(sum_weights != 0, sum_weighted_prices / sum_weights, 0)
np.save(f'jet_price_ave{name}.npy', jet_prices_average)

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

jet_price_stats = calculate_stats_uncertainty(jet_prices_average)

df = pd.DataFrame(jet_price_stats)

# Sort the DataFrame by the 'mean' column
df_sorted = df.sort_values(by='mean').reset_index()

df_price = df_sorted

#%% Load arrays for CI calculations

ethanol_unit_transp_GHG_each_jet = np.load(f'Ethanol_transport_CI{name}.npy') # in kgCO2eq/Mg ethanol
ethanol_GHG_kgCO2_per_kg = np.load(f'ethanol_GHG_kgCO2_per_kg{name}.npy') # in kgCO2e/kg ethanol

#%%
ethanol_unit_transp_GHG_each_jet_per_kg = ethanol_unit_transp_GHG_each_jet / 1000 # in in kgCO2eq/kg ethanol 

ethanol_delivered_GHG_each_jet = ethanol_unit_transp_GHG_each_jet_per_kg + np.expand_dims(np.array(ethanol_GHG_kgCO2_per_kg), axis=-1)

#%%

def load_all_state_data(file_path, GHG_type, sizes, start_col=48):
    data_dict = {}
    
    # Loop through each size and load the corresponding Excel sheet
    for size in sizes:
        sheet_name = f'GHG{GHG_type}_{size}'
        try:
            # Load the specific sheet into a DataFrame
            df = pd.read_excel(file_path, sheet_name=sheet_name, header=0)
            
            # Convert column titles to lowercase
            df.columns = df.columns.str.lower()
            
            # Select only the columns starting from start_col
            df_selected = df.iloc[:, start_col:]
            
            # Store the selected DataFrame in the dictionary with size as the key
            data_dict[size] = df_selected
        except Exception as e:
            print(f"Error loading sheet {sheet_name}: {e}")
    
    return data_dict

# Load the entire Excel file into a dictionary of DataFrames


file_name = 'Jet_fuel_GHG_for_python.xlsx'
sizes = [5, 14, 31, 32, 44, 45, 62, 71, 80]  # Refinery sizes in MM Gal/yr
file_path = os.path.join(input_data_folder, file_name)

# Load all state data once
y_values_GHG_0 = load_all_state_data(file_path, GHG_type=0, sizes=sizes)
y_values_GHG_19 = load_all_state_data(file_path, GHG_type=19, sizes=sizes)

#%%
# Load all necessary sheets into memory
sheet_names = {
    226367.: ('a1_80', 'b1_80'),
    210909.: ('a1_71', 'b1_71'),
    184545.: ('a1_62', 'b1_62'),
    134094.: ('a1_45', 'b1_45'),
    131818.: ('a1_44', 'b1_44'),
    94549.: ('a1_32', 'b1_32'),
    92273.: ('a1_31', 'b1_31'),
    41822.: ('a1_14', 'b1_14'),
    15458.: ('a1_5', 'b1_5'),
}

# Dictionary to store loaded DataFrames
data_frames = {}

for key, (a1_sheet, b1_sheet) in sheet_names.items():
    data_frames[a1_sheet] = pd.read_excel(file_path, sheet_name=a1_sheet)
    data_frames[b1_sheet] = pd.read_excel(file_path, sheet_name=b1_sheet)

# Initialize the jet_prices array
jet_CIs = np.zeros((1000, 1000, 2))

# Dictionary to store previously calculated a1 and b1 values
calculated_values = {}

# Start processing
for loc in range(num_points):  # num_points = 1000
    for jet in range(2):  # Two refineries
        start_time = time.time() # Start timing
        state_name = jet_producers_states_supplied[loc][jet].lower()
        #size = np.round(sent_ethanol_array[0], 0)[loc][jet]
        size = np.round(sent_ethanol_no_uncertainty,0)[loc][jet]

        if size == 0.:
            a1 = np.zeros(1000)
            b1 = np.zeros(1000)

        elif size in calculated_values and state_name in calculated_values[size]: # verificar
            # Retrieve cached values
            a1, b1 = calculated_values[size][state_name]
            
        elif size in sheet_names:
            # Access the preloaded data from the dictionary
            a1_sheet = data_frames[sheet_names[size][0]]
            b1_sheet = data_frames[sheet_names[size][1]]
            
            # Ensure the state name matches case-insensitively
            matching_columns_a1 = [col for col in a1_sheet.columns if col.lower() == state_name]
            matching_columns_b1 = [col for col in b1_sheet.columns if col.lower() == state_name]

            if matching_columns_a1 and matching_columns_b1:
                a1 = np.array(a1_sheet[matching_columns_a1[0]])
                b1 = np.array(b1_sheet[matching_columns_b1[0]])
            else:
                raise ValueError(f"State '{state_name}' not found in sheets for size {size}.")
                
        
        # Generate jet CIs for this location and refinery
        x = ethanol_delivered_GHG_each_jet[:, loc, jet]  # Shape (1000, )
        jet_CIs[:, loc, jet] = a1 * x + b1

        elapsed_time = time.time() - start_time  # End timing
        print(f"Elapsed time for location {loc}, jet {jet}: {elapsed_time:.2f} seconds")
        

# Save results
np.save(f'jet_GHG_separated{name}.npy', jet_CIs)

#%% compute weighted average CI for each jet producer
# Compute the weighted CIs
weighted_CIs = jet_CIs * sent_ethanol_no_uncertainty

# Sum the weighted CI and weights
sum_weighted_CIs = np.sum(weighted_CIs, axis=2)
sum_weights = np.sum(sent_ethanol_no_uncertainty, axis=1)

# Compute the weighted average CIs
# Avoid division by zero by using np.where to handle cases where sum_weights is zero
jet_CIs_average = np.where(sum_weights != 0, sum_weighted_CIs / sum_weights, 0)

np.save(f'jet_CI_ave{name}.npy', jet_CIs_average)

jet_CIs_stats = calculate_stats_uncertainty(jet_CIs_average)

df = pd.DataFrame(jet_CIs_stats)

# Sort the DataFrame by the 'mean' column
df_sorted = df.sort_values(by='mean').reset_index()

df_CIs = df_sorted

#%% CI vs Costs Pareto Front

# Sort the dataframes by the 'index' column before plotting
df_price_sorted = df_price.sort_values(by='index')
df_CIs_sorted = df_CIs.sort_values(by='index')


# Convert to numpy arrays for efficient processing
costs = df_price_sorted['mean'].to_numpy()
carbon_intensity = df_CIs_sorted['mean'].to_numpy()

# Initialize a list to store Pareto frontier points
pareto_front = []

# Sort by costs (ascending)
sorted_indices = np.argsort(costs)
sorted_costs = costs[sorted_indices]
sorted_carbon_intensity = carbon_intensity[sorted_indices]

# Start with the first point (lowest cost), which is Pareto-efficient
pareto_front.append((sorted_costs[0], sorted_carbon_intensity[0]))

# Traverse the sorted points and find Pareto-efficient points
min_ci = sorted_carbon_intensity[0]  # The smallest carbon intensity seen so far
for i in range(1, len(sorted_costs)):
    if sorted_carbon_intensity[i] < min_ci:
        pareto_front.append((sorted_costs[i], sorted_carbon_intensity[i]))
        min_ci = sorted_carbon_intensity[i]

# Convert Pareto front to a numpy array for plotting
pareto_front = np.array(pareto_front)


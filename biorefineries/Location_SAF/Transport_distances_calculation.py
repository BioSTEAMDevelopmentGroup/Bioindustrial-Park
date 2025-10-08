#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 28 14:21:54 2025

@author: bianco3
"""

from Feedstock_Transport_Class import *


model1 = FeedstockTransportModel(feedstock = 'switchgrass', 
                                 num_points = 100,
                                 samples_for_uncertainty = 100,
                                 sizes_of_biorefineries = 80,
                                 blending_capacity_set=0.1
                                 )
model1.calculate_location_parameters()

model1.plot_candidate_locations()

results1 = model1.results

#%%

model2 = FeedstockTransportModel(feedstock = 'switchgrass', 
                                 num_points = 100,
                                 samples_for_uncertainty = 100,
                                 sizes_of_biorefineries = 80,
                                 blending_capacity_set=0.2
                                 )
model2.calculate_location_parameters()

model2.plot_candidate_locations()

results2 = model2.results

#%%

model1_1000 = FeedstockTransportModel(feedstock = 'switchgrass', 
                                 num_points = 1000,
                                 samples_for_uncertainty = 100,
                                 sizes_of_biorefineries = 80,
                                 blending_capacity_set=0.1
                                 )
model1_1000.calculate_location_parameters()

#model1_1000.plot_candidate_locations()

results1_1000 = model1_1000.results

#%%

model2_1000 = FeedstockTransportModel(feedstock = 'miscanthus', 
                                 num_points = 1000,
                                 samples_for_uncertainty = 100,
                                 sizes_of_biorefineries = 80,
                                 blending_capacity_set=0.1
                                 )
model2_1000.calculate_location_parameters()

#model2_1000.plot_candidate_locations()

results2_1000 = model2_1000.results

#%%
dist_feedstock_res2 = np.array(results2_1000["farm_to_bio_mean_distances"])


numerator = np.array(results2_1000["ethanol_unit_transp_cost_each_jet"])
denominator = np.array(results2_1000["ethanol_unit_transp_cost_each_jet"])

where_to_send_model2_1000 = np.divide(
    numerator, 
    denominator, 
    out=np.zeros_like(numerator, dtype=float),  # fill with 0 by default
    where=denominator != 0                      # only divide where denominator ≠ 0
)

distances_ethanol_model2 = np.array(results2_1000["ethanol_transport_distances"][:,:3])*where_to_send_model2_1000

sent_ethanol_model2 = np.array(results2_1000["sent_ethanol"])

weighted_avg = np.average(distances_ethanol_model2, axis=1, weights=sent_ethanol_model2)


trans_cost_biomass_km_model2 = np.array(results2_1000["trans_cost_biomass_km"])

trans_GHG_model2 =np.array(results2_1000["trans_GHG"])
#%%


def print_array_stats(data, name="Array"):
    """
    Prints mean, median, and selected percentiles (5, 25, 75, 95) for a NumPy array.
    
    Parameters:
        data (array-like): Input array.
        name (str): Optional name for the dataset to include in the output.
    """
    # Calculate stats
    mean_val = np.mean(data)
    percentiles = np.percentile(data, [5, 25, 50, 75, 95])
    p5, p25, median_val, p75, p95 = percentiles
    
    # Print results
    print(f"\nStatistics for {name}:")
    print("-" * (15 + len(name)))
    print(f"Mean:       {mean_val:.4f}")
    print(f"Median:     {median_val:.4f}")
    print(f"5th pct:    {p5:.4f}")
    print(f"25th pct:   {p25:.4f}")
    print(f"75th pct:   {p75:.4f}")
    print(f"95th pct:   {p95:.4f}")
    print("-" * (15 + len(name)))
    
#%%


import geopandas as gpd
import numpy as np
import pandas as pd

def sample_farms_nearest_refinery(farms_gdf, refineries_gdf, n_samples=10000, random_state=None):
    """
    Randomly sample farms and calculate distance to nearest refinery.
    
    Returns a DataFrame with farm index, geometry, nearest refinery index, and distance (m and km)
    """
    rng = np.random.default_rng(random_state)
    
    # Sample farm indices
    n_farms = len(farms_gdf)
    sampled_indices = rng.choice(n_farms, size=n_samples, replace=False)
    sampled_farms = farms_gdf.iloc[sampled_indices].copy()
    
    # Use GeoPandas sjoin_nearest to find nearest refinery
    nearest = gpd.sjoin_nearest(
        sampled_farms,
        refineries_gdf,
        how='left',
        distance_col='distance_m'
    )
    
    # Add distance in km
    nearest['distance_km'] = nearest['distance_m'] / 1000
    
    # Keep relevant columns
    result = nearest.reset_index()[[
        'index',            # sampled farm index
        'geometry',         # farm geometry
        'index_right',      # nearest refinery index
        'distance_m',
        'distance_km'
    ]]
    
    result.rename(columns={'index': 'farm_index', 'index_right': 'refinery_index'}, inplace=True)
    
    return result

# Example usage:
distances_df = sample_farms_nearest_refinery(model1_1000.sw_data, model1_1000.jet_producers, n_samples=10000, random_state=42)
print(distances_df.head())

#%%
# Example usage:
distances_df2 = sample_farms_nearest_refinery(model2_1000.mis_data, model2_1000.jet_producers, n_samples=10000, random_state=2)



#%%


    

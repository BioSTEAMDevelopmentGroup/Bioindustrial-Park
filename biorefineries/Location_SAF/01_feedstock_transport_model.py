#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 10:50:15 2024

@author: bianco3

Results from this model:
    - Feedstock delivered price
    - Feedstock delivered CI
    - Ethanol unit cost of transportation to each jet producer (if biorefinery supplies to more than one)
    - Ethanol unit CI of transportation to each jet producer
    - Sent ethanol flow to each jet producer
    
To run model without uncertainty (only baseline values):
    Line 478-479 (choose feedstock and number of candidate locations)
    
To run model with uncertainty:
    Line 918-924 (choose feedstock, number of candidate locations, and number of samples for uncertainty)
    This will save results as npy files that will be used in the Biorefinery and SAF models as inputs
    
Next model to run: 02_ethanol_refinery_model

"""

# Import packages
import pandas as pd
import geopandas as gpd
import numpy as np
from scipy.spatial.distance import cdist
import libpysal as ps
import matplotlib.pyplot as plt
from pyproj import Transformer
import math
import matplotlib as mpl
import os
import ast
from scipy.optimize import curve_fit
import chaospy as cp
from SALib.sample import latin
from SALib.analyze import sobol

#%% Import files we are going to use

# Get directory for input files
# Get folder path of this file
folder = os.path.dirname(__file__)
# Join the folder path with the folder name where GIS data is located
input_data_folder = os.path.join(folder, 'GIS_data')

map_shapefile_path = os.path.join(input_data_folder, 'US_map_for_boundary.shp')
map_USA = gpd.read_file(map_shapefile_path) # Read the shapefile into a GeoDataFrame
map_USA = map_USA.to_crs('EPSG:5070')

map_shapefile_path = os.path.join(input_data_folder, 'USA-rainfedStates.shp')
USA_rainfed = gpd.read_file(map_shapefile_path) # Read the shapefile into a GeoDataFrame
USA_rainfed = USA_rainfed.to_crs('EPSG:5070')


# Import Jet producers, ethanol producers, yield and GHG for each crop, transportation costs
# Use same crs for all, in this project I will use a crs able to give distances in meters 'EPSG:5070'

# Import jet fuel producers shp file with capacity in Jet_Mbpd (Thousand barrels per day)
file_path = os.path.join(input_data_folder, 'Jet-fuel-producers.shp')
# Read the shapefile into a GeoDataFrame
jet_producers_gpd = gpd.read_file(file_path)
# Transform to same crs
jet_producers_gpd = jet_producers_gpd.to_crs('EPSG:5070')
# distances are in meters

# Ethanol biorefineries
# Import ethanol plants shp file with capacity in Mmgal per year (Millions of gallons per year)
file_path = os.path.join(input_data_folder, 'Ethanol_Plant.shp')
# Read the shapefile into a GeoDataFrame
ethanol_plants = gpd.read_file(file_path)
ethanol_plants = ethanol_plants.to_crs('EPSG:5070') # change crs to one that has units in meters for distance calculations


# Import switchgrass yield as pandas dataframe (yield in dry Mg/ha) Source: Fan et al. 2024
# File includes CI (GHG emissions) in kg CO2 eq / Mg feedstock
file_path = os.path.join(input_data_folder, 'switchgrass_data.csv')
# Read csv
sw_pd = pd.read_csv(file_path)
# Transform pandas dataframe into geopandas
sw_data = gpd.GeoDataFrame(sw_pd, crs = "EPSG:4326", geometry = gpd.points_from_xy(x = sw_pd.long, y = sw_pd.lat)).to_crs('EPSG:5070')


# Import miscanthus yield as pandas dataframe (yield in dry Mg/ha) Source: Fan et al. 2024
# File includes CI (GHG emissions) in kg CO2 eq / Mg feedstock
file_path = os.path.join(input_data_folder, 'miscanthus_data.csv')
# Read csv
mis_pd = pd.read_csv(file_path)
# Transform pandas dataframe into geopandas
mis_data = gpd.GeoDataFrame(mis_pd, crs = "EPSG:4326", geometry = gpd.points_from_xy(x = mis_pd.long, y = mis_pd.lat)).to_crs('EPSG:5070')

# Read layer of transportation costs
file_path = os.path.join(input_data_folder, 'States_with_transport.shp')
# Read the shapefile into a GeoDataFrame
transport_costs_gpd = gpd.read_file(file_path)
# Rename columns
# transport_costs_gpd = transport_costs_gpd.rename(columns={'Cost_trans': 'Flat_bed', 'Cost_tra_1': 'Tanker_truck'})
transport_costs_gpd = transport_costs_gpd.rename(columns={'Tanker': 'Tanker_truck'})
# Flat bed costs are in $/km/Mg of feedstock
# Tanker truck costs are in $/km/Mg of ethanol
# Transform crs to match the project
transport_costs_gpd = transport_costs_gpd.to_crs('EPSG:5070')


# Tortuosity factors per crop (per possible location, calculated in separate code - using seed = 1)
# For switchgrass
file_path = os.path.join(input_data_folder, 'TF_per_poss_locations.xlsx')
TF_sw = pd.read_excel(file_path)
mean_TF_sw = np.array(TF_sw["Mean"])
# For miscanthus
# file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Tortuosity/TF_per_poss_locations_mis.xlsx'
file_path = os.path.join(input_data_folder, 'TF_per_poss_locations_mis.xlsx')
TF_mis = pd.read_excel(file_path)
mean_TF_mis = np.array(TF_mis["Mean"])

# Real distances for ethanol transport
# For switchgrass
file_path = os.path.join(input_data_folder, 'results_TF_ethanol_sw.xlsx')
TF_ethanol_sw = pd.read_excel(file_path)
# Parse the string values into lists
TF_ethanol_sw['Tortuosity Factors'] = TF_ethanol_sw['Tortuosity Factors'].apply(lambda x: eval(x))
# Calculate the average for each row
TF_ethanol_sw['Average Tortuosity'] = TF_ethanol_sw['Tortuosity Factors'].apply(np.mean)
mean_TF_ethanol_sw = np.array(TF_ethanol_sw['Average Tortuosity'])
real_dist_ethanol_sw = TF_ethanol_sw['Real distance'].apply(ast.literal_eval)

# For miscanthus
file_path = os.path.join(input_data_folder, 'results_TF_ethanol_mis_complete.xlsx')
TF_ethanol_mis = pd.read_excel(file_path)
# Parse the string values into lists
TF_ethanol_mis['Tortuosity Factors'] = TF_ethanol_mis['Tortuosity Factors'].apply(lambda x: eval(x))
# Calculate the average for each row
TF_ethanol_mis['Average Tortuosity'] = TF_ethanol_mis['Tortuosity Factors'].apply(np.mean)
mean_TF_ethanol_mis = np.array(TF_ethanol_mis['Average Tortuosity'])
real_dist_ethanol_mis = TF_ethanol_mis['Real distance'].apply(ast.literal_eval)

#%% Conversion factors

working_days = 350 # Same as used in cellulosic BioSTEAM refinery
gal_per_barrel = 42 # gal/barrel of jet fuel
ethanol_to_jet = 0.31554808068174084 # from BioSTEAM ATJ simulation gal of jet fuel/gal of ethanol # up to 0.56 can be found in literature 
liters_per_gal = 3.785 # L/gal
liters_per_m3 = 1000 # L/m3
density_of_ethanol = 747.58 # from BioSTEAM simulation (kg/m3)
kg_per_ton = 1000 # kg/Mg (Metric tonnes)

# Capacity of trucks (Assumed same capacity transporting miscanthus and switchgrass as bales)
truck_capacity = 20 # Mg of feedstock per trip
tanker_capacity = 26.88 # Mg of ethanol per trips

# Number of candidate locations for biorefinery:
num_points = 1000 # This number can be changed if we want to evaluate more or less possible locations

# Area of each farm in ha
area = 4*4 #km2 Resolution of yield data 
area = area * 100 # ha
area = area *0.2 # assume 20% of area is planted with feedstock (either switchgrass or miscanthus)

# Conversion of feedstock to ethanol
mis_to_ethanol_gallon = 111.47 #  from BioSTEAM, in Gal of ethanol per dry ton of feedstock
sw_to_ethanol_gallon = 99.45 #  from BioSTEAM, in Gal of ethanol per dry ton of feedstock

#%% Function definitions
def sort_by_min_distance(first_mat, sorted_mat, num_points):
    """ Function to sort the first matrix according to indexes of the second matrix.
    This is performed for every row in num_points rows, returning the sorted matrix.
    """
    sorted_mat = np.array([first_mat[i][sorted_mat[i]] for i in range(num_points)])
    return sorted_mat


def calculate_farm_totals(sorted_mat, farm_num, cost_or_GHG_sorted, num_points):
    """Function to calculte total costs or emissions for farms. 
    It takes the sorted matrix up to farm_num (quantity of farms that supply the biorefinery according to its size)
    and multiplies the sorted_mat (farm_size_sorted) by the corresponding costs or GHG matrix up to farm_num.
    And sums all the results, performing this operation for the num_points rows (for every potential candidate
    biorefinery location).
    """
    result = np.array([np.sum(sorted_mat[i][:farm_num[i]]*cost_or_GHG_sorted[i][:farm_num[i]]) for i in range(num_points)])
    return result

def calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_or_GHG, num_points):
    """Function to calculte total costs or emissions for feedstock transport. 
    It takes the sorted distance up to farm_num (quantity of farms that supply the biorefinery according to its size)
    and multiplies the dist_sorted by the size of the farms. And sums all the results, multiplying the result by the 
    unit cost or GHG emissions of transporting one ton one kilometer, performing this operation for the num_points rows 
    (for every potential candidate biorefinery location).
    """
    # km * ton * $/km-ton = $ total
    result = np.array([np.sum(dist_sorted[i][:farm_num[i]]*farm_sizes_sorted[i][:farm_num[i]])*trans_cost_or_GHG[i] for i in range(num_points)])
    return result

def calculate_ethanol_transport_totals(dist_sorted, jet_ref_num, jet_ref_sizes_sorted, jet_ref_count_mat, trans_cost_or_GHG, num_points):
    """Function to calculte total costs or emissions for biofuel transport. 
    It takes the sorted distance up to jet_ref_num (quantity of jet fuel refineries that are supplied by the biorefinery according to 
    its size and their capacity * blending percentage capacity), and multiplies the dist_sorted by the size of the jet fuel producers. 
    And sums all the results, multiplying the result by the unit cost or GHG emissions of transporting one ton one kilometer, 
    performing this operation for the num_points rows (for every potential candidate biorefinery location).
    """
    result = np.array([np.sum(dist_sorted[i][:jet_ref_num[i]]*jet_ref_sizes_sorted[i][:jet_ref_num[i]]*jet_ref_count_mat[i][:jet_ref_num[i]])*trans_cost_or_GHG[i] for i in range(num_points)])
    return result


#%% Refineries capacity definitions

sizes_of_biorefineries = 80 # In Million Gallons of ethanol per year
# Median size of existing ethanol facilities

blending_capacity_set = 0.2 # Maximum capacity that jet producers are willing to divert to SAF production from their current capacity
# This will determine the max amount of ethanol to be delivered to each jet producer according to their capacity. 

# Conversion to mass of ethanol in Metric tonnes to calculate transportation costs
mass_of_biorefineries = sizes_of_biorefineries * liters_per_gal # In Million liters of ethanol per year
mass_of_biorefineries = mass_of_biorefineries *1000000 / liters_per_m3 # In m3 ethanol per year 
mass_of_biorefineries = mass_of_biorefineries * density_of_ethanol / kg_per_ton # In tonnes of ethanol per year

# To consider capacity of jet fuel producers as maximum to deliver
Capacity_jet_producers = np.array(jet_producers_gpd.Jet_Mbpd) # Thousands of barrels of jet fuel per day
# Conversion of Jet fuel capacity to tons of ethanol per year
Capacity_jet_producers *= working_days # Capacity in thousand of barrels of jet fuel per year

Capacity_jet_producers *= gal_per_barrel*1000 # Capacity in gallons of jet fuel per year
Capacity_jet_producers /= ethanol_to_jet # Capacity in gallons of ethanol per year
Capacity_jet_producers *= liters_per_gal # Capacity in L of ethanol per year
Capacity_jet_producers /= liters_per_m3 # Capacity in m3 of ethanol per year
Capacity_jet_producers *= density_of_ethanol # Capacity in kg of ethanol per year
Capacity_jet_producers /= kg_per_ton # Capacity in tonnes of ethanol per year (to supply 100% target)

#%% Function to calculate all location parameters

def calculate_location_parameters(feedstock, num_points):
    if feedstock == 'switchgrass':
        crop_data = sw_data
        crop_to_ethanol_gallon = sw_to_ethanol_gallon
        
        Tortuosity_reshaped = mean_TF_sw[:, np.newaxis] # Needs to be (1000,1) to match the shape of dist for multiplication 
        perc_5th_TF = np.array(TF_sw["5th percentile"])
        perc_5th_TF =perc_5th_TF[:, np.newaxis]
        perc_95th_TF = np.array(TF_sw["95th percentile"])
        perc_95th_TF =perc_95th_TF[:, np.newaxis]
        
        mean_TF_ethanol = mean_TF_ethanol_sw
        
        real_dist_ethanol = real_dist_ethanol_sw
    
    elif feedstock == 'miscanthus':
        crop_data = mis_data
        crop_to_ethanol_gallon = mis_to_ethanol_gallon
        Tortuosity_reshaped = mean_TF_mis[:, np.newaxis]
        perc_5th_TF = np.array(TF_mis["5th percentile"])
        perc_5th_TF =perc_5th_TF[:, np.newaxis]
        perc_95th_TF = np.array(TF_mis["95th percentile"])
        perc_95th_TF =perc_95th_TF[:, np.newaxis]

        mean_TF_ethanol = mean_TF_ethanol_mis
        real_dist_ethanol = real_dist_ethanol_mis

    feedstock_TF_D = []
    for i in range(num_points):
        list = [perc_5th_TF[i], Tortuosity_reshaped[i], perc_95th_TF[i]]
        list_of_floats = [float(arr.item()) for arr in list]
        feedstock_TF_D.append(list_of_floats)

    # Extract yield (in dry Mg/ha) / 0.8 to get wet tons (assuming 20% moisture content)
    yield_crop = np.array(crop_data.y_data) /0.8
    # Extract GHG (in kgCO2e/dry Mg) * 0.8 to get wet tons (assuming 20% moisture content)
    GHG_crop = np.array(crop_data.Net_CI) * 0.8
    # Extract Breakeven price of feedstock at farm gate in $/dry Mg * 0.8 to get wet tons (assuming 20% moisture content)
    Price_crop = np.array(crop_data.Farm_price) * 0.8

    # Calculate Farm costs per ton of feedstock 
    inflation_rate = 0.0319
    years = 2023-2016
    
    crop_data['Adjusted_Price'] = (crop_data['Farm_price']*0.8) * (1 + inflation_rate) ** years # in $/ wet Mg


    size_farms = area * yield_crop * np.ones(num_points)[..., None] 
    # yield in wet tons/ha * area in ha = total wet tons per farm

    probability_crop = yield_crop/np.sum(yield_crop)

    # Demand of crop
    sizes_of_crop_demand = sizes_of_biorefineries * 1000000 / (crop_to_ethanol_gallon *0.8 ) # In wet tonnes of crop of demand


    np.random.seed(1) 
    # Choose locations random, with higher probability for the ones with higher yields
    chosen_ind = np.random.choice(np.arange(len(yield_crop)),num_points, replace = False, p = probability_crop)   
    # This gives me the filtered yields for only the 1000 indices that were chosen before
    # We are going to use these 1000 locations as possible locations for ethanol biorefineries
    mask = crop_data.index.isin(chosen_ind)
    possible_locations = crop_data[mask]
    
    # Add  transportation costs
    possible_locations = gpd.sjoin(possible_locations, transport_costs_gpd[['Flat_bed', 'Tanker_truck', 'GHG_truck','geometry']], how='left', predicate='within')
        
    # Extract transportation costs for biomass in $/km/Mg and transportation costs for ethanol in $/km/Mg of ethanol
    # I need $/km for uncertainty so I multiply for capacity assumed, and then divide for the capacity considered with uncertainty
    # Different for both crops only because the probabilities are different so we might have different 1000 locations as candidates
    trans_cost_biomass_km = (np.array(possible_locations.Flat_bed)).astype(float) * 20 # It was uploaded as string so we have to transform to float
    trans_cost_biomass = trans_cost_biomass_km/truck_capacity
    trans_cost_ethanol_km = (np.array(possible_locations.Tanker_truck)).astype(float) *26.88
    trans_cost_ethanol = trans_cost_ethanol_km/tanker_capacity
    
    # Extract transportation GHG emissions for trucks in kgCO2e/km. Same unit emissions are assumed for transporting ethanol vs switchgrass,
    # The difference is in the capacity of each
    # Also difference in crops but only because we might be using different candidate locations for biorefineries
    trans_GHG = np.array(possible_locations.GHG_truck) #kgCO2e/km
    trans_GHG_truck = trans_GHG/truck_capacity #kgCO2e/km/Mg
    trans_GHG_tanker = trans_GHG/tanker_capacity #kgCO2e/km/Mg
    
    # We are calculating here the distance matrix between possible locations of biorefineries and 
    # total farms (each point in the yield data) = dji
    dist_mat = possible_locations.geometry.apply(lambda g: crop_data.distance(g)) # distance in meters
    # Transform to numpy array to perform operations faster
    dist = np.array(dist_mat)/1000 * Tortuosity_reshaped # /1000 to get values in km
    
    # This part is to sort the farms from min distance to max in rows
    sorted_farms = np.argsort(dist, axis = 1) # same for both crops
    
    # We order the farm sizes in the same way we ordered the farms (by minimum distance indices)
    farm_sizes_sorted = sort_by_min_distance(size_farms, sorted_farms, num_points) 
    # Here we calculate the accumulated sum of the farm size (in tons produced)
    farm_sizes_cum = np.cumsum(farm_sizes_sorted, axis = 1) 
    # We rearrange the distance matrix to respect the indices of sorted farms
    dist_sorted = sort_by_min_distance(dist, sorted_farms, num_points)
    
    # This calculates the distance matrix between the possible locations of ethanol biorefineries
    dist_2 = np.array(possible_locations.geometry.apply(lambda g: jet_producers_gpd.distance(g)))/1000 # divided by 1000 to get values in km
    # This part is to sort the jet_producers from min distance to max in rows 
    # This gives me indices
    sorted_jet_ref = np.argsort(dist_2, axis = 1)
    # We rearrange the distance matrix to respect the indices of sorted jet_fuel_refineries
    dist_2_sorted = sort_by_min_distance(dist_2, sorted_jet_ref, num_points)
    
    # This code calculates the number of farms required to fill the biorefinery demand in tonnes of feedstock
    farm_num=np.array([np.sum(farm_sizes_cum[i]<sizes_of_crop_demand) for i in range(num_points)]) 

    # We calculate our metric for the transportation of feedstock from the farms to the ethanol biorefinery
    # The transportation is given in total $
    cost_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_biomass, num_points) 
    
    # GHG emissions from farm to bio (transportation)
    GHG_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_GHG_truck, num_points) 

    size_jet_ref = Capacity_jet_producers * blending_capacity_set * np.ones(num_points)[..., None] # Same for all feedstocks
    
    # We order the jet refinery sizes in the same way we ordered the farms (by minimum distance indices)
    jet_ref_sizes_sorted = sort_by_min_distance(size_jet_ref, sorted_jet_ref, num_points) # Same for all feedstocks
    
    # Here we calculate the accumulated sum of the jet refinery size (in tons that can be received)
    jet_ref_sizes_cum = np.cumsum(jet_ref_sizes_sorted, axis = 1) # Same for all feedstocks
    
    # This code calculates the number of jet refineries required to deliver all the biorefinery production
    jet_ref_num=np.array([np.sum(jet_ref_sizes_cum[i]<mass_of_biorefineries) for i in range(num_points)]) +1 
    
    # This code calculates the number of jet refineries required to deliver all the biorefinery production
    jet_ref_count = np.zeros(num_points)# Initialize an empty array to store the number of jet_refs needed for each row
    for i in range(num_points):# Iterate over each row
        for j in range(54):
            if j == 0:
                if jet_ref_sizes_cum[i][j]>= mass_of_biorefineries:
                    jet_ref_count[i] = mass_of_biorefineries/jet_ref_sizes_cum[i][j] 
                    break
            else:
                if jet_ref_sizes_cum[i][j]>= mass_of_biorefineries:
                    jet_ref_count[i] = j + (mass_of_biorefineries- jet_ref_sizes_cum[i][j-1])/(jet_ref_sizes_cum[i][j] - jet_ref_sizes_cum[i][j-1])
                    break
    jet_ref_count_mat = np.zeros((num_points, num_points))
    for p, val in enumerate(jet_ref_count):
        integer_part = int(val)
        decimal_part = val - integer_part
        jet_ref_count_mat[p, :integer_part] = 1 # Fill the row with 1s up to the integer part of the value
        jet_ref_count_mat[p, integer_part] = decimal_part # Fill the rest of the row with the decimal part
    
    
    
    farm_cost_all_loc = np.array([np.mean(np.array((crop_data.iloc[sorted_farms[i][:farm_num[i]]]).Adjusted_Price)) for i in range(num_points)])
    # mean cost of feedstock for all farms supplying to each possible location (1000,0)
    
    # do the same for CI 
    farm_GHG_all_loc = np.array([np.mean(np.array((crop_data.iloc[sorted_farms[i][:farm_num[i]]]).Net_CI*0.8)) for i in range(num_points)])
    # in kgCO2e/wet ton crop
    # mean GHG of feedstock for all farms supplying to each possible location (1000,0)
    
    # Already calculated GHG from farm to bio
    # GHG_transp_farm_to_bio.shape # in total kgCO2e
    
    # Already calculated costs from farm to bio
    # cost_transp_farm_to_bio.shape # in total $
    
    # Total kilograms of feedstock produced (for biosteam price has to be in $/kg and GHG in kgCO2e/kg)
    kg_crop = sizes_of_crop_demand *1000 
    
    # Feedstock delivered price (biosteam assumes price of feedstock includes transportation)
    # farm cost in $/ton / 1000 kg/ton + transport cost in total $ / total kg of crop
    feedstock_delivered_price = farm_cost_all_loc/1000 + cost_transp_farm_to_bio/kg_crop # in $/ wet kg
    
    # farm GHG in kgCO2e/ton crop / 1000 kg/ton + transport GHG in total kgCO2e /total kg crop
    feedstock_delivered_GHG = farm_GHG_all_loc/1000 + GHG_transp_farm_to_bio/kg_crop # in kgCO2e/ wet kg crop
    
    # sizes_jet_1: sizes of jet producers we deliver in ton ethanol /year (10% of their total capacity)
    sizes_jet_1 = np.zeros((num_points, np.max(jet_ref_num)))
    for i in range(num_points):
        sizes_jet_1[i][:jet_ref_num[i]] = jet_ref_sizes_sorted[i][:jet_ref_num[i]]
        
    # sent_ethanol: actual amount of ethanol sent to jet producers (last one might be less than what it can receive)
    sent_ethanol = sizes_jet_1.copy()
    for i in range(num_points):    
        sent_ethanol[i][jet_ref_num[i]-1] = sizes_jet_1[i][jet_ref_num[i]-1]* jet_ref_count_mat[i][jet_ref_num[i]-1]
    # jet_indexes: to know which jet producer we supply to
    jet_indexes = np.zeros((num_points, np.max(jet_ref_num)))
    for i in range(num_points):
        jet_indexes[i][:jet_ref_num[i]] = sorted_jet_ref[i][:math.ceil(jet_ref_count[i])]
    # Now we combine the jet_indexes with the sent_ethanol to get tuples with which jet producer is being supplied and quantity
    combined = []
    for i in range(jet_indexes.shape[0]):
        row_i = []
        for j in range(jet_indexes.shape[1]):
            row_i.append((jet_indexes[i, j], sent_ethanol[i, j]))
        combined.append(row_i)
    

    if sizes_of_biorefineries == 80 and blending_capacity_set == 0.2:
        real_dist_ethanol_array = np.zeros((num_points, 2))
        for i in range(num_points):
            for j in range(2):
                try:
                    real_dist_ethanol_array[i][j] = real_dist_ethanol[i][j]
                except:
                    real_dist_ethanol_array[i][j] = 0
        dist_2_sorted_real = real_dist_ethanol_array/1000
    
    else:
        Tortuosity_ethanol_reshaped = mean_TF_ethanol[:, np.newaxis] 
        dist_2_sorted_real = dist_2_sorted * Tortuosity_ethanol_reshaped


    # ethanol_transp_cost_each_jet: in total $/year to each jet producer
    ethanol_transp_cost_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
    # This is in total $/year to each jet producer
    for i in range(num_points):
        ethanol_transp_cost_each_jet[i][:jet_ref_num[i]] = dist_2_sorted_real[i][:jet_ref_num[i]]*jet_ref_sizes_sorted[i][:jet_ref_num[i]]*jet_ref_count_mat[i][:jet_ref_num[i]] *trans_cost_ethanol[i]
    
    
    # ethanol_unit_transp_cost_each_jet: in $/Mg ethanol for each jet producer
    ethanol_unit_transp_cost_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
    for i in range(num_points):    
        ethanol_unit_transp_cost_each_jet[i][:jet_ref_num[i]] = ethanol_transp_cost_each_jet[i][:jet_ref_num[i]]/ sent_ethanol[i][:jet_ref_num[i]] # in $/Mg ethanol
    
    
    # ethanol_transp_GHG_each_jet: in total kgCO2e/year to each jet producer
    ethanol_transp_GHG_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
    for i in range(num_points):
        ethanol_transp_GHG_each_jet[i][:jet_ref_num[i]] = dist_2_sorted_real[i][:jet_ref_num[i]]*jet_ref_sizes_sorted[i][:jet_ref_num[i]]*jet_ref_count_mat[i][:jet_ref_num[i]] *trans_GHG_tanker[i]
    # This is in total kgCO2e/year
    
    
    # ethanol_unit_transp_GHG_each_jet: in kgCO2e/Mg ethanol for each jet producer
    # Ethanol transport GHG to each SAF producer delivered
    ethanol_unit_transp_GHG_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
    for i in range(num_points):    
        ethanol_unit_transp_GHG_each_jet[i][:jet_ref_num[i]] = ethanol_transp_GHG_each_jet[i][:jet_ref_num[i]]/ sent_ethanol[i][:jet_ref_num[i]] # in kgCO2e/Mg ethanol

    return (feedstock_delivered_price, feedstock_delivered_GHG, ethanol_unit_transp_cost_each_jet, 
            ethanol_unit_transp_GHG_each_jet, sent_ethanol, crop_to_ethanol_gallon, feedstock_TF_D, 
            Price_crop, yield_crop, GHG_crop, trans_cost_biomass_km, trans_cost_ethanol_km,
           trans_GHG, possible_locations, jet_indexes)

#%% Run the model without uncertainty (for  baseline values and num_points number of locations)

# choose feedstock and number of candidate locations
feedstock ='miscanthus'
num_points = 1000

# Run without uncertainty
feedstock_delivered_price, feedstock_delivered_GHG, ethanol_unit_transp_cost_each_jet, \
ethanol_unit_transp_GHG_each_jet, sent_ethanol, crop_to_ethanol_gallon, feedstock_TF_D, \
Price_crop, yield_crop, GHG_crop, trans_cost_biomass_km, trans_cost_ethanol_km, \
trans_GHG, possible_locations, jet_indexes = calculate_location_parameters(feedstock = feedstock, num_points = num_points )

possible_locations = possible_locations.drop(columns=['index_right'])
possible_locations_with_state = gpd.sjoin(possible_locations, USA_rainfed, how="left", predicate='intersects')

States_list = possible_locations_with_state.NAME
States_list = States_list.astype('category')

jet_producers_states_supplied = []
for i in range(len(jet_indexes)):
    list_for_each_possible_loc = []
    for j in range(len(jet_indexes[0])):
        list_for_each_possible_loc.append(jet_producers_gpd.iloc[int(jet_indexes[i][j])]['State'])
    jet_producers_states_supplied.append(list_for_each_possible_loc)


if feedstock =='switchgrass':
    np.save(sent_ethanol, 'sent_ethanol_no_uncertainty.npy')
    np.save(jet_producers_states_supplied, 'jet_producers_states_supplied.npy')
    possible_locations_with_state.to_file("possible_locations_with_state.shp")
elif feedstock == 'miscanthus':
    np.save(sent_ethanol, 'sent_ethanol_no_uncertainty_mis.npy')
    possible_locations_with_state.to_file("possible_locations_with_state_mis.shp")
    np.save(jet_producers_states_supplied, 'jet_producers_states_supplied_mis.npy')
    

#%% Definition of function to calculate samples from distribution of uncertain parameters

# number of samples
N = 1000
# number of possible locations
num_points = 1000
def calculate_distributions(N, num_points, crop_to_ethanol_gallon,
                           feedstock_TF_D, Price_crop, yield_crop,
                           GHG_crop, trans_cost_biomass_km, 
                           trans_cost_ethanol_km, trans_GHG):
    np.random.seed(seed=1234)
    #================ONE VALUE FOR UNCERTAINTY=====================================
    # Define distributions using chaospy
    eth_conv_rate_D = cp.Triangle(0.85*crop_to_ethanol_gallon, crop_to_ethanol_gallon, 1.15*crop_to_ethanol_gallon)
    size_farms_D = cp.Triangle(80, 320, 480)
    truck_capacity_D = cp.Uniform(17.5, 22.5)
    tanker_capacity_D = cp.Uniform(23.89, 29.87)
    
    # Create a list of distribution objects
    distributions = [eth_conv_rate_D, size_farms_D, truck_capacity_D, tanker_capacity_D]
    
    # Generate LHS samples using SALib
    num_samples = N
    
    # Define the problem with distributions
    problem = {
        'num_vars': 4,
        'names': ['ethanol conversion rate', 'size of farms', 'truck capacity', 'tanker capacity'],
        'bounds': [[0, 1]*4] # All bounds are [0, 1] because LHS samples are uniformly distributed between 0 and 1.
    }
    
    lhs_samples = latin.sample(problem, num_samples)
    
    # Convert LHS samples to actual parameter values using chaospy
    samples = np.zeros((num_samples, len(distributions)))
    for i, dist in enumerate(distributions):
        samples[:, i] = dist.ppf(lhs_samples[:, i])
    
    # Convert to a DataFrame for easy manipulation
    df_samples = pd.DataFrame(samples, columns=problem['names'])

    #================1000 VALUES (ONE VALUE PER POSSIBLE LOCATION) FOR UNCERTAINTY=====================================
    # possible locations = num_vars
    # definition to generate samples for 1000 candidate locations # this number should be changed if possible locations are more than 1000
    #
    def samples_per_possible_location(list_of_distributions, N):
        # Define your problem (assuming you have 1000 variables)
        num_vars = len(list_of_distributions)
        num_samples = N
        problem = {
            'num_vars': num_vars,
            'names': ['x' + str(i) for i in range(num_vars)],
            'bounds': [(0, 1)] * num_vars
        }
        
        # Generate N samples using LHS
        lhs_samples = latin.sample(problem, num_samples)
        
        # Convert the LHS samples to match the triangular distributions
        samples_from_distributions = np.array([
            [list_of_distributions[i].ppf(lhs_samples[j, i]) for i in range(num_vars)]
            for j in range(num_samples)
        ])
        
        # Each row in samples_from_distributions corresponds to one sample, 
        # with each column representing a draw from the corresponding triangular distribution
        return samples_from_distributions
    
    TF_feedstock_D = []
    # Define distributions using chaospy
    for i in range(num_points):
        a = cp.Triangle(feedstock_TF_D[i][0], feedstock_TF_D[i][1], feedstock_TF_D[i][2])
        TF_feedstock_D.append(a)

    Price_crop_D = []
    # Define distributions using chaospy
    for i in range(1):
        a = cp.Uniform(0.9, 1.1)
        Price_crop_D.append(a)

    yield_crop_D = []
    # Define distributions using chaospy
    for i in range(1): 
        a = cp.Normal(mu=0, sigma=1)
        yield_crop_D.append(a)

    GHG_crop_D = []
    # Define distributions using chaospy
    for i in range(1): 
        a = cp.Normal(mu=0, sigma=1)
        GHG_crop_D.append(a)

    trans_cost_biomass_km_D = []
    # Define distributions using chaospy
    for i in range(num_points):
        a = cp.Triangle(trans_cost_biomass_km[i]*0.9, trans_cost_biomass_km[i], trans_cost_biomass_km[i]*1.1)
        trans_cost_biomass_km_D.append(a)

    trans_cost_ethanol_km_D = []
    # Define distributions using chaospy
    for i in range(num_points):
        a = cp.Triangle(trans_cost_ethanol_km[i]*0.9, trans_cost_ethanol_km[i], trans_cost_ethanol_km[i]*1.1)
        trans_cost_ethanol_km_D.append(a)

    trans_GHG_D = []
    # Define distributions using chaospy
    for i in range(num_points):
        a = cp.Uniform(trans_GHG[i]*0.85, trans_GHG[i]*1.15)
        trans_GHG_D.append(a)

    dist_factor_ethanol_D = []
    # Define distributions using chaospy
    for i in range(num_points):
        a = cp.Triangle(0.9, 1, 1.1)
        dist_factor_ethanol_D.append(a)
        
    samples_from_TF_D = samples_per_possible_location(TF_feedstock_D, N)
    samples_from_Price_crop_D = samples_per_possible_location(Price_crop_D, N)
    samples_from_yield_crop_D = samples_per_possible_location(yield_crop_D, N)
    samples_from_GHG_crop_D = samples_per_possible_location(GHG_crop_D, N)
    samples_from_trans_cost_biomass = samples_per_possible_location(trans_cost_biomass_km_D, N)
    samples_from_trans_cost_ethanol = samples_per_possible_location(trans_cost_ethanol_km_D, N)
    samples_from_trans_GHG_D = samples_per_possible_location(trans_GHG_D, N)
    samples_from_dist_factor_ethanol_D = samples_per_possible_location(dist_factor_ethanol_D, N)

    return(df_samples, samples_from_TF_D, samples_from_Price_crop_D, samples_from_yield_crop_D,
    samples_from_GHG_crop_D, samples_from_trans_cost_biomass, samples_from_trans_cost_ethanol, 
    samples_from_trans_GHG_D, samples_from_dist_factor_ethanol_D)

#%% Function to calculate uncertainty for feedstock transport model and ethanol transport model

# number of samples
N = 1000
# number of possible locations
num_points = 1000
# function to do uncertainty
def calculate_uncertainty(N, num_points, feedstock):
    # with N being the number of samples
    # ===========function to calculate model 1st time, with chosen feedstock===============================
    # function calculate_location_parameters is outside this one
    feedstock_delivered_price, feedstock_delivered_GHG, ethanol_unit_transp_cost_each_jet, \
    ethanol_unit_transp_GHG_each_jet, sent_ethanol, crop_to_ethanol_gallon, feedstock_TF_D, \
    Price_crop, yield_crop, GHG_crop, trans_cost_biomass_km, trans_cost_ethanol_km, \
    trans_GHG = calculate_location_parameters(feedstock, num_points)
    # ==========functions to create distributions and create distributions====================
    # function calculate_distributions defined outside this one
    df_samples, samples_from_TF_D, samples_from_Price_crop_D, samples_from_yield_crop_D,\
    samples_from_GHG_crop_D, samples_from_trans_cost_biomass, samples_from_trans_cost_ethanol, \
    samples_from_trans_GHG_D, samples_from_dist_factor_ethanol_D = calculate_distributions(N = N, num_points = num_points,
                                                                                           crop_to_ethanol_gallon = crop_to_ethanol_gallon,
                                                                                           feedstock_TF_D = feedstock_TF_D, 
                                                                                           Price_crop=Price_crop, yield_crop=yield_crop,
                                                                                           GHG_crop=GHG_crop, 
                                                                                           trans_cost_biomass_km=trans_cost_biomass_km,
                                                                                           trans_cost_ethanol_km=trans_cost_ethanol_km,
                                                                                           trans_GHG=trans_GHG)
                                                                           
    #create lists of metrics
    feedstock_delivered_price_list =[]
    feedstock_delivered_GHG_list = []
    ethanol_unit_transport_cost_each_jet_list = []
    ethanol_unit_transport_GHG_each_jet_list = []
    sent_ethanol_list = []

    for sample in range(N): 
        
        # define samples for the variables that don't vary by possible location
        tanker_capacity = df_samples['tanker capacity'][sample]
        truck_capacity = df_samples['truck capacity'][sample]
        area = df_samples['size of farms'][sample]
        crop_to_ethanol_gallon = df_samples['ethanol conversion rate'][sample]
        Tortuosity_reshaped = samples_from_TF_D[sample][:, np.newaxis]
        
        if feedstock == 'switchgrass':
            crop_data = sw_data

            
            mean_TF_ethanol = mean_TF_ethanol_sw
            # no lo voy a variar en incertidumbre por ahora
            
            real_dist_ethanol = real_dist_ethanol_sw
            ## UNCERTAINTY - +/- 10% a este valor
        
        elif feedstock == 'miscanthus':
            crop_data = mis_data
            mean_TF_ethanol = mean_TF_ethanol_mis
            real_dist_ethanol = real_dist_ethanol_mis

               
        #  Yield in wet Mg/ha
        yield_crop = np.maximum(0, (np.array(crop_data.y_data) /0.8) + samples_from_yield_crop_D[sample]*3.5) # to transform from standard normal to N(Mean,3.5)
        # GHG (in kgCO2e/wet Mg) 
        GHG_crop = (np.array(crop_data.Net_CI) * 0.8) + samples_from_GHG_crop_D[sample]*78
        # Breakeven price of feedstock at farm gate in $/wet Mg 
        Price_crop = (np.array(crop_data.Farm_price) * 0.8) * samples_from_Price_crop_D[sample]
    
        # Calculate Farm costs per ton of feedstock 
        inflation_rate = 0.0319
        years = 2023-2016
        
        
        crop_data['Adjusted_Price'] = (crop_data['Farm_price']) * (1 + inflation_rate) ** years # in $/ wet Mg
    
    
        size_farms = area * yield_crop * np.ones(num_points)[..., None] 
        # yield in wet tons/ha * area in ha = total wet tons per farm
    
        probability_crop = yield_crop/np.sum(yield_crop)
    
        # Demand of crop
        sizes_of_crop_demand = sizes_of_biorefineries * 1000000 / (crop_to_ethanol_gallon *0.8 ) # In wet tonnes of crop of demand
    
    
        np.random.seed(1) 
        chosen_ind = np.random.choice(np.arange(len(yield_crop)),num_points, replace = False, p = probability_crop)   
        # This gives me the filtered yields for only the 1000 indices that were chosen before
        # We are going to use these 1000 locations as possible locations for ethanol biorefineries
        mask = crop_data.index.isin(chosen_ind)
        possible_locations = crop_data[mask] # In this case, possible_locations is the same for both crops (since the data is taken from the same points)
        
        # Add  transportation costs
        possible_locations = gpd.sjoin(possible_locations, transport_costs_gpd[['Flat_bed', 'Tanker_truck', 'GHG_truck','geometry']], how='left', predicate='within')
            

        trans_cost_biomass_km = samples_from_trans_cost_biomass[sample]
        trans_cost_biomass = trans_cost_biomass_km/truck_capacity # $/km/Mg 
        trans_cost_ethanol_km = samples_from_trans_cost_ethanol[sample]
        trans_cost_ethanol = trans_cost_ethanol_km/tanker_capacity # $/km/Mg ethanol
        
        # Extract transportation GHG emissions for trucks in kgCO2e/km. Same unit emissions are assumed for transporting ethanol vs switchgrass,
        # The difference is in the capacity of each
        # Also difference in crops but only because we might be using different candidate locations for biorefineries
        trans_GHG = samples_from_trans_GHG_D[sample] #kgCO2e/km
        trans_GHG_truck = trans_GHG/truck_capacity #kgCO2e/km/Mg
        trans_GHG_tanker = trans_GHG/tanker_capacity #kgCO2e/km/Mg
        
        # We are calculating here the distance matrix between possible locations of biorefineries and 
        # total farms (each point in the yield data) = dji
        dist_mat = possible_locations.geometry.apply(lambda g: crop_data.distance(g)) # distance in meters, same for both crops
        # Transform to numpy array to perform operations faster
        dist = np.array(dist_mat)/1000 * Tortuosity_reshaped # /1000 to get values in km
        
        # This part is to sort the farms from min distance to max in rows
        sorted_farms = np.argsort(dist, axis = 1) # same for both crops
        
        # We order the farm sizes in the same way we ordered the farms (by minimum distance indices)
        farm_sizes_sorted = sort_by_min_distance(size_farms, sorted_farms, num_points) 
        
        # Here we calculate the accumulated sum of the farm size (in tons produced)
        farm_sizes_cum = np.cumsum(farm_sizes_sorted, axis = 1) 
        # We rearrange the distance matrix to respect the indices of sorted farms
        dist_sorted = sort_by_min_distance(dist, sorted_farms, num_points)
        
        # This calculates the distance matrix between the possible locations of ethanol biorefineries
        dist_2 = np.array(possible_locations.geometry.apply(lambda g: jet_producers_gpd.distance(g)))/1000 # divided by 1000 to get values in km
        # This part is to sort the jet_producers from min distance to max in rows 
        # This gives me indices
        sorted_jet_ref = np.argsort(dist_2, axis = 1)
        # We rearrange the distance matrix to respect the indices of sorted jet_fuel_refineries
        dist_2_sorted = sort_by_min_distance(dist_2, sorted_jet_ref, num_points)
        
        # This code calculates the number of farms required to fill the biorefinery demand in tonnes of feedstock
        farm_num=np.array([np.sum(farm_sizes_cum[i]<sizes_of_crop_demand) for i in range(num_points)]) 
        
        # We calculate our metric for the transportation of feedstock from the farms to the ethanol biorefinery
        # The transportation is given in total $
        cost_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_biomass, num_points) 
        
        # GHG emissions from farm to bio (transportation)
        GHG_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_GHG_truck, num_points) 
    
        size_jet_ref = Capacity_jet_producers * blending_capacity_set * np.ones(num_points)[..., None] # Same for all feedstocks
        
        # We order the jet refinery sizes in the same way we ordered the farms (by minimum distance indices)
        jet_ref_sizes_sorted = sort_by_min_distance(size_jet_ref, sorted_jet_ref, num_points) # Same for all feedstocks
        
        # Here we calculate the accumulated sum of the jet refinery size (in tons that can be received)
        jet_ref_sizes_cum = np.cumsum(jet_ref_sizes_sorted, axis = 1) # Same for all feedstocks
        
        # This code calculates the number of jet refineries required to deliver all the biorefinery production
        jet_ref_num=np.array([np.sum(jet_ref_sizes_cum[i]<mass_of_biorefineries) for i in range(num_points)]) +1 
        
        # This code calculates the number of jet refineries required to deliver all the biorefinery production
        jet_ref_count = np.zeros(num_points)# Initialize an empty array to store the number of jet_refs needed for each row
        for i in range(num_points):# Iterate over each row
            for j in range(54):
                if j == 0:
                    if jet_ref_sizes_cum[i][j]>= mass_of_biorefineries:
                        jet_ref_count[i] = mass_of_biorefineries/jet_ref_sizes_cum[i][j] 
                        break
                else:
                    if jet_ref_sizes_cum[i][j]>= mass_of_biorefineries:
                        jet_ref_count[i] = j + (mass_of_biorefineries- jet_ref_sizes_cum[i][j-1])/(jet_ref_sizes_cum[i][j] - jet_ref_sizes_cum[i][j-1])
                        break
        jet_ref_count_mat = np.zeros((num_points, num_points))
        for p, val in enumerate(jet_ref_count):
            integer_part = int(val)
            decimal_part = val - integer_part
            jet_ref_count_mat[p, :integer_part] = 1 # Fill the row with 1s up to the integer part of the value
            jet_ref_count_mat[p, integer_part] = decimal_part # Fill the rest of the row with the decimal part
        
    
        farm_cost_all_loc = np.array([np.mean(np.array((crop_data.iloc[sorted_farms[i][:farm_num[i]]]).Adjusted_Price*0.8)) for i in range(num_points)])
        # mean cost of feedstock for all farms supplying to each possible location (1000,0)
        
        # do the same for CI 
        farm_GHG_all_loc = np.array([np.mean(np.array((crop_data.iloc[sorted_farms[i][:farm_num[i]]]).Net_CI*0.8)) for i in range(num_points)])
        # in kgCO2e/wet ton crop
        # mean GHG of feedstock for all farms supplying to each possible location (1000,0)
        
        # Already calculated GHG from farm to bio
        # GHG_transp_farm_to_bio.shape # in total kgCO2e
        
        # Already calculated costs from farm to bio
        # cost_transp_farm_to_bio.shape # in total $
        
        # Total kilograms of feedstock produced (for biosteam price has to be in $/kg and GHG in kgCO2e/kg)
        kg_crop = sizes_of_crop_demand *1000 
        
        # Feedstock delivered price (biosteam assumes price of feedstock includes transportation)
        # farm cost in $/ton / 1000 kg/ton + transport cost in total $ / total kg of crop
        feedstock_delivered_price = farm_cost_all_loc/1000 + cost_transp_farm_to_bio/kg_crop # in $/ wet kg
        
        # farm GHG in kgCO2e/ton crop / 1000 kg/ton + transport GHG in total kgCO2e /total kg crop
        feedstock_delivered_GHG = farm_GHG_all_loc/1000 + GHG_transp_farm_to_bio/kg_crop # in kgCO2e/ wet kg crop
        
        # sizes_jet_1: sizes of jet producers we deliver in ton ethanol /year (10% of their total capacity)
        sizes_jet_1 = np.zeros((num_points, np.max(jet_ref_num)))
        for i in range(num_points):
            sizes_jet_1[i][:jet_ref_num[i]] = jet_ref_sizes_sorted[i][:jet_ref_num[i]]
            
        # sent_ethanol: actual amount of ethanol sent to jet producers (last one might be less than what it can receive)
        sent_ethanol = sizes_jet_1.copy()
        for i in range(num_points):    
            sent_ethanol[i][jet_ref_num[i]-1] = sizes_jet_1[i][jet_ref_num[i]-1]* jet_ref_count_mat[i][jet_ref_num[i]-1]
        # jet_indexes: to know which jet producer we supply to
        jet_indexes = np.zeros((num_points, np.max(jet_ref_num)))
        for i in range(num_points):
            jet_indexes[i][:jet_ref_num[i]] = sorted_jet_ref[i][:math.ceil(jet_ref_count[i])]
        # Now we combine the jet_indexes with the sent_ethanol to get tuples with which jet producer is being supplied and quantity
        combined = []
        for i in range(jet_indexes.shape[0]):
            row_i = []
            for j in range(jet_indexes.shape[1]):
                row_i.append((jet_indexes[i, j], sent_ethanol[i, j]))
            combined.append(row_i)
    
        if sizes_of_biorefineries == 80 and blending_capacity_set == 0.2:
            real_dist_ethanol_array = np.zeros((num_points, 2))
            for i in range(num_points):
                for j in range(2):
                    try:
                        real_dist_ethanol_array[i][j] = real_dist_ethanol[i][j]
                    except:
                        real_dist_ethanol_array[i][j] = 0

            dist_2_sorted_real = real_dist_ethanol_array/1000 * samples_from_dist_factor_ethanol_D[sample][:, np.newaxis]
        
        else:
            Tortuosity_ethanol_reshaped = mean_TF_ethanol[:, np.newaxis] 
            
            dist_2_sorted_real = dist_2_sorted * Tortuosity_ethanol_reshaped * samples_from_dist_factor_ethanol_D[sample][:, np.newaxis]
    
    
        # ethanol_transp_cost_each_jet: in total $/year to each jet producer
        ethanol_transp_cost_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
        # This is in total $/year to each jet producer
        for i in range(num_points):
            ethanol_transp_cost_each_jet[i][:jet_ref_num[i]] = dist_2_sorted_real[i][:jet_ref_num[i]]*jet_ref_sizes_sorted[i][:jet_ref_num[i]]*jet_ref_count_mat[i][:jet_ref_num[i]] *trans_cost_ethanol[i]
        
        
        # ethanol_unit_transp_cost_each_jet: in $/Mg ethanol for each jet producer
        ethanol_unit_transp_cost_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
        for i in range(num_points):    
            ethanol_unit_transp_cost_each_jet[i][:jet_ref_num[i]] = ethanol_transp_cost_each_jet[i][:jet_ref_num[i]]/ sent_ethanol[i][:jet_ref_num[i]] # in $/Mg ethanol
        
        
        # ethanol_transp_GHG_each_jet: in total kgCO2e/year to each jet producer
        ethanol_transp_GHG_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
        for i in range(num_points):
            ethanol_transp_GHG_each_jet[i][:jet_ref_num[i]] = dist_2_sorted_real[i][:jet_ref_num[i]]*jet_ref_sizes_sorted[i][:jet_ref_num[i]]*jet_ref_count_mat[i][:jet_ref_num[i]] *trans_GHG_tanker[i]
        # This is in total kgCO2e/year
        
        
        # ethanol_unit_transp_GHG_each_jet: in kgCO2e/Mg ethanol for each jet producer
        # Ethanol transport GHG to each SAF producer delivered
        ethanol_unit_transp_GHG_each_jet = np.zeros((num_points, np.max(jet_ref_num)))
        for i in range(num_points):    
            ethanol_unit_transp_GHG_each_jet[i][:jet_ref_num[i]] = ethanol_transp_GHG_each_jet[i][:jet_ref_num[i]]/ sent_ethanol[i][:jet_ref_num[i]] # in kgCO2e/Mg ethanol
    
        feedstock_delivered_price_list.append(feedstock_delivered_price)
        feedstock_delivered_GHG_list.append(feedstock_delivered_GHG)
        ethanol_unit_transport_cost_each_jet_list.append(ethanol_unit_transp_cost_each_jet)
        ethanol_unit_transport_GHG_each_jet_list.append(ethanol_unit_transp_GHG_each_jet)
        sent_ethanol_list.append(sent_ethanol)
        

    return (feedstock_delivered_price_list, feedstock_delivered_GHG_list,
            ethanol_unit_transport_cost_each_jet_list, ethanol_unit_transport_GHG_each_jet_list,
            sent_ethanol_list)

#%% Run feedstock transport and ethanol transport models with uncertainty


# This part takes over 12 hours to run on a MacBook Pro with an Apple M2 Max processor, 64 GB RAM, macOS Sonoma 14.3.1.

# define feedstock
feedstock = 'miscanthus'

# define number of candidate locations
num_points = 1000

# Define number of samples for uncertainty
N=1000

# Wait 12 hours if num_points and N are = 1000
feedstock_delivered_price_list, feedstock_delivered_GHG_list,\
ethanol_unit_transp_cost_each_jet_list, ethanol_unit_transp_GHG_each_jet_list,\
sent_ethanol_list = calculate_uncertainty(N=N, num_points=num_points, feedstock=feedstock)

# Tranform to numpy arrays 
feedstock_delivered_price_array = np.array(feedstock_delivered_price_list)
feedstock_delivered_GHG_array = np.array(feedstock_delivered_GHG_list)
ethanol_unit_transp_cost_each_jet_array = np.array(ethanol_unit_transp_cost_each_jet_list)
ethanol_unit_transp_GHG_each_jet_array = np.array(ethanol_unit_transp_GHG_each_jet_list)
sent_ethanol_array = np.array(sent_ethanol_list)

# Save array to a .npy file
if feedstock == 'switchgrass':
    name = ''
else: 
    name = '_mis'
    
np.save(f'Feedstock_price{name}.npy', feedstock_delivered_price_array)
np.save(f'Feedstock_GHG{name}.npy', feedstock_delivered_GHG_array)
np.save(f'Ethanol_transport_cost{name}.npy', ethanol_unit_transp_cost_each_jet_array)
np.save(f'Ethanol_transport_GHG{name}.npy', ethanol_unit_transp_GHG_each_jet_array)
np.save(f'Sent_ethanol{name}.npy', sent_ethanol_array)

#%% Next model to run

# Run BioSTEAM file for biorefinery simulation '02_ethanol_refinery_model'









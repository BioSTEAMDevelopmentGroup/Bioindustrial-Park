#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 13:12:45 2024

@author: bianco3

Location model Switchgrass
Biorefinery of 80 million gal ethanol per year

To couple with BioSTEAM bioethanol refinery with delivered cost of feedstock and delivered GHG
(Includes cultivation, harvesting, and transportation)

Missing: delivered price of ethanol, and delivered GHG to couple with 
ATJ model (include transportaion in BioSTEAM's price and GWP of ethanol')

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
# import re

#%% Directory path for readable files and for results to be saved at

folder = os.path.dirname(__file__)
# st_data_file = os.path.join(folder, 'state_scenarios_for_import.xlsx')

input_data_folder = os.path.join(folder, 'input_data')

results_folder = os.path.join(folder, 'results')

#%% Read shapefiles and pandas data frames for inputs and maps
# TODO: save all files I need in same folder as this code (new folder called input_data), so I can use os instead of the complete path

# map_shapefile_path = os.path.join(input_data_folder, 'US_map_for_boundary.shp')
map_shapefile_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos descargados -varios/US states/US_map_for_boundary.shp' 
map_USA = gpd.read_file(map_shapefile_path) # Read the shapefile into a GeoDataFrame
map_USA = map_USA.to_crs('EPSG:5070')

#map_shapefile_path = os.path.join(input_data_folder, 'USA-rainfedStates.shp')
map_shapefile_path ='/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos descargados -varios/US states/USA-rainfedStates.shp' 
USA_rainfed = gpd.read_file(map_shapefile_path) # Read the shapefile into a GeoDataFrame
USA_rainfed = USA_rainfed.to_crs('EPSG:5070')

# Import all the shp and pd files we are going to use:
# Jet producers, ethanol producers, yield and GHG for each crop, transportation costs
# Use same crs for all, in this project I will use a crs able to give distances in meters 'EPSG:5070'

# Import jet fuel producers shp file with capacity in Jet_Mbpd (Thousand barrels per day)
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/GIS/Jet-fuel-producers.shp'
# file_path = os.path.join(input_data_folder, 'Jet-fuel-producers.shp')
# Read the shapefile into a GeoDataFrame
jet_producers_gpd = gpd.read_file(file_path)
# Transform to same crs
jet_producers_gpd = jet_producers_gpd.to_crs('EPSG:5070')
# distances are in meters

# Ethanol biorefineries
# Import ethanol plants shp file with capacity in Mmgal per year (Millions of gallons per year)
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos descargados -varios/Ethanol_Plants_US_EIA_-6086905423768615929/Ethanol_Plant.shp'
# file_path = os.path.join(input_data_folder, 'Ethanol_Plant.shp')
# Read the shapefile into a GeoDataFrame
ethanol_plants = gpd.read_file(file_path)
ethanol_plants = ethanol_plants.to_crs('EPSG:5070') # change crs to one that has units in meters for distance calculations


# Import switchgrass yield as pandas dataframe (yield in dry Mg/ha) Xinxin's data
# File includes CI (GHG emissions) in kg CO2 eq / Mg feedstock
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos GHG Xinxin/switchgrass_data.csv'
# file_path = os.path.join(input_data_folder, 'switchgrass_data.csv')
# Read csv
sw_pd = pd.read_csv(file_path)
# Transform pandas dataframe into geopandas
sw_data = gpd.GeoDataFrame(sw_pd, crs = "EPSG:4326", geometry = gpd.points_from_xy(x = sw_pd.long, y = sw_pd.lat)).to_crs('EPSG:5070')
# Extract yield (in dry Mg/ha)
# yield_sw = np.array(sw_data.y_data)
# # Extract GHG (in kgCO2e/Mg)
# GHG_sw = np.array(sw_data.Net_CI)
# # Extract Breakeven price of feedstock at farm gate in $/Mg
# Price_sw = np.array(sw_data.Farm_price)


# Import miscanthus yield as pandas dataframe (yield in dry Mg/ha) Xinxin's data
# File includes CI (GHG emissions) in kg CO2 eq / Mg feedstock
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos GHG Xinxin/miscanthus_data.csv'
# file_path = os.path.join(input_data_folder, 'miscanthus_data.csv')
# Read csv
mis_pd = pd.read_csv(file_path)
# Transform pandas dataframe into geopandas
mis_data = gpd.GeoDataFrame(mis_pd, crs = "EPSG:4326", geometry = gpd.points_from_xy(x = mis_pd.long, y = mis_pd.lat)).to_crs('EPSG:5070')
# # Extract yield (in dry Mg/ha)
# yield_mis = np.array(mis_data.y_data)
# # Extract GHG (in kgCO2e/Mg)
# GHG_mis = np.array(mis_data.Net_CI)
# # Extract Breakeven price of feedstock at farm gate in $/Mg
# Price_mis = np.array(mis_data.Farm_price)



# Read layer of transportation costs

file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Datos descargados -varios/US states/States_with_transport.shp'
# file_path = os.path.join(input_data_folder, 'States_with_transport.shp')
# Read the shapefile into a GeoDataFrame
transport_costs_gpd = gpd.read_file(file_path)
# Rename columns
# transport_costs_gpd = transport_costs_gpd.rename(columns={'Cost_trans': 'Flat_bed', 'Cost_tra_1': 'Tanker_truck'})
transport_costs_gpd = transport_costs_gpd.rename(columns={'Tanker': 'Tanker_truck'})
# Flat bed costs are in $/km/Mg of switchgrass
# Tanker truck costs are in $/km/Mg of ethanol
# Transform crs to match the project
transport_costs_gpd = transport_costs_gpd.to_crs('EPSG:5070')


# Tortuosity factors per crop (per possible location, calculated in separate code - using seed = 1)
# For switchgrass
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Tortuosity/TF_per_poss_locations.xlsx'
# file_path = os.path.join(input_data_folder, 'TF_per_poss_locations.xlsx')
TF_sw = pd.read_excel(file_path)
mean_TF_sw = np.array(TF_sw["Mean"])
# For miscanthus
## ADD here
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Tortuosity/TF_per_poss_locations_mis.xlsx'
# file_path = os.path.join(input_data_folder, 'TF_per_poss_locations_mis.xlsx')
TF_mis = pd.read_excel(file_path)
mean_TF_mis = np.array(TF_mis["Mean"])

# TODO: import here real distance from PR 
# For switchgrass
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Tortuosity/results_TF_ethanol_sw.xlsx'
# file_path = os.path.join(input_data_folder, 'results_TF_ethanol_sw.xlsx')
TF_ethanol_sw = pd.read_excel(file_path)
# Step 2: Parse the string values into lists
TF_ethanol_sw['Tortuosity Factors'] = TF_ethanol_sw['Tortuosity Factors'].apply(lambda x: eval(x))
# Step 3: Calculate the average for each row
TF_ethanol_sw['Average Tortuosity'] = TF_ethanol_sw['Tortuosity Factors'].apply(np.mean)
mean_TF_ethanol_sw = np.array(TF_ethanol_sw['Average Tortuosity'])
real_dist_ethanol_sw = TF_ethanol_sw['Real distance'].apply(ast.literal_eval)

# For miscanthus
file_path = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/Location model/SAF/Tortuosity/results_TF_ethanol_mis_complete.xlsx'
# file_path = os.path.join(input_data_folder, 'results_TF_ethanol_mis_complete.xlsx')
TF_ethanol_mis = pd.read_excel(file_path)
# Step 2: Parse the string values into lists
TF_ethanol_mis['Tortuosity Factors'] = TF_ethanol_mis['Tortuosity Factors'].apply(lambda x: eval(x))
# Step 3: Calculate the average for each row
TF_ethanol_mis['Average Tortuosity'] = TF_ethanol_mis['Tortuosity Factors'].apply(np.mean)
mean_TF_ethanol_mis = np.array(TF_ethanol_mis['Average Tortuosity'])
real_dist_ethanol_mis = TF_ethanol_mis['Real distance'].apply(ast.literal_eval)



#%% Conversion factors 
# (Make sure we are using the same for both feedstocks, except conversion to ethanol yield)
working_days = 350 # This number is the same as used in cellulosic biosteam refinery
gal_per_barrel = 42 # gal/barrel of jet fuel
ethanol_to_jet = 0.31554808068174084 # from BioSTEAM #0.56 from literature # gal of jet fuel/gal of ethanol
liters_per_gal = 3.785 # L/gal
liters_per_m3 = 1000 # L/m3
density_of_ethanol = 747.58 # from biosteam # 789 from internet # kg/m3
kg_per_ton = 1000 # kg/Mg (Metric tonnes)

# Capacity of trucks (Assumed same capacity transporting miscanthus and switchgrass as bales)
truck_capacity = 20 # tons of switchgrass per trip
tanker_capacity = 26.88 # tons of ethanol per trips

# Number of candidate locations for biorefinery:
num_points = 1000 # Check with more values to see if it slows a lot or not

# Area of each farm in ha
# Resolution of yield data will give me the area of each farm is 4 km2
area = 4*4 #km2
area = area * 100 # ha
area = area *0.2 # assume 20% of area is planted with feedstock (either switchgrass or miscanthus)

# CHECK WITH BIOSTEAM AND ADJUST NUMBER 
# 98 gallons of ethanol per U.S. ton of miscanthus to ethanol
mis_to_ethanol_gallon = 111.47 # # BIOSTEAM number in Gal of ethanol per dry ton
# sw_to_ethanol_gallon = 83.98 # gallons of ethanol per tons of switchgrass. NUMBER TAKEN FROM XINXIN'S SI 
sw_to_ethanol_gallon = 99.45 # 139.552145234607(SAF model) # BIOSTEAM number in Gal of ethanol per dry ton
# Xinxin's number is dry tons according to her

#%% CHANGE FEEDSTOCKS

feedstock = 'switchgrass' # change here to miscanthus

if feedstock == 'switchgrass':
    crop_data = sw_data
    crop_to_ethanol_gallon = sw_to_ethanol_gallon
    Tortuosity_reshaped = mean_TF_sw[:, np.newaxis] # Needs to be (1000,1) to match the shape of dist for multiplication 
    mean_TF_ethanol = mean_TF_ethanol_sw
    real_dist_ethanol = real_dist_ethanol_sw
elif feedstock == 'miscanthus':
    crop_data = mis_data
    crop_to_ethanol_gallon = mis_to_ethanol_gallon
    Tortuosity_reshaped = mean_TF_mis[:, np.newaxis]
    mean_TF_ethanol = mean_TF_ethanol_mis
    real_dist_ethanol = real_dist_ethanol_mis
    
# TODO: add here Tortuosity factors per crop

# Extract yield (in dry Mg/ha) / 0.8 to get wet tons (assuming 20% moisture content)
yield_crop = np.array(crop_data.y_data) /0.8
# Extract GHG (in kgCO2e/dry Mg) * 0.8 to get wet tons (assuming 20% moisture content)
GHG_crop = np.array(crop_data.Net_CI) * 0.8
# Extract Breakeven price of feedstock at farm gate in $/dry Mg * 0.8 to get wet tons (assuming 20% moisture content)
Price_crop = np.array(crop_data.Farm_price) * 0.8



#%% # FARM COSTS

# Calculate Farm costs per ton of feedstock 
inflation_rate = 0.0319
years = 2023-2016

Farm_cost_per_ton = Price_crop * (1+inflation_rate)**years # in $/ wet Mg (Price_crop was modified to wet tons)
# From Xinxin's data prices in 2016 dollars, discount rate 10%, actualized to 2022 by inflation

# Average inflation rate = 3.19% per year between 2016 and 2023 (price actualized to 2023 us dollars)
# Source of inflation rate: https://www.officialdata.org/us/inflation/2016?endYear=2022&amount=100
# Accessed on June 5, 2024


# We repeat the same for num_points number of rows, so each row represents a potential biorefinery to compare
Farm_cost_per_ton_mat = Farm_cost_per_ton * np.ones(num_points)[..., None]

crop_data['Adjusted_Price'] = (crop_data['Farm_price']) * (1 + inflation_rate) ** years # in $/ wet Mg

#%% # GHG EMISSIONS

# Array representing the GHG emissions of the field in kgCO2e/Mg of switchgrass (wet tons, since GHG_crop was modified to wet tons) 
# We repeat the same for num_points number of rows, so each row represents a potential biorefinery location to compare
GHG_crop_data = GHG_crop * np.ones(num_points)[..., None]


#%% # SIZE OF FARMS

# The size of the farms is going to be in tons produced and we repeat the 121784 farms in 1000 rows so 
# each row represents a potential biorefinery location to compare to.
size_farms = area * yield_crop * np.ones(num_points)[..., None] 
# yield in wet tons/ha * area in ha = total wet tons per farm


#%% # In this part, we divide the yield data by the total yield to get a number
# This number is used as a probability of choosing that location as a potential biorefinery candidate
# Locations with higher yield will have higher probability
probability_crop = yield_crop/np.sum(yield_crop)


#%% SIZE OF BIOREFINERY
# Modify if several sizes need to be included
# I was using a 200 Million Gallon ethanol per year refinery (mean size), but I switched to 80 MMgal (median size)

sizes_of_biorefineries = 80 # In Million Gallons of ethanol per year
# Median size of existing ethanol facilities
# The name is the same as the one used for many sizes

blending_capacity_set = 0.2 # I'm using 10% blending capacity at jet producers (Set if I want to consider many possibilities)
# INCREASED TO 20% TO SEE IF WE CAN DELIVER TO LESS JET PRODUCERS PER BIOREFINERY (22 JUL)
# This will determine the max amount of ethanol to be delivered to each jet
# producer according to their capacity. 
# TODO: check if it's betther to separate in 2 modules, feedstock and ethanol

# Conversion to mass of ethanol in Metric tonnes to calculate transportation costs
mass_of_biorefineries = sizes_of_biorefineries * liters_per_gal # In Million liters of ethanol per year
mass_of_biorefineries = mass_of_biorefineries *1000000 / liters_per_m3 # In m3 ethanol per year 
mass_of_biorefineries = mass_of_biorefineries * density_of_ethanol / kg_per_ton # In tonnes of ethanol per year

# To consider capacity of jet fuel producers as maximum to deliver
Capacity_jet_producers = np.array(jet_producers_gpd.Jet_Mbpd) # Thousands of barrels of jet fuel per day
# Conversion of Jet fuel capacity to tons of ethanol per year
Capacity_jet_producers *= working_days # Capacity in thousand of barrels of jet fuel per year
# TODO: see if I have to include an uptime

Capacity_jet_producers *= gal_per_barrel*1000 # Capacity in gallons of jet fuel per year
Capacity_jet_producers /= ethanol_to_jet # Capacity in gallons of ethanol per year
Capacity_jet_producers *= liters_per_gal # Capacity in L of ethanol per year
Capacity_jet_producers /= liters_per_m3 # Capacity in m3 of ethanol per year
Capacity_jet_producers *= density_of_ethanol # Capacity in kg of ethanol per year
Capacity_jet_producers /= kg_per_ton # Capacity in tonnes of ethanol per year (to supply 100% target)



# Demand of crop
sizes_of_crop_demand = sizes_of_biorefineries * 1000000 / (crop_to_ethanol_gallon *0.8 ) # In wet tonnes of crop of demand


#%% # MULTI-OBJECTIVE SCENARIO PARAMETER

# For multi-objective scenario, I define an alpha that represents the value of offsetting each ton of CO2eq
alpha = np.linspace(0,1,100) # Logaritmic scale? Not Based on VCM price because it favors GHG emissions, since they are larger
                       #(0.4,380,100)         # Based on Voluntary Carbon Market price and projections for limiting temperature rises
                                 # See Data-with-sources, sheet 'Carbon-price' for sources
# Still don't know if I will be doing this

#%% # Function definitions
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

#%% Run model part 1

np.random.seed(1)

# This line is for choosing a random set of 1000 indices from the yield array.
# In this run I will do only Switchgrass
# replace = False is used to ensure the chosen indices are unique. 
chosen_ind = np.random.choice(np.arange(len(yield_crop)),num_points, replace = False, p = probability_crop)


# This gives me the filtered yields for only the 1000 indices that were chosen before
# We are going to use these 1000 locations as possible locations for ethanol biorefineries
mask = crop_data.index.isin(chosen_ind)
possible_locations = crop_data[mask] # In this case, possible_locations is the same for both crops (since the data is taken from the same points)


# Add  transportation costs
possible_locations = gpd.sjoin(possible_locations, transport_costs_gpd[['Flat_bed', 'Tanker_truck', 'GHG_truck','geometry']], how='left', predicate='within')
    
# Extract transportation costs for biomass in $/km/Mg and transportation costs for ethanol in $/km/Mg of ethanol
# Different for both crops only because the probabilities are different so we might have different 1000 locations as candidates
trans_cost_biomass = (np.array(possible_locations.Flat_bed)).astype(float) # It was uploaded as string so we have to transform to float
trans_cost_ethanol = (np.array(possible_locations.Tanker_truck)).astype(float)
# TODO: check if trans_cost_ethanol should be in a different module

# Extract transportation GHG emissions for trucks in kgCO2e/km. Same unit emissions are assumed for transporting ethanol vs switchgrass,
# The difference is in the capacity of each
# Also difference in crops but only because we might be using different candidate locations for biorefineries
trans_GHG_truck = np.array(possible_locations.GHG_truck)/truck_capacity
trans_GHG_tanker = np.array(possible_locations.GHG_truck)/tanker_capacity
# TODO: check if trans_GHG_tanker should be in a different module

# We are calculating here the distance matrix between possible locations of biorefineries and 
# total farms (each point in the yield data) = dji
# It can take a while
dist_mat = possible_locations.geometry.apply(lambda g: crop_data.distance(g)) # distance in meters, same for both crops
# Transform to numpy array to perform operations faster
dist = np.array(dist_mat)/1000 * Tortuosity_reshaped # /1000 to get values in km

# This part is to sort the farms from min distance to max in rows
sorted_farms = np.argsort(dist, axis = 1) # same for both crops

# We order the farm sizes in the same way we ordered the farms (by minimum distance indices)
farm_sizes_sorted = sort_by_min_distance(size_farms, sorted_farms, num_points) 


### We do the same for farm costs and we multiply by farm size so we get total costs per farm
Farm_cost_per_ton_crop_sorted = sort_by_min_distance(Farm_cost_per_ton_mat, sorted_farms, num_points) 


# We order the GHG emissions in the same way we order the farms (by minimum distance indices)
GHG_crop_sorted = sort_by_min_distance(GHG_crop_data, sorted_farms, num_points) 

# Here we calculate the accumulated sum of the farm size (in tons produced)
farm_sizes_cum = np.cumsum(farm_sizes_sorted, axis = 1) 

# We rearrange the distance matrix to respect the indices of sorted farms
dist_sorted = sort_by_min_distance(dist, sorted_farms, num_points)

# This calculates the distance matrix between the possible locations of ethanol biorefineries
# and the jet fuel producers
# Same for both feedstocks
dist_2 = np.array(possible_locations.geometry.apply(lambda g: jet_producers_gpd.distance(g)))/1000 # divided by 1000 to get values in km
# TODO: check to see if it is better to put in a different module
# This part is to sort the jet_producers from min distance to max in rows 
# This gives me indices
sorted_jet_ref = np.argsort(dist_2, axis = 1)
# We rearrange the distance matrix to respect the indices of sorted jet_fuel_refineries
dist_2_sorted = sort_by_min_distance(dist_2, sorted_jet_ref, num_points)

#%% # Run model part 2

# This code calculates the number of farms required to fill the biorefinery demand in tonnes of feedstock
farm_num=np.array([np.sum(farm_sizes_cum[i]<sizes_of_crop_demand) for i in range(num_points)]) 
# TODO: sizes_of_sw_demand[b] if we have many sizes of biorefineries


# Field GHG emissions can be calculated first
Field_GHG_emissions = calculate_farm_totals(farm_sizes_sorted, farm_num,  GHG_crop_sorted, num_points) 

### We can calculate Farm costs the same way ($/ton * ton/farm = $/farm, and sum for all the supplying farms)
Cost_farm = calculate_farm_totals(farm_sizes_sorted, farm_num, Farm_cost_per_ton_crop_sorted, num_points) 

# We calculate our metric for the transportation of feedstock from the farms to the ethanol biorefinery
# The transportation is given in total $
cost_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_biomass, num_points) 

# GHG emissions from farm to bio (transportation)
GHG_transp_farm_to_bio = calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_GHG_truck, num_points) 

#%% Run model part 3 - jet fuel producers (Check to see if a different module would be better)
# TODO: check if this should be included in a different module
# The size of the jet refineries is going to be in tons tons of ethanol per year (to supply target%) 
# and we repeat the 54 jet refineries in 1000 rows so 
# each row represents a potential biorefinery location to compare to.
size_jet_ref = Capacity_jet_producers * blending_capacity_set * np.ones(num_points)[..., None] # Same for all feedstocks
# TODO: blending_capacity_set[t] if many capacities are considered

# We order the jet refinery sizes in the same way we ordered the farms (by minimum distance indices)
jet_ref_sizes_sorted = sort_by_min_distance(size_jet_ref, sorted_jet_ref, num_points) # Same for all feedstocks

# Here we calculate the accumulated sum of the jet refinery size (in tons that can be received)
jet_ref_sizes_cum = np.cumsum(jet_ref_sizes_sorted, axis = 1) # Same for all feedstocks

# This code calculates the number of jet refineries required to deliver all the biorefinery production
jet_ref_num=np.array([np.sum(jet_ref_sizes_cum[i]<mass_of_biorefineries) for i in range(num_points)]) +1 
# TODO: mass_of_biorefineries[b] if many are included (same in the for loop below)

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

# We calculate our metric for the transportation of ethanol from the biorefinery to the jet fuel producers
# The transportation is given in tons*meters (demand*distance)
cost_transp_bio_to_jet_2 = calculate_ethanol_transport_totals(dist_2_sorted, jet_ref_num, jet_ref_sizes_sorted, jet_ref_count_mat, trans_cost_ethanol, num_points)

# Here we calculate Transportation GHG emissions for ethanol to the jet fuel producers
GHG_transp_bio_to_jet = calculate_ethanol_transport_totals(dist_2_sorted, jet_ref_num, jet_ref_sizes_sorted, jet_ref_count_mat, trans_GHG_tanker, num_points)

#%% Totals

# Total transportation costs
total_cost_transp_crop = cost_transp_farm_to_bio + cost_transp_bio_to_jet_2 

### Total costs including farm operations and transportation
total_cost_2_crop = total_cost_transp_crop + Cost_farm 

# Total transportation GHG emissions
total_transp_GHG_crop = GHG_transp_farm_to_bio + GHG_transp_bio_to_jet 

# Total GHG emissions = Transportation + Field -> Used to minimize
total_GHG_field_transp_crop = total_transp_GHG_crop + Field_GHG_emissions 

# Multi-objective
Total_crop = total_cost_2_crop + alpha[...,None]*total_GHG_field_transp_crop 

# No minimization is required so far, we are going to minimize MJSP and jet GWP final values
#%% INFO TO INPUT INTO BIOSTEAM

# THIS IS IF WE WANTED THE MINIMUM WHICH WE DONT  
# farms_min_cost = sw_data.iloc[sorted_farms[np.argmin(total_cost_2_sw)][:farm_num_sw[np.argmin(total_cost_2_sw)]]]

# Cost of farm for all 1000 locations
# *1.25 as operating margin of farms 25%?
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
kg_crop_hr = kg_crop/(working_days)/24 # flow in kg/hr (wet kg)


# Feedstock delivered price (biosteam assumes price of feedstock includes transportation)
# farm cost in $/ton / 1000 kg/ton + transport cost in total $ / total kg of crop
feedstock_delivered_price = farm_cost_all_loc/1000 + cost_transp_farm_to_bio/kg_crop # in $/ wet kg

# farm GHG in kgCO2e/ton crop / 1000 kg/ton + transport GHG in total kgCO2e /total kg crop
feedstock_delivered_GHG = farm_GHG_all_loc/1000 + GHG_transp_farm_to_bio/kg_crop # in kgCO2e/ wet kg crop

#### FALTA: VER ETANOL, PRECIO DE TRANSPORTE Y GHG. A CUANTAS REFINERIAS DE PETROLEO LE ENTREGA Y CUANTO A CADA UNA
### buscar en SAF-location-sw-and_mis-multiple_runs-07-Combined-cost-GHG
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
# Finally, we get the unique tuples to remove repeated values
# To get unique tuples, flatten the list of lists and then get unique tuples
combined_flattened = [item for sublist in combined for item in sublist]
unique_tuples = list(set(combined_flattened))
# And this is to sort them by jet_index
# Sort unique tuples by the first element
sorted_unique_tuples = sorted(unique_tuples, key=lambda x: x[0])

#%% transportation for ethanol

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


#%% STATES list

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


#%%

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#data = feedstock_delivered_price # in $/ wet kg
data = feedstock_delivered_GHG
data_2 = ethanol_unit_transp_cost_each_jet[ethanol_unit_transp_cost_each_jet != 0]


# Calculate statistics
min_value = np.min(data)
max_value = np.max(data)
median_value = np.median(data)
mean_value = np.mean(data)

# Print statistics
print('Feedstock delivered GHG stats')
print(f'Min: {min_value}')
print(f'Max: {max_value}')
print(f'Median: {median_value}')
print(f'Mean: {mean_value}')

# Calculate statistics
min_value_2 = np.min(data_2)
max_value_2 = np.max(data_2)
median_value_2 = np.median(data_2)
mean_value_2 = np.mean(data_2)

# Print statistics
print('Ethanol transport price stats (in $/Mg ethanol)')
print(f'Min: {min_value_2}')
print(f'Max: {max_value_2}')
print(f'Median: {median_value_2}')
print(f'Mean: {mean_value_2}')

# # Plot the distribution
# # plt.figure(figsize=(10, 6))
# # sns.histplot(data, kde=True)
# # plt.title('Distribution of Data')
# # plt.xlabel('Value')
# # plt.ylabel('Frequency')
# # plt.show()

# # boxplot
# plt.figure(figsize=(10, 6))
# plt.boxplot(data, vert=False)
# plt.title('Boxplot of Data')
# plt.xlabel('Value')
# plt.show()

#%% Filter for Illinois

States_list = np.array(States_list)
data = np.array(data)

# Filter the data for "Illinois"
illinois_data = data[States_list == "Illinois"]

# Calculate the average value
average_value = np.mean(illinois_data)

print("Average value for Illinois:", average_value)


# Filter the data for "Kansas"
kansas_data = data[States_list == "Kansas"]

# Calculate the average value
average_value = np.mean(kansas_data)

print("Average value for Kansas:", average_value)

# Filter the data for "Alabama"
alabama_data = data[States_list == "Alabama"]

# Calculate the average value
average_value = np.mean(alabama_data)

print("Average value for Alabama:", average_value)

# Filter the data for "Arkansas"
arkansas_data = data[States_list == "Arkansas"]

# Calculate the average value
average_value = np.mean(arkansas_data)

print("Average value for Arkansas:", average_value)

#%% For ethanol delivered price stats

ethanol_delivered_price_each_jet = np.load('ethanol_delivered_price_each_jet.npy') # in $/gal ethanol taken from code file "Ethanol_price_per_location"

sent_ethanol_no_uncertainty = np.load('sent_ethanol_no_uncertainty.npy')


# Compute the weighted prices
weighted_prices = ethanol_delivered_price_each_jet * sent_ethanol_no_uncertainty

# Sum the weighted prices and weights
sum_weighted_prices = np.sum(weighted_prices, axis=2)
sum_weights = np.sum(sent_ethanol_no_uncertainty, axis=1)

# Compute the weighted average prices
# Avoid division by zero by using np.where to handle cases where sum_weights is zero
ethanol_prices_average = np.where(sum_weights != 0, sum_weighted_prices / sum_weights, 0)

mean_value_ethanol_price = np.mean(ethanol_prices_average)
print("Average ethanol delivered price for all states is:", mean_value_ethanol_price, " in USD per Gal")

arkansas_ethanol_data = ethanol_prices_average[States_list == "Arkansas"]
mean_value_ethanol_price_Arkansas = np.mean(arkansas_ethanol_data)
print("Average ethanol delivered price for Arkansas is:", mean_value_ethanol_price_Arkansas, " in USD per Gal")



#%% USEFUL FUNCTIONS

# def get_file_name(name):
#     return os.path.join(results_folder, name)

# model.table.to_excel(get_file_name('IP.xlsx'))



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 10:46:05 2025

@author: bianco3

Feedstock Transport model class

This class can be used to get feedstock delivered price and CI, and ethanol unit transport price and CI
for a biorefinery supplying ethanol to existing petroleum refineries (up to blending_capacity_set (20% default) of their capacity)
to produce SAF using switchgrass or miscanthus as feedstocks.

Simplifications are made to calculate distances when the class runs for a number of candidate locations different than 1000 (num_points).

For uncertainty, it is not advised to go over N=1000 samples because the computational time can be very long - (around 12 hours for 1000 candidate sites and 1000 uncertain scenarios)


"""

import os
import pandas as pd
import geopandas as gpd
import numpy as np
import ast
import math
import chaospy as cp
from SALib.sample import latin
import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point

class FeedstockTransportModel:
    
    # Class-level default conversion factors 
    default_conversion_factors = {
        'working_days': (350, 'days/year'),
        'gal_per_barrel': (42, 'gal/barrel'),
        'ethanol_to_jet': (0.31554808068174084, 'gal jet/gal ethanol'),
        'liters_per_gal': (3.785, 'L/gal'),
        'liters_per_m3': (1000, 'L/m³'),
        'density_of_ethanol': (747.58, 'kg/m³'),
        'kg_per_ton': (1000, 'kg/Mg'),
        'truck_capacity': (20, 'Mg/trip'),
        'tanker_capacity': (26.88, 'Mg/trip'),
        'farm_area_ha': (4 * 4 * 100 * 0.2, 'ha'),  # 4x4 km² → ha → 20% used
        'mis_to_ethanol_gallon': (111.47, 'gal ethanol/dry ton'),
        'sw_to_ethanol_gallon': (99.45, 'gal ethanol/dry ton'),
    }
    
    def __init__(self, feedstock = 'switchgrass', base_folder=None, num_points = 1000,
                 samples_for_uncertainty = 1000,
                 sizes_of_biorefineries=80,  # MMgal/year
                 blending_capacity_set=0.2,# Fraction (e.g. 0.2 = 20%)
                 bounding_box=None, # min_lon, min_lat, max_lon, max_lat) in EPSG:4326.
                 seed = 1, # used to choose candidate locations at random
                 refinery_locations = None): # (lat, lon) in EPSG:4326
        """
        Parameters
        ----------
        feedstock : str
            Feedstock should be either switchgrass or miscanthus. Default is switchgrass.
        base_folder : str or None
            Path to the folder where GIS and input data is stored. Defaults to a 'GIS_data' subfolder.
        num_points : int
            Number of candidate biorefinery locations or simulations. Default is 1000.
        sizes_of_biorefineries : float
            Size of biorefineries in MMgal/year (Million Gallons per year). Default is 80.
        blending_capacity_set : float
            Maximum share of petroleum refinery capacity that can be derived for SAF production. Default is 0.2.
        bounding_box: tuple or list of (min_lon, min_lat, max_lon, max_lat) in EPSG:4326.
            If None, all data is loaded. If provided, filters miscanthus and switchgrass data.
        seed: int 
            Seed used for the algorithm to choose candidate locations (with higher probability points with higher yields). Default = 1.
        refinery_locations: tuple or list of (lat, lon) in EPSG: 4326
            Candidate location to be analyzed. If provided, these locations will be analyzed as candidate refinery locations, num_points and seed will be ignored.
            Default is None, default candidate locations are chosen at random based on locations with higher yields (with seed) and number of points desired (num_points).
        """
        
        if base_folder is None:
            try:
                base_folder = os.path.join(os.path.dirname(__file__), 'GIS_data')
            except NameError:
                base_folder = os.path.join(os.getcwd(), 'GIS_data')
        
        self.base_folder = base_folder
        self._set_default_paths()
        
        self.feedstock = feedstock
        self.N = samples_for_uncertainty
        self.seed = seed
        self.refinery_locations = refinery_locations
        
        # Set num_points based on refinery_locations, if provided
        if self.refinery_locations is not None:
            self.num_points = len(self.refinery_locations)
        else:
            self.num_points = num_points
       
        
        self.conversion_factors = self.default_conversion_factors.copy() # Instance-level copy to allow individual modification
        
        self.sizes_of_biorefineries = sizes_of_biorefineries  # MM gal/year
        self.blending_capacity_set = blending_capacity_set
        self.bounding_box = bounding_box
        
        # Notify users if default values change
        if sizes_of_biorefineries != 80:
            print(f"[INFO] Biorefinery size set to {sizes_of_biorefineries} MMgal/year.")
        if blending_capacity_set != 0.2:
            print(f"[INFO] Blending capacity set to {blending_capacity_set*100:.1f}%.")
        
        self._load_all_data()
        
        self.results = {}  # store results from calculation for access later
        
    def _set_default_paths(self):
        """Set default paths to all files."""
        self.paths = {
            'map_shapefile': os.path.join(self.base_folder, 'US_map_for_boundary.shp'),
            'rainfed_shapefile': os.path.join(self.base_folder, 'USA-rainfedStates.shp'),
            'jet_producers': os.path.join(self.base_folder, 'Jet-fuel-producers.shp'),
            'ethanol_plants': os.path.join(self.base_folder, 'Ethanol_Plant.shp'),
            'switchgrass_data': os.path.join(self.base_folder, 'switchgrass_data.csv'),
            'miscanthus_data': os.path.join(self.base_folder, 'miscanthus_data.csv'),
            'transport_costs': os.path.join(self.base_folder, 'States_with_transport.shp'),
            'TF_sw': os.path.join(self.base_folder, 'TF_per_poss_locations.xlsx'),
            'TF_mis': os.path.join(self.base_folder, 'TF_per_poss_locations_mis.xlsx'),
            'TF_ethanol_sw': os.path.join(self.base_folder, 'results_TF_ethanol_sw.xlsx'),
            'TF_ethanol_mis': os.path.join(self.base_folder, 'results_TF_ethanol_mis_complete.xlsx'),
            'TF_per_state': os.path.join(self.base_folder, 'Tortuosity_per_state.xlsx'),
        }

    def update_path(self, key, new_path):
        """Update a specific file path (e.g., 'switchgrass_data', 'jet_producers')."""
        if key in self.paths:
            self.paths[key] = new_path
            print(f"Updated path for {key}.")
            self._load_all_data()
        else:
            raise KeyError(f"{key} is not a recognized data file key.")

    def _load_all_data(self):
        """Load all files into memory using geopandas/pandas."""
        self.map_USA = gpd.read_file(self.paths['map_shapefile']).to_crs('EPSG:5070')
        self.USA_rainfed = gpd.read_file(self.paths['rainfed_shapefile']).to_crs('EPSG:5070')
        self.jet_producers = gpd.read_file(self.paths['jet_producers']).to_crs('EPSG:5070')
        self.ethanol_plants = gpd.read_file(self.paths['ethanol_plants']).to_crs('EPSG:5070')

        # Optional spatial filter in EPSG:4326
        if self.bounding_box is not None:
            from shapely.geometry import box
            min_lon, min_lat, max_lon, max_lat = self.bounding_box
            bbox_poly = box(min_lon, min_lat, max_lon, max_lat)
            
        # # Switchgrass yield
        # sw_df = pd.read_csv(self.paths['switchgrass_data'])
        # self.sw_data = gpd.GeoDataFrame(
        #     sw_df, crs='EPSG:4326',
        #     geometry=gpd.points_from_xy(sw_df['long'], sw_df['lat'])
        # ).to_crs('EPSG:5070')
        
        # Switchgrass yield
        sw_df = pd.read_csv(self.paths['switchgrass_data'])
        sw_gdf = gpd.GeoDataFrame(
            sw_df, crs='EPSG:4326',
            geometry=gpd.points_from_xy(sw_df['long'], sw_df['lat'])
        )
    
        if self.bounding_box is not None:
            sw_gdf = sw_gdf[sw_gdf.geometry.within(bbox_poly)]
    
        self.sw_data = sw_gdf.to_crs('EPSG:5070')

        # # Miscanthus yield
        # mis_df = pd.read_csv(self.paths['miscanthus_data'])
        # self.mis_data = gpd.GeoDataFrame(
        #     mis_df, crs='EPSG:4326',
        #     geometry=gpd.points_from_xy(mis_df['long'], mis_df['lat'])
        # ).to_crs('EPSG:5070')
        
        # Miscanthus yield
        mis_df = pd.read_csv(self.paths['miscanthus_data'])
        mis_gdf = gpd.GeoDataFrame(
            mis_df, crs='EPSG:4326',
            geometry=gpd.points_from_xy(mis_df['long'], mis_df['lat'])
        )
    
        if self.bounding_box is not None:
            mis_gdf = mis_gdf[mis_gdf.geometry.within(bbox_poly)]
    
        self.mis_data = mis_gdf.to_crs('EPSG:5070')

        # Transport costs
        transport = gpd.read_file(self.paths['transport_costs']).to_crs('EPSG:5070')
        self.transport_costs = transport.rename(columns={'Tanker': 'Tanker_truck'})
        
        if self.num_points == 1000:
            # Tortuosity Factors
            self.mean_TF_sw = np.array(pd.read_excel(self.paths['TF_sw'])["Mean"])
            self.mean_TF_mis = np.array(pd.read_excel(self.paths['TF_mis'])["Mean"])
    
            # Ethanol transport distances and tortuosities
            self.real_dist_ethanol_sw, self.mean_TF_ethanol_sw = self._parse_ethanol_file(self.paths['TF_ethanol_sw'])
            self.real_dist_ethanol_mis, self.mean_TF_ethanol_mis = self._parse_ethanol_file(self.paths['TF_ethanol_mis'])
            self.calculate_eth_distance = False
        else:
            # Fill with default value but ensure correct shape
            default_tortuosity = 1.25
            self.mean_TF_sw = np.full((self.num_points,), default_tortuosity)
            self.mean_TF_mis = np.full((self.num_points,), default_tortuosity)
            self.mean_TF_ethanol_sw = np.full((self.num_points,), default_tortuosity)
            self.mean_TF_ethanol_mis = np.full((self.num_points,), default_tortuosity)
            self.calculate_eth_distance = True # calculate distances between points and jet producers
            self.TF_per_state = pd.read_excel(self.paths['TF_per_state'])


    def _parse_ethanol_file(self, path):
        df = pd.read_excel(path)
        df['Tortuosity Factors'] = df['Tortuosity Factors'].apply(lambda x: eval(x))
        df['Average Tortuosity'] = df['Tortuosity Factors'].apply(np.mean)
        distances = df['Real distance'].apply(ast.literal_eval)
        return distances, np.array(df['Average Tortuosity'])

        
    def print_conversion_factors(self):
        """Print all conversion factors with units."""
        print("Conversion Factors:")
        for key, (val, unit) in self.conversion_factors.items():
            print(f"  {key}: {val} {unit}")
            
    def set_conversion_factor(self, name, value, unit=None):
        """Update or add a conversion factor."""
        if unit is None and name in self.conversion_factors:
            unit = self.conversion_factors[name][1]  # Keep original unit
        elif unit is None:
            unit = ''
        self.conversion_factors[name] = (value, unit)
    
    
    ########## functions definition ########## They will be used for calculations inside the class
    @staticmethod
    def sort_by_min_distance(first_mat, sorted_mat, num_points):
        """Sort first_mat based on the indices in sorted_mat for each row."""
        return np.array([first_mat[i][sorted_mat[i]] for i in range(num_points)])

    @staticmethod
    def calculate_farm_totals(sorted_mat, farm_num, cost_or_GHG_sorted, num_points):
        """Calculate total cost or GHG for farms."""
        return np.array([
            np.sum(sorted_mat[i][:farm_num[i]] * cost_or_GHG_sorted[i][:farm_num[i]])
            for i in range(num_points)
        ])

    @staticmethod
    def calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_or_GHG, num_points):
        """Calculate total feedstock transport cost or GHG emissions."""
        return np.array([
            np.sum(dist_sorted[i][:farm_num[i]] * farm_sizes_sorted[i][:farm_num[i]]) * trans_cost_or_GHG[i]
            for i in range(num_points)
        ])

    @staticmethod
    def calculate_ethanol_transport_totals(dist_sorted, jet_ref_num, jet_ref_sizes_sorted, jet_ref_count_mat, trans_cost_or_GHG, num_points):
        """Calculate total ethanol transport cost or GHG emissions."""
        return np.array([
            np.sum(
                dist_sorted[i][:jet_ref_num[i]] *
                jet_ref_sizes_sorted[i][:jet_ref_num[i]] *
                jet_ref_count_mat[i][:jet_ref_num[i]]
            ) * trans_cost_or_GHG[i]
            for i in range(num_points)
        ])
    
    ########## end of functions definition ##########
    
    ######### functions to use inside the class but they are not static methods because they use instance attributes like self.sizes_of_biorefineries #########
    def calculate_mass_of_biorefineries(self, size_in_MG=None):
        if size_in_MG is None:
            size_in_MG = self.sizes_of_biorefineries
    
        cf = self.conversion_factors
        mass_ethanol = size_in_MG * cf['liters_per_gal'][0] * 1e6
        mass_ethanol /= cf['liters_per_m3'][0]
        mass_ethanol *= cf['density_of_ethanol'][0]
        mass_ethanol /= cf['kg_per_ton'][0]
        self.mass_of_biorefineries = mass_ethanol

    def calculate_capacity_of_jet_producers(self, jet_producers_gpd):
        if jet_producers_gpd is None:
            self.Capacity_jet_producers = None
            return
    
        cf = self.conversion_factors
        capacity = np.array(jet_producers_gpd.Jet_Mbpd)
        capacity *= cf['working_days'][0]
        capacity *= cf['gal_per_barrel'][0] * 1e3
        capacity /= cf['ethanol_to_jet'][0]
        capacity *= cf['liters_per_gal'][0] / cf['liters_per_m3'][0]
        capacity *= cf['density_of_ethanol'][0]
        capacity /= cf['kg_per_ton'][0]
        self.Capacity_jet_producers = capacity
        
    ######################################################
    
    def calculate_location_parameters(self):
        """
        If model = FeedstockTransportClass()
        Run model.results to see all results, or obtain each results by running model.get{result_name}() (e.g., model.get_feedstock_delivered_price())
        """
        num_points = self.num_points
        cf = self.conversion_factors
        
        if self.feedstock == 'switchgrass':
            crop_data = self.sw_data
        
            crop_to_ethanol_gallon = cf['sw_to_ethanol_gallon'][0]
            
            Tortuosity_reshaped = self.mean_TF_sw[:, np.newaxis] # Needs to be (1000,1) to match the shape of dist for multiplication 
            TF_sw = pd.read_excel(self.paths['TF_sw'])
            perc_5th_TF = np.array(TF_sw["5th percentile"])
            perc_5th_TF =perc_5th_TF[:, np.newaxis]
            perc_95th_TF = np.array(TF_sw["95th percentile"])
            perc_95th_TF =perc_95th_TF[:, np.newaxis]
            
            mean_TF_ethanol = self.mean_TF_ethanol_sw
            
            if self.calculate_eth_distance == False:
                real_dist_ethanol = self.real_dist_ethanol_sw
        
        elif self.feedstock == 'miscanthus':
            crop_data = self.mis_data
            crop_to_ethanol_gallon = cf['mis_to_ethanol_gallon'][0]
            Tortuosity_reshaped = self.mean_TF_mis[:, np.newaxis]
            TF_mis = pd.read_excel(self.paths['TF_mis'])
            perc_5th_TF = np.array(TF_mis["5th percentile"])
            perc_5th_TF =perc_5th_TF[:, np.newaxis]
            perc_95th_TF = np.array(TF_mis["95th percentile"])
            perc_95th_TF =perc_95th_TF[:, np.newaxis]

            mean_TF_ethanol = self.mean_TF_ethanol_mis
            
            if self.calculate_eth_distance == False:
                real_dist_ethanol = self.real_dist_ethanol_mis
            
        else:
            raise ValueError(f'Unsupported feedstock: {self.feedstock}')
        
        
        
        feedstock_TF_D = []
        for i in range(num_points):
            list_1 = [perc_5th_TF[i], Tortuosity_reshaped[i], perc_95th_TF[i]]
            list_of_floats = [float(arr.item()) for arr in list_1]
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

        area = cf['farm_area_ha'][0]
        size_farms = area * yield_crop * np.ones(num_points)[..., None] 
        
        # yield in wet tons/ha * area in ha = total wet tons per farm

        probability_crop = yield_crop/np.sum(yield_crop)

        # Demand of crop
        sizes_of_crop_demand = self.sizes_of_biorefineries * 1000000 / (crop_to_ethanol_gallon *0.8 ) # In wet tonnes of crop of demand

        if self.refinery_locations is not None:
            # Convert input list to GeoDataFrame in EPSG:4326
            gdf_input = gpd.GeoDataFrame(geometry=[Point(lon, lat) for lat, lon in self.refinery_locations], crs='EPSG:4326')
            
            gdf_input_5070 = gdf_input.to_crs('EPSG:5070')# Transform to EPSG:5070
            
            if crop_data.crs != 'EPSG:5070':
                crop_data = crop_data.to_crs('EPSG:5070') # Ensure crop_data is also in EPSG:5070
            
            from scipy.spatial import cKDTree # Find nearest crop_data point for each refinery location
        
            # Create arrays of coordinates
            crop_coords = np.array(list(zip(crop_data.geometry.x, crop_data.geometry.y)))
            refinery_coords = np.array(list(zip(gdf_input_5070.geometry.x, gdf_input_5070.geometry.y)))
            
            # Use KDTree for efficient nearest neighbor search
            tree = cKDTree(crop_coords)
            _, nearest_idx = tree.query(refinery_coords, k=1)
        
            # Get the rows from crop_data corresponding to nearest points
            possible_locations = crop_data.iloc[nearest_idx]
            self.num_points = len(possible_locations)
        else:
            # Original behavior: choose locations at random with yield-based probabilities
            np.random.seed(self.seed)
            chosen_ind = np.random.choice(len(yield_crop), self.num_points, replace=False, p=probability_crop)
            possible_locations = crop_data.iloc[chosen_ind]
            
        # np.random.seed(self.seed) 
        # # # Choose locations random, with higher probability for the ones with higher yields
        # # chosen_ind = np.random.choice(np.arange(len(yield_crop)),num_points, replace = False, p = probability_crop)  
        # # print('chosen_ind', chosen_ind)
        # # # This gives me the filtered yields for only the 1000 indices that were chosen before
        # # # We are going to use these 1000 locations as possible locations for ethanol biorefineries
        # # mask = crop_data.index.isin(chosen_ind)
        # # possible_locations = crop_data[mask]
        # # print('possible_locations', possible_locations)
        # chosen_ind = np.random.choice(len(yield_crop), num_points, replace=False, p=probability_crop)
        # possible_locations = crop_data.iloc[chosen_ind]

        
        # Add  transportation costs
        possible_locations = gpd.sjoin(possible_locations, self.transport_costs[['Flat_bed', 'Tanker_truck', 'GHG_truck','geometry']], how='left', predicate='within')
        
        possible_locations_copy = possible_locations.copy()# Make a copy before modification
        if 'index_right' in possible_locations_copy.columns:
            possible_locations_copy = possible_locations_copy.drop(columns=['index_right'])

        possible_locations_with_state = gpd.sjoin(possible_locations_copy, self.USA_rainfed, how="left", predicate='intersects') # Spatial join with state information

        self.candidate_locations = possible_locations
        
        if self.num_points != 1000:
            self.To_get_TF = possible_locations_with_state.merge(self.TF_per_state, left_on='NAME',right_on='State',  how='left')
            if self.feedstock == 'switchgrass':
                mean_TF = self.To_get_TF['Mean sw'].values
                perc_5th_TF = np.array(self.To_get_TF["5th percentile sw"])
                perc_5th_TF =perc_5th_TF[:, np.newaxis]
                perc_95th_TF = np.array(self.To_get_TF["95th percentile sw"])
                perc_95th_TF =perc_95th_TF[:, np.newaxis]
            elif self.feedstock == 'miscanthus':
                mean_TF = self.To_get_TF['Mean mis'].values
                perc_5th_TF = np.array(self.To_get_TF["5th percentile mis"])
                perc_5th_TF =perc_5th_TF[:, np.newaxis]
                perc_95th_TF = np.array(self.To_get_TF["95th percentile mis"])
                perc_95th_TF =perc_95th_TF[:, np.newaxis]
                
            Tortuosity_reshaped = mean_TF.reshape(-1, 1)
            
                
        # Extract transportation costs for biomass in $/km/Mg and transportation costs for ethanol in $/km/Mg of ethanol
        # I need $/km for uncertainty so I multiply for capacity assumed, and then divide for the capacity considered with uncertainty
        # Different for both crops only because the probabilities are different so we might have different 1000 locations as candidates
        truck_capacity = cf['truck_capacity'][0]
        tanker_capacity = cf['tanker_capacity'][0]
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
        farm_sizes_sorted = FeedstockTransportModel.sort_by_min_distance(size_farms, sorted_farms, num_points) 
        # Here we calculate the accumulated sum of the farm size (in tons produced)
        farm_sizes_cum = np.cumsum(farm_sizes_sorted, axis = 1) 
        # We rearrange the distance matrix to respect the indices of sorted farms
        dist_sorted = FeedstockTransportModel.sort_by_min_distance(dist, sorted_farms, num_points)
        
        # This calculates the distance matrix between the possible locations of ethanol biorefineries
        dist_2 = np.array(possible_locations.geometry.apply(lambda g: self.jet_producers.distance(g)))/1000 # divided by 1000 to get values in km
        # This part is to sort the jet_producers from min distance to max in rows 
        # This gives me indices
        sorted_jet_ref = np.argsort(dist_2, axis = 1)
        # We rearrange the distance matrix to respect the indices of sorted jet_fuel_refineries
        dist_2_sorted = FeedstockTransportModel.sort_by_min_distance(dist_2, sorted_jet_ref, num_points)
        
        # This code calculates the number of farms required to fill the biorefinery demand in tonnes of feedstock
        farm_num=np.array([np.sum(farm_sizes_cum[i]<sizes_of_crop_demand) for i in range(num_points)]) 

        # We calculate our metric for the transportation of feedstock from the farms to the ethanol biorefinery
        # The transportation is given in total $
        cost_transp_farm_to_bio = FeedstockTransportModel.calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_biomass, self.num_points) 
        
        # GHG emissions from farm to bio (transportation)
        GHG_transp_farm_to_bio = FeedstockTransportModel.calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_GHG_truck, self.num_points) 

        self.calculate_capacity_of_jet_producers(jet_producers_gpd = self.jet_producers)
        size_jet_ref = self.Capacity_jet_producers * self.blending_capacity_set * np.ones(num_points)[..., None] # Same for all feedstocks
        
        # We order the jet refinery sizes in the same way we ordered the farms (by minimum distance indices)
        jet_ref_sizes_sorted = FeedstockTransportModel.sort_by_min_distance(size_jet_ref, sorted_jet_ref, num_points) # Same for all feedstocks
        
        # Here we calculate the accumulated sum of the jet refinery size (in tons that can be received)
        jet_ref_sizes_cum = np.cumsum(jet_ref_sizes_sorted, axis = 1) # Same for all feedstocks
        
        self.calculate_mass_of_biorefineries()
        # This code calculates the number of jet refineries required to deliver all the biorefinery production
        jet_ref_num=np.array([np.sum(jet_ref_sizes_cum[i]<self.mass_of_biorefineries) for i in range(num_points)]) +1 
        
        # This code calculates the number of jet refineries required to deliver all the biorefinery production
        jet_ref_count = np.zeros(num_points)# Initialize an empty array to store the number of jet_refs needed for each row
        for i in range(num_points):# Iterate over each row
            for j in range(54):
                if j == 0:
                    if jet_ref_sizes_cum[i][j]>= self.mass_of_biorefineries:
                        jet_ref_count[i] = self.mass_of_biorefineries/jet_ref_sizes_cum[i][j] 
                        break
                else:
                    if jet_ref_sizes_cum[i][j]>= self.mass_of_biorefineries:
                        jet_ref_count[i] = j + (self.mass_of_biorefineries- jet_ref_sizes_cum[i][j-1])/(jet_ref_sizes_cum[i][j] - jet_ref_sizes_cum[i][j-1])
                        break
        # Get the maximum number of jet refineries needed (integer part + 1)
        max_jet_producers = int(np.ceil(np.max(jet_ref_count)))
        # Create matrix sized appropriately: rows = num_points, columns = max_jet_producers
        jet_ref_count_mat = np.zeros((num_points, max_jet_producers))
        for p, val in enumerate(jet_ref_count):
            integer_part = int(val)
            decimal_part = val - integer_part
            if integer_part < max_jet_producers:
                jet_ref_count_mat[p, :integer_part] = 1  # Fully filled refineries
                if integer_part < max_jet_producers:     # Add remaining fractional refinery
                    jet_ref_count_mat[p, integer_part] = decimal_part
            else:
                # Edge case: val is exactly equal to max_jet_producers
                jet_ref_count_mat[p, :max_jet_producers] = 1
        
        
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

        if self.sizes_of_biorefineries == 80 and self.blending_capacity_set == 0.2:
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
            
        if self.num_points != 1000:
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
    
        self.results['feedstock_delivered_price'] = feedstock_delivered_price
        self.results['feedstock_delivered_GHG'] = feedstock_delivered_GHG
        self.results['ethanol_unit_transp_cost_each_jet'] = ethanol_unit_transp_cost_each_jet
        self.results['ethanol_unit_transp_GHG_each_jet'] = ethanol_unit_transp_GHG_each_jet
        self.results['sent_ethanol'] = sent_ethanol
        self.results['crop_to_ethanol_gallon'] = crop_to_ethanol_gallon
        self.results['feedstock_TF_D'] = feedstock_TF_D
        self.results['Price_crop'] = Price_crop
        self.results['yield_crop'] = yield_crop
        self.results['GHG_crop'] = GHG_crop
        self.results['trans_cost_biomass_km'] = trans_cost_biomass_km
        self.results['trans_cost_ethanol_km'] = trans_cost_ethanol_km
        self.results['trans_GHG'] = trans_GHG
        self.results['possible_locations'] = possible_locations        
        self.results['jet_indexes'] = jet_indexes
        self.results['ethanol_transport_distances'] = dist_2_sorted_real

        
    # Functions to get the results from the above method (calculate_location_parameters)
    def get_feedstock_delivered_price(self):
        """
        Returns an array shape (num_points,) of the feedstock delivered price for each point (candidate location)
        in $/ wet kg (units used to input in BioSTEAM's refinery')
        """
        return self.results.get('feedstock_delivered_price', None)
    
    def get_feedstock_delivered_GHG(self):
        """
        Returns an array shape (num_points,) of the feedstock delivered CI for each point (candidate location)
        in kgCO2e/ wet kg crop (units used to input in BioSTEAM's refinery')
        """
        return self.results.get('feedstock_delivered_GHG', None)
    
    def get_ethanol_unit_transp_cost_each_jet(self):
        """
        Returns an array shape (num_points, number of jet producers supplied) of the ethanol unit transport cost to each petroleum
        refinery location being supplied (2 in the base case) for each candidate location (num_points)
        in $/Mg ethanol for each jet producer
        """
        return self.results.get('ethanol_unit_transp_cost_each_jet', None)

    def get_ethanol_unit_transp_GHG_each_jet(self):
        """
        Returns an array shape (num_points, number of jet producers supplied) of the ethanol unit transport CI to each petroleum
        refinery location being supplied (2 in the base case) for each candidate location (num_points)
        in kgCO2e/Mg ethanol for each jet producer
        """
        return self.results.get('ethanol_unit_transp_GHG_each_jet', None)

    def get_sent_ethanol(self):
        """
        Returns an array shape (num_points, number of jet producers supplied) of the actual amount of ethanol
        being supplied to each jet producer (20% of their capacity is the upper limit)
        in ton ethanol /year 
        """
        return self.results.get('sent_ethanol', None)

    def get_crop_to_ethanol_gallon(self):
        """
        Returns the value (float) of the conversion of crop to ethanol in gal ethanol/dry ton.
        It should be the same as the default conversion factor for the desired crop.

        """
        return self.results.get('crop_to_ethanol_gallon', None)
    
    def get_feedstock_TF_D(self):
        return self.results.get('feedstock_TF_D', None)
    
    def get_Price_crop(self):
        return self.results.get('Price_crop', None)
    
    def get_yield_crop(self):
        return self.results.get('yield_crop', None)
    
    def get_GHG_crop(self):
        return self.results.get('GHG_crop', None)
    
    def get_trans_cost_biomass_km(self):
        return self.results.get('trans_cost_biomass_km', None)
    
    def get_trans_cost_ethanol_km(self):
        return self.results.get('trans_cost_ethanol_km', None)
    
    def get_trans_GHG(self):
        return self.results.get('trans_GHG', None)
    
    def get_possible_locations(self):
        """
        Returns a gpd of candidate locations.
        """
        return self.results.get('possible_locations', None)
    
    def get_jet_indexes(self):
        """
        Returns the indexes of the jet producers being supplied by each candidate location.
        """
        return self.results.get('jet_indexes', None)
    
    ########## Function to calculate samples from distribution of uncertain parameters ########3
    def calculate_distributions(self,  crop_to_ethanol_gallon, 
                               feedstock_TF_D, Price_crop, yield_crop,
                               GHG_crop, trans_cost_biomass_km, 
                               trans_cost_ethanol_km, trans_GHG,
                               N , num_points):
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


    def calculate_uncertainty(self, N = None):
        if N == None:
            N = self.N
        feedstock = self.feedstock
        num_points = self.num_points
        
        
        # with N being the number of samples
        # ===========function to calculate model 1st time, with chosen feedstock===============================
        # function calculate_location_parameters is outside this one
        # feedstock_delivered_price, feedstock_delivered_GHG, ethanol_unit_transp_cost_each_jet, \
        # ethanol_unit_transp_GHG_each_jet, sent_ethanol, crop_to_ethanol_gallon, feedstock_TF_D, \
        # Price_crop, yield_crop, GHG_crop, trans_cost_biomass_km, trans_cost_ethanol_km, \
        # trans_GHG = self.calculate_location_parameters()
        self.calculate_location_parameters()
        feedstock_delivered_price = self.get_feedstock_delivered_price()
        feedstock_delivered_GHG = self.get_feedstock_delivered_GHG()
        ethanol_unit_transp_cost_each_jet = self.get_ethanol_unit_transp_cost_each_jet()
        ethanol_unit_transp_GHG_each_jet = self.get_ethanol_unit_transp_GHG_each_jet()
        sent_ethanol = self.get_sent_ethanol()
        crop_to_ethanol_gallon = self.get_crop_to_ethanol_gallon()
        feedstock_TF_D = self.get_feedstock_TF_D()
        Price_crop = self.get_Price_crop()
        yield_crop = self.get_yield_crop()
        GHG_crop = self.get_GHG_crop()
        trans_cost_biomass_km = self.get_trans_cost_biomass_km()
        trans_cost_ethanol_km = self.get_trans_cost_ethanol_km()
        trans_GHG = self.get_trans_GHG()
        # ==========functions to create distributions and create distributions====================
        # function calculate_distributions defined outside this one
        df_samples, samples_from_TF_D, samples_from_Price_crop_D, samples_from_yield_crop_D,\
        samples_from_GHG_crop_D, samples_from_trans_cost_biomass, samples_from_trans_cost_ethanol, \
        samples_from_trans_GHG_D, samples_from_dist_factor_ethanol_D = self.calculate_distributions(N = N, num_points = num_points,
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
                crop_data = self.sw_data

                
                mean_TF_ethanol = self.mean_TF_ethanol_sw
                # no lo voy a variar en incertidumbre por ahora
                
                if self.num_points == 1000:
                    real_dist_ethanol = self.real_dist_ethanol_sw
                    ## UNCERTAINTY - +/- 10% a este valor
                else:
                    real_dist_ethanol = self.results.get('ethanol_transport_distances', None)
            
            elif feedstock == 'miscanthus':
                crop_data = self.mis_data
                mean_TF_ethanol = self.mean_TF_ethanol_mis
                
                if self.num_points == 1000:
                    real_dist_ethanol = self.real_dist_ethanol_mis
                else:
                    real_dist_ethanol = self.results.get('ethanol_transport_distances', None)
                   
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
            sizes_of_crop_demand = self.sizes_of_biorefineries * 1000000 / (crop_to_ethanol_gallon *0.8 ) # In wet tonnes of crop of demand
        
            if self.refinery_locations is not None:
                # Convert input list to GeoDataFrame in EPSG:4326
                gdf_input = gpd.GeoDataFrame(geometry=[Point(lon, lat) for lat, lon in self.refinery_locations], crs='EPSG:4326')
                
                gdf_input_5070 = gdf_input.to_crs('EPSG:5070')# Transform to EPSG:5070
                
                if crop_data.crs != 'EPSG:5070':
                    crop_data = crop_data.to_crs('EPSG:5070') # Ensure crop_data is also in EPSG:5070
                
                from scipy.spatial import cKDTree # Find nearest crop_data point for each refinery location
            
                # Create arrays of coordinates
                crop_coords = np.array(list(zip(crop_data.geometry.x, crop_data.geometry.y)))
                refinery_coords = np.array(list(zip(gdf_input_5070.geometry.x, gdf_input_5070.geometry.y)))
                
                # Use KDTree for efficient nearest neighbor search
                tree = cKDTree(crop_coords)
                _, nearest_idx = tree.query(refinery_coords, k=1)
            
                # Get the rows from crop_data corresponding to nearest points
                possible_locations = crop_data.iloc[nearest_idx]
                self.num_points = len(possible_locations)
            else:
                # Original behavior: choose locations at random with yield-based probabilities
                np.random.seed(self.seed)
                chosen_ind = np.random.choice(len(yield_crop), self.num_points, replace=False, p=probability_crop)
                possible_locations = crop_data.iloc[chosen_ind]
                
            # np.random.seed(self.seed) 
            # chosen_ind = np.random.choice(len(yield_crop), num_points, replace=False, p=probability_crop)
            # possible_locations = crop_data.iloc[chosen_ind]
            
            # Add  transportation costs
            possible_locations = gpd.sjoin(possible_locations, self.transport_costs[['Flat_bed', 'Tanker_truck', 'GHG_truck','geometry']], how='left', predicate='within')
                

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
            farm_sizes_sorted = FeedstockTransportModel.sort_by_min_distance(size_farms, sorted_farms, num_points) 
            
            # Here we calculate the accumulated sum of the farm size (in tons produced)
            farm_sizes_cum = np.cumsum(farm_sizes_sorted, axis = 1) 
            # We rearrange the distance matrix to respect the indices of sorted farms
            dist_sorted = FeedstockTransportModel.sort_by_min_distance(dist, sorted_farms, num_points)
            
            # This calculates the distance matrix between the possible locations of ethanol biorefineries
            dist_2 = np.array(possible_locations.geometry.apply(lambda g: self.jet_producers.distance(g)))/1000 # divided by 1000 to get values in km
            # This part is to sort the jet_producers from min distance to max in rows 
            # This gives me indices
            sorted_jet_ref = np.argsort(dist_2, axis = 1)
            # We rearrange the distance matrix to respect the indices of sorted jet_fuel_refineries
            dist_2_sorted = FeedstockTransportModel.sort_by_min_distance(dist_2, sorted_jet_ref, num_points)
            
            # This code calculates the number of farms required to fill the biorefinery demand in tonnes of feedstock
            farm_num=np.array([np.sum(farm_sizes_cum[i]<sizes_of_crop_demand) for i in range(num_points)]) 
            
            # We calculate our metric for the transportation of feedstock from the farms to the ethanol biorefinery
            # The transportation is given in total $
            cost_transp_farm_to_bio = FeedstockTransportModel.calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_cost_biomass, num_points) 
            
            # GHG emissions from farm to bio (transportation)
            GHG_transp_farm_to_bio = FeedstockTransportModel.calculate_feedstock_transport_totals(dist_sorted, farm_num, farm_sizes_sorted, trans_GHG_truck, num_points) 
        
            size_jet_ref = self.Capacity_jet_producers * self.blending_capacity_set * np.ones(num_points)[..., None] # Same for all feedstocks
            
            # We order the jet refinery sizes in the same way we ordered the farms (by minimum distance indices)
            jet_ref_sizes_sorted = FeedstockTransportModel.sort_by_min_distance(size_jet_ref, sorted_jet_ref, num_points) # Same for all feedstocks
            
            # Here we calculate the accumulated sum of the jet refinery size (in tons that can be received)
            jet_ref_sizes_cum = np.cumsum(jet_ref_sizes_sorted, axis = 1) # Same for all feedstocks
            
            # This code calculates the number of jet refineries required to deliver all the biorefinery production
            jet_ref_num=np.array([np.sum(jet_ref_sizes_cum[i]<self.mass_of_biorefineries) for i in range(num_points)]) +1 
            
            # This code calculates the number of jet refineries required to deliver all the biorefinery production
            jet_ref_count = np.zeros(num_points)# Initialize an empty array to store the number of jet_refs needed for each row
            for i in range(num_points):# Iterate over each row
                for j in range(54):
                    if j == 0:
                        if jet_ref_sizes_cum[i][j]>= self.mass_of_biorefineries:
                            jet_ref_count[i] = self.mass_of_biorefineries/jet_ref_sizes_cum[i][j] 
                            break
                    else:
                        if jet_ref_sizes_cum[i][j]>= self.mass_of_biorefineries:
                            jet_ref_count[i] = j + (self.mass_of_biorefineries- jet_ref_sizes_cum[i][j-1])/(jet_ref_sizes_cum[i][j] - jet_ref_sizes_cum[i][j-1])
                            break

            # Get the maximum number of jet refineries needed (integer part + 1)
            max_jet_producers = int(np.ceil(np.max(jet_ref_count)))
            # Create matrix sized appropriately: rows = num_points, columns = max_jet_producers
            jet_ref_count_mat = np.zeros((num_points, max_jet_producers))
            for p, val in enumerate(jet_ref_count):
                integer_part = int(val)
                decimal_part = val - integer_part
                if integer_part < max_jet_producers:
                    jet_ref_count_mat[p, :integer_part] = 1  # Fully filled refineries
                    if integer_part < max_jet_producers:     # Add remaining fractional refinery
                        jet_ref_count_mat[p, integer_part] = decimal_part
                else:
                    # Edge case: val is exactly equal to max_jet_producers
                    jet_ref_count_mat[p, :max_jet_producers] = 1
        
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
        
            if self.sizes_of_biorefineries == 80 and self.blending_capacity_set == 0.2:
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
            
        self.uncertain_feedstock_delivered_price = np.array(feedstock_delivered_price_list)
        self.uncertain_feedstock_delivered_GHG = np.array(feedstock_delivered_GHG_list)
        self.uncertain_ethanol_unit_transport_cost_each_jet = np.array(ethanol_unit_transport_cost_each_jet_list)
        self.uncertain_ethanol_unit_transport_GHG_each_jet = np.array(ethanol_unit_transport_GHG_each_jet_list)
        self.uncertain_sent_ethanol = np.array(sent_ethanol_list)
    
    def get_uncertain_feedstock_delivered_price(self):
        """
        Returns an array shape (num_samples, num_points) of the feedstock delivered price for each point (candidate location)
        with uncertainty 
        in $/ wet kg (units used to input in BioSTEAM's refinery')
        and saves the file as 'Feedstock_price.npy' if it is switchgrass and 'Feedstock_price_mis.npy' if it is miscanthus
        """
        # Save array to a .npy file
        if self.feedstock == 'switchgrass':
            name = ''
        else: 
            name = '_mis'
        
        np.save(f'Feedstock_price{name}.npy', self.uncertain_feedstock_delivered_price)
        
        return self.uncertain_feedstock_delivered_price
    
    def get_uncertain_feedstock_delivered_GHG(self):
        """
        Returns an array shape (num_samples, num_points) of the feedstock delivered price for each point (candidate location)
        with uncertainty 
        in kgCO2e/ wet kg crop 
        and saves the file as 'Feedstock_GHG.npy' if it is switchgrass and 'Feedstock_GHG_mis.npy' if it is miscanthus
        """
        # Save array to a .npy file
        if self.feedstock == 'switchgrass':
            name = ''
        else: 
            name = '_mis'
        
        np.save(f'Feedstock_GHG{name}.npy', self.uncertain_feedstock_delivered_GHG)
        
        return self.uncertain_feedstock_delivered_GHG
    
    def get_uncertain_ethanol_transport_cost(self):
        """
        Returns an array shape (num_samples, num_points, num jet producers supplied) of the ethanol transport unit cost for each point (candidate location) to each jet fuel producer
        with uncertainty 
        in $/Mg ethanol
        and saves the file as 'Ethanol_transport_cost.npy' if it is switchgrass and 'Ethanol_transport_cost_mis.npy' if it is miscanthus
        """
        # Save array to a .npy file
        if self.feedstock == 'switchgrass':
            name = ''
        else: 
            name = '_mis'
        
        np.save(f'Ethanol_transport_cost{name}.npy', self.uncertain_ethanol_unit_transport_cost_each_jet)
        
        return self.uncertain_ethanol_unit_transport_cost_each_jet
    
    def get_uncertain_ethanol_transport_GHG(self):
        """
        Returns an array shape (num_samples, num_points, num jet producers supplied) of the ethanol transport unit GHG for each point (candidate location) to each jet fuel producer
        with uncertainty 
        in kgCO2e/Mg ethanol
        and saves the file as 'Ethanol_transport_GHG.npy' if it is switchgrass and 'Ethanol_transport_GHG_mis.npy' if it is miscanthus
        """
        # Save array to a .npy file
        if self.feedstock == 'switchgrass':
            name = ''
        else: 
            name = '_mis'
        
        np.save(f'Ethanol_transport_cost{name}.npy', self.uncertain_ethanol_unit_transport_GHG_each_jet)
        
        return self.uncertain_ethanol_unit_transport_GHG_each_jet
    
    def get_uncertain_sent_ethanol(self):
        """
        Returns an array shape (num_samples, num_points, num jet producers supplied) of the amount of ethanol sent from each point (candidate location) to each jet fuel producer
        with uncertainty 
        in Mg ethanol/year
        and saves the file as 'Sent_ethanol.npy' if it is switchgrass and 'Sent_ethanol_mis.npy' if it is miscanthus
        """
        # Save array to a .npy file
        if self.feedstock == 'switchgrass':
            name = ''
        else: 
            name = '_mis'
        
        np.save(f'Sent_ethanol{name}.npy', self.uncertain_sent_ethanol)
        
        return self.uncertain_sent_ethanol
    

    def plot_candidate_locations(self, color_metric=None, cmap='viridis',
                                 vmax=None, vmin=None, markersize = 30,
                                 crop_map = False, save_fig = False):
        """
        Plot candidate locations over the USA rainfed boundary, color-coded by a specified metric.
    
        Parameters:
        -----------
        color_metric : str or None
            Column name in the associated GeoDataFrame (e.g., TF) to use for color-coding.
            If None, plots the candidate locations as uniform markers.
        cmap : str
            Colormap to use for point coloring.
        vmax : float, optional
            Maximum value for colormap scaling.
        vmin : float, optional
            Minimum value for colormap scaling.
        markersize : int, optional
            Size of marker for points. Default is 10.
        crop_map : bool, optional
            Whether to crop the plot to a bounding box. Default is False.
        save_fig : bool, optional
            If True, saves the figure in high resolution. Default is False.
        """
        # Make sure candidate_locations and USA_rainfed exist
        if not hasattr(self, 'candidate_locations') or not hasattr(self, 'USA_rainfed'):
            print("Missing required data: 'candidate_locations' or 'USA_rainfed'")
            return
        gdf = self.candidate_locations.copy().to_crs(epsg=4326)
        rainfed = self.USA_rainfed.copy().to_crs(epsg=4326)
        
        if crop_map:
            minx, miny, maxx, maxy = self.bounding_box
            gdf = gdf.cx[minx:maxx, miny:maxy]
            rainfed = rainfed.cx[minx:maxx, miny:maxy]
        
        # Define mapping from abbreviation to actual column name
        metric_map = {
        'TF': {'col': 'TF', 'unit': ''},
        'fdp': {'col': 'feedstock_delivered_price', 'unit': 'USD/wet kg'},
        'fdg': {'col': 'feedstock_delivered_GHG', 'unit': 'kg CO₂e/wet kg'},
        'tcb_km': {'col': 'trans_cost_biomass_km', 'unit': 'USD/km'},
        'tce_km': {'col': 'trans_cost_ethanol_km', 'unit': 'USD/km'},
        'transport_CI': {'col': 'trans_GHG', 'unit': 'kg CO₂e/km'}
    }
        
        if color_metric is not None:
            if color_metric == 'TF':
                if self.num_points == 1000:
                    if self.feedstock == 'switchgrass':
                        gdf['TF'] = self.mean_TF_sw
                    elif self.feedstock == 'miscanthus':
                        gdf['TF'] = self.mean_TF_mis
                elif hasattr(self, 'To_get_TF'):
                    if self.feedstock == 'switchgrass':
                        gdf['TF'] = self.To_get_TF['Mean sw'].values
                    elif self.feedstock == 'miscanthus':
                        gdf['TF'] = self.To_get_TF['Mean mis'].values
                else:
                    print("Missing GeoDataFrame 'To_get_TF' with the color metric")
                    return
        
            elif color_metric in metric_map and hasattr(self, 'results'):
                full_col = metric_map[color_metric]['col']
                if full_col in self.results:
                    gdf[full_col] = self.results[full_col]
                else:
                    print(f"Metric '{full_col}' not found in self.results.")
                    return
            else:
                print(f"Unrecognized or unsupported color metric: {color_metric}")
                return
    
        fig, ax = plt.subplots(figsize=(12, 8))
        rainfed.boundary.plot(ax=ax, edgecolor='black', linewidth=1, alpha = 0.4)
        
        col_key = metric_map.get(color_metric, {'col': color_metric})['col']
        unit = metric_map.get(color_metric, {'unit': ''})['unit']
        
        if color_metric is not None and col_key in gdf.columns:
            gdf.plot(ax=ax, column=col_key, cmap=cmap,
                     legend=True, markersize=markersize, vmin=vmin, vmax=vmax)
            col_name = metric_map.get(color_metric, {'col': color_metric})['col']
            ax.set_title(f"Candidate Locations Colored by {col_name} ({unit})", fontsize=14)
        else:
            gdf.plot(ax=ax, color='teal', markersize=markersize)
            ax.set_title("Candidate Locations", fontsize=14)

        
        ax.axis('off')
        if save_fig:
            if color_metric is not None and col_key in gdf.columns:
                plt.savefig(f"Candidate Locations Colored by {col_name} ({unit})", dpi= 1200, bbox_inches='tight', transparent = True)
            else:
                plt.savefig("Candidate Locations", dpi= 1200, bbox_inches='tight', transparent = True)
        plt.show()
        
    def plot_ethanol_supply(self, markersize=30, line_color='gray', line_alpha=0.5, figsize=(12, 8),
                            save_fig = False):
        """
        Plot candidate locations and the jet producers they supply with connecting lines.
    
        Parameters:
        -----------
        markersize : int
            Size of the candidate and jet producer markers.
        line_color : str
            Color of the connecting lines.
        line_alpha : float
            Transparency of the connecting lines.
        figsize : tuple
            Size of the figure.
        save_fig : bool, optional
            If True, saves the figure in high resolution. Default is False.
        """
    
        # Basic checks
        if not hasattr(self, 'candidate_locations') or not hasattr(self, 'jet_producers'):
            print("Missing required GeoDataFrames: 'candidate_locations' or 'jet_producers'")
            return
        if 'jet_indexes' not in self.results:
            print("Missing 'jet_indexes' in results.")
            return
    
        gdf_candidates = self.candidate_locations.copy()
        gdf_jets = self.jet_producers.copy()
        rainfed = self.USA_rainfed.copy()
        jet_indexes = self.results['jet_indexes']  # shape (n_candidates, 2)
        sent_array = self.results['sent_ethanol']  # shape (n_candidates, 2)
        
    
        fig, ax = plt.subplots(figsize=figsize)
    
        # Optional: Plot USA rainfed background
        if hasattr(self, 'USA_rainfed'):
            rainfed.boundary.plot(ax=ax, edgecolor='black', linewidth=1)
    
        # Plot candidate locations and jet producers
        gdf_candidates.plot(ax=ax, color='#00a996', markersize=markersize, label='Candidate Locations')
        gdf_jets.plot(ax=ax, color='tab:red', markersize=markersize, label='Jet Producers')
    
        # Normalize sent ethanol to control linewidth scaling
        max_sent = sent_array.max()
        min_width = 0.5
        max_width = 2
    
        for idx, candidate in enumerate(gdf_candidates.itertuples()):
            candidate_point = candidate.geometry
            connected_jet_ids = jet_indexes[idx]  # Should be array/list of 2 jet indexes
    
            for j, jet_id in enumerate(connected_jet_ids):
                sent_amount = sent_array[idx, j]
    
                if sent_amount > 0 and 0 <= jet_id < len(gdf_jets):
                    jet_point = gdf_jets.loc[jet_id].geometry
                    line = LineString([candidate_point, jet_point])
                    normalized_width = min_width + (sent_amount / max_sent) * (max_width - min_width)
    
                    gpd.GeoSeries([line]).plot(
                        ax=ax,
                        color=line_color,
                        alpha=line_alpha,
                        linewidth=normalized_width
                    )
    
        ax.set_title("Candidate Locations and Connections to Supplied Jet Producers", fontsize=14)
        ax.legend()
        ax.axis('off')
        if save_fig:
            plt.savefig('Ethanol_supply.png', dpi= 1200, bbox_inches='tight', transparent = True)
        plt.show()



#%% Example of usage

# First import class 
# from Feedstock_Transport_Class import FeedstockTransportModel

# Create an instance of the class with default file paths
# model = FeedstockTransportModel()


# Example to update file
# model.update_path(key = 'switchgrass_data', 'path/to/file')
# keys are:
# map_shapefile' - US map
# 'rainfed_shapefile' - US rainfed states
# 'jet_producers' - Jet fuel producers 
# 'ethanol_plants' - Ethanol plants
# 'switchgrass_data'
# 'miscanthus_data'
# 'transport_costs' - shape file 
# 'TF_sw' - Tortuosity factor switchgrass
# 'TF_mis' - Tortuosity factors miscanthus
# 'TF_ethanol_sw' - Tortuosity factors for ethanol transport - switchgrass
# 'TF_ethanol_mis' - Tortuosity factors for ethanol transport - miscanthus


# Example to print data
# print(model.sw_data.head())

# Example to check conversion factor and constant of the model
# model.default_conversion_factors
# model.conversion_factors
# or model.print_conversion_factors()

# Example to change a conversion factor/constant
# model.set_conversion_factor('working_days', 360)
# Warning: the set_conversion_factor() method allows to input a unit value but the model itself is not able to adapt to the new unit,
# please use the same units as the default conversion factors


## Example to run the model without uncertainty (for base values)
# model = FeedstockTransportModel()
# model.calculate_location_parameters()
## and different arrays of results can be obtained with their correspondent method
# feedstock_delivered_GHG = model.get_feedstock_delivered_GHG()
# feedstock_delivered_price = model.get_feedstock_delivered_price()
# ethanol_unit_transp_cost_each_jet = model.get_ethanol_unit_transp_cost_each_jet()
# ethanol_unit_transp_GHG_each_jet = model.get_ethanol_unit_transp_GHG_each_jet()
# sent_ethanol = model.get_sent_ethanol() # amount of ethanol sent by each candidate location to each jet producer
# candidate_locations = model.get_possible_locations()
# jet_producers_supplied = model.get_jet_indexes() # indexes of jet producers being supplied by each candidate location


# Example to run the model with uncertainty
# model = FeedstockTransportModel(samples_for_uncertainty = 3)
# model.calculate_uncertainty()
## get arrays with uncertainty
# uncertain_feedstock_delivered_price = model.get_uncertain_feedstock_delivered_price()
# uncertain_feedstock_delivered_GHG = model.get_uncertain_feedstock_delivered_GHG()
# uncertain_ethanol_transport_cost = model.get_uncertain_ethanol_transport_cost()
# uncertain_ethanol_transport_GHG = model.get_uncertain_ethanol_transport_GHG()
# uncertain_sent_ethanol = model.get_uncertain_sent_ethanol()

# Example to plot candidate locations with metrics
# model.plot_candidate_locations(color_metric = 'fdp') # for feedstock delivered price

# Example to plot ethanol supplied
# model1.plot_ethanol_supply()

# Example to use not in all rainfed region but only in one state (e.g. = Illinois)
# Define bounding box for state
# bbox_illinois = (-91.5131, 36.9703, -87.4948, 42.5083)
# Run model with bounding box
# model3 = FeedstockTransportModel(num_points = 2, bounding_box = bbox_illinois)
# Calculate location parameters and uncertainty
# model3.calculate_location_parameters()
# model3.calculate_uncertainty()
# get results and visualize plots same as before
# model3.plot_candidate_locations()
# model3.plot_ethanol_supply()




import os
import numpy as np
import geopandas as gpd

# Load training inputs
folder_path = "/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/01_Location model/02_Machine Learning/01_Code_clean_ML/100_RBF_interpolation"

file_path = os.path.join(folder_path, "inputs_rand_500k_Rs_2perc.npy")
inputs_2perc = np.load(file_path, allow_pickle=True).T


# Load geospatial data
input_data_folder = '/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/01_Location model/02_Machine Learning/01_Code_clean_ML/00_Sample_generation/Input_data/GIS data'

file_name = 'USA-rainfedStates.shp'
file_path = os.path.join(input_data_folder, file_name)
USA_rainfed = gpd.read_file(file_path)
USA_rainfed = USA_rainfed.to_crs('EPSG:5070')# distances are in meters
# Great_Lakes
file_name = 'Great_Lakes/GL250515_lam.shp'
file_path = os.path.join(input_data_folder, file_name)
GL = gpd.read_file(file_path) 
GL = GL.to_crs('EPSG:5070')

# USA without great lakes
USA_rainfed = gpd.overlay(USA_rainfed, GL, how="difference")

# Export variables for import
__all__ = ['USA_rainfed', 'inputs_2perc']

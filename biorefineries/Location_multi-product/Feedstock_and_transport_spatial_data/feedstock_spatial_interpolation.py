## This code generates raster maps of feedstock yields, feedstock costs, and transportation unit costs across the study area. 
## It uses spatial interpolation techniques to estimate values at unsampled locations based on the known values at sampled locations (farms). 
## Author: M.F. Bianco, 2026

#%% Import necessary libraries
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import geopandas as gpd
import os
from shapely.geometry import Point, LineString
from scipy.interpolate import griddata
import matplotlib as mpl

#%% Import paths
from project.paths import MISCANTHUS_DATA, US_RAINFED, GREAT_LAKES, TRANSPORT_COSTS

# Import miscanthus yield as pandas dataframe (yield in dry Mg/ha) Xinxin Fan et al. 2024 data
# File includes CI (GHG emissions) in kg CO2 eq / Mg feedstock
mis_pd = pd.read_csv(MISCANTHUS_DATA) # Read csv
mis_data = gpd.GeoDataFrame(mis_pd, crs = "EPSG:4326", geometry = gpd.points_from_xy(x = mis_pd.long, y = mis_pd.lat)).to_crs('EPSG:5070') # Transform pandas dataframe into geopandas

# Read layer of transportation costs
transport_costs_gpd = gpd.read_file(TRANSPORT_COSTS)
transport_costs_gpd = transport_costs_gpd.rename(columns={'Tanker': 'Tanker_truck'}) # Rename columns
# Flat bed costs are in $/km/Mg of miscanthus
# Tanker truck costs are in $/km/Mg of ethanol
transport_costs_gpd = transport_costs_gpd.to_crs('EPSG:5070') # Transform crs to match the project

# USA_rainfed
USA_rainfed = gpd.read_file(US_RAINFED) 
USA_rainfed = USA_rainfed.to_crs('EPSG:5070') # distances are in meters

# Great_Lakes
GL = gpd.read_file(GREAT_LAKES) 
GL = GL.to_crs('EPSG:5070')

# USA without great lakes
USA_rainfed = gpd.overlay(USA_rainfed, GL, how="difference")

#%% Actualize costs to 2023 dollars and convert yields to wet Mg/ha, CI to kgCO2e/wet Mg

feedstock = 'miscanthus' # change here to miscanthus
if feedstock == 'miscanthus':
    crop_data = mis_data

# Calculate Farm costs per ton of feedstock
inflation_rate = 0.0319
years = 2023-2016

# From Xinxin's data prices in 2016 dollars, discount rate 10%, actualized to 2023 by inflation
# Average inflation rate = 3.19% per year between 2016 and 2023 (price actualized to 2023 us dollars)
# Source of inflation rate: https://www.officialdata.org/us/inflation/2016?endYear=2022&amount=100
# Accessed on June 5, 2024

crop_data['Farm_price'] = (crop_data['Farm_price']*0.8) * (1 + inflation_rate) ** years # in $/ wet Mg

# 20% moisture content assumed
crop_data['y_data'] = (crop_data['y_data']/0.8)  # in wet Mg/ha
crop_data['Net_CI'] = (crop_data['Net_CI']*0.8)  # in kgCO2e/wet Mg

# Perform the spatial join
crop_data = gpd.sjoin(crop_data, transport_costs_gpd, how="left", predicate="within")
# Flat bed costs are in $/km/Mg of switchgrass
# Tanker truck costs are in $/km/Mg of ethanol


#%% Create grid for interpolation

xmin, ymin, xmax, ymax = USA_rainfed.total_bounds

# 4km resolution for each cell
grid_x = np.linspace(xmin, xmax, np.int32((xmax-xmin)/4000))
grid_y = np.linspace(ymin, ymax, np.int32((ymax-ymin)/4000))

#%% Interpolate farm prices
# Convert data to NumPy arrays
points = np.vstack((crop_data.geometry.x, crop_data.geometry.y)).T
values = crop_data.Farm_price * crop_data.y_data

# Query grid points
grid_points = np.array([(x, y) for x in grid_x for y in grid_y])
interpolated_prices = griddata(points, values, grid_points, method='nearest') # instead of linear method

# Reshape to grid
Farm_prices_interpolated_map = interpolated_prices.reshape(len(grid_x), len(grid_y))


# save file
from project.paths import OUTPUT
maps_folder = OUTPUT / "Interpolated_maps"
maps_folder.mkdir(parents=True, exist_ok=True)  # creates folder if it doesn't exist

file_name = 'Farm_prices_interpolated_map_mis.npy'
output_path = maps_folder / file_name

# Save the array
np.save(output_path, Farm_prices_interpolated_map)
print(f"Saved file to: {output_path}")

#%% Interpolate truck prices
# Convert data to NumPy arrays
points = np.vstack((crop_data.geometry.x, crop_data.geometry.y)).T
values = crop_data.Flat_bed * crop_data.y_data

# Query grid points
grid_points = np.array([(x, y) for x in grid_x for y in grid_y])
interpolated_prices = griddata(points, values, grid_points, method='nearest')

# Reshape to grid
Truck_prices_interpolated_map = interpolated_prices.reshape(len(grid_x), len(grid_y))

# save file
file_name = 'Truck_prices_interpolated_map.npy'
output_path = maps_folder / file_name
np.save(output_path, Truck_prices_interpolated_map)
print(f"Saved file to: {output_path}")

#%% Interpolate crop yields
# Convert data to NumPy arrays
points = np.vstack((crop_data.geometry.x, crop_data.geometry.y)).T
values = crop_data.y_data

# Query grid points
grid_points = np.array([(x, y) for x in grid_x for y in grid_y])
interpolated_yields = griddata(points, values, grid_points, method='nearest')

# Reshape to grid
crop_yields_interpolated_map = interpolated_yields.reshape(len(grid_x), len(grid_y))

file_name = 'Crop_yields_interpolated_map_mis.npy'
output_path = maps_folder / file_name
np.save(output_path, crop_yields_interpolated_map)
print(f"Saved file to: {output_path}")



#%% function to plot interpolated maps
def plot_interpolated_map(
    grid_x,
    grid_y,
    interpolated_map,
    USA_boundary_gdf,
    cmap="viridis",
    colorbar_label="",
    num_ticks=8,
    save_path=None
):
    """
    General function to plot an interpolated 2D map with US boundary and colorbar.

    Parameters
    ----------
    grid_x, grid_y : 2D arrays
        Grid coordinates for pcolormesh.
    interpolated_map : 2D array
        Values to plot.
    USA_boundary_gdf : GeoDataFrame
        US boundary for overlay.
    cmap : str
        Matplotlib colormap.
    colorbar_label : str
        Label for the colorbar.
    num_ticks : int
        Number of ticks on the colorbar.
    save_path : str or Path, optional
        File path to save the figure. If None, figure is only shown.
    """
    
    plt.rcParams['font.family'] = 'Arial'
    
    vmin = interpolated_map.min()
    vmax = interpolated_map.max()
    
    plt.figure(figsize=(10, 6))
    
    mesh = plt.pcolormesh(
        grid_x, grid_y,
        interpolated_map.T,
        shading='auto',
        cmap=cmap,
        vmin=vmin,
        vmax=vmax
    )
    
    cbar = plt.colorbar(mesh)
    tick_values = np.linspace(vmin, vmax, num_ticks)
    cbar.set_ticks(tick_values)
    cbar.set_ticklabels([f"{t:.1f}" for t in tick_values])
    cbar.set_label(colorbar_label, fontsize=16)
    cbar.ax.tick_params(labelsize=13)
    
    USA_boundary_gdf.boundary.plot(ax=plt.gca(), color='black', linewidth=0.5)
    
    plt.axis('off')
    plt.tight_layout()
    
    if save_path is not None:
        plt.savefig(save_path, dpi=800, transparent=True)
    
    plt.show()

# Example usage:

# Farm prices
# plot_interpolated_map(
#     grid_x, grid_y,
#     Farm_prices_interpolated_map,
#     USA_rainfed,
#     cmap="RdYlGn_r",
#     colorbar_label="Interpolated farm price [USD·ha$^{-1}$]",
#     num_ticks=7,
#     save_path=PROJECT_ROOT / "figures" / "Farm_price_interpolated.png"
# )

# # Truck prices
# plot_interpolated_map(
#     grid_x, grid_y,
#     Truck_prices_interpolated_map,
#     USA_rainfed,
#     cmap="viridis",
#     colorbar_label="Interpolated truck price [USD·km$^{-1}$·ha$^{-1}$]",
#     num_ticks=8,
#     save_path=PROJECT_ROOT / "figures" / "Truck_price_interpolated.png"
# )

# # Crop yields
# plot_interpolated_map(
#     grid_x, grid_y,
#     crop_yields_interpolated_map,
#     USA_rainfed,
#     cmap="viridis",
#     colorbar_label="Interpolated yields [Mg·ha$^{-1}$]",
#     num_ticks=8,
#     save_path=PROJECT_ROOT / "figures" / "Crop_yields_interpolated.png"
# )
## This file generates num_samples random samples of farms
## for N_rand number of refineries across the US rainfed region,
## between irregular areas formed by 4 radii (defined with a 
## minimum and maximum value according to the feedstock desity scenario)
## The outputs are numpy files (.npy) of x and y coordinates of the 
## farms and the refineries.
## Author: M.F. Bianco - 2026-01-30

#%% Import libraries
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import nearest_points
from scipy.spatial import cKDTree  # Fast lookup
from scipy.interpolate import RegularGridInterpolator
import os
import joblib 
import re


class RefineryFarmSampler:
    """Class to generate spatial samples of collection areas and farms around refineries based on irregular radii.
    
    Parameters:
    - scenario: str, one of "2perc", "5perc", "10perc" or custom (if custom, provide max_rad, and keep format "Xperc")
    - x_coords: np.array, x coordinates of points where refineries can be sampled from. For this project, we use the 
                suitable land for perennial grass cultivation at 4x4km resulution across the US rainfed region.
    - y_coords: np.array, y coordinates of points where refineries can be sampled from.
    - max_rad: float, maximum radius for sampling farms around refineries (in meters). 
                If scenario is one of the predefined ones, this is set automatically.
    - min_rad: float, minimum radius for sampling farms around refineries (in meters). Default is 1km.
    - n_interp_points: int, number of points to define the irregular area around each refinery, equally spaced in angle (default is 4).
    - n_samples_per_refinery: int, number of farm samples to generate for each refinery (default is 100).
    - N_rand: int, number of random refinery locations to sample (default is 500,000).
    - seed: int, random seed for reproducibility (default is 123).  
    """

    SCENARIO_RADII = {
        "2perc": 270_000,
        "5perc": 170_000,
        "10perc": 120_000
    } # default radii in meters

    def __init__(self,
                 scenario,
                 x_coords = None,
                 y_coords = None,
                 max_rad=None,
                 min_rad=1_000,
                 n_interp_points=4,
                 n_samples_per_refinery=100,
                 N_rand=500_000,
                 seed=123):
        
        match = re.search(r"\d+perc", scenario)
        if match:
            self.scenario = match.group()
        else:
            self.scenario = scenario

        if x_coords == None or y_coords == None:
            # Read suitable land for perennial grass cultivation coordinates
            from project.paths import SUITABLE_LAND_X, SUITABLE_LAND_Y
            self.x_coords = np.load(SUITABLE_LAND_X)
            self.y_coords = np.load(SUITABLE_LAND_Y)
        else:
            self.x_coords = x_coords
            self.y_coords = y_coords

        self.min_rad = min_rad
        self.n_interp_points = n_interp_points
        self.n_samples_per_refinery = n_samples_per_refinery
        self.N_rand = N_rand
        self.seed = seed

        # Handle max radius
        if self.scenario in self.SCENARIO_RADII:
            self.max_rad = self.SCENARIO_RADII[self.scenario]
        elif max_rad is not None:
            self.max_rad = max_rad
        else:
            raise ValueError("Provide max_rad if scenario not recognized.")

        # Precompute angles
        self.angles_interp = np.linspace(0, 2*np.pi, n_interp_points+1)[:-1]
        self.angles_plot = np.linspace(0, 2*np.pi, 100)

        np.random.seed(self.seed)
        
    # -------------------------------------------------
    # MAIN SAMPLING
    # -------------------------------------------------
    def generate_samples(self):

        # Random refinery selection
        self.index_rand = np.random.choice(
            np.arange(len(self.x_coords)),
            size=self.N_rand,
            replace=True
        )

        # Radii sampling
        interp_vals_rand = np.random.uniform(low=self.min_rad, high=self.max_rad, 
                                             size=(self.n_interp_points, self.N_rand)
                                             )

        self.inputs_rand = np.concatenate((self.x_coords[self.index_rand][None, :],
            self.y_coords[self.index_rand][None, :],
            interp_vals_rand
        ))

        # Sample farms
        self._sample_farms()

        # Compute areas
        self.areas = self.calculate_areas() / 1e6

        return self._build_output()

    # -------------------------------------------------
    # FARM SAMPLING
    # -------------------------------------------------
    def _sample_farms(self):

        random_theta = []
        random_r = []

        for i in range(self.N_rand):

            radii_interp = np.interp(self.angles_plot, self.angles_interp, self.inputs_rand[2:, i], period=2*np.pi)

            Ftheta = np.cumsum(radii_interp) / np.sum(radii_interp)

            theta_rand = np.interp(
                np.random.uniform(size=self.n_samples_per_refinery),
                Ftheta,
                self.angles_plot
            )

            R_theta = np.interp(
                theta_rand,
                self.angles_interp,
                self.inputs_rand[2:, i],
                period=2*np.pi
            )

            r_rand = np.sqrt(np.random.uniform(size=self.n_samples_per_refinery)) * R_theta

            random_theta.extend(theta_rand)
            random_r.extend(r_rand)

        random_theta = np.degrees(np.array(random_theta)) # from radians to degrees
        random_r = np.array(random_r)

        # Convert to Cartesian
        x_rep = np.repeat(self.x_coords[self.index_rand], self.n_samples_per_refinery)
        y_rep = np.repeat(self.y_coords[self.index_rand], self.n_samples_per_refinery)

        self.x_farms = x_rep[None,:] + random_r[None,:] * np.cos(np.radians(random_theta))
        self.y_farms = y_rep[None,:] + random_r[None,:] * np.sin(np.radians(random_theta))

        self.x_bioref = self.x_coords[self.index_rand]
        self.y_bioref = self.y_coords[self.index_rand]

    # -------------------------------------------------
    # AREA CALCULATION
    # -------------------------------------------------
    def calculate_areas(self):

        R = self.inputs_rand.T[:, -self.n_interp_points:]
        Rk = np.roll(R, -1, axis=1)

        term1 = np.pi / self.n_interp_points * R**2
        term2 = np.pi / self.n_interp_points * R * (Rk - R)
        term3 = np.pi / (3 * self.n_interp_points) * (Rk - R)**2

        return np.sum(term1 + term2 + term3, axis=1)

    # -------------------------------------------------
    # OUTPUT FORMAT
    # -------------------------------------------------
    def _build_output(self):

        random_bio_to_farm = np.vstack([
            np.repeat(self.x_bioref, self.n_samples_per_refinery),
            np.repeat(self.y_bioref, self.n_samples_per_refinery),
            self.x_farms,
            self.y_farms
        ]).T

        return {
            "random_bio_to_farm": random_bio_to_farm,
            "inputs_rand": self.inputs_rand,
            "areas_rand": self.areas,
            "x_coords_bioref": self.x_bioref,
            "y_coords_bioref": self.y_bioref,
            "x_coords_farms_sampled": self.x_farms,
            "y_coords_farms_sampled": self.y_farms
        }

    # -------------------------------------------------
    # SAVE
    # -------------------------------------------------
    def save(self, output_folder, label):

        os.makedirs(output_folder, exist_ok=True)

        data = self._build_output()

        for name, arr in data.items():
            np.save(os.path.join(output_folder, f"{name}.npy"), arr)

        print(f"✅ Saved scenario: {label}")

    # -------------------------------------------------
    # STATIC LOADER
    # -------------------------------------------------
    @staticmethod
    def load(folder):

        data = {}
        for file in os.listdir(folder):
            if file.endswith(".npy"):
                data[file.replace(".npy", "")] = np.load(os.path.join(folder, file))

        return data
    
    def plot_irregular_areas(self, num_refineries=3, save_path=None):
        """Visualize the irregular areas around a few sampled refineries."""
        
        n_interp_points=self.n_interp_points

        angles_interp = np.linspace(0, 2*np.pi, n_interp_points, endpoint=False)
        angles_plot = np.linspace(0, 2*np.pi, 200)

        colors = [
            "#79bf82", "#ED586F", "#f3c354",
            "#60c1cf", "#90918e", "#a280b9", "#f98f60"
        ]

        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(8, 8))

        ax.set_theta_zero_location("E")
        ax.set_theta_direction(1)
        ax.grid(color='lightgray', linestyle='--', linewidth=0.8)

        for i in range(num_refineries):
            radii_interp = self.inputs_rand[2:, i] / 1000

            radii = np.interp(
                angles_plot,
                angles_interp,
                radii_interp,
                period=2*np.pi
            )

            ax.fill(
                angles_plot,
                radii,
                color=colors[i % len(colors)],
                alpha=0.4,
                edgecolor="k",
                linewidth=0.7,
                label=f"Refinery {i+1}"
            )

        ax.scatter(0, 0, color='black', s=30, zorder=5)
        ax.legend(frameon=False)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, transparent=True)

        plt.show()
    
    def plot_biorefinery_farms(self, bio_index = 1, save_path=None):
        """Plot a specific biorefinery index and its sampled farms with the irregular area boundary."""
        
        num_farms=self.n_samples_per_refinery
        n_interp_points=self.n_interp_points

        # --- Coordinates (convert to km)
        bio_x = self.x_bioref[bio_index] / 1000
        bio_y = self.y_bioref[bio_index] / 1000
        
        x_farms_reshaped = self.x_farms.reshape(self.N_rand, self.n_samples_per_refinery)
        y_farms_reshaped = self.y_farms.reshape(self.N_rand, self.n_samples_per_refinery)

        farm_x = x_farms_reshaped[bio_index, :num_farms] / 1000
        farm_y = y_farms_reshaped[bio_index, :num_farms] / 1000
        

        # Shift to origin
        farm_x -= bio_x
        farm_y -= bio_y

        # --- Angles
        angles_interp = np.linspace(0, 2*np.pi, n_interp_points, endpoint=False)
        angles_plot = np.linspace(0, 2*np.pi, 200)

        # --- Radii
        radii_interp = self.inputs_rand[2:, bio_index] / 1000
        radii = np.interp(angles_plot, angles_interp, radii_interp, period=2*np.pi)

        # Convert to Cartesian
        area_x = radii * np.cos(angles_plot)
        area_y = radii * np.sin(angles_plot)

        # --- Plot
        plt.figure(figsize=(8, 6))

        plt.fill(area_x, area_y, color='#79bf82', alpha=0.25,
                edgecolor='k', linewidth=0.7, label='Collection area')

        for fx, fy in zip(farm_x, farm_y):
            plt.plot([0, fx], [0, fy], 'k--', alpha=0.3)

        plt.scatter(farm_x, farm_y, c='#60c1cf', s=60, edgecolor='k', label='Farms')
        plt.scatter(0, 0, c='#ED586F', s=120, edgecolor='k', label='Biorefinery')

        plt.xlabel("X [km]")
        plt.ylabel("Y [km]")
        plt.legend(frameon=False)
        plt.axis('equal')
        plt.grid(alpha=0.4)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, transparent=True)

        plt.show()
    


# Example usage:
if __name__ == "__main__":
    sampler = RefineryFarmSampler(
        scenario="5perc",
        n_samples_per_refinery=100,
        N_rand=500_000
    )

    data = sampler.generate_samples()

    sampler.save("Outputs/500k_5perc", label="500k_5perc")
    
    sampler.plot_irregular_areas(num_refineries=3, save_path="Outputs/500k_5perc/irregular_areas.png")
    sampler.plot_biorefinery_farms(bio_index=1, save_path="Outputs/500k_5perc/biorefinery_farms.png")

    # To load later:
    loaded_data = RefineryFarmSampler.load("Outputs/500k_5perc")
    print(loaded_data.keys())
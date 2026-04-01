# This script generates the total logisticscosts per year and 
# total metric tons of biomass per year delivered to each biorefinery for 
# the samples generated in the "sample_generation" script. 
# It uses the distances calculated in the "calculate_distances_with_NN" script and 
# it needs a farm density scenario as input."
# @author: bianco3

from pathlib import Path
import re
import numpy as np
import geopandas as gpd
import joblib
from shapely.geometry import Point
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt

from project.paths import US_RAINFED, GREAT_LAKES, FARM_DENSITY, OUTPUT


class BioRefineryAnalysis:

    def __init__(self, scenario, farm_density_scenario=None):
        self.scenario = scenario
        self.base_folder = Path(OUTPUT)
        self.scenario_folder = self.base_folder / scenario

        if not self.scenario_folder.exists():
            raise FileNotFoundError(f"Scenario folder '{scenario}' not found")

        # 🔹 Auto-detect density if not provided
        if farm_density_scenario is None:
            self.farm_density = self._infer_density_from_scenario()
        else:
            self.farm_density = farm_density_scenario

        print(f"Using farm density: {self.farm_density}")

    # =========================================================
    # 🔹 Helpers
    # =========================================================
    def _infer_density_from_scenario(self):
        import re
        match = re.search(r'(\d+)perc', self.scenario)
        if match:
            return int(match.group(1)) / 100
        else:   
            raise ValueError("Could not infer farm density from scenario name")

    def _load_npy(self, name):
        path = self.scenario_folder / name
        if not path.exists():
            raise FileNotFoundError(f"{name} not found in scenario folder")
        return np.load(path)

    def _fill_nans(self, array):
        nan_mask = np.isnan(array)
        dist, indices = distance_transform_edt(nan_mask, return_indices=True)
        array[nan_mask] = array[tuple(indices[i][nan_mask] for i in range(array.ndim))]
        return array
    

    # =========================================================
    # 🔹 Load data
    # =========================================================
    def load_data(self):
        self.inputs_rand = self._load_npy('inputs_rand.npy')
        self.random_bio_to_farm = self._load_npy('random_bio_to_farm.npy')

        self.x_farms = self._load_npy('x_coords_farms_sampled.npy').flatten()
        self.y_farms = self._load_npy('y_coords_farms_sampled.npy').flatten()

        self.x_bio = self._load_npy('x_coords_bioref.npy')
        self.y_bio = self._load_npy('y_coords_bioref.npy')

        self.areas = self._load_npy('areas_rand.npy')

        self.N_rand = self.inputs_rand.shape[1]
        self.n_points = int(len(self.random_bio_to_farm) / self.N_rand)

    # =========================================================
    # 🔹 Load interpolated maps
    # =========================================================
    def load_maps(self):
        folder_maps = self.base_folder / "Interpolated_maps"
        if not folder_maps.exists():
            raise FileNotFoundError("Run interpolation script first")

        self.truck_map = self._fill_nans(np.load(folder_maps / 'Truck_prices_interpolated_map.npy'))
        self.farm_map = self._fill_nans(np.load(folder_maps / 'Farm_prices_interpolated_map_mis.npy'))
        self.yield_map = self._fill_nans(np.load(folder_maps / 'Crop_yields_interpolated_map_mis.npy'))

    # =========================================================
    # 🔹 Interpolation
    # =========================================================
    def interpolate(self):

        xmin, ymin, xmax, ymax = -1000350.6182353649, 259072.13554903926, 2263785.931181358, 3013035.67700148337

        grid_x = np.linspace(xmin, xmax, int((xmax - xmin) / 4000))
        grid_y = np.linspace(ymin, ymax, int((ymax - ymin) / 4000))

        points = np.column_stack((self.x_farms, self.y_farms))

        def interp(map_data):
            f = RegularGridInterpolator(
                (grid_x, grid_y),
                map_data,
                method='nearest',
                bounds_error=False,
                fill_value=None
            )
            return f(points)

        self.prices_farm = interp(self.farm_map)
        self.prices_truck = interp(self.truck_map)
        self.yields = interp(self.yield_map)

    # =========================================================
    # 🔹 Load NN distances
    # =========================================================
    def load_distances(self):
        path = self.scenario_folder / f'distances_bio_to_farm_{self.scenario}.npy'
        if not path.exists():
            raise FileNotFoundError("Run NN distance script first")

        self.distances = np.load(path)

    # =========================================================
    # 🔹 Function to apply penalty to points outside
    # =========================================================

    def build_us_penalty_grid(self, resolution=10000):

        USA = gpd.read_file(US_RAINFED).to_crs('EPSG:5070')
        GL = gpd.read_file(GREAT_LAKES).to_crs('EPSG:5070')

        USA = gpd.overlay(USA, GL, how="difference")

        xmin, ymin, xmax, ymax = USA.total_bounds

        x_vals = np.arange(xmin, xmax, resolution)
        y_vals = np.arange(ymin, ymax, resolution)

        grid_x, grid_y = np.meshgrid(x_vals, y_vals)
        points = np.c_[grid_x.ravel(), grid_y.ravel()]

        boundary = USA.geometry.unary_union

        distances = np.array([Point(x, y).distance(boundary) for x, y in points])
        distance_grid = distances.reshape(grid_x.shape)

        kd_tree = cKDTree(points)

        self.grid_x = grid_x
        self.grid_y = grid_y
        self.distance_grid = distance_grid
        self.kd_tree = kd_tree
    
    
    def compute_penalty(self):

        coords = np.column_stack((self.x_farms, self.y_farms))

        _, idx = self.kd_tree.query(coords)
        nearest = self.kd_tree.data[idx]

        grid_x_unique = np.unique(self.grid_x)
        grid_y_unique = np.unique(self.grid_y)

        interpolator = RegularGridInterpolator(
            (grid_x_unique, grid_y_unique),
            self.distance_grid.T,
            bounds_error=False,
            fill_value=np.max(self.distance_grid)
        )

        distances = interpolator((nearest[:, 0], nearest[:, 1]))

        self.penalty = np.where(distances > 8000., 0.0, 1.0)

    # =========================================================
    # 🔹 Compute costs
    # =========================================================
    def compute_costs(self):

        # Base cost
        total = (self.prices_truck * self.distances.flatten()) * 2 + self.prices_farm
        
        # print(total[:10])

        # Apply US penalty
        yields_penalized = self.yields * self.penalty
        total_penalized = total * self.penalty

        # Replace zero-cost (outside US) with large penalty
        penalty_value = np.nanmax(total_penalized)
        total_penalized[total_penalized == 0] = penalty_value

        # Load density
        density_gdf = gpd.read_file(FARM_DENSITY).to_crs('EPSG:5070')

        points_gdf = gpd.GeoDataFrame(
            geometry=[Point(x, y) for x, y in zip(self.x_farms, self.y_farms)],
            crs="EPSG:5070"
        )

        joined = gpd.sjoin(points_gdf, density_gdf, how="left", predicate="within")
        density = np.array(joined["percentage"].fillna(0))

        # Apply density + % planted
        yields_final = yields_penalized * density / 100 * self.farm_density # tons/ha
        costs_final = total_penalized * density / 100 * self.farm_density #in $/ha


        # Reshape
        yields_reshaped = yields_final.reshape(self.N_rand, self.n_points)
        costs_reshaped = costs_final.reshape(self.N_rand, self.n_points)
        

        mean_yields = np.mean(yields_reshaped, axis=1)
        mean_costs = np.mean(costs_reshaped, axis=1)

        # Totals
        production = yields_reshaped * self.areas[:,None] * 100
        total_costs = costs_reshaped * self.areas[:,None] * 100
        

        # Save outputs
        np.save(self.scenario_folder / 'production.npy', production) # in total tons per farm assuming all farms are miscanthus
        np.save(self.scenario_folder / 'Costs_biorefinery.npy', total_costs) # total $ 

        print("✅ Costs and production computed")

        print("✅ Analysis complete")

    # =========================================================
    # 🔹 Run all
    # =========================================================
    def run(self):
        print("Starting load_data")
        self.load_data()
        print("Starting load_maps")
        self.load_maps()
        print("Starting interpolate")
        self.interpolate()
        print("Starting load_distances")
        self.load_distances()
        print("Starting build_us_penalty_grid")
        self.build_us_penalty_grid()  
        print("Starting compute_penalty")
        self.compute_penalty()         
        print("Starting compute_costs")
        self.compute_costs()
        print("Done")


#%% Example usage
if __name__ == "__main__":
    print("Script started - initializing BioRefineryAnalysis")
    analysis = BioRefineryAnalysis("500k_5perc")
    print("Running analysis...")
    analysis.run()
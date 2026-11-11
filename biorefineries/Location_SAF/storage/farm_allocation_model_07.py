import os
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString, Point
from scipy.spatial import distance
import pandas as pd


# =============================================================================
# MODULE-LEVEL STATIC-DATA CACHE
#
# FeedstockAllocationModel is rebuilt from scratch many times over the course
# of a Monte Carlo run (once per sample x combo), but the underlying GIS
# reference layers (US boundary shapefile, rainfed-states shapefile,
# switchgrass/miscanthus point data) are identical every time. Loading and
# reprojecting them from disk on every instantiation is the single biggest
# avoidable cost in a large run. These helpers load each file once per
# process and hand back a cheap .copy() after that.
# =============================================================================
_STATIC_SHAPEFILE_CACHE = {}
_STATIC_CSV_POINTS_CACHE = {}


def _load_static_shapefile(path, target_crs='EPSG:5070'):
    """Load + reproject a shapefile once per process; reuse after that."""
    key = (path, target_crs)
    if key not in _STATIC_SHAPEFILE_CACHE:
        _STATIC_SHAPEFILE_CACHE[key] = gpd.read_file(path).to_crs(target_crs)
    return _STATIC_SHAPEFILE_CACHE[key].copy()


def _load_static_csv_points(path):
    """Load a lat/long CSV into a GeoDataFrame (EPSG:4326) once per process;
    reuse after that. Left in EPSG:4326 (not reprojected) so callers can
    apply a bounding_box filter in the same CRS as the original code before
    reprojecting themselves."""
    if path not in _STATIC_CSV_POINTS_CACHE:
        df = pd.read_csv(path)
        gdf = gpd.GeoDataFrame(
            df,
            crs='EPSG:4326',
            geometry=gpd.points_from_xy(df['long'], df['lat']),
        )
        _STATIC_CSV_POINTS_CACHE[path] = gdf

    return _STATIC_CSV_POINTS_CACHE[path]


class FeedstockAllocationModel:

    default_conversion_factors = {
        'working_days': (350, 'days/year'),
        'liters_per_gal': (3.785, 'L/gal'),
        'liters_per_m3': (1000, 'L/m³'),
        'density_of_ethanol': (747.58, 'kg/m³'),
        'kg_per_ton': (1000, 'kg/Mg'),
        'farm_area_ha': (4 * 4 * 100, 'ha'),
        'mis_to_ethanol_gallon': (113.9, 'gal ethanol/dry ton'),
        'sw_to_ethanol_gallon': (99.45, 'gal ethanol/dry ton'),
    }

    def __init__(
        self,
        feedstock='miscanthus',
        base_folder=None,
        num_points=1,
        samples_for_uncertainty=10,
        sizes_of_biorefineries=80,
        suitable_land_density=0.1,
        penalty_for_unmet_demand=0.3,
        delivered_biomass=0,
        harvested_biomass=0,
        storage_costs=0,
        bounding_box=None,
        seed=1,
        refinery_locations=None,
        moisture_content=0.2,
        yield_bonus=0.0,
        conversion_efficiency=113.9,
        extra_biomass_with_yield_bonus=True,
    ):
        """
        Parameters
        ----------
        extra_biomass_with_yield_bonus : bool, default True
            Controls how ``yield_bonus`` is applied.

            False:
                Yield bonus is applied before farm selection. Therefore,
                fewer farms may be contracted to meet the target biomass.

            True:
                Farms are selected/contracted using baseline yield, so the
                contracted farms are unchanged by the yield bonus. The yield
                bonus then changes the actual biomass delivered by those
                contracted farms.
                In our simulations, extra_biomass_with_yield_bonus is always True.

                IMPORTANT:
                The refinery pays for ALL actual biomass delivered, including
                biomass resulting from a positive yield bonus.
        """

        if base_folder is None:
            try:
                base_folder = os.path.join(
                    os.path.dirname(__file__),
                    'GIS_data'
                )
            except NameError:
                base_folder = os.path.join(
                    os.getcwd(),
                    'GIS_data'
                )

        self.base_folder = base_folder
        self._set_default_paths()

        self.feedstock = feedstock
        self.N = samples_for_uncertainty
        self.seed = seed
        self.refinery_locations = refinery_locations
        self.suitable_land_density = suitable_land_density
        self.penalty_for_unmet_demand = penalty_for_unmet_demand
        self.delivered_biomass = delivered_biomass
        self.harvested_biomass = harvested_biomass
        self.storage_costs = storage_costs
        self.moisture_content = moisture_content
        self.yield_bonus = yield_bonus
        self.conversion_efficiency = conversion_efficiency
        self.extra_biomass_with_yield_bonus = extra_biomass_with_yield_bonus

        if self.refinery_locations is not None:
            self.num_points = len(self.refinery_locations)
        else:
            self.num_points = num_points

        self.conversion_factors = self.default_conversion_factors.copy()
        self.sizes_of_biorefineries = sizes_of_biorefineries
        self.bounding_box = bounding_box

        if sizes_of_biorefineries != 80:
            print(
                f"[INFO] Biorefinery size set to "
                f"{sizes_of_biorefineries} MMgal/year."
            )

        self._load_all_data()
        self.results = {}

        # Cache for repeated MC runs
        self._distance_cache = None

    def _set_default_paths(self):
        self.paths = {
            'map_shapefile':
                os.path.join(
                    self.base_folder,
                    'US_map_for_boundary.shp'
                ),
            'rainfed_shapefile':
                os.path.join(
                    self.base_folder,
                    'USA-rainfedStates.shp'
                ),
            'switchgrass_data':
                os.path.join(
                    self.base_folder,
                    'switchgrass_data.csv'
                ),
            'miscanthus_data':
                os.path.join(
                    self.base_folder,
                    'miscanthus_data.csv'
                ),
            'transport_costs':
                os.path.join(
                    self.base_folder,
                    'States_with_transport.shp'
                ),
        }

    def update_path(self, key, new_path):
        if key in self.paths:
            self.paths[key] = new_path
            self._load_all_data()
        else:
            raise KeyError(
                f"{key} is not a recognized data file key."
            )

    def _load_all_data(self):
        # Static reference layers
        self.map_USA = _load_static_shapefile(
            self.paths['map_shapefile']
        )

        self.USA_rainfed = _load_static_shapefile(
            self.paths['rainfed_shapefile']
        )

        if self.bounding_box is not None:
            from shapely.geometry import box

            min_lon, min_lat, max_lon, max_lat = self.bounding_box
            bbox_poly = box(
                min_lon,
                min_lat,
                max_lon,
                max_lat
            )

        # Switchgrass
        sw_gdf = _load_static_csv_points(
            self.paths['switchgrass_data']
        ).copy()

        if self.bounding_box is not None:
            sw_gdf = sw_gdf[
                sw_gdf.geometry.within(bbox_poly)
            ]

        self.sw_data = sw_gdf.to_crs('EPSG:5070')

        # Miscanthus
        mis_gdf = _load_static_csv_points(
            self.paths['miscanthus_data']
        ).copy()

        if self.bounding_box is not None:
            mis_gdf = mis_gdf[
                mis_gdf.geometry.within(bbox_poly)
            ]

        self.mis_data = mis_gdf.to_crs('EPSG:5070')

    def print_conversion_factors(self):
        for key, (val, unit) in self.conversion_factors.items():
            print(f"  {key}: {val} {unit}")

    def set_conversion_factor(self, name, value, unit=None):
        if unit is None and name in self.conversion_factors:
            unit = self.conversion_factors[name][1]
        elif unit is None:
            unit = ''

        self.conversion_factors[name] = (value, unit)

    @staticmethod
    def sort_by_min_distance(first_mat, sorted_mat, num_points):
        return np.take_along_axis(
            first_mat,
            sorted_mat,
            axis=1
        )

    def _sample_active_farms(self, crop_data):
        np.random.seed(self.seed)

        active_mask = (
            np.random.rand(len(crop_data))
            < self.suitable_land_density
        )

        return crop_data[active_mask].copy()

    def _calculate_extra_biomass_costs(
        self,
        unit_price,
        actual_delivered_biomass_Mg=None,
    ):
        """
        Calculate the amount and cost of biomass that must be purchased
        beyond the biomass actually delivered by contracted farms.

        The refinery demand is calculated from refinery size and conversion
        efficiency. Any shortfall is purchased at the farm price multiplied
        by the unmet-demand penalty.

        Positive yield bonuses are already included in
        ``actual_delivered_biomass_Mg`` and therefore do NOT get added again.
        """

        biomass_need_total_dry = (
            self.sizes_of_biorefineries
            * 1e6
            / self.conversion_efficiency
        )

        biomass_need_total_wet = (
            biomass_need_total_dry
            / (1 - self.moisture_content)
        )

        if actual_delivered_biomass_Mg is not None:
            supplied_biomass = actual_delivered_biomass_Mg
        else:
            supplied_biomass = self.delivered_biomass

        extra_biomass_needed = np.maximum(
            biomass_need_total_wet - supplied_biomass,
            0.0
        )

        extra_cost_total = (
            extra_biomass_needed
            * unit_price
            * (1 + self.penalty_for_unmet_demand)
        )

        return (
            extra_biomass_needed,
            extra_cost_total
        )

    def _precompute_distance_cache(
        self,
        crop_data_active,
        possible_locations
    ):
        """Vectorized distance matrix precomputation using scipy's cdist."""
        from scipy.spatial.distance import cdist

        loc_coords = np.column_stack((
            possible_locations.geometry.x,
            possible_locations.geometry.y
        ))

        farm_coords = np.column_stack((
            crop_data_active.geometry.x,
            crop_data_active.geometry.y
        ))

        dist = cdist(
            loc_coords,
            farm_coords
        ) / 1000.0  # km

        sorted_farms = np.argsort(
            dist,
            axis=1
        )

        dist_sorted = np.take_along_axis(
            dist,
            sorted_farms,
            axis=1
        )

        return {
            "dist": dist,
            "sorted_farms": sorted_farms,
            "dist_sorted": dist_sorted,
        }

    def calculate_location_parameters(self):

        num_points = self.num_points
        cf = self.conversion_factors

        # ------------------------------------------------------------------
        # Select feedstock
        # ------------------------------------------------------------------
        if self.feedstock == 'switchgrass':
            crop_data = self.sw_data

        elif self.feedstock == 'miscanthus':
            crop_data = self.mis_data

        else:
            raise ValueError(
                f'Unsupported feedstock: {self.feedstock}'
            )

        crop_data_active = self._sample_active_farms(
            crop_data
        )

        factor_MC = 1 - self.moisture_content

        # ------------------------------------------------------------------
        # Yield used for contracting / farm selection
        #
        # extra_biomass_with_yield_bonus=False
        #     Bonus is applied before farm selection.
        #
        # extra_biomass_with_yield_bonus=True
        #     Farms are selected using baseline yield. The yield bonus
        #     changes actual biomass delivered by those same farms.
        # ------------------------------------------------------------------
        if self.extra_biomass_with_yield_bonus:

            # Baseline yield used for contracting
            yield_crop_contract = (
                np.array(crop_data_active.y_data)
                / factor_MC
            )

        else:

            # Original behavior: bonus affects farm selection
            yield_crop_contract = (
                np.array(crop_data_active.y_data)
                * (1 + self.yield_bonus)
                / factor_MC
            )

        # Actual yield used for biomass delivered
        yield_crop = (
            np.array(crop_data_active.y_data)
            * (1 + self.yield_bonus)
            / factor_MC
        )

        GHG_crop = (
            np.array(crop_data_active.Net_CI)
            * factor_MC
        )

        Price_crop = (
            np.array(crop_data_active.Farm_price)
            * factor_MC
        )

        crop_data_active['Adjusted_Price'] = (
            crop_data_active['Farm_price']
            * factor_MC
        )

        area = cf['farm_area_ha'][0]

        # ------------------------------------------------------------------
        # Farm biomass used for CONTRACTING / FARM SELECTION
        #
        # This remains based on baseline yield when
        # extra_biomass_with_yield_bonus=True.
        # ------------------------------------------------------------------
        size_farms_contract = (
            area
            * yield_crop_contract
            * np.ones(num_points)[..., None]
        )

        probability_crop = (
            yield_crop_contract
            / np.sum(yield_crop_contract)
        )

        sizes_of_crop_demand = self.harvested_biomass

        # ------------------------------------------------------------------
        # Refinery location -> nearest crop grid point
        # ------------------------------------------------------------------
        if self.refinery_locations is not None:

            gdf_input = gpd.GeoDataFrame(
                geometry=[
                    Point(lon, lat)
                    for lat, lon in self.refinery_locations
                ],
                crs='EPSG:4326'
            )

            gdf_input_5070 = gdf_input.to_crs(
                'EPSG:5070'
            )

            if crop_data.crs != 'EPSG:5070':
                crop_data = crop_data.to_crs(
                    'EPSG:5070'
                )

            from scipy.spatial import cKDTree

            crop_coords = np.column_stack((
                crop_data.geometry.x,
                crop_data.geometry.y
            ))

            refinery_coords = np.column_stack((
                gdf_input_5070.geometry.x,
                gdf_input_5070.geometry.y
            ))

            tree = cKDTree(crop_coords)

            _, nearest_idx = tree.query(
                refinery_coords,
                k=1
            )

            possible_locations = crop_data.iloc[
                nearest_idx
            ]

            self.num_points = len(
                possible_locations
            )

        else:

            np.random.seed(self.seed)

            chosen_ind = np.random.choice(
                len(yield_crop_contract),
                self.num_points,
                replace=False,
                p=probability_crop
            )

            possible_locations = crop_data.iloc[
                chosen_ind
            ]

        possible_locations = (
            possible_locations
            .sort_values("siteid")
            .reset_index(drop=True)
        )

        self.candidate_locations = possible_locations

        # ------------------------------------------------------------------
        # Distance calculations
        # ------------------------------------------------------------------
        if self._distance_cache is None:

            self._distance_cache = (
                self._precompute_distance_cache(
                    crop_data_active,
                    possible_locations
                )
            )

        dist = self._distance_cache["dist"]

        sorted_farms = self._distance_cache[
            "sorted_farms"
        ]

        dist_sorted = self._distance_cache[
            "dist_sorted"
        ]

        # ------------------------------------------------------------------
        # Sort CONTRACTED farm sizes by distance
        # ------------------------------------------------------------------
        farm_sizes_sorted = (
            FeedstockAllocationModel.sort_by_min_distance(
                size_farms_contract,
                sorted_farms,
                num_points
            )
        )

        farm_sizes_cum = np.cumsum(
            farm_sizes_sorted,
            axis=1
        )

        prices_by_sorted = (
            crop_data_active['Adjusted_Price'].values[
                sorted_farms
            ]
        )

        # ------------------------------------------------------------------
        # Number of farms needed to meet the contracted biomass target
        # ------------------------------------------------------------------
        if np.any(farm_sizes_cum[:, -1] < sizes_of_crop_demand - 1e-8):
            raise ValueError('Insufficient farm biomass for the requested contract')
        contracted_by_farm = np.minimum(farm_sizes_sorted, np.maximum(
            sizes_of_crop_demand - (farm_sizes_cum - farm_sizes_sorted), 0.0))
        contract_fraction = np.divide(contracted_by_farm, farm_sizes_sorted,
            out=np.zeros_like(farm_sizes_sorted), where=farm_sizes_sorted > 0)
        farm_num = np.sum(contracted_by_farm > 0, axis=1)

        # ------------------------------------------------------------------
        # Contracted farm mask
        # ------------------------------------------------------------------
        col_idx = np.arange(
            farm_sizes_sorted.shape[1]
        )

        mask = (
            col_idx[None, :]
            < farm_num[:, None]
        )

        # ------------------------------------------------------------------
        # Actual biomass produced by contracted farms
        #
        # IMPORTANT:
        # This uses yield_crop, which includes the yield bonus.
        #
        # Therefore:
        #
        #   negative yield -> less biomass delivered
        #   positive yield -> more biomass delivered
        #
        # and ALL of this actual biomass is paid for.
        # ------------------------------------------------------------------
        actual_biomass_all_sorted = (
            area
            * yield_crop
            * np.ones(num_points)[..., None]
        )

        actual_biomass_sorted = (
            FeedstockAllocationModel.sort_by_min_distance(
                actual_biomass_all_sorted,
                sorted_farms,
                num_points
            )
        )

        actual_biomass_harvested = np.where(
            mask,
            actual_biomass_sorted * contract_fraction,
            0.0
        )

        total_actual_Mg = (
            actual_biomass_harvested.sum(axis=1)
        )

        # ------------------------------------------------------------------
        # Baseline biomass represented by the contracted farms
        # ------------------------------------------------------------------
        contracted_baseline_biomass = np.where(
            mask,
            contracted_by_farm,
            0.0
        ).sum(axis=1)

        # ------------------------------------------------------------------
        # Bonus biomass
        #
        # This is now purely a reporting quantity.
        # It is INCLUDED in total_actual_Mg and therefore INCLUDED
        # in total_feedstock_cost.
        # ------------------------------------------------------------------
        bonus_biomass_Mg = np.maximum(
            total_actual_Mg
            - contracted_baseline_biomass,
            0.0
        )

        # ------------------------------------------------------------------
        # Populate harvested_farms
        #
        # biomass_supplied_Mg now represents ACTUAL biomass delivered,
        # not baseline contracted biomass.
        # ------------------------------------------------------------------
        harvested_farms = []

        for i in range(num_points):

            n_farms = farm_num[i]

            harvested_idx = (
                sorted_farms[i][:n_farms]
            )

            harvested_subset = (
                crop_data_active
                .iloc[harvested_idx]
                .copy()
            )

            harvested_subset[
                'biomass_supplied_Mg'
            ] = actual_biomass_harvested[
                i
            ][:n_farms]

            harvested_subset[
                'distance_to_refinery_km'
            ] = dist_sorted[
                i
            ][:n_farms]

            harvested_farms.append(
                harvested_subset
            )

        self.harvested_farms = harvested_farms

        # ------------------------------------------------------------------
        # FEEDSTOCK COST
        #
        # The refinery pays for ALL actual biomass delivered by contracted
        # farms, including any positive yield bonus.
        # ------------------------------------------------------------------
        total_feedstock_cost = (
            actual_biomass_harvested
            * prices_by_sorted
        ).sum(axis=1)

        farm_cost_all_loc = np.where(
            total_actual_Mg > 0,
            total_feedstock_cost / total_actual_Mg,
            0.0
        )

        # ------------------------------------------------------------------
        # GHG
        # ------------------------------------------------------------------
        CI_by_sorted = (
            crop_data_active['Net_CI'].values[
                sorted_farms
            ]
            * 0.8
        )

        sum_CI = (CI_by_sorted * actual_biomass_harvested).sum(axis=1)

        farm_num_safe = np.where(
            farm_num > 0,
            farm_num,
            1
        )

        farm_GHG_all_loc = (
            np.divide(sum_CI, total_actual_Mg, out=np.zeros_like(sum_CI), where=total_actual_Mg > 0)
        )

        # ------------------------------------------------------------------
        # Mean transport distance
        #
        # Weight by ACTUAL biomass delivered.
        # ------------------------------------------------------------------
        with np.errstate(
            divide='ignore',
            invalid='ignore'
        ):

            mean_distances_bio = (
                np.sum(
                    dist_sorted
                    * actual_biomass_harvested,
                    axis=1
                )
                / total_actual_Mg
            )

            mean_distances_bio = np.where(
                total_actual_Mg > 0,
                mean_distances_bio,
                0.0
            ).tolist()

        # ------------------------------------------------------------------
        # EXTRA BIOMASS
        #
        # Compare refinery demand against ACTUAL biomass delivered.
        #
        # Therefore a positive yield bonus reduces the amount of extra
        # biomass that needs to be purchased.
        # ------------------------------------------------------------------
        extra_biomass_amount, extra_biomass_cost_total = (
            self._calculate_extra_biomass_costs(
                farm_cost_all_loc,
                actual_delivered_biomass_Mg=total_actual_Mg
            )
        )

        # ------------------------------------------------------------------
        # Store results
        # ------------------------------------------------------------------

        # Average price of actual biomass delivered
        self.results[
            'feedstock_unit_price_contract'
        ] = farm_cost_all_loc / 1000

        # Cost of ALL actual biomass harvested
        self.results[
            'Total_contracted_feedstock_cost'
        ] = total_feedstock_cost

        # # Cost based on actual harvested biomass
        # self.results[
        #     'Total_supplied_feedstock_costs'
        # ] = (
        #     farm_cost_all_loc
        #     * total_actual_Mg
        # )

        self.results[
            'feedstock_unit_GHG'
        ] = farm_GHG_all_loc / 1000

        self.results[
            'Price_crop'
        ] = Price_crop

        self.results[
            'yield_crop'
        ] = yield_crop

        self.results[
            'GHG_crop'
        ] = GHG_crop

        self.results[
            'possible_locations'
        ] = possible_locations

        self.results[
            'farm_to_bio_mean_distances'
        ] = mean_distances_bio

        self.results[
            'extra_biomass_amount'
        ] = extra_biomass_amount

        self.results[
            'extra_biomass_cost_total'
        ] = extra_biomass_cost_total

        # Baseline biomass represented by contracted farms
        self.results[
            'contracted_biomass_Mg'
        ] = contracted_baseline_biomass

        # Actual biomass delivered, including yield bonus or yield loss
        self.results[
            'delivered_biomass_Mg'
        ] = total_actual_Mg

        # Positive difference between actual and baseline production
        self.results[
            'bonus_biomass_Mg'
        ] = bonus_biomass_Mg

        self.results[
            'farm_count'
        ] = farm_num

        self.results[
            'yield_bonus_applied'
        ] = self.yield_bonus

        # Active farms
        self.active_farms = crop_data_active

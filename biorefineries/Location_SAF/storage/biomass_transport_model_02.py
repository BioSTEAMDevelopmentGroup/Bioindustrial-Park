## Code to calculate transport costs
## Can choose between road network calculations (may take a lot of time, only use for small number of points),
## using neural network predictions to estimate distances,
## or straight-line distance + tortuosity adjustment

"""
BiomassTransport
================
A class for calculating transport distances and costs for biomass supply chains.

Coordinate Reference Systems
-----------------------------
- Farm 'lat' / 'long' columns        : EPSG:4326 decimal degrees
- Farm 'geometry' column             : EPSG:5070 Albers metres
- Farm 'distance_to_refinery_km'     : straight-line distance pre-computed by the
                                       upstream model in EPSG:5070 — used directly,
                                       no reprojection or haversine needed.
- Refinery coordinates               : EPSG:4326 decimal degrees

Distance Methods
----------------
- "straight_line"  : uses 'distance_to_refinery_km' from the farm DataFrame as-is.
- "tortuosity"     : distance_to_refinery_km × tortuosity factor.
                     If use_per_state_tortuosity=True (default), the factor is read
                     from Tortuosity_per_state.xlsx ('Mean mis' column) per farm state.
                     Falls back to tortuosity_factor constant if the state value is NA.
- "road_network"   : shortest road path via osmnx + NetworkX.
                     Uses farm 'lat'/'long' (EPSG:4326) and refinery coords.
                     Downloads the OSM graph automatically on first call.

State Resolution
----------------
If the farm DataFrame has a 'state' column it is used directly.
Otherwise, states are resolved via a spatial join against a US boundary
shapefile (GIS_data/US_map_for_boundary.shp). The shapefile is loaded once
at initialisation and reprojected to EPSG:4326 to match the farm coordinates.
The state name is read from the 'NAME' column (configurable via state_col).

Install dependencies
--------------------
    pip install geopandas                 # state lookup via shapefile
    pip install osmnx networkx            # only for road_network method
"""

import math
import warnings
import numpy as np
import pandas as pd
from typing import Union


def load_transport_costs_csv(csv_path: str) -> dict:
    """Standalone version of BiomassTransport._load_transport_costs so a
    driver script can load it once up front and pass the resulting dict in
    via cost_lookup=, instead of re-reading the CSV on every instantiation."""
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()

    state_col = next((c for c in df.columns if "state" in c), None)
    cost_col  = next((c for c in df.columns if "cost"  in c), None)

    if state_col is None or cost_col is None:
        raise ValueError(
            "transport_costs_csv must contain a 'state' column and a cost column "
            "(e.g. 'cost_per_km_per_trip')."
        )
    return dict(zip(df[state_col].str.strip().str.lower(), df[cost_col]))


class BiomassTransport:
    """
    Calculate biomass transport distances and costs for one or more refineries.

    Parameters
    ----------
    refinery_locations : list of (lat, lon) tuples
        EPSG:4326 coordinates of each refinery.
    harvested_farms : list of pd.DataFrame
        One DataFrame per refinery (pre-assigned by the upstream model, i.e.
        harvested_farms[i] belongs to refinery_locations[i]).
        Required columns:
            'lat', 'long'               — EPSG:4326 decimal degrees
            'distance_to_refinery_km'   — straight-line distance (EPSG:5070-based)
            'biomass_supplied_Mg'       — biomass at this farm in Mg
        Optional columns:
            'state'                     — if present, used directly for cost lookup
                                          (skips reverse geocoding)
    transport_costs_csv : str
        Path to a CSV with at least two columns:
            state                  — US state name (e.g. "Oklahoma")
            cost_per_km_per_trip   — freight cost in $/km per truck trip
    transport_capacity_Mg : float
        Maximum biomass per truck trip in Mg.
        Used to convert trip-based costs to $/Mg/km, and to count trips.
    amount_transported_Mg : float or list of float
        Total biomass to be transported per refinery.
        A single float is broadcast to all refineries.
    distance_method : str
        One of 'straight_line', 'tortuosity', 'road_network'.
        Default: 'tortuosity'.
    tortuosity_factor : float
        Multiplier for the 'tortuosity' method.  Default: 1.3.
    shapefile_path : str
        Path to the US boundary shapefile (e.g. 'GIS_data/US_map_for_boundary.shp').
        Used to resolve state names when the farm DataFrame has no 'state' column.
        The shapefile is loaded and reprojected to EPSG:4326 once at initialisation.
    state_col : str
        Column name in the shapefile that contains state names.  Default: 'NAME'.
    use_per_state_tortuosity : bool
        If True (default) and distance_method == 'tortuosity', per-state tortuosity
        factors are loaded from tortuosity_xlsx and applied per farm.  The constant
        tortuosity_factor is used as a fallback for states with a missing ('Mean mis')
        value.  If False, the constant tortuosity_factor is applied to all farms.
    tortuosity_xlsx : str
        Path to the Excel file with per-state tortuosity factors.
        Expected columns: 'State' and 'Mean mis'.
        Default: 'GIS_data/Tortuosity_per_state.xlsx'.
    default_cost_per_km : float
        Fallback cost ($/km/trip) when a state is not found in the CSV.
        Emits a warning when used.  Default: 0.0.
    """

    VALID_METHODS = {"straight_line", "tortuosity", "road_network", "neural_network"}

    def __init__(
        self,
        refinery_locations: list,
        harvested_farms: list,
        transport_costs_csv: str,
        # transport_capacity_Mg: float,
        amount_transported_Mg: Union[float, list],
        bale_capacity_per_truck: Union[int, dict] = 42, # Round bale = 42, Rectangular bale = 36
        bale_mass_Mg: Union[float, dict] = 0.27, # Round bale = 0.27, Rectangular bale = 0.35
        shapefile_path: str = "GIS_data/US_map_for_boundary.shp",
        state_col: str = "NAME",
        use_per_state_tortuosity: bool = True,
        tortuosity_xlsx: str = "GIS_data/Tortuosity_per_state.xlsx",
        distance_method: str = "tortuosity",
        tortuosity_factor: float = 1.3,
        default_cost_per_km: float = 0.0,
        nn_model_path: str = 'single_bio_to_farm_distance_model.pth',
        nn_scaler_path: str = 'scaler.pkl',
        nn_device: str = "cpu",
        us_states_gdf=None,          # optional pre-loaded GeoDataFrame — skips
                                      # gpd.read_file()+to_crs() on every call
        cost_lookup: dict = None,    # optional pre-loaded {state: cost} dict —
                                      # skips re-reading transport_costs_csv
        tortuosity_lookup: dict = None,  # optional pre-loaded tortuosity dict

    ):
        # ── validate alignment ────────────────────────────────────────────
        if len(refinery_locations) != len(harvested_farms):
            raise ValueError(
                f"refinery_locations has {len(refinery_locations)} entries but "
                f"harvested_farms has {len(harvested_farms)} DataFrames — must match."
            )

        self.refinery_locations = refinery_locations
        self.harvested_farms    = [df.copy() for df in harvested_farms]

        # ── validate required columns ─────────────────────────────────────
        required = {"lat", "long", "distance_to_refinery_km", "biomass_supplied_Mg"}
        for i, df in enumerate(self.harvested_farms):
            missing = required - set(df.columns)
            if missing:
                raise ValueError(
                    f"harvested_farms[{i}] is missing required columns: {missing}"
                )

        # ── transport capacity & amounts ──────────────────────────────────
        # if transport_capacity_Mg <= 0:
        #     raise ValueError("transport_capacity_Mg must be > 0.")
        # self.transport_capacity_Mg = transport_capacity_Mg

        self.bale_capacity_per_truck = bale_capacity_per_truck
        self.bale_mass_Mg = bale_mass_Mg

        if isinstance(amount_transported_Mg, (int, float)):
            self.amount_transported_Mg = (
                [float(amount_transported_Mg)] * len(refinery_locations)
            )
        else:
            if len(amount_transported_Mg) != len(refinery_locations):
                raise ValueError(
                    "amount_transported_Mg list length must match number of refineries."
                )
            self.amount_transported_Mg = [float(v) for v in amount_transported_Mg]

        # ── distance method ───────────────────────────────────────────────
        if distance_method not in self.VALID_METHODS:
            raise ValueError(f"distance_method must be one of {self.VALID_METHODS}.")
        self.distance_method          = distance_method
        self.tortuosity_factor        = tortuosity_factor
        self.use_per_state_tortuosity = use_per_state_tortuosity
        if tortuosity_lookup is not None:
            self._tortuosity_lookup = tortuosity_lookup
        else:
            self._tortuosity_lookup = (
                self._load_tortuosity(tortuosity_xlsx)
                if use_per_state_tortuosity and distance_method == "tortuosity"
                else {}
            )

        # ── transport costs lookup ────────────────────────────────────────
        # Pass cost_lookup= to skip re-reading transport_costs_csv on every
        # BiomassTransport() instantiation (it's the same file every call in
        # a Monte Carlo loop).
        self.default_cost_per_km = default_cost_per_km
        self._cost_lookup = (
            cost_lookup if cost_lookup is not None
            else self._load_transport_costs(transport_costs_csv)
        )

        # ── US boundary shapefile (for state lookup) ──────────────────────
        # Pass us_states_gdf= to skip gpd.read_file()+to_crs() on every call.
        self.state_col = state_col
        self._us_states_gdf = (
            us_states_gdf if us_states_gdf is not None
            else self._load_shapefile(shapefile_path)
        )

        # ── lazy-loaded OSM road graph ────────────────────────────────────
        self._road_graph = None

        # ── neural-network distance model ──────────────────────────────
        self.nn_model_path = nn_model_path
        self.nn_scaler_path = nn_scaler_path
        self.nn_device = nn_device
        self._nn_model = None
        self._nn_scaler = None


        if self.distance_method == "neural_network":
            if self.nn_model_path is None:
                raise ValueError(
                    "nn_model_path must be provided when "
                    "distance_method='neural_network'."
                )
            if self.nn_scaler_path is None:
                raise ValueError(
                    "nn_scaler_path must be provided when "
                    "distance_method='neural_network'."
                )
            self._load_nn_model()
            self._load_nn_scaler()

        # ── results ───────────────────────────────────────────────────────
        self.results = None


    # ======================================================================
    # Load neural network function
    # ======================================================================

    def _load_nn_model(self):
        """
        Load trained PyTorch neural network for distance prediction.
        """
        try:
            import torch
            import torch.nn as nn
        except ImportError:
            raise ImportError(
                "PyTorch is required for neural_network distance method.\n"
                "Install with: pip install torch"
            )

        class NeuralNetwork(nn.Module):
            def __init__(self, input_size):
                super().__init__()

                self.model = nn.Sequential(
                    nn.Linear(input_size, 128),
                    nn.ReLU(),

                    nn.Linear(128, 128),
                    nn.ReLU(),

                    nn.Linear(128, 128),
                    nn.ReLU(),

                    nn.Linear(128, 128),
                    nn.ReLU(),

                    nn.Linear(128, 1)
                )

            def forward(self, x):
                return self.model(x)

        model = NeuralNetwork(input_size=4)

        model.load_state_dict(
            torch.load(
                self.nn_model_path,
                map_location=self.nn_device
            )
        )

        model.to(self.nn_device)
        model.eval()

        self._nn_model = model

    def _to_epsg5070(self, lats, lons):
        """
        Convert EPSG:4326 → EPSG:5070
        """
        from pyproj import Transformer

        transformer = Transformer.from_crs(
            "EPSG:4326",
            "EPSG:5070",
            always_xy=True
        )

        x, y = transformer.transform(lons, lats)

        return np.column_stack([x, y])

    def _load_nn_scaler(self):
        """
        Load sklearn scaler used during neural-network training.
        """

        import joblib

        self._nn_scaler = joblib.load(self.nn_scaler_path)

    # ======================================================================
    # Bale type per farm helpers
    # ======================================================================

    def _get_bale_params(self, row: pd.Series):
        """
        Returns (bale_capacity, bale_mass_Mg)
        """

        # default scalar case
        if isinstance(self.bale_capacity_per_truck, (int, float)):
            return float(self.bale_capacity_per_truck), float(self.bale_mass_Mg)

        # optional dict-based extension (future-proof)
        state = str(row.get("state_used", "unknown")).lower()

        cap = self.bale_capacity_per_truck.get(state, 42)
        mass = self.bale_mass_Mg.get(state, 0.25)

        return float(cap), float(mass)

    # ======================================================================
    # Cost helpers
    # ======================================================================

    def _load_transport_costs(self, csv_path: str) -> dict:
        """Load state → cost_per_km_per_trip from CSV. Column names are flexible."""
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip().str.lower()

        state_col = next((c for c in df.columns if "state" in c), None)
        cost_col  = next((c for c in df.columns if "cost"  in c), None)

        if state_col is None or cost_col is None:
            raise ValueError(
                "transport_costs_csv must contain a 'state' column and a cost column "
                "(e.g. 'cost_per_km_per_trip')."
            )
        return dict(zip(df[state_col].str.strip().str.lower(), df[cost_col]))

    def _get_cost_for_state(self, state: str) -> float:
        """Return cost per km per trip for the given state name."""
        key = str(state).strip().lower()
        if key not in self._cost_lookup:
            warnings.warn(
                f"State '{state}' not found in transport costs CSV. "
                f"Using default: {self.default_cost_per_km} $/km/trip."
            )
            return self.default_cost_per_km
        return self._cost_lookup[key]

    # def _trips_needed(self, biomass_Mg: float) -> int:
    #     """Ceiling number of truck trips for a given biomass quantity."""
    #     return math.ceil(biomass_Mg / self.transport_capacity_Mg)

    def _trips_needed(self, biomass_Mg: float, row: pd.Series) -> int:
        bale_cap, bale_mass = self._get_bale_params(row)

        bales = biomass_Mg / bale_mass
        return math.ceil(bales / bale_cap)

    # ======================================================================
    # State resolution via shapefile spatial join
    # ======================================================================

    def _load_shapefile(self, shapefile_path: str):
        """
        Load and reproject the US boundary shapefile to EPSG:4326.
        Keeps only the state name column to minimise memory use.
        Called once at __init__ time.
        """
        try:
            import geopandas as gpd
        except ImportError:
            raise ImportError(
                "geopandas is required for state lookup via shapefile.\n"
                "Install with:  pip install geopandas\n"
                "Or add a 'state' column to your farm DataFrames to skip this step."
            )

        gdf = gpd.read_file(shapefile_path)

        if self.state_col not in gdf.columns:
            raise ValueError(
                f"Column '{self.state_col}' not found in shapefile. "
                f"Available columns: {list(gdf.columns)}"
            )

        # Keep only what we need and ensure EPSG:4326 to match farm lat/long
        gdf = gdf[[self.state_col, "geometry"]].to_crs(epsg=4326)
        return gdf

    def _lookup_states_from_shapefile(self, df: pd.DataFrame) -> pd.Series:
        """
        Spatial join of farm points (EPSG:4326) against the US boundary GeoDataFrame.

        All farms in df are joined in a single vectorised call.
        Returns a Series of state name strings aligned to df's index.
        Falls back to 'unknown' for any point that does not intersect a polygon
        (e.g. farms very close to a state border that fall outside due to projection).
        """
        import geopandas as gpd
        from shapely.geometry import Point

        # Build a GeoDataFrame from farm lat/long (already EPSG:4326)
        farm_gdf = gpd.GeoDataFrame(
            df[[]],                                          # no extra columns needed
            geometry=gpd.points_from_xy(df["long"], df["lat"]),
            crs="EPSG:4326",
        )

        # Spatial join: each farm point → state polygon
        joined = gpd.sjoin(
            farm_gdf,
            self._us_states_gdf,
            how="left",
            predicate="within",
        )

        # Fill any unmatched points (e.g. on a border) with 'unknown'
        states = joined[self.state_col].fillna("unknown")
        states.index = df.index   # restore original index alignment
        return states

    # ======================================================================
    # Distance methods
    # ======================================================================

    def _load_tortuosity(self, xlsx_path: str) -> dict:
        """
        Load per-state tortuosity factors from Excel.

        Expected columns: 'State', 'Mean mis'.
        NA values are kept as NaN so _distance_tortuosity can fall back to the
        constant tortuosity_factor for those states.
        Returns a dict: lowercase state name → float (or NaN).
        """
        df = pd.read_excel(xlsx_path)
        df.columns = df.columns.str.strip()

        if "State" not in df.columns or "Mean mis" not in df.columns:
            raise ValueError(
                "Tortuosity_per_state.xlsx must contain 'State' and 'Mean mis' columns. "
                f"Found: {list(df.columns)}"
            )

        lookup = {}
        for _, row in df.iterrows():
            key   = str(row["State"]).strip().lower()
            value = row["Mean mis"]
            lookup[key] = float(value) if pd.notna(value) else float("nan")
        return lookup

    def _get_tortuosity_factor(self, state: str) -> float:
        """
        Return the tortuosity factor for a given state.

        Priority:
          1. Per-state value from the xlsx (if use_per_state_tortuosity=True and not NA)
          2. Constant self.tortuosity_factor as fallback
        """
        if not self.use_per_state_tortuosity or not self._tortuosity_lookup:
            return self.tortuosity_factor

        key    = str(state).strip().lower()
        factor = self._tortuosity_lookup.get(key, float("nan"))

        if math.isnan(factor):
            warnings.warn(
                f"No valid tortuosity factor for state '{state}' in xlsx "
                f"(missing or NA). Using constant fallback: {self.tortuosity_factor}."
            )
            return self.tortuosity_factor

        return factor

    def _distance_straight_line(self, farm_row: pd.Series) -> float:
        """Return the precomputed straight-line distance (EPSG:5070-based, in km)."""
        return float(farm_row["distance_to_refinery_km"])

    def _distance_tortuosity(self, farm_row: pd.Series) -> float:
        """
        Precomputed straight-line distance × tortuosity factor.
        Uses per-state factor from xlsx if use_per_state_tortuosity=True,
        falling back to the constant tortuosity_factor for missing states.
        """
        state  = farm_row.get("state_used", "unknown")
        factor = self._get_tortuosity_factor(state)
        return float(farm_row["distance_to_refinery_km"]) * factor

    def _load_road_graph(self):
        """
        Lazily download one OSM driving graph covering all refineries and farms.
        Called only when distance_method == 'road_network'.
        """
        try:
            import osmnx as ox
        except ImportError:
            raise ImportError(
                "osmnx and networkx are required for distance_method='road_network'.\n"
                "Install with:  pip install osmnx networkx"
            )

        if self._road_graph is not None:
            return self._road_graph

        print("Downloading OSM road network (this may take a moment)…")

        all_lats, all_lons = [], []
        for i, (rlat, rlon) in enumerate(self.refinery_locations):
            all_lats.append(rlat);  all_lons.append(rlon)
            all_lats.extend(self.harvested_farms[i]["lat"].tolist())
            all_lons.extend(self.harvested_farms[i]["long"].tolist())

        center_lat = float(np.mean(all_lats))
        center_lon = float(np.mean(all_lons))

        # Haversine needed only here to size the OSM download radius
        def _hav_km(lat1, lon1, lat2, lon2):
            R = 6371.0
            a = (math.sin(math.radians(lat2 - lat1) / 2) ** 2
                 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
                 * math.sin(math.radians(lon2 - lon1) / 2) ** 2)
            return R * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

        max_km = max(
            _hav_km(center_lat, center_lon, lat, lon)
            for lat, lon in zip(all_lats, all_lons)
        )
        radius_m = (max_km + 10) * 1000

        self._road_graph = ox.graph_from_point(
            (center_lat, center_lon), dist=radius_m, network_type="drive"
        )
        print(
            f"Road network ready: "
            f"{len(self._road_graph.nodes):,} nodes, "
            f"{len(self._road_graph.edges):,} edges."
        )
        return self._road_graph

    def _distance_road_network(self, farm_row: pd.Series, refinery_idx: int) -> float:
        """
        Shortest road distance in km via osmnx + NetworkX.
        Falls back to the precomputed straight-line distance if no path exists.
        """
        import osmnx as ox
        import networkx as nx

        G = self._load_road_graph()
        rlat, rlon = self.refinery_locations[refinery_idx]

        orig = ox.nearest_nodes(G, X=farm_row["long"], Y=farm_row["lat"])
        dest = ox.nearest_nodes(G, X=rlon,             Y=rlat)

        try:
            return nx.shortest_path_length(G, orig, dest, weight="length") / 1000.0
        except nx.NetworkXNoPath:
            warnings.warn(
                f"No road path for farm at "
                f"({farm_row['lat']:.5f}, {farm_row['long']:.5f}). "
                "Falling back to precomputed straight-line distance."
            )
            return float(farm_row["distance_to_refinery_km"])

    def _get_distance(self, farm_row: pd.Series, refinery_idx: int) -> float:
        """Dispatch to the selected distance method."""
        if self.distance_method == "straight_line":
            return self._distance_straight_line(farm_row)
        elif self.distance_method == "tortuosity":
            return self._distance_tortuosity(farm_row)
        elif self.distance_method == "road_network":
            return self._distance_road_network(farm_row, refinery_idx)
        elif self.distance_method == "neural_network":
            raise RuntimeError(
                "Neural-network distances should be computed in batch."
            )

     # Neural network prediction
    def _distance_neural_network_batch(
        self,
        farms_df: pd.DataFrame,
        refinery_idx: int,
    ):
        """
        Predict transport distances using trained neural network.

        Returns
        -------
        np.ndarray
            Predicted distances in km.
        """

        import torch

        rlat, rlon = self.refinery_locations[refinery_idx]

        # Convert farm coordinates to EPSG:5070
        farm_xy = self._to_epsg5070(
            farms_df["lat"].values,
            farms_df["long"].values,
        )

        # Convert refinery coordinates to EPSG:5070
        refinery_xy = self._to_epsg5070(
            np.full(len(farms_df), rlat),
            np.full(len(farms_df), rlon),
        )

        # Build feature matrix
        # [farm_x, farm_y, refinery_x, refinery_y]
        X = np.column_stack([
            farm_xy[:, 0],
            farm_xy[:, 1],
            refinery_xy[:, 0],
            refinery_xy[:, 1],
        ])

        # Apply SAME scaler used during training
        X_scaled = self._nn_scaler.transform(X)

        # Convert to tensor
        X_tensor = torch.tensor(
            X_scaled,
            dtype=torch.float32,
            device=self.nn_device,
        )

        # Predict
        with torch.no_grad():
            preds = self._nn_model(X_tensor)

        preds = preds.cpu().numpy().flatten()

        # Optional safety clamp
        # Avoid negative predicted distances
        preds = np.maximum(preds, 0)

        return preds
    # ======================================================================
    # Public API
    # ======================================================================

    def compute(self) -> list:
        """
        Run the full pipeline for all refineries.

        Adds the following columns to each farm DataFrame:
            state_used              — state name used for cost lookup
            computed_distance_km    — distance per selected method
            cost_per_km_per_trip    — raw value from the CSV ($/km/trip)
            cost_per_Mg_per_km      — cost_per_km_per_trip / transport_capacity_Mg
            transport_cost_total    — cost_per_Mg_per_km × biomass × distance  ($)
            trips_needed            — ceil(biomass / capacity)

        Returns
        -------
        list of pd.DataFrame, one per refinery (also stored as self.results).
        """
        results = []

        for i, (farms_df, (rlat, rlon)) in enumerate(
            zip(self.harvested_farms, self.refinery_locations)
        ):
            df = farms_df.copy()

            # ── state resolution ──────────────────────────────────────────
            if "state" in df.columns:
                df["state_used"] = df["state"].astype(str)
            else:
                # print(
                #     f"Refinery {i}: no 'state' column — "
                #     f"reverse-geocoding {len(df)} farm locations…"
                # )
                df["state_used"] = self._lookup_states_from_shapefile(df)

            # ── distances ─────────────────────────────────────────────────
            if self.distance_method == "neural_network":
                df["computed_distance_km"] = (
                    self._distance_neural_network_batch(df, i)
                )
            else:
                df["computed_distance_km"] = df.apply(
                    lambda row: self._get_distance(row, i),
                    axis=1,
                )

            # Record the tortuosity factor actually used per farm (informational)
            if self.distance_method == "tortuosity":
                df["tortuosity_factor_used"] = df["state_used"].apply(
                    self._get_tortuosity_factor
                )

            # ── costs ─────────────────────────────────────────────────────
            #   cost_per_km_per_trip  →  divide by capacity  →  $/Mg/km
            #   transport_cost_total  =  $/Mg/km  ×  Mg  ×  km
            df["cost_per_km_per_trip"] = df["state_used"].apply(self._get_cost_for_state)

            self.transport_capacity_Mg = self.bale_capacity_per_truck * self.bale_mass_Mg

            df["cost_per_Mg_per_km"]   = (
                df["cost_per_km_per_trip"] / self.transport_capacity_Mg
            )

            df["trips_needed"] = df.apply(
                lambda row: self._trips_needed(row["biomass_supplied_Mg"], row),
                axis=1
            )

            df["transport_cost_total"] = (
                df["cost_per_Mg_per_km"]
                * df["biomass_supplied_Mg"]
                * df["computed_distance_km"]
                *2 #round-trip
            )

            # df["transport_cost_total"] = (
            #     df["trips_needed"]
            #     * df["cost_per_km_per_trip"]
            #     * df["computed_distance_km"]
            #     *2 #round-trip
            # )

            df["cost_per_Mg"] = (
                df["transport_cost_total"] / df["biomass_supplied_Mg"]
            )

            # ── trips ─────────────────────────────────────────────────────
            # df["trips_needed"] = df["biomass_supplied_Mg"].apply(self._trips_needed)

            results.append(df)

        self.results = results
        return results

    def summary(self) -> pd.DataFrame:
        """
        Per-refinery summary of distances, costs, and trip counts.

        Raises RuntimeError if compute() has not been called yet.
        """
        if self.results is None:
            raise RuntimeError("Call compute() before summary().")

        rows = []
        for i, (df, (rlat, rlon)) in enumerate(
            zip(self.results, self.refinery_locations)
        ):
            total_biomass = df["biomass_supplied_Mg"].sum()
            total_cost    = df["transport_cost_total"].sum()

            rows.append({
                "refinery_idx"              : i,
                "refinery_lat"              : rlat,
                "refinery_lon"              : rlon,
                "num_farms"                 : len(df),
                "total_biomass_supplied_Mg" : round(total_biomass, 3),
                "amount_to_transport_Mg"    : self.amount_transported_Mg[i],
                "trips_for_total_amount"    : math.ceil((self.amount_transported_Mg[i] / self.bale_mass_Mg)/ self.bale_capacity_per_truck),
                "mean_distance_km"          : round(df["computed_distance_km"].mean(), 3),
                "max_distance_km"           : round(df["computed_distance_km"].max(), 3),
                "total_transport_cost_USD"  : round(total_cost, 2),
                "cost_per_Mg_USD"           : (
                    round(total_cost / total_biomass, 4) if total_biomass else 0
                ),
                "distance_method"           : self.distance_method,
            })

        return pd.DataFrame(rows)

    def get_farm_results(self, refinery_idx: int = 0) -> pd.DataFrame:
        """
        Detailed per-farm result DataFrame for a given refinery index.

        Parameters
        ----------
        refinery_idx : int
            Index into harvested_farms / refinery_locations.

        Raises RuntimeError if compute() has not been called.
        Raises IndexError if refinery_idx is out of range.
        """
        if self.results is None:
            raise RuntimeError("Call compute() before get_farm_results().")
        if refinery_idx >= len(self.results):
            raise IndexError(
                f"refinery_idx {refinery_idx} out of range "
                f"(have {len(self.results)} refineries)."
            )
        return self.results[refinery_idx]

    def get_total_transport_costs(self) -> list:
        """
        Returns total transport cost per refinery as a list.

        Returns
        -------
        list of float
            Total transport costs [USD] for each refinery.
        """

        if self.results is None:
            raise RuntimeError(
                "Call compute() before get_total_transport_costs()."
            )

        return [
            df["transport_cost_total"].sum()
            for df in self.results
        ]

    def __repr__(self) -> str:
        status = "computed" if self.results is not None else "not computed"
        return (
            f"BiomassTransport("
            f"refineries={len(self.refinery_locations)}, "
            f"method='{self.distance_method}', "
            # f"capacity={self.transport_capacity_Mg} Mg, "
            f"capacity={self.bale_capacity_per_truck} bales, "
            f"status={status})"
        )

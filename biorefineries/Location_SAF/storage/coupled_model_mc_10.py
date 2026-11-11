from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pub_style import apply_pub_style, style_axes, sync_twins, save_fig
from biomass_transport_model_02 import BiomassTransport, load_transport_costs_csv
from farm_allocation_model_07 import FeedstockAllocationModel
from refinery_model import RefineryModel
from storage_model_05 import (Plant, StorageTech, Storage, StorageScenario, Simulation,
                                SimulationRunner, StorageConfig)

# =============================================================================
# CONFIGURATION (Constants preserved exactly)
# =============================================================================
REFINERY_SIZE     = 80          # MMgal/yr
YIELD_ACCOUNTING_VERSION = "fixed_contract_and_planned_capacity_v1"
MANUSCRIPT_CORRECTIONS_VERSION = "partial_last_farm_no_barn_storage_conversion_v1"
REFINERY_LOCATION = [(35.4676, -97.5164)]       # OK
MOISTURE_CONTENT  = 0.20
CONVERSION        = 113.97      # gal / dry Mg miscanthus

_ROOT = Path(__file__).resolve().parent
TRANSPORT_COSTS_CSV = str(_ROOT / "Truck_cost_per_state.csv")
SHAPEFILE_PATH      = str(_ROOT / "GIS_data" / "US_map_for_boundary.shp")
TORTUOSITY_XLSX     = str(_ROOT / "GIS_data" / "Tortuosity_per_state.xlsx")
DISTANCE_METHOD     = "tortuosity"
TORTUOSITY_FACTOR   = 1.25


MIN_PENALTY = 0.0
MID_PENALTY = 0.3
MAX_PENALTY = 1.0

MIN_LAND_DENSITY = 0.02
MID_LAND_DENSITY = 0.10
MAX_LAND_DENSITY = 0.20

SALVAGE_PRICE = 0 # price to sell extra purchased biomass (between 0 and feedstock price) in USD_per_Mg
# =============================================================================
# STATIC LAYERS FOR BiomassTransport, LOADED ONCE PER PROCESS
#
# BiomassTransport() used to re-read+reproject SHAPEFILE_PATH and re-read
# TRANSPORT_COSTS_CSV on every single instantiation (thousands of times
# across a Monte Carlo run). Load them once here and pass them in via
# us_states_gdf= / cost_lookup= instead.
# =============================================================================
_US_STATES_GDF = gpd.read_file(SHAPEFILE_PATH)[["NAME", "geometry"]].to_crs(epsg=4326)
_COST_LOOKUP   = load_transport_costs_csv(TRANSPORT_COSTS_CSV)

START_HARVEST_DATE = "2025-11-01"
HARVEST_DAYS = 150

YIELD_BONUS = {
    10: {"min": 0.00, "peak": 0.15, "max": 0.30},
    11: {"min": 0.00, "peak": 0.13, "max": 0.25},
    12: {"min": 0.00, "peak": 0.07, "max": 0.20},
     1: {"min": 0.00, "peak": 0.05, "max": 0.20},
     2: {"min": 0.00, "peak": 0.02, "max": 0.10},
     3: {"min": 0.00, "peak": 0.00, "max": 0.00},
}

MIN_YIELD_VARIATION = 0.75
MAX_YIELD_VARIATION = 1.25
MID_YIELD_VARIATION = 1.00

BALE_PARAMS = {
    "round":  {"bale_capacity_per_truck": 42, "bale_mass_Mg": 0.27},
    "square": {"bale_capacity_per_truck": 36, "bale_mass_Mg": 0.35},
}

DML_RANGES   = {"uncovered": (0.15, 0.45), "tarped": (0.05, 0.15), "barn": (0.02, 0.08)}
ALPHA_RANGES = {"uncovered": (1.00, 1.25), "tarped": (1.00, 1.15), "barn": (1.00, 1.10)}

STORAGE_UNIT_COSTS = {
    "uncovered": {"round": {"fixed": 0.000, "variable": 0.000}, "square": {"fixed": 0.000, "variable": 0.000}},
    "tarped":    {"round": {"fixed": 0.0066, "variable": 0.0006}, "square": {"fixed": 0.0037, "variable": 0.0004}},
    "barn":      {"round": {"fixed": 0.10,   "variable": 0.000},  "square": {"fixed": 0.09,  "variable": 0.000}},
}

STORAGE_COST_FACTOR_RANGE = (0.70, 1.30)   # ±30% uncertainty

Storage_control_costs = {
    "uncovered": {"round": {"fixed": 0.000, "variable": 0.000}, "square": {"fixed": 0.000, "variable": 0.000}},
    "tarped":    {"round": {"fixed": 0.0,   "variable": 0.0},   "square": {"fixed": 0.0,   "variable": 0.0}},
    "barn":      {"round": {"fixed": 0.0,   "variable": 0.000}, "square": {"fixed": 0.00,  "variable": 0.003}},
    # $/kg/day, applied only to barn/square early-harvest bales for first 42 days
}

HARVEST_DELTA_MEAN = 110.31 - 115.51
HARVEST_DELTA_STD  = np.sqrt(2.20**2 + 2.24**2) / np.sqrt(2)
EARLY_HARVEST_MONTHS = {10, 11}

MSP_SLOPE     = -0.0108
MSP_RESID_STD =  0.0161

STORAGE_CONV_DELTA_MEAN = -1.97
STORAGE_CONV_DELTA_STD  =  3.55

INCLUDE_STORAGE_CONVERSION_EFFECTS = True
INCLUDE_HARVEST_TIMING_EFFECTS     = True

COMBOS = [
    ("tarped", "square", "FIFO"),
    ("tarped", "square", "LIFO"),
    ("tarped", "round",  "FIFO"),
    ("tarped", "round",  "LIFO"),
]

CACHE_RESOLUTION_KG = 1000

### ==== helper
def _sample_triangular(rng, left, mode, right, N, decimals=2):
    if left == mode == right:
        return np.full(N, left)

    samples = rng.triangular(left, mode, right, N)
    return np.round(samples, decimals)


def _sample_uniform(rng, low, high, N, decimals=2):
    if low == high:
        return np.full(N, low)

    samples = rng.uniform(low, high, N)
    return np.round(samples, decimals)

### ===================================
# Random inputs - master file

def generate_random_inputs(
    N=500,
    seed=42,
):
    """
    Generate one master set of Monte Carlo random inputs.

    These inputs can be saved and reused across different scenarios so that
    scenario comparisons are paired by Monte Carlo sample.
    """

    rng = np.random.default_rng(seed)

    # ---------------------------------------------------------
    # Common inputs
    # ---------------------------------------------------------

    inputs = pd.DataFrame({
        "sample": np.arange(N),

        "penalty": _sample_triangular(
            rng,
            MIN_PENALTY,
            MID_PENALTY,
            MAX_PENALTY,
            N
        ),

        "land_density": _sample_triangular(
            rng,
            MIN_LAND_DENSITY,
            MID_LAND_DENSITY,
            MAX_LAND_DENSITY,
            N
        ),

        "storage_cost_factor": _sample_uniform(
            rng,
            *STORAGE_COST_FACTOR_RANGE,
            N
        ),

        "storage_conv_delta": rng.normal(
            STORAGE_CONV_DELTA_MEAN,
            STORAGE_CONV_DELTA_STD,
            N
        ),

        "msp_resid": rng.normal(
            0,
            MSP_RESID_STD,
            N
        ),

        "yield_variation":_sample_triangular(
            rng,
            MIN_YIELD_VARIATION,
            MID_YIELD_VARIATION,
            MAX_YIELD_VARIATION,
            N
        ),
    })

    # ---------------------------------------------------------
    # DML and alpha
    # Only generate storage types actually used in COMBOS
    # ---------------------------------------------------------

    unique_types = list(dict.fromkeys(c[0] for c in COMBOS))

    for storage_type in unique_types:

        inputs[f"dml_{storage_type}"] = _sample_uniform(
            rng,
            *DML_RANGES[storage_type],
            N
        )

        inputs[f"alpha_{storage_type}"] = _sample_uniform(
            rng,
            *ALPHA_RANGES[storage_type],
            N
        )

    return inputs

## HARVEST INPUTS
def add_harvest_inputs(
    random_inputs,
    start_date,
    harvest_days,
    seed=1234,
):
    """
    Add harvest-window-specific random inputs to a master random-input
    DataFrame.

    The original DataFrame is not modified.
    """

    rng = np.random.default_rng(seed)

    N = len(random_inputs)

    harvest_conv_delta, yield_bonus, month_weights = (
        sample_harvest_window_vectorized(
            start_date,
            harvest_days,
            N,
            rng
        )
    )

    out = random_inputs.copy()

    out["harvest_conv_delta"] = harvest_conv_delta
    out["yield_bonus"] = yield_bonus

    for m in [10, 11, 12, 1, 2, 3, 4, 5]:
        out[f"month_{m}_frac"] = [
            mw.get(m, 0.0)
            for mw in month_weights
        ]

    return out

# =============================================================================
# HIGH-SPEED VECTORIZED HELPERS
# =============================================================================

def _make_plant(harvest_days, conversion):
    demand_kg = (REFINERY_SIZE * 1e6 / conversion) / (1.0 - MOISTURE_CONTENT) * 1000.0
    return Plant(annual_demand_kg=demand_kg, operating_days=350, harvest_days=harvest_days)


def _make_storage_scenario(storage_type, bale_shape, rule, dml_max, alpha, max_days, storage_cost_factor=1.0, yield_multiplier=1.0):
    cfg = STORAGE_UNIT_COSTS[storage_type][bale_shape]
    return StorageScenario(
        name=f"{storage_type}_{bale_shape}_{rule}",
        dml_max=dml_max,
        ql_max=0.0,
        max_days=int(max_days),
        alpha=alpha,
        yield_multiplier=yield_multiplier,
        rule=rule,
        bale_shape=bale_shape,
        storage_type=storage_type,
        fixed_cost_per_kg_capacity=cfg["fixed"] * storage_cost_factor,
        variable_cost_per_kg_throughput=cfg["variable"] * storage_cost_factor,
    )


# def _run_transport(farm_model, bale_shape, supplied_biomass_kg, cache):
#     supply_bucket = round(supplied_biomass_kg / CACHE_RESOLUTION_KG) * CACHE_RESOLUTION_KG
#     cache_key = (bale_shape, supply_bucket)
#     if cache_key in cache:
#         return cache[cache_key]

#     bp = BALE_PARAMS[bale_shape]
#     bt = BiomassTransport(
#         refinery_locations       = REFINERY_LOCATION,
#         harvested_farms          = farm_model.harvested_farms,
#         transport_costs_csv      = TRANSPORT_COSTS_CSV,
#         amount_transported_Mg    = supplied_biomass_kg / 1000.0,
#         bale_capacity_per_truck  = bp["bale_capacity_per_truck"],
#         bale_mass_Mg             = bp["bale_mass_Mg"],
#         shapefile_path           = SHAPEFILE_PATH,
#         state_col                = "NAME",
#         use_per_state_tortuosity = False,
#         tortuosity_xlsx          = TORTUOSITY_XLSX,
#         distance_method          = DISTANCE_METHOD,
#         tortuosity_factor        = TORTUOSITY_FACTOR,
#     )
#     bt.compute()
#     raw = bt.get_total_transport_costs()
#     result = float(sum(raw) if isinstance(raw, (list, tuple, np.ndarray)) else raw)
#     cache[cache_key] = result
#     return result

# # =============================================================================
# # AFTER (fixed): cache key includes the inputs that actually determine
# # which farms get selected, so cross-sample collisions can't happen
# # =============================================================================
# def _run_transport(farm_model, bale_shape, supplied_biomass_kg, land, penalty, cache):
#     supply_bucket = round(supplied_biomass_kg / CACHE_RESOLUTION_KG) * CACHE_RESOLUTION_KG
#     # land_density and penalty are what FeedstockAllocationModel uses to choose
#     # which farms are contracted, so they must be part of the cache identity.
#     # Rounding to 6 dp just guards against float noise from array indexing;
#     # it does NOT create cross-sample collisions since land/penalty are drawn
#     # from continuous distributions (rng.triangular) and are effectively
#     # unique per sample.
#     cache_key = (bale_shape, round(land, 6), round(penalty, 6), supply_bucket)
#     if cache_key in cache:
#         return cache[cache_key]

#     bp = BALE_PARAMS[bale_shape]
#     bt = BiomassTransport(
#         refinery_locations       = REFINERY_LOCATION,
#         harvested_farms          = farm_model.harvested_farms,
#         transport_costs_csv      = TRANSPORT_COSTS_CSV,   # unused when cost_lookup is passed; kept for the required-arg signature
#         us_states_gdf            = _US_STATES_GDF,        # skip re-reading/reprojecting the shapefile every call
#         cost_lookup              = _COST_LOOKUP,           # skip re-reading the cost CSV every call
#         amount_transported_Mg    = supplied_biomass_kg / 1000.0,
#         bale_capacity_per_truck  = bp["bale_capacity_per_truck"],
#         bale_mass_Mg             = bp["bale_mass_Mg"],
#         shapefile_path           = SHAPEFILE_PATH,
#         state_col                = "NAME",
#         use_per_state_tortuosity = False,
#         tortuosity_xlsx          = TORTUOSITY_XLSX,
#         distance_method          = DISTANCE_METHOD,
#         tortuosity_factor        = TORTUOSITY_FACTOR,
#     )
#     bt.compute()
#     raw = bt.get_total_transport_costs()
#     result = float(sum(raw) if isinstance(raw, (list, tuple, np.ndarray)) else raw)
#     cache[cache_key] = result
#     return result

def _run_transport(farm_model, bale_shape, supplied_biomass_kg, land, penalty, cache):
    bp = BALE_PARAMS[bale_shape]

    bt = BiomassTransport(
        refinery_locations=REFINERY_LOCATION,
        harvested_farms=farm_model.harvested_farms,
        transport_costs_csv=TRANSPORT_COSTS_CSV,
        us_states_gdf=_US_STATES_GDF,
        cost_lookup=_COST_LOOKUP,
        amount_transported_Mg=supplied_biomass_kg / 1000.0,
        bale_capacity_per_truck=bp["bale_capacity_per_truck"],
        bale_mass_Mg=bp["bale_mass_Mg"],
        shapefile_path=SHAPEFILE_PATH,
        state_col="NAME",
        use_per_state_tortuosity=False,
        tortuosity_xlsx=TORTUOSITY_XLSX,
        distance_method=DISTANCE_METHOD,
        tortuosity_factor=TORTUOSITY_FACTOR,
    )

    bt.compute()
    raw = bt.get_total_transport_costs()

    return float(
        sum(raw) if isinstance(raw, (list, tuple, np.ndarray)) else raw
    )

def sample_harvest_window_vectorized(start_date, harvest_days, N, rng):
    """
    Vectorized harvest window calculations. Evaluates the 150-day window matrix
    simultaneously instead of iterating row-by-row.
    """
    start = pd.Timestamp(start_date)
    dates = [start + pd.Timedelta(days=d) for d in range(harvest_days)]
    months = np.array([d.month for d in dates], dtype=np.int32)

    unique_months, counts = np.unique(months, return_counts=True)
    month_weights = {int(m): float(c) / harvest_days for m, c in zip(unique_months, counts)}

    means = np.zeros_like(months, dtype=np.float64)
    means[np.isin(months, list(EARLY_HARVEST_MONTHS))] = HARVEST_DELTA_MEAN
    means[months == 12] = HARVEST_DELTA_MEAN / 2.0

    # Generate the same number of draws as before, preserving the
    # random-number sequence used subsequently for yield bonuses.
    random_deltas = rng.normal(
        means[:, np.newaxis],
        HARVEST_DELTA_STD,
        size=(harvest_days, N),
    )

    # Retain the original nonpositive conversion effects for
    # October–November and the transitional December period.
    random_deltas = np.minimum(random_deltas, 0.0)

    # Late harvest is the reference: no harvest-conversion adjustment.
    # This includes all months whose specified mean effect is zero.
    reference_days = means == 0.0
    random_deltas[reference_days, :] = 0.0

    mean_deltas = random_deltas.mean(axis=0)

    random_bonuses = np.zeros((harvest_days, N))
    for m in unique_months:
        if m in YIELD_BONUS:
            v = YIELD_BONUS[m]
            if v["min"] != v["max"]:
                idx = np.where(months == m)[0]
                random_bonuses[idx, :] = rng.triangular(v["min"], v["peak"], v["max"], size=(len(idx), N))

    mean_bonuses = random_bonuses.mean(axis=0)

    return mean_deltas, mean_bonuses, [month_weights] * N


# =============================================================================
# COUPLED MONTE CARLO RUNNER
# =============================================================================

def run_coupled_mc(
    N=1000,
    seed=42,
    verbose=True,
    include_storage_conversion_effects=INCLUDE_STORAGE_CONVERSION_EFFECTS,
    include_harvest_timing_effects=INCLUDE_HARVEST_TIMING_EFFECTS,
    max_storage_days=200,
    random_inputs=None,
    scenario="storage",
):
    rng = np.random.default_rng(seed)

    if random_inputs is None:

        # penalty_samples = np.full(N, MAX_PENALTY) if MIN_PENALTY == MAX_PENALTY else rng.triangular(MIN_PENALTY, MID_PENALTY, MAX_PENALTY, N)
        # land_samples    = rng.triangular(MIN_LAND_DENSITY, MID_LAND_DENSITY, MAX_LAND_DENSITY, N)

        # storage_cost_factor_samples = rng.uniform(*STORAGE_COST_FACTOR_RANGE,N,)

        # harvest_conv_delta_raw, yield_bonus_raw, month_weights_list = sample_harvest_window_vectorized(
        #     START_HARVEST_DATE, HARVEST_DAYS, N, rng
        # )

        # unique_types  = list(dict.fromkeys(c[0] for c in COMBOS))
        # dml_samples   = {st: rng.uniform(*DML_RANGES[st],   N) for st in unique_types}
        # alpha_samples = {st: rng.uniform(*ALPHA_RANGES[st], N) for st in unique_types}

        # storage_conv_delta_raw = rng.normal(STORAGE_CONV_DELTA_MEAN, STORAGE_CONV_DELTA_STD, N)
        # msp_resid_raw          = rng.normal(0, MSP_RESID_STD, N)

        # Generate master random inputs
        random_inputs = generate_random_inputs(N=N, seed=seed)

        if include_harvest_timing_effects:
            random_inputs = add_harvest_inputs(random_inputs, START_HARVEST_DATE, HARVEST_DAYS, seed=seed + 1000)

        else:
            random_inputs["harvest_conv_delta"] = 0.0

            random_inputs["yield_bonus"] = 0.0

            for m in [10, 11, 12, 1, 2, 3, 4, 5,]:
                random_inputs[f"month_{m}_frac"] = 0.0

    else:
        # Use supplied random inputs
        random_inputs = (random_inputs.copy().reset_index(drop=True))
        # N comes directly from the DataFrame
        N = len(random_inputs)
        if N == 0:
            raise ValueError(
                "random_inputs contains zero rows.")

    ## VALIDATE INPUTS

    required_columns = [
        "sample",
        "penalty",
        "land_density",
        "storage_cost_factor",
        "storage_conv_delta",
        "msp_resid",
        "yield_variation",
    ]

    if include_harvest_timing_effects:

        required_columns += [
            "harvest_conv_delta",
            "yield_bonus",
        ]

    unique_types = list(dict.fromkeys(c[0]for c in COMBOS))

    if scenario != "base":
        for storage_type in unique_types:
            required_columns += [
                f"dml_{storage_type}",
                f"alpha_{storage_type}",
            ]
    missing = [c for c in required_columns if c not in random_inputs.columns]
    if missing:
        raise ValueError("random_inputs is missing required columns: "f"{missing}")

    # CONVERT INPUTS TO NUMPY ARRAYS
    sample_ids = (random_inputs[ "sample"].to_numpy())

    penalty_samples = (random_inputs["penalty"].to_numpy())

    land_samples = (random_inputs["land_density"].to_numpy())

    storage_cost_factor_samples = (random_inputs["storage_cost_factor"].to_numpy())

    storage_conv_delta_raw = (random_inputs["storage_conv_delta"].to_numpy())

    msp_resid_raw = (random_inputs["msp_resid"].to_numpy())


    # HARVEST INPUTS
    if include_harvest_timing_effects:
        harvest_conv_delta_raw = (random_inputs["harvest_conv_delta"].to_numpy())
        yield_bonus_raw = (random_inputs["yield_bonus"].to_numpy())
    else:
        harvest_conv_delta_raw = (np.zeros(N))
        yield_bonus_raw = (np.zeros(N))

    # DML AND ALPHA INPUTS
    if scenario == "base":
        # BASE = no storage degradation
        dml_samples = {st: np.zeros(N)for st in unique_types}
        alpha_samples = {st: np.ones(N)for st in unique_types}
        yield_variation_samples = np.ones(N)

    else:
        dml_samples = {st: random_inputs[f"dml_{st}"].to_numpy() for st in unique_types}

        alpha_samples = {st: random_inputs[f"alpha_{st}"].to_numpy()for st in unique_types}

        yield_variation_samples = random_inputs["yield_variation"].to_numpy()

    if not np.isfinite(yield_variation_samples).all() or (yield_variation_samples < 0).any():
        raise ValueError("yield_variation must be finite and nonnegative")



    # MONTH WEIGHTS
    month_list = [10, 11, 12, 1, 2,3, 4,5,]
    month_cols = [f"month_{m}_frac"for m in month_list]

    if all(
        c in random_inputs.columns
        for c in month_cols
    ):

        month_weights_list = [
            {
                m: row[
                    f"month_{m}_frac"
                ]

                for m in month_list
            }

            for _, row
            in random_inputs.iterrows()
        ]

    else:

        month_weights_list = [
            {}
            for _ in range(N)
        ]

    _transport_cache = {}
    cache_hits = 0
    records = []

    for i in range(N):
        penalty = penalty_samples[i]
        land    = land_samples[i]
        storage_cost_factor = storage_cost_factor_samples[i]

        yield_bonus      = float(yield_bonus_raw[i])       if include_harvest_timing_effects else 0.0
        harvest_cv_delta = float(harvest_conv_delta_raw[i]) if include_harvest_timing_effects else 0.0
        conversion_sample = CONVERSION + harvest_cv_delta

        plant = _make_plant(HARVEST_DAYS, conversion_sample)

        # Pass start_date and harvest_days so the runner can compute the
        # early-harvest cutoff for barn/square scenarios automatically.
        runner = SimulationRunner(
            plant,
            start_date=START_HARVEST_DATE,
            harvest_days=HARVEST_DAYS,
        )

        # ── Farm allocation model, built ONCE per sample ────────────────────
        # land/penalty/seed are identical across all 12 combos for this sample,
        # so a single FeedstockAllocationModel instance can be reused for all
        # of them: calculate_location_parameters() overwrites self.results and
        # self.harvested_farms in place each call and reuses self._distance_cache
        # (populated on the first call) instead of recomputing the farm-distance
        # matrix from scratch 12x per sample.
        farm_model = FeedstockAllocationModel(
            feedstock                = "miscanthus",
            sizes_of_biorefineries   = REFINERY_SIZE,
            suitable_land_density    = land,
            refinery_locations       = REFINERY_LOCATION,
            seed                     = 1,
            moisture_content         = MOISTURE_CONTENT,
            harvested_biomass        = 0.0,   # set per-combo below
            delivered_biomass        = 0.0,   # set per-combo below
            penalty_for_unmet_demand = penalty,
            yield_bonus              = 0.0,   # set per-combo below
            conversion_efficiency    = CONVERSION,
        )

        for storage_type, bale_shape, rule in COMBOS:
            dml   = dml_samples[storage_type][i]
            alpha = alpha_samples[storage_type][i]
            annual_yield_multiplier = float(yield_variation_samples[i])

            # 1. Storage simulation
            sc   = _make_storage_scenario(storage_type, bale_shape, rule, dml, alpha, max_storage_days, storage_cost_factor,
                                          yield_multiplier=annual_yield_multiplier)
            stor = runner.run(sc)

            # 2. Storage conversion penalty
            # Outdoor composition response only; barn keeps physical DML/costs
            # but has no storage-related conversion adjustment.
            if include_storage_conversion_effects and storage_type != "barn":
                mean_age         = stor["mean_dispatch_age"]
                storage_cv_delta = storage_conv_delta_raw[i] * (1.0 - np.exp(-mean_age / 270.0))
            else:
                storage_cv_delta = 0.0

            total_conversion = conversion_sample + storage_cv_delta
            delta_conv       = total_conversion - CONVERSION

            # 3. Farm allocation — update this combo's inputs on the shared
            # instance and re-run (see comment above the instance creation).
            # Contract farms and pay for the planned amount before the annual
            # yield shock. Seasonal harvest assumptions remain unchanged.
            farm_model.harvested_biomass = stor["contracted_biomass"] / 1000.0
            farm_model.delivered_biomass = stor["supplied_biomass"]  / 1000.0
            farm_model.yield_bonus       = yield_bonus
            farm_model.calculate_location_parameters()
            farm = farm_model.results

            # 4. Allocate storage biomass to the actual refinery demand.
            # ``supplied_biomass`` is the regular annual dispatch from storage.
            # End-of-year salvage plus surviving bonus biomass is a separate,
            # finite inventory pool: it must be used at most once.
            real_demand_wet_Mg     = (REFINERY_SIZE * 1e6 / total_conversion) / (1.0 - MOISTURE_CONTENT)
            baseline_demand_wet_Mg = (REFINERY_SIZE * 1e6 / CONVERSION)       / (1.0 - MOISTURE_CONTENT)
            conv_gap_Mg            = real_demand_wet_Mg - baseline_demand_wet_Mg
            bonus_yield_Mg         = float(np.asarray(farm["bonus_biomass_Mg"] * annual_yield_multiplier * (1 - stor["avg_dml"])).item())
            harvested_Mg           = float(np.asarray(farm["delivered_biomass_Mg"]).item()) * annual_yield_multiplier
            base_storage_supply_Mg = float(stor["supplied_biomass"] / 1000.0)
            extra_stored_Mg        = float(stor["salvage_biomass_Mg"] + bonus_yield_Mg)

            # In a late-harvest case the conversion gap can be negative.  Do
            # not use max(conv_gap, 0): allocate against *real* demand instead.
            base_to_refinery_Mg = min(base_storage_supply_Mg, real_demand_wet_Mg)
            discarded_Mg        = max(base_storage_supply_Mg - real_demand_wet_Mg, 0.0)
            remaining_demand_Mg = real_demand_wet_Mg - base_to_refinery_Mg
            extra_used_Mg       = min(extra_stored_Mg, remaining_demand_Mg)
            supplied_Mg         = base_to_refinery_Mg + extra_used_Mg
            extra_spot_Mg       = remaining_demand_Mg - extra_used_Mg
            salvaged_biomass_Mg = extra_stored_Mg - extra_used_Mg
            overpurchase_Mg        = (alpha-1) * baseline_demand_wet_Mg
            overpurchase_Mg        = float(np.asarray(overpurchase_Mg).item())

            # Farm-to-storage balance.  Bonus DML is already captured here:
            # it must not be added a second time.
            DML_Mg = harvested_Mg - supplied_Mg - salvaged_biomass_Mg - discarded_Mg

            # 5. Transport (cached): transport exactly the biomass supplied
            # from storage, including any end-of-year inventory that is used.
            before = len(_transport_cache)
            transport_cost = _run_transport(farm_model, bale_shape, supplied_Mg * 1000.0,
                                            land, penalty,
                                            _transport_cache)
            if len(_transport_cache) == before:
                cache_hits += 1



            unit_price_raw    = farm["feedstock_unit_price_contract"]
            unit_price_per_Mg = float(unit_price_raw.item() if hasattr(unit_price_raw, "item") else unit_price_raw) * 1000.0


            # 5b. Refinery model pricing
            # extra_cost_raw    = farm["extra_biomass_cost_total"]
            # extra_cost_M      = float(extra_cost_raw.item() if hasattr(extra_cost_raw, "item") else extra_cost_raw)

            extra_cost = extra_spot_Mg * unit_price_per_Mg * (1.0 + penalty)
            extra_cost = float(extra_cost.item() if hasattr(extra_cost, "item") else extra_cost)

            # Reporting only: a negative conversion gap has no positive
            # procurement cost.  Spot procurement can cover both a baseline
            # shortage and a positive conversion requirement.
            conv_gap_cost = min(extra_spot_Mg, max(conv_gap_Mg, 0.0)) * unit_price_per_Mg * (1.0 + penalty)

            # salvage price of discarded biomass at the end of the year
            salvage_price = SALVAGE_PRICE * salvaged_biomass_Mg

            contracted_cost_raw = farm["Total_contracted_feedstock_cost"] # includes yield variation, and bonus yield for harvest timing
            contracted_cost   = float(contracted_cost_raw.item() if hasattr(contracted_cost_raw, "item") else contracted_cost_raw)

            # 5c. Control cost: barn/square early-harvest bales, first 42 days only.
            # The storage model returns early_controlled_kg_days = 0 for all other
            # combos (tracking is disabled), so this expression is always safe.
            if storage_type == "barn" and bale_shape == "square" and include_harvest_timing_effects:
                ctrl_var            = Storage_control_costs["barn"]["square"]["variable"]
                control_extra_cost  = stor["early_controlled_kg_days"] * ctrl_var
            else:
                control_extra_cost  = 0.0

            control_extra_cost = float(control_extra_cost.item() if hasattr(control_extra_cost, "item") else control_extra_cost)

            feedstock_cost = (contracted_cost + extra_cost - salvage_price)  #+ conv_gap_cost
            feedstock_cost = float(feedstock_cost.item() if hasattr(feedstock_cost, "item") else feedstock_cost)
            storage_cost   = float(stor["total_storage_cost"]) + control_extra_cost

            refinery = RefineryModel(
                refinery_locations    = REFINERY_LOCATION,
                size_of_refineries    = REFINERY_SIZE,
                feedstock_costs       = feedstock_cost,
                storage_costs         = storage_cost,
                transport_costs       = transport_cost,
                product               = "Ethanol",
                moisture_content      = MOISTURE_CONTENT,
                conversion_efficiency = CONVERSION,
            )
            refinery.compute()

            msp_raw = refinery.results["MSP"]
            msp     = float(msp_raw.item() if hasattr(msp_raw, "item") else msp_raw)
            msp_adj = msp + MSP_SLOPE * delta_conv + msp_resid_raw[i]

            # 6. Accumulate record
            total_system      = storage_cost + feedstock_cost + transport_cost
            # storage_demand_Mg = (stor["supplied_biomass"] + stor["unmet_demand"]) / 1000.0

            # supplied_biomass_total = stor["supplied_biomass"]
            # unmet_demand_total     = stor["unmet_demand"]
            # sl_denom               = supplied_biomass_total + unmet_demand_total
            # Storage service level excludes spot purchases and does not count
            # unused inventory twice.  With spot purchases, system service is
            # 100% whenever real demand is positive.
            service_level_val      = 100.0 * supplied_Mg / real_demand_wet_Mg if real_demand_wet_Mg > 0 else 0.0

            contracted_biomass_raw = farm["contracted_biomass_Mg"]
            contracted_Mg_val      = float(contracted_biomass_raw.item() if hasattr(contracted_biomass_raw, "item") else contracted_biomass_raw)

            farm_count_raw = farm["farm_count"]
            farm_count_val = int(farm_count_raw.item() if hasattr(farm_count_raw, "item") else farm_count_raw)

            avg_dml = stor["avg_dml"]
            avg_qual = stor["avg_quality"]


            records.append({
                "sample": sample_ids[i], "scenario": scenario,
                "storage_type": storage_type, "bale_shape": bale_shape, "rule": rule,
                "combo": f"{storage_type}/{bale_shape}/{rule}",
                "flag_storage_conv": include_storage_conversion_effects,
                "flag_harvest_eff":  include_harvest_timing_effects,
                "storage_cost_factor": storage_cost_factor,
                "dml_max": dml, "alpha": alpha, "max_days": max_storage_days,
                "alpha_nominal": alpha,
                "alpha_effective": stor["alpha_effective"],
                "annual_yield_multiplier": annual_yield_multiplier,
                "yield_accounting_version": YIELD_ACCOUNTING_VERSION,
                "planned_contract_Mg": stor["contracted_biomass"] / 1000.0,
                "storage_harvested_Mg": stor["harvested_biomass"] / 1000.0,
                "fixed_storage_capacity_Mg": stor["fixed_cost_capacity_kg"] / 1000.0,
                "realized_peak_storage_Mg": stor["max_storage_required"] / 1000.0,
                "storage_capacity_excess_Mg": stor["capacity_excess_kg"] / 1000.0,
                "fixed_storage_cost_M": stor["fixed_storage_cost"] / 1e6,
                "variable_storage_cost_M": stor["variable_storage_cost"] / 1e6,
                "avg_dml": avg_dml, "avg_quality": avg_qual,
                "penalty": penalty, "land_density": land,
                "yield_bonus": yield_bonus, "harvest_cv_delta": harvest_cv_delta,
                "storage_cv_delta": storage_cv_delta,
                "conversion_sample": conversion_sample, "conversion_total": total_conversion,
                "delta_conversion": delta_conv,
                "storage_cost_M":    storage_cost / 1e6,
                "control_extra_cost_M": control_extra_cost / 1e6,   # if the square bale needs to be controled (fire ignition risk)
                "early_controlled_kg_days": stor["early_controlled_kg_days"],  #  (raw, for QA)
                "service_level": service_level_val,
                "feedstock_unit_price_Mg": unit_price_per_Mg, # feedstock cost per Mg at farm-gate
                "unmet_Mg":       stor["unmet_demand"]     / 1000.0,
                "harvested_Mg":   harvested_Mg, # Actual biomass, including yield bonus or yield loss
                "bonus_biomass_Mg": bonus_yield_Mg, # bonus yield due to harvest timing effects with decay
                "extra_stored_Mg": extra_stored_Mg,# biomass in storage at the end of the year
                "extra_stored_used_Mg": extra_used_Mg,
                "base_storage_supply_Mg": base_storage_supply_Mg,
                # "supplied_Mg":    supplied_biomass_total / 1000.0,
                "supplied_Mg":    supplied_Mg  , # amount delivered to the refinery after DML
                "overpurchase_Mg": overpurchase_Mg, # alpha * baseline demand
                "mean_dispatch_age": stor["mean_dispatch_age"],
                "conv_gap_Mg": conv_gap_Mg, # only for reporting (this is included in the extra_Mg
                "conv_gap_cost_M": conv_gap_cost / 1e6,
                "extra_spot_Mg": extra_spot_Mg, # extra biomass purchased on Spot market
                "extra_cost_M":      extra_cost      / 1e6,
                "real_demand_wet_Mg":real_demand_wet_Mg, # consumed biomass by the refinery
                "baseline_demand_wet_Mg": baseline_demand_wet_Mg, # baseline demand with baseline conversion efficiency
                "discarded_Mg": discarded_Mg, # amount discarded due to overpurchase
                "contracted_cost_M": contracted_cost / 1e6,
                "contracted_Mg": contracted_Mg_val,
                "feedstock_cost_M": feedstock_cost/ 1e6, # cost of feedstock only, without storage = contracted + extra - salvage
                "salvage_price_M": salvage_price / 1e6,
                "salvaged_biomass_Mg": salvaged_biomass_Mg, # to sell to an alternative market
                "DML_Mg": DML_Mg, # effective dry matter loss in supply chain
                "farm_count": farm_count_val,
                "transport_cost_M": transport_cost / 1e6,
                "msp": msp, "msp_adj": msp_adj,
                "feedstock_delivered_cost_M":    total_system / 1e6, # includes farm + extra + storage + transport
                # "total_system_demand_Mg": storage_demand_Mg + conv_gap_Mg,
                # **{f"month_{m}_frac": month_weights_list[i].get(m, 0.0)
                #    for m in [10, 11, 12, 1, 2, 3, 4, 5]},

            })

        if verbose and (i + 1) % 100 == 0:
            total_calls = (i + 1) * len(COMBOS)
            print(f"  {i+1}/{N}  |  transport cache hits: {cache_hits}/{total_calls} ({100*cache_hits/total_calls:.0f}%)")

    return pd.DataFrame(records)


# =============================================================================
# SUMMARY TABLE & VISUALIZATION
# =============================================================================

def _print_summary(df):
    print("\n" + "=" * 80)
    print("COUPLED MC SUMMARY  —  total system cost ($M)  [median, p10–p90]")
    print("=" * 80)
    hdr = f"  {'Combo':<30}  {'Median':>8}  {'p10':>8}  {'p90':>8}  {'SL%':>6}  {'MSP':>8}  {'MSP_adj':>9}  {'Ctrl$M':>8}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for combo, grp in df.groupby("combo", sort=False):
        vals = grp["feedstock_delivered_cost_M"].dropna().values.astype(float)
        p10, med, p90 = np.percentile(vals, [10, 50, 90])
        sl   = grp["service_level"].median()
        ctrl = grp["control_extra_cost_M"].median()
        print(
            f"  {combo:<30}  ${med:>7.2f}  ${p10:>7.2f}  ${p90:>7.2f}"
            f"  {sl:>5.1f}%  ${grp['msp'].median():>7.4f}"
            f"  ${grp['msp_adj'].median():>8.4f}  ${ctrl:>7.4f}"
        )


def plot_coupled_mc(df, save_path="coupled_mc.png"):
    apply_pub_style()
    combos    = df["combo"].unique().tolist()
    n_combos  = len(combos)
    BASE_COLORS  = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"]
    short_labels = [c.replace("tarped","tar").replace("square","sq").replace("round","rd") for c in combos]

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    plt.subplots_adjust(hspace=0.44, wspace=0.36)
    twins = {key: style_axes(ax, grid=True) for key, ax in np.ndenumerate(axes)}

    # (a) Total System Cost
    ax = axes[0, 0]
    data  = [df[df.combo == c]["feedstock_delivered_cost_M"].values for c in combos]
    parts = ax.violinplot(data, showmedians=True)
    for j, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(BASE_COLORS[j % len(BASE_COLORS)])
        pc.set_edgecolor("black"); pc.set_alpha(0.75)
    parts["cmedians"].set_color("black")
    ax.set_xticks(range(1, n_combos + 1))
    ax.set_xticklabels(short_labels, rotation=20, ha="right")
    ax.set_title("(a) Total System Cost"); ax.set_ylabel("Cost ($M)")
    sync_twins(ax, *twins[(0, 0)])

    # (b) Cost Components (median) — now includes control cost
    ax  = axes[0, 1]
    x   = np.arange(n_combos)
    s_m  = [df[df.combo==c]["storage_cost_M"].median()    for c in combos]
    c_m  = [df[df.combo==c]["contracted_cost_M"].median() for c in combos]
    e_m  = [df[df.combo==c]["extra_cost_M"].median()      for c in combos]
    tr_m = [df[df.combo==c]["transport_cost_M"].median()  for c in combos]
    ct_m = [df[df.combo==c]["control_extra_cost_M"].median() for c in combos]

    ax.bar(x, s_m,  color="#1f77b4", label="Storage",   edgecolor="white", lw=0.5)
    ax.bar(x, c_m,  bottom=s_m, color="#2ca02c", label="Contracted", edgecolor="white", lw=0.5)
    bot2 = [a+b for a,b in zip(s_m, c_m)]
    ax.bar(x, e_m,  bottom=bot2, color="#d62728", label="Extra",      edgecolor="white", lw=0.5)
    bot3 = [a+b for a,b in zip(bot2, e_m)]
    ax.bar(x, tr_m, bottom=bot3, color="#ff7f0e", label="Transport",  edgecolor="white", lw=0.5)
    bot4 = [a+b for a,b in zip(bot3, tr_m)]
    ax.bar(x, ct_m, bottom=bot4, color="#9467bd", label="Control",    edgecolor="white", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(short_labels, rotation=20, ha="right")
    ax.set_title("(b) Cost Components (median)"); ax.set_ylabel("Cost ($M)")
    ax.legend(frameon=True, fontsize=9)
    sync_twins(ax, *twins[(0, 1)])

    # (c) MSP adjusted
    ax = axes[0, 2]
    data_msp  = [df[df.combo == c]["msp_adj"].values for c in combos]
    parts_msp = ax.violinplot(data_msp, showmedians=True)
    for j, pc in enumerate(parts_msp["bodies"]):
        pc.set_facecolor(BASE_COLORS[j % len(BASE_COLORS)])
        pc.set_edgecolor("black"); pc.set_alpha(0.75)
    parts_msp["cmedians"].set_color("black")
    ax.set_xticks(range(1, n_combos + 1))
    ax.set_xticklabels(short_labels, rotation=20, ha="right")
    ax.set_title("(c) Min. Selling Price (adjusted)"); ax.set_ylabel("MSP ($/gal)")
    sync_twins(ax, *twins[(0, 2)])

    # (d–f) Scatter plots
    metrics  = [
        ("penalty",   "Unmet Demand Penalty"),
        ("land_density", "Suitable Land Density"),
        ("unmet_Mg",  "Unmet Demand (Mg wet)"),
    ]
    y_cols   = ["feedstock_delivered_cost_M", "feedstock_delivered_cost_M", "transport_cost_M"]
    titles   = ["(d) Total Cost vs Penalty", "(e) Total Cost vs Land Density", "(f) Transport Cost vs Unmet Demand"]
    y_labels = ["Total Cost ($M)", "Total Cost ($M)", "Transport Cost ($M)"]

    for idx, (x_col, x_label) in enumerate(metrics):
        ax = axes[1, idx]
        ax.set_title(titles[idx])
        for j, c in enumerate(combos):
            sub = df[df.combo == c]
            ax.scatter(sub[x_col], sub[y_cols[idx]], c=BASE_COLORS[j],
                       s=6, alpha=0.35, label=short_labels[j], edgecolors="none")
        ax.set_xlabel(x_label); ax.set_ylabel(y_labels[idx])
        ax.legend(frameon=True, fontsize=8, markerscale=2)
        sync_twins(ax, *twins[(1, idx)])

    fig.suptitle(
        f"Coupled Storage × Farm × Transport × Refinery MC (N={df['sample'].nunique()})",
        fontsize=12, y=1.02
    )
    save_fig(fig, save_path, dpi=300)
    plt.show()
    plt.close()


if __name__ == "__main__":
    print("Running optimized coupled MC (N=1000) …")
    df_full = run_coupled_mc(N=1000, seed=42, verbose=True)
    _print_summary(df_full)
    df_full.to_csv("coupled_mc_results_full.csv", index=False)
    plot_coupled_mc(df_full)

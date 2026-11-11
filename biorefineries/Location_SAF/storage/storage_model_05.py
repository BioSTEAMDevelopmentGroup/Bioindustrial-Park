import numpy as np
import pandas as pd
from dataclasses import dataclass, replace


# =============================================================================
# WEIGHTED PERCENTILE HELPER
# =============================================================================

def _weighted_percentile(values, weights, percentile):
    """
    Mass-weighted percentile of `values`, weighted by `weights`, without
    materializing a repeated array (see call site in Simulation.run() for
    why that matters: weights here are dispatch masses in kg, so a
    np.repeat-based approach builds an array with hundreds of millions of
    elements for a single year's dispatch history).

    Sorts once (O(n log n) in the number of dispatch records, not in the
    total weight), then walks the cumulative-weight curve to find the value
    at which the target fraction of total weight has been reached —
    equivalent to np.percentile(np.repeat(values, weights), percentile)
    but without ever building that repeated array.
    """
    values  = np.asarray(values,  dtype=float)
    weights = np.asarray(weights, dtype=float)

    order          = np.argsort(values)
    values_sorted  = values[order]
    weights_sorted = np.maximum(weights[order], 0.0)

    cum_weights = np.cumsum(weights_sorted)
    total       = cum_weights[-1]
    if total <= 0:
        return 0.0

    target = (percentile / 100.0) * total
    idx    = np.searchsorted(cum_weights, target)
    idx    = min(idx, len(values_sorted) - 1)
    return float(values_sorted[idx])


# =============================================================================
# EARLY HARVEST CUTOFF HELPER
# =============================================================================

# EARLY_HARVEST_CUTOFF_MONTH = 12
# EARLY_HARVEST_CUTOFF_DOM   = 15   # mid-December


# def _early_harvest_cutoff(start_date, harvest_days):
#     """
#     Return the last 0-based simulation day that counts as early harvest
#     (Oct through mid-December). Returns None if no early-harvest days exist
#     in the harvest window.
#     """
#     start = pd.Timestamp(start_date)
#     cutoff = None
#     for d in range(harvest_days):
#         dt = start + pd.Timedelta(days=d)
#         if dt.month == EARLY_HARVEST_CUTOFF_MONTH and dt.day >= EARLY_HARVEST_CUTOFF_DOM:
#             cutoff = d - 1   # last day before mid-Dec
#             break
#         if dt.month > EARLY_HARVEST_CUTOFF_MONTH:
#             cutoff = d - 1
#             break
#     # If the entire harvest window falls before mid-Dec, all days are early
#     return cutoff if cutoff is not None else harvest_days - 1

EARLY_HARVEST_MONTHS = {10, 11}
EARLY_HARVEST_CUTOFF_MONTH = 12
EARLY_HARVEST_CUTOFF_DOM = 15   # December 15


def _early_harvest_cutoff(start_date, harvest_days):
    """
    Return the last 0-based simulation day that counts as early harvest:
    October, November, and December 1–14.

    Returns None if the harvest window contains no early-harvest days.
    """
    start = pd.Timestamp(start_date)

    early_days = []

    for d in range(harvest_days):
        dt = start + pd.Timedelta(days=d)

        is_early = (
            dt.month in EARLY_HARVEST_MONTHS
            or (
                dt.month == EARLY_HARVEST_CUTOFF_MONTH
                and dt.day < EARLY_HARVEST_CUTOFF_DOM
            )
        )

        if is_early:
            early_days.append(d)

    if not early_days:
        return None

    return max(early_days)

# =============================================================================
# DATACLASSES
# =============================================================================

@dataclass
class StorageConfig:
    storage_type: str
    bale_shape: str
    fixed_cost_per_kg_capacity: float
    variable_cost_per_kg_throughput: float


@dataclass
class StorageScenario:
    name: str = "default"
    dml_max: float = 0.10
    ql_max: float = 0.05
    max_days: int = 365
    rule: str = "FIFO"
    alpha: float = 1.0
    storage_type: str = "tarped"
    bale_shape: str = "round"
    fixed_cost_per_kg_capacity: float = 0.0
    variable_cost_per_kg_throughput: float = 0.0
    # Annual physical-yield multiplier; alpha remains the procurement decision.
    yield_multiplier: float = 1.0


# =============================================================================
# PLANT
# =============================================================================

class Plant:
    def __init__(self, annual_demand_kg, operating_days=350, harvest_days=150):
        self.annual_demand = annual_demand_kg
        self.operating_days = operating_days
        self.harvest_days = harvest_days
        self.daily_demand = annual_demand_kg / operating_days


# =============================================================================
# STORAGE TECH
# =============================================================================

class StorageTech:
    def __init__(self, name, dml_max, ql_max, max_days):
        self.name = name
        self.dml_max = dml_max
        self.ql_max = ql_max
        self.max_days = max_days


# =============================================================================
# STORAGE
# =============================================================================

class Storage:
    def __init__(self, tech, storage_days=365 * 3, rule="FIFO", config=None):
        self.tech = tech
        self.storage_days = storage_days
        self.rule = rule
        self.config = config

        self.mass = np.zeros(storage_days)
        self.quality = np.ones(storage_days)

        self.early_harvest_mass = np.zeros(storage_days)

        self.cumulative_stored_gross = 0.0
        self.max_inventory_gross = 0.0

        # Pre-calculate static arrays
        self._ages = np.arange(self.storage_days)
        self._loss_factor = (
            1.0 - self.tech.dml_max * (1.0 - np.exp(-self._ages / self.tech.max_days))
        )
        self._effective_quality = (
            1.0 - self.tech.ql_max * (1.0 - np.exp(-self._ages / self.tech.max_days))
        )

    def add_biomass(self, amount, is_early_harvest=False):
        self.mass[0] += amount

        if is_early_harvest:
            self.early_harvest_mass[0] += amount

        self.quality[0] = 1.0
        self.cumulative_stored_gross += amount

        current_total = self.mass.sum()
        if current_total > self.max_inventory_gross:
            self.max_inventory_gross = current_total

    def age_storage(self):
        self.mass[1:] = self.mass[:-1]
        self.mass[0] = 0.0

        self.early_harvest_mass[1:] = self.early_harvest_mass[:-1]
        self.early_harvest_mass[0] = 0.0

        self.quality[1:] = self.quality[:-1]
        self.quality[0] = 1.0

    def remove_biomass(self, demand):
        # 1. Calculate available biomass today across all bins
        effective_mass = self.mass * self._loss_factor
        total_available = effective_mass.sum()

        if total_available <= 0 or demand <= 0:
            return 0.0, 0.0, 0.0, np.zeros(self.storage_days), 0.0

        take_total = min(total_available, demand)

        # 2. Vectorized cumulative sum allocation based on rule
        if self.rule == "FIFO":
            # FIFO: Take from oldest bins first (highest index down to 0)
            cum_available = np.cumsum(effective_mass[::-1])[::-1]
            cum_next = np.zeros_like(cum_available)
            cum_next[:-1] = cum_available[1:]
            take_effective = np.where(
                cum_available <= take_total,
                effective_mass,
                np.where(cum_next < take_total, take_total - cum_next, 0.0)
            )
        else:
            # LIFO: Take from newest bins first (index 0 upwards)
            cum_available = np.cumsum(effective_mass)
            cum_next = np.zeros_like(cum_available)
            cum_next[1:] = cum_available[:-1]
            take_effective = np.where(
                cum_available <= take_total,
                effective_mass,
                np.where(cum_next < take_total, take_total - cum_next, 0.0)
            )

        # 3. Update the true underlying mass array
        nonzero_mask = effective_mass > 0

        removal_fraction = np.where(
            nonzero_mask,
            take_effective / np.where(nonzero_mask, effective_mass, 1.0),
            0.0
        )

        # Update total biomass
        self.mass = np.where(
            nonzero_mask,
            self.mass * (1.0 - removal_fraction),
            0.0
        )
        self.mass[self.mass < 1e-9] = 0.0

        # Update early-harvest biomass using the same removal fraction
        self.early_harvest_mass *= (1.0 - removal_fraction)
        self.early_harvest_mass[self.early_harvest_mass < 1e-9] = 0.0

        # 4. Compute metrics
        avg_quality = np.sum(take_effective * self._effective_quality) / take_total
        avg_age     = np.sum(take_effective * self._ages) / take_total
        avg_dml = np.sum(take_effective * (1-self._loss_factor)) / take_total

        # Return take_effective so callers can do per-bin accounting
        return take_total, avg_quality, avg_age, take_effective, avg_dml

    def total_inventory(self):
        return self.mass.sum()

    def inventory_with_loss(self):
        return (self.mass * self._loss_factor).sum()


# =============================================================================
# SIMULATION
# =============================================================================

CONTROL_MAX_DAYS = 42   # first 6 weeks of storage subject to control cost


class Simulation:
    def __init__(self, plant, storage, alpha=1.0, early_harvest_cutoff_day=None,
                 yield_multiplier=1.0, fixed_cost_capacity_kg=None):
        self.plant = plant
        self.storage = storage
        self.alpha = alpha
        if not np.isfinite(yield_multiplier) or yield_multiplier < 0:
            raise ValueError("yield_multiplier must be finite and nonnegative")
        self.yield_multiplier = float(yield_multiplier)
        self.fixed_cost_capacity_kg = fixed_cost_capacity_kg
        # 0-based simulation day index of the last early-harvest day.
        # None means control cost tracking is disabled for this scenario.
        self.early_harvest_cutoff_day = early_harvest_cutoff_day

    def run(self):
        days = 365
        total_harvested = 0.0
        total_supplied  = 0.0
        unmet_demand    = 0.0

        quality_history      = []
        mean_age_history     = []
        old_inventory_history = []
        max_age_history      = []
        inventory_history    = []
        dispatch_records     = []
        dml_history = []

        # Accumulator for the control-cost calculation (barn/square/early only)
        early_controlled_kg_days = 0.0
        track_control = self.early_harvest_cutoff_day is not None

        storage        = self.storage
        daily_demand   = self.plant.daily_demand
        harvest_days   = self.plant.harvest_days
        operating_days = self.plant.operating_days
        harvest_amount = self.alpha * self.yield_multiplier * (self.plant.annual_demand / harvest_days)
        ages           = storage._ages
        loss_factor    = storage._loss_factor

        for day in range(days):
            storage.age_storage()

            if day < harvest_days:
                is_early_harvest = (
                    track_control
                    and day <= self.early_harvest_cutoff_day
                )

                storage.add_biomass(
                    harvest_amount,
                    is_early_harvest=is_early_harvest
                )

                total_harvested += harvest_amount

            # --- Control labor cost: early-harvest inventory within first 42 days ---
            if track_control:
                early_mask = ages <= CONTROL_MAX_DAYS

                early_controlled_kg_days += np.sum(
                    storage.early_harvest_mass[early_mask]
                )


            if day < operating_days:
                removed, avg_q, avg_age, take_effective, avg_dml = storage.remove_biomass(daily_demand)
                total_supplied += removed

                if removed < daily_demand:
                    unmet_demand += (daily_demand - removed)

                if removed > 0:
                    quality_history.append(avg_q)
                    dml_history.append(avg_dml)
                    dispatch_records.append((removed, avg_age, avg_dml))

            total_inventory = storage.mass.sum()
            if total_inventory > 0:
                mean_age = np.sum(ages * storage.mass) / total_inventory
                nonzero_indices = np.nonzero(storage.mass)[0]
                max_age = nonzero_indices[-1] if nonzero_indices.size > 0 else 0
            else:
                mean_age = 0.0
                max_age  = 0.0

            mean_age_history.append(mean_age)
            max_age_history.append(max_age)
            old_inventory_history.append(storage.mass[365:].sum())
            # Same value as the old mass_history_matrix[day,:]=storage.mass
            # followed by (mass_history_matrix * loss_factor).sum(axis=1) at
            # the end — computed inline per-day instead, so the full
            # (365 x storage_days) matrix never needs to be allocated.
            inventory_history.append(float(np.dot(storage.mass, loss_factor)))

        if len(dispatch_records) > 0:
            arr_records      = np.array(dispatch_records)
            mean_dispatch_age = np.sum(arr_records[:, 0] * arr_records[:, 1]) / np.sum(arr_records[:, 0])
            mean_dispatch_dml = (np.sum(arr_records[:,0] * arr_records[:,2])/ np.sum(arr_records[:,0]))
        else:
            mean_dispatch_age = 0.0
            mean_dispatch_dml = 0.0

        cfg          = storage.config
        fixed_capacity = (storage.max_inventory_gross if self.fixed_cost_capacity_kg is None
                          else self.fixed_cost_capacity_kg)
        fixed_cost   = cfg.fixed_cost_per_kg_capacity * fixed_capacity
        variable_cost = cfg.variable_cost_per_kg_throughput * storage.cumulative_stored_gross

        # new 7-Aug #
        inventory_array = np.array(inventory_history)

        average_inventory = inventory_array.mean()
        peak_inventory = inventory_array.max()
        days_empty = np.sum(inventory_array < 1e-6)

        if dispatch_records:
            dispatch_mass = np.array([x[0] for x in dispatch_records])
            dispatch_age = np.array([x[1] for x in dispatch_records])

            # Mass-weighted percentile without materializing a repeated array.
            # The previous np.repeat(dispatch_age, dispatch_mass.astype(int))
            # approach built an array of ~sum(dispatch_mass) elements — with
            # dispatch_mass in kg (millions/day), that's hundreds of millions
            # of elements per call, which is what was dominating runtime
            # (np.percentile's internal sort/partition over that array).
            # This computes the same weighted percentile in O(n log n) with
            # n = number of dispatch days (at most ~350), by sorting once and
            # walking the cumulative-mass curve to the target quantile.
            age_p90 = _weighted_percentile(dispatch_age, dispatch_mass, 90)
            age_p95 = _weighted_percentile(dispatch_age, dispatch_mass, 95)

        else:
            age_p90 = 0
            age_p95 = 0

        #### ---####

        # =========================================================================
        # END-OF-YEAR SALVAGE INVENTORY
        # =========================================================================
        # storage.mass = gross physical biomass remaining in each age bin.
        # Apply the age-dependent loss factor to obtain recoverable biomass.

        salvage_inventory_gross = storage.mass.sum()

        salvage_inventory_after_dml = np.sum(
            storage.mass * loss_factor
        )

        # Biomass already lost to DML and therefore NOT available for salvage
        salvage_dml_loss = (
            salvage_inventory_gross
            - salvage_inventory_after_dml
        )

        return {
            "alpha":                    self.alpha,
            "yield_multiplier":         self.yield_multiplier,
            "alpha_effective":          self.alpha * self.yield_multiplier,
            "contracted_biomass":       self.alpha * self.plant.annual_demand,
            "fixed_cost_capacity_kg":    fixed_capacity,
            "capacity_excess_kg":        max(storage.max_inventory_gross - fixed_capacity, 0.0),
            "harvested_biomass":        total_harvested,
            "supplied_biomass":         total_supplied,
            "unmet_demand":             unmet_demand,
            "inventory":                inventory_history,
            "avg_quality":              np.mean(quality_history) if quality_history else 0.0,
            "avg_dml":              np.mean(dml_history) if dml_history else 0.0,
            "mean_age_history":         mean_age_history,
            "old_inventory_history":    old_inventory_history,
            "max_age_history":          max_age_history,
            "max_storage_required":     storage.max_inventory_gross,
            "total_stored_mass_gross":  storage.cumulative_stored_gross,
            "fixed_storage_cost":       fixed_cost,
            "variable_storage_cost":    variable_cost,
            "total_storage_cost":       fixed_cost + variable_cost,
            "mean_dispatch_age":        mean_dispatch_age,
            "mean_dispatch_dml": mean_dispatch_dml,
            "dispatch_records":         dispatch_records,
            # kg·days for early-harvested bales within their first 42 days.
            # Zero for all non-barn/non-square scenarios (tracking disabled).
            "early_controlled_kg_days": early_controlled_kg_days,
            "average_inventory": average_inventory,
            "peak_inventory": peak_inventory,
            "days_empty": days_empty,
            "dispatch_age_p90": age_p90,
            "dispatch_age_p95": age_p95,
            "max_dispatch_age": max(dispatch_age) if dispatch_records else 0,
            "salvage_biomass_gross_Mg": salvage_inventory_gross / 1000.0,
            "salvage_biomass_Mg": salvage_inventory_after_dml / 1000.0,
            "salvage_dml_loss_Mg": salvage_dml_loss / 1000.0,
        }


# =============================================================================
# SIMULATION RUNNER
# =============================================================================

class SimulationRunner:
    def __init__(self, plant, start_date=None, harvest_days=None):
        self.plant        = plant
        self.start_date   = start_date    # e.g. "2025-11-01" or pd.Timestamp
        self.harvest_days = harvest_days  # e.g. 150

    def run(self, scenario: StorageScenario):
        if not np.isfinite(scenario.yield_multiplier) or scenario.yield_multiplier < 0:
            raise ValueError("yield_multiplier must be finite and nonnegative")
        # Plan capacity at normal annual yield, with the SAME alpha, harvest
        # window, dispatch rule, and decay draw. Do not resize capital after
        # observing the annual yield shock. This is a cost basis, not a hard
        # inventory cap; capacity_excess_kg reports any physical exceedance.
        fixed_capacity = None
        if scenario.yield_multiplier != 1.0:
            planned = self.run(replace(scenario, yield_multiplier=1.0))
            fixed_capacity = planned["max_storage_required"]
        tech = StorageTech(
            name=scenario.name,
            dml_max=scenario.dml_max,
            ql_max=scenario.ql_max,
            max_days=scenario.max_days,
        )
        config = StorageConfig(
            storage_type=scenario.storage_type,
            bale_shape=scenario.bale_shape,
            fixed_cost_per_kg_capacity=scenario.fixed_cost_per_kg_capacity,
            variable_cost_per_kg_throughput=scenario.variable_cost_per_kg_throughput,
        )
        storage = Storage(tech=tech, rule=scenario.rule, config=config)

        # Compute early-harvest cutoff day only for barn/square combos.
        # For everything else, cutoff_day stays None → tracking is disabled,
        # early_controlled_kg_days will be 0 in the output.
        cutoff_day = None
        if (scenario.storage_type == "barn"
                and scenario.bale_shape == "square"
                and self.start_date is not None
                and self.harvest_days is not None):
            cutoff_day = _early_harvest_cutoff(self.start_date, self.harvest_days)

        sim = Simulation(
            plant=self.plant,
            storage=storage,
            alpha=scenario.alpha,
            early_harvest_cutoff_day=cutoff_day,
            yield_multiplier=scenario.yield_multiplier,
            fixed_cost_capacity_kg=fixed_capacity,
        )



        return sim.run()

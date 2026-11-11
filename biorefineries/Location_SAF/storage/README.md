# Coupled biomass storage, farm allocation, transport, and refinery model

This workflow couples biomass harvest
and storage simulation with spatial farm allocation, truck transportation, and
location-dependent refinery minimum selling price (MSP) estimates.

This is a storage study. It is separate from the cellulosic ethanol-to-SAF siting
workflow in the parent `Location_SAF` directory. The six model files alone are
not sufficient to run: the plotting helper and input datasets listed below are
also required.

## Files

| File | Purpose |
| --- | --- |
| `coupled_model_mc_10.py` | Coupled Monte Carlo simulation and diagnostic plots |
| `storage_model_05.py` | Harvest, inventory, degradation, dispatch, and storage costs |
| `farm_allocation_model_07.py` | Spatial feedstock allocation and farm costs |
| `biomass_transport_model_02.py` | Truck transportation distances and costs |
| `refinery_model.py` | Refinery MSP surrogate calculation |
| `bc_functions.py` | Location-dependent interpolation of refinery coefficients |
| `pub_style.py` | Plot styling used by the coupled model |
| `run_example.py` | Small configurable command-line example |
| `requirements.txt` | Dependency versions from the author's environment |
| `DATA.md` | Input inventory and provenance |
| `prepare_inputs.py`, `INPUTS.json` | Prepare and verify existing repository inputs |

## Using the modules separately

The farm allocation, storage, transport, and refinery modules can be used
individually by supplying their required inputs, or linked through the coupled
Monte Carlo model. Run the input preparation step below before using modules that
read GIS or refinery coefficient files. For standalone use, run from this folder
or provide explicit input paths.

### Farm allocation — `farm_allocation_model_07.py`

`FeedstockAllocationModel` allocates farms to a refinery based on their distance
and available biomass. Inputs include farm coordinates, crop yields and farm-gate
prices, refinery location and production capacity, feedstock type (miscanthus or
switchgrass), moisture content, and biomass-to-product conversion efficiency.
The prepared feedstock tables and GIS boundaries provide the spatial inputs;
`base_folder` can point to another compatible input directory.

Adjust `suitable_land_density` to change the fraction of land available for
feedstock production, `yield_bonus` to represent yield changes, and
`penalty_for_unmet_demand` to change the price premium for additional biomass.
`extra_biomass_with_yield_bonus` determines whether yield changes affect initial
farm selection or the actual harvest from farms contracted at baseline yield.
Call `calculate_location_parameters()` to populate farm allocation and cost
results, including selected farms, biomass quantities, farm counts, and distances.

### Storage — `storage_model_05.py`

The storage module simulates daily biomass inventory and dispatch over one year.
`Plant` specifies annual biomass demand (kg), operating days, and harvest duration.
`StorageTech`, `StorageConfig`, and `Storage` describe degradation, costs, and
inventory dispatch; `Simulation` runs the daily balance. It tracks biomass supplied,
unmet demand, remaining inventory, dry-matter loss, quality, dispatch age, storage
capacity requirements, and costs.

Scenarios can represent **tarped, barn, or uncovered storage**, **square or round
bales**, and **FIFO** (oldest inventory first) or **LIFO** (newest inventory first)
dispatch. Supply the corresponding parameters for each scenario: storage-type and
bale-shape labels alone do not select cost or degradation values in the standalone
module. Adjust fixed cost per kg of capacity, variable cost per kg of throughput,
maximum dry-matter and quality loss (`dml_max`, `ql_max`), and the degradation time
scale (`max_days`). Here, `max_days` is used in the decay curve, not as a mandatory
inventory disposal age.

The procurement multiplier `alpha` controls planned buffer purchases: for example,
`alpha=1.10` plans 10% more biomass than baseline annual demand. The separate
`yield_multiplier` changes realized harvest. The standalone storage model reports
shortfalls; the coupled model adds spot purchases and their costs. Uncertainty
ranges for decay and buffer purchases are defined in the coupled module rather
than sampled automatically by the standalone storage classes.

### Biomass transport — `biomass_transport_model_02.py`

`BiomassTransport` estimates truck trips, farm-to-refinery distances, and transport
costs. Provide refinery coordinates, allocated farm tables (including `lat`,
`long`, `distance_to_refinery_km`, and `biomass_supplied_Mg`), the biomass quantity
to transport, and state-level costs per kilometre. Bale mass and the number of
bales per truck can be adjusted for round or square bales.

The default distance method applies state-specific tortuosity factors to
straight-line distances using the boundary shapefile and tortuosity spreadsheet.
A constant tortuosity factor and straight-line distance are also available;
road-network and neural-network methods require the extra dependencies described
below. Call `compute()` to calculate results, then use `summary()`,
`get_farm_results()`, or `get_total_transport_costs()` to inspect them.

### Refinery — `refinery_model.py`

`RefineryModel` estimates minimum selling price from refinery location, capacity,
and annual feedstock, storage, and transport costs. Moisture content and conversion
efficiency determine biomass demand. Its `compute()` method returns MSP,
feedstock price, biomass demand, and separate refinery, feedstock, storage, and
transport contributions.

This module uses a location-dependent surrogate rather than running a full
BioSTEAM process simulation. It can be supplied with costs from other models or
user-defined scenarios, but requires `bc_functions.py` and the prepared coefficient
arrays. The example workflow uses ethanol; other product coefficient functions
are retained but are not validated by the example.

### Refinery coefficients — `bc_functions.py`

This helper interpolates the location-dependent `b` and `c` coefficients used by
the refinery model. Supply location coordinates in the same order as the refinery
model, `(latitude, longitude)`. `b_c_vectorized()` provides ethanol coefficients;
additional functions provide acrylic acid, succinic acid, potassium sorbate, and
lactic acid coefficients. It can be used separately to inspect coefficients at
candidate sites, provided the input arrays have been prepared.

### Coupled Monte Carlo model — `coupled_model_mc_10.py`

`run_coupled_mc()` links harvest and storage behavior with farm allocation,
transport, and refinery pricing. It samples uncertain inputs and evaluates the
storage combinations listed in `COMBOS`. Adjust the sample count and seed, refinery
location and size, harvest window, `DML_RANGES`, `ALPHA_RANGES`, storage unit costs,
transport assumptions, and harvest/storage conversion effects to define scenarios.
It accounts for additional spot purchases when biomass supply is insufficient
and returns a table of sampled inputs and costs for each sample and combination.

`pub_style.py` supplies plotting styles; it does not change model calculations.
`run_example.py` provides a command-line entry point to the coupled workflow.

## Installation

The author's environment uses Python 3.9.21. Create a separate environment instead
of installing these versions into an existing BioSTEAM environment. This storage
workflow uses saved refinery surrogate arrays; it does not import BioSTEAM.

From this directory:

```bash
python3.9 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python prepare_inputs.py
```

On Windows, activate with `.venv\Scripts\activate` instead. In VS Code, open this
folder and choose the environment with **Python: Select Interpreter**. Run the
commands below in its integrated terminal.

The setup command verifies and copies existing repository data into local input
folders. These generated copies are ignored by Git. Keep the complete repository
checkout (including the source paths in `INPUTS.json`); downloading only this
`storage` folder is not sufficient. Run `python prepare_inputs.py --check` to
verify the prepared inputs.

## Run a small example first

```bash
python run_example.py --samples 2 --seed 42 --output outputs/example.csv
```

This produces one row per Monte Carlo draw and storage combination: 8 rows with
the four default combinations. It is a functionality check, not a converged
uncertainty analysis. Importing `bc_functions.py` constructs five interpolators
from million-row arrays, so startup can take several minutes and substantial RAM.

For the default 1,000-draw run:

```bash
python run_example.py --samples 1000 --seed 42 --output outputs/coupled_mc.csv --plot
```

For a noninteractive machine, set `MPLBACKEND=Agg` before running. The output CSV
path is relative to the working directory; the required input paths in the
coupled entry point are relative to this code directory.

You can also use the original entry point, `python coupled_model_mc_10.py`, which
runs 1,000 draws, writes `coupled_mc_results_full.csv`, and creates a plot.

## Default configuration

Study constants and scenario definitions are near the top of
`coupled_model_mc_10.py`. The shipped configuration uses:

- An 80 million gallon/year ethanol refinery at (35.4676, -97.5164)
  (latitude, longitude), Oklahoma.
- Miscanthus, 20% moisture content, and baseline conversion of 113.97 gal/dry Mg.
- Harvest start date 2025-11-01 and a 150-day harvest window.
- Tarped round and square bales, with FIFO and LIFO dispatch (four combinations).
- State-specific tortuosity for transportation distances.
- Harvest timing and storage conversion effects enabled.
- `run_coupled_mc` default maximum storage duration of 200 days.

Other storage technologies have parameter definitions but are not included in
the default `COMBOS`. Reproducing a specific manuscript analysis requires its
scenario settings, sample count, and seed; the default run does not claim to
reproduce every manuscript table or figure.

## Outputs and units

Each row identifies `sample`, `scenario`, `combo`, and sampled parameters.
Key outputs include `msp` and `msp_adj` (USD/gal), `feedstock_delivered_cost_M`,
`storage_cost_M`, and `transport_cost_M` (million USD), biomass quantities ending
in `_Mg` (metric tonnes; check the column's wet/dry basis in the source),
`service_level` (percent), and `mean_dispatch_age` (days). Full field definitions and
accounting are in the records dictionary returned by `run_coupled_mc`.
The random seed controls sampling; use the same dependency versions for comparisons.

## Scope and limitations

The default tortuosity workflow uses local files and does not need network access.
Optional `road_network` and `neural_network` distance methods in the transport
class are outside this example. They need additional packages (OSMnx/NetworkX or
PyTorch/joblib) and, for the neural network, model/scaler files not bundled here.
The other product interpolators are retained from the supplied code; this release
checks the ethanol workflow only.

The publication copy updates input path resolution in `bc_functions.py` and
`coupled_model_mc_10.py`, and corrects stale cost-column references in the summary
and plotting functions to `feedstock_delivered_cost_M`. Model equations and
default scientific parameters are preserved. The original local working files are not modified.

## Citation

The manuscript is in preparation.

# Coupled biomass storage, farm allocation, transport, and refinery model

Research code for the biomass storage project by Maria Florencia Bianco at the
University of Illinois Urbana-Champaign. This workflow couples biomass harvest
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

## Citation and licensing

Cite the associated storage study when its final bibliographic details are
available, and record the Git commit used for an analysis. Do not cite the SAF
siting paper as though it describes this storage model. For an archival code DOI,
a versioned repository release can later be deposited in Zenodo or another archive.
See the repository's top-level `LICENSE.txt` for its software license. Input data
provenance is documented separately in `DATA.md`; third-party terms remain applicable.

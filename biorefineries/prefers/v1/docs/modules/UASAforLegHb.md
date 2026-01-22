# Uncertainty & Sensitivity Analysis (`LegHb`)

**Directory:** `v1/LegHb/analyses/`  
**Data Script:** `gen_data.py`  
**Figures Script:** `gen_figures.py`  
**Orchestrator:** `UA_SA.py`

This module provides comprehensive tools for quantifying risk and parameter sensitivity in the LegHemoglobin biorefinery through Monte Carlo simulations. The pipeline is now split into a cold-data generation step and a separate plotting step to avoid redundant computation.

## 1. Analysis Strategy

The analysis supports **two main scenarios** with flexible configuration:

### Scenario 1: Fixed Scale Monte Carlo (`exclude_production_scale=True`)
- Keeps production scale constant (default: 275 kg/hr)
- Isolates biological and process uncertainties (yield, titer, productivity)
- **Outputs:** Spearman correlations, KDE plots, contour plots

### Scenario 2: Variable Scale Monte Carlo (`exclude_production_scale=False`)
- Varies production scale alongside fermentation parameters
- Studies economies of scale effects
- Useful for capacity planning decisions

## 2. Architecture

### Split Pipeline (Cold Data + Figures)

1. **Cold Data Generation** (`gen_data.py`)
    - Runs baseline evaluation, Monte Carlo (UA), and single-point sensitivity (SA)
    - Saves all raw results to `analyses/results_{config}_{timestamp}/data/`
    - Outputs robust, portable files (CSV, Excel, Pickle)

2. **Figure Generation** (`gen_figures.py`)
    - Loads the precomputed datasets from `/data`
    - Produces figures in `/figure` without re-running simulations

3. **Orchestrator** (`UA_SA.py`)
    - Optional convenience wrapper for running data and figures in sequence

### Utility Modules Used

The plotting and reporting helpers are centralized in `v1/utils/` and are used by `gen_data.py` and `gen_figures.py`:

- `utils/plot_utils.py`: common chart builders and styling
- `utils/sankey_utils.py`: Sankey preparation and export
- `utils/style.py`: global matplotlib style settings
- `utils/plots.py`: composite plotting helpers

### Core Functions

#### `evaluate_single_sample(sample_index_and_data, baseline_production_kg_hr, exclude_production_scale=False, config='config1')`
Worker function for parallel Monte Carlo evaluation. Evaluates a single sample point:
- Creates fresh model instance per sample
- Sets parameter values from sample array
- Simulates system and extracts metrics
- Returns results dict or None on failure

#### `run_monte_carlo(model, N_target, baseline_production_kg_hr, exclude_production_scale=False, batch_size=100, scenario_name="", config='config1')`
Robust parallel Monte Carlo driver:
- Generates Latin Hypercube samples
- Distributes work across CPU cores
- Collects and aggregates results
- Returns DataFrame with all valid samples

#### `generate_2d_contour_plots(results_table, param_indices, metric_indices, timestamp)`
Creates scatter plots with color-coded z-values showing joint impact of fermentation parameters on economic/environmental metrics.

#### `calculate_percentiles_by_bin(x_values, y_values, bins)`
Statistical utility for binned percentile analysis.

## 3. Model Parameters

The `_models.py` module defines the `create_model()` function which builds a `bst.Model` with uncertain parameters:

### Fermentation Performance Parameters

| Parameter | Distribution | Baseline | Units |
|-----------|--------------|----------|-------|
| `Titer` | Triangular/Uniform | 7.27 | g/L |
| `Productivity` | Triangular/Uniform | 0.101 | g/L/hr |
| `Yield` | Triangular/Uniform | 0.0224 | g product/g glucose |

### Operating Cost Parameters

| Parameter | Distribution | Baseline | Units |
|-----------|--------------|----------|-------|
| `Glucose Price` | Uniform | varies | $/kg |
| `Electricity Price` | Uniform | varies | $/kWh |

### Optional Scale Parameter

| Parameter | Distribution | Baseline | Units |
|-----------|--------------|----------|-------|
| `Production Rate` | Triangular | 275 | kg/hr |

## 4. Model Metrics

The analysis tracks Key Performance Indicators (KPIs) defined in `_models.py`:

### Economic Metrics

| Metric | Function | Units |
|--------|----------|-------|
| `MSP` | `get_MSP()` | $/kg LegHb |
| `TCI` | `get_TCI()` | $ Million |
| `AOC` | `get_AOC()` | $ Million/yr |
| `Specific CapEx` | `get_specific_capex()` | $/kg annual capacity |

### Environmental Metrics

| Metric | Function | Units |
|--------|----------|-------|
| `GWP` | `get_GWP()` | kg CO2-eq/kg LegHb |

**GWP Calculation Method:** Uses BioSTEAM's LCA displacement allocation table, extracting total GWP (sum of inputs - outputs + process impacts) per kg of product.

### Technical Metrics

| Metric | Function | Units |
|--------|----------|-------|
| `LegHb Content` | `get_LegHb_content()` | mass % |
| `Protein Purity` | `get_protein_purity()` | % of total protein |
| `Actual Production` | `get_actual_production()` | kg/hr |
| `Annual Production` | `get_annual_production()` | metric tons/yr |

## 5. Outputs

Running `gen_data.py` generates results in a timestamped subfolder `results_{config}_{YYYYMMDD_HHMM}/`:

### Excel Datasets

| File | Contents |
|------|----------|
| `baseline_metrics.xlsx` | Baseline KPI values |
| `monte_carlo_results.pkl` | Full MC results with MultiIndex (for plotting) |
| `monte_carlo_results.xlsx/csv` | Flat MC results for portability |
| `tornado_sensitivity.xlsx/csv` | Single-point sensitivity data |
| `analysis_metadata.json` | Run configuration + metric/parameter indices |

### Visualizations

Generated by `gen_figures.py` using the cold data in `/data`:

| Type | Description |
|------|-------------|
| **Spearman Heatmap** | Parameter ranking via Spearman correlation |
| **Joint Marginal Plot** | MSP vs GWP trade-off |
| **Tornado Plot** | Single-point MSP sensitivity |
| **Distribution Plots** | Uncertainty distributions for MSP and GWP |

## 6. Configuration Options

The analysis supports both `config1` and `config2` process configurations:

### Example Commands

```bash
# 1) Generate cold data only
python v1/LegHb/analyses/gen_data.py --config config1 --samples 500

# 2) Generate figures from existing results
python v1/LegHb/analyses/gen_figures.py --results-dir v1/LegHb/analyses/results_config1_YYYYMMDD_HHMM

# 3) One-shot (data + figures)
python v1/LegHb/analyses/UA_SA.py --config config1 --samples 500
```

## 7. Verification

The `_models.py` module includes `verify_model_integration()` function that tests:
- System creation and simulation
- Design specification mode (production rate setting)
- TEA integration and MSP calculation
- Metric extraction for all defined metrics

Run directly:
```bash
python -m biorefineries.prefers.v1.LegHb._models
```

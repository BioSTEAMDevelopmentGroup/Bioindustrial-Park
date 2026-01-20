# Uncertainty & Sensitivity Analysis (`LegHb`)

**Directory:** `v1/LegHb/analyses/`  
**Main Script:** `uncertainty_and_sensitivity.py`

This module provides comprehensive tools for quantifying risk and parameter sensitivity in the LegHemoglobin biorefinery through parallel Monte Carlo simulations.

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

Running `uncertainty_and_sensitivity.py` generates results in a timestamped subfolder `results_YYYY.MM.DD-HH.MM/`:

### Excel Datasets

| File Pattern | Contents |
|--------------|----------|
| `LegHb_MC_no_scale_*.xlsx` | Results with fixed scale |
| `LegHb_MC_with_scale_*.xlsx` | Results with variable scale |

### Visualizations

| Type | Description |
|------|-------------|
| **Scatter Plots (2D)** | Color-coded scatter plots (e.g., MSP vs Titer/Yield) showing joint parameter impacts |
| **KDE Plots** | 1D and 2D Kernel Density Estimates showing probability distributions of MSP and GWP |
| **Spearman Plots** | Tornado plots ranking parameters by Spearman correlation coefficient on MSP and GWP |
| **Contour Plots** | 2D contour visualizations of response surfaces |

## 6. Configuration Options

The analysis supports both `config1` and `config2` process configurations:

```python
# Example: Run fixed-scale analysis with config1
model = create_model(baseline_production_kg_hr=275, config='config1')
results, stats = run_monte_carlo(
    model, 
    N_target=500,
    baseline_production_kg_hr=275,
    exclude_production_scale=True,
    config='config1'
)
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

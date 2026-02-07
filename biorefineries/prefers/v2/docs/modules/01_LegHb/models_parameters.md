# LegHb Uncertainty Model: Parameters & Metrics (v2)

**Source:** `v2/LegHb/_models.py`  
**Model Factory:** `create_model(baseline_production_kg_hr=150, config='config1', verbose=True)`  
**Organism:** *Pichia pastoris*

## Overview

The LegHb uncertainty model quantifies risk and sensitivity by sampling **20 parameters** across 5 categories (Fermentation, DSP, Facilities, Economics, GWP) and evaluating **10 metrics** (4 economic/environmental + 6 composition).

### Changes from v1

Distribution shapes were upgraded in v2 to better reflect the statistical nature of each parameter category:

| Category | v1 Shape | v2 Shape | Rationale |
|----------|----------|----------|-----------|
| Titer | Triangle | Trunc(LogNormal) | Biological titers are right-skewed; PERT intended but `scipy.special.btdtri` unavailable in current env |
| Tau | Triangle | Triangle | Residence time uncertainty is symmetric |
| Yields | Triangle | TruncNormal | Measurement errors around established values are approximately normal |
| DSP efficiencies | Uniform | TruncNormal | Process efficiencies cluster around design targets (Beta intended for cell disruption but unavailable) |
| BT efficiency | Uniform | Uniform | Facility parameter, not DSP — unchanged |
| Raw material prices | Triangle | Trunc(LogNormal) | Commodity prices exhibit multiplicative (right-skewed) uncertainty |
| Operating days / IRR / Tax | Triangle | Triangle | Non-commodity parameters — unchanged |
| GWP factors | Triangle/Uniform | Trunc(LogNormal) | Emission factors follow multiplicative uncertainty similar to prices |

**Critical constraint:** All distributions preserve the same lower/upper bounds as v1, ensuring Baseline, Lower-Bound, and Upper-Bound evaluations are identical.

### Sigma Constants

```python
_ln_sigma_titer  = (np.log(1.5) - np.log(0.5)) / 4  # ≈ 0.2747
_ln_sigma_10pct  = (np.log(1.1) - np.log(0.9)) / 4  # ≈ 0.0502
```

- **`_ln_sigma_titer`**: For titer range [0.5×, 1.5×] baseline — ~95% of un-truncated LogNormal falls within bounds
- **`_ln_sigma_10pct`**: For ±10% range — ~95% of un-truncated LogNormal falls within bounds
- **TruncNormal sigma**: `half_range / 2` so ~95% of un-truncated Normal falls within bounds

---

## Model Parameters

### 1. Fermentation Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 1.1 | Fermentation titer | ~5.0 | g/L | **Trunc(LogNormal)** | 0.5× baseline | 1.5× baseline | `mu=ln(baseline)`, `σ=0.2747` |
| 1.2 | Fermentation tau | ~72 | hr | Triangle | 0.8× baseline | 1.2× baseline | Unchanged from v1 |
| 1.3 | Product yield | ~3.71 | % | **TruncNormal** | 0.9× baseline | 1.1× baseline | `mu=baseline`, `σ=baseline×0.05/2` |
| 1.4 | Biomass yield | ~53 | % | **TruncNormal** | 0.9× baseline | 1.1× baseline | `mu=baseline`, `σ=baseline×0.05/2` |

### 2. DSP Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 2.1 | Centrifuge split | 1.0 | multiplier | **TruncNormal** | 0.95 | 1.02 | `mu=1.0`, `σ=0.0175` |
| 2.2 | Cell disruption eff | 0.87 | — | **TruncNormal** | 0.82 | 0.92 | `mu=0.87`, `σ=0.025` |
| 2.3 | Filtration capture | 0.85 | — | **TruncNormal** | 0.75 | 0.95 | `mu=0.85`, `σ=0.05` |
| 2.4 | DF retention | 0.95 | — | **TruncNormal** | 0.90 | 0.99 | `mu=0.95`, `σ=0.0225` |

### 3. Facilities Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 3.1 | BT efficiency | ~0.85 | — | Uniform | 0.80 | 0.90 | Unchanged from v1 |

### 4. Economic Parameters (Isolated)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 4.1 | Electricity price | ~0.030 | $/kWh | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 4.2 | Glucose price | ~0.420 | $/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 4.3 | Ammonia price | ~0.115 | $/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 4.4 | Buffer/Seed cost | 1.0 | multiplier | **Trunc(LogNormal)** | 0.9 | 1.1 | `mu=0`, `σ=0.0502` |
| 4.5 | Operating days | 333 | days/yr | Triangle | 300 | 350 | Unchanged from v1 |
| 4.6 | IRR | 0.18 | fraction | Triangle | 0.10 | 0.25 | Unchanged from v1 |
| 4.7 | Income tax | 0.17 | fraction | Triangle | 0.10 | 0.25 | Unchanged from v1 |

### 5. GWP Parameters (Isolated)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 5.1 | Electricity GWP | ~0.550 | kg CO₂/kWh | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.2 | Glucose GWP | ~1.610 | kg CO₂/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.3 | Ammonia GWP | ~0.710 | kg CO₂/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.4 | Buffer/Seed GWP | 1.0 | multiplier | **Trunc(LogNormal)** | 0.9 | 1.1 | `mu=0`, `σ=0.0502` |

---

## Model Specification

The model specification function is called before each evaluation:
1. Calls `run_titer_convergence()` which iteratively:
   - Adjusts glucose for target titer (`adjust_glucose_for_titer`)
   - Optimizes NH₃ loading (`optimize_NH3_loading`)
   - Sets production rate (`set_production_rate`)

---

## Metrics

### Economic Metrics

| Metric | Units | Element | Description |
|--------|-------|---------|-------------|
| MSP | $/kg | PreFerS | Minimum Selling Price via `tea.solve_price(LegHb_3)` |
| TCI | 10⁶ $ | PreFerS | Total Capital Investment / 10⁶ |
| AOC | 10⁶ $/yr | PreFerS | Annual Operating Cost / 10⁶ |
| GWP | kg CO₂-eq/kg | PreFerS | LCA displacement allocation (total per kg product) |

### Composition Metrics

| Metric | Units | Element | Description |
|--------|-------|---------|-------------|
| Fat Content | wt% | Composition | OleicAcid / total mass × 100 |
| Carbohydrates | wt% | Composition | (Glucan + Glucose + Chitin) / total mass × 100 |
| Product Content | wt% | Composition | Leghemoglobin / total mass × 100 |
| Total Solids | wt% | Composition | (total − H₂O) / total × 100 |
| Protein Purity | % | Composition | LegHb / (LegHb + Globin + Mannoprotein) × 100 |
| Heme Equivalent | wt% | Composition | Based on LegHb molar content: `mol/763 × MW_Heme_b` |

---

## Summary Statistics

| Property | Value |
|----------|-------|
| Total parameters | 20 (4 ferm + 4 DSP + 1 facilities + 7 econ + 4 GWP) |
| Total metrics | 10 (4 economic/environmental + 6 composition) |
| Coupled parameters | 9 (require system re-simulation) |
| Isolated parameters | 11 (only affect TEA/LCA calculations) |
| Trunc(LogNormal) params | 9 (titer + 4 prices + 4 GWP) |
| TruncNormal params | 6 (2 yields + 4 DSP efficiencies) |
| Triangle params | 4 (tau + op days + IRR + tax) |
| Uniform params | 1 (BT efficiency) |

### Verified Baselines (v2)

| Scenario | MSP ($/kg) | GWP (kg CO₂-eq/kg) |
|----------|-----------|---------------------|
| Baseline | 11.81 | 17.21 |
| Lower Bound | 14.10 | — |
| Upper Bound | 10.07 | — |

## Last Updated
2026-02-07

# HemDx Uncertainty Model: Parameters & Metrics (v2)

**Source:** `v2/HemDx/_models.py`  
**Model Factory:** `create_model(baseline_production_kg_hr=150, config='config1', verbose=True)`  
**Organism:** *Corynebacterium glutamicum*  
**Product:** N-HemoDextrin (Nicotinamide-Stabilized Hemodextrin)

## Overview

The HemDx uncertainty model supports **3 configurations** (config1/2/3) with up to **23 parameters** across 5 categories and **10 metrics**. Parameter count varies by config as some DSP units are absent in config2 (no S404) and config3 (no S402, C403, S405).

### Changes from v1

Distribution shapes were upgraded in v2 to better reflect the statistical nature of each parameter category:

| Category | v1 Shape | v2 Shape | Rationale |
|----------|----------|----------|-----------|
| Titer | Triangle | Trunc(LogNormal) | Biological titers are right-skewed; PERT intended but `scipy.special.btdtri` unavailable in current env |
| Tau | Triangle | Triangle | Residence time uncertainty is symmetric |
| Yields | Triangle | TruncNormal | Measurement errors around established values are approximately normal |
| Secretion fraction | Uniform | Uniform | Config-dependent range too wide for bell-curve assumption — unchanged |
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
| 1.1 | Fermentation titer | `R302.titer` | g/L | **Trunc(LogNormal)** | 0.5× baseline | 1.5× baseline | `mu=ln(baseline)`, `σ=0.2747` |
| 1.2 | Fermentation tau | `R302.tau` | hr | Triangle | 0.8× baseline | 1.2× baseline | Unchanged from v1 |
| 1.3 | Product yield | `R302.target_yield×100` | % | **TruncNormal** | 0.9× baseline | 1.1× baseline | `mu=baseline`, `σ=baseline×0.05/2` |
| 1.4 | Biomass yield | CG2 growth yield×100 | % | **TruncNormal** | 0.9× baseline | 1.1× baseline | `mu=baseline`, `σ=baseline×0.05/2` |
| 1.5 | Secretion fraction | Config-dependent | fraction | Uniform | Config-dependent | Config-dependent | Unchanged from v1 |

#### Secretion Fraction Distribution by Config

| Config | Baseline SF | Distribution | Lower | Upper | Rationale |
|--------|------------|-------------|-------|-------|-----------|
| config1 (Base) | 0.45 | Uniform | 0.2 | 0.8 | Mixed secretion |
| config2 (Intracellular) | 0.10 | Uniform | 0.0 | 0.3 | Low secretion |
| config3 (Extracellular) | 0.90 | Uniform | 0.7 | 1.0 | High secretion |

### 2. DSP Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Configs | Notes |
|---|-----------|----------|-------|-------------|-------|-------|---------|-------|
| 2.1 | Centrifuge split | 1.0 | multiplier | **TruncNormal** | 0.95 | 1.02 | All | `mu=1.0`, `σ=0.0175` |
| 2.2 | Cell disruption eff | 0.87 | — | **TruncNormal** | 0.82 | 0.92 | 1, 2 | `mu=0.87`, `σ=0.025` |
| 2.3 | Filtration capture | 0.95 | — | **TruncNormal** | 0.90 | 0.99 | Varies | `mu=0.95`, `σ=0.0225` |
| 2.4 | DF retention | 0.95 | — | **TruncNormal** | 0.90 | 0.99 | All | `mu=0.95`, `σ=0.0225` |
| 2.5 | Resin yield | ~0.97 | — | **TruncNormal** | 0.95 | 0.99 | All | `mu=baseline`, `σ=0.01` |
| 2.6 | CSTR conversion | 0.95 | — | **TruncNormal** | 0.87 | 0.97 | All | `mu=0.95`, `σ=0.025` |

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
| 5.1 | Electricity GWP | ~0.550 | kg CO₂-eq/kWh | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.2 | Glucose GWP | ~1.610 | kg CO₂-eq/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.3 | Ammonia GWP | ~0.710 | kg CO₂-eq/kg | **Trunc(LogNormal)** | 0.9× base | 1.1× base | `mu=ln(base)`, `σ=0.0502` |
| 5.4 | Buffer/Seed GWP | 1.0 | multiplier | **Trunc(LogNormal)** | 0.9 | 1.1 | `mu=0`, `σ=0.0502` |

---

## Model Specification

The model specification function (called before each evaluation):
1. `adjust_glucose_for_titer()` — scales yield based on `R302.titer`
2. `optimize_NH3_loading()` — adjusts ammonia for nitrogen balance
3. `set_production_rate()` — scales system to target production
4. All wrapped in `run_titer_convergence()` for iterative convergence

---

## Metrics

### Economic Metrics

| Metric | Units | Element | Description |
|--------|-------|---------|-------------|
| MSP | $/kg | PreFerS | Minimum Selling Price via `tea.solve_price(NHemDx_Product)` |
| TCI | 10⁶ $ | PreFerS | Total Capital Investment / 10⁶ |
| AOC | 10⁶ $/yr | PreFerS | Annual Operating Cost / 10⁶ |
| GWP | kg CO₂-eq/kg | PreFerS | LCA displacement allocation (total per kg product) |

### Composition Metrics

| Metric | Units | Element | Description |
|--------|-------|---------|-------------|
| Salt Content | wt% | Composition | (NaCl + NaOH) / total mass × 100 |
| Residual CD | wt% | Composition | γ-Cyclodextrin / total mass × 100 |
| Residual Nicotinamide | wt% | Composition | Nicotinamide / total mass × 100 |
| Intermediate HemDx | wt% | Composition | HemoDextrin (intermediate) / total mass × 100 |
| N-HemoDextrin | wt% | Composition | N-HemoDextrin (final product) / total mass × 100 |
| Heme Equivalent | wt% | Composition | `(mol_NHemDx + mol_HemDx) × 0.0014711 × 616.487` / total × 100 |

---

## Parameter Count by Config

| Config | Fermentation | DSP | Facilities | Economics | GWP | Total |
|--------|-------------|-----|-----------|-----------|-----|-------|
| config1 | 5 | 6 | 1 | 7 | 4 | 23 |
| config2 | 5 | 5¹ | 1 | 7 | 4 | 22 |
| config3 | 5 | 4² | 1 | 7 | 4 | 21 |

¹ No S404 → filtration capture may have fewer units  
² No S402 (cell disruption), C403, S405 → fewer DSP params

---

## Summary Statistics (config1)

| Property | Value |
|----------|-------|
| Total parameters | 23 |
| Total metrics | 10 |
| Coupled parameters | 12 |
| Isolated parameters | 11 |
| Trunc(LogNormal) params | 9 (titer + 4 prices + 4 GWP) |
| TruncNormal params | 6 (2 yields + 4 DSP efficiencies + resin yield + CSTR) |
| Triangle params | 4 (tau + op days + IRR + tax) |
| Uniform params | 3 (secretion fraction + BT efficiency) |

### Verified Baselines (v2 config1)

| Scenario | MSP ($/kg) | GWP (kg CO₂-eq/kg) |
|----------|-----------|---------------------|
| Baseline | 8.86 | 8.29 |
| Lower Bound | 15.09 | — |
| Upper Bound | 7.38 | — |

## Last Updated
2026-02-07

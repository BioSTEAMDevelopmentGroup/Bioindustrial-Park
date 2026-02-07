# HemDx Uncertainty Model: Parameters & Metrics

**Source:** `v1/HemDx/_models.py`  
**Model Factory:** `create_model(baseline_production_kg_hr=150, config='config1', verbose=True)`  
**Organism:** *Corynebacterium glutamicum*  
**Product:** N-HemoDextrin (Nicotinamide-Stabilized Hemodextrin)

## Overview

The HemDx uncertainty model supports 3 configurations (config1/2/3) with up to 21 parameters across 5 categories and 10 metrics. Parameter count varies by config as some DSP units are absent in config2 (no S404) and config3 (no S402, C403, S405).

---

## Model Parameters

### 1. Fermentation Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 1.1 | Fermentation titer | `R302.titer` | g/L | Triangle | 0.5× baseline | 1.5× baseline | Sets `R302.titer` and `R302.target_titer` |
| 1.2 | Fermentation tau | `R302.tau` | hr | Triangle | 0.8× baseline | 1.2× baseline | Updates `R302.target_productivity` |
| 1.3 | Product yield | `R302.target_yield×100` | % | Triangle | 0.9× baseline | 1.1× baseline | Sets `R302.target_yield` |
| 1.4 | Biomass yield | CG2 growth yield×100 | % | Triangle | 0.9× baseline | 1.1× baseline | Updates both CG1 and CG2 reactions |
| 1.5 | Secretion fraction | 0.45 (config-dependent) | fraction | Uniform | Config-dependent | Config-dependent | Updates `R302.reaction_params['SF']` and redefines R302 spec |

#### Secretion Fraction Distribution by Config

| Config | Baseline SF | Distribution | Lower | Upper | Rationale |
|--------|------------|-------------|-------|-------|-----------|
| config1 (Base) | 0.45 | Uniform | 0.2 | 0.8 | Mixed secretion |
| config2 (Intracellular) | 0.10 | Uniform | 0.0 | 0.3 | Low secretion |
| config3 (Extracellular) | 0.90 | Uniform | 0.7 | 1.0 | High secretion |

### 2. DSP Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Configs | Notes |
|---|-----------|----------|-------|-------------|-------|-------|---------|-------|
| 2.1 | Centrifuge split | 1.0 | multiplier | Uniform | 0.95 | 1.0 | All | C401, C402, C403 (if present) |
| 2.2 | Cell disruption efficiency | 0.87 | - | Uniform | 0.82 | 0.92 | 1, 2 | S402 only (absent in config3) |
| 2.3 | Filtration capture | 0.95 | - | Uniform | 0.90 | 0.99 | Varies | S404, S405 (config-dependent) |
| 2.4 | Diafiltration retention | 0.95 | - | Uniform | 0.90 | 0.99 | All | U601, U801 |
| 2.5 | Resin yield | `U501.TargetProduct_Yield` (~0.99) | - | Uniform | 0.90 | 0.99 | All | U501 adsorption column |
| 2.6 | CSTR conversion | 0.95 | - | Uniform | 0.87 | 0.97 | All | R702 complexation reactions |

### 3. Facilities Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 3.1 | Turbogenerator efficiency | `BT.turbogenerator_efficiency` (~0.85) | - | Uniform | 0.80 | 0.90 | BT power generation |

### 4. Economic Parameters (Isolated)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 4.1 | Electricity price | `PowerUtility.price` | $/kWh | Triangle | 0.9× base | 1.1× base | ±10% |
| 4.2 | Glucose price | `glucose.price` | $/kg | Triangle | 0.9× base | 1.1× base | ±10% |
| 4.3 | Ammonia price | `ammonia.price` | $/kg | Triangle | 0.9× base | 1.1× base | ±10% |
| 4.4 | Buffer/Seed cost | 1.0 | multiplier | Uniform | 0.9 | 1.1 | Scales all buffer stream prices |
| 4.5 | Operating days | 333 | days/yr | Triangle | 300 | 350 | Also updates `operating_hours` |
| 4.6 | IRR | 0.18 | fraction | Triangle | 0.10 | 0.25 | Internal rate of return |
| 4.7 | Income tax | 0.17 | fraction | Triangle | 0.10 | 0.25 | Corporate income tax rate |

### 5. GWP Parameters (Isolated)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 5.1 | Electricity GWP | Grid emission factor | kg CO₂-eq/kWh | Triangle | 0.9× base | 1.1× base | Set as (val, val) tuple |
| 5.2 | Glucose GWP | `glucose.CF['GWP']` | kg CO₂-eq/kg | Triangle | 0.9× base | 1.1× base | |
| 5.3 | Ammonia GWP | `ammonia.CF['GWP']` | kg CO₂-eq/kg | Triangle | 0.9× base | 1.1× base | |
| 5.4 | Buffer/Seed GWP | 1.0 | multiplier | Uniform | 0.9 | 1.1 | Scales all buffer GWP CFs |

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

## Last Updated
2026-02-07

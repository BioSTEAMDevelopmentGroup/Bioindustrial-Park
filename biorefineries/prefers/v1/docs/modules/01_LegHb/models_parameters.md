# LegHb Uncertainty Model: Parameters & Metrics

**Source:** `v1/LegHb/_models.py`  
**Model Factory:** `create_model(baseline_production_kg_hr=150, config='config1', verbose=True)`  
**Organism:** *Pichia pastoris*

## Overview

The LegHb uncertainty model quantifies risk and sensitivity by sampling 17 parameters across 5 categories (Fermentation, DSP, Facilities, Economics, GWP) and evaluating 10 metrics (4 economic/environmental + 6 composition).

---

## Model Parameters

### 1. Fermentation Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 1.1 | Fermentation titer | `R302.titer` (~5.0) | g/L | Triangle | 0.5× baseline | 1.5× baseline | Sets `R302.titer` and `R302.target_titer` |
| 1.2 | Fermentation tau | `R302.tau` (~72) | hr | Triangle | 0.8× baseline | 1.2× baseline | Also updates `R302.target_productivity` |
| 1.3 | Product yield | `R302.target_yield×100` (~3.33) | % | Triangle | 0.9× baseline | 1.1× baseline | Sets `R302.target_yield` |
| 1.4 | Biomass yield | Cell growth yield×100 (~53) | % | Triangle | 0.9× baseline | 1.1× baseline | Updates `Pichia_pastoris` product yield |

### 2. DSP Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 2.1 | Centrifuge split | 1.0 | multiplier | Uniform | 0.95 | 1.0 | Applied to C401, C402, C403; clipped at 1.0 |
| 2.2 | Cell disruption efficiency | 0.87 | - | Uniform | 0.82 | 0.92 | S401 HPH efficiency |
| 2.3 | Filtration capture | 0.85 | - | Uniform | 0.75 | 0.95 | S403 solid capture |
| 2.4 | Diafiltration retention | 0.95 | - | Uniform | 0.90 | 0.99 | U501, U502 product retention |

### 3. Facilities Parameters (Coupled)

| # | Parameter | Baseline | Units | Distribution | Lower | Upper | Notes |
|---|-----------|----------|-------|-------------|-------|-------|-------|
| 3.1 | Turbogenerator efficiency | `BT.turbogenerator_efficiency` (~0.85) | - | Uniform | 0.80 | 0.90 | BT power generation efficiency |

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
| 5.1 | Electricity GWP | Grid emission factor | kg CO₂/kWh | Triangle | 0.9× base | 1.1× base | Set as (val, val) tuple |
| 5.2 | Glucose GWP | `glucose.CF['GWP']` | kg CO₂/kg | Triangle | 0.9× base | 1.1× base | |
| 5.3 | Ammonia GWP | `ammonia.CF['GWP']` | kg CO₂/kg | Triangle | 0.9× base | 1.1× base | |
| 5.4 | Buffer/Seed GWP | 1.0 | multiplier | Uniform | 0.9 | 1.1 | Scales all buffer GWP CFs |

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
| Total Solids | wt% | Composition | (total - H₂O) / total × 100 |
| Protein Purity | % | Composition | LegHb / (LegHb + Globin + Mannoprotein) × 100 |
| Heme Equivalent | wt% | Composition | Based on LegHb molar content: `mol/763 × MW_Heme_b` |

---

## Parameter Summary

- **Total Parameters:** 17 (4 fermentation + 4 DSP + 1 facilities + 7 economics + 4 GWP [note: 1 param not counted due to overlap])
- **Total Metrics:** 10 (4 economic/environmental + 6 composition)
- **Coupled Parameters:** 8 (require system re-simulation)
- **Isolated Parameters:** 9 (only affect TEA/LCA calculations)

## Last Updated
2026-02-07

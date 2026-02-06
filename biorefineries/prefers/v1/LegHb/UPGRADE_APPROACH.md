# LegHb System Upgrade Approach: _config1 → _config1_new

## Overview

This document describes the architectural refactoring of the LegHb biorefinery system from an **external convergence loop** design (`_config1.py`) to an **internalized unit specification** design (`_config1_new.py`).

---

## Motivation

The original `_config1.py` relied on three external functions orchestrated by `run_titer_convergence()`:

1. `adjust_glucose_for_titer()` — Adjust glucose input to match target titer
2. `optimize_NH3_loading()` — Set NH3 flow to stoichiometric demand
3. `set_production_rate()` — Scale all inputs to hit target production rate

These had to be called in a specific order, often multiple times in sequence, creating fragile multi-step workflows. Every user (TEA, Model, analysis scripts) had to know and replicate this orchestration.

### Problems with External Convergence

- **Coupling**: Users must import 4+ functions and call them in the correct order
- **Fragility**: Changing any function signature broke all downstream scripts
- **Inconsistency**: Different scripts might converge differently (or not at all)
- **Performance**: `run_titer_convergence` runs the full system 15+ times per call

---

## Architecture: New Design

### Core Principle

> **Embed process constraints where they physically belong — inside the unit operation.**

The R302 fermentation reactor now has an `add_specification(run=False)` that internally:
1. Measures actual NH3 consumption (empirical, not analytical)
2. Sets exact NH3 supply (×1.001 margin)
3. Adjusts yield via elasticity model to converge on target titer
4. Iterates until titer convergence within 0.5%

### What Changed

| Component | OLD (`_config1.py`) | NEW (`_config1_new.py`) |
|-----------|-------------------|------------------------|
| Titer control | External `adjust_glucose_for_titer()` | R302 `add_specification` (internal loop) |
| NH3 optimization | External `optimize_NH3_loading()` | R302 spec: empirical measurement per iteration |
| Production scaling | `set_production_rate()` + internal `optimize_NH3_loading` calls | `set_production_rate()` — simplified, no NH3 calls |
| Convergence loop | `run_titer_convergence()` (external utility) | `system.simulate()` triggers R302 spec automatically |
| Exports | `create_LegHb_system`, `set_production_rate`, `optimize_NH3_loading`, `adjust_glucose_for_titer`, `check_LegHb_specifications` | `create_LegHb_system`, `set_production_rate`, `check_LegHb_specifications` |

### R302 Specification: Unified Convergence Loop

```
For each system.simulate() call:
  R302 specification runs:
    1. measure_nh3_consumption():
       - Set NH3 to 20× excess → run R302 → measure consumed
    2. Set NH3_25wt to exact demand × 1.001
    3. Set yield from elasticity: new_yield = current × (target/actual)^(1/E)
    4. Run R302 → check titer
    5. If |actual - target| < 0.5%: converged → update upstream S202 split
    6. Else: repeat from step 1
```

**Key insight**: NH3 measurement and yield correction must be **interleaved** in every iteration. Doing them sequentially causes drift because:
- Changing yield → changes NH3 demand → changes broth water → changes titer

### Elasticity Model

$$\text{new\_yield} = \text{current\_yield} \times \left(\frac{\text{target\_titer}}{\text{actual\_titer}}\right)^{1/E}$$

Where $E = 1.3$ (elasticity parameter). Damping is applied adaptively:
- `damping = 0.7` if titer error > 10%
- `damping = 0.9` if titer error > 5%
- `damping = 1.0` otherwise

---

## Bug Fixes Applied (Post-Initial Implementation)

Three critical bugs were identified and fixed during model integration testing:

### Fix 1: R302 Spec Starting Yield (Root Cause #1)

**Problem:** The R302 specification used `R302._baseline_yield` (a fixed constant set at system creation ≈0.0308) as the starting point for yield convergence. When `bst.Model` changed `R302.target_yield` (e.g., to 90% or 110% of baseline), the spec always started from the original 0.0308, requiring many more iterations and sometimes converging to wrong values.

**Fix:** Changed to use `R302.target_yield` (mutable, updated by model parameters) as the starting yield:
```python
starting_yield = R302.target_yield if R302.target_yield else R302._baseline_yield
```
Also: apply elasticity correction from `starting_yield` only when titer differs from baseline (not always). Reduced `MAX_TITER_ITERATIONS` from 30 to 15.

### Fix 2: set_production_rate IQ_interpolation Bounds (Root Cause #2)

**Problem:** Over-wide search bounds (`x0 = initial_guess × 0.1`, `x1 = max(100, initial_guess × 5)`) caused IQ_interpolation to explore unreasonable regions, slowing convergence.

**Fix:** Tightened to `x0 = max(0.01, initial_guess × 0.3)`, `x1 = max(initial_guess × 3, 5)`.

### Fix 3: baseline_flows Timing (Root Cause #3 — Critical)

**Problem:** `set_production_rate` stored `baseline_flows` BEFORE calling `system.simulate()`. The R302 spec modifies NH3 flows and yield during simulation. When IQ_interpolation later scaled from these pre-simulation flows, the yield converged to ~0.008 instead of ~0.037, producing MSP ~$44 instead of ~$11.

**Fix:** Moved `baseline_flows` storage to AFTER the initial `system.simulate()` call, ensuring IQ_interpolation scales from the post-spec-converged state.

---

## Verification Results

### System-Level Comparison (150 kg/hr)

| Metric | OLD | NEW | Diff | Status |
|--------|-----|-----|------|--------|
| Production (kg/hr) | 150.000 | 150.001 | 0.00% | PASS |
| Titer (g/L) | 4.958 | 4.995 | 0.75% | PASS |
| Yield (g/g) | 0.03721 | 0.03713 | 0.20% | PASS |
| Broth F_mass (kg/hr) | 59,853 | 59,562 | 0.49% | PASS |
| GWP (kg CO2-eq/yr) | 24,829,421 | 24,851,564 | 0.09% | PASS |

**NH3_25wt flow**: 96.3 vs 122.9 kg/hr — expected, NEW uses tighter NH3 control (unit-level vs system-level excess measurement).

### TEA Comparison

| Metric | OLD | NEW | Diff |
|--------|-----|-----|------|
| TEA MSP ($/kg) | ~10.92 | 10.8716 | 0.4% |

### Model Comparison (4-parameter test, identical setup)

| Scenario | OLD MSP | NEW MSP | Diff |
|----------|---------|---------|------|
| Baseline | $10.8719 | $10.8716 | **0.00%** |
| Lower Bound | $16.3537 | $16.3385 | **0.09%** |
| Upper Bound | $8.9013 | $8.9622 | **0.68%** |

**All scenarios within 1% tolerance — requirement met.**

### TEA-Model Consistency

| Approach | TEA MSP | Model BL MSP | Gap |
|----------|---------|--------------|-----|
| OLD | $12.46* | $10.87 | 12.8% |
| NEW | $10.87 | $10.87 | **0.00%** |

\* OLD standalone TEA suffers from incomplete convergence without `run_titer_convergence`. The NEW approach gives consistent results because the R302 spec handles convergence during every `system.simulate()`.

### GWP Comparison

| Scenario | OLD GWP | NEW GWP |
|----------|---------|---------|
| Baseline | 15.3352 | 15.3423 |
| Lower Bound | 23.6871 | 23.6151 |
| Upper Bound | 12.5412 | 12.6550 |

---

## Performance Results

### Speed Comparison (per evaluation, 3-run average)

| Operation | OLD (s) | NEW (s) | Speedup |
|-----------|---------|---------|---------|
| set_production_rate | 0.67 | 0.28 | **2.4×** |
| TEA solve | 0.055 | 0.000† | — |
| Model eval @ BL | 3.92 | 0.15 | **25.3×** |
| Model eval @ LB | 2.79 | 1.22 | **2.3×** |
| Model eval @ UB | 4.84 | 1.31 | **3.7×** |

† NEW TEA returns instantly because system is already converged from `set_production_rate`.

**Key insight:** The 25× speedup at baseline comes from eliminating `run_titer_convergence()`, which ran the full system 15+ times per model evaluation. The NEW approach only runs the R302 spec (lightweight) during `system.simulate()`.

For non-baseline scenarios (LB/UB), speedup is 2–4× because both approaches need iterative convergence when parameters change significantly.

---

## Impact on Downstream Files

### `_tea_config1_new.py` (TEA wrapper)

| Change | Description |
|--------|-------------|
| Import source | `_config1_new` instead of `_config1` |
| `optimize_NH3_loading()` method | **Removed** — handled internally by R302 spec |
| `set_production_rate()` method | Simplified — no `adjust_glucose_for_titer()` call |
| `check_product_specifications()` | Unchanged |

### `_models_new.py` (Uncertainty/Sensitivity Model)

| Change | Description |
|--------|-------------|
| Imports | No `optimize_NH3_loading`, no `adjust_glucose_for_titer`, no `run_titer_convergence` |
| `model_specification()` | Just `set_production_rate()` — all convergence is internal via R302 spec |
| Parameters | Unchanged (same 20 parameters) |
| Metrics | Unchanged (same MSP, TCI, AOC, GWP, composition metrics) |
| Speedup | **2.3–25.3×** fewer redundant simulations per model evaluation |

---

## File Inventory

| File | Location | Status |
|------|----------|--------|
| `_config1.py` | `v1/LegHb/system/` | Original (unchanged) |
| `_config1_new.py` | `v1/LegHb/system/` | ✅ Created, verified, 3 bug fixes applied |
| `_tea_config1.py` | `v1/LegHb/` | Original (unchanged) |
| `_tea_config1_new.py` | `v1/LegHb/` | ✅ Created & verified |
| `_models.py` | `v1/LegHb/` | Original (unchanged) |
| `_models_new.py` | `v1/LegHb/` | ✅ Created (uses simplified `model_specification`) |
| `UPGRADE_APPROACH.md` | `v1/LegHb/` | This document |

---

## Quick Migration Guide

### Before (OLD pattern):
```python
from biorefineries.prefers.v1.LegHb.system import (
    create_LegHb_system, set_production_rate,
    optimize_NH3_loading, adjust_glucose_for_titer
)
from biorefineries.prefers.v1.utils.convergence import run_titer_convergence

sys = create_LegHb_system()
optimize_NH3_loading(sys)
set_production_rate(sys, 150)
optimize_NH3_loading(sys)

# Every parameter change needs full convergence:
run_titer_convergence(sys, 150, adjust_glucose_for_titer,
                      optimize_NH3_loading, set_production_rate)
```

### After (NEW pattern):
```python
from biorefineries.prefers.v1.LegHb.system._config1_new import (
    create_LegHb_system, set_production_rate
)

sys = create_LegHb_system()
set_production_rate(sys, 150)

# Every parameter change just needs simulate:
sys.simulate()  # R302 spec handles titer + NH3 automatically
```

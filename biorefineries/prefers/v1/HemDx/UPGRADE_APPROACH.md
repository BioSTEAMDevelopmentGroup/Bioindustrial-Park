# HemDx Upgrade Approach: Internalized Specifications (_config_new)

## Overview

This document describes the upgrade from the OLD external convergence approach to the NEW internalized `add_specification` pattern for the HemDx N-HemoDextrin production system.

## Problem Statement

The OLD approach used three separate external functions called in a convergence loop:
1. `adjust_glucose_for_titer()` — adjusted yield based on titer target
2. `optimize_NH3_loading()` — set NH3 flow and S202 split to match demand
3. `run_titer_convergence()` — outer loop calling both + `set_production_rate()`

This required **12+ outer-loop iterations** per evaluation, each calling `set_production_rate` which itself ran `optimize_NH3_loading` (2 system.simulate() calls per NH3 optimization). Result: **~8-10 seconds per baseline evaluation**.

## Solution: Internalized R302 Specification

All titer control and NH3 optimization logic is now embedded in a single `@R302.add_specification(run=False)` that executes **during** `system.simulate()`:

```python
@R302.add_specification(run=False)
def update_fermentation_and_nh3():
    # 1. Measure NH3 demand empirically (20x excess → measure consumption)
    # 2. Set exact NH3 on R302's inlet
    # 3. Run R302 and check titer
    # 4. Adjust yield using elasticity correction (ELASTICITY=2.0)
    # 5. Repeat until convergence (TITER_TOL=0.005, max 15 iterations)
    # 6. Update upstream NH3_25wt source and S202 split
```

### Key Design Decisions (HemDx-specific)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| ELASTICITY | 2.0 | Higher than LegHb (1.3) — HemDx oscillates more due to lower yields |
| TITER_TOL | 0.005 | 0.5% relative tolerance (sufficient for MSP <1% match) |
| Yield clamps | [0.0001, 0.05] | Much lower than LegHb [0.005, 0.15] due to Y_Heme ≈ 0.00175 |
| SF handling | Dynamic via `R302.reaction_params['SF']` | Model parameter can change SF without replacing spec |
| 4-reaction update | All yields split by SF | Heme_b(SF), Heme_b_In(1-SF), ProtoporphyrinIX(SF*Y_pp), ProtoporphyrinIX_In((1-SF)*Y_pp) |

### Three Critical Fixes (from LegHb session, applied to HemDx)

1. **Starting yield = R302.target_yield** (mutable, respects model parameters) — not fixed `_baseline_yield`
2. **IQ_interpolation bounds tightened**: x0=0.3×guess, x1=3×guess (was 0.1× and 10×)
3. **baseline_flows stored AFTER system.simulate()** — captures post-spec stream state

## Files Created

| File | Purpose |
|------|---------|
| `system/_config1_new.py` | Config 1 (SF=0.45) with internalized specs |
| `system/_config2_new.py` | Config 2 (SF=0.1, intracellular) with internalized specs |
| `system/_config3_new.py` | Config 3 (SF=0.9, extracellular) with internalized specs |
| `_tea_config1_new.py` | TEA wrapper for config1_new |
| `_tea_config2_new.py` | TEA wrapper for config2_new |
| `_tea_config3_new.py` | TEA wrapper for config3_new |
| `_models_new.py` | Unified model supporting all 3 configs via _new modules |

## Verification Results

### Config1 TEA Baseline (150 kg/hr)
| Metric | OLD | NEW | Diff |
|--------|-----|-----|------|
| MSP | $9.1903/kg | $9.1958/kg | 0.06% |
| Time | 9.16s | 0.52s | **17.6x faster** |

### Config1 Model (BL/LB/UB)
| Scenario | OLD MSP | NEW MSP | Diff | Speedup |
|----------|---------|---------|------|---------|
| Baseline | $8.7596 | $8.7658 | 0.07% | **36.1x** |
| LB | $15.6798 | $15.6995 | 0.13% | **6.4x** |
| UB | $7.3873 | $7.3679 | 0.26% | **6.0x** |

### Config2/3 Model Baseline
| Config | OLD MSP | NEW MSP | Diff | Speedup |
|--------|---------|---------|------|---------|
| Config2 | $8.7530 | $8.7636 | 0.12% | **33.2x** |
| Config3 | $7.2022 | $7.2069 | 0.07% | **33.9x** |

**All outputs within 1% of OLD. All configs pass.**

## Architecture Comparison

### OLD Flow
```
model_specification() →
  run_titer_convergence() →
    for i in range(15):
      adjust_glucose_for_titer()    ← separate function
      optimize_NH3_loading()        ← 2x system.simulate() per call
      set_production_rate()         ← calls optimize_NH3_loading again inside
    check convergence
```

### NEW Flow
```
model_specification() →
  set_production_rate() →
    system.simulate() →
      R302.update_fermentation_and_nh3()  ← all-in-one spec
        measure_nh3(excess) → set_exact → run → check_titer → adjust_yield
```

## Backward Compatibility

- OLD files (`_config1.py`, `_tea_config1.py`, `_models.py`) remain untouched
- NEW files are additive (`_config1_new.py`, etc.)
- `_models_new.py` supports config switching: `create_model(config='config1'|'config2'|'config3')`

## Date
2026-02-03

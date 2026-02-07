# System Optimization & Convergence Logic

This document details the optimization methodologies used in the PreFerS v2 biorefinery models (LegHb and HemDx) to ensure physical feasibility and meet design specifications.

## Overview: Internalized Specification Architecture

v2 replaces the external `convergence.py` utility with **internalized unit specifications** attached directly to the fermentation reactor (R302) via `@R302.add_specification(run=False)`. This means titer convergence and NH3 optimization execute automatically whenever BioSTEAM simulates the system — no external wrapper function is needed.

**Key advantages over the v1 approach:**
- **3–45x faster** simulation (no redundant outer-loop system simulates)
- **Simpler model specification**: `model_specification()` just calls `set_production_rate()` + `system.simulate()`
- **No external convergence imports** needed in config or model files
- **Self-contained**: each config file is fully standalone

## Titer Convergence (R302 Specification)

The standard design specification requires achieving a specific **target titer** (g/L) and **target production rate** (kg/hr). Because titer, yield, and production volume are coupled non-linearly (due to cell growth, density changes, and downstream recoveries), a simple single-step calculation is insufficient.

### Architecture

The convergence logic is embedded as `@R302.add_specification(run=False)`, which means:
1. BioSTEAM calls this spec *instead of* `R302._run()` during simulation
2. The spec internally calls `R302._run()` within its iteration loop
3. No external `run_titer_convergence()` call is needed

### Logic Flow

The `update_fermentation_and_nh3()` specification (attached to R302) executes the following loop until the actual titer matches the target titer or max iterations (15) are reached:

1.  **Initialize Yield from Elasticity Correction**:
    *   If `target_titer` differs from `baseline_titer`, applies an initial correction:
    *   $Yield_{init} = Yield_{starting} \times (\frac{Titer_{target}}{Titer_{baseline}})^{\frac{1}{Elasticity}}$
    *   LegHb: $Elasticity = 1.3$, $Tolerance = 0.1\%$
    *   HemDx: $Elasticity = 2.0$, $Tolerance = 0.5\%$

2.  **Measure NH3 Demand Empirically** (per iteration):
    *   Runs R302 with 20× excess NH3 at the current yield
    *   Measures actual NH3 consumption: $consumed = \sum NH3_{in} - \sum NH3_{out}$
    *   This replaces the old stoichiometric `optimize_NH3_loading` function

3.  **Set Exact NH3 and Run R302**:
    *   Sets the ammonia inlet to 100.1% of measured consumption (minimal safety margin)
    *   Runs R302 with correct yield and NH3

4.  **Check Titer Convergence**:
    *   Calculates relative error: $error = |Titer_{target} - Titer_{actual}| / Titer_{target}$
    *   If $error < Tolerance$, loop terminates successfully

5.  **Adjust Yield with Adaptive Damping**:
    *   Coarse mode (error > 10%): full elasticity step ($damping = 1.0$)
    *   Fine mode (error ≤ 10%): half step ($damping = 0.5$)
    *   $step = (\frac{Titer_{target}}{Titer_{actual}})^{\frac{1}{Elasticity} \times damping}$
    *   $Yield_{new} = Yield_{current} \times step$, clamped to $[0.005, 0.15]$ (LegHb) or $[0.0001, 0.05]$ (HemDx)

6.  **Update Upstream NH3 Source**:
    *   After convergence, updates the S202 splitter and NH3_25wt source stream to reflect total NH3 demand (R301 + R302) for downstream consistency

## Production Rate Scaling

The `set_production_rate(system, target_kg_hr)` function uses `flexsolve.IQ_interpolation` (root-finding) to find the exact input stream scaling factor that achieves the target production rate in the final product stream. This is necessary because non-linear separations (membrane rejection, centrifuge recovery) make simple linear scaling inaccurate.

## Usage in System Configs

The internalized specification is implemented in:
*   `LegHb/system/_config1.py` — `update_fermentation_and_nh3()` on R302
*   `HemDx/system/_config1.py` — `update_fermentation_and_nh3()` on R302
*   `HemDx/system/_config2.py` — `update_fermentation_and_nh3()` on R302
*   `HemDx/system/_config3.py` — `update_fermentation_and_nh3()` on R302

Each config is fully self-contained with no imports from `utils/convergence.py` (which has been removed in v2).

## Usage in Models (Uncertainty Analysis)

In `_models.py`, the `model_specification` function is simplified to:

```python
def model_specification():
    set_production_rate(system, baseline_production_kg_hr, verbose=False)
```

The R302 specification handles titer/NH3 convergence automatically during `system.simulate()`. This guarantees that for *every* Monte Carlo sample — where parameters like yield, efficiency, and prices change — the system physically converges to the correct titer and production rate before metrics are calculated.

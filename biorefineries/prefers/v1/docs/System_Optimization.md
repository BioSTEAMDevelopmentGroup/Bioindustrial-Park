# System Optimization & Convergence Logic

This document details the iterative optimization methodologies used in the PreFerS biorefinery models (LegHb and HemDx) to ensure physical feasibility and meet design specifications.

## Titer Convergence Loop

The standard design specification for `PreFerS` systems requires achieving a specific **target titer** (g/L) and **target production rate** (kg/hr). Because titer, yield, and production volume are coupled non-linearly (due to cell growth, density changes, and downstream recoveries), a simple single-step calculation is insufficient.

We employ a cohesive **Iterative Titer Convergence** logic, encapsulated in `biorefineries/prefers/v1/utils/convergence.py`.

### Logic Flow

The `run_titer_convergence` utility executes the following loop until the actual titer matches the target titer (within 1% tolerance) or max iterations (15) are reached:

1.  **Adjust Titer (Yield) Correction**:
    *   Compares the `actual_titer` from the previous simulation against the `target_titer`.
    *   Uses an elasticity-based correction model to adjust the fermentation **yield** parameter.
    *   Equation: $Yield_{new} = Yield_{old} \times (\frac{Titer_{target}}{Titer_{actual}})^{\frac{1}{Elasticity}}$

2.  **Optimize Nutrient Loading (Ammonia)**:
    *   Calculates the required stoichiometric Nitrogen for the new biomass and product yield.
    *   Adjusts the `NH3` feed rate to ensure nitrogen is not limiting, nor excessively over-supplied (which would impact wastewater treatment).
    *   *Note*: This optimization is critical to perform *before* scaling flow rates.

3.  **Set Production Rate**:
    *   Scales all input streams proportionally to match the `Target Production Rate` (kg/hr) of the final product stream.
    *   Uses a `flexsolve` root-finding algorithm to find the exact scaling factor, as non-linear separations (like membrane rejection) can make linear scaling inaccurate.

4.  **Re-Optimize Nutrients**:
    *   After scaling, minor adjustments to Ammonia may be needed to ensure precise stoichiometry at the new flow rates.

5.  **Simulate System**:
    *   Runs a full mass-and-energy balance simulation (`sys.simulate()`).

6.  **Check Convergence**:
    *   Calculates relative error of Titer.
    *   If error < 1%, loop terminates successfully.
    *   If not, repeats from Step 1 with the newly simulated data.

## Usage in System Configs

This logic is implemented in:
*   `LegHb/system/_config1.py`
*   `HemDx/system/_config1.py`
*   `HemDx/system/_config2.py`
*   `HemDx/system/_config3.py`

And utilized in both checking scripts (`if name == main`) and `_models.py` for uncertainty analysis.

## Usage in Models (Uncertainty Analysis)

In `_models.py`, the `model_specification` function wraps this convergence loop. This guarantees that for *every* Monte Carlo sample—where parameters like yield, efficiency, and prices change—the system physically converges to the correct titer and production rate before metrics are calculated. This prevents "physical infeasibility" noise in sensitivity analysis results.

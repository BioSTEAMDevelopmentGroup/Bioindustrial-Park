# Uncertainty & Sensitivity Analysis (`LegH`)

**Directory:** `v1/LegH/analyses/`  
**Main Script:** `uncertainty_and_sensitivity.py`

This module provides comprehensive tools for quantifying risk and parameter sensitivity in the LegHemoglobin biorefinery.

## 1. Analysis Strategy

The analysis is divided into three scenarios, handled by the `run_monte_carlo` driver:
1.  **Fixed Scale Monte Carlo:** Runs simulations while keeping the production scale fixed/constant. Used to isolate the effect of biological and process uncertainties (yield, titer) on economics and sustainability.
    *   *Outputs:* Spearman Correlations, KDE Plots.
2.  **Variable Scale Monte Carlo:** Varies the production scale alongside other parameters to study economies of scale.
3.  **Single-Point Sensitivity:** (Planned) Generates Tornado diagrams.

## 2. Key Uncertain Parameters

The model (`bst.Model`) samples from defined distributions (Uniform, Triangular, etc.). Key parameters include:

*   **Fermentation Performance:**
    *   `Titer` [g/L]
    *   `Yield` [g product / g glucose]
    *   `Productivity` [g/L/hr]
*   **Operating Costs:**
    *   `Electricity Price` [$/kWh]
    *   `Substrate Price` (Glucose) [$/kg]
*   **CapEx Factors:**
    *   `Lang Factor` (Installation cost multiplier)

## 3. Metrics Evaluated

The analysis tracks the following Key Performance Indicators (KPIs):

*   **Economic:**
    *   `MSP`: Minimum Selling Price [$/kg LegH]
    *   `TCI`: Total Capital Investment [$ Million]
    *   `AOC`: Annual Operating Cost [$ Million/yr]
*   **Environmental:**
    *   `GWP`: Global Warming Potential [kg CO2-eq/kg LegH]
*   **Technical:**
    *   `LegH Content`: Final product purity [%]

## 4. Outputs

Running `uncertainty_and_sensitivity.py` generates:

1.  **Excel Datasets:**
    *   `LegH_MC_no_scale_*.xlsx`: Results with fixed scale.
    *   `LegH_MC_with_scale_*.xlsx`: Results with variable scale.

2.  **Visualizations:**
    *   **Scatter Plots (2D):** Color-coded scatter plots (e.g., MSP vs Titer/Yield).
    *   **KDE Plots:** 1D and 2D Kernel Density Estimates showing probability distributions of MSP and GWP.
    *   **Spearman Plots:** Tornado plots ranking parameters by their correlation coefficient (impact) on MSP and GWP.

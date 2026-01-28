# Project Status: Biorefinery Simulation Platform (PREFERS v1)

> **Date:** 2026-01-28  
> **Version:** v1.0  
> **Focus:** Leghemoglobin (LegHb) & Hemedextrin (HemDx)

---

## 1. Project Goal
*   **Objective**: Develop rigorous techno-economic analysis (TEA) and life cycle assessment (LCA) models for recombinant protein production in *Pichia pastoris*.
*   **Target Products**:
    1.  **Leghemoglobin (LegHb)**: Heme-protein for plant-based meat alternatives.
    2.  **Hemedextrin (HemDx)**: Soluble, protein-free heme-iron supplement.
*   **Technology Stack**: Python-based simulations using [BioSTEAM](https://biosteam.readthedocs.io).

## 2. Current Architecture
A modular, object-oriented framework is established at `biorefineries/prefers/v1/`:

*   **Unit Library (`v1/_units.py`)**: Custom `bst.Unit` subclasses for detailed modeling of:
    *   *Perfusion/Fed-Batch Fermentation*
    *   *High-Pressure Homogenization* (Cell Disruption)
    *   *Tangential Flow Filtration* (Ultrafiltration/Diafiltration)
*   **Data Interface (`v1/_process_settings.py`)**: "Single Source of Truth" for:
    *   Utility Prices (Electricity, Natural Gas determined by region).
    *   LCA characterization factors (GWP100).
*   **Generalized TEA (`PreFerSTEA`)**: Flexible economic model allowing rapid switching between different biorefinery configurations.

## 3. Progress Update (Current Stage)

### ✅ Completed Milestones
*   **LegHb Implementation**:
    *   **Config 1 (Food-Grade)**: HTST pasteurization + Ultrafiltration flowsheets verified.
    *   **Model Integration**: `optimize_NH3_loading` function successfully coupled with system design.
    *   **Refinements**: Updated water management and stream naming conventions for clarity.
*   **HemDx Implementation**:
    *   **Feasibility Fix**: Resolved specific mass balance issues (FeSO4 limitation) at 1500 kg/hr scale via wastewater treatment supplementation.
    *   **Base Topology**: Established baseline fermentation and recovery flowsheet.
*   **Documentation**:
    *   Synchronized `v1/docs` with latest script changes to ensure architectural accuracy.

### 🔄 In Progress / Refinement
*   **Model Robustness**: Final verifications of system integration scripts (`test_models_integration.py`).
*   **Visualization**: Refining TEA breakdown plots and Spearman correlation heatmaps in `gen_figure.py`.

## 4. Next Steps (Immediate Plan)

1.  **Verification Run**:
    *   Execute `verify_all_scaling.py` to confirm mass balance and TEA convergence across all target production capacities (500 - 10,000 kg/yr).
2.  **Uncertainty Analysis**:
    *   Run `LegHb/analyses/gen_data_mc.py` to generate Monte Carlo usage data.
    *   Produce Tornado plots and Sensitivity Heatmaps to identify key cost drivers.
3.  **New Design Phase (N-HemDx)**:
    *   Begin detailed implementation of the **"Split-Stream Heme Recovery"** topology.
    *   Scale-up design based on the "Basis of Design" document.
4.  **Reporting**:
    *   Compile generated figures and sensitivity results into the final deck.

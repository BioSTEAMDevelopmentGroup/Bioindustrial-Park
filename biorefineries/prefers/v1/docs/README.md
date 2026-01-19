# Developer Handbook: PREFERS Source Code (v1)

**Scope:** `biorefineries/prefers/v1/`  
**Version:** v1.0 (Production)

## Project Overview
This codebase implements the "Process Evaluation and Feedstock Engineering for Renewable Sustainability" (PREFERS) platform using **BioSTEAM**. It models biorefinery flowsheets for products like Leghemoglobin (`LegHb`) and Hemedextrin (`HemDx`).

## Directory Map

| Directory                 | Description                                                                                                   |
| :------------------------ | :------------------------------------------------------------------------------------------------------------ |
| `v1/LegHb/`                | **Leghemoglobin Module**. Contains `_system.py` (flowsheet), `_chemicals.py`, `_tea.py`, and `analyses/`.     |
| `v1/HemDx/`              | **HemDx Module**. Protein-free Heme B production logic.                                                      |
| `v1/_units.py`            | **Process Library**. Shared library of custom Unit Operations (SeedTrain, Disruption, Ultrafiltration, etc.). |
| `v1/_process_settings.py` | **Data Interface**. Global Truth for Utility Prices and LCA Factors.                                          |

## Documentation Index

1.  [Process Flow & Control (`LegHb`)](modules/01_LegHb_Process.md) - Detailed sequence of unit operations and control loops.
2.  [Unit Operations Catalog](Unit_Operations_Catalog.md) - Design parameters, sizing logic, and costing models for all custom units.
3.  [UASA Analysis Tools](modules/UASAforLegHb.md) - Guide to Uncertainty and Sensitivity Analysis scripts.

## Quick Start

### 1. Run the LegHb Flowsheet
To execute the Leghemoglobin simulation and see the results:

```bash
python v1/LegHb/_system.py
```

### 2. Run the TEA (Techno-Economic Analysis)
To run the full economic analysis with design specifications:

```bash
python v1/LegHb/_tea.py
```

### 3. Run Uncertainty Analysis
To perform Monte Carlo simulations:

```bash
python v1/LegHb/analyses/uncertainty_and_sensitivity.py
```

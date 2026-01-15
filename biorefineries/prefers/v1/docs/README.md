# Developer Handbook: PREFERS Source Code (v1)

**Scope:** `biorefineries/prefers/v1/`
**Version:** v1.0 (Production)

## Project Overview
This codebase implements the "Process Evaluation and Feedstock Engineering for Renewable Sustainability" (PREFERS) platform using **BioSTEAM**. It allows for modular construction of biorefinery flowsheets for products like Leghemoglobin (LegH) and Hemedextrin (HemeIn).

## Directory Map

| Directory                | Description                                                                                  |
| :----------------------- | :------------------------------------------------------------------------------------------- |
| `v1/LegH/`               | **Leghemoglobin Module**. Contains `system.py` (flowsheet), `chemicals.py`, and `tea.py`.    |
| `v1/HemeIn/`             | **HemeIn Module**. Protein-free Heme B production logic.                                     |
| `v1/units.py`            | **Process Library**. Shared library of custom Unit Operations (SeedTrain, Disruption, etc.). |
| `v1/process_settings.py` | **Data Interface**. Global Truth for Utility Prices and LCA Factors.                         |
| `v1/docs/`               | This documentation.                                                                          |

## Quick Start

### 1. Run the LegH Flowsheet
To execute the Leghemoglobin simulation and see the results:

```bash
python v1/LegH/system.py
```

### 2. Run the HemeIn Flowsheet
To execute the HemeIn simulation:

```bash
python v1/HemeIn/system.py
```

### 3. Understanding the Output
The scripts will print:
- Variable Operating Costs (VOC)
- CAPEX Breakdown
- Unit Operation Results
- TEA Summary (MPSP, NPV)

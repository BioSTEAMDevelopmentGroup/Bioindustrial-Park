# Developer Handbook: PREFERS Source Code (v2)

**Scope:** `biorefineries/prefers/v2/`  
**Version:** v2 (Production)

## Project Overview
This codebase implements the "Process Evaluation and Feedstock Engineering for Renewable Sustainability" (PREFERS) platform using **BioSTEAM**. It models biorefinery flowsheets for products like Leghemoglobin (`LegHb`) and Hemedextrin (`HemDx`).

## Directory Map

| Directory                 | Description                                                                                                                     |
| :------------------------ | :------------------------------------------------------------------------------------------------------------------------------ |
| `LegHb/`                  | **Leghemoglobin Module**. Contains `system/`, `_chemicals.py`, `_tea_config1.py`, `_models.py`, and `analyses/`.               |
| `LegHb/system/_config1.py`| Config 1 flowsheet (internalized titer + NH3 specs).                                                                           |
| `LegHb/analyses/`         | Uncertainty and sensitivity analysis scripts (split data + figures).                                                           |
| `HemDx/`                  | **HemDx Module**. N-HemoDextrin production logic with configs 1-3.                                                             |
| `HemDx/system/`           | Configurations: `_config1.py` (base), `_config2.py` (intracellular), `_config3.py` (extracellular).                            |
| `_units.py`               | **Process Library**. Shared custom Unit Operations (SeedTrain, Disruption, Diafiltration, ResinColumn, etc.).                 |
| `_units_adv.py`           | **Advanced Units**. Mechanistic Resin/Diafiltration/Filtration models.                                                         |
| `_process_settings.py`    | **Data Interface**. Global prices, utilities, and LCA factors.                                                                  |
| `_tea.py`                 | **TEA Base Class**. `PreFerSTEA` with IRAS depreciation schedules.                                                              |
| `utils/`                  | **Utilities & Plotting**. Plotting styles, helper functions, and file I/O (includes original `PreFerS` palette + `PreFerS_softlight` suffix variant). |

## Documentation Index

### LegHb Module
1.  [Process Flow Config 1](modules/01_LegHb/config1.md) - Intracellular recovery with HTST pasteurization
2.  [Uncertainty & Sensitivity Analysis](modules/UASAforLegHb.md) - Monte Carlo workflow (split data + figures)
3.  [Model Parameters](modules/01_LegHb/models_parameters.md) - Uncertainty model parameters and metrics

### HemDx Module
4.  [Process Flow Config 1](modules/02_HemDx/config1.md) - Integrated full system (base secretion)
5.  [Process Flow Config 2](modules/02_HemDx/config2.md) - Intracellular variant (low secretion)
6.  [Process Flow Config 3](modules/02_HemDx/config3.md) - Extracellular variant (high secretion)
7.  [Model Parameters](modules/02_HemDx/models_parameters.md) - Uncertainty model parameters and metrics

### Shared Resources
8.  [Unit Operations Catalog](Unit_Operations_Catalog.md) - Design parameters, sizing logic, and costing models for all custom units
9.  [System Architecture](System_Architecture.md) - High-level system design
10. [System Optimization](System_Optimization.md) - Internalized titer + NH3 convergence
11. [Data Interface](Data_Interface.md) - Process settings and LCA factors
12. [TEA Framework](modules/_tea.md) - `PreFerSTEA` base class and reporting
13. [Utilities & Plotting](../utils/README.md) - Plotting suite and helper functions

## Quick Start

### 1. Run the LegHb Flowsheet (Config 1 - Default)
To execute the Leghemoglobin simulation with food-grade configuration:

```bash
python -c "from biorefineries.prefers.v2.LegHb.system._config1 import create_LegHb_system; sys = create_LegHb_system(); sys.simulate(); sys.show()"
```

Or run the config script directly:
```bash
python biorefineries/prefers/v2/LegHb/system/_config1.py
```

### 2. Run the HemDx Flowsheet (Config 1)
```bash
python biorefineries/prefers/v2/HemDx/system/_config1.py
```

### 3. Run the TEA (Techno-Economic Analysis)
To run the full economic analysis with design specifications:

```bash
python -m biorefineries.prefers.v2.LegHb._tea_config1 --production 150
```

HemDx TEA:
```bash
python -m biorefineries.prefers.v2.HemDx._tea_config1 --production 150
```

### 4. Run Uncertainty & Sensitivity Analysis
The UA/SA workflow is split into data generation and plotting:

```bash
# Generate baseline + Monte Carlo data
python biorefineries/prefers/v2/LegHb/analyses/gen_data_base.py --config config1
python biorefineries/prefers/v2/LegHb/analyses/gen_data_mc.py --config config1 --samples 500

# Generate figures from saved data
python biorefineries/prefers/v2/LegHb/analyses/gen_figure.py --results-dir biorefineries/prefers/v2/LegHb/analyses/results_config1_YYYYMMDD_HHMM
```

## Configuration Notes

- LegHb ships with Config 1 in v2 (internalized titer + NH3 specs in `R302`).
- HemDx supports Configs 1-3 (base, intracellular, extracellular).

## Configuration Comparison

| Feature             | Config 1 (v2)                |
| ------------------- | ---------------------------- |
| Primary Separation  | Biomass harvest + cell wash  |
| Chromatography      | None                         |
| Polishing           | HTST Pasteurization (72°C)   |
| Final Formulation   | Antioxidant + dilution       |
| Product Temperature | 4°C                          |
| Debris Handling     | Boiler feed (via BT)         |

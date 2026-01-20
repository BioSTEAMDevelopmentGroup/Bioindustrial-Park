# Developer Handbook: PREFERS Source Code (v1)

**Scope:** `biorefineries/prefers/v1/`  
**Version:** v1.0 (Production)

## Project Overview
This codebase implements the "Process Evaluation and Feedstock Engineering for Renewable Sustainability" (PREFERS) platform using **BioSTEAM**. It models biorefinery flowsheets for products like Leghemoglobin (`LegHb`) and Hemedextrin (`HemDx`).

## Directory Map

| Directory                 | Description                                                                                                   |
| :------------------------ | :------------------------------------------------------------------------------------------------------------ |
| `v1/LegHb/`                | **Leghemoglobin Module**. Contains `system/` (flowsheet configs), `_chemicals.py`, `_tea.py`, `_models.py`, and `analyses/`. |
| `v1/LegHb/system/`         | System configurations: `_config1.py` (food-grade HTST), `_config2.py` (research-grade with IX). |
| `v1/LegHb/analyses/`       | Uncertainty and sensitivity analysis scripts. |
| `v1/HemDx/`              | **HemDx Module**. Protein-free Heme B production logic.                                                      |
| `v1/_units.py`            | **Process Library**. Shared library of custom Unit Operations (SeedTrain, Disruption, Ultrafiltration, etc.). |
| `v1/_process_settings.py` | **Data Interface**. Global Truth for Utility Prices and LCA Factors.                                          |

## Documentation Index

### LegHb Module
1.  [Process Flow Config 1 (Food-Grade)](modules/01_LegHb/config1.md) - Intracellular recovery with HTST pasteurization
2.  [Process Flow Config 2 (Research-Grade)](modules/01_LegHb/config2.md) - Traditional downstream with Ion Exchange
3.  [Uncertainty & Sensitivity Analysis](modules/UASAforLegHb.md) - Monte Carlo simulation tools

### Shared Resources
4.  [Unit Operations Catalog](Unit_Operations_Catalog.md) - Design parameters, sizing logic, and costing models for all custom units
5.  [System Architecture](System_Architecture.md) - High-level system design
6.  [Data Interface](Data_Interface.md) - Process settings and LCA factors

## Quick Start

### 1. Run the LegHb Flowsheet (Config 1 - Default)
To execute the Leghemoglobin simulation with food-grade configuration:

```bash
python -c "from biorefineries.prefers.v1.LegHb.system import create_LegHb_system; sys = create_LegHb_system(); sys.simulate(); sys.show()"
```

Or run the config script directly:
```bash
python v1/LegHb/system/_config1.py
```

### 2. Run with Config 2 (Research-Grade with IX)
```bash
python v1/LegHb/system/_config2.py
```

### 3. Run the TEA (Techno-Economic Analysis)
To run the full economic analysis with design specifications:

```bash
python v1/LegHb/_tea.py
```

Optional arguments for config selection:
```bash
python v1/LegHb/_tea.py --config config2 --production 500
```

### 4. Run Uncertainty Analysis
To perform Monte Carlo simulations:

```bash
python v1/LegHb/analyses/uncertainty_and_sensitivity.py
```

## Configuration Comparison

| Feature | Config 1 (Food-Grade) | Config 2 (Research-Grade) |
|---------|----------------------|---------------------------|
| Primary Separation | Biomass harvest + cell wash | Direct cell disruption |
| Chromatography | None | Ion Exchange (ResinColumn2) |
| Polishing | HTST Pasteurization (72°C) | Nanofiltration DF |
| Final Formulation | Antioxidant + dilution | Direct concentration |
| Product Temperature | 4°C | 0°C |
| Debris Handling | Disposal | Boiler feed |

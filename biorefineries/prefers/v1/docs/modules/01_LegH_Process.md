# LegH Process Module

**Directory:** `v1/LegH/`

## Process Overview
This module simulates the production of **Leghemoglobin** using an aerobic fermentation of *Pichia pastoris*.

## Process Flow Summary
1.  **Preparation**: Glucose and Ammonia feed preparation.
2.  **Seed Train**: 5-stage expansion (`R301`).
3.  **Fermentation**: Fed-batch production (`R302`).
4.  **Harvest**: Centrifugation (`S401`) to separate biomass.
    *   *Note: In LegH system, downstream is currently simplified or uses standard library units.*

## Configuration & Compatibility Check

### System Overrides
The `LegH/system.py` script overrides several defaults from `units.py`.

**1. SeedTrain (R301)**
*   **Temperature:** Explicitly set to `32 + 273.15` K (Matches default).

**2. AeratedFermentation (R302)**
*   **Temperature:** Set to `T-operation` (30°C).
    *   *Default*: Undefined in base class, effectively controlled by cooling loop.
*   **Efficiency**: `compressor_isentropic_efficiency` = 0.85.

### Chemicals (`LegH/chemicals.py`)
*   **Key Components**: `Leghemoglobin`, `Pichia_pastoris`, `Glucose`, `H2SO4`, `NH3`.

### TEA Settings (`LegH/tea.py`)
Uses `PreFerSTEA` wrapper.
*   **Depreciation**: IRAS 6-year (Singapore).
*   **Tax**: 17% (Singapore Corporate Tax).

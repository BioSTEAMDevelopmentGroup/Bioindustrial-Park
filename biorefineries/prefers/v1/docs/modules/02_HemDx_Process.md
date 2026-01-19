# HemDx Process Module

**Directory:** `v1/HemDx/`

## Process Overview
Produces **Protein-Free Hemedextrin (Heme B)** using *Corynebacterium glutamicum*.

## Process Differences (vs LegHb)
This is NOT a clone. It contains unique biology and downstream logic.

1.  **Biology**:
    *   Organism: *Corynebacterium glutamicum* (vs *Pichia*).
    *   Product: Heme B (Intracellular & Extracellular).
2.  **Fermentation**:
    *   2-Stage logic: Growth Phase (72h) + Production Phase (90h).
    *   Operating Temp: 30°C.
    *   DO Setpoint: 45%.

## Downstream Specifics
**Custom Unit Logic in `system.py`:**

### S402 CellDisruption
*Uses `units.CellDisruption` but with custom biology mapping.*
*   **Target**: `Corynebacterium_glutamicum`.
*   **Component Fractions**:
    *   Protein: 0.45
    *   Cellulose: 0.22
    *   Xylan: 0.15
    *   OleicAcid: 0.08
    *   RNA: 0.10

### U401 Resin Adsorption
*Custom unit `units.ResinAdsorption` (Verified in `HemDx/system.py`)*
*   **Target Yield**: 98%.
*   **Impurity Removal**: 99.9%.

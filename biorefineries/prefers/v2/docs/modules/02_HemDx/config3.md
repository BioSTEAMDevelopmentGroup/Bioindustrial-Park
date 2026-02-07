# N-HemDx Full Production System - Config 3 (Extracellular)

## Overview
- **Product:** N-HemDx
- **Organism:** *Corynebacterium glutamicum*
- **Feedstock:** Glucose
- **Key Feature:** Extracellular production (High Secretion).
- **Secretion Factor (SF):** 0.9 (Baseline)

## Process Description
Config 3 represents a scenario where the product is efficiently secreted into the broth, minimizing the need for cell disruption.

### Key Modifications (vs Config 1)
- **Area 400 (Recovery):**
    - **Cell Disruption Train REMOVED:** Units S402 (Homogenizer), H401, S403 (Debris Centrifuge), S405 (Lysate Filter), and M404 are removed.
    - **Stream Routing:**
        - **WashedCellCream (C402-0):** Sent directly to S406 (Screw Press) for dewatering.
        - **SupernatantCake (S404-0):** Sent to S406.
    - **M405:** Receives only Filtered Supernatant (S404-1).

### Key Parameters
- **Secretion Fraction:** 0.7 - 1.0 (Uniform Distribution)
- **Titer Control:** Elasticity-based yield adjustment solved with `flexsolve.IQ_interpolation` plus one-step correction; model titer bounds tightened to ±1% for $<0.01$ g/L deviation.

## Unit Operations
See `config1.md` for shared units.
- **S401:** Centrifuge
- **S404:** Supernatant Microfiltration (Retained for extracellular product recovery)
- **S406:** Screw Press (Handles all solids: washed cells + filter cake)

## Related Files
- System Config: `biorefineries/prefers/v2/HemDx/system/_config3.py`
- TEA Config: `biorefineries/prefers/v2/HemDx/_tea_config3.py`

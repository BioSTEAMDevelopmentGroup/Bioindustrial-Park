# N-HemDx Full Production System - Config 2 (Intracellular)

## Overview
- **Product:** N-HemDx
- **Organism:** *Corynebacterium glutamicum*
- **Feedstock:** Glucose
- **Key Feature:** Intracellular production (Low Secretion).
- **Secretion Factor (SF):** 0.1 (Baseline)

## Process Description
Config 2 represents a scenario where the product is primarily retained intracellularly or secretion is very low.

### Key Modifications (vs Config 1)
- **Area 400 (Recovery):**
    - **S404 Removed:** Supernatant microfiltration is removed as most product is in the cells.
    - **M504 Added:** Supernatant from S401 is routed directly to Wastewater Treatment (via M504).
    - **M405:** Receives only Clarified Lysate (S405-1).

### Key Parameters
- **Secretion Fraction:** 0.0 - 0.3 (Uniform Distribution)

## Unit Operations
See `config1.md` for shared units.
- **S401:** Centrifuge (Separates CellCream and Supernatant)
- **M504:** Supernatant Mixer -> WWT
- **S402:** Cell Disruption (Standard)
- **S403/S405:** Lysate Clarification

## Related Files
- System Config: `biorefineries/prefers/v1/HemDx/system/_config2.py`
- TEA Config: `biorefineries/prefers/v1/HemDx/_tea_config2.py`

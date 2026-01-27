# Centrifuge Unit Specification

**Class:** `Centrifuge`
**Base Class:** `bst.SolidsCentrifuge`
**Source:** `v1/_units.py`

## Purpose
Extended centrifuge model suitable for pilot-scale or low-solids operations, preventing warnings where standard BioSTEAM models impose strict lower bounds.

## Modifications

1. **Solids Loading Range**:
   - **Extended Lower Bound**: 0.1 ton/hr (vs. standard 1-2 ton/hr).
   - Allows simulation of pilot production rates without `lb_warning`.

2. **Cost Extrapolation**:
   - Uses `min_solids_loading_cost` (0.1 ton/hr) as a floor for costing calculations.
   - Maintains standard scaling exponents:
     - Reciprocating Pusher: $n=0.5$
     - Scroll Solid Bowl: $n=0.3$

## References
- Humbird et al. (2011) - NREL Biofuels Baseline.
- Seider et al. (2017) - Product and Process Design Principles.

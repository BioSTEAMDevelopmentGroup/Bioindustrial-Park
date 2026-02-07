# NeutralizationTank1 Unit Specification

**Class:** `NeutralizationTank1`
**Type:** CSTR / MixTank with Reaction and Cooling
**Source:** `v1/_units.py`

## Purpose
Performs pH adjustment (neutralization) of acidic/basic streams, specifically handling heat of neutralization. used in Area 500/900 for waste treatment.

## Design Basis

| Parameter          | Default | Unit  | Description               |
| :----------------- | :------ | :---- | :------------------------ |
| **Residence Time** | 2.0     | hr    | Tank hold time            |
| **Agitation**      | 0.5     | kW/m³ | Mixing power density      |
| **Temperature**    | 25      | °C    | Target outlet temperature |

## Sub-Units
1. **Reactor**:
   - Handles H2SO4 + NaOH -> Na2SO4 reactions.
   - Calculates adiabatic heat release.
2. **Cooler (`HXutility`)**:
   - Removes heat of neutralization to maintain setpoint `T`.
3. **Agitator**:
   - Sized based on liquid volume.

## Cost Model

Includes three components:
1. **Tank**: $96,000 \cdot (V/10)^{0.7}$
2. **Agitator**: $31,800 \cdot (P/11.3)^{0.5}$
3. **Cooler**: $50,000 \cdot (Q/10^6)^{0.6}$ (Simplified duty-based cost)

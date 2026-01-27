# AeratedFermentation Unit Specification

**Class:** `AeratedFermentation`
**Base Class:** `bst.AeratedBioreactor`
**Source:** `v1/_units.py`

## Purpose
Main production fermenter for LegHb. Optimized for aerobic yeast fermentation (*K. marxianus* / *P. pastoris*).

## Design Specification

| Parameter          | Default | Unit     | Description                                    |
| :----------------- | :------ | :------- | :--------------------------------------------- |
| **V_max_default**  | 500     | m³       | Maximum vessel working volume                  |
| **Agitation**      | 0.985   | kW/m³    | Power input for oxygen transfer (Stirred Tank) |
| **Cooler dP**      | 20,684  | Pa       | Pressure drop across cooling loop              |
| **Compressor Eff** | 0.85    | -        | Isentropic efficiency for air supply           |
| **Q_O2**           | -110    | kcal/mol | Heat of reaction per mole O2 consumed          |

## Logic Enhancements

1. **Reaction Safeguards**:
   - Prevents negative flows.
   - Checks glucose availability (< 1e-6 threshold).
   - Forces water constraint (non-negative).
2. **Vent Management**:
   - Separates non-condensables (CO2, O2, N2) to vent stream.
   - Enforces correct phase behavior.

## Usage
Used in Area 300 (Conversion) with `create_fermentation_reactions` supplying the reaction kinetics.

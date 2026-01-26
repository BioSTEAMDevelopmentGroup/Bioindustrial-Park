# Walkthrough - ResinColumn Mass Balance Fix

I have successfully refactored the `ResinColumn` unit operation to strictly enforce mass balance in **Adsorption (Flow-Through) Mode**.

## 1. Problem
The previous `ResinColumn` implementation in `_units.py` had a logic error in `_run_adsorption` where auxiliary buffer streams (Wash, Elution, Regeneration) were accounted for in the inputs but ignored in the outputs, leading to a massive mass balance violation.

## 2. Changes Implemented
- **Refactored `_run_adsorption`**:
    - Modeled the unit as a **Negative Chromatography** step (Flow-Through for Product, Bind for Impurities).
    - **Stream Routing**:
        - `ins[0]` (Feed) &rarr; `outs[0]` (Product Flow-Through) + `outs[2]` (Adsorbed Impurities).
        - `ins[1]` (Wash Buffer) &rarr; `outs[3]` (Spent Wash).
        - `ins[2]` (Elute Buffer) + `ins[3]` (Regen Buffer) &rarr; `outs[2]` (Combined Chemical Waste).
    - **Mass Balance**: Enforced strict conservation of mass for all components.
- **Updated Design Parameters**:
    - Added literature citations for `EBCT` (5-30 min), `Superficial Velocity` (5-20 m/h), and `Adsorbent Density` to the unit's parameters.

## 3. Verification Results
Executed `_config1.py` and inspected `U501`:

| Stream  | Flow (kmol/hr)     | Content                          | Status       |
| :------ | :----------------- | :------------------------------- | :----------- |
| **In**  | **CombinedLysate** | ~276 H2O + Solutes               | Feed         |
| **In**  | **Wash Buffer**    | ~829 H2O + NaCl                  | Aux          |
| **In**  | **Elute + Regen**  | ~18 H2O + Chemicals              | Aux          |
| **Out** | **ResinTreated**   | ~276 H2O + Heme                  | **Balanced** |
| **Out** | **ResinWash**      | ~829 H2O + NaCl                  | **Balanced** |
| **Out** | **ResinWaste**     | ~18 H2O + Chemicals + Impurities | **Balanced** |

The unit now executes without error and correctly separates impurities (Protein, DNA) from the Heme product while accounting for all buffer flows.

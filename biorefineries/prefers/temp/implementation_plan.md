# ResinColumn Upgrade & Mass Balance Fix Plan

## Goal Description
Fix the mass balance violation in `ResinColumn` (U501) where Wash, Elution, and Regeneration buffers vanish during simulation. Align the unit's logic with a **Flow-Through (Negative) Adsorption** strategy, where impurities bind to the resin and the target Heme product flows through. Additionally, validate design parameters against literature where possible.

## User Review Required
> [!IMPORTANT]
> **Adsorption Mode Logic Change**: The unit will be explicitly modeled as **Negative Chromatography**:
> - **ResinTreated (Out 0)**: Contains the Target Product (Heme) and Flow-Through components.
> - **ResinWash (Out 3)**: Contains the Wash Buffer (NaCl).
> - **ResinWaste (Out 2)**: Contains the Regeneration Buffer (Ethanol/NaOH) + **Desorbed Impurities**.
> - **ResinFlowthrough (Out 1)**: Will be empty (unused) in this configuration.

## Proposed Changes

### [biorefineries/prefers/v1/_units.py]
#### [MODIFY] [ResinColumn](file:///c:/Programming/PREFERS/Bioindustrial-Park/biorefineries/prefers/v1/_units.py)
- **Refactor `_run_adsorption`**:
    - Assign `outs[3]` (Wash) -> Copy `ins[1]` (Wash Buffer).
    - Assign `outs[2]` (Waste) -> Copy `ins[2]` (Elution is likely not used in negative mode, or combined with Regen) + `ins[3]` (Regen).
        - *Correction*: The config has 3 buffers: Wash (NaCl), Elute (NaOH), Regen (Ethanol).
        - In Negative Mode, "Elution" step might be the "Strip Impurities" step. "Regen" step might be "Sanitize/Re-equilibrate".
        - I will sum all "Cleaning" buffers (Elute + Regen) into `ResinWaste` (Out 2) along with the Desorbed Impurities.
    - Ensure `ResinTreated` (Out 0) captures only the flow-through mass (Feed - Impurities).
- **Update Parameters**:
    - Add literature citations for `EBCT` (Empty Bed Contact Time), `velocity`, and `loading` if available from recent search.

### [biorefineries/prefers/v1/HemDx/system/_config1.py]
#### [MODIFY] [_config1.py](file:///c:/Programming/PREFERS/Bioindustrial-Park/biorefineries/prefers/v1/HemDx/system/_config1.py)
- Verify stream connections match the new logic.
- (Optional) Adjust `TargetProduct_Yield` (Removal Efficiency) if literature suggests different values for protein removal on hydrophobic resins.

## Verification Plan

### Automated Tests
1. **Mass Balance Check**:
    - Run `_config1.py`.
    - Retrieve `U501`.
    - Assert `abs(U501.mass_in - U501.mass_out) < 1e-5` (Strict Mass Balance).
2. **Logic Check**:
    - Verify `ResinTreated` contains high Heme and low Protein.
    - Verify `ResinWaste` contains the Protein + Ethanol/NaOH.
    - Verify `ResinWash` contains NaCl/Water.

### Manual Verification
- Execute `_config1.py` and inspect the `U501.show()` output at the end of the script.

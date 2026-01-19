# LegHemoglobin Process Flow (`LegHb`)

**Source:** `v1/LegHb/system/_config1.py`  
**System ID:** `LegHb_config1`

This module defines the complete production process for LegHemoglobin, including upstream fermentation, downstream purification, and facility support systems.

## 1. Upstream Process (Fermentation)

The upstream section focuses on preparing the inoculum and running the main fermentation to produce LegHemoglobin and Heme B.

### Sequence
*   **M301 (MixTank):** Preparation of Seed Solution 1.
*   **M302 (MixTank):** Preparation of Seed Solution 2 and Culture medium.
*   **M303 (SeedHoldTank):** Holding tank combining seed streams.
*   **R301 (SeedTrain):** 5-stage seed train for inoculum expansion.
    *   *Input:* Seed mix.
    *   *Output:* Inoculum (to R302).
*   **M304 (MixTank):** Glucose feed preparation.
*   **T301 (StorageTank):** Glucose feed holding tank.
*   **T302 (AmmoniaStorageTank):** Ammonia feed storage.
*   **R302 (AeratedFermentation):** Main production fermenter (Fed-batch).
    *   *Inputs:* Inoculum, Glucose feed, Ammonia, Air.
    *   *Outputs:* Fermentation Broth, Vent gases.

### Control Specifications
*   **R302 Specification (`update_reaction_time_and_yield`):**
    *   Calculates residence time (`tau`): `target_titer / target_productivity`.
    *   Updates reaction product yield based on `target_yield`.

## 2. Downstream Process (Purification)

The downstream section isolates and purifies Leghemoglobin through cell disruption, filtration, and chromatography.

### Sequence
*   **S401 (CellDisruption):** High-pressure homogenizer to lyse cells and release product.
*   **S402 (Centrifuge):** Separates cell debris (solid) from the liquid supernatant containing product.
    *   *Deposit:* To ScrewPress (S403).
    *   *Supernatant:* To Purification (M401/U401).
*   **S403 (ScrewPress):** Dewaters cell mass waste.
*   **S404 (Splitter):** Residual cell mass removal (Mass balance correction).
*   **M401 (MixTank) & H401 (HXutility):** Preparation and cooling of Ultrafiltration Buffer.
*   **U401 (Diafiltration):** Primary ultrafiltration for buffer exchange and initial concentration.
    *   *Target:* Retain LegHemoglobin.
*   **M402, M403, M404 (MixTanks):** Preparation of IX buffers (Equilibration, Elution, Regeneration).
*   **H402, H403, H404 (HXutility):** Cooling of IX buffers.
*   **U402 (IonExchangeCycle):** Chromatography capture step.
    *   *Inputs:* U401 Retentate + Buffers.
    *   *Outputs:* Product Eluate, Waste streams.
*   **M405 (MixTank) & H405 (HXutility):** Preparation and cooling of Nanofiltration Buffer.
*   **U403 (Diafiltration):** Secondary filtration (Polishing/Concentration).
*   **U404 (Ultrafiltration):** Final concentration step.
*   **H406 (HXutility):** Final product cooling to 0°C.

### Control Specifications
*   **U404 Specification (`U404_adjust_water_recovery`):**
    *   **Goal:** Achieve target total solids content (sets to 12.0%).
    *   **Variable:** `FeedWater_Recovery_to_Permeate` (Water removal rate).
    *   **Method:** Uses `flexsolve.IQ_interpolation` to optimize water recovery dynamically based on feed composition.

## 3. Waste Treatment & Facilities

*   **T501 (SulfuricAcidStorageTank) & M502 (NeutralizationTank1):** Neutralization of base waste from IX regeneration.
*   **S501, S503 (ReverseOsmosis):** Water recovery from waste streams.
*   **Facilities:**
    *   `CT` (Cooling Tower)
    *   `CWP` (Chilled Water Package)
    *   `BT` (Boiler Turbogenerator): Burns sludge/biomass for energy.
    *   `PWC` (Process Water Center): Manages water recycling and makeup.

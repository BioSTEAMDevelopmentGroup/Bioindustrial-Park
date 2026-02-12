# LegHemoglobin Process Flow (`LegHb`) - Config 1

**Source:** `biorefineries/prefers/v2/LegHb/system/_config1.py`  
**System ID:** `LegHb_sys`  
**Architecture:** Modular Process Areas

This configuration implements an **intracellular Leghemoglobin recovery** process optimized for minimal purification steps while meeting product specifications. The design emphasizes biomass harvest, high-pressure homogenization, and UF/DF purification with HTST pasteurization.

## Module Architecture

The system is organized into discrete Process Area functions following standard BioSTEAM conventions:

```
create_LegHb_system()
├── create_area_200_media_prep()      # Media Preparation
├── create_area_300_conversion()      # Fermentation
├── create_area_400_recovery()        # Harvest & Cell Disruption
├── create_area_500_purification()    # UF/DF Membrane Separation
├── create_area_600_formulation()     # Thermal Treatment & Product
└── create_area_900_facilities()      # Utilities & Wastewater
```

### Shared Configuration Functions

| Function                          | Purpose                                       |
| --------------------------------- | --------------------------------------------- |
| `get_fermentation_parameters()`   | Returns validated fermentation parameters     |
| `create_fermentation_reactions()` | Creates reaction systems for LegHb production |

---

## Area 200: Media Preparation

**Function:** `create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)`

Purpose: Prepare all feed solutions for fermentation.

### Unit Operations

| Unit ID | Type               | Description                                              |
| ------- | ------------------ | -------------------------------------------------------- |
| M201    | MixTank            | Seed Solution 1 preparation (τ = 16 hr)                  |
| M202    | MixTank            | Seed Solution 2 & Culture medium preparation (τ = 16 hr) |
| M203    | SeedHoldTank       | Combines seed streams                                    |
| M204    | MixTank            | Glucose solution preparation (50% dilution)              |
| T201    | StorageTank        | Glucose storage tank (tau = 16*4 + 72 hr)                |
| T202    | AmmoniaStorageTank | Ammonia storage tank (25 wt% aqueous ammonia)            |
| S202    | Splitter           | Splits Ammonia between seed train and main fermenter     |

### Outputs

| Key           | Stream  | Description                         |
| ------------- | ------- | ----------------------------------- |
| `seed_out`    | M203Out | Combined seed solution → Area 300   |
| `glucose_out` | T201Out | Diluted glucose solution → Area 300 |
| `ammonia_out` | T202Out | Ammonia solution → Area 300         |

---

## Area 300: Conversion (Fermentation)

**Function:** `create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, reactions, params)`

Purpose: Convert glucose to LegHb through microbial fermentation.

### Unit Operations

| Unit ID | Type                | Description                                          |
| ------- | ------------------- | ---------------------------------------------------- |
| R301    | SeedTrain           | 5-stage seed train for inoculum expansion (T = 32°C) |
| R302    | AeratedFermentation | Main production fermenter (Fed-batch)                |

### Fermentation Parameters

| Parameter             | Value         | Units        |
| --------------------- | ------------- | ------------ |
| `titer_LegHb`         | 5.0           | g/L          |
| `productivity_LegHb`  | 5/72          | g/L/hr       |
| `yield_LegHb (Y_p)`   | ~0.0308       | wt/wt        |
| `Y_b (biomass yield)` | 0.53          | wt/wt        |
| `theta_O2`            | 0.5           | % saturation |
| `agitation_power`     | 0.985         | kW/m³        |
| `T_operation`         | 305.15 (32°C) | K            |
| `V_max`               | 500           | m³           |

### Control Specifications (Internalized in R302)

v2 internalizes titer convergence and NH3 optimization as a single `@R302.add_specification(run=False)` — no external `convergence.py`, `optimize_NH3_loading`, or `adjust_glucose_for_titer` functions needed.

- **Combined Specification (`update_fermentation_and_nh3`):**
  - Empirically measures NH3 demand by running R302 with 20× excess, then sets exact NH3 (replaces `optimize_NH3_loading`).
  - Iteratively adjusts fermentation yield to hit target titer using elasticity-based correction ($Elasticity = 1.3$, $Tolerance = 0.1\%$) with adaptive damping (replaces `adjust_glucose_for_titer`).
  - Updates residence time (`tau`): `target_titer / target_productivity`.
  - Clamps yield to $[0.005, 0.15]$ to preserve physical feasibility.
  - Updates upstream `NH3_25wt` source and `S202` split ratio for consistency.

### Outputs

| Key         | Stream | Description                |
| ----------- | ------ | -------------------------- |
| `broth_out` | Broth  | Cell suspension → Area 400 |

---

## Area 400: Recovery (Biomass Harvest & Cell Disruption)

**Function:** `create_area_400_recovery(broth_in, DfUltraBuffer2)`

Purpose: Release and clarify intracellular LegHb from yeast cells.

### Step 1: Biomass Harvest (Primary Separation)

| Unit ID | Type             | Description                                                 |
| ------- | ---------------- | ----------------------------------------------------------- |
| C401    | SolidsCentrifuge | Primary harvest centrifuge (moisture_content = 0.40)        |
| M401    | MixTank          | Wash buffer preparation (DfUltraBuffer2 + Water4)           |
| H402    | HXutility        | Cool wash buffer to 10°C                                    |
| M402    | MixTank          | Cell wash mixer                                             |
| C402    | SolidsCentrifuge | Washed cell separation (moisture_content = 0.55)            |

- **Cell capture:** From `c.chemical_groups['SolidsCentrifuge']`

### Step 2: Cell Disruption (High-Pressure Homogenization)

| Unit ID | Type           | Description                                    |
| ------- | -------------- | ---------------------------------------------- |
| S401    | CellDisruption | High-pressure homogenizer at 1000 bar          |
| H401    | HXutility      | Immediate post-valve cooling to <15°C          |

- **Cell disruption efficiency:** 0.87
- **Operating pressure:** 1000 bar (100 MPa)

### Step 3: Lysate Clarification (Debris Removal)

| Unit ID | Type                   | Description                                  |
| ------- | ---------------------- | -------------------------------------------- |
| C403    | SolidsCentrifuge       | Debris removal (moisture_content = 0.20)     |
| S403    | FiltrationAdv (MF)     | Depth filtration (solid_capture_eff = 0.85)  |
| S404    | ScrewPress             | Debris dewatering (split = 0.99)             |

- **Solid capture efficiency:** 0.85 (FiltrationAdv)

### Outputs

| Key                 | Stream           | Description                           |
| ------------------- | ---------------- | ------------------------------------- |
| `clarified_lysate`  | ClarifiedLysate  | Cell-free protein solution → Area 500 |
| `spent_media`       | SpentMedia       | To wastewater (Area 900)              |
| `wash_effluent`     | WashEffluent     | To wastewater (Area 900)              |
| `dehydrated_debris` | DehydratedDebris | To disposal (Area 900)                |
| `press_liquor`      | PressLiquor      | To wastewater (Area 900)              |

---

## Area 500: Purification (UF/DF Membrane Separation)

**Function:** `create_area_500_purification(clarified_lysate, DfUltraBuffer1)`

Purpose: Concentrate LegHb and remove impurities using tangential flow filtration.

### Unit Operations

| Unit ID | Type                      | Description                                     |
| ------- | ------------------------- | ----------------------------------------------- |
| M501    | MixTank                   | DF buffer preparation (DfUltraBuffer1 + Water5) |
| H501    | HXutility                 | Cool DF buffer to 5°C                           |
| U501    | DiafiltrationAdv (UF)      | TFF with 3-10 kDa MWCO membrane                 |
| U502    | DiafiltrationAdv           | Final concentration step                        |

### Process Parameters

| Parameter      | Value                   |
| -------------- | ----------------------- |
| LegH retention | 95%                     |
| Salt retention | 5% (salts wash through) |
| Diavolumes     | 6 (buffer water = 6× feed water) |

### Control Specifications

- **U502 Specification (`U502_adjust_water_recovery`):**
  - Fixed `Salt_Retention = 0.05` and run with configured water recovery

### Outputs

| Key                      | Stream           | Description               |
| ------------------------ | ---------------- | ------------------------- |
| `concentrated_product`   | ConcentratedLegH | Purified LegHb → Area 600 |
| `uf_permeate`            | U501.outs[1]     | To wastewater (Area 900)  |
| `concentration_permeate` | U502.outs[1]     | To wastewater (Area 900)  |

---

## Area 600: Formulation (Thermal Treatment & Final Product)

**Function:** `create_area_600_formulation(concentrated_product, LegHb_3)`

Purpose: Ensure microbiological safety and achieve final product specification.

### Step 1: Thermal Stabilization (HTST Pasteurization)

| Unit ID | Type      | Description                  |
| ------- | --------- | ---------------------------- |
| H603    | HXutility | Heat to 72°C (70-75°C range) |
| T601    | MixTank   | Hold time tank (~30 seconds) |

- Inactivates pathogens
- Precipitates unstable host proteins
- LegHb (stable heme-protein) remains soluble

### Step 2: Final Formulation & Cooling

| Unit ID | Type      | Description                                           |
| ------- | --------- | ----------------------------------------------------- |
| M604    | MixTank   | Formulation mixer with antioxidant and dilution water |
| H606    | HXutility | Final rapid cooling to 4°C                            |

### Formulation Targets

| Specification    | Target   | Range   |
| ---------------- | -------- | ------- |
| Total solids     | 20%      | max 24% |
| LegHb content    | 7.5%     | 6-9%    |
| Sodium ascorbate | 0.1% w/w | -       |

### Control Specifications

- **M604 Specification (`update_formulation`):**
  - Calculates dilution water from total solids and LegHb targets.
  - Sets sodium ascorbate at 0.1% w/w of final product (with 10% w/w solution water).

### Outputs

| Key       | Stream  | Description          |
| --------- | ------- | -------------------- |
| `product` | LegHb_3 | Final product stream |

---

## Area 900: Facilities & Utilities

**Function:** `create_area_900_facilities(area_400, area_500, effluent1, use_area_convention)`

Purpose: Provide utilities and manage wastewater.

### Wastewater Treatment

Wastewater is routed to BioSTEAM’s wastewater treatment system rather than direct RO. This avoids combining high-phosphate and large-molecule streams before biological treatment.

| Unit/System ID           | Type                               | Description                              |
| ------------------------ | ---------------------------------- | ---------------------------------------- |
| wastewater_treatment_sys | create_wastewater_treatment_system | Centralized WWT (biological + RO)        |
| (commented) M901, S901   | MixTank + ReverseOsmosis           | Legacy direct-RO path (kept as comments) |

### Utility Systems

| Unit ID | Type                 | Description                                    |
| ------- | -------------------- | ---------------------------------------------- |
| CT      | CoolingTower         | Cooling water system                           |
| CWP     | ChilledWaterPackage  | Chilled water for process cooling              |
| M902    | Mixer                | Combines dehydrated debris with WWT sludge     |
| BT      | BoilerTurbogenerator | Biomass + biogas combustion and power          |
| PWC     | ProcessWaterCenter   | Water recycling and makeup management (WWT RO) |

**Note:** Debris stream routes to the boiler via the prefers-local `BoilerTurbogenerator` subclass with guarded emissions enthalpy updates and glucose-based thermo fallbacks for large biomolecules.

---

## Exported Functions

### System Creation

```python
create_LegHb_system(ins, outs, use_area_convention=False)
```
Factory function creating the complete production system using modular Process Area functions.

### Area Functions

| Function                         | Description                   |
| -------------------------------- | ----------------------------- |
| `create_area_200_media_prep()`   | Create Media Preparation area |
| `create_area_300_conversion()`   | Create Fermentation area      |
| `create_area_400_recovery()`     | Create Recovery area          |
| `create_area_500_purification()` | Create Purification area      |
| `create_area_600_formulation()`  | Create Formulation area       |
| `create_area_900_facilities()`   | Create Facilities area        |

### Design Specifications

```python
set_production_rate(system, target_production_rate_kg_hr, verbose=True)
```
Adjusts system inputs to achieve target LegHb_3 production rate using global scaling factor with flexsolve IQ_interpolation. NH3 optimization is handled automatically by R302's internalized specification during simulation.

```python
check_LegHb_specifications(product_stream)
```
Verifies LegHb_3 product stream meets composition specifications:
- Fat (OleicAcid): 0-2%
- Carbohydrates: 0-4%
- Leghemoglobin: 6-9%
- Total Solids: 0-24%
- Protein Purity: ≥65%

### Heme Equivalent Reporting
For comparison with Heme-rich products (e.g. HemDx), the system reports:
- **Heme Equivalent Yield**: Mass flow of Heme B prosthetic group contained within the Leghemoglobin product.
  - Calculation: `LegHb_moles / 763 * MW_Heme_B` (accounting for C1-normalized LegHb definition)


---

## Change History

| Date       | Version | Description                                                            |
| ---------- | ------- | ---------------------------------------------------------------------- |
| 2026-02-07 | 2.4     | Internalized titer + NH3 spec in R302 with empirical NH3 measurement   |
| 2026-01-23 | 2.2     | Replaced direct RO with BioSTEAM wastewater treatment system           |
| 2026-01-21 | 2.0     | Modular refactoring into Process Areas                                 |

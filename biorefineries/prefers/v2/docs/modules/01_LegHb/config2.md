# LegHemoglobin Process Flow (`LegHb`) - Config 2

**Source:** `v1/LegHb/system/_config2.py`  
**System ID:** `LegHb_sys`

This configuration implements a **traditional downstream purification** process with ion exchange chromatography. The design includes intracellular recovery (aligned with Config 1), ultrafiltration/diafiltration, ion exchange chromatography, nanofiltration diafiltration, and final ultrafiltration concentration.

## Module Architecture

The system is organized into discrete Process Area functions following the Config 1 coding style:

```
create_LegHb_system()
├── create_area_200_media_prep()      # Media Preparation
├── create_area_300_conversion()      # Fermentation
├── create_area_400_recovery()        # Harvest & Cell Disruption
├── create_area_500_purification()    # UF/DF + IEX + NF + UF
├── create_area_600_formulation()     # Final cooling
└── create_area_900_facilities()      # Utilities & Wastewater
```

## 1. Upstream Process (Fermentation)

Identical to Config 1 - prepares inoculum and runs main fermentation.

### Unit Operations

| Unit ID | Type                | Description                                              |
| ------- | ------------------- | -------------------------------------------------------- |
| M301    | MixTank             | Seed Solution 1 preparation (τ = 16 hr)                  |
| M302    | MixTank             | Seed Solution 2 & Culture medium preparation (τ = 16 hr) |
| M303    | SeedHoldTank        | Combines seed streams                                    |
| R301    | SeedTrain           | 5-stage seed train for inoculum expansion (T = 32°C)     |
| M304    | MixTank             | Glucose feed preparation (τ = 16 hr)                     |
| T301    | StorageTank         | Glucose feed holding (τ = 136 hr)                        |
| T302    | AmmoniaStorageTank  | Ammonia feed storage                                     |
| R302    | AeratedFermentation | Main production fermenter (Fed-batch)                    |

### Fermentation Parameters
Same as Config 1:

| Parameter             | Value         | Units  |
| --------------------- | ------------- | ------ |
| `titer_LegHb`         | 5.0           | g/L    |
| `productivity_LegHb`  | 0.069         | g/L/hr |
| `yield_LegHb (Y_p)`   | 0.0333        | wt/wt  |
| `Y_b (biomass yield)` | 0.53          | wt/wt  |
| `T_operation`         | 305.15 (32°C) | K      |
| `V_max`               | 500           | m³     |

### Control Specifications
- **R302 Specification (`update_reaction_time_and_yield`):**
  - Calculates residence time (`tau`): `target_titer / target_productivity`
  - Updates reaction product yield based on `target_yield`

## 2. Downstream Process (Traditional Purification)

### Cell Disruption & Solid-Liquid Separation
Aligned with Config 1 recovery train before downstream purification.

| Unit ID | Type                   | Description                                     |
| ------- | ---------------------- | ----------------------------------------------- |
| C401    | Centrifuge             | Primary harvest                                 |
| M401/H402/M402/C402 | MixTank + HX + Centrifuge | Wash buffer prep and cell washing         |
| S401    | CellDisruption         | High-pressure homogenizer                       |
| H401    | HXutility              | Post-valve cooling                              |
| S402    | SolidsCentrifuge       | Debris removal                                  |
| S403    | Filtration (MF preset) | Lysate clarification                            |
| S404    | ScrewPress             | Debris dewatering                               |

### Primary Ultrafiltration/Diafiltration
| Unit ID | Type                      | Description                               |
| ------- | ------------------------- | ----------------------------------------- |
| M401    | MixTank                   | DfUltraBuffer preparation (4× feed water) |
| H401    | HXutility                 | Cool to 5°C                               |
| U401    | Diafiltration (UF preset) | Primary buffer exchange and concentration |

**U401 Parameters:**
- Target Product: Leghemoglobin
- Salt retention: Low (salts wash through)

### Ion Exchange Chromatography
| Unit ID   | Type         | Description                                    |
| --------- | ------------ | ---------------------------------------------- |
| M402      | MixTank      | IXEquilibriumBuffer preparation                |
| M403      | MixTank      | IXElutionBuffer preparation                    |
| M404      | MixTank      | IXRegenerationSolution preparation             |
| H402-H404 | HXutility    | Cool buffers to 5°C                            |
| U402      | ResinColumn2 | Ion exchange capture step (IonExchange preset) |

**U402 Outputs:**
1. Product eluate (U402-0)
2. Flowthrough waste (U402-1)
3. Wash waste (U402-2)
4. Regeneration waste (U402-3)

**Buffer Scaling Specifications:**
- Equilibration: `wash_CV × feed water volume`
- Elution: `elution_CV × feed water volume`
- Regeneration: `regeneration_CV × feed water volume`

### Secondary Diafiltration (Nanofiltration)
| Unit ID | Type                      | Description                                   |
| ------- | ------------------------- | --------------------------------------------- |
| M405    | MixTank                   | DfNanoBuffer preparation (2× IX eluate water) |
| H405    | HXutility                 | Cool to 5°C                                   |
| U403    | Diafiltration (NF preset) | Polishing/concentration step                  |

### Final Concentration
| Unit ID | Type            | Description                 |
| ------- | --------------- | --------------------------- |
| U404    | Ultrafiltration | Final product concentration |

**U404 Control Specification (`U404_adjust_water_recovery`):**
- **Target:** LegHb content = 7.5% (range 6-9%)
- **Method:** Uses `flexsolve.IQ_interpolation` to optimize water recovery
- **Salt_Retention:** Dynamically set to 0.05 to wash out salts
- **Search range:** 50-99% water removal

## 3. Waste Treatment & Neutralization

### Waste Collection
| Unit/System ID           | Type                     | Description                                                     |
| ------------------------ | ------------------------ | --------------------------------------------------------------- |
| T501                     | SulfuricAcidStorageTank  | Acid for IX regeneration neutralization                         |
| M502                     | NeutralizationTank1      | Neutralizes IX regeneration waste (T = 20°C)                    |
| M904                     | Mixer                    | Routes neutralized regeneration waste to `effluent2`            |
| M905                     | MixTank                  | Combines WWT feeds and neutralizes H2SO4 in-stream              |
| S904                     | Splitter                 | Sends any residual H2SO4 to disposal (protects WWT)             |
| wastewater_treatment_sys | create_wastewater_treatment_system | Centralized WWT (biological + RO)                     |
| S902/S903                | Splitter                 | Splits sludge/biogas to boiler vs. system outlets               |
| (commented) M501, S501-S503 | MixTank + RO/Splitter | Legacy direct-RO path (kept as comments)                         |

**Acid Requirement:** Stoichiometric 1:2 ratio with NaOH plus 0.1% excess

### Water Recovery
Water recovery is handled within the wastewater treatment system (biological treatment + RO), with treated water routed to the Process Water Center.

### Final Product Cooling
| Unit ID | Type      | Description         |
| ------- | --------- | ------------------- |
| H406    | HXutility | Cool product to 0°C |

## 4. Facilities

| Unit ID | Type                 | Description                           |
| ------- | -------------------- | ------------------------------------- |
| CT      | CoolingTower         | Cooling water system                  |
| CWP     | ChilledWaterPackage  | Chilled water for process cooling     |
| BT      | BoilerTurbogenerator | Steam & power from cell mass waste    |
| PWC     | ProcessWaterCenter   | Water recycling and makeup management |

**Makeup Water Streams:** cooling tower, M301-M304 waters, buffer waters (Water1-8), boiler makeup

**Process Water Streams:** WWT treated water plus makeup streams

**BT Note:** Uses prefers-local `BoilerTurbogenerator` with guarded emissions enthalpy update and glucose-based thermo fallbacks for large biomolecules.

## 5. Exported Functions

### `create_LegHb_system(ins, outs, use_area_convention=False)`
Factory function creating the complete production system.

### `set_production_rate(system, target_production_rate_kg_hr)`
Adjusts system inputs to achieve target LegHb_3 production rate using global scaling factor.
- Search bounds: 0.1× to 5.0× baseline
- Uses flexsolve IQ_interpolation

### `check_LegHb_specifications(product_stream)`
Verifies LegHb_3 product stream meets composition specifications:
- Fat (OleicAcid): 0-2%
- Carbohydrates: 0-4%
- Leghemoglobin: 6-9%
- Total Solids: 0-24%
- Protein Purity: ≥65%

## 6. Key Differences from Config 1

| Aspect             | Config 1                      | Config 2                    |
| ------------------ | ----------------------------- | --------------------------- |
| Primary Separation | Biomass harvest centrifuges   | Direct cell disruption      |
| Chromatography     | None (simplified process)     | Ion Exchange (ResinColumn2) |
| Polishing          | HTST Pasteurization           | Nanofiltration DF           |
| Final Formulation  | Antioxidant + dilution mixing | Direct concentration        |
| Product Temp       | 4°C                           | 0°C                         |
| Waste Treatment    | Debris to disposal            | Cell mass to boiler         |

## Change History

| Date       | Version | Description                                           |
| ---------- | ------- | ----------------------------------------------------- |
| 2026-01-23 | 1.3     | Refactored to modular area functions and aligned recovery with Config 1 |
| 2026-01-23 | 1.2     | Replaced direct RO with BioSTEAM wastewater treatment system |
| 2026-01-23 | 1.1     | Verified script runs successfully with BT/CWP updates |
| 2025-06-04 | 1.0     | Initial implementation                                |


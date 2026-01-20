# LegHemoglobin Process Flow (`LegHb`) - Config 2

**Source:** `v1/LegHb/system/_config2.py`  
**System ID:** `LegHb_sys`

This configuration implements a **traditional downstream purification** process with ion exchange chromatography. The design includes cell disruption, microfiltration, ultrafiltration/diafiltration, ion exchange chromatography, nanofiltration diafiltration, and final ultrafiltration concentration.

## 1. Upstream Process (Fermentation)

Identical to Config 1 - prepares inoculum and runs main fermentation.

### Unit Operations

| Unit ID | Type | Description |
|---------|------|-------------|
| M301 | MixTank | Seed Solution 1 preparation (τ = 16 hr) |
| M302 | MixTank | Seed Solution 2 & Culture medium preparation (τ = 16 hr) |
| M303 | SeedHoldTank | Combines seed streams |
| R301 | SeedTrain | 5-stage seed train for inoculum expansion (T = 32°C) |
| M304 | MixTank | Glucose feed preparation (τ = 16 hr) |
| T301 | StorageTank | Glucose feed holding (τ = 136 hr) |
| T302 | AmmoniaStorageTank | Ammonia feed storage |
| R302 | AeratedFermentation | Main production fermenter (Fed-batch) |

### Fermentation Parameters
Same as Config 1:

| Parameter | Value | Units |
|-----------|-------|-------|
| `titer_LegHb` | 7.27 | g/L |
| `productivity_LegHb` | 0.101 | g/L/hr |
| `yield_LegHb (Y_p)` | 0.0224 | wt/wt |
| `Y_b (biomass yield)` | 0.43 | wt/wt |
| `T_operation` | 305.15 (32°C) | K |
| `V_max` | 500 | m³ |

### Control Specifications
- **R302 Specification (`update_reaction_time_and_yield`):**
  - Calculates residence time (`tau`): `target_titer / target_productivity`
  - Updates reaction product yield based on `target_yield`

## 2. Downstream Process (Traditional Purification)

### Cell Disruption & Solid-Liquid Separation
| Unit ID | Type | Description |
|---------|------|-------------|
| S401 | CellDisruption | High-pressure homogenizer |
| S402 | Filtration (MF preset) | Solid/liquid separation |
| S403 | ScrewPress | Cell mass dewatering (99.9% split, 1% moisture) |
| S404 | Splitter | Residual cell mass removal |

### Primary Ultrafiltration/Diafiltration
| Unit ID | Type | Description |
|---------|------|-------------|
| M401 | MixTank | DfUltraBuffer preparation (4× feed water) |
| H401 | HXutility | Cool to 5°C |
| U401 | Diafiltration (UF preset) | Primary buffer exchange and concentration |

**U401 Parameters:**
- Target Product: Leghemoglobin
- Salt retention: Low (salts wash through)

### Ion Exchange Chromatography
| Unit ID | Type | Description |
|---------|------|-------------|
| M402 | MixTank | IXEquilibriumBuffer preparation |
| M403 | MixTank | IXElutionBuffer preparation |
| M404 | MixTank | IXRegenerationSolution preparation |
| H402-H404 | HXutility | Cool buffers to 5°C |
| U402 | ResinColumn2 | Ion exchange capture step (IonExchange preset) |

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
| Unit ID | Type | Description |
|---------|------|-------------|
| M405 | MixTank | DfNanoBuffer preparation (2× IX eluate water) |
| H405 | HXutility | Cool to 5°C |
| U403 | Diafiltration (NF preset) | Polishing/concentration step |

### Final Concentration
| Unit ID | Type | Description |
|---------|------|-------------|
| U404 | Ultrafiltration | Final product concentration |

**U404 Control Specification (`U404_adjust_water_recovery`):**
- **Target:** LegHb content = 7.5% (range 6-9%)
- **Method:** Uses `flexsolve.IQ_interpolation` to optimize water recovery
- **Salt_Retention:** Fixed at 0.95 (keep salts with product)
- **Search range:** 50-99% water removal

## 3. Waste Treatment & Neutralization

### Waste Collection
| Unit ID | Type | Description |
|---------|------|-------------|
| M501 | MixTank | Combines permeate streams from S403, U401, U402, U403, U404 |
| T501 | SulfuricAcidStorageTank | Acid for IX regeneration neutralization |
| M502 | NeutralizationTank1 | Neutralizes IX regeneration waste (T = 20°C) |

**Acid Requirement:** Stoichiometric 1:2 ratio with NaOH plus 0.1% excess

### Water Recovery
| Unit ID | Type | Description |
|---------|------|-------------|
| S501 | ReverseOsmosis | Primary wastewater treatment |
| S502 | Splitter | Salt waste separation (99% split) |
| S503 | ReverseOsmosis | Secondary water recovery |

### Final Product Cooling
| Unit ID | Type | Description |
|---------|------|-------------|
| H406 | HXutility | Cool product to 0°C |

## 4. Facilities

| Unit ID | Type | Description |
|---------|------|-------------|
| CT | CoolingTower | Cooling water system |
| CWP | ChilledWaterPackage | Chilled water for process cooling |
| BT | BoilerTurbogenerator | Steam & power from cell mass waste |
| PWC | ProcessWaterCenter | Water recycling and makeup management |

**Makeup Water Streams:** cooling tower, M301-M304 waters, buffer waters (Water1-8), boiler makeup

**Process Water Streams:** RO permeates (S501, S503) plus makeup streams

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

| Aspect | Config 1 | Config 2 |
|--------|----------|----------|
| Primary Separation | Biomass harvest centrifuges | Direct cell disruption |
| Chromatography | None (simplified process) | Ion Exchange (ResinColumn2) |
| Polishing | HTST Pasteurization | Nanofiltration DF |
| Final Formulation | Antioxidant + dilution mixing | Direct concentration |
| Product Temp | 4°C | 0°C |
| Waste Treatment | Debris to disposal | Cell mass to boiler |

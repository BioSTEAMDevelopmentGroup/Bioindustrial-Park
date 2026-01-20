# LegHemoglobin Process Flow (`LegHb`) - Config 1

**Source:** `v1/LegHb/system/_config1.py`  
**System ID:** `LegHb_sys`

This configuration implements an **intracellular Leghemoglobin recovery** process optimized for minimal purification steps while meeting product specifications. The design emphasizes biomass harvest, high-pressure homogenization, and UF/DF purification with HTST pasteurization.

## 1. Upstream Process (Fermentation)

The upstream section prepares inoculum and runs main fermentation to produce intracellular Leghemoglobin.

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

| Parameter | Value | Units |
|-----------|-------|-------|
| `titer_LegHb` | 7.27 | g/L |
| `productivity_LegHb` | 0.101 | g/L/hr |
| `yield_LegHb (Y_p)` | 0.0224 | wt/wt |
| `Y_b (biomass yield)` | 0.43 | wt/wt |
| `theta_O2` | 0.5 | % saturation |
| `agitation_power` | 0.985 | kW/m³ |
| `T_operation` | 305.15 (32°C) | K |
| `V_max` | 500 | m³ |

### Control Specifications
- **R302 Specification (`update_reaction_time_and_yield`):**
  - Calculates residence time (`tau`): `target_titer / target_productivity`
  - Updates reaction product yield based on `target_yield`

## 2. Downstream Process (Intracellular LegHb Recovery)

This configuration uses a **6-step intracellular recovery process**:

### Step 1: Biomass Harvest (Primary Separation)
| Unit ID | Type | Description |
|---------|------|-------------|
| C401 | SolidsCentrifuge | Disk stack centrifuge densifying biomass to ~45% dry solids |
| M401 | MixTank | Wash buffer preparation (1.5 diavolumes) |
| H401_wash | HXutility | Cool wash buffer to 10°C |
| M402_wash | MixTank | Cell wash mixer |
| C402 | SolidsCentrifuge | Washed cell separation |

- **Cell capture:** 98%
- **Moisture content:** 55%

### Step 2: Cell Disruption (High-Pressure Homogenization)
| Unit ID | Type | Description |
|---------|------|-------------|
| S401 | CellDisruption | Multi-pass HPH at 1000 bar (800-1200 bar spec) |
| H401 | HXutility | Immediate post-valve cooling to <15°C |

- **Cell disruption efficiency:** 90%
- **Operating pressure:** 1000 bar (100 MPa)

### Step 3: Lysate Clarification (Debris Removal)
| Unit ID | Type | Description |
|---------|------|-------------|
| S402 | SolidsCentrifuge | High-speed centrifuge for bulk solids |
| S403 | Filtration (MF preset) | Depth filtration series (30µm → 5µm → 0.5µm) |
| S404 | ScrewPress | Combined debris dewatering |

- **Solid capture efficiency:** 95%
- **Key goal:** Remove Mannoprotein for ≥65% protein purity

### Step 4: UF/DF Purification
| Unit ID | Type | Description |
|---------|------|-------------|
| M403 | MixTank | DF buffer preparation (6 diavolumes) |
| H402 | HXutility | Cool DF buffer to 5°C |
| U401 | Diafiltration (UF preset) | TFF with 3-10 kDa MWCO membrane |
| U404 | Ultrafiltration | Final concentration step |

- **LegH retention:** 99%
- **Salt retention:** 5% (salts wash through)
- **VCF:** 5-10X
- **Diavolumes:** 6 (range 5-7)

### Step 5: Thermal Stabilization (HTST Pasteurization)
| Unit ID | Type | Description |
|---------|------|-------------|
| H403 | HXutility | Heat to 72°C (70-75°C range) |
| T401 | MixTank | Hold time tank (~30 seconds) |

- Inactivates pathogens
- Precipitates unstable host proteins
- LegH (stable heme-protein) remains soluble

### Step 6: Final Formulation & Cooling
| Unit ID | Type | Description |
|---------|------|-------------|
| M404 | MixTank | Formulation mixer with antioxidant and dilution water |
| H406 | HXutility | Final rapid cooling to 4°C |

**Formulation Targets:**
- Total solids: 20% (max 24%)
- LegHb content: 7.5% (range 6-9%)
- Sodium ascorbate: 0.1% w/w

### Control Specifications
- **U404 Specification (`U404_adjust_water_recovery`):**
  - Fixed Salt_Retention = 0.95 for concentration
- **M404 Specification (`update_formulation`):**
  - Calculates dilution water from total solids target
  - Sets antioxidant at 0.1% w/w of final product

## 3. Wastewater Treatment & Facilities

### Waste Collection
| Unit ID | Type | Description |
|---------|------|-------------|
| M501 | MixTank | Combines all liquid waste streams |
| S501 | ReverseOsmosis | Water recovery from wastewater |
| M502_compat | Mixer | IX buffer compatibility (zero flow) |
| S503 | ReverseOsmosis | Overflow RO (placeholder) |
| M503_debris | Mixer | Debris disposal route |

### Facilities
| Unit ID | Type | Description |
|---------|------|-------------|
| CT | CoolingTower | Cooling water system |
| CWP | ChilledWaterPackage | Chilled water for process cooling |
| BT | BoilerTurbogenerator | Steam & power (80% boiler, 85% turbogenerator) |
| PWC | ProcessWaterCenter | Water recycling and makeup management |

**Note:** Debris stream (S404-0) routes to disposal, not boiler, due to protein thermodynamic property limitations.

## 4. Exported Functions

### `create_LegHb_system(ins, outs, use_area_convention=False)`
Factory function creating the complete production system.

### `set_production_rate(system, target_production_rate_kg_hr, verbose=True)`
Adjusts system inputs to achieve target LegHb_3 production rate using global scaling factor with flexsolve IQ_interpolation.

### `check_LegHb_specifications(product_stream)`
Verifies LegHb_3 product stream meets composition specifications:
- Fat (OleicAcid): 0-2%
- Carbohydrates: 0-4%
- Leghemoglobin: 6-9%
- Total Solids: 0-24%
- Protein Purity: ≥65%

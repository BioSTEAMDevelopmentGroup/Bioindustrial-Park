# LegHemoglobin Process Flow (`LegHb`) - Config 1

**Source:** `v1/LegHb/system/_config1.py`  
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

| Function | Purpose |
|----------|---------|
| `get_fermentation_parameters()` | Returns validated fermentation parameters |
| `create_fermentation_reactions()` | Creates reaction systems for LegHb production |

---

## Area 200: Media Preparation

**Function:** `create_area_200_media_prep(SeedIn1, SeedIn2, CultureIn, Glucose, NH3_25wt)`

Purpose: Prepare all feed solutions for fermentation.

### Unit Operations

| Unit ID | Type | Description |
|---------|------|-------------|
| M301 | MixTank | Seed Solution 1 preparation (τ = 16 hr) |
| M302 | MixTank | Seed Solution 2 & Culture medium preparation (τ = 16 hr) |
| M303 | SeedHoldTank | Combines seed streams |
| M304 | MixTank | Glucose feed preparation (τ = 16 hr) |
| T301 | StorageTank | Glucose feed holding (τ = 136 hr) |
| T302 | AmmoniaStorageTank | Ammonia feed storage |

### Outputs

| Key | Stream | Description |
|-----|--------|-------------|
| `seed_out` | M303Out | Combined seed solution → Area 300 |
| `glucose_out` | T301Out | Diluted glucose solution → Area 300 |
| `ammonia_out` | T302Out | Ammonia solution → Area 300 |

---

## Area 300: Conversion (Fermentation)

**Function:** `create_area_300_conversion(seed_in, glucose_in, ammonia_in, vent1, vent2, reactions, params)`

Purpose: Convert glucose to LegHb through microbial fermentation.

### Unit Operations

| Unit ID | Type | Description |
|---------|------|-------------|
| R301 | SeedTrain | 5-stage seed train for inoculum expansion (T = 32°C) |
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

### Outputs

| Key | Stream | Description |
|-----|--------|-------------|
| `broth_out` | Broth | Cell suspension → Area 400 |

---

## Area 400: Recovery (Biomass Harvest & Cell Disruption)

**Function:** `create_area_400_recovery(broth_in, DfUltraBuffer_wash)`

Purpose: Release and clarify intracellular LegHb from yeast cells.

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

### Outputs

| Key | Stream | Description |
|-----|--------|-------------|
| `clarified_lysate` | ClarifiedLysate | Cell-free protein solution → Area 500 |
| `spent_media` | SpentMedia | To wastewater (Area 900) |
| `wash_effluent` | WashEffluent | To wastewater (Area 900) |
| `dehydrated_debris` | DehydratedDebris | To disposal (Area 900) |
| `press_liquor` | PressLiquor | To wastewater (Area 900) |

---

## Area 500: Purification (UF/DF Membrane Separation)

**Function:** `create_area_500_purification(clarified_lysate, DfUltraBuffer)`

Purpose: Concentrate LegHb and remove impurities using tangential flow filtration.

### Unit Operations

| Unit ID | Type | Description |
|---------|------|-------------|
| M403 | MixTank | DF buffer preparation (6 diavolumes) |
| H402 | HXutility | Cool DF buffer to 5°C |
| U401 | Diafiltration (UF preset) | TFF with 3-10 kDa MWCO membrane |
| U404 | Diafiltration | Final concentration step |

### Process Parameters

| Parameter | Value |
|-----------|-------|
| LegH retention | 99% |
| Salt retention | 5% (salts wash through) |
| VCF | 5-10X |
| Diavolumes | 6 (range 5-7) |

### Control Specifications

- **U404 Specification (`U404_adjust_water_recovery`):**
  - Fixed Salt_Retention = 0.95 for concentration

### Outputs

| Key | Stream | Description |
|-----|--------|-------------|
| `concentrated_product` | ConcentratedLegH | Purified LegHb → Area 600 |
| `uf_permeate` | UFPermeate | To wastewater (Area 900) |
| `concentration_permeate` | ConcentrationPermeate | To wastewater (Area 900) |

---

## Area 600: Formulation (Thermal Treatment & Final Product)

**Function:** `create_area_600_formulation(concentrated_product, LegHb_3)`

Purpose: Ensure microbiological safety and achieve final product specification.

### Step 1: Thermal Stabilization (HTST Pasteurization)

| Unit ID | Type | Description |
|---------|------|-------------|
| H403 | HXutility | Heat to 72°C (70-75°C range) |
| T401 | MixTank | Hold time tank (~30 seconds) |

- Inactivates pathogens
- Precipitates unstable host proteins
- LegHb (stable heme-protein) remains soluble

### Step 2: Final Formulation & Cooling

| Unit ID | Type | Description |
|---------|------|-------------|
| M404 | MixTank | Formulation mixer with antioxidant and dilution water |
| H406 | HXutility | Final rapid cooling to 4°C |

### Formulation Targets

| Specification | Target | Range |
|---------------|--------|-------|
| Total solids | 20% | max 24% |
| LegHb content | 7.5% | 6-9% |
| Sodium ascorbate | 0.1% w/w | - |

### Control Specifications

- **M404 Specification (`update_formulation`):**
  - Calculates dilution water from total solids target
  - Sets antioxidant at 0.1% w/w of final product

### Outputs

| Key | Stream | Description |
|-----|--------|-------------|
| `product` | LegHb_3 | Final product stream |

---

## Area 900: Facilities & Utilities

**Function:** `create_area_900_facilities(area_400, area_500, effluent1, use_area_convention)`

Purpose: Provide utilities and manage wastewater.

### Wastewater Treatment

| Unit ID | Type | Description |
|---------|------|-------------|
| M501 | MixTank | Combines all liquid waste streams |
| S501 | ReverseOsmosis | Water recovery from wastewater |
| M503_debris | Mixer | Debris disposal route |

### Utility Systems

| Unit ID | Type | Description |
|---------|------|-------------|
| CT | CoolingTower | Cooling water system |
| CWP | ChilledWaterPackage | Chilled water for process cooling |
| PWC | ProcessWaterCenter | Water recycling and makeup management |

**Note:** Debris stream routes to the boiler via the prefers-local `BoilerTurbogenerator` subclass with guarded emissions enthalpy updates and glucose-based thermo fallbacks for large biomolecules.

---

## Exported Functions

### System Creation

```python
create_LegHb_system(ins, outs, use_area_convention=False)
```
Factory function creating the complete production system using modular Process Area functions.

### Area Functions

| Function | Description |
|----------|-------------|
| `create_area_200_media_prep()` | Create Media Preparation area |
| `create_area_300_conversion()` | Create Fermentation area |
| `create_area_400_recovery()` | Create Recovery area |
| `create_area_500_purification()` | Create Purification area |
| `create_area_600_formulation()` | Create Formulation area |
| `create_area_900_facilities()` | Create Facilities area |

### Design Specifications

```python
set_production_rate(system, target_production_rate_kg_hr, verbose=True)
```
Adjusts system inputs to achieve target LegHb_3 production rate using global scaling factor with flexsolve IQ_interpolation.

```python
check_LegHb_specifications(product_stream)
```
Verifies LegHb_3 product stream meets composition specifications:
- Fat (OleicAcid): 0-2%
- Carbohydrates: 0-4%
- Leghemoglobin: 6-9%
- Total Solids: 0-24%
- Protein Purity: ≥65%

---

## Change History

| Date | Version | Description |
|------|---------|-------------|
| 2026-01-21 | 2.0 | Modular refactoring into Process Areas |
| 2025-06-04 | 1.0 | Initial monolithic implementation |

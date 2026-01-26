# N-HemDx Process Module (Integrated Full System)

**Source:** `v1/HemDx/system/_config1.py`  
**System ID:** `NHemDx_sys`  
**Product:** N-HemoDextrin (Nicotinamide-Stabilized Hemodextrin)

## Process Overview
Complete production system for **N-HemDx** using *Corynebacterium glutamicum*, integrating:
- Upstream media preparation
- Fermentation (2-phase: Growth + Production)
- Split-stream downstream processing (DSP)
- γ-Cyclodextrin complexation
- Nicotinamide stabilization
- Full wastewater treatment and utility systems

**Note:** This configuration uses a single **Adsorption** ResinColumn with realistic impurity distribution parameters on the combined stream. Buffer prep uses solute feeds with CV-based water addition.

## Key Performance Metrics

| Metric                         | Value            |
| ------------------------------ | ---------------- |
| N-HemoDextrin                  | 23.81 kg/hr      |
| HemoDextrin                    | 1.25 kg/hr       |
| Free Heme (Heme_b + Heme_b_In) | 0.047 kg/hr      |
| Total Product Stream           | 360.09 kg/hr     |
| Annual Production              | 2,880.75 MT/year |

## Process Areas

### Area 200: Media Preparation
- M201-M204: Seed and glucose solution mixing
- T201-T202: Storage tanks (glucose, ammonia)
- R301: Seed train

### Area 300: Conversion (Fermentation)
- **R302**: Aerated Fermentation (*Corynebacterium glutamicum*)
  - 2-Stage logic: Growth (72h) + Production (90h)
  - Operating Temp: 30°C
  - DO Setpoint: 45%
  - Secretion Fraction (SF): 45% extracellular

### Area 400: Clarification & Cell Disruption (Split-Stream)
- **S401**: Primary Centrifuge → Supernatant + CellCream
- **S404**: Microfiltration (MF preset) → FilteredSupernatant
- **Cell Cream Processing (LegHb sequence)**
  - **M401**: Wash buffer preparation
  - **H402**: Cool wash buffer to 10°C
  - **M402**: Cell wash mixer
  - **C402**: Washed cell centrifuge
  - **S402**: Cell Disruption (HPH, 1000 bar)
    - Target: `Corynebacterium_glutamicum`
    - Disruption efficiency: 55%
    - Component fractions: Protein (0.45), Cellulose (0.22), Xylan (0.15), OleicAcid (0.08), RNA (0.10)
  - **H401**: Post-valve cooling to 15°C
  - **S403**: Debris centrifuge → CrudeLysate + CellDebris
  - **S405**: Microfiltration (MF preset) → ClarifiedLysate
- **M404**: Debris mixer → CellDebrisRaw
- **S406**: ScrewPress for debris dewatering
  - Outputs: DehydratedDebris → Boiler (M902)
  - PressLiquor → Wastewater Treatment
- **M405**: Mixer → FilteredSupernatant + ClarifiedLysate → CombinedLysate

### Area 500: Capture & Purification (Single Adsorption Column)
Single ResinColumn (Adsorption preset) treating the combined stream:

**U501 (Combined Path)**:
- Preset: Adsorption
- Feed: CombinedLysate (FilteredSupernatant + ClarifiedLysate)
- Purpose: Capture heme products with realistic impurity distribution
- **Impurity Distribution Parameters** (2026-01-26 upgrade):
  - `NonTarget_Removal`: 0.99 (99% to flowthrough)
  - `Wash_Impurity_Carryover`: 0.02 (2% to wash stream)
  - `Regen_Impurity_Carryover`: 0.01 (1% to regen stream)
- Buffer prep: M501/M502/M503 with NaCl, NaOH, and Ethanol solute feeds
- CV-based water addition uses wash/elution/regeneration volumes from `U501`

### Area 600: Concentration
- **U601**: Diafiltration (NF preset)
  - Heme retention: 98%
  - Concentration factor: 5×
- **M601**: DF buffer preparation (DfUltraBuffer solute + water)

### Area 700: Formulation
Creates the final N-HemoDextrin complex:

**Solution Prep:**
- **M701**: GammaCyclodextrin solution preparation
- **M703**: Nicotinamide solution preparation

**R702 (Complexation + Stabilization)**:
```
Heme_b + GammaCyclodextrin → HemoDextrin (95%)
Heme_b_In + GammaCyclodextrin → HemoDextrin (95%)
HemoDextrin + Nicotinamide → N-HemoDextrin (95%)
```

Dosing:
- γ-CD: 13.45× molar ratio (based on stoichiometry)
- Nicotinamide: 2× molar with 2.5% excess

### Area 800: Final Product
- **H801**: Cold storage conditioning (4°C)
- Output: NHemDx_Product

### Area 900: Wastewater Treatment & Facilities (LegHb Pattern)

**Wastewater Treatment System** (`bst.create_wastewater_treatment_system`):
- Collects: ResinFlowthrough, ResinWash, ResinRegen, NFPermeate, WashEffluent, PressLiquor
- Outputs: biogas, sludge, RO_treated_water, ProcessWaste
- `mockup=True`, `area=500`

**Utility Systems:**
- **CT**: CoolingTower
- **CWP**: ChilledWaterPackage
- **PWC**: ProcessWaterCenter (with makeup water streams)

**Solids Handling:**
- **M902**: Mixer → DehydratedDebris + WWT sludge → SolidsToBoiler
- **BT**: BoilerTurbogenerator
  - Inputs: SolidsToBoiler, biogas, makeup water, natural gas, lime, chems
  - Outputs: emissions, rejected_water, ash_disposal
  - Boiler efficiency: 80%, Turbogenerator efficiency: 85%

## Chemical Additions

| Chemical                 | Role                    | Stoichiometry         |
| ------------------------ | ----------------------- | --------------------- |
| γ-Cyclodextrin (MW 1297) | Carrier molecule        | 13.45:1 molar w/ heme |
| Nicotinamide (MW 122)    | Axial ligand stabilizer | 2:1 molar w/ 2.5% XS  |
| HemoDextrin (MW ~25)     | Intermediate complex    | Reaction product      |
| N-HemoDextrin (MW ~25)   | Final product           | Reaction product      |

## Related Files
- `_config1_bp.py`: Backup of original upstream-only config
- `_config1_DSP_Draft.py`: DSP draft (standalone, for reference)
- `_chemicals.py`: Chemical definitions including HemDx formulation chemicals (with Glucose model copying for WWT compatibility)
- `_streams.py`: Input stream definitions

## Last Updated
2026-01-26 - Upgraded ResinColumn impurity distribution, added ScrewPress, full WWT system with CT/CWP/PWC


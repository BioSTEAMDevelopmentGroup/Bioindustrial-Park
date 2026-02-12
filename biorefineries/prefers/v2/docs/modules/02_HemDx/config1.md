# N-HemDx Process Module (Integrated Full System)

**Source:** `biorefineries/prefers/v2/HemDx/system/_config1.py`  
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

This configuration uses **ResinColumnAdv (Adsorption)** with impurity distribution parameters on the combined lysate stream. Buffer prep uses CV-based water addition.

## Process Areas

### Area 200: Media Preparation
- M201-M204: Seed and glucose solution mixing
- T201-T202: Storage tanks (glucose, ammonia)
  - **S202**: Ammonia Splitter (Seed/Fermentation split control)
- R301: Seed train

### Area 300: Conversion (Fermentation)
- **R302**: Aerated Fermentation (*Corynebacterium glutamicum*)
  - 2-stage logic: Growth (72 hr) + Production (90 hr)
  - Operating temp: 30°C, DO setpoint: 45%
  - Secretion Fraction (SF): 0.45 baseline
  - **Internalized specification:** Empirical NH3 measurement and elasticity-based titer control
    - $Elasticity = 2.0$, $Tolerance = 0.5\%$
    - Yield clamped to $[0.0001, 0.05]$
    - Updates NH3 source stream and S202 split for consistency

### Area 400: Clarification & Cell Disruption (Split-Stream)
- **C401**: Centrifuge → CellCream + Supernatant (moisture_content = 0.40)
- **S404**: FiltrationAdv (MF) → FilteredSupernatant (solid_capture_eff = 0.85, cake_moisture = 0.30)
- **Cell Cream Processing**
  - **M401/H402/M402/C402**: Wash buffer prep + cooling + wash (C402 moisture_content = 0.55)
  - **S402**: Cell Disruption (HPH, 1000 bar, efficiency = 0.87)
    - Component fractions: Protein 0.45, Cellulose 0.22, Xylan 0.15, OleicAcid 0.08, RNA 0.10
  - **H401/C403**: Cooling + debris centrifuge (moisture_content = 0.20)
  - **S405**: FiltrationAdv (MF) → ClarifiedLysate
- **M404/S406**: Debris mixer + ScrewPress (split = 0.999, moisture_content = 0.001)
- **M405**: Combine FilteredSupernatant + ClarifiedLysate → CombinedLysate

### Area 500: Capture & Purification (Single Adsorption Column)
Single ResinColumn (Adsorption preset) treating the combined stream:

**U501 (Combined Path)**:
- Unit: `ResinColumnAdv` (preset = Adsorption)
- TargetProduct_IDs: Heme_b, Heme_b_In, ProtoporphyrinIX, ProtoporphyrinIX_In
- TargetProduct_Yield = 0.97, NonTarget_Removal = 0.97
- wash_CV = 3, elution_CV = 0.05, regeneration_CV = 0.05
- Wash/Regen impurity carryover: 0.02 / 0.01
- Buffer prep: M501/M502/M503 with NaCl, NaOH, Ethanol feeds
- CV-based water addition is driven by U501 CVs

### Area 600: Concentration
- **U601**: DiafiltrationAdv (NF preset)
  - TargetProduct_Retention = 0.95
  - Salt_Retention = 0.10
  - diavolumes = 5.0
- **M601**: DF buffer prep (DfUltraBuffer1 + Water8), water = 2× feed water

### Area 700: Formulation
Creates the final N-HemoDextrin complex:

**Solution Prep:**
- **M701**: GammaCyclodextrin solution preparation
- **M703**: Nicotinamide solution preparation

**R702 (Complexation + Stabilization)**:
```
0.0014711 Heme_b + 0.0197913 GammaCyclodextrin → HemoDextrin (95%)
HemoDextrin + 0.0029422 Nicotinamide → N-HemoDextrin (95%)
```

Dosing:
- γ-CD: molar ratio set by reaction stoichiometry; water = 10× γ-CD mass
- Nicotinamide: 2× molar with 2.5% excess; water = 5× nicotinamide mass

### Area 800: Final Product
- **U801**: DiafiltrationAdv (UF) for final concentration and salt removal
  - TargetProduct_Retention = 0.99, Salt_Retention = 0.10, diavolumes = 5
  - Spec targets 7.5 wt% N-HemoDextrin
- **H802**: HTST pasteurization at 74°C
- **M802**: Formulation MixTank (Antioxidant + Water12)
- **H803**: Final cooling to 4°C
- **T801**: Storage tank (1 week)
- Output: NHemDx_Product

### Area 900: Wastewater Treatment & Facilities (LegHb Pattern)

**Wastewater Treatment System** (`bst.create_wastewater_treatment_system`):
- Collects: PressLiquor, WashEffluent, ResinFlowthrough, ResinRegenWaste, ResinWash, NFPermeate, FinalPermeate
- Supplemental nutrients: `SupplementalNH4SO4` (NH3 = 0.5 kg/hr, (NH4)2SO4 = 50 kg/hr) and `SupplementalFeSO4`
- Outputs: biogas, sludge, RO_treated_water, ProcessWaste
- `mockup=True`, `area=500`

**Utility Systems:**
**Utility Systems:**
- CT: CoolingTower
- CWP: ChilledWaterPackage
- PWC: ProcessWaterCenter (recycled + makeup water)

**Solids Handling:**
- M902: Mixer → DehydratedDebris + WWT sludge → SolidsToBoiler
- BT: BoilerTurbogenerator (boiler_efficiency = 0.80, turbogenerator_efficiency = 0.85)

## Product Specifications (QA)

`check_HemDx_specifications()` enforces:
- Salt < 2.0 wt%
- Residual Cyclodextrin < 4.0 wt%
- Residual Nicotinamide < 2.0 wt%
- Intermediate HemDx < 2.0 wt%
- N-HemoDextrin 6.5-8.5 wt%

## Related Files
- `_chemicals.py`: HemDx chemical definitions
- `_streams.py`: Input stream definitions

### System Creation

```python
create_NHemDx_system(ins, outs)
```
Factory function creating the complete production system.

### Design Specifications

```python
set_production_rate(system, target_production_rate_kg_hr, verbose=True)
```
Adjusts system inputs to achieve target production rate. NH3 optimization is handled automatically by R302's internalized specification during simulation.

## Last Updated
2026-02-07 - Internalized titer + NH3 specification and advanced filtration/diafiltration units

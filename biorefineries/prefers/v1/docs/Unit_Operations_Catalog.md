# Unit Operations Catalog

**Source:** `v1/_units.py`  
**Abstract:** This catalog documents all custom unit operations in the shared library with comprehensive details on design parameters, sizing logic, and costing models.

---

## 1. SeedTrain

**Class:** `SeedTrain(bst.Unit)`  
**Description:** A 5-stage serial reactor train for inoculum expansion with geometric scale-up (10x volume per stage).

### 1.1 Design Parameters (Defaults)

| Parameter   | Default Value | Unit | Description                        |
| :---------- | :------------ | :--- | :--------------------------------- |
| `N_stages`  | 5             | -    | Number of reactor stages in series |
| `N_trains`  | 2             | -    | Number of parallel seed trains     |
| `tau_batch` | 16            | hr   | Cycle time for each batch          |
| `T`         | 305.15        | K    | Operating temperature (32°C)       |

### 1.2 Sizing Logic (Hardcoded Constants)

* **Geometric Scaling:** Each stage is 10x the volume of the previous stage (`vol *= 10`)
* **Total Volume:** Calculated as `outs[1].F_vol * tau_turnover`
* **Turnover Time:** `tau_batch / N_trains`

### 1.3 Costing Model

Uses decorator-based costing with power-law scaling:

| Equipment          | Base Cost ($) | Base Size (S)        | Exponent (n) | CE Index | BM Factor | Power (kW) |
| :----------------- | :------------ | :------------------- | :----------- | :------- | :-------- | :--------- |
| Pumps              | 24,800        | 43,149 kg/hr         | 0.8          | 522      | 2.3       | 40         |
| Stage #1 reactors  | 37,700        | 20 gal (0.076 m³)    | 0.7          | 522      | 1.8       | -          |
| Stage #2 reactors  | 58,300        | 200 gal (0.76 m³)    | 0.7          | 522      | 1.8       | -          |
| Stage #3 reactors  | 78,800        | 2,000 gal (7.6 m³)   | 0.7          | 522      | 1.8       | -          |
| Stage #4 reactors  | 176,000       | 20,000 gal (76 m³)   | 0.7          | 522      | 1.8       | -          |
| Stage #4 agitators | 13,000        | 20,000 gal (76 m³)   | 0.5          | 522      | 1.5       | 7.5        |
| Stage #5 reactors  | 590,000       | 200,000 gal (757 m³) | 0.7          | 522      | 1.8       | -          |
| Stage #5 agitators | 21,500        | 200,000 gal (757 m³) | 0.5          | 522      | 1.5       | 10         |

**Equation:** `Cost = N_trains × (CE/522) × Base_Cost × (S/S_base)^n`

---

## 2. AeratedFermentation

**Class:** `AeratedFermentation(bst.AeratedBioreactor)`  
**Description:** Main production fermenter extending BioSTEAM's AeratedBioreactor with custom reaction handling.

### 2.1 Design Parameters (Defaults)

| Parameter          | Default Value | Unit    | Description                                          |
| :----------------- | :------------ | :------ | :--------------------------------------------------- |
| `V_max_default`    | 500           | m³      | Maximum vessel volume                                |
| `dT_hx_loop`       | 8             | °C      | Coolant temperature difference in heat exchange loop |
| `Q_O2_consumption` | -460,240      | kJ/kmol | Heat of fermentation (oxygen basis)                  |
| `optimize_power`   | True          | -       | Automatically optimizes agitation power              |
| `batch`            | True          | -       | Batch operation mode                                 |

### 2.2 Sizing Logic

* Inherits sizing from `bst.AeratedBioreactor`
* Handles four reaction types: fermentation, cell growth, respiration, neutralization
* Vent streams: CO₂, O₂, N₂

### 2.3 Costing Model

* Inherits costing from `bst.AeratedBioreactor`

---

## 3. CellDisruption

**Class:** `CellDisruption(bst.Unit)`  
**Description:** High-pressure homogenizer for cell lysis. Converts biomass into constituent components using manual mass balance. Models power consumption and thermal effects.

### 3.1 Design Parameters (Defaults)

| Parameter                    | Default Value     | Unit | Description                                  |
| :--------------------------- | :---------------- | :--- | :------------------------------------------- |
| `Cell_ID`                    | 'Pichia_pastoris' | -    | Target cell type for disruption              |
| `cell_disruption_efficiency` | 0.55              | -    | Fraction of cells disrupted (50-60% typical) |
| `P_high`                     | 150×10⁵           | Pa   | High operating pressure (150 bar)            |
| `P_low`                      | 101,325           | Pa   | Outlet pressure (1 atm)                      |

### 3.2 Default Component Fractions

When `component_fractions=None`, disrupted cell mass splits into:

| Component    | Fraction |
| :----------- | :------- |
| Mannoprotein | 40%      |
| Glucan       | 50%      |
| OleicAcid    | 6%       |
| Chitin       | 3%       |
| RNA          | 1%       |

### 3.3 Sizing Logic (Hardcoded Constants)

* **Compression Ratio Per Stage:** 4.0
* **Number of Stages:** `ceil(log(P_high / P_inlet) / log(4.0))`
* **Power Calculation:**
  ```
  power_kW = (F_vol × 1000) × dP / (36000 × 0.85)
  ```
  where `dP = (P_high - P_low) / 1e5` and efficiency = 85%

### 3.4 Costing Model

* **Equation:**
  ```
  Cost = 90,000 × (F_vol × 1000)^0.5 × (dP/1000)^1.5
  ```
* **Bare Module Factor (BM):** 3.5

---

## 4. ReverseOsmosis

**Class:** `ReverseOsmosis(bst.Unit)`  
**Description:** Reverse osmosis unit for water recovery from brine based on fractional water recovery.

### 4.1 Design Parameters (Defaults)

| Parameter                 | Default Value | Unit    | Description                         |
| :------------------------ | :------------ | :------ | :---------------------------------- |
| `water_recovery`          | 0.987         | -       | Fraction of water recovered (98.7%) |
| `membrane_flux`           | 40            | L/m²/hr | Membrane flux rate                  |
| `membrane_cost_per_m2`    | 300           | $/m²    | Membrane cost                       |
| `membrane_lifetime_years` | 3             | years   | Membrane replacement cycle          |
| `plant_lifetime_years`    | 20            | years   | Total plant lifetime                |
| `operating_pressure_bar`  | 30            | bar     | Operating pressure                  |

### 4.2 Sizing Logic

* **Membrane Area:** `(F_vol × 1000) / membrane_flux`
* **Power:** Calculated via internal `bst.Pump` at operating pressure

### 4.3 Costing Model

* **Pump:** Standard `bst.Pump` costing
* **Membrane Replacement (OPEX):**
  ```
  Cost = (plant_lifetime / membrane_lifetime) × Area × cost_per_m2
  ```
* **Bare Module Factors:**
  - Pump: 2.3
  - Membrane replacement: 1.65
  - Hardware: 2.0

---

## 5. Diafiltration

**Class:** `Diafiltration(bst.Unit)`  
**Description:** Diafiltration unit for size-based separation, retaining large molecules (proteins) while allowing smaller ones (salts, water) through the permeate.

### 5.1 Design Parameters (Defaults)

| Parameter                        | Default Value | Unit    | Description                               |
| :------------------------------- | :------------ | :------ | :---------------------------------------- |
| `TargetProduct_Retention`        | 0.98          | -       | Retention of target product (98%)         |
| `Salt_Retention`                 | 0.05          | -       | Retention of salts (5%)                   |
| `OtherLargeMolecules_Retention`  | 0.98          | -       | Retention of other large molecules        |
| `DefaultSolutes_Retention`       | 0.08          | -       | Default retention for unlisted solutes    |
| `FeedWater_Recovery_to_Permeate` | 0.75          | -       | Feed water going to permeate (75%)        |
| `membrane_flux_LMH`              | 40.0          | L/m²/hr | Membrane flux (Ultra: 20-80, Nano: 10-40) |
| `TMP_bar1`                       | 2.0           | bar     | Transmembrane pressure (feed pump)        |
| `TMP_bar2`                       | 2.0           | bar     | Transmembrane pressure (recirc pump)      |
| `membrane_cost_USD_per_m2`       | 150.0         | $/m²    | Membrane cost                             |
| `membrane_lifetime_years`        | 2.0           | years   | Membrane replacement cycle                |
| `equipment_lifetime_years`       | 20.0          | years   | Equipment lifetime                        |
| `module_cost_factor`             | 25,000        | $       | Base module cost factor                   |
| `module_cost_exponent`           | 0.7           | -       | Scale exponent                            |
| `base_CEPCI`                     | 500.0         | -       | Base CE index for cost adjustment         |
| `recirculation_ratio`            | 10.0          | -       | Recirculation ratio (5-20 typical)        |

### 5.2 Sizing Logic

* **Membrane Area:** `permeate_vol_L_hr / membrane_flux_LMH`
* **Pump Efficiency:** 85%
* **Power:** Sum of feed pump and recirculation pump power

### 5.3 Costing Model

* **Membrane System (CAPEX):**
  ```
  Cost = module_cost_factor × Area^module_cost_exponent × (CE / base_CEPCI)
  ```
* **Membrane Replacement:**
  ```
  Cost = (equipment_lifetime / membrane_lifetime) × Area × cost_per_m2
  ```
* **Pumps:** `pump1.cost + pump2.cost`
* **Bare Module Factors:**
  - Membrane System: 1.65
  - Membrane replacement: 1.65
  - Pump: 1.89

---

## 6. Ultrafiltration

**Class:** `Ultrafiltration(bst.Unit)`  
**Description:** Single-pass ultrafiltration unit. Inherits retention logic from Diafiltration but without diafiltration buffer stream.

### 6.1 Design Parameters

Inherits all defaults from `Diafiltration`, with:

| Parameter | Default Value | Unit | Description                         |
| :-------- | :------------ | :--- | :---------------------------------- |
| `TMP_bar` | 2.0           | bar  | Single TMP value (no recirculation) |

### 6.2 Sizing & Costing

* Same as Diafiltration, but with single feed pump only

---

## 7. IonExchangeCycle

**Class:** `IonExchangeCycle(bst.Unit)`  
**Description:** Complete Ion Exchange (IEX) chromatography cycle as a steady-state equivalent. Models loading, washing, elution, and regeneration steps.

### 7.1 Design Parameters (Defaults)

| Parameter                       | Default Value      | Unit  | Description                        |
| :------------------------------ | :----------------- | :---- | :--------------------------------- |
| `cycle_time_hr`                 | 4.0                | hr    | Total cycle duration               |
| `equilibration_CV`              | 5.0                | CV    | Equilibration buffer volume        |
| `wash_CV`                       | 5.0                | CV    | Wash buffer volume                 |
| `elution_CV`                    | 3.0                | CV    | Elution buffer volume              |
| `regeneration_CV`               | 5.0                | CV    | Regeneration solution volume       |
| `resin_DBC_g_L`                 | 50.0               | g/L   | Dynamic Binding Capacity           |
| `load_safety_factor`            | 0.8                | -     | Loading safety factor (80% of DBC) |
| `TargetProduct_IDs`             | ('Leghemoglobin',) | -     | Target product chemical IDs        |
| `TargetProduct_Yield`           | 0.95               | -     | Product recovery (95%)             |
| `BoundImpurity_IDs`             | ('Heme_b',)        | -     | Bound impurity IDs                 |
| `BoundImpurity_Removal`         | 0.93               | -     | Impurity removal (93%)             |
| `NonBinding_Carryover`          | 0.04               | -     | Non-binding carryover (4%)         |
| `resin_cost_USD_per_L`          | 30.0               | $/L   | Resin cost                         |
| `resin_lifetime_years`          | 5.0                | years | Resin replacement cycle            |
| `column_hardware_cost_factor`   | 30,000             | $     | Hardware cost factor               |
| `column_hardware_cost_exponent` | 0.6                | -     | Hardware scale exponent            |

### 7.2 Sizing Logic

* **Resin Volume:**
  ```
  resin_volume_L = (target_mass_kg × 1000) / (DBC × safety_factor)
  ```
* **Pump Pressure:** 4 bar (400 kPa)

### 7.3 Costing Model

* **Column Hardware:**
  ```
  Cost = (CE/500) × factor × resin_volume^exponent
  ```
* **Resin:**
  ```
  Cost = cost_per_L × resin_volume × (plant_life / resin_life)
  ```
* **Pump:** Standard `bst.Pump` costing
* **Bare Module Factors:**
  - IEX Column Hardware: 2.5
  - IEX Resin: 1.5
  - Pump: 1.89

---

## 8. AmmoniaStorageTank

**Class:** `AmmoniaStorageTank(bst.StorageTank)`  
**Description:** Storage tank for ammonia solution with decorator-based costing.

### 8.1 Costing Model

| Equipment | Base Cost ($) | Base Size (S) | Unit  | Exponent (n) | CE Index | BM Factor |
| :-------- | :------------ | :------------ | :---- | :----------- | :------- | :-------- |
| Tank      | 196,000       | 1,171         | kg/hr | 0.7          | 522      | 2.0       |

---

## 9. SulfuricAcidStorageTank

**Class:** `SulfuricAcidStorageTank(bst.StorageTank)`  
**Description:** Storage tank for sulfuric acid with pump.

### 9.1 Costing Model

| Equipment | Base Cost ($) | Base Size (S) | Unit  | Exponent (n) | CE Index | BM Factor | Power (kW) |
| :-------- | :------------ | :------------ | :---- | :----------- | :------- | :-------- | :--------- |
| Tank      | 96,000        | 1,981         | kg/hr | 0.7          | 522      | 1.5       | -          |
| Pump      | 7,493         | 1,981         | kg/hr | 0.8          | 522      | 2.3       | 0.5        |

---

## 10. SeedHoldTank

**Class:** `SeedHoldTank(bst.Mixer)`  
**Description:** Hold tank for seed solution before fermentation.

### 10.1 Costing Model

| Equipment | Base Cost ($) | Base Size (S) | Unit  | Exponent (n) | CE Index | BM Factor | Power (kW) |
| :-------- | :------------ | :------------ | :---- | :----------- | :------- | :-------- | :--------- |
| Pump      | 8,200         | 43,149        | kg/hr | 0.8          | 522      | 2.3       | 7.457      |
| Tank      | 439,000       | 40,414        | kg/hr | 0.7          | 522      | 1.8       | -          |
| Agitator  | 31,800        | 40,414        | kg/hr | 0.5          | 522      | 1.5       | 11.32      |

---

## 11. NeutralizationTank1

**Class:** `NeutralizationTank1(bst.Unit)`  
**Description:** Neutralization reactor for acid-base reactions with cooling.

### 11.1 Design Parameters (Defaults)

| Parameter            | Default Value | Unit  | Description               |
| :------------------- | :------------ | :---- | :------------------------ |
| `T`                  | 298.15        | K     | Target temperature (25°C) |
| `residence_time`     | 2.0           | hr    | Tank residence time       |
| `agitator_kW_per_m3` | 0.5           | kW/m³ | Agitator power intensity  |

### 11.2 Default Reactions

1. `H2SO4 + 2 NaOH -> Na2SO4 + 2 H2O` (X=1.0)
2. `H2SO4 + Na2SO4 -> 2 NaHSO4` (X=1.0)

### 11.3 Sizing Logic

* **Tank Volume:** `(water_mass_flow / 1000) × residence_time`
* **Agitator Power:** `tank_volume × agitator_kW_per_m3`
* **Cooling:** Via auxiliary `bst.HXutility`

### 11.4 Costing Model

* **Tank:**
  ```
  Cost = (CE/522) × 96,000 × (volume/10)^0.7
  ```
* **Agitator:**
  ```
  Cost = (CE/522) × 31,800 × (power/11.3)^0.5
  ```
* **Cooler:**
  ```
  Cost = (CE/522) × 50,000 × (duty/1e6)^0.6
  ```
* **Bare Module Factors:**
  - Tank: 2.0
  - Agitator: 1.5
  - Cooler: 2.3

---

## Module Constants

| Constant    | Value     | Description             |
| :---------- | :-------- | :---------------------- |
| `_gal2m3`   | 0.003785  | Gallons to cubic meters |
| `_gpm2m3hr` | 0.227124  | GPM to m³/hr            |
| `_hp2kW`    | 0.7457    | Horsepower to kilowatts |
| `_Gcal2kJ`  | 4,184,000 | Gcal to kJ              |

---

*Last Updated: 2026-01-15*

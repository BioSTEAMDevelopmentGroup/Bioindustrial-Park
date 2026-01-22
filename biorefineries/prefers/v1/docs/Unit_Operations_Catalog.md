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

## 4. ResinColumn
**Class:** `ResinColumn`
**Detailed Specification:** [ResinColumn_spec.md](./units/ResinColumn_spec.md)

Polymorphic unit for ion exchange and adsorption processes.

| Preset        | Key Parameter | Default | Status    |
| :------------ | :------------ | :------ | :-------- |
| `IonExchange` | DBC           | 50 g/L  | Validated |
| `Adsorption`  | EBCT          | 5 min   | Validated |

## 5. ReverseOsmosis

**Class:** `ReverseOsmosis`  
**Description:** Reverse osmosis unit for water recovery with literature-validated parameters.
**Detailed Specification:** [ReverseOsmosis_spec.md](./units/ReverseOsmosis_spec.md)

### 4.1 Summary of Validated Parameters

| Parameter            | Default | Unit   | Status             |
| :------------------- | :------ | :----- | :----------------- |
| `water_recovery`     | 0.85    | -      | 75-95% (Ind. Std.) |
| `membrane_flux`      | 40      | LMH    | Validated          |
| `operating_pressure` | 25      | bar    | Conservative       |
| `specific_energy`    | 3.0     | kWh/m³ | New Metric         |

> For full design basis, equations, and references, see the [Detailed Specification](./units/ReverseOsmosis_spec.md).

---


## 6. Diafiltration

**Class:** `Diafiltration(bst.Unit)`  
**Description:** Verified Diafiltration unit for size-based separation and buffer exchange.
**Detailed Specification:** [Diafiltration_spec.md](./units/Diafiltration_spec.md)

### 5.1 Factory Presets (`Diafiltration.from_preset(name, ...)`)

| Preset | Application                                        | Flux (LMH) | Pressure (bar) | Cost ($/m²) | Life (yr) |
| :----- | :------------------------------------------------- | :--------- | :------------- | :---------- | :-------- |
| `'UF'` | Protein Concentration (Ultrafiltration)            | 50.0       | 2.0            | 150         | 3.0       |
| `'NF'` | Buffer Exchange / Fine Separation (Nanofiltration) | 25.0       | 8.0            | 250         | 2.0       |

### 5.2 Summary of Validated Parameters

| Parameter                  | Default | Unit  | Description                  |
| :------------------------- | :------ | :---- | :--------------------------- |
| `membrane_flux_LMH`        | 40.0    | LMH   | Permeate flux (Validated)    |
| `membrane_cost_USD_per_m2` | 150.0   | $/m²  | Membrane replacement cost    |
| `membrane_lifetime_years`  | 2.0     | years | Replacement frequency        |
| `TargetProduct_Retention`  | 0.99    | -     | Product retention efficiency |

> See [Detailed Specification](./units/Diafiltration_spec.md) for full design basis.

---

## 7. Ultrafiltration

**Class:** `Ultrafiltration(bst.Unit)`  
**Description:** Single-pass ultrafiltration unit. Inherits retention logic from Diafiltration but without diafiltration buffer stream.

### 6.1 Factory Presets

Inherits presets from `Diafiltration`. Use `Ultrafiltration.from_preset('UF', ...)` or `('NF', ...)`.

### 6.2 Design Parameters

Inherits all defaults from `Diafiltration`, with:

| Parameter | Default Value | Unit | Description                         |
| :-------- | :------------ | :--- | :---------------------------------- |
| `TMP_bar` | 2.0           | bar  | Single TMP value (no recirculation) |

### 6.3 Sizing & Costing

* Same as Diafiltration, but with single feed pump only

---

## 8. BoilerTurbogenerator

**Class:** `BoilerTurbogenerator(bst.BoilerTurbogenerator)`  
**Description:** Prefers-local subclass with guarded emissions enthalpy update to handle robust combustion/energy balance when large biomolecules lack complete vapor-phase thermo.  
**Detailed Specification:** [BoilerTurbogenerator_spec.md](./units/BoilerTurbogenerator_spec.md)

---

## 9. AmmoniaStorageTank

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

## 12. VacuumPSA

**Class:** `VacuumPSA(bst.Unit)`  
**Description:** Vacuum Pressure Swing Adsorption (VPSA) unit for gas separation. Designed for separating H₂, CO, and C₂H₄ from syngas using zeolite 13X columns with cyclic adsorption/desorption.

### 12.1 Design Parameters (Defaults)

| Parameter                | Default Value | Unit   | Description                          |
| :----------------------- | :------------ | :----- | :----------------------------------- |
| `P_ads`                  | 6×10⁵         | Pa     | Adsorption pressure (6 bar)          |
| `P_des`                  | 1×10⁴         | Pa     | Desorption pressure (0.1 bar vacuum) |
| `cycle_time`             | 600           | s      | Total cycle time                     |
| `N_beds`                 | 2             | -      | Number of adsorbent beds             |
| `adsorbent_loading`      | 2.0           | mol/kg | Adsorbent capacity                   |
| `adsorbent_bulk_density` | 650           | kg/m³  | Bulk density of zeolite 13X          |
| `vacuum_efficiency`      | 0.70          | -      | Vacuum pump efficiency               |
| `adsorbent_cost`         | 5.0           | $/kg   | Cost of adsorbent material           |

### 12.2 Default Split Factors

Based on Zeolite 13X selectivity:
* **Product (H₂-rich):** 95% H₂, 30% CO, 20% C₂H₄, 10% CO₂, 85% N₂
* **Purge (CO/C₂H₄-rich):** Remainder

### 12.3 Sizing Logic

* **Adsorbent Mass:** 
  ```
  m_ads = (F_mol_feed × cycle_time) / loading
  ```
* **Bed Volume:** `m_ads / bulk_density`
* **Vacuum Power:**
  ```
  Power = F_purge × R × T × ln(P_ads/P_des) / efficiency
  ```

### 12.4 Costing Model

Uses power-law scaling:

| Equipment        | Base Cost ($)  | Base Size | Exponent (n) | BM Factor |
| :--------------- | :------------- | :-------- | :----------- | :-------- |
| Pressure Vessels | 50,000         | 10 m³     | 0.6          | 2.5       |
| Vacuum Pump      | 30,000         | 100 kW    | 0.7          | 1.8       |
| Adsorbent        | -              | -         | 1.0 (Linear) | 1.0       |
| Piping & Valves  | 15% of Vessels | -         | 1.0          | 1.0       |

---

## 13. Filtration

**Class:** `Filtration(bst.Unit)`  
**Description:** Continuous Rotary Drum Vacuum Filter (RDVF) for high-volume solid-liquid separation.
**Detailed Specification:** [Filtration_spec.md](./units/Filtration_spec.md)

### 13.1 Factory Presets (`Filtration.from_preset(name, ...)`)

| Preset | Application                                         | Loading (kg/m²/hr) | Pressure (bar) | Cost ($/m²) | Life (yr) |
| :----- | :-------------------------------------------------- | :----------------- | :------------- | :---------- | :-------- |
| `'MF'` | Cell Mass / Primary Clarification (Microfiltration) | 30.0               | 1.5            | 80          | 4.0       |
| `'UF'` | Virus Removal / Fine Separation (Ultrafiltration)   | 15.0               | 3.0            | 150         | 2.5       |

### 13.2 Summary of Validated Parameters

| Parameter                  | Default | Unit     | Description             |
| :------------------------- | :------ | :------- | :---------------------- |
| `solids_loading`           | 20.0    | kg/m²/hr | Filter surface loading  |
| `cake_moisture_content`    | 0.20    | -        | Residual moisture (20%) |
| `solid_capture_efficiency` | 0.99    | -        | Solids capture (99%)    |
| `power_per_m2`             | 1.0     | kW/m²    | Vacuum/drive power      |

---

*Last Updated: 2026-01-18*

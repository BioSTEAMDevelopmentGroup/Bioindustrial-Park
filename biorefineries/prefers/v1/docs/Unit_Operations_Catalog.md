# Unit Operations Catalog

**Source:** `v1/units.py`
**Abstract:** This catalog documents the custom unit operations available in the shared library. It lists **hardcoded default values** found in `__init__` and the **costing equations** used in `_cost`.

---

## 1. SeedTrain

**Class:** `SeedTrain(bst.Unit)`
**Description:** A 5-stage serial reactor train for inoculum expansion.

### Hardcoded Defaults
| Parameter   | Default Value   | Source Code |
| :---------- | :-------------- | :---------- |
| `N_stages`  | 5               | Line 92     |
| `N_trains`  | 2               | Line 95     |
| `tau_batch` | 16 hr           | Line 98     |
| `T`         | 32°C (305.15 K) | Line 106    |

### Costing Model
Uses decorator-based costing (`@cost`).
*   **Reactors (5 Stages):** Exponential scale-up ($37k \to \$590k$).
*   **Agitators (Stages 4 & 5):** Added separately ($13k \to \$21.5k$).
*   **Pumps:** Base cost $24,800.

---

## 2. AeratedFermentation

**Class:** `AeratedFermentation(bst.AeratedBioreactor)`
**Description:** Main production fermenter. Extends BioSTEAM's `AeratedBioreactor`.

### Hardcoded Defaults
| Parameter          | Default Value    | Notes                               |
| :----------------- | :--------------- | :---------------------------------- |
| `V_max_default`    | 500 m³           | Line 203                            |
| `dT_hx_loop`       | 8°C              | Coolant temperature difference      |
| `Q_O2_consumption` | -460,240 kJ/kmol | Heat of fermentation (Oxygen basis) |
| `optimize_power`   | True             | Automatically adjusts agitation     |

---

## 3. CellDisruption

**Class:** `CellDisruption(bst.Unit)`
**Description:** High-pressure homogenizer simulation.

### Hardcoded Defaults
| Parameter                    | Default Value    | Notes                       |
| :--------------------------- | :--------------- | :-------------------------- |
| `cell_disruption_efficiency` | 0.55 (55%)       | Fraction of cells disrupted |
| `P_high`                     | 150 bar (15 MPa) | Operating Pressure          |
| `P_low`                      | 1.01325 bar      | Outlet Pressure             |

### Default Component Fractions
If not provided, the cell mass splits into:
*   Mannoprotein: 40%
*   Glucan: 50%
*   OleicAcid: 6%
*   Chitin: 3%
*   RNA: 1%

### Costing Model
Scaling based on flow rate and pressure drop.
*   **Base Cost:** $90,000
*   **Equation:**
    $$Cost = 90000 \cdot (F_{vol, m3/hr} \cdot 1000)^{0.5} \cdot \left(\frac{P_{high} - P_{low}}{1000 \text{ bar?}}\right)^{1.5}$$
    *(Note: The code uses `dP/1000` where dP is in bar. 150 bar / 1000 = 0.15)*

---

## 4. ReverseOsmosis

**Class:** `ReverseOsmosis(bst.Unit)`

### Hardcoded Defaults
| Parameter                 | Default Value |
| :------------------------ | :------------ |
| `water_recovery`          | 0.987 (98.7%) |
| `membrane_flux`           | 40 L/m²/hr    |
| `operating_pressure_bar`  | 30 bar        |
| `membrane_cost_per_m2`    | $300          |
| `membrane_lifetime_years` | 3             |

### Costing Model
*   **Membrane Replacement:** Operating Cost (OPEX).
*   **Pump:** `bst.Pump` standard costing.

---

## 5. Diafiltration

**Class:** `Diafiltration(bst.Unit)`

### Hardcoded Defaults
| Parameter                  | Default Value |
| :------------------------- | :------------ |
| `TargetProduct_Retention`  | 0.98          |
| `Salt_Retention`           | 0.05          |
| `FeedWater_Recovery`       | 0.75          |
| `membrane_flux_LMH`        | 40.0          |
| `membrane_cost_USD_per_m2` | 150.0         |
| `module_cost_factor`       | 25,000        |

### Costing Model
*   **Membrane System (Capital):**
    $$Cost = 25000 \cdot (Area)^{0.7}$$
*   **Pumps:** Calculated separately for feed and recirculation.

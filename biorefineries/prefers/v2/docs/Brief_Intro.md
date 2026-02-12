**Project Brief: Comparative Techno-Economic Analysis (TEA) and Life Cycle Assessment (LCA) of Heme-Functionalized Ingredients for Meat Analogues**

**Date:** February 12, 2026  
**Prepared By:** Senior Process Engineer, Bioprocess Division  
**Subject:** Technical Overview of Leghemoglobin (LegHb) vs. N-HemoDextrin (HemDx) Production Routes — PREFERS v2

---

### 1. Executive Summary

This project evaluates the industrial viability and environmental footprint of two distinct bioprocess routes for manufacturing heme-based organoleptic agents used in meat analogues. The objective is to provide a rigorous engineering basis for process selection, focusing on:

* **LegHb** — recombinant leghemoglobin protein (~16 kDa) produced via *Pichia pastoris* fed-batch fermentation.
* **N-HemoDextrin (HemDx)** — free heme B (616.5 Da) produced by *Corynebacterium glutamicum*, stabilized by γ-Cyclodextrin inclusion complexation and Nicotinamide co-stabilization.

The analysis integrates stoichiometry, mass transfer limitations, downstream processing (DSP) complexity, capital intensity, and stochastic risk quantification. **PREFERS v2** models both routes using BioSTEAM with internalized design specifications for production targets and nitrogen balance, Singapore-specific IRAS depreciation schedules, and Monte Carlo uncertainty analysis.

---

### 2. Process Definitions and Topology

#### Route A: Leghemoglobin (LegHb) — The "Protein Route"

**Concept:** Production of a monomeric heme-binding protein (~16 kDa) structurally similar to myoglobin.  
**Reference Basis:** *Shao et al. (2022)*; *Tian et al. (2024)*; *Fraser et al. (2018)*.  
**v2 Source:** `LegHb/system/_config1.py`

* **Upstream Processing (USP):**
  * **Host:** Methylotrophic yeast (*Pichia pastoris*) or *Kluyveromyces marxianus*.
  * **Fermentation Mode:** High-cell-density fed-batch (R302, `AeratedFermentation`).
  * **Baseline Titer:** 5.0 g/L; productivity ≈ 0.069 g/L/hr (72 hr batch).
  * **Operating Conditions:** T = 32 °C, V_max = 500 m³, DO setpoint = 50% saturation.
  * **Key Constraint:** Heme biosynthesis flux must match globin translation rates to prevent accumulation of apo-leghemoglobin or toxic free heme.
  * **v2 Control:** Internalized specification on R302 — empirical NH₃ demand measurement and elasticity-based yield adjustment (Elasticity = 1.3, Tolerance = 0.1%).

* **Downstream Processing (DSP):**
  * **Target:** Intracellular soluble protein.
  * **Area 400 — Recovery:** Disc-stack centrifuge (C401, moisture 0.40) → cell wash (buffer + cooling to 10 °C) → washed centrifuge (C402, moisture 0.55) → high-pressure homogenization (S401, 1000 bar, efficiency 0.87) → post-valve cooling → debris centrifuge → microfiltration.
  * **Area 500 — Purification:** Diafiltration/Ultrafiltration (UF/DF) for concentration and buffer exchange (product retention ≥ 0.95).
  * **Area 600 — Formulation:** HTST pasteurization at 72 °C → antioxidant addition (sodium ascorbate) → cooling to 4 °C → storage (liquid concentrate, **no spray drying**).

**v2 Configurations:** Config 1 only (intracellular recovery).

---

#### Route B: N-HemoDextrin (HemDx) — The "Small Molecule Route"

**Concept:** Biosynthesis of free heme B (iron-protoporphyrin IX), followed by *ex-situ* stabilization via inclusion complexation with **γ-Cyclodextrin** and co-stabilization with **Nicotinamide**.  
**Reference Basis:** *Ko et al. (2021)*; *Yang et al. (2024)*; *Li et al. (2025)*.  
**v2 Source:** `HemDx/system/_config1.py` (base), `_config2.py` (intracellular), `_config3.py` (extracellular)

* **Upstream Processing (USP):**
  * **Host:** *Corynebacterium glutamicum* (GRAS status).
  * **Metabolic Engineering:** Overexpression of C5 pathway (glutamate → ALA) and downstream porphyrin pathway; deletion of heme-degrading enzymes.
  * **Fermentation Mode:** 2-phase fed-batch (R302, `AeratedFermentation`) — 72 hr growth + 90 hr production, T = 30 °C, DO = 45%.
  * **Key Constraint:** **Cytotoxicity.** Free heme is toxic to the host. The **Secretion Fraction (SF)** governs whether heme exits the cell or remains intracellular, driving the choice of DSP configuration.
  * **v2 Control:** Internalized specification on R302 — empirical NH₃ demand measurement and elasticity-based yield adjustment (Elasticity = 2.0, Tolerance = 0.5%, yield clamp [0.0001, 0.05]).

* **Downstream Processing (DSP) — Config 1 (Base, SF = 0.45):**
  * **Area 400 — Split-Stream Clarification:** Disc-stack centrifuge (C401) splits cell cream from supernatant. **Supernatant path:** FiltrationAdv (MF, S404). **Cell cream path:** wash → HPH (S402, 1000 bar) → debris centrifuge → FiltrationAdv (MF, S405). Both paths recombine (M405 → CombinedLysate).
  * **Area 500 — Capture:** ResinColumnAdv (Adsorption, U501) — captures heme and protoporphyrin species (yield 0.97, non-target removal 0.97).
  * **Area 600 — Concentration:** DiafiltrationAdv (NF preset, U601) — product retention 0.95, 5 diavolumes.
  * **Area 700 — Formulation:** γ-Cyclodextrin complexation (R702, 95% conversion) + Nicotinamide co-stabilization (95% conversion).
  * **Area 800 — Final Product:** DiafiltrationAdv (UF, U801, target 7.5 wt% N-HemoDextrin) → HTST at 74 °C → antioxidant formulation → cooling to 4 °C → storage (liquid concentrate, **no freeze drying**).

* **Config Variants:**

| Config | Secretion Fraction | Key DSP Change |
| ------ | ------------------ | -------------- |
| **Config 1** (Base) | 0.45 (range 0.2–0.8) | Full split-stream: supernatant filtration + cell disruption |
| **Config 2** (Intracellular) | 0.10 (range 0.0–0.3) | Supernatant → WWT (M504); only lysate path to capture |
| **Config 3** (Extracellular) | 0.90 (range 0.7–1.0) | Cell disruption train removed; only supernatant filtration to capture |

---

### 3. Techno-Economic Analysis (TEA) Framework

PREFERS v2 implements `PreFerSTEA` (inheriting `bst.TEA`) with **IRAS depreciation schedules** for Singapore-based facility economics. Key parameters: IRR = 18%, income tax = 17%, operating days = 333/yr.

#### 3.1 Capital Expenditure (CAPEX)

* **LegHb (Route A):** High CAPEX driven by bioreactor cooling requirements (high metabolic heat from aerobic yeast fermentation, `Q_O2` = −460 kJ/kmol O₂) and membrane filtration units (UF/DF) for protein concentration.
* **HemDx (Route B):** Moderate CAPEX. Resin adsorption columns (ResinColumnAdv) and membrane concentration units (DiafiltrationAdv) replace traditional chromatography skids. The complexation reactor and Nicotinamide dosing add formulation-specific equipment.

#### 3.2 Operating Expenditure (OPEX)

* **Titer Sensitivity:**
  * *LegHb:* v2 baseline titer = 5.0 g/L with uncertainty range [2.5, 7.5] g/L (Trunc-LogNormal). Higher titers reduce fermenter volume and improve MSP.
  * *HemDx:* Generally lower titers than protein routes. Viability is sensitive to secretion fraction — higher SF reduces DSP cost by eliminating cell disruption (Config 3).

* **Raw Materials:**
  * *LegHb:* Glucose ($0.42/kg), ammonia ($0.46/kg), trace metals (iron is critical).
  * *HemDx:* Glucose ($0.42/kg), ammonia ($0.115/kg), and crucially **γ-Cyclodextrin** ($6.0/kg) and **Nicotinamide** ($9.0/kg) for the formulation step — significant cost adders vs. simple protein buffers.

#### 3.3 Production Rate Scaling

The `set_production_rate()` function uses `flexsolve.IQ_interpolation` (root-finding) to find the exact input stream scaling factor that achieves the target production rate in the final product stream. Non-linear separations (membrane rejection, centrifuge recovery) make simple linear scaling inaccurate.

---

### 4. Life Cycle Assessment (LCA) Framework

#### 4.1 System Boundaries

* **Cradle-to-Gate:** From raw material extraction (glucose corn syrup, minerals) to the formulated liquid product at the factory gate.

#### 4.2 Key GWP Characterization Factors

| Input | GWP (kg CO₂-eq/unit) | Source |
| ----- | -------------------- | ------ |
| Electricity (SG) | 0.55 /kWh | Ecoinvent 3.11 |
| Glucose | 1.61 /kg | Ecoinvent 3.11 (GLO) |
| Ammonia (SEA) | 2.84 /kg | Ecoinvent 3.11 (SEA) |
| NaOH | 1.41 /kg | Ecoinvent 3.11 (RoW) |
| γ-Cyclodextrin | 4.75 /kg | Assumed |
| Nicotinamide | 7.26 /kg | FineChem2 prediction |

#### 4.3 Environmental Hotspots

* **Energy Consumption (GWP):**
  * **Bioreactor Aeration:** Heme synthesis is oxygen-intensive. High oxygen demand drives electricity usage for agitation and aeration — the largest single GWP contributor for both routes.

* **Water Usage:**
  * *LegHb:* Diafiltration requires multiple buffer volumes (diavolumes), generating significant wastewater routed to the integrated WWT system.
  * *HemDx:* Resin column buffer washes (CV-based) and NF/UF diavolumes are the primary water consumers. v2 uses aqueous-only DSP — **no organic solvents**.

* **Eutrophication:**
  * Spent fermentation broth containing residual nitrogen and phosphorus is treated by the integrated BioSTEAM wastewater treatment system. *C. glutamicum* broths are often cleaner than yeast broths but still require rigorous treatment.

---

### 5. Critical Engineering Recommendations

1. **Mass Transfer Strategy:** For both routes, dissolved oxygen is rate-limiting. Pressurized fermentation (moderate backpressure) is recommended to increase Oxygen Transfer Rate (OTR).
2. **Toxicity Management (HemDx):** Prioritize strains with engineered efflux pumps (e.g., *HrtBA*) to secrete heme into the medium, decoupling production from biomass volume. Config 3 (SF = 0.9) models this scenario.
3. **Formulation (LegHb):** HTST pasteurization at 72 °C preserves the heme-protein coordination bond. Temperatures above 80 °C may denature the globin, leading to color loss.
4. **Formulation (HemDx):** The γ-Cyclodextrin + Nicotinamide dual-stabilization step (R702) is critical for product stability. Ensure stoichiometric γ-CD and 2× molar Nicotinamide dosing with 2.5% excess.
5. **Regulatory (Safety):**
   * *LegHb:* Must demonstrate bio-equivalency to soy root leghemoglobin (GRAS Notice 737).
   * *HemDx:* Must validate that the cyclodextrin complex prevents free radical generation in the gut (heme toxicity mitigation) as suggested by *Li et al. (2025)*.

---

### 6. Block Flow Diagram (Simplified Comparison)

| Step | LegHb (Protein Route) | HemDx (Small Molecule Route) |
| --- | --- | --- |
| **Fermentation** | Yeast fed-batch (aerobic, 32 °C, 72 hr) | Bacteria 2-phase fed-batch (aerobic, 30 °C, 72+90 hr) |
| **Primary Separation** | Centrifugation → cell wash → washed centrifuge | Centrifugation → split-stream (supernatant MF + cell cream) |
| **Cell Processing** | HPH at 1000 bar → debris removal → MF | HPH at 1000 bar → debris removal → MF (Configs 1, 2 only) |
| **Capture / Purification** | UF/DF membrane concentration | ResinColumnAdv (Adsorption) + NF concentration |
| **Formulation** | HTST 72 °C + antioxidant | γ-CD complexation + Nicotinamide stabilization |
| **Final Product** | Liquid concentrate at 4 °C | UF polishing → HTST 74 °C → liquid concentrate at 4 °C |

---

### 7. Uncertainty Quantification

Both routes are evaluated under stochastic uncertainty using BioSTEAM's `Model` framework with Latin Hypercube sampling:

| Aspect | LegHb | HemDx |
| ------ | ----- | ----- |
| **Parameters** | 20 | Up to 23 (config-dependent) |
| **Key Distributions** | Trunc(LogNormal) for titer, prices, GWP; TruncNormal for yields/efficiencies | Same + Uniform for secretion fraction |
| **Titer Range** | [0.5×, 1.5×] baseline | [0.5×, 1.5×] baseline |
| **Key Metrics** | MSP ($/kg), TCI ($M), AOC ($M/yr), GWP (kg CO₂-eq/kg), composition (6 quality metrics) | MSP, TCI, AOC, GWP, composition (6 quality metrics) |
| **Pipeline** | `gen_data_base.py` → `gen_data_mc.py` → `gen_figure.py` | Same structure |

---

### 8. Product Specification Contracts

| Quality Parameter | LegHb Target | HemDx Target |
| ----------------- | ------------ | ------------ |
| Active Ingredient | LegHb 6–9 wt% | N-HemoDextrin 6.5–8.5 wt% |
| Fat (OleicAcid) | < 2 wt% | — |
| Carbohydrates | < 4 wt% | — |
| Protein Purity | ≥ 65% | — |
| Residual Salt | — | < 2.0 wt% |
| Residual γ-CD | — | < 4.0 wt% |
| Residual Nicotinamide | — | < 2.0 wt% |
| Total Solids | < 24 wt% | — |

---

*End of Brief*

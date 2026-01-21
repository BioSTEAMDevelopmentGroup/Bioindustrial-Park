# Biorefinery Simulation Optimization Master Plan

**Goal:** Transition from monolithic scripts to a modular, "Process Area" architecture (BioSTEAM Standard) while ensuring data integrity and documentation synchronization across all configurations.

---

## Phase 0: Codebase Hygiene & Flowsheet Pruning

**Objective:** Clean the unit library and remove legacy code to prevent technical debt from interfering with refactoring.

1.  **Refactor `v1\units.py`:**
    *   **Delete Legacy:** Remove old `ReverseOsmosis` and `Ultrafiltration` classes.
    *   **Rename & Standardize:** Rename `RO2` -> `ReverseOsmosis` and `ResinColumn2` -> `ResinColumn`. Do not use aliases; rename the classes directly.
    *   **Update Imports:** Scan all scripts and update references to use the new names.

2.  **Prune `LegHb_config1`:**
    *   **Remove unused units:** Delete `S503`, `M503_debris`, and `M502` (if specific to Config 2).
    *   **Route Solid Waste:** Ensure the debris stream from Area 400 is routed explicitly to `BT` (Boiler) for electricity generation.

3.  **Verification:** Run `LegHb_config1` to ensure no import errors or mass balance tears occur.

---

## Phase 1: Baseline Establishment & Benchmarking

**Objective:** Capture the performance of the *current* working state to create a validation standard.

1.  **Execution:** Run all current system configurations (`LegHb_config1`, `LegHb_config2`, `HemDx_config1`).
2.  **Data Logging:** Create a summary file `baseline_metrics.txt` recording:
    *   **Techno-Economic:** MSP, TCI.
    *   **Process:** Final Product Titer, Yield, Purity.
    *   **LCA:** GWP (Global Warming Potential).
3.  **Outcome:** A frozen benchmark to validate Phase 2 against.

---

## Phase 2: Modular Refactoring, Documentation & Verification

**Objective:** Refactor code into standard Process Areas (Subsystems), update documentation simultaneously, and verify correctness.

### Step 2a: Pilot Refactoring (`LegHb_config1`)

*   **Action:** Break the monolithic script into discrete functions following the **Area 100-900** convention. Use the specific hierarchy below (enforce this in Markdown docs as well):
    *   **Area 200 (Media Prep):** Mixing, Sterilization, Heat Exchange. (Distinct from conversion).
    *   **Area 300 (Conversion):** Fermentation, Seed Trains, Air/Gas handling.
    *   **Area 400 (Recovery/Concentration):** Centrifugation, Cell Disruption, Filtration (Primary).
    *   **Area 500 (Purification):** Resin columns (Chromatography), Ultrafiltration/Diafiltration.
    *   **Area 600 (Formulation):** Stabilization, Formulation tanks, Packaging.
    *   **Area 900 (Facilities):** bst.CoolingTower, bst.BoilerTurbineGenerator (CHP), bst.ChilledWaterPackage, Wastewater Treatment (if applicable).
*   **Constraint:** Ensure `bst.System` objects are linked correctly between areas.

### Step 2b: Documentation Alignment (New Constraint)

*   **Action:** As code is refactored, immediately update the corresponding Markdown files in `\docs`.
*   **Requirement:** The "Mock Management" (Process Area definitions) in the documentation must match the new Python class/function structure exactly.
    *   *Example:* If Python defines `create_media_area()`, the Markdown must describe "Area 200: Media Preparation".

### Step 2c: The Fixing-Loop (Verification)

*   **Action:** Run the new modular `LegHb_config1`.
*   **Logic:**
    1.  **Code Check:** Compare new TEA/LCA results with `baseline_metrics.txt`.
    2.  **Doc Check:** Verify that the logic described in `\docs` matches the actual code execution.
    3.  **If Error/Deviation:** Debug code OR fix documentation, then re-run.
    4.  **If Match:** Mark Step 2a as "Verified".

### Step 2d: Propagation

*   **Action:** Apply the Verified Modular Architecture (and Documentation updates) to `LegHb_config2` and `HemDx_config1`.
*   **Requirement:** Perform the **Fixing-Loop** for these configurations as well.

---

## GATEWAY: Pre-Phase 3 Check

**STOP and VERIFY:**

*   Are all configs (`LegHb`, `HemDx`) running without errors?
*   Do the TEA/LCA results match the Baseline?
*   **Is the Markdown documentation fully synced with the new Modular Scripts?**

*Only proceed to Phase 3 if all answers are YES.*

---

## Phase 3: "BioSTEAM-Mock" Skill Generation

**Objective:** Encapsulate the verified architecture into a reusable "Skill" for future projects.

1.  **Abstraction:** Create a `BioprocessFactory` (or `MockManager`) class that can dynamically assemble the verified Areas from Phase 2.
2.  **Intelligence:**
    *   The Skill must be able to recognize missing units.
    *   **Trigger:** If a future process requires a unit not in the library, the Agent will search Zotero/Web for design equations and generate a new `bst.Unit` subclass automatically.
3.  **Constraint:** The Skill must enforce the Area 100-900 ID standard automatically for any new process generated.

---

### Consultant's Note:

I have structured the **Gateway** to ensure that Phase 3 (which relies on the AI learning from your code) does not start until the "Training Data" (your scripts and markdown) are perfectly aligned. This will prevent the agent from "learning" the wrong architecture.

# BoilerTurbogenerator Unit Specification

**Class:** `BoilerTurbogenerator`
**Base Class:** `bst.BoilerTurbogenerator`
**Source:** `v1/_units.py`

## Purpose
Combined Heat and Power (CHP) unit that burns biomass residues and supplementary natural gas to generate steam and electricity. 

## Design Basis

1. **Efficiency**:
   - Boiler Efficiency ($\eta_B$): 0.80 (LHV basis).
   - Turbogenerator Efficiency ($\eta_{TG}$): 0.85 (Enthalpy-to-Work).
2. **Emissions Control**:
   - **SO2 Scrubbing**: Uses Lime ($Ca(OH)_2$ or $CaO$) with 20% stoichiometric excess if SO2 is present.
   - **NOx/Particulates**: Modeled via standard emission factors (derived from composition).
   - **Ammonia**: Routed to emissions (thermal destruction assumed or accounted in balance).
3. **Steam Balance**:
   - Satisfies process steam demand first.
   - Excess heat generates electricity.
   - If process steam > biomass heat, natural gas is supplemented.

## Logic Flow

1. **Fuel Evaluation**:
   - Calculates LHV/HHV of all feed solids and gas.
   - Determines if supplementary CH4 is needed.
2. **Combustion**:
   - complete combustion assumed.
   - Enthalpy of emissions updated with `Guarded` logic (robust to missing thermo params).
3. **Electricity**:
   - $Work = (H_{combustion} \cdot \eta_B - H_{steam,process}) \cdot \eta_{TG}$
   - Handles "Electricity Demand" satisfaction mode (optional).

## Cost Model

Standard NREL/BioSTEAM scaling:
- **Boiler**: Scaled on steam flow.
- **Turbogenerator**: Scaled on Work output (kW).
- **Baghouse/Scrubber**: Included in Boiler factor or separate (param dependent).

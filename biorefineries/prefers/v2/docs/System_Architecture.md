# System Architecture

## The System Factory Pattern
PREFERS uses BioSTEAM's `@SystemFactory` decorator to encapsulate process logic. This allows the system to be instantiated with flexible configurations (default config is used unless overridden).

### Construction Flow
1. **Chemicals Loaded**: `create_chemicals_[Module]` initializes thermodynamics.
2. **Streams Defined**: Static input streams are declared in each module's `_streams.py`.
3. **Process Definition**:
   - Units are instantiated from `biorefineries/prefers/v2/_units.py` and `_units_adv.py`.
   - Streams are connected (`ins` and `outs`).
   - **Specs**: Control logic is attached via `@unit.add_specification` (e.g., R302 internalized titer + NH3 convergence).
4. **System Creation**: `bst.System(...)` compiles the connectivity graph and registers the flowsheet.

## Internalized Specifications (v2)

The main fermenter (`R302`) embeds titer/NH3 convergence as a unit specification:
- Empirical NH3 demand measurement with excess-NH3 runs.
- Elasticity-based yield adjustment to hit target titer.
- Updates upstream NH3 source and S202 split for consistency.

This removes the need for external convergence utilities and keeps each config fully self-contained.

## TEA & LCA Engines

### TEA (Techno-Economic Analysis)
Economics are calculated by `PreFerSTEA` (inheriting from `bst.TEA`) in `biorefineries/prefers/v2/_tea.py`.
- **Inputs**: The `System` object (CAPEX/VOC and utilities).
- **Parameters**: IRR, depreciation schedule (IRAS/MACRS), tax rate, operating days.
- **Output**: Minimum Selling Price (MPSP) and standard CAPEX/FOC tables.

Project-specific subclasses live in:
- `biorefineries/prefers/v2/LegHb/_tea_config1.py`
- `biorefineries/prefers/v2/HemDx/_tea_config1.py` (and config2/3 equivalents)

### LCA (Life Cycle Assessment)
LCA factors are injected at stream level:
- `biorefineries/prefers/v2/_process_settings.py` defines `GWP_CFs`.
- `set_GWPCF` and `set_GWPCF_Multi` attach CFs to streams.
- BioSTEAM aggregates mass flow * CFs to compute GWP metrics.

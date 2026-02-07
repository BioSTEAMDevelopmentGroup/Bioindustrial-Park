# System Architecture

## The "System Factory" Pattern
PREFERS uses BioSTEAM's `@SystemFactory` decorator to encapsulate process logic. This allows the system to be instantiated with flexible configurations (though usually run with defaults).

### Construction Flow
1.  **Chemicals Loaded**: `create_chemicals_[Module]` is called to set up thermodynamics.
2.  **Streams Defined**: Static streams (inputs) are initialized.
3.  **Process Definition**:
    - Units are instantiated from `v1/units.py`.
    - Streams are connected (`ins` and `outs`).
    - **Specs**: Logic (like finding titer-based residence time) is added via `@unit.add_specification`.
4.  **System Creation**: `bst.System(...)` compiles the connectivity graph.

## TEA & LCA Engines

### TEA (Techno-Economic Analysis)
The economics are calculated by a custom `PreFerSTEA` class (inheriting from `bst.TEA`) found in `[Module]/tea.py`.
- **Inputs**: The `System` object (for CAPEX/VOC).
- **Parameters**: Interest rate, depreciation schedule (IRAS), tax rate.
- **Output**: Minimum Selling Price (MPSP).

### LCA (Life Cycle Assessment)
Links directly to the flowsheet streams.
- **Mechanism**: `process_settings.py` loads `GWP_CFs` (Global Warming Potential Characterization Factors).
- **Injection**: These factors are attached to stream objects. BioSTEAM then tracks mass flows * factors to calculate total GWP.

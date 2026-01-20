# LegHb Downstream Process Redesign Implementation Plan

## Goal

Redesign the downstream process section in `_config1.py` to implement a more realistic intracellular Leghemoglobin recovery process, replacing the current Ion Exchange-based purification with a food-grade streamlined approach.

---

## User Review Required

> [!IMPORTANT]
> **Upstream Process Unchanged**: All upstream units (M301–M304, R301–R302, T301–T302) remain exactly as-is. Any changes to upstream will require explicit user approval.

> [!WARNING]
> **Ion Exchange Removal**: The current design includes `ResinColumn2` (U402) for ion exchange chromatography. The new specification **removes** this step entirely, replacing it with:
> - Sequential Centrifugation + Depth Filtration for clarification
> - Single UF/DF stage for purification
> - HTST Pasteurization for thermal stabilization
>
> This is a significant simplification. Please confirm this is the intended design direction.

> [!CAUTION]
> **Potential Breaking Changes**: 
> - Input streams `IXEquilibriumBuffer`, `IXElutionBuffer`, `IXRegenerationSolution` will become unused
> - Output stream `RegenerationWaste` and associated neutralization (M502) may be removed
> - May affect TEA parameters and LCA inventory

---

## Current vs. New Process Flow Comparison

| Step | Current Design | New Design |
|------|----------------|------------|
| **1. From Fermentation** | R302 Broth → S401 | R302 Broth → Centrifuge (new) |
| **2. Biomass Harvest** | S401 (CellDisruption) | Disk Stack Centrifuge + Wash |
| **3. Cell Disruption** | *(disruption first)* | High-Pressure Homogenizer (after harvest) |
| **4. Lysate Clarification** | S402 (Filtration MF) → S403 (ScrewPress) | Centrifuge → Depth Filtration (series) |
| **5. Purification** | U401 (UF Diafiltration) → U402 (Ion Exchange) → U403 (NF Diafiltration) → U404 (UF concentration) | UF/DF only (combine U401+U404) |
| **6. Thermal Treatment** | None | HTST Pasteurization (new) |
| **7. Final Formulation** | H406 (cooling to 0°C) | Mixing + Rapid Chilling |

---

## Proposed Changes

### Component 1: Biomass Harvest (New Step Before Disruption)

> [!NOTE]
> Currently, S401 (CellDisruption) immediately processes fermentation broth. The new design separates cells first, then disrupts.

#### [MODIFY] Current S401 (CellDisruption) Position

**Current location**: Line 232-236  
**Change**: Move cell disruption AFTER biomass separation

**New unit sequence**:
```
R302 Broth
    ↓
S401_new: Disk Stack Centrifuge (Biomass Harvest)
    ├─→ Washed Cell Cream (to disruption)
    └─→ Spent Supernatant (to wastewater)
    ↓
S402_new: High-Pressure Homogenizer (relabel from old S401)
    ↓
(continue to clarification)
```

**Implementation**:
- **S401**: Replace with `bst.SolidsCentrifuge` or existing `Filtration.from_preset('MF')` for cell concentration
- Add cell wash step using buffer stream (1-2 diavolumes)
- Target: 40-50% WCW (Wet Cell Weight)

---

### Component 2: Cell Disruption (Keep, Re-position)

#### [KEEP] `_units.py` `CellDisruption` class

**Current**: Lines 232-236 in `_config1.py`  
**Change**: 
- Receives washed cell cream instead of raw broth
- Verify pressure parameters (800-1200 bar spec vs current 150 bar)
- Add post-valve cooling (HXutility) to return lysate to <15°C

```python
# Current:
S401 = u.CellDisruption('S401', ins=R302-1, outs='DisruptedBroth')

# New concept:
C401 = <CentrifugeSeparator>('C401', ins=R302-1, outs=('CellCream', 'SpentMedia'))
S401 = u.CellDisruption('S401', ins=C401-0, outs='DisruptedBroth', 
                         P_high=1000e5)  # 1000 bar = 100 MPa
H401_cool = bst.HXutility('H401_cool', ins=S401-0, T=15+273.15, cool_only=True)
```

---

### Component 3: Lysate Clarification (Replace MF + ScrewPress)

#### [MODIFY] S402, S403, S404

**Current**: 
- S402: `Filtration.from_preset('MF')` - Lines 238-244
- S403: `bst.ScrewPress` - Line 247  
- S404: `bst.Splitter` - Line 249

**New Design** (per spec):
1. High-Speed Centrifuge: Remove bulk solids (cell wall ghosts)
2. Depth Filtration: Series filtration (30 µm → 5 µm → 0.5 µm)

**Implementation Options**:
- Option A: Use existing `Filtration.from_preset('MF')` for centrifuge-equivalent, add depth filter as second `Filtration` unit
- Option B: Use `bst.SolidsCentrifuge` + series `Filtration` units

**Recommended** (minimal new units):
```python
# Replace S402-S404 with:
S402 = bst.SolidsCentrifuge('S402', ins=H401_cool-0, 
                            outs=('CellDebris', 'CrudeLysate'),
                            split={'cellmass': 0.95, 'Glucan': 0.90, ...})
S403 = u.Filtration.from_preset('MF', 'S403', ins=S402-1, 
                                 outs=('DebrisWaste', 'ClarifiedLysate'))
```

---

### Component 4: Ultrafiltration & Diafiltration (Simplify)

#### [MODIFY] U401, [DELETE] U402, [DELETE] U403, [MODIFY] U404

**Current flow**:
```
S404 → U401 (UF Diafiltration) → U402 (Ion Exchange) → U403 (NF Diafiltration) → U404 (UF concentration)
```

**New flow** (per spec - no Ion Exchange):
```
S403 → U401 (UF/DF combined) → U404 (concentration to final spec)
```

**Implementation**:
- **U401**: Keep `Diafiltration.from_preset('UF')`, update parameters:
  - MWCO: 3-10 kDa (current is appropriate for 16 kDa LegH)
  - VCF: 5X-10X concentration
  - Diafiltration: 5-7 diavolumes
  
- **U402, U403**: **DELETE** - Ion Exchange step removed
  - Remove associated buffer tanks (M402, M403, M404, H402, H403, H404)
  
- **U404**: Keep for final concentration, update target:
  - Target: 6-10% w/v protein content
  - Current spec already aims for 7.5% LegHb

---

### Component 5: Thermal Stabilization (NEW)

#### [NEW] HTST Pasteurizer

**Add after U404**, before final cooling:

```python
H407 = bst.HXutility('H407', ins=U404-0, T=72+273.15, heat_only=True)  # Heat to 72°C
# Hold time: Built into HX or add MixTank with tau=45/3600 (45 sec)
T401 = bst.MixTank('T401', ins=H407-0, outs='PasteurizedConcentrate', tau=45/3600)
```

**Alternative**: Use plate heat exchanger if available in BioSTEAM

---

### Component 6: Final Formulation & Cooling

#### [MODIFY] H406

**Current**: `bst.HXutility` cooling to 0°C (lines 501-507)

**New Design**:
- Add antioxidant mixing (ascorbate/erythorbate)
- Rapid chill to <4°C

```python
M501_form = bst.MixTank('M501_form', ins=(T401-0, 'AntioxidantStream'), 
                         outs='FormulatedProduct')
H406 = bst.HXutility('H406', ins=M501_form-0, outs=LegHb_3, 
                     T=4+273.15, cool_only=True)  # 4°C not 0°C
```

---

### Component 7: Water Recycle & Heat Integration

#### [MODIFY] Wastewater handling

**Current wastewater units**: M501, S501, S502, M502, S503

**New per spec**:
- UF/DF permeate → RO polishing → recycle to biomass wash or upstream
- Pasteurizer heat recovery (regenerative exchange)

**Implementation**:
- Keep existing RO units (S501, S503)
- Add heat integration between pasteurizer inlet/outlet (HXN or manual)

---

### Component 8: Facilities (Minimal Changes)

#### [KEEP] CT, CWP, BT, PWC

No changes to facilities expected. Update `makeup_water_streams` if new water streams added.

---

## File Modification Summary

| File | Scope | Details |
|------|-------|---------|
| `_config1.py` | **MAJOR** | Downstream section lines 228-546 |
| `_streams.py` | MINOR | May remove IX buffer streams, add antioxidant stream |
| `_chemicals.py` | MINOR | May need to add `SodiumErythorbate` if not using `SodiumAscorbate` |
| `_tea.py` | MINOR | Adjust if equipment list changes significantly |
| `_units.py` | NONE | Use existing units only |

---

## Verification Plan

### Automated Test
```bash
# Run from project root after modifications
cd c:\Programming\PreFerS\Bioindustrial-Park
python -c "from biorefineries.prefers.v1.LegHb.system._config1 import create_LegHb_system; sys = create_LegHb_system(); sys.simulate(); print('SUCCESS:', sys.flowsheet.stream.LegHb_3.F_mass, 'kg/hr')"
```

### Manual Verification
1. **Convergence Check**: System simulates without errors
2. **Mass Balance**: LegHb_3 output shows reasonable production rate (~275 kg/hr at target scale)
3. **Product Spec**: Run `check_LegHb_specifications(ss.LegHb_3)` - should pass all specs
4. **Diagram Review**: Generate `sys.diagram(format='html')` and visually verify flowsheet matches new design

---

## Questions for User

1. **Biomass Harvest Method**: Should we use `bst.SolidsCentrifuge` (disk stack) or adapt existing `Filtration.from_preset('MF')` (crossflow MF)?

2. **Ion Exchange Removal Confirmed?**: The new spec removes IEX entirely - this significantly simplifies the process but removes a purification step. Is this intentional?

3. **Pasteurization Parameters**: Spec says 70-75°C for 15-45 seconds. Use 72°C and 30 seconds as midpoint?

4. **Antioxidant Addition**: Current chemicals include `SodiumAscorbate`. Add mixing tank for this, or assume it's already in formulation buffer?

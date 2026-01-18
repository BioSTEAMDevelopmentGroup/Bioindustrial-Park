# VacuumPSA Specification

**Date:** 2026-01-18  
**Author:** PreFerS Development Team  
**Target Module:** `prefers/v1/_units.py`

---

## 1. Functional Description

**Purpose:** Separate gas mixture components (H₂, CO, C₂H₄) using vacuum pressure swing adsorption technology.

**Mechanism:** Cyclic adsorption process using zeolite 13X adsorbent. At high pressure (adsorption step), heavier gases (CO, C₂H₄) are preferentially adsorbed while H₂ passes through. At low/vacuum pressure (desorption step), adsorbed gases are released.

**Application:** Gas separation in syngas processing, particularly for hydrogen purification and olefin recovery.

---

## 2. Operating Conditions

| Parameter | Value | Unit | Source |
|:----------|:------|:-----|:-------|
| Adsorption Pressure | 6 × 10⁵ | Pa | Typical industrial PSA |
| Desorption Pressure | 1 × 10⁴ | Pa | VPSA vacuum regeneration |
| Temperature | 298.15 | K | Ambient operation |
| Cycle Time | 600 | s | 10-minute cycle typical |
| Number of Beds | 2 | - | Continuous operation |

---

## 3. Stream Configuration

### Inlets

| Index | Name | Description |
|:------|:-----|:------------|
| 0 | feed | Mixed gas feed (H₂, CO, C₂H₄, etc.) |

### Outlets

| Index | Name | Description |
|:------|:-----|:------------|
| 0 | product | H₂-rich product stream (raffinate) |
| 1 | purge | CO/C₂H₄-rich tail gas (extract) |

---

## 4. Mass Balance (`_run`)

### Separation Mechanism

The unit uses component-specific split factors based on adsorption selectivity:

| Component | Split to Product | Split to Purge | Basis |
|:----------|:-----------------|:---------------|:------|
| H₂ | 0.95 | 0.05 | Low adsorption affinity |
| CO | 0.30 | 0.70 | Moderate-high affinity |
| C₂H₄ | 0.20 | 0.80 | High affinity (π-bonding) |
| CO₂ | 0.10 | 0.90 | Strong quadrupole interaction |
| N₂ | 0.85 | 0.15 | Low affinity |
| CH₄ | 0.60 | 0.40 | Moderate affinity |

Default splits are based on zeolite 13X selectivity at 6 bar, 25°C from literature.

### Balance Equations

```python
for component in feed.chemicals.IDs:
    split = split_factors.get(component, 0.5)  # Default 50% if unknown
    product.imol[component] = feed.imol[component] * split
    purge.imol[component] = feed.imol[component] * (1 - split)
```

### Phase Behavior

- [x] Single phase (gas)
- [ ] VLE equilibrium
- [ ] LLE equilibrium
- [ ] Reactions involved

---

## 5. Sizing Logic (`_design`)

### Primary Sizing Parameters

| Parameter | Equation | Units |
|:----------|:---------|:------|
| Adsorbent Mass | m_ads = F_feed × τ_cycle / q_loading | kg |
| Bed Volume | V_bed = m_ads / ρ_bulk | m³ |
| Vacuum Power | P_vac = F_purge × ln(P_ads/P_des) × R × T / η | kW |

### Design Constants

| Constant | Value | Unit | Source |
|:---------|:------|:-----|:-------|
| Adsorbent loading (q) | 2.0 | mol/kg | Zeolite 13X for CO |
| Bulk density (ρ_bulk) | 650 | kg/m³ | Packed bed zeolite |
| Vacuum efficiency (η) | 0.70 | - | Typical vacuum pump |
| Compression ratio | 60 | - | P_ads/P_des = 6e5/1e4 |

### Utility Requirements

| Utility | Calculation | Notes |
|:--------|:------------|:------|
| Electricity | Vacuum pump + control | Primary utility |
| Cooling | Minor blower heating | Often negligible |

---

## 6. Cost Correlations (`_cost`)

### Equipment Cost Model

| Equipment | Base Cost ($) | Base Size | Exponent | CE Index | BM Factor |
|:----------|:--------------|:----------|:---------|:---------|:----------|
| Pressure Vessels | 50,000 | 10 m³ | 0.6 | 567 | 2.5 |
| Vacuum Pump | 30,000 | 100 kW | 0.7 | 567 | 1.8 |
| Adsorbent (initial) | 5 | $/kg | 1.0 | - | 1.0 |
| Piping & Valves | 15% of vessels | - | - | - | 1.0 |

### Cost Equations

```python
# Pressure vessels (2 beds)
C_vessels = 50000 * (V_bed / 10) ** 0.6 * N_beds

# Vacuum pump
C_vacuum = 30000 * (P_vac / 100) ** 0.7

# Adsorbent fill (consumable, treated as initial investment)
C_adsorbent = 5 * m_ads * N_beds

# Piping and valves
C_piping = 0.15 * C_vessels
```

### Additional Costs

| Item | Method | Notes |
|:-----|:-------|:------|
| Installation | BM factor | Included |
| Adsorbent replacement | Every 3-5 years | Operating cost |
| Instrumentation | 10% of equipment | Assumed included |

---

## 7. References

> [!IMPORTANT]
> **Source Quality Check:**
> - [x] Primary design data (Split factors): Ruthven (1994), Yang (1997)
> - [x] Cost data: Heuristic estimates scaled from standard pressure vessel costs (Peters & Timmerhaus)
> - [x] Application context: Naquash et al. (2021)

### Primary Sources (Design & Mechanism)
1. Ruthven, D.M., Farooq, S., Knaebel, K.S. *Pressure Swing Adsorption*, VCH Publishers, 1994. (Adsorption isotherms, cycle logic)
2. Yang, R.T. *Gas Separation by Adsorption Processes*, Imperial College Press, 1997. (Zeolite 13X selectivity)

### Secondary Sources (Application)
3. Naquash et al., "Hydrogen purification from syngas", *ICAE*, 2021. DOI: 10.46855/energy-proceedings. (Process configuration)

### Costing Basis
4. Peters, M.S., Timmerhaus, K.D., West, R.E. *Plant Design and Economics for Chemical Engineers*, 5th Ed., McGraw-Hill, 2003. (Pressure vessel and pump power laws)
5. NREL/TP-5100-47764 (Humbird et al. 2011). (Standard bare module factors for fluid processing equipment)

## 8. Implementation Notes

### Assumptions

- [x] Steady-state approximation of cyclic process
- [x] Ideal mixing in gas phase
- [x] Isothermal operation
- [x] Constant split factors (no pressure dependence)

### Limitations

- Valid for: P_ads = 3-10 bar, T = 15-40°C
- Not suitable for: Trace contaminant removal (<100 ppm)
- Model fidelity: Screening-level TEA/LCA (not detailed design)

### Future Improvements

- [ ] Pressure-dependent selectivity
- [ ] Temperature swing hybrid (VTSA)
- [ ] Multi-bed kinetic model
- [ ] Adsorbent degradation over time

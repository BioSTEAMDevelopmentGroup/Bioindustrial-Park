# Filtration Unit Specification

**Class:** `Filtration` (formerly `Filtration2`)
**Type:** Continuous Rotary Drum Vacuum Filter (RDVF)

> [!NOTE]
> Replaces generic Centrifuge models for rigorous biomass separation modeling.

---

## Design Basis

| Parameter | Default | Unit | Literature Range | Reference |
|:---|:---|:---|:---|:---|
| **Solids Loading** | 20.0 | kg/m²/hr | 10-50 kg/m²/hr | [1] Perry's Handbook (Table 18-5) |
| **Cake Moisture** | 20% | wt% | 20-80% | [2] Industrial benchmarks |
| **Capture Efficiency** | 99% | - | 95-99.9% | [3] Filtration Standards |
| **Power Intensity** | 1.0 | kW/m² | 0.5-2.0 kW/m² | [4] Perry's Handbook (Vacuum+Drive) |

---

## Key Equations

### 1. Sizing (Filter Area)
$$ A [m^2] = \frac{M_{solids} [kg/hr]}{L [kg/m^2/hr]} $$
- $M_{solids}$: Mass flow of dry solids in cake.
- $L$: Specific solids loading rate (20 kg/m²/hr default).

### 2. Utility (Power)
$$ P [kW] = A [m^2] \times P_{spec} [kW/m^2] $$
- Represents total power for drum drive and vacuum pump.

### 3. Costing
$$ C_p = 65,000 \times \left( \frac{A}{10} \right)^{0.6} $$
- Base cost ~$65k for 10 m² unit (Representative RDVF cost).

---

## References

1. **Perry's Chemical Engineers' Handbook, 8th Ed.** *Section 18: Liquid-Solid Operations.*
2. **Techno-economic comparison of dewatering technologies.**
3. **General Industrial Filtration Standards.**

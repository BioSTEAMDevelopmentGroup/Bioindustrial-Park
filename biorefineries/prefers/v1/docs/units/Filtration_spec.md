# Filtration Unit Specification

**Class:** `Filtration`
**Type:** Continuous Rotary Drum Vacuum Filter (RDVF)
**Last Updated:** 2026-01-20

> [!NOTE]
> Replaces generic Centrifuge models for rigorous biomass separation modeling.
> Parameters validated for yeast/mycelia filtration.

---

## Design Basis

| Parameter | Default | Unit | Literature Range | Reference |
|:---|:---|:---|:---|:---|
| **Solids Loading** | 20.0 | kg/m²/hr | 10-50 kg/m²/hr | [1] Perry's Handbook (Table 18-5) |
| **Cake Moisture** | 20% | wt% | 20-80% | [2] Industrial benchmarks |
| **Capture Efficiency** | 99% | - | 95-99.9% | [3] Filtration Standards |
| **Power Intensity** | 1.0 | kW/m² | 0.5-2.0 kW/m² | [1] Perry's Handbook |
| **TMP** | 2.0 | bar | 1.5-3.0 bar | RDVF typical |
| **Membrane Cost** | 100.0 | USD/m² | $50-$200/m² | Filter cloth |
| **Membrane Lifetime** | 3.0 | Years | 2-5 years | Industrial |

### Presets

| Preset | Solids Loading | TMP (bar) | Power (kW/m²) | Application |
|:---|:---:|:---:|:---:|:---|
| **MF** | 30 | 1.5 | 0.8 | Cell separation, primary clarification |
| **UF** | 15 | 3.0 | 1.2 | Finer separation, virus removal |

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

## Application Notes for LegHb Production

For P. pastoris / K. marxianus biomass:
- RDVF suitable for primary clarification [1]
- Cell disruption creates smaller debris requiring downstream polishing
- Typical solids loading 10-30 kg/m²/hr for yeast [4]
- Consider centrifugation for high-density cell harvest

---

## References

1. **Perry's Chemical Engineers' Handbook, 8th Ed.** *Section 18: Liquid-Solid Operations.*
2. **Techno-economic comparison of dewatering technologies.**
3. **General Industrial Filtration Standards.**
4. **Yao et al.** (2025). "Recent advances in microbial fermentation for the production of leghemoglobin." *World J Microbiol Biotechnol* 41:404.

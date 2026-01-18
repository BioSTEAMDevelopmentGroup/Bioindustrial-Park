# Filtration2 (Rotary Vacuum Filter)

**Class:** `Filtration2`
**Type:** Continuous Rotary Drum Vacuum Filter (RDVF)

## Description
A continuous filtration unit operation modeling a Rotary Drum Vacuum Filter. It is designed for high-volume solid-liquid separation in biorefineries, particularly for separating biomass (cells) from fermentation broth. The model allows for specifying cake dryness (moisture content) and solid capture efficiency.

## Design Basis

| Parameter | Value | Units | Source |
|:---|:---|:---|:---|
| **Specific Loading** | 20.0 | kg-solids/m²/hr | [1] Perry's Handbook, Table 18-5 (Biomass/Mycelia range 10-50) |
| **Cake Moisture** | 20-50% | wt% | [2] Techno-comparison of dewatering (RVF typical 50-80% moisture, optimized 20%) |
| **Capture Eff.** | 99.0% | - | [3] Typical industrial filtration target |
| **Power Demand** | 1.0 | kW/m² | [4] Perry's Handbook (Vacuum pump + drive) |
| **Filter Cost** | - | USD | [5] NREL Template / Matches CEPCI scaling |

## Equations

### 1. Sizing
Filter Area ($A$) is determined by the solids throughput ($M_s$) and specific loading rate ($L$):
$$ A [m^2] = \frac{M_s [kg/hr]}{L [kg/m^2/hr]} $$

### 2. Utility
Power is estimated based on filter area:
$$ P [kW] = A [m^2] \times 1.0 [kW/m^2] $$

### 3. Costing
Purchase cost ($C_p$) is calculated using a size-based power law:
$$ C_p = C_{base} \times \left( \frac{A}{A_{base}} \right)^{0.6} $$

## References
[1] Perry's Chemical Engineers' Handbook, 8th Ed.
[2] "Techno-economic comparison of dewatering technologies."
[3] Industrial filtration standards.
[4] Typical drive and vacuum power for RDVF.
[5] NREL Process Design and Economics Guidelines.

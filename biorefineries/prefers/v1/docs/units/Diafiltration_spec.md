# Diafiltration Unit Specification

**Class:** `Diafiltration` (formerly `Diafiltration2`)
**Type:** Membrane Separation (Ultrafiltration/Nanofiltration) with Buffer Exchange

> [!NOTE]
> Validated against industrial literature for protein purification applications.

---

## Design Basis

| Parameter | Default | Unit | Literature Range | Reference |
|:---|:---|:---|:---|:---|
| **Membrane Flux** | 40.0 | L/m²/hr | 20-80 LMH (UF) | [1] Membranes.com |
| **Membrane Cost** | 150.00 | USD/m² | $80-$350/m² | [2] SNS Insider |
| **Membrane Lifetime** | 2.0 | Years | 1-3 Years | [3] DuPont |
| **Product Retention** | 99% | - | - | User Specified |
| **Op. Pressure** | 3.0 | bar | 2-10 bar | [5] Synder Filtration |

---

## Key Equations

### 1. Mass Balance (Solute Retention)
$$ M_{permeate,i} = M_{in,i} \times (1 - R_i) $$
$$ M_{retentate,i} = M_{in,i} \times R_i $$
- $R_i$: Retention coefficient for component $i$.

### 2. Sizing (Membrane Area)
$$ A [m^2] = \frac{Q_p [L/hr]}{J [L/m^2/hr]} $$
- $Q_p$: Permeate flow rate.
- $J$: Membrane Flux (40 LMH default).

### 3. Costing
$$ C_{sys} = C_{base} \times \left( \frac{A}{100} \right)^{0.7} $$
- $C_{base}$: ~$25,000 for 100 m² module (illustrative module factor).

---

## References

1. **Membranes.com.** *Ultrafiltration Membrane Flux Rates.*
2. **SNS Insider.** *Membrane Market Analysis 2023.*
3. **DuPont Water Solutions.** *Membrane Life and Maintenance.*
4. **Humbird et al.** *NREL Technical Report TP-5100-47764.* (2011)
5. **Synder Filtration.** *UF Operating Parameters.*

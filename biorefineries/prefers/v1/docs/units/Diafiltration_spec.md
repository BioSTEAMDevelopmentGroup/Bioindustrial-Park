# Diafiltration Unit Specification

**Class:** `Diafiltration`
**Type:** Membrane Separation (Ultrafiltration/Nanofiltration) with Buffer Exchange
**Last Updated:** 2026-01-20

> [!NOTE]
> Validated against industrial literature for protein purification applications.
> Parameters verified for LegHb production per Yao et al. (2025).

---

## Design Basis

| Parameter | Default | Unit | Literature Range | Reference |
|:---|:---|:---|:---|:---|
| **Membrane Flux** | 40.0 | L/m²/hr | 20-80 LMH (UF), 10-40 LMH (NF) | [1] Membranes.com |
| **Membrane Cost** | 150.00 | USD/m² | $80-$350/m² (industrial) | [2] SNS Insider 2023 |
| **Membrane Cost (Pharma)** | 1500-5000 | USD/m² | $1500-$5000/m² (TFF cassettes) | [2] SNS Insider 2023 |
| **Membrane Lifetime** | 2.0 | Years | 1-3 Years | [3] DuPont Water Solutions |
| **Product Retention** | 99% | - | >95% typical for UF | User/Application Specific |
| **Op. Pressure (TMP)** | 2.0-4.0 | bar | 2-10 bar (UF), 10-25 bar (NF) | [5] Synder Filtration |
| **Recirculation Ratio** | 10.0 | - | 5-20 | Fouling control |
| **Equipment Lifetime** | 20.0 | Years | - | [4] NREL Guidelines |

### Presets

| Preset | Flux (LMH) | TMP (bar) | Cost ($/m²) | Lifetime (yr) | Application |
|:---|:---:|:---:|:---:|:---:|:---|
| **UF** | 50 | 2.0/1.5 | 150 | 3 | Protein concentration |
| **NF** | 25 | 8.0/5.0 | 250 | 2 | Buffer exchange, fine separation |

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

## Application Notes for LegHb Purification

For Leghemoglobin (~16 kDa monomer) purification:
- Use 3-10 kDa MWCO UF membranes for high retention
- 5-7 diavolumes recommended for buffer exchange [6]
- TMP 2-4 bar to avoid protein denaturation
- Temperature <10°C to preserve activity

---

## References

1. **Membranes.com.** *Ultrafiltration Membrane Flux Rates.*
2. **SNS Insider.** *Membrane Market Analysis 2023.*
3. **DuPont Water Solutions.** *Membrane Life and Maintenance.*
4. **Humbird et al.** *NREL Technical Report TP-5100-47764.* (2011)
5. **Synder Filtration.** *UF Operating Parameters.*
6. **Yao et al.** (2025). "Recent advances in microbial fermentation for the production of leghemoglobin." *World J Microbiol Biotechnol* 41:404.

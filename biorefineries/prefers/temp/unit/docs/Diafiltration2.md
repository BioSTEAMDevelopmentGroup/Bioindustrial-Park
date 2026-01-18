# Diafiltration2

**Class:** `Diafiltration2`
**Type:** Membrane Separation (Ultrafiltration/Nanofiltration) with Buffer Exchange

## Description
A rigorous Diafiltration unit operation designed for protein purification and buffer exchange. It models the separation of molecules based on size (retention) while continuously adding wash buffer to remove permeable impurities (salts, small molecules). The model includes cost estimation for membrane modules, pump power, and replacement costs based on industrial literature.

## Design Basis

| Parameter | Value | Units | Source |
|:---|:---|:---|:---|
| **Flux** | 40.0 | L/m²/hr | [1] Membranes.com (Typical UF range 20-80 LMH) |
| **Membrane Cost** | 150.00 | USD/m² | [2] SNS Insider, Samcotech (Range \$80-\$350/m²) |
| **Membrane Lifetime** | 2.0 | Years | [3] DuPont, STERIS (Typical 1-3 years for polymeric) |
| **Equipment Life** | 20.0 | Years | [4] NREL TEA Guidelines |
| **Op. Pressure (P1)** | 3.0 | bar | [5] Synder Filtration (UF typical 2-10 bar) |
| **Recirculation** | 10x | - | [5] Engineering heuristic for cross-flow velocity |

## Equations

### 1. Mass Balance
The unit employs a rigorous solute flux model:
$$ M_{permeate,i} = M_{in,i} \times (1 - R_i) $$
$$ M_{retentate,i} = M_{in,i} \times R_i $$
Where $R_i$ is the retention coefficient for component $i$.

### 2. Sizing
Membrane area ($A$) is calculated from the permeate volumetric flow ($Q_p$) and flux ($J$):
$$ A [m^2] = \frac{Q_p [L/hr]}{J [L/m^2/hr]} $$

### 3. Costing
Capital cost ($C_{sys}$) is estimated using a power law:
$$ C_{sys} = C_{base} \times \left( \frac{A}{A_{base}} \right)^{n} $$
Where $C_{base} = \$25,000$ (module factor) and $n = 0.7$.

Membrane replacement cost is annualized or calculated as total NPV over plant life.

## References
[1] Membranes.com. "Ultrafiltration Membrane Flux Rates."
[2] SNS Insider. "Membrane Market Analysis 2023."
[3] DuPont Water Solutions. "Membrane Life and Maintenance."
[4] Humbird et al. (2011). NREL Technical Report TP-5100-47764.
[5] Synder Filtration. "UF Operating Parameters."

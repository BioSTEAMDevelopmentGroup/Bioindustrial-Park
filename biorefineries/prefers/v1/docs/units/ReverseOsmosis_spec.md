# ReverseOsmosis (RO2) Unit Documentation

Reverse Osmosis unit with literature-validated parameters for biorefinery wastewater treatment.

> [!NOTE]
> As of v1.1, `ReverseOsmosis` is an alias to `RO2`. The legacy implementation is available as `ReverseOsmosis_Legacy` but is deprecated.

---

## Design Basis

| Parameter | Default | Unit | Literature Range | Reference |
|:----------|:--------|:-----|:-----------------|:----------|
| Water Recovery | 0.85 (85%) | - | 75-95% | Axeon Water Technologies |
| Membrane Flux | 40 | L/m²/hr | 17-40 LMH | membranes.com |
| Membrane Cost | 200 | USD/m² | $80-$350/m² | SNS Insider 2023 |
| Membrane Lifetime | 3 | years | 3-5 years | DuPont/DOW |
| Plant Lifetime | 20 | years | 15-25 years | NREL TEA |
| Operating Pressure | 25 | bar | 15-30 bar | Axeon Water |
| Specific Energy | 3.0 | kWh/m³ | 2-6 kWh/m³ | VSEP, DuPont |

---

## Key Equations

### Membrane Area
```
A_membrane = Q_permeate / J_flux
```
Where:
- `A_membrane` = Required membrane area (m²)
- `Q_permeate` = Permeate flow rate (L/hr) = Feed × Recovery × 1000
- `J_flux` = Membrane flux (L/m²/hr)

### Power Consumption
```
P = max(P_SEC, P_pump)
```
Where:
- `P_SEC` = Q_permeate × SEC (specific energy consumption)
- `P_pump` = Pump power from BioSTEAM Pump unit

### Costing
```
C_pump = f(P) from BioSTEAM Pump correlation
C_membrane = A × Cost_per_m² × (Plant_life / Membrane_life)
C_hardware = 50000 × (A / 100)^0.7
```

---

## References

1. **Axeon Water Technologies.** *Understanding RO Water Recovery Rates.* [axeonwater.com](https://www.axeonwater.com)

2. **DuPont Water Solutions.** *FilmTec Reverse Osmosis Technical Manual.* (2023)

3. **SNS Insider.** *Reverse Osmosis Membrane Market Report.* (2023) - Industrial membrane pricing data.

4. **NREL.** *Process Design and Economics for Biochemical Conversion of Lignocellulosic Biomass to Ethanol.* Technical Report NREL/TP-5100-47764. (2011)

5. **Peters, M.S. & Timmerhaus, K.D.** *Plant Design and Economics for Chemical Engineers.* 5th Ed. McGraw-Hill.

6. **membranes.com.** *Industrial RO Membrane Specifications.* - Flux ranges for wastewater applications.

---

## Usage Example

```python
import biosteam as bst
from biorefineries.prefers.v1._units import ReverseOsmosis  # Now aliases to RO2

# Create feed stream
feed = bst.Stream('RO_feed', H2O=1000, NaCl=10, units='kg/hr')

# Create RO unit with default validated parameters
RO = ReverseOsmosis('RO1', ins=feed)

# Or with custom parameters for high-recovery application
RO_high = ReverseOsmosis('RO2', ins=feed, 
              water_recovery=0.90,
              operating_pressure_bar=35)

# Simulate
RO.simulate()
print(RO.design_results)
print(RO.purchase_costs)
```

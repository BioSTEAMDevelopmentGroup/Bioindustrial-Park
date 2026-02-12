# FiltrationAdv Specification

**Advanced Filtration Unit** with capacity-based sizing and bidirectional parameter modes.

**Module**: `biorefineries.prefers.v2._units_adv`

## Overview
`FiltrationAdv` implements a dead-end filtration model with capacity-based area sizing:
- **Backward** (default): User sets capture efficiency, area calculated from capacity
- **Forward**: User sets filter area, efficiency determined by loading

## Mechanistic Model

### Dead-End Filtration Sizing
```
Area_needed = V_batch / (V_max × η_safety)
```

Where:
- **V_max** = Maximum capacity before filter clogs (L/m²)
- **η_safety** = Safety factor (default 0.8)

### Solids Separation
- Solid capture efficiency applied to identified solids
- Solubles follow water with 1% entrainment
- Cake moisture controlled by mass balance

## Parameters

| Parameter                  | Default      | Description                    |
| -------------------------- | ------------ | ------------------------------ |
| `parameter_mode`           | `'backward'` | `'backward'` or `'forward'`    |
| `solid_capture_efficiency` | 0.98         | Fraction of solids captured    |
| `cake_moisture_content`    | 0.25         | Moisture in cake (wt fraction) |
| `V_max_L_per_m2`           | 200.0        | Max capacity (L/m²)            |
| `safety_factor`            | 0.8          | η_safety                       |
| `TMP_bar`                  | 1.5          | Transmembrane pressure (bar)   |
| `TMP_ref`                  | 1.5          | Reference TMP for capacity scaling |
| `power_per_m2`             | 0.8          | Specific power (kW/m²)         |
| `efficiency`               | 0.95         | Global efficiency scalar       |

## Presets

| Preset | V_max (L/m²) | TMP (bar) | Capture Eff. | Description                        |
| ------ | ------------ | --------- | ------------ | ---------------------------------- |
| `'MF'` | 200.0        | 1.5       | 0.98         | Microfiltration (cell separation)  |
| `'UF'` | 150.0        | 3.0       | 0.995        | Ultrafiltration (finer separation) |

## Usage Example
```python
from biorefineries.prefers.v2._units_adv import FiltrationAdv

S404 = FiltrationAdv.from_preset(
    'MF', 'S404',
    ins=supernatant,
    outs=('FilterCake', 'Filtrate'),
    solid_capture_efficiency=0.95,
    parameter_mode='backward',
)
```

## Reference
- User model: `filtration_and_diafiltation_adv.md`
- Perry's Handbook 8th Ed, Table 18-5 (RDVF parameters)

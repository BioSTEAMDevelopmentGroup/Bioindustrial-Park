# DiafiltrationAdv Specification

**Advanced Diafiltration Unit** with mechanistic modeling and bidirectional parameter modes.

**Module**: `biorefineries.prefers.v1._units_adv`

## Overview
`DiafiltrationAdv` implements a **Constant Volume Diafiltration (CVD)** model assuming **Perfect Mixing** for diafiltration with two operating modes:
- **Backward** (default): User sets retention targets directly
- **Forward**: User sets diavolumes, retention is calculated via washout equation

## Mechanistic Model

### Washout Equation
After N diavolumes:
```
C_out / C_in = exp(-k × N_DV × (1 - σ))
```

Where:
- **N_DV** = diavolumes (V_buffer / V_retentate)
- **σ** = rejection coefficient (0 = free pass, 1 = full retention)
- **k** = calibration constant (0.92, calibrated for 99% removal at 5 DV for σ=0)

### Backward Calculation
To achieve target retention:
```
N_DV = -ln(Retention) / (k × (1 - σ))
```

## Parameters

| Parameter                 | Default      | Description                       |
| ------------------------- | ------------ | --------------------------------- |
| `parameter_mode`          | `'backward'` | `'backward'` or `'forward'`       |
| `TargetProduct_Retention` | 0.99         | Fraction of product retained      |
| `Salt_Retention`          | 0.05         | Fraction of salts retained        |
| `diavolumes`              | 5.0          | Number of diavolumes              |
| `sigma_target`            | 0.99         | Rejection coefficient for product |
| `sigma_salt`              | 0.05         | Rejection coefficient for salts   |
| `membrane_flux_LMH`       | 50.0         | Membrane flux (L/m²/hr)           |
| `efficiency`              | 0.95         | Global efficiency scalar          |

## Presets

| Preset | Flux (LMH) | TMP (bar) | Description     |
| ------ | ---------- | --------- | --------------- |
| `'UF'` | 50.0       | 2.0       | Ultrafiltration |
| `'NF'` | 25.0       | 8.0       | Nanofiltration  |

## Usage Example
```python
from biorefineries.prefers.v1._units_adv import DiafiltrationAdv

U601 = DiafiltrationAdv.from_preset(
    'NF', 'U601',
    ins=(eluate_in, buffer),
    outs=('Retentate', 'Permeate'),
    TargetProduct_Retention=0.98,
    Salt_Retention=0.10,
    parameter_mode='backward',
)
```

## Reference
- User model: `filtration_and_diafiltation_adv.md`
- Standard TFF theory: Constant Volume Diafiltration (CVD) with Perfect Mixing

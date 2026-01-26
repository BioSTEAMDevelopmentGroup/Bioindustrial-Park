# N-HemDx DSP (Downstream Processing) Design Draft

## Overview

This document summarizes the **preliminary DSP design** for N-HemDx (Nicotinamide-Stabilized Hemodextrin) production, implemented in `_config1_DSP_Draft.py`.

## Process Topology: Split-Stream Heme Recovery

```
Fermentation Broth
        │
        ▼
  ┌─────────────┐
  │ S101 Cent.  │ Area 100: Clarification
  └──────┬──────┘
         │
    ┌────┴─────┐
    │          │
    ▼          ▼
Supernatant  Wet Pellet
 (Path A)    (Path B)
    │          │
    │     ┌────┴────┐
    │     │ T201    │ Area 200: Recovery
    │     │ U202 HPH│
    │     │ S203    │
    │     │ S204    │
    │     └────┬────┘
    │          │
    └────┬─────┘
         │
         ▼
  ┌─────────────┐
  │ M301 Merge  │ Area 300: Capture
  │ U302 Resin  │ (IonExchange)
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │ U401 NF/DF  │ Area 400: Concentration
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │ M501 γ-CD   │ Area 500: Formulation
  │ R502 Cplex  │ Heme + γ-CD → HemoDextrin
  │ M503 Nic    │ HemoDextrin + Nic → N-HemDx
  │ R504 Stab   │
  └──────┬──────┘
         │
         ▼
  ┌─────────────┐
  │ S601 Cool   │ Area 600: Product
  └──────┬──────┘
         │
         ▼
   N-HemDx Product
```

## Key Process Areas

### Area 100: Clarification
- **S101**: Centrifuge separating broth into supernatant (extracellular heme) and pellet (intracellular)
- Split ratio configurable via `extracellular_fraction` parameter

### Area 200: Intracellular Recovery
- **T201**: Resuspension tank with buffer
- **U202**: High-pressure homogenizer (1000 bar, 55% disruption)
- **S203**: Debris centrifuge
- **S204**: Guard filter (MF)

### Area 300: Capture & Purification
- **M301**: Stream merger (Path A + Path B)
- **U302**: ResinColumn (IonExchange preset)
  - Target: Heme_b, ProtoporphyrinIX (90% yield)
  - Remove: Salts, Glucose, Glycine

### Area 400: Concentration
- **U401**: Diafiltration (NF preset)
  - Heme retention: 98%
  - Salt removal: 90%
  - 5× concentration factor

### Area 500: Formulation
- **M501 + R502**: γ-Cyclodextrin addition and complexation
  - `Heme_b + GammaCyclodextrin → HemoDextrin` (95%)
- **M503 + R504**: Nicotinamide addition and stabilization
  - `HemoDextrin + Nicotinamide → N-HemoDextrin` (90%)

### Area 600: Final Product
- **S601**: Cold storage conditioning (4°C)

## Chemicals Added

| Chemical       | MW (g/mol) | Role                    |
| -------------- | ---------- | ----------------------- |
| γ-Cyclodextrin | 1297.12    | Carrier molecule        |
| Nicotinamide   | 122.12     | Axial ligand stabilizer |
| HemoDextrin    | ~1913      | Intermediate complex    |
| N-HemoDextrin  | ~2035      | Final product           |

## Verification Status

✅ Simulation completed successfully
- N-HemoDextrin: 0.077 kg/hr
- HemoDextrin: 0.009 kg/hr
- Free Heme: 0.110 kg/hr

## Integration Notes

This DSP module needs to be integrated with:
1. Upstream media preparation (Area 200 in full system)
2. Fermentation (Area 300 in full system)

The fermentation broth output should connect to S101 input for the split-stream processing.

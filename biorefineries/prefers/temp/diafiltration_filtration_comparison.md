# Diafiltration & Filtration: Old vs New Comparison

## 1. Diafiltration

| Feature | Old `Diafiltration` | New `Diafiltration` (Upgraded) | Improvement |
|:---|:---|:---|:---|
| **Flux** | Not explicitly validated (Default implicit) | **40 LMH** (Validated) | Based on industrial UF data (Membranes.com, Synder) |
| **Membrane Cost** | Generic defaults | **$150 / m²** | Updated 2023 market pricing (SNS Insider) |
| **Replacement** | Generic defaults | **2 years** | Aligned with DuPont technical manuals |
| **Solute Logic** | Mixed logic (if/else blocks) | Unified `retention_map` | Cleaner, more robust mass balance |
| **Validation** | None provided | **Full Citations** in docstring | Traceable design basis |

## 2. Filtration (Replacing Centrifuge)

| Feature | Old `Centrifuge` (Standard) | New `Filtration` (Rotary Vacuum) | Improvement |
|:---|:---|:---|:---|
| **Type** | Black-box Splitter | **Physics-based Filtration** | Models surface area and cake formation |
| **Sizing Basis** | Throughput (Flow Rate) | **Solids Loading (kg/m²/hr)** | Explicit scaling with solid content (Perry's Handbook) |
| **Capture Eff.** | Fixed Split (0.995) | **99% (Configurable)** | Consistent with RDVF performance |
| **Power** | Generic correlation | **1.0 kW/m²** | Physics-based power demand (Vacuum + Drive) |
| **Costing** | Generic Centrifuge Curve | **RDVF Cost Curve** | More accurate for large-scale biomass separation |

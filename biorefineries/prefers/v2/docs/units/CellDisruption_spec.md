# CellDisruption Unit Specification

**Class:** `CellDisruption`
**Type:** High-Pressure Homogenizer (HPH)
**Source:** `v1/_units.py`

## Purpose
Releases intracellular products (LegHb, Heme) by lysing yeast cells under high pressure.

## Process Physics

1. **Mass Balance**:
   - Converts `imass[Cell_ID]` to component fractions (Mannoprotein, Glucan, RNA, etc.) based on `cell_disruption_efficiency`.
   - Releases intracellular species (`_In` suffix) to extracellular phase.
   - Default Efficiency: 0.55 (0.90 typically used in config).
2. **Pressure Dynamics**:
   - Simulates multi-stage pumping to reach `P_high` (e.g., 1000 bar).
   - Limits per-stage compression ratio (<4.0) and head (<3000 ft) to conservatively model realistic pump constraints and prevent BioSTEAM warnings.
3. **Thermal Effect**:
   - Modeled via `bst.IsenthalpicValve` from `P_high` to `P_low`.
   - Captures temperature rise due to energy dissipation (Joule-Thomson effect).

## Cost Model

**Correlation**:
$$ Cost = 90,000 \cdot (Q \cdot 1000)^{0.5} \cdot (dP/1000)^{1.5} $$

- **Q**: Flow rate (m³/hr)
- **dP**: Pressure drop (bar)
- **Ref**: Inguva et al. (2024), Chem Eng Res Des.

**Power**:
$$ P_{kW} = \frac{Q_{L/hr} \cdot dP_{bar}}{36000 \cdot \eta} $$
- $\eta$: 0.85 (Efficiency)

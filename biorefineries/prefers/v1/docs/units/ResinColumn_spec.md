# ResinColumn Specification

**Class:** `ResinColumn`
**Type:** Polymorphic Separation Unit
**Presets:** 'IonExchange', 'Adsorption'
**Last Updated:** 2026-01-20

> [!NOTE]
> Parameters validated against literature. See Walch & Jungbauer (2017) for continuous chromatography design principles.

## Overview
The `ResinColumn2` unit is a versatile packed-bed column model that can simulate two distinct physical processes based on the selected `preset`:

1.  **Ion Exchange (Chromatography):** Models a bind-and-elute cycle where a target molecule binds to the resin and is subsequently eluted. Key design parameter is Dynamic Binding Capacity (DBC).
2.  **Adsorption (e.g., Activated Carbon):** Models a flow-through or capture process where contaminants are removed from the stream. Key design parameter is Empty Bed Contact Time (EBCT).

## 1. Ion Exchange Preset

**Process:** Bind-Elute
**Key Physics:** Resin volume is calculated to bind the total mass of target product per cycle, subject to a safety factor on the Dynamic Binding Capacity (DBC).

### Validated Parameters

| Parameter              | Default | Unit | Source             | Notes                                          |
| :--------------------- | :------ | :--- | :----------------- | :--------------------------------------------- |
| `resin_DBC_g_L`        | 50.0    | g/L  | Typical Industrial | Range 30-100 g/L depending on resin/target     |
| `cycle_time_hr`        | 4.0     | hr   | Operation Schedule | Load + Wash + Elute + Regen                    |
| `equilibration_CV`     | 5.0     | CV   | Heuristic          |                                                |
| `wash_CV`              | 5.0     | CV   | Heuristic          |                                                |
| `elution_CV`           | 3.0     | CV   | Heuristic          |                                                |
| `regeneration_CV`      | 5.0     | CV   | Heuristic          |                                                |
| `resin_lifetime_years` | 5.0     | yr   | Vendor Data        | 1-10 years typical                             |
| `resin_cost_USD_per_L` | 30.0    | $/L  | Vendor Quotes      | Strong Cation: \$10-50, Strong Anion: \$30-100 |

### Equations
$$ V_{resin} = \frac{M_{target, cycle}}{DBC \times SafetyFactor} $$
$$ Cost_{resin} = V_{resin} \times Cost_{per\_L} $$

---

## 2. Adsorption Preset

**Process:** Flow-Through / Capture
**Key Physics:** Bed volume is calculated based on the required Empty Bed Contact Time (EBCT) to achieve separation. Column dimensions are constrained by superficial velocity.

### Validated Parameters

| Parameter                   | Default | Unit  | Source               | Notes                                   |
| :-------------------------- | :------ | :---- | :------------------- | :-------------------------------------- |
| `EBCT_min`                  | 5.0     | min   | urbansaqua.com       | 2-10 min typical for liquids (Cl, VOCs) |
| `superficial_velocity_m_h`  | 10.0    | m/h   | aquaenergyexpo.com   | 5-20 m/h typical                        |
| `adsorbent_bulk_density`    | 450     | kg/m³ | calgoncarbon.com     | GAC: 350-550 kg/m³                      |
| `adsorbent_cost_USD_per_kg` | 5.0     | $/kg  | Market Data          | Activated Carbon: ~$2-8/kg              |
| `adsorbent_lifetime_years`  | 3.0     | yr    | Operation Experience | Depends on regeneration frequency       |

### Equations
$$ V_{bed} = Q_{vol} \times EBCT $$
$$ A_{min} = \frac{Q_{vol}}{v_{superficial}} $$
$$ Mass_{adsorbent} = V_{bed} \times \rho_{bulk} $$

## Usage Examples

```python
# Ion Exchange Mode
R201 = ResinColumn2('R201', ins=[feed, buffer_A, buffer_B, regen], 
                    preset='IonExchange',
                    resin_DBC_g_L=60)

# Adsorption Mode (e.g. Carbon Filter)
AC301 = ResinColumn2('AC301', ins=[feed, None, None, regen_waste_dest],
                     preset='Adsorption',
                     EBCT_min=10.0)
```

---

## LegHb Application Notes

For leghemoglobin production, the IonExchange preset is typically used for final polishing:
- **Target binding:** LegHb at pI ~4.5 binds well to strong anion resins above pH 6.0
- **Recommended DBC:** 40-60 g/L for recombinant LegHb (lower than native proteins due to size/conformation)
- **Elution:** Salt gradient (0.1-0.5 M NaCl) or pH step

Adsorption mode used for decolorization/off-flavor removal in final stages.

## References

1. Walch, N., & Jungbauer, A. (2017). *Continuous chromatography for protein purification: A review.* Journal of Chromatography A, 1498, 1-11.
2. Thömmes, J., & Etzel, M. (2007). *Alternatives to chromatographic separations.* Biotechnology Progress, 23(1), 42-45.
3. Vendor documentation: GE Healthcare (Cytiva), Bio-Rad, Purolite ion exchange resins.
4. Calgon Carbon Corporation (2023). *Granular Activated Carbon Technical Specifications.*

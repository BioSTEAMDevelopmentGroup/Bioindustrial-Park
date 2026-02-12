# Unit Operations Catalog

**Source:** `biorefineries/prefers/v2/_units.py`, `biorefineries/prefers/v2/_units_adv.py`

This catalog lists the custom unit operations used in PREFERS v2. Detailed design and costing are documented in the unit specs under `docs/units/`.

## Catalog (v2)

| Unit | Class | Module | Detail |
| :--- | :---- | :----- | :----- |
| Seed Train | `SeedTrain` | `_units.py` | [SeedTrain_spec.md](./units/SeedTrain_spec.md) |
| Fermenter | `AeratedFermentation` | `_units.py` | [AeratedFermentation_spec.md](./units/AeratedFermentation_spec.md) |
| Cell Disruption | `CellDisruption` | `_units.py` | [CellDisruption_spec.md](./units/CellDisruption_spec.md) |
| Centrifuge | `Centrifuge` | `_units.py` | [Centrifuge_spec.md](./units/Centrifuge_spec.md) |
| Diafiltration | `Diafiltration` | `_units.py` | [Diafiltration_spec.md](./units/Diafiltration_spec.md) |
| Diafiltration (Adv) | `DiafiltrationAdv` | `_units_adv.py` | [DiafiltrationAdv_spec.md](./units/DiafiltrationAdv_spec.md) |
| Filtration | `Filtration` | `_units.py` | [Filtration_spec.md](./units/Filtration_spec.md) |
| Filtration (Adv) | `FiltrationAdv` | `_units_adv.py` | [FiltrationAdv_spec.md](./units/FiltrationAdv_spec.md) |
| Resin Column | `ResinColumn` | `_units.py` | [ResinColumn_spec.md](./units/ResinColumn_spec.md) |
| Resin Column (Adv) | `ResinColumnAdv` | `_units_adv.py` | [ResinColumnAdv_spec.md](./units/ResinColumnAdv_spec.md) |
| Reverse Osmosis | `ReverseOsmosis` | `_units.py` | [ReverseOsmosis_spec.md](./units/ReverseOsmosis_spec.md) |
| Neutralization | `NeutralizationTank1` | `_units.py` | [NeutralizationTank1_spec.md](./units/NeutralizationTank1_spec.md) |
| Boiler/TG | `BoilerTurbogenerator` | `_units.py` | [BoilerTurbogenerator_spec.md](./units/BoilerTurbogenerator_spec.md) |
| VPSA | `VacuumPSA` | `_units.py` | [vpsa_spec.md](vpsa_spec.md) |
| HemDx Reactor | `HemDxCSTR` | `_units.py` | Formulation reactor (no spec doc) |

## Notes

- Advanced units are re-exported from `_units.py` for convenience (`u.DiafiltrationAdv`, `u.FiltrationAdv`, `u.ResinColumnAdv`).
- Filtration and diafiltration presets define MF/UF/NF operating defaults; config files override as needed.
- Centrifuge auto-distributes soluble species based on `moisture_content` for realistic liquid following.

## Module Constants

| Constant    | Value     | Description             |
| :---------- | :-------- | :---------------------- |
| `_gal2m3`   | 0.003785  | Gallons to cubic meters |
| `_gpm2m3hr` | 0.227124  | GPM to m³/hr            |
| `_hp2kW`    | 0.7457    | Horsepower to kilowatts |
| `_Gcal2kJ`  | 4,184,000 | Gcal to kJ              |

*Last Updated: 2026-02-11*

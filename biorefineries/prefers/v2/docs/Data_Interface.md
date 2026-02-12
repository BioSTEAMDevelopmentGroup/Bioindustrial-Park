# Data Interface

**Source File:** `biorefineries/prefers/v2/_process_settings.py`

This file is the single source of truth for economic and environmental parameters external to the process physics.

## Global Utility Prices

`load_process_settings()` injects utility prices and agent settings into BioSTEAM:

| Utility                   | Price | Unit  | Notes |
| :------------------------ | :---- | :---- | :---- |
| **Electricity (SG)**      | 0.03  | $/kWh | SP Group tariff (Singapore) |
| **Electricity (Global)**  | 0.065 | $/kWh | Global average fallback |
| **Process Water**         | 0.0002 | $/kg | ($0.2/m3) |
| **High Pressure Steam**   | ~0.977 | $/kg | `17.6/1e3/18.01528` in `load_process_settings()` |
| **Medium Pressure Steam** | ~0.849 | $/kg | `15.3/1e3/18.01528` in `load_process_settings()` |
| **Low Pressure Steam**    | ~0.733 | $/kg | `13.2/1e3/18.01528` in `load_process_settings()` |
| **Cooling Water**         | ~0.0015 | $/kg | `0.027/1e3/18.01528` in `load_process_settings()` |

## Material Prices (TEA)

Defined in the `price` dictionary. Selected entries used by v2 systems:

| Material | Price ($/kg) |
| :------- | :----------- |
| Glucose | 0.42 |
| NH3 | 0.46 |
| NH3_25wt | 0.46 |
| NH4OH | 0.40 |
| H2SO4 | 0.05 |
| NaOH | 0.28 |
| KH2PO4 | 1.31 |
| NaCl | 0.078 |
| Glycerol | 0.45 |
| SodiumAscorbate | 4.35 |
| Ethanol | 0.66 |
| GammaCyclodextrin | 6.0 |
| Nicotinamide | 9.0 |

## GWP Factors (LCA)

Defined in `GWP_CFs` and applied via `set_GWPCF`/`set_GWPCF_Multi`.

| Factor | Value (kg CO2e/kg) | Source |
| :----- | :----------------- | :----- |
| Electricity (SG) | 0.55 | Ecoinvent 3.11 |
| Glucose | 1.61 | Ecoinvent 3.11 (GLO) |
| Sulfuric Acid | 0.165 | Ecoinvent 3.11 (RoW) |
| Ammonia (SEA) | 2.84 | Ecoinvent 3.11 (SEA) |
| Ammonia (CN) | 5.07 | Ecoinvent 3.11 (CN) |
| NaOH | 1.41 | Ecoinvent 3.11 (RoW) |
| NaCl | 0.268 | Ecoinvent 3.11 (GLO) |
| Ethanol | 0.90 | Ecoinvent 3.11 (RoW) |
| GammaCyclodextrin | 4.75 | Assumed |
| Nicotinamide | 7.26 | FineChem2 prediction |

## Helper Functions

| Function | Purpose |
| :------- | :------ |
| `load_process_settings()` | Inject utilities, CEPCI, and electricity pricing into BioSTEAM |
| `set_GWPCF(obj, name, dilution=None)` | Set single GWP factor on streams/utilities |
| `set_GWPCF_Multi(obj, names, dilutions=None)` | Combine multiple GWP factors for composite feeds |

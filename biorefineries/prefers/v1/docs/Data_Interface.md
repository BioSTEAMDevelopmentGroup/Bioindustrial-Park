# Data Interface

**Source File:** `v1/process_settings.py`

This file is the "Single Source of Truth" for economic and environmental parameters external to the process physics.

## Global Utility Prices

These prices are injected into the BioSTEAM settings object when `load_process_settings()` is called.

| Utility                   | Price    | Unit  | Notes                       |
| :------------------------ | :------- | :---- | :-------------------------- |
| **Electricity (SG)**      | $0.03    | $/kWh | SP Group Tariff (Singapore) |
| **Electricity (Global)**  | $0.065   | $/kWh | Global Average fallback     |
| **Process Water**         | $0.0002  | $/kg  | ($0.2/m3)                   |
| **High Pressure Steam**   | ~$0.977  | $/kg  | ($17.6/kmol)                |
| **Medium Pressure Steam** | ~$0.849  | $/kg  | ($15.3/kmol)                |
| **Low Pressure Steam**    | ~$0.733  | $/kg  | ($13.2/kmol)                |
| **Cooling Water**         | ~$0.0015 | $/kg  | ($0.027/kmol)               |

## Material Prices (TEA)
Defined in `price` dictionary.

| Material      | Price ($/kg) |
| :------------ | :----------- |
| Glucose       | 0.42         |
| Ammonia       | 0.46         |
| Sulfuric Acid | 0.05         |
| NaOH          | 0.28         |
| KH2PO4        | 1.31         |
| Glycerol      | 0.45         |

## GWP 100 Factors (LCA)
Defined in `GWP_CFs` dictionary.

| Factor               | Value (kg CO2e/kg) | Source         |
| :------------------- | :----------------- | :------------- |
| **Electricity (SG)** | 0.55               | Ecoinvent 3.11 |
| **Glucose**          | 1.61               | Global         |
| **Sulfuric Acid**    | 0.165              | RoW            |
| **Ammonia (SEA)**    | 2.84               | SEA            |
| **Ammonia (CN)**     | 5.07               | CN             |

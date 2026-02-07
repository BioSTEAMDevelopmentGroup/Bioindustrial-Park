# SeedTrain Unit Specification

**Class:** `SeedTrain`
**Type:** Fermentation Seed Train (Multi-stage)
**Source:** `v1/_units.py`

## Purpose
Prepares inoculum for main fermentation by growing biomass in successively larger reactor stages (scale-up ratio ~10x per stage).

## Design Basis

| Parameter         | Default | Unit | Description                         |
| :---------------- | :------ | :--- | :---------------------------------- |
| **N_stages**      | 5       | -    | Number of reactors in series        |
| **N_trains**      | 2       | -    | Number of parallel seed trains      |
| **Cycle Time**    | 16      | hr   | Batch time per train (tau_batch)    |
| **Turnover Time** | 8       | hr   | tau_batch / N_trains (tau_turnover) |
| **Temperature**   | 32      | °C   | Operating temperature (305.15 K)    |

## Process Architecture

1. **Staging Logic**:
   - Volume scales by 10x per stage.
   - Stage N volume = `total_vol * 10^-N_stages * 10^(stage_num)`.
2. **Reactions**:
   - Competing reactions: Glucose -> Ethanol, Biomass, Glycerol, Succinic Acid.
   - Validated against *Z. mobilis* and yeast kinetics (default set).
3. **Control**:
   - Isothermal operation (cooling duty calculated).
   - Vents CO2/O2 gases.

## Cost Model

Costing is applied per stage using BioSTEAM's decorator pattern.

| Stage    | Volume basis | Cost Basis (2010$) | Reference                  |
| :------- | :----------- | :----------------- | :------------------------- |
| Stage #1 | 20 gal       | \$37,700           | NREL Humbird et al. (2011) |
| Stage #2 | 200 gal      | \$58,300           | "                          |
| Stage #3 | 2,000 gal    | \$78,800           | "                          |
| Stage #4 | 20,000 gal   | \$176,000          | "                          |
| Stage #5 | 200,000 gal  | \$590,000          | "                          |

*Note: Agitators are costed separately for stages 4 and 5.*

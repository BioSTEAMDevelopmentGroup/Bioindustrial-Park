# Leghemoglobin Fermentation Parameters

**Date:** 2026-01-22
**Purpose:** Define critical fermentation parameters (Y_LegHb, Y_b, Respiration, Carbon Consumption) for Leghemoglobin production in *Pichia pastoris* (or alternative hosts) based on Zotero library and internal code review.

## 1. Key Performance Indicators (KPIs)

| Parameter               | Value (Current Code) | Value (Lit. - Pichia) | Value (Lit. - K. marxianus) | Notes                                                                                                                          |
| :---------------------- | :------------------- | :-------------------- | :-------------------------- | :----------------------------------------------------------------------------------------------------------------------------- |
| **Titer (P)**           | **7.27 g/L**         | **3.5 g/L**           | **7.27 g/L**                | Current code cites Tian et al. (2024) for *K. marxianus*. Shao et al. (2022) reports 3.5 g/L for *P. pastoris*.                |
| **Productivity**        | **0.10 g/L/h**       | ~0.05 g/L/h           | High                        | Based on 72h batch time.                                                                                                       |
| **Biomass Yield (Y_b)** | **0.43 g/g**         | **0.40 - 0.60 g/g**   | -                           | 0.43 is a conservative estimate for *Pichia* on glycerol/methanol mixes. Growth on glycerol alone is typically ~0.55-0.65 g/g. |
| **Product Yield (Y_p)** | **~2.2 wt%**         | ~1-3.5% (est)         | -                           | Calculated from titer and total substrate feed.                                                                                |
| **Heme Binding**        | -                    | **93%**               | -                           | Critical quality attribute (Shao et al. 2022).                                                                                 |

## 2. Stoichiometry & Reaction Splitting
Based on analysis of `_config1.py` and literature:

### A. Carbon Consumption
**Does all carbon feedstock (glucose, glycerol) get consumed?**
*   **Literature**: In highly optimized fed-batch, residual carbon is kept near zero (<1 g/L) to prevent Crabtree effect (ethanol production) or repression. "All" is effectively consumed (>99.5%).
*   **Current Model**: Assumes **100% conversion** split between Biomass, Product, and Respiration.
    *   Specific Reaction: `Glucose -> Product` (Yield = $Y_{LegHb}$)
    *   Growth Reaction: `Glucose -> Biomass` (Yield = $Y_b \times (1 - Y_{LegHb\_correction})$)
    *   Respiration: `Glucose + O2 -> CO2 + H2O` (Consumes **remainder** of carbon).

### B. Respiration (CO2 Evolution)
*   **Range**: Respiration typically accounts for **40-60%** of carbon flux in high-density heterologous protein production.
*   **Recommendation**: 
    *   If $Y_b \approx 0.45$, Product $\approx 0.02$, then Respiration $\approx 0.53$ (53%).
    *   Current code calculates this dynamically: `1 - cell_growth_X - fermentation_X`. This is physically correct.

## 3. Literature References (Zotero)

### [1] High-level secretory production of leghemoglobin in Pichia pastoris...
*   **Authors**: Shao et al. (2022), *Bioresource Technology*
*   **Key Data**:
    *   Titer: **3.5 g/L** (Highest reported for Pichia secretory).
    *   Strategy: Heme pathway engineering + gene dosage.
    *   Heme Binding: 93%.
    *   Fold Increase: 83-fold over wild type.

### [2] Recent advances in microbial fermentation...
*   **Authors**: Yao et al. (2025), *World J Microbiol Biotechnol*
*   **Key Data**: Review paper. Likely source of the "7.27 g/L" figure (citing Tian et al.).

### [3] Safety Evaluation of Soy Leghemoglobin...
*   **Authors**: Fraser et al. (2018), *Int. J. Toxicol.* (Impossible Foods)
*   **Key Data**:
    *   Confirms *Pichia pastoris* host.
    *   Focus on safety/purity (No mutagenicity).
    *   implies highly efficient scale-up (commercial scale).

## 4. Recommendations for Bioreactor Setup
To increase precision in `_config1.py`:
1.  **Split Yields by Phase**: 
    *   **Growth Phase (Glycerol)**: $Y_{b} \approx 0.6$ g/g. No product.
    *   **Induction Phase (Methanol/Glucose)**: $Y_{b} \approx 0.3-0.4$ g/g. High product formation ($Y_{p}$ peak). High respiration.
2.  **Sensitivity Analysis Range**:
    *   **$Y_b$**: Range **0.35 - 0.55**.
    *   **$Y_{LegHb}$**: Range **1.0 - 5.0 g/L** (titer check) or **0.5% - 4.0%** (yield check).
    *   **Residual Carbon**: Assume 0% (efficient) to 5% (inefficient/overflow metabolism).


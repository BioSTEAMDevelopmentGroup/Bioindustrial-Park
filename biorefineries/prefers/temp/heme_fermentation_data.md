# Heme Precision Fermentation Parameters

This document summarizes the key fermentation parameters for **free heme** production (non-protein bound) based on recent literature found in the Zotero library.

## Summary Table

| Host Strain       | Titer (g/L) | Productivity (mg/L/h) | Yield (g/g glucose) | Extracellular %  | Reference         |
| :---------------- | :---------- | :-------------------- | :------------------ | :--------------- | :---------------- |
| **C. glutamicum** | 1.59        | N/A                   | N/A                 | 45.5%            | Yang et al., 2024 |
| **E. coli**       | 1.03        | 21.5                  | N/A                 | High (Secretory) | Choi et al., 2022 |
| **S. cerevisiae** | 0.38        | 4.2                   | N/A                 | N/A              | Guo et al., 2024  |
| **C. glutamicum** | 0.31        | 6.44                  | ~0.0021             | 79%              | Ko et al., 2021   |
| **B. subtilis**   | 0.25        | ~1.72                 | N/A                 | 89%              | Yang et al., 2023 |

> **Note**: Respiration yield and Carbon Conversion percentages were not explicitly reported in the abstracts or accessible full text for most papers. Standard aerobic respiration typically consumes ∼40-60% of carbon source.

## Detailed Data by Reference

### [1] Multi-modular metabolic engineering of heme synthesis in C. glutamicum (2024)
*   **Reference**: Yang et al. (Zotero: `QY8VNLGG`)
*   **Strain**: *Corynebacterium glutamicum* HS12
*   **Titer**: 1592 mg/L (Total Iron-porphyrin derivatives)
*   **Secretion**: 45.5% extracellular
*   **Pathway**: Siroheme-dependent (SHD) pathway identified as optimal.
*   **Modifications**: RBS engineering, Heme oxygenase knockout (prevent degradation).

### [2] Improved production of heme using metabolically engineered E. coli (2022)
*   **Reference**: Choi et al. (Zotero: `G9EF7CRT`)
*   **Strain**: *Escherichia coli* HAEM7
*   **Titer**: 1.03 g/L
*   **Productivity**: 21.5 mg/L/h
*   **Strategy**: Optimization of carbon sources, iron concentration, pH, induction, and feeding.
*   **Carbon Source**: "Different carbon sources examined" (likely Glucose/Glycerol).

### [3] Animal-free heme production... in C. glutamicum (2021)
*   **Reference**: Ko et al. (Zotero: `RDW3FQ3W`)
*   **Strain**: *C. glutamicum* (Engineered)
*   **Titer**: 309.18 mg/L (Total), 242.95 mg/L (Secreted)
*   **Yield**: 0.61 mmol/mol glucose (~0.0021 g/g glucose)
*   **Productivity**: 6.44 mg/L/h
*   **Feedstock**: Glucose.

### [4] Improved biosynthesis of heme in B. subtilis (2023)
*   **Reference**: Yang et al. (Zotero: `U9HM7W7V`)
*   **Strain**: *Bacillus subtilis* BSH11
*   **Titer**: 248.26 mg/L (at 144 h)
*   **Extracellular**: 221.83 mg/L (89%)
*   **Carbon Source**: Glucose (initial) then Sucrose (feed).
*   **Key Insight**: Secretion is very efficient in B. subtilis (89%).

### [5] Multidimensional engineering of S. cerevisiae (2024)
*   **Reference**: Guo et al. (Zotero: `GU9TRNLS`)
*   **Strain**: *Saccharomyces cerevisiae* R5-M
*   **Titer**: 380.5 mg/L
*   **Productivity**: 4.2 mg/L/h
*   **Tolerance**: Engineered for heme tolerance and ROS quenching.

## Carbon & Respiration Notes
*   **Respiration Yield**: Specific data (e.g., "50% carbon to CO2") was not found in the texts.
*   **Carbon Consumption**: In *B. subtilis*, 400 g/L sucrose feed was used, suggesting high carbon consumption for relatively low product titer (<0.3 g/L), implying a low mass yield and significant carbon loss to biomass/respiration or by-products.
*   **By-products**:
    *   *B. subtilis*: ALA (precursor) accumulation/outflow observed.
    *   *General*: Mixed acid fermentation products (acetate, lactate) are common in E. coli/B. subtilis if overflow metabolism occurs, though fed-batch usually aims to minimize this.


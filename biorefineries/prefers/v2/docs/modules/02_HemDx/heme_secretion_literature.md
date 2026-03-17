# Heme b Fermentation Performance: Literature Reference Table

**Purpose:** Reference data for the Secretion Fraction (SF) parameter used in the NHemDx process model (config1 baseline: SF = 0.45 for *C. glutamicum* R302 fed-batch fermentation).

**Secretion Fraction definition:** SF = Extracellular heme / Total heme × 100%

---

## Table 1. Reported Heme b Production Performance Across Engineered Microbial Strains

| # | Host | Strain / Key Genotype | Mode | Total Titer (mg/L) | Extra­cellular (mg/L) | Intra­cellular (mg/L) | **SF (%)** | Productivity (mg/L/h) | Yield | Ref |
|---|------|-----------------------|------|-------------------|----------------------|----------------------|:----------:|----------------------|-------|-----|
| 1 | *E. coli* BL21(DE3) | HAEM7 (ΔldhA ΔyfeX Δpta; C5 pathway; ccmABC↑) | Fed-batch (+L-Glu feed) | 239.2 ± 7.2 | 151.4 ± 5.6 | 87.7 ± 1.7 | **63.3** | 3.32 | 0.14 mg/g glucose^a | 1 |
| 2 | *E. coli* BL21(DE3) | HAEM7-opt (Round H24; glycerol-based) | Fed-batch (glycerol) | 1034.3 ± N/R | 471.1 | 563.2 | **45.6** | 21.5 | 4.14 mg/g glycerol | 2 |
| 3 | *C. glutamicum* | ΔSAT::ALSDtE-YHQBA(c) (DtxR↑; heme exporters↑; heme-binding ΔKO; cell-wall mod.) | Fed-batch | 309.18 ± 16.43 | 242.95 ± 11.45 | 66.23 | **78.6** | 6.44 | 2.09 mg/g glucose^b | 3 |
| 4 | *C. glutamicum* | HS12 (SHD pathway; RBS-tuned; Δheme oxygenase) | Fed-batch | 98.9^c | ~53.0 | ~45.9 | **46.4**^c | 3.53^c | N/R | 4 |
| 5 | *B. subtilis* | BSH11 (C5 + SHD pathway; heme exporters↑; 10 L) | Fed-batch (10 L) | 248.26 | 221.83 | 26.43 | **89.3** | 1.72^d | N/R | 5 |
| 6 | *S. cerevisiae* CENPK.2-1C | R5-M (ALA C4+C5; HEM2/3/4/12/13↑; ΔHMX1 ΔHAP1 ΔROX1; CARM-evolved; 5 L) | Fed-batch (5 L) | 380.5 | — | 380.5 | Intracellular only^e | 4.2 | N/R | 6 |
| 7 | *S. cerevisiae* D452-2 | ΔHMX1/H3&12+PUG1 (HEM3/12↑; PUG1↑; ΔHMX1) | Fed-batch | 28 | — | 28 | Intracellular only^e | 0.42^d | 0.091 mg/g glucose | 7 |
| 8 | *S. cerevisiae* KCCM 12638 | ΔHMX1_H2/3/12/13 (HEM2/3/12/13↑; ΔHMX1) | Fed-batch | 67 | — | 67 | Intracellular only^e | N/R | N/R | 8 |

> **Abbreviations:** SF = Secretion Fraction; Extra = extracellular heme; Intra = intracellular heme; — = not detected / not secreted; N/R = not reported.

---

## Notes

**^a** Zhao et al. 2018 (Ref. 1): Yield of 0.000138 mol heme/mol glucose reported for HAEM7 without L-Glu supplementation (glucose-only fed-batch: 115.5 mg/L total, 63.5% SF). The maximal result (239.2 mg/L, 63.3% SF) is from the L-Glu-supplemented fed-batch (72 h). Calculated 0.14 mg/g glucose using 0.000138 mol/mol × 616.5 g/mol ÷ 180 g/mol × 1000.

**^b** Ko et al. 2021 (Ref. 3): Reported yield = 0.61 mmol heme/mol glucose; converted: 0.61 × 10⁻³ mol × 616.5 g/mol ÷ 180 g/mol = 2.09 mg/g glucose.

**^c** Yang et al. 2024 (Ref. 4): The HS12 strain accumulates primarily **iron-containing porphyrin derivatives (ICPDs)**; total ICPDs = 1592 mg/L with SF = 45.5% at 28 h (productivity 56.86 mg/L/h). Pure heme b by HPLC = 98.9 mg/L with SF = 46.4%. Values reported here use heme b (HPLC) data for consistency. The SF of ~46% is consistent regardless of whether ICPDs or pure heme are used.

**^d** Productivity calculated from reported titer and fermentation time (Yang 2023: 248.26 mg/L ÷ 144 h ≈ 1.72 mg/L/h; Lee 2024: 28 mg/L ÷ 66 h ≈ 0.42 mg/L/h).

**^e** *S. cerevisiae* produces heme exclusively as an **intracellular product**: heme is almost entirely bound to hemoproteins inside the cell, and extracellular free heme was not detected in any of the three yeast studies (Refs 6–8). The PUG1 heme excretion mechanism was confirmed to be negligible under normal fermentation conditions (Ref. 7). No SF is calculable for these strains.

---

## Contextual Summary for NHemDx Model

The config1 model baseline of **SF = 0.45** (45%) for *C. glutamicum* R302 is well-supported by literature:

- **Ko 2021 (*C. glutamicum*):** SF = 78.6% — highest reported for *C. glutamicum*; reflects optimized membrane engineering and exporter overexpression.
- **Yang 2024 (*C. glutamicum* HS12):** SF = 46.4% (heme) / 45.5% (ICPDs) — closely matches the model baseline of 0.45.
- **Zhao 2018 (*E. coli*):** SF = 63.3% — *E. coli* with explicit heme exporter (CcmABC) overexpression.
- **Choi 2022 (*E. coli*):** SF = 45.6% — similar to model baseline despite higher absolute titers.
- **Yang 2023 (*B. subtilis*):** SF = 89.3% — *B. subtilis* uniquely high, likely due to natural heme export mechanisms.
- **S. cerevisiae (Refs 6–8):** Intracellular product only — yeast is not suitable for secretory heme production without dedicated engineering; no free extracellular heme detected.

The SF = 0.45 baseline in config1 represents a **conservative mid-range** for bacterial heme production, appropriate for *C. glutamicum* without extreme exporter engineering.

---

## References

1. Zhao, X. R., Choi, K. R. & Lee, S. Y. Metabolic engineering of *Escherichia coli* for secretory production of free haem. *Nat. Catal.* **1**, 720–728 (2018). https://doi.org/10.1038/s41929-018-0126-1

2. Choi, K. R., Yu, H. E., Lee, H. & Lee, S. Y. Improved production of heme using metabolically engineered *Escherichia coli*. *Biotechnol. Bioeng.* **119**, 3178–3193 (2022). https://doi.org/10.1002/bit.28194

3. Ko, Y. J. *et al.* Animal-free heme production for artificial meat in *Corynebacterium glutamicum* via systems metabolic and membrane engineering. *Metab. Eng.* **66**, 217–228 (2021). https://doi.org/10.1016/j.ymben.2021.04.013

4. Yang, Q., Sun, X., Wang, H., Chen, T. & Wang, Z. Multi-modular metabolic engineering of heme synthesis in *Corynebacterium glutamicum*. *Synth. Syst. Biotechnol.* **9**, 285–293 (2024). https://doi.org/10.1016/j.synbio.2024.02.008

5. Yang, S. *et al.* Improved biosynthesis of heme in *Bacillus subtilis* through metabolic engineering assisted fed-batch fermentation. *Microb. Cell Fact.* **22**, 102 (2023). https://doi.org/10.1186/s12934-023-02077-3

6. Guo, Q. *et al.* Multidimensional engineering of *Saccharomyces cerevisiae* for the efficient production of heme by exploring the cytotoxicity and tolerance of heme. *Metab. Eng.* **85**, 46–60 (2024). https://doi.org/10.1016/j.ymben.2024.07.007

7. Lee, H. J., Shin, D. J., Nho, S. B., Lee, K. W. & Kim, S. K. Metabolic engineering of *Saccharomyces cerevisiae* for fermentative production of heme. *Biotechnol. J.* **19**, e202400351 (2024). https://doi.org/10.1002/biot.202400351

8. Lee, H. J. *et al.* Enhanced heme production in industrial *Saccharomyces cerevisiae* through metabolic engineering. *npj Sci. Food* **9**, 244 (2025). https://doi.org/10.1038/s41538-025-00618-1

---

*Last updated: 2025 | Sources: Zotero library – "Heme Molecule literatures"*

# Task 2.4 — TEA, LCA, and Optimization for Deployment (2026 Update)

> **Task Director:** Jeremy Guest (UIUC)  
> **Task Investigators:** Yalin Li (Rutgers), Vijay Singh (UIUC), Wen Shan Yew (NUS), Won Jae Choi (SIT), Na Wei (UIUC)  
> **Updated:** February 2026 (Year 2 progress review)

---

## Background and Motivation

Precision fermentation offers a promising route to produce high-value food-grade nutrients (heme, lactoferrin, vitamins) at industrial scale. Task 2.4 develops **open-source process simulation, TEA, and LCA** to identify the most promising production pathways, quantify their economic and environmental performance, and guide R&D prioritization across the PreFerS research center.

Using the **BioSTEAM** platform (created by PI's group; Cortes-Peña et al., 2020), we build rigorous process models that capture downstream processing physics, propagate uncertainty, and benchmark precision-fermented products against existing alternatives.

---

## Year 1–2 Accomplishments

### Delivered (as of Feb 2026)

| Deliverable | Detail |
|---|---|
| **Unit operation library** | 15 custom units for precision fermentation DSP (fermentation, separation, membranes, chromatography, formulation) |
| **Mechanistic models** | 5 physics-based units: FiltrationAdv, DiafiltrationAdv, ResinColumnAdv (flux equations, column-volume breakthrough, MWCO-dependent rejection) |
| **LegHb TEA/LCA** | Leghemoglobin (*P. pastoris*) — 1 config, intracellular recovery. MSP: **$11.53/kg**, GWP: **16.72 kg CO₂-eq/kg** |
| **HemDx TEA/LCA** | N-HemoDextrin (*C. glutamicum*) — 3 configs (base/intra/extra). MSP: **$7.28–8.85/kg**, GWP: **8.08–9.54 kg CO₂-eq/kg** |
| **Uncertainty analysis** | 500-sample Latin Hypercube with 20+ uncertain parameters per route |
| **Platform extension** | Extended BioSTEAM from cellulosic biorefineries to precision fermentation for food-grade nutrients — a new application domain |

### Strategic Significance

1. **First-mover TEA/LCA for heme ingredients:** LegHb and HemDx are not currently modeled elsewhere in the research center or published literature at this level of process detail. Our work provides a meaningful techno-economic preview for these representative precision fermentation products.
2. **Reusable unit operation infrastructure:** ResinColumnAdv and DiafiltrationAdv are workhorse separation operations across **all** precision fermentation downstream processes — not just heme. These mechanistic models transfer directly to Lactoferrin, VB12, and other future targets.
3. **LegHb as Lactoferrin template:** The LegHb process route (yeast fermentation → centrifuge → cell disruption → UF/DF → pasteurization) shares ~70–80% of its DSP structure with Lactoferrin production. The LegHb flowsheet will serve as a starting template for the Lactoferrin model.

---

## Updated Yearly Plan (2026 Revision)

> **Simple format** — designed for direct copy-paste into presentation slides.  
> **Design principle:** Use generic categories (not specific product names) so the slide ages well if center priorities shift. Specific products belong in speaker notes and this reference document.

| Year | Focus |
|---|---|
| **Year 1** ✅ | Core unit operation library for precision fermentation downstream processing |
| **Year 2** ✅★ | TEA/LCA of representative precision fermentation products; uncertainty analysis |
| **Year 3** | TEA/LCA of prioritized precision fermentation products; MCDA framework development |
| **Year 4** | ML-enhanced unit operation models; MCDA of all products; deployment scenarios |
| **Year 5** | MCDA under uncertainty → deployment recommendations; AI-assisted process design |

> ★ **Current:** Year 2 complete (Feb 2026)
>
> **Reading key:** "Representative" = products chosen to demonstrate the platform (currently: heme LegHb & HemDx). "Prioritized" = products the experimental team currently focuses on (currently: Lactoferrin, VB12). The plan extends to future precision fermentation targets as the center's portfolio evolves.

---

## Year 3–5 Detailed Approaches

### 1. Product Expansion: Lactoferrin (Year 3 — Primary Target)

**Rationale:** Lactoferrin is the research center's current experimental priority. An ~80 kDa iron-binding glycoprotein, typically produced in *P. pastoris* (same host as LegHb), it shares substantial DSP overlap with the LegHb route.

**Approach:**
- Start from the LegHb flowsheet template — modify fermentation kinetics, protein size (80 vs 16 kDa), and adjust membrane MWCO settings
- If secreted expression: Remove cell disruption (similar to HemDx Config 3 pathway concept)
- If intracellular: Same DSP structure as LegHb, adjust for larger protein
- May add cation-exchange chromatography step (ResinColumnAdv already available)
- Full TEA/LCA with uncertainty, benchmarked against animal-derived lactoferrin

**Feasibility: ✅ HIGH** — ~70–80% of unit operations transfer directly from LegHb. Main effort is fermentation kinetics and chromatography configuration.

### 2. Product Expansion: VB12 (Year 3 Scoping → Year 4 Complete)

**Rationale:** Vitamin B12 (cobalamin, 1355 Da) is another priority target at the research center. Produced by *C. glutamicum* or *Pseudomonas denitrificans*, it has a distinct DSP from heme proteins.

**Approach (Year 3 scoping):**
- Literature review of VB12 production routes and published process data
- Identify DSP gaps: likely need new unit operations for chemical conversion (hydroxocobalamin → cyanocobalamin) and crystallization
- Begin flowsheet architecture design

**Approach (Year 4 complete):**
- Develop 2–3 new unit operations (chemical conversion, crystallization, possibly solvent extraction)
- Complete VB12 flowsheet, TEA/LCA with uncertainty
- Compare with chemical synthesis route

**Feasibility: ⚠️ MODERATE** — ~50–60% unit operation transfer. Needs 2–3 new units. Staggered timeline (scope Year 3, deliver Year 4) is realistic.

### 3. MCDA Framework (Year 3 Framework → Year 4–5 Execution)

**Rationale:** As the product portfolio grows (heme, lactoferrin, VB12, others), we need a systematic framework to compare alternatives across multiple criteria.

**Approach:**
- Year 3: Define MCDA criteria (MSP, GWP, water, capital intensity, technology readiness), select methodology (TOPSIS or similar), implement Python framework
- Year 4: Run MCDA across all products at baseline; design deployment scenarios (Singapore, ASEAN)
- Year 5: MCDA under uncertainty with Monte Carlo–propagated scenario analysis → deployment recommendations

**Feasibility: ✅ HIGH** — Well-established methodology. PI has extensive MCDA experience. The framework builds on BioSTEAM's existing uncertainty engine.

### 4. ML-Enhanced Unit Operations (Year 3 Pilot → Year 4 Integration)

**Rationale:** As experimental collaborators generate laboratory data (fermentation kinetics, membrane flux, resin capacity), machine learning can refine our unit operation models beyond literature-based empirical correlations.

**Approach:**
- Year 3: Pilot — train ML models on available experimental data for 1–2 priority unit operations (e.g., membrane fouling prediction, fermentation kinetic fitting)
- Year 4: Integrate validated ML models into BioSTEAM unit operations; demonstrate improved prediction accuracy vs. empirical baselines
- Data sources: Experimental data from collaborating tasks within PreFerS

**Feasibility: ⚠️ MODERATE** — Depends on experimental data availability. Conceptually sound and well-aligned with BioSTEAM's modular architecture. Frame as "data-driven refinement" rather than "ML replacement."

**Risk:** Insufficient experimental data in Year 3 → limited training sets. **Mitigation:** Start with literature data augmentation; focus on parameters with largest sensitivity impact.

### 5. AI-Assisted Process Design (Year 5 — Exploratory)

**Rationale:** As the product portfolio and unit operation library grow, manual process design becomes a bottleneck. AI agents could accelerate discovery by reasoning over literature databases and proposing candidate process configurations for BioSTEAM evaluation.

**Concept:**
- LLM agent receives: target product specifications, available organisms, literature database of DSP options
- Agent uses structured prompts to propose 3–5 candidate process configurations
- BioSTEAM automatically evaluates each → TEA/LCA comparison
- Agent iterates to optimize

**Approach:**
- Year 4: Develop structured prompt templates and curate literature database of process configurations
- Year 5: Proof-of-concept demonstration — AI-proposed vs. manually-designed process configurations

**Feasibility: ⚠️ LOW-MODERATE** — Genuinely novel research direction. No mature frameworks exist. LLMs risk thermodynamic inconsistency. Frame as "exploratory proof-of-concept." **Potential impact is HIGH if successful** — could become a differentiating capability for BioSTEAM.

**Risk:** LLM hallucination in process design. **Mitigation:** All AI-proposed configurations must pass BioSTEAM simulation (thermodynamic consistency check) before evaluation. Human review required.

---

## Feasibility Summary

| Direction | Readiness | Timeline | Confidence |
|---|---|---|---|
| Lactoferrin TEA/LCA | Template from LegHb; shared yeast platform | Year 3 | 🟢 HIGH |
| VB12 TEA/LCA | Needs 2–3 new unit ops; different DSP | Year 3 scope → Year 4 | 🟡 MODERATE |
| MCDA framework | Established methodology; PI expertise | Year 3 → Year 5 | 🟢 HIGH |
| ML-enhanced unit ops | Depends on experimental data availability | Year 3 pilot → Year 4 | 🟡 MODERATE |
| AI-assisted process design | Novel; no mature frameworks | Year 4–5 exploratory | 🟠 LOW-MODERATE |

---

## Synergy with Other Tasks

- **Inputs from experimental tasks:** Fermentation TRYs, separation efficiency data, product specs for Lactoferrin and VB12 → feed into our process models and ML training
- **Outputs to experimental tasks:** TEA/LCA results set quantitative performance targets (e.g., "achieve titer >X g/L to reach MSP <$Y/kg") → guide R&D prioritization
- **Cross-cutting:** Reusable unit operations (ResinColumnAdv, DiafiltrationAdv, FiltrationAdv) benefit all product teams, not just heme

---

## References

- Cortes-Peña, Y. et al. (2020). BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. *ACS Sustainable Chemistry & Engineering*, 8(8), 3302–3310.
- ISO 14040:2006. Environmental management — Life cycle assessment — Principles and framework.
- Saltelli, A. et al. (2009). *Global Sensitivity Analysis: The Primer*. John Wiley & Sons.
- GFI/SPINS (2024). US Plant-Based Market Overview.
- Ahmad, M. et al. (2023). Plant-based meat analogues: A review. *Trends in Food Science & Technology*, 141, 104199.

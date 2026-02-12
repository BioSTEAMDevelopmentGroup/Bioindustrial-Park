# Task 2.4 — TEA, LCA, and Optimization for Deployment

> **Task Director:** Jeremy Guest (UIUC)  
> **Task Investigators:** Yalin Li (Rutgers), Vijay Singh (UIUC), Wen Shan Yew (NUS), Won Jae Choi (SIT), Na Wei (UIUC)

---

## Background and Motivation

Considering the expansive landscape of precision fermentation development pathways, it is critical to **identify the most promising pathways at early stages** to prioritize research and development.  Given the complex interactions of these technologies with existing infrastructure, these pathways should be optimized at the **systems level** to consider potential synergies and tradeoffs that exceed the boundaries of a single system.  Further, to expedite technology deployment, systems should be designed, evaluated, and refined considering a range of **future scenarios** to enhance resilience in the context of a rapidly changing environment.

---

## Approaches

### 1. Systematic Product Selection

- Perform a **systematic literature review** to narrow candidate products.
- Leverage the Vitamin and Mineral Nutrition Information System (VMNIS; WHO) to identify **dietary micronutrient deficiencies** in Singapore and neighboring regions.
- Perform **preliminary TEA/LCA** to shortlist candidates of large demand, high market competitiveness, and low environmental impacts, guiding experimental R&D.
- Include contextual considerations (e.g., consumers' preferences) for **Singapore and the ASEAN region**.

### 2. Open-Source Process Simulation in BioSTEAM

- Develop **open-source process simulation modules** in the [BioSTEAM](https://biosteam.readthedocs.io) platform (Cortes-Peña et al., 2020).
- Model **all biological and downstream processes**, as well as interconnected systems.
- BioSTEAM: open-source, Python-based platform for:
  - Mass/energy balances
  - TEA (including DCFROR analysis)
  - LCA (ISO 14040/14044 guidelines)
  - **Uncertainty analysis** (Monte Carlo, Latin Hypercube)

### 3. Benchmarking Metrics

Technologies will be benchmarked against existing (non-precision-fermentation) production routes using:

| Metric | Unit |
|---|---|
| Minimum Product Selling Price (MSP) | $/kg |
| Total Capital Investment (TCI) | $M |
| Life-Cycle GHG Intensity (GWP) | kg CO₂-eq/kg |
| Water Usage/Consumption | kg water/kg product |

### 4. Sensitivity & Optimization

- **Global sensitivity methods:** Morris method, Sobol's variance-based sensitivity indices (Saltelli et al., 2009).
- Prioritize parameters the system is most sensitive to for:
  - R&D (technological parameters)
  - Scenario analysis (contextual parameters)
  - Optimization (design and operation decisions)
- TEA/LCA modules quantify implications of technological advancements (e.g., increases in TRYs, new separation approaches).

### 5. AI-Assisted Process Design

- Apply a **superstructure-based process synthesis** approach.
- Leverage **heuristically guided search algorithms** and ML principles to automate generation, pruning, and optimization of superstructures.
- Use ML to **predict missing physicochemical properties** (e.g., boiling point, solubility) of reaction intermediates and products.

### 6. Multi-Criteria Decision Analysis (MCDA)

- Develop a **Python-based MCDA framework** leveraging key indicators (MSP, GWP, etc.).
- Robust methodologies (e.g., TOPSIS — Technique for Order of Preference by Similarity).
- Evaluate **multiple aspects:** technical and legislative feasibility, financial viability, environmental impacts.
- Design **distinct and representative scenarios** reflecting:
  - Design and operating decision combinations (stakeholder preferences)
  - Technological parameters (R&D levels)
  - Contextual parameters (deployment context)

---

## Yearly Research Plan

| Year | Focus |
|---|---|
| **Year 1** | Development of unit operations for prioritized precision fermentation processes |
| **Year 2** ★ | Development of unit operations for up-/downstream processes of priority precision fermentation processes |
| **Year 3** | Development of additional processes of interest; identification of MCDA criteria, indicators, and ranking methodologies |
| **Year 4** | Initial MCDA of all technology alternatives at baseline; design of deployment scenarios |
| **Year 5** | MCDA across scenarios and under uncertainty |

> ★ **Current progress:** Year 2 (as of Feb 2026)

---

## Expected Outcomes and Deliverables

1. **Rigorous simulation modules** for technology alternatives (precision fermentation, up/downstream processes, and interconnected systems)
2. **Robust MCDA framework** and open-source MCDA tool in Python
3. **Optimized systems design** under varying deployment scenarios

---

## Potential Pitfalls and Mitigation

| Risk | Mitigation |
|---|---|
| High parametric uncertainty due to early-stage nature of precision fermentation technologies | Leverage **screening techniques** (e.g., Morris method) to identify high-sensitivity parameters and restrict uncertainty analysis to those parameters |
| Computation-intensive simulations with many uncertain inputs | Focus Monte Carlo/Latin Hypercube on the most impactful parameters |

---

## Synergy with Other Tasks/Thrusts

- **Inputs:** Advancements from other tasks provide the technology basis for system simulation (fermentation TRYs, separation data, product specs).
- **Outputs:** Results from Task 2.4 can be leveraged to **prioritize technology R&D** and **set specific performance targets** for experimental teams.

---

## Open-Source Sustainability

- Long-term operation and stability ensured by complementary resources (e.g., through CABBI supported by US DOE).
- BioSTEAM platform has gained considerable momentum with users/contributors from academia, research institutions, and industries.
- Expected to evolve into an **open-source ecosystem** with far-reaching impacts.

---

## References

- Cortes-Peña, Y. et al. (2020). BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. *ACS Sustainable Chemistry & Engineering*, 8(8), 3302–3310.
- ISO 14040:2006. Environmental management — Life cycle assessment — Principles and framework.
- Saltelli, A. et al. (2009). *Global Sensitivity Analysis: The Primer*. John Wiley & Sons.

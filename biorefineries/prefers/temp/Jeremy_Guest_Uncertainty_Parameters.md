# Uncertainty and Sensitivity Analysis Parameters in Biorefinery Design
Based on the methodologies of the Guest Group (BioSTEAM frameworks), particularly referencing *sustainable lactic acid production*, *BioSTEAM-LCA*, and *quantitative sustainable design (QSD)* tutorial papers.

## 1. Financial & Economic Parameters
These parameters capture market volatility and capital estimation uncertainty.
*   **Plant Uptime / Operating Days**: Variation in annual operating hours (e.g., 90% - 96% availability).
*   **Total Capital Investment (TCI) Ratio**: A factor to account for uncertainty in capital cost estimates (e.g., ±25-30% of baseline).
*   **Discount Rate / IRR**: Often fixed for baseline MPSP (e.g., 10%), but can be a parameter in broader financial risk analysis.
*   **Tax Rate**: Sometimes varied if relevant to specific policy contexts, though often fixed.

## 2. Feedstock & Utilities (Market & Composition)
Parameters related to the cost and properties of inputs.
*   **Feedstock Price**: Price per dry ton (highly volatile).
*   **Utility Prices**:
    *   Electricity price ($/kWh)
    *   Natural gas price ($/MMBtu or $/kg)
*   **Chemical/Material Prices**:
    *   Sulfuric Acid ($/kg)
    *   Lime / Calcium Hydroxide ($/kg)
    *   Ammonia / CSL / Nutrients pricing
    *   Gypsum (or other byproduct) selling price (often modeled with potential to be negative/disposal cost).

## 3. Technical Process Parameters
These capture the performance uncertainty of the technology.

### A. Pretreatment
*   **Solids Loading**: The ratio of solids to liquid (e.g., 25% - 40%).
*   **Chemical Loading**: Acid/Base loading per unit biomass (e.g., mg H2SO4/g dry biomass).
*   **Conversion Efficiencies**:
    *   Glucan-to-Pretreated-Glucan/Glucose conversion.
    *   Xylan-to-Xylose conversion.

### B. Enzymatic Hydrolysis & Saccharification
*   **Enzyme Loading**: mg protein/g cellulose (critical cost driver).
*   **Hydrolysis Time** (Residence Time): Impacting capital cost (tank size).
*   **Solids Loading**: In the hydrolysis reactor.
*   **Saccharification Yield**: Glucan-to-Glucose conversion percentage.

### C. Fermentation
*   **Yield**: Product yield (g product / g sugar) - often the most sensitive parameter.
*   **Titer**: Final concentration (g/L) - impacts downstream separation energy.
*   **Productivity**: Rate (g/L/hr) - impacts reactor volume/CapEx.
*   **Microbial/Nutrient Requirements**: e.g., CSL loading, seed train inoculum ratio.
*   **Byproduct Formation**: Yield of impurities (e.g., acetic acid, glycerol) which burden downstream separations.

### D. Downstream Processing & Separations
*   **Separation Efficiency / Splits**: Recovery rates of the product in key separation units (e.g., chromatography, distillation).
*   **Conversion Factors**: If reactive separations (e.g., esterification/hydrolysis) are involved.
*   **Residence Times**: For specific crystallization or reaction steps.

### E. Facilities
*   **Boiler Efficiency**: Efficiency of steam generation/combustion.
*   **Wastewater Treatment**: Chemical Oxygen Demand (COD) removal efficiency or biogas yield (sometimes varied).

## 4. Life Cycle Assessment (LCA) Parameters
Parameters affecting environmental impact indicators (GWP, FEC).
*   **GWP Characterization Factors**: CO2 equivalent emissions associated with material production (e.g., kg CO2-eq/kg H2SO4).
*   **Electricity Grid Intensity**: GHG intensity of the grid (kg CO2-eq/kWh).
*   **Direct Emissions**: Process-specific emissions (e.g., from fermentation or combustion).

## Summary Note
In a typical **Monte Carlo** simulation (N=1000+ runs), these parameters are sampled simultaneously from defined distributions (Uniform, Triangle, etc.).
In a **Sensitivity Analysis** (e.g., Tornado plot or Spearman Rank Correlation), the impact of each individual parameter on key metrics (MPSP, GWP) is evaluated to identify "bottlenecks" and research priorities.

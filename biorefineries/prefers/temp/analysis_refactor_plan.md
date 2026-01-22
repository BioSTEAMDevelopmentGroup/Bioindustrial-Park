# Analysis Refactor Implementation Plan

This plan outlines the creation of analysis tools and figure generation scripts for the PreFerS v1 codebase. The goal is to separate data generation from figure plotting and establish a reusable utility library.

## 1. Directory Structure & File Organization

### Analyses Directory (`biorefineries/prefers/v1/LegHb/analyses/`)
We will create three clean, single-purpose scripts:
1.  **`gen_data_base.py`**:
    *   **Purpose**: Runs the baseline simulation for the system.
    *   **Outputs**: Generates detailed stream data (mass/energy flows) and unit design parameters.
    *   **Usage**: Primary data source for Task A (Sankey) and Task D (Block Flow).
2.  **`gen_data_mc.py`**:
    *   **Purpose**: Runs Monte Carlo simulations.
    *   **Features**: Configurable parameters (N_targets, batch_size, production scale). Supports multi-config execution (starting with Config 1).
    *   **Outputs**: MC results dataframe (MSP, GWP, etc.) and breakdown of costs.
    *   **Usage**: Primary data source for Task B (Economics) and Task C (Uncertainty).
3.  **`gen_figure.py`**:
    *   **Purpose**: Loads generated data and calls utility plotting functions to create figures.
    *   **Outputs**: Saves `.png/.svg` files to `results_YYYY.MM.DD-HH.MM`.

### Utils Directory (`biorefineries/prefers/v1/utils/`)
We will expand the utility library to support these specific figure types:
*   `sankey_utils.py`: Logic for hierarchical node aggregation and Sankey generation.
*   `plot_utils.py`: (Extensions to existing) Specialized plot styles for split-axis boxplots and KDEs.
*   `diagram_utils.py`: Tools for generating block flow diagrams from system objects or docs.

---

## 2. Detailed Task Breakdown

### Task A: Hierarchical Sankey Diagrams
*   **Target**: Two hierarchical Sankey diagrams (Carbon Mass Flow, Energy Flow).
*   **Data Generation (`gen_data_base.py`)**:
    *   Implement logic to traverse the system.
    *   **Hierarchical Aggregation**: Group units into "Steps" (e.g., Feedstock Handling, Fermentation, Separation) similar to the Mock System step division.
    *   **Carbon Tracking**: Calculate C-flow in/out for each step.
    *   **Energy Tracking**: Calculate LHV/Process Heat flow for each step.
*   **Plotting (`gen_figure.py` + `utils`)**:
    *   Use `sankey_utils.py` (refining draft from `temp/sankey.py`).
    *   Figure A1: Step-wise Carbon Flow.
    *   Figure A2: Step-wise Energy Flow.

### Task B: Economics Summary (MSP & Cost Breakdown)
*   **Target**: Combined plot with 3 Columns (Config 1, 2, 3 - currently implementing Config 1).
    *   **Upper Panel**: MSP Box Plot.
    *   **Lower Panel**: Installed Equipment Cost (ISBL/OSBL) Stacked Bar.
*   **Data Generation (`gen_data_mc.py`)**:
    *   Run MC for Config 1.
    *   Save MSP distribution.
    *   Calculate and save Cost Breakdown per "Step" (matching Task A groups).
    *   *Note*: Script should be written to easily loop over multiple Configs in the future.
*   **Plotting (`gen_figure.py` + `utils`)**:
    *   Create a reusable "Economic Panel" plotting function.
    *   **Dynamic Columns**: Auto-detect number of configs.
    *   **Breakdown**: Use Step groups for the bar chart segments.

### Task C: Uncertainty Analysis (KDE + Box Plot)
*   **Target**: Joint Kernel Density Estimate (KDE) with marginal Box Plots.
*   **Style**: Use defined `prefers` color palette (Reference `v1/utils/style.py`).
*   **Data Generation (`gen_data_mc.py`)**:
    *   Use MC results (MSP vs GWP, or similar key metrics).
*   **Plotting (`gen_figure.py` + `utils`)**:
    *   Implement "Seaborn jointplot" style wrapper using Matplotlib/Scipy (to avoid heavy dependencies if preferred, or use Seaborn if available).
    *   Ensure consistency with Project Fonts and Colors.

### Task D: Block Flow Diagram (BFD)
*   **Target**: Simplified Block Flow Diagram separating "Conversion" and "Separation".
*   **Style**: "Picture1" style (legacy style) - distinct colors for Conversion vs Separation.
*   **Logic**:
    *   **Source**: Instead of purely hardcoding, attempt to read structure from `v1/docs` or traverse the `system` object to identify "Key Units".
    *   **Separation**: Group units into "Conversion" (Upstream) and "Separation" (Downstream).
*   **Plotting**:
    *   Generate a static schematic or use a graphing tool (like Graphviz logic or simple patch drawing in Matplotlib) to represent the flow.
    *   *Note*: If dynamic generation is too complex, this may fall back to a curated drawing script, but the goal is to track key units automatically.

---

## 3. Implementation Steps

1.  **Refactor Utils**:
    *   Move/Refactor `temp/sankey.py` -> `v1/utils/sankey_utils.py`.
    *   Verify `v1/utils/style.py` has correct colors.
2.  **Create Data Scripts**:
    *   Implement `gen_data_base.py` for System traversal.
    *   Implement `gen_data_mc.py` for Monte Carlo execution.
3.  **Create Figure Script**:
    *   Implement `gen_figure.py` to orchestrate the plots.
4.  **Verification**:
    *   Run Base -> Check Sankey & BFD.
    *   Run MC (small batch) -> Check Box/KDE plots.

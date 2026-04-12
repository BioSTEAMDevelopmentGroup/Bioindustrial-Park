# PREFERS v2 One-Page Architecture

Last updated: 2026-04-01
Scope: Shared v2 core + LegHb + HemDx (config1/2/3)
Detailed Mermaid diagrams: ARCHITECTURE_MERMAID.md

## 1) Control Plane (Who runs what)

The top-level orchestrator is `_RUN.py`.

Execution order:
1. LegHb config1: `gen_data_base.py` + `gen_data_mc.py` (parallel subprocesses)
2. Pause
3. LegHb config2: `gen_data_base.py` + `gen_data_mc.py` (parallel subprocesses)
4. Pause
5. HemDx config1: `gen_data_base.py` + `gen_data_mc.py` (parallel subprocesses)
6. Pause
7. HemDx config2: `gen_data_base.py` + `gen_data_mc.py` (parallel subprocesses)
8. Pause
9. HemDx config3: `gen_data_base.py` + `gen_data_mc.py` (parallel subprocesses)
10. Pause
11. Figure stage: 5 parallel `gen_figure.py` runs (LegHb c1/c2, HemDx c1/c2/c3)

## 2) Dependency Spine (How scripts connect)

```text
_RUN.py
  -> LegHb/analyses/gen_data_base.py
      -> LegHb/_models.py:create_model(config)
          -> LegHb/system/_config1.py or _config2.py:create_LegHb_system()
          -> LegHb/_tea_config1.py or _tea_config2.py:PreFerSTEA
      -> utils/report_generator.py + utils/sankey_utils.py

  -> LegHb/analyses/gen_data_mc.py
      -> LegHb/_models.py:create_model(config)
      -> multiprocessing worker evaluation (batched MC)
      -> writes MC datasets for figures

  -> LegHb/analyses/gen_figure.py
      -> reads baseline + MC outputs
      -> renders summary/uncertainty/sensitivity figures

  -> HemDx/analyses/gen_data_base.py
      -> HemDx/_models.py:create_model(config)
          -> HemDx/system/_config1.py|_config2.py|_config3.py:create_NHemDx_system()
          -> HemDx/_tea_config1.py|_tea_config2.py|_tea_config3.py:PreFerSTEA
      -> utils/report_generator.py + utils/sankey_utils.py

  -> HemDx/analyses/gen_data_mc.py
      -> HemDx/_models.py:create_model(config)
      -> multiprocessing worker evaluation (batched MC)
      -> writes MC datasets for figures

  -> HemDx/analyses/gen_figure.py
      -> reads baseline + MC outputs
      -> renders summary/uncertainty/sensitivity figures
```

## 3) Shared Core Layer

- `_process_settings.py`
  - Central prices, utility settings, and GWP characterization factors.
  - Called early by model factories via `load_process_settings()`.
- `_units.py`, `_units_adv.py`
  - Custom BioSTEAM units used by both product lines.
- `_tea.py`
  - Base `PreFerSTEA` class used by config-specific TEA wrappers.
- `utils/`
  - `utils.py`: output directory and IO helpers.
  - `report_generator.py`: TEA/LCA breakdown export.
  - `sankey_utils.py`, `plots.py`, `plot_utils.py`, `style.py`: visualization stack.

## 4) Runtime Behavior by Script Type

- Baseline (`gen_data_base.py`)
  - Builds model/system at chosen config.
  - Runs baseline simulation + sensitivity summary.
  - Exports baseline metrics, stream/unit summaries, breakdown workbook, sankey JSON.

- Monte Carlo (`gen_data_mc.py`)
  - Samples model parameters and evaluates in worker processes.
  - Stores scenario datasets (fixed/varied scale and filtered groups).
  - Produces statistical inputs for figure scripts.

- Figures (`gen_figure.py`)
  - Loads latest config result directory.
  - Builds distributions, correlation/sensitivity, and contribution figures.

## 5) Output Contract (per config run)

```text
<module>/analyses/results_<config>_<timestamp>/
  README.txt
  data/
    analysis_manifest.json
    analysis_metadata.json
    baseline_metrics.*
    baseline_streams.*
    baseline_unit_design.*
    breakdown and MC result tables/files
    sankey_*.json
  figure/
    distribution/correlation/sensitivity/breakdown plots
```

## 6) Key Cross-Module Architecture Risks

1. `_RUN.py` uses a hardcoded Python executable path, reducing portability across machines/environments.
2. Parallel subprocess output is consumed process-by-process; heavy logs can block a non-drained process pipe.
3. System routers expose multi-config availability but package-level exports are primarily config1-bound.
4. HemDx baseline script imports `PreFerSTEA` from config1 directly; config2/3 TEA routing should be explicit for consistency.

## 7) Mental Model

Treat PREFERS v2 as a 3-layer system:
- Layer A (Control): `_RUN.py` coordinates all runs.
- Layer B (Model/System): `<Module>/_models.py` selects config-specific `system/_configX.py` + `_tea_configX.py`.
- Layer C (Analysis IO): `analyses/*.py` and `utils/*` transform simulation outputs into datasets, reports, and figures.

This separation is what enables one runner to drive both LegHb and HemDx with the same analysis pipeline pattern.

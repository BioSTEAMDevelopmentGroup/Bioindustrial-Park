<!--
Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>

This module is under the UIUC open-source license. See
github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
for license details.
-->

# Sabre pricing notes

This file merges two pieces of documentation-only pricing analysis for
sabre (nothing in the codebase loads either):

1. A cross-check of sabre's flat per-kg wastewater/solids disposal
   assumptions against `biorefineries/wwt`'s COD-based pricing method.
2. An audit of economic assumptions used by `legacy/analyses/*.py` that
   are not sourced from `data/*.yaml` (compiled 2026-07-17, originally
   `legacy_analyses_assumptions_audit.yaml`).

## COD-based wastewater pricing cross-check

Applied to sabre's actual effluent streams (`pelagic` feedstock, default
simulation).

### Method (from `biorefineries/wwt/_utils.py`)

1. `compute_stream_COD(stream)` derives each chemical's theoretical oxygen
   demand from its elemental composition (C, H, O, N, S, P), sums over the
   stream, and returns a concentration in kg-O2/m3.
2. COD load = `compute_stream_COD(stream) * stream.F_vol` (kg-COD/hr).
3. Unit disposal rate: **-$0.3676/kg-COD**, averaged from Schueller (2020)
   *Municipal Residential Wastewater Rates* excess-COD surcharge range
   ($0.127-$0.2065/lb).
4. `stream.price = rate * COD_load / stream.F_mass` ($/kg stream), the form
   BioSTEAM's TEA expects.

wwt module could not be imported as a package as-is in this checkout (its
`__init__.py` chain uses `from ... import Stream, Unit`, which assumes a
different repo layout than this one and raises `ImportError: attempted
relative import beyond top-level package`). `_chemicals.py` and `_utils.py`
were loaded directly via `importlib` with a stub `biorefineries.wwt`
package registered in `sys.modules` to satisfy their internal relative
imports, bypassing `__init__.py`. No wwt files were modified.

### Results

Sabre's actual streams were simulated (`create_ad_fermentation_system`,
`create_ad_biomethane_system`, feedstock='pelagic') and
`compute_stream_COD` applied directly to the resulting `Stream` objects.

| Stream | Flow (kg/hr) | COD conc. (mg/L) | COD-based price | sabre's current flat price |
|---|---:|---:|---:|---:|
| `wastewater` (VFA fermentation pathway) | 471,561 | 24,289 | -$0.00896/kg (-$8.93/m3) | -$0.005/kg |
| `liquid_digestate` (biomethane pathway) | 450,460 | 18,307 | -$0.00662/kg (-$6.73/m3) | -$0.002/kg |
| `acidogenic_residual_solids` (VFA AD solids cake) | 38,569 | 88,626 | -$0.02846/kg | -$0.04/kg |
| `soil_amendment` (biomethane solids cake) | 53,857 | 112,154 | -$0.03691/kg | -$0.02/kg |

sabre's current flat rates (`FERM_WASTEWATER_DISPOSAL_USD_PER_KG`,
`LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG`, `SOLIDS_DISPOSAL_USD_PER_KG`,
`SOIL_AMENDMENT_DISPOSAL_USD_PER_KG`) are set in
`analyses/msp_comparison.py` (and duplicated across
`legacy/analyses/integrated_tea.py`, `legacy/analyses/vfa_fermentation_tea.py`,
`legacy/analyses/ad_tea_final.py`).

For the `wastewater` stream, COD is dominated by residual VFAs and
unconverted seaweed polysaccharides (via `get_COD_breakdown`): acetic acid
34%, alginate 21%, fucoidan 16%, propionic acid 14%, valeric+butyric acids
16%, mannitol 5%, microbial oil <1%.

### Takeaways

- For the two **liquid** effluent streams, the COD-based method prices
  disposal roughly **1.8x (`wastewater`) and 3.3x (`liquid_digestate`)**
  higher than sabre's current flat assumptions -- at 8,000 h/yr, about
  -$33.8M/yr and -$23.9M/yr vs. sabre's current ~-$18.9M/yr and ~-$7.2M/yr.
- For the **solids cakes**, the COD-based numbers land in the same
  ballpark as sabre's existing biosolids assumptions (within ~30-85%).

### Caveat

wwt's $0.3676/kg-COD rate is a municipal-wastewater "excess COD" surcharge
calibrated for conventional, much more dilute sewage. Sabre's streams run
18,000-112,000 mg/L COD -- roughly 10-100x typical municipal strength --
so a flat linear $/kg-COD extrapolation likely overstates true disposal
cost at this concentration; streams this strong are usually candidates for
on-site anaerobic treatment (which sabre's AD pathways already partly
provide) or resource recovery rather than municipal-surcharge billing.
Treat this as an order-of-magnitude cross-check against sabre's
literature-based flat rates, not a rigorous replacement for them.

## Legacy economics audit

Assumptions used by `biorefineries/sabre/legacy/analyses/*.py` that are
**not** sourced from `data/*.yaml`. None of `data/feedstock.yaml`,
`data/preprocessing.yaml`, or `data/vfa_fermentation.yaml` cover
economics, so these were never "migrated" from anywhere -- this was the
first time they were written down in one place. Candidate seed list for a
future `data/economics.yaml`. Compiled 2026-07-17; kept verbatim (as YAML)
from the original standalone audit file rather than reformatted, since
reformatting dense structured data into prose risks introducing errors.

The `disposal_costs_usd_per_kg` block below is the direct source for the
"sabre's current flat price" column in the COD cross-check table above.

```yaml
legacy_only_economics:
  # None of these have a data/*.yaml counterpart. Grouped by category;
  # each entry lists every legacy_analyses script that defines it and
  # flags where the same concept appears under different names/values.

  nutrient_and_reagent_prices:
    mgso4_usd_per_kg: {value: 0.40, note: "industrial grade", sources: [vfa_fermentation_tea.py]}
    ammonia_usd_per_kg: {value: 0.50, note: "liquid ammonia, N source", sources: [vfa_fermentation_tea.py]}
    phosphate_usd_per_kg: {value: 0.80, note: "KH2PO4, P source", sources: [vfa_fermentation_tea.py]}
    naoh_usd_per_kg: {value: 0.35, note: "pH base", sources: [vfa_fermentation_tea.py]}
    oil_extraction_reagent_usd_per_kg_oil:
      value: 0.50
      sources: [vfa_fermentation_tea.py, integrated_tea.py]
      consistency: "same value in both files"
      sensitivity_range: [0.50, 1.00, 1.50]  # vfa_fermentation_tea.py REAGENT_SENSITIVITY

  disposal_costs_usd_per_kg:
    vfa_microfilter_retentate: {value: -0.005, sources: [vfa_fermentation_tea.py, integrated_tea.py]}
    fermentation_wastewater: {value: -0.005, sources: [vfa_fermentation_tea.py, integrated_tea.py]}
    acidogenic_residual_solids:
      value: -0.04
      sources: [vfa_fermentation_tea.py, integrated_tea.py]
      basis: "standard biosolids, $30-80/dry ton"
      sensitivity_range: [-0.02, -0.04, -0.10]  # vfa_fermentation_tea.py SOLIDS_SENSITIVITY
    liquid_digestate:
      value: -0.002
      basis: "$2/m3 wastewater disposal"
      sources: [ad_feed_tea.py, ad_tea_final.py, ad_biostimulant_price.py,
                plot_two_pretreatment_figures_from_tea.py, integrated_tea.py]
      also_present_as: "ad_heatmap.py LIQUID_DIGESTATE_USD_PER_KG (same value, different name)"
      consistency: "identical value across all 6 sources"
    solid_digestate:
      value: -0.02
      basis: "standard biosolids; heavy metals in pelagic Sargassum digestate prevent land application"
      sources: [ad_feed_tea.py, ad_tea_final.py, ad_biostimulant_price.py,
                plot_two_pretreatment_figures_from_tea.py, integrated_tea.py]
      also_present_as: "ad_heatmap.py SOLIDS_DIGESTATE_USD_PER_KG (same value, different name); integrated_tea.py SOLID_DIGESTATE_DISPOSAL_USD_PER_KG (singular 'SOLID', same value)"
      consistency: "identical value across all sources, naming is inconsistent (SOLIDS_ vs SOLID_)"

  biostimulant_pricing_usd_per_kg:
    base_case_no_revenue:
      value: 0.00
      sources: [vfa_fermentation_tea.py, ad_feed_tea.py, ad_tea_final.py,
                plot_two_pretreatment_figures_from_tea.py]
      note: "base-case assumption: no biostimulant co-product revenue"
    sensitivity_range: [0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00]
      # union of BIOSTIMULANT_SENSITIVITY_USD_PER_KG ([0,0.5,1,2] in
      # vfa_fermentation_tea.py/ad_tea_final.py) and BIOSTIMULANT_PRICES /
      # BIOSTIMULANT_PRICE_CASES ([0,0.5,1,1.5,2,2.5,3] in ad_biostimulant_price.py,
      # vfa_price_scenarios.py) -- ranges are NOT identical across scripts,
      # flagged for awareness, not necessarily a bug (different sweep resolutions).
    integrated_pathway_market_price:
      value: 0.50
      source: integrated_tea.py (BIOSTIMULANT_USD_PER_KG)
      note: >
        DIFFERENT parameter from base_case_no_revenue above -- this is a
        market-price assumption used for an NPV-at-market-price sensitivity,
        not the base-case revenue assumption. Same "$0.50/kg standard
        assumption" language as the sensitivity range's second point, but
        serves a different modeling purpose in this one script.

  feedstock_price_usd_per_kg_wet:
    standard_collection_cost:
      value: 0.02
      sources: [ad_biostimulant_price.py, plot_two_pretreatment_figures_from_tea.py,
                vfa_price_scenarios.py (FEED_PRICE_BASE), vfa_yield_sensitivities.py,
                integrated_tea.py (FEED_PRICE_BASE)]
    near_zero_baseline:
      value: 0.00
      source: hauke_stream_tables.py (FEED_PRICE)
      note: >
        hauke_stream_tables.py uses $0.00/kg where every other single-scenario
        script above uses $0.02/kg "standard collection cost" as its baseline.
        Flagged as an inconsistency worth checking -- may be intentional
        (this script's purpose is stream-table generation, not TEA), not
        verified against sabre_ad or a citation.
    tiered_scenario_cases:
      # FEED_PRICE_CASES in vfa_fermentation_tea.py and ad_tea_final.py
      # (identical in both):
      tipping_fee: -0.02
      near_zero: 0.00
      low_cost_collect: 0.02
      beach_midpoint: 0.05
      source_citation: "Rodriguez-Martinez et al. 2023 (Sci Total Environ); cleanup costs $19-85/m3 wet = $0.024-0.106/kg wet"
      note: >
        This matches data/feedstock.yaml's feed_price_sensitivity section
        almost exactly, EXCEPT feedstock.yaml renamed beach_midpoint to
        mid_cost_collect_usd_per_kg_wet and added a 5th point,
        high_cost_collect_usd_per_kg_wet: 0.10, that no legacy script has.
        feedstock.yaml's feed_price_sensitivity section is itself unread by
        any code (confirmed separately) -- both copies are documentation only.
    sweep_grids:
      vfa_heatmap_and_ad_heatmap: [-0.02, 0.00, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10]
      ad_feed_tea_sensitivity: [-0.02, 0.00, 0.02, 0.05]
      note: "different grids for different sweep purposes, not a discrepancy"

  market_reference_prices:
    biomethane_usd_per_mmbtu:
      ad_tea_final_sensitivity_points: [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
        # basis: US Henry Hub $2-8/MMBtu typical range
      integrated_tea_single_benchmark: 10.00
        # basis: European TTF benchmark
      integrated_tea_sensitivity_points: [3.0, 10.0, 14.0]
        # Henry Hub, TTF, JKM
      naming_collision: >
        ad_tea_final.py and integrated_tea.py both define a constant named
        BIOMETHANE_MARKET_MMBTU, but ad_tea_final.py's is a list of
        sensitivity points and integrated_tea.py's is a single scalar
        benchmark price -- same name, different type and meaning, in two
        different (never co-imported) modules. Not a live bug, but
        confusing if read side by side.
    crude_oil_market_usd_per_kg:
      value: 1.00
      source: integrated_tea.py (OIL_MARKET_USD_PER_KG)
      note: "midpoint of soybean oil $0.62-1.50/kg range"
    soybean_oil_reference_range_usd_per_kg:
      low: 0.62
      high: 1.50
      sources: [vfa_heatmap.py, vfa_price_scenarios.py]

  vfa_ids_constant:
    value: ["AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"]
    sources: [vfa_fermentation_tea.py (module-level VFA_IDS)]
    note: >
      Duplicates data/vfa_fermentation.yaml's vfa_fermentation.vfa_IDs list
      (now the actual code default via _VFA_FERM["vfa_IDs"] in
      _ad_fermentation_system.py) and analyses/msp_comparison.py's own
      VFA_IDS constant. Three independent copies of the same list; only the
      yaml one is load-bearing for simulation now.
```

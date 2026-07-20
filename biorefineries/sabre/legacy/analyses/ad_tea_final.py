# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Methanogenic AD TEA
=====================================
Runs all five pretreatment cases
  1. Base case (press_mill_only) with full diagnostics
  2. Pretreatment case comparison table (MSP vs case)
  3. Feed price sensitivity (same 4 scenarios as VFA TEA)
  4. Biomethane market price sensitivity (NPV at target $/MMBtu)
  5. Biostimulant price sensitivity

Product: biomethane
Co-product: biostimulant liquid concentrate (pressate concentrator output) --> not sold in base case
MSP reported as $/kg CH4 and $/MMBtu

Sources:
  - Feed price scenarios
  - Biomethane market: US natural gas Henry Hub ~$3-6/MMBtu typical range
  - AD maintenance: ADBC spreadsheet, $10/m3/yr
  - Biogas upgrading maintenance: IEA Bioenergy Task 37 (2014), 3.5%/yr
  - Biostimulant price: Market.us Seaweed Extract Biostimulant Market (2024)
      Global market $1.1B (2024), brown seaweed 64.7% share, liquid 76.8%.
      Wholesale bulk liquid concentrate: $0.50-2.00/kg.
  - Digestate disposal: heavy metal contamination (As, Cd, Pb) from
      pelagic Sargassum prevents land application as fertilizer.
"""

import sys
from pathlib import Path

import biosteam as bst

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.systems import create_ad_biomethane_system
from biorefineries.sabre._tea import create_tea, solve_biomethane_msp

OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# -------------------------
# Pretreatment cases
# -------------------------
PRETREATMENT_CASES = [
    "press_mill_only",
    "enzymatic",
    "peroxide",
    "combined_PE",
    "combined_PTE",
]

# -------------------------
# Feed price scenarios
# Source: Rodriguez-Martinez et al. 2023 (Sci Total Environ)
# -------------------------
FEED_PRICE_CASES = [
    ("tipping_fee",      -0.02),
    ("near_zero",         0.00),
    ("low_cost_collect",  0.02),
    ("beach_midpoint",    0.05),
]

# -------------------------
# Biomethane market prices for NPV sensitivity
# Basis: US Henry Hub $2-8/MMBtu typical range
# -------------------------
BIOMETHANE_MARKET_MMBTU = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0]

# -------------------------
# Disposal costs
# Digestate and solids --> modeled as waste
# -------------------------
LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG = -0.002  # $2/m3 wastewater disposal
SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG = -0.02   # standard biosolids

# -------------------------
# Biostimulant co-product revenue
# Pressate concentrator output: membrane concentrate of Alginate, Fucoidan,
# Mannitol, Protein — the bioactive compounds in brown seaweed biostimulants.
# Market: $1.1B global (2024), 10.5% CAGR; brown seaweed dominant.
# Wholesale bulk liquid: $0.50-2.00/kg.
# -------------------------
BIOSTIMULANT_PRICE_BASE_USD_PER_KG = 0.00  # base case with no biostimulant revenue
BIOSTIMULANT_SENSITIVITY_USD_PER_KG = [0.00, 0.50, 1.00, 2.00]


def _apply_stream_economics(
    sys,
    biostimulant_price: float = BIOSTIMULANT_PRICE_BASE_USD_PER_KG,
):
    """
    Assign prices to all outlet streams.
    Positive price -> revenue (reduces MSP).
    Negative price -> disposal cost (increases MSP).
    """
    annual_hours = 330.0 * 24.0
    summary = {}

    # Pressate permeate — zero-cost discharge (floating biorefinery assumption)
    for sid in ("pressate_permeate",):
        try:
            s = sys.flowsheet.stream[sid]
            s.price = 0.0
        except Exception:
            pass

    # Biostimulant product (BiostimulantEvaporator outlet -- the actual
    # system-boundary product stream; biostimulant_membrane_concentrate is
    # now an internal stream consumed by BiostimulantEvaporator, so pricing
    # it has no effect on TEA revenue)
    for sid in ("biostimulant_product",):
        try:
            s = sys.flowsheet.stream[sid]
            s.price = biostimulant_price
            summary[sid] = s.price * s.F_mass * annual_hours
        except Exception:
            pass

    # Solid digestate disposal
    for sid in ("soil_amendment",):
        try:
            s = sys.flowsheet.stream[sid]
            s.price = SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG
            summary[sid] = s.price * s.F_mass * annual_hours
        except Exception:
            pass

    # Liquid digestate disposal
    for sid in ("liquid_digestate",):
        try:
            s = sys.flowsheet.stream[sid]
            s.price = LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG
            summary[sid] = s.price * s.F_mass * annual_hours
        except Exception:
            pass

    return summary


def build_and_simulate(
    pretreatment_case: str,
    feed_price_per_kg_wet: float,
    feedstock_type: str = "pelagic",
):
    bst.main_flowsheet.clear()
    set_thermo()
    sys = create_ad_biomethane_system(
        feedstock=feedstock_type,
        pretreatment_case=pretreatment_case,
    )
    sys.feeds[0].price = feed_price_per_kg_wet
    sys.simulate()
    return sys

def print_unit_cost_breakdown(sys):
    annual_hours = 330.0 * 24.0
    print("\n=== UNIT COST BREAKDOWN ===")
    print(f"  {'Unit':<8} {'Purchase Cost ($)':>18} {'add_OPEX ($/yr)':>18} {'Utility ($/yr)':>16}")
    print("  " + "-" * 64)
    for u in sys.units:
        capex = sum(u.baseline_purchase_costs.values()) if u.baseline_purchase_costs else 0.0
        opex = sum(u.add_OPEX.values()) * annual_hours if hasattr(u, 'add_OPEX') and u.add_OPEX else 0.0
        utility = u.power_utility.cost * annual_hours if hasattr(u, 'power_utility') else 0.0
        if capex > 0 or opex > 0 or utility > 0:
            print(f"  {u.ID:<8} ${capex:>17,.0f} ${opex:>17,.0f} ${utility:>15,.0f}")

def run_unit_capex_opex_table():
    annual_hours = 330.0 * 24.0
    print("\n" + "=" * 90)
    print("UNIT CAPEX AND OPEX BY PRETREATMENT CASE (near_zero feed, biostimulant=$0)")
    print("=" * 90)

    for case in PRETREATMENT_CASES:
        sys = build_and_simulate(case, feed_price_per_kg_wet=0.00)
        _apply_stream_economics(sys, biostimulant_price=0.00)

        print(f"\n  --- {case} ---")
        print(f"  {'Unit':<8} {'Purchase Cost ($)':>18} {'add_OPEX ($/yr)':>18} {'Utility ($/yr)':>16} {'Total OPEX ($/yr)':>18}")
        print("  " + "-" * 82)

        for u in sys.units:
            capex   = sum(u.baseline_purchase_costs.values()) if u.baseline_purchase_costs else 0.0
            add_op  = sum(u.add_OPEX.values()) * annual_hours if hasattr(u, 'add_OPEX') and u.add_OPEX else 0.0
            utility = u.power_utility.cost * annual_hours if hasattr(u, 'power_utility') else 0.0
            total_op = add_op + utility
            if capex > 0 or total_op > 0:
                print(f"  {u.ID:<8} ${capex:>17,.0f} ${add_op:>17,.0f} ${utility:>15,.0f} ${total_op:>17,.0f}")

def run_case(
    pretreatment_case: str,
    feed_price_per_kg_wet: float,
    case_label: str = "",
    run_diagnostics: bool = True,
    silent: bool = False,
    biostimulant_price: float = BIOSTIMULANT_PRICE_BASE_USD_PER_KG,
    target_biomethane_prices_mmbtu: list | None = None,
):
    sys = build_and_simulate(pretreatment_case, feed_price_per_kg_wet)

    stream_economics = _apply_stream_economics(sys, biostimulant_price)

    tea = create_tea(sys)
    msp = solve_biomethane_msp(tea, sys.flowsheet.stream.biomethane)

    # NPV at target market prices
    npv_results = {}
    if target_biomethane_prices_mmbtu:
        biomethane = sys.flowsheet.stream.biomethane
        ch4_mmbtu_per_kg = 0.0526
        ch4_mass = float(biomethane.imass["Methane"])
        total_mass = float(biomethane.F_mass)
        ch4_frac = ch4_mass / total_mass if total_mass > 0 else 0.0

        for target_mmbtu in target_biomethane_prices_mmbtu:
            target_usd_per_kg_stream = target_mmbtu * ch4_mmbtu_per_kg * ch4_frac
            old_price = biomethane.price
            biomethane.price = target_usd_per_kg_stream
            try:
                npv_results[target_mmbtu] = tea.NPV
            finally:
                biomethane.price = old_price

    if silent:
        return tea, msp, sys, npv_results

    # ======================================================
    # PRINT RESULTS
    # ======================================================
    label_str = f" [{case_label}]" if case_label else ""
    print("\n" + "=" * 65)
    print(f"FEED PRICE = {feed_price_per_kg_wet:.4f} $/kg wet Sargassum{label_str}")
    print(f"Pretreatment: {pretreatment_case} | Biostimulant: ${biostimulant_price:.2f}/kg")
    print("=" * 65)

    if run_diagnostics:
        print("\n" + "-" * 65)
        print("DIAGNOSTIC: AD MASS BALANCE & DESIGN")
        print("-" * 65)

        feed = sys.feeds[0]
        print(f"\n  [sargassum_feed]")
        print(f"    Total F_mass:  {feed.F_mass:>12.2f} kg/hr")
        try:
            dry = feed.F_mass - float(feed.imass["Water"])
            print(f"    Dry mass:      {dry:>12.2f} kg/hr")
            print(f"    Moisture:      {float(feed.imass['Water'])/feed.F_mass*100:>11.1f}%")
        except Exception:
            pass

        try:
            ad = sys.flowsheet.unit.AD
            print(f"\n  [AD Design Results]")
            for k in [
                "AD inlet TS (kg/hr)", "AD inlet VS (kg/hr)",
                "Biodegradable pool (kg/hr)", "Biodegradable destroyed (kg/hr)",
                "Methane production (kg/hr)", "CH4 yield (kg/kg VS fed)",
                "Total digester volume (m3)", "Number of digesters",
                "Digester volume each (m3)", "Heating duty (kJ/h)",
            ]:
                v = ad.design_results.get(k)
                if v is not None:
                    print(f"    {k:<45} {v:.3f}")
        except Exception as e:
            print(f"  Could not read AD: {e}")

        try:
            up = sys.flowsheet.unit.UP
            bm = sys.flowsheet.stream.biomethane
            print(f"\n  [Biomethane]")
            print(f"    Total F_mass:      {bm.F_mass:>10.2f} kg/hr")
            print(f"    Methane:           {float(bm.imass['Methane']):>10.2f} kg/hr")
            print(f"    CH4 mol%:          {up.design_results.get('Biomethane CH4 mol%', 0):.1f}%")
            print(f"    Raw biogas Nm3/h:  {up.design_results.get('Raw biogas flow (Nm3/h)', 0):.0f}")
        except Exception:
            pass

        try:
            pc = sys.flowsheet.unit.PC
            conc = pc.outs[0]
            print(f"\n  [Biostimulant concentrate (PC output)]")
            print(f"    Total F_mass:   {conc.F_mass:>10.2f} kg/hr")
            for cid in ["Alginate", "Fucoidan", "Mannitol", "Protein", "OtherSolids", "Water"]:
                if cid in conc.chemicals:
                    m = float(conc.imass[cid])
                    if m > 0:
                        print(f"    {cid:<20} {m:>10.2f} kg/hr")
            annual_revenue = biostimulant_price * conc.F_mass * 330 * 24
            print(f"    Annual revenue @ ${biostimulant_price:.2f}/kg: ${annual_revenue/1e6:.2f}M/yr")
        except Exception:
            pass

        try:
            sa = sys.flowsheet.stream.soil_amendment
            ld = sys.flowsheet.stream.liquid_digestate
            print(f"\n  [Digestate]")
            print(f"    soil_amendment:   {sa.F_mass:>10.2f} kg/hr")
            print(f"    liquid_digestate: {ld.F_mass:>10.2f} kg/hr")
            print(f"    Note: heavy metal content prevents fertilizer credit")
        except Exception:
            pass

        print("\n" + "-" * 65)

    # Stream economics summary
    print("\n=== STREAM ECONOMICS (revenue +, cost −) ===")
    for name, val in stream_economics.items():
        tag = "revenue" if val > 0 else "cost"
        print(f"  {name:<40} ${val/1e6:>+8.2f}M/yr  ({tag})")

    print("\n=== FEEDS SEEN BY SYSTEM ===")
    for s in sys.feeds:
        if s.F_mass == 0.0 and s.source is None and s.ID.startswith("s"):
            continue
        print(f"  {s.ID:<30} {round(s.F_mass, 3):>12.3f} kg/hr   price={s.price:.4f} $/kg")

    print("\n=== TEA SUMMARY ===")
    print(f"  TCI ($):              {tea.TCI:>18,.0f}")
    print(f"  FCI ($):              {tea.FCI:>18,.0f}")
    print(f"  VOC ($/yr):           {tea.VOC:>18,.0f}")
    print(f"  FOC ($/yr):           {tea.FOC:>18,.0f}")
    print(f"  Material cost ($/yr): {tea.material_cost:>18,.0f}")
    print(f"  Utility cost ($/yr):  {tea.utility_cost:>18,.0f}")

    print("\n=== MSP ===")
    for k, v in msp.items():
        print(f"  {k:<35} {v}")

    if npv_results:
        print("\n=== NPV AT TARGET BIOMETHANE PRICES ===")
        print(f"  {'Market ($/MMBtu)':<25} {'NPV ($M)':<20} {'Viable?'}")
        print("  " + "-" * 50)
        for price, npv in npv_results.items():
            viable = "✓" if npv > 0 else "✗"
            print(f"  ${price:<24.1f} ${npv/1e6:<19.2f} {viable}")

    print("\n=== TOP UNIT PURCHASE COSTS ===")
    purchase = [
        (u.ID, sum(u.baseline_purchase_costs.values()) if u.baseline_purchase_costs else 0.0)
        for u in sys.units
    ]
    purchase.sort(key=lambda x: x[1], reverse=True)
    for uid, cost in purchase[:12]:
        print(f"  {uid:<20} ${cost:>15,.0f}")

    return tea, msp, sys, npv_results


# =============================================================
# Comparison and sensitivity tables
# =============================================================

def run_pretreatment_comparison(
    feed_price: float = 0.02,
    biostimulant_price: float = BIOSTIMULANT_PRICE_BASE_USD_PER_KG,
):
    print("\n" + "=" * 80)
    print(
        f"PRETREATMENT CASE COMPARISON — feed=${feed_price:.3f}/kg  "
        f"biostimulant=${biostimulant_price:.2f}/kg"
    )
    print("=" * 80)
    print(
        f"  {'Case':<20} {'MSP ($/MMBtu)':<18} {'MSP ($/kg CH4)':<18} "
        f"{'TCI ($M)':<12} {'VOC ($M/yr)'}"
    )
    print("  " + "-" * 78)

    for case in PRETREATMENT_CASES:
        tea, msp, sys, _ = run_case(
            pretreatment_case=case,
            feed_price_per_kg_wet=feed_price,
            biostimulant_price=biostimulant_price,
            silent=True,
        )
        msp_mmbtu = msp.get("usd_per_mmbtu", float("nan"))
        msp_ch4 = msp.get("usd_per_kg_ch4", float("nan"))
        print(
            f"  {case:<20} ${msp_mmbtu:<17.3f} ${msp_ch4:<17.4f} "
            f"${tea.TCI/1e6:<11.1f} ${tea.VOC/1e6:.1f}M"
        )


def run_feed_price_sensitivity(pretreatment_case: str = "press_mill_only"):
    print("\n" + "=" * 65)
    print(f"FEED PRICE SENSITIVITY — pretreatment: {pretreatment_case}")
    print("=" * 65)
    print(
        f"  {'Scenario':<20} {'Feed ($/kg)':<14} "
        f"{'MSP ($/MMBtu)':<18} {'MSP ($/kg CH4)'}"
    )
    print("  " + "-" * 65)

    for label, price in FEED_PRICE_CASES:
        tea, msp, sys, _ = run_case(
            pretreatment_case=pretreatment_case,
            feed_price_per_kg_wet=price,
            silent=True,
        )
        msp_mmbtu = msp.get("usd_per_mmbtu", float("nan"))
        msp_ch4 = msp.get("usd_per_kg_ch4", float("nan"))
        print(
            f"  {label:<20} ${price:<13.3f} "
            f"${msp_mmbtu:<17.3f} ${msp_ch4:.4f}"
        )


def run_biomethane_npv_sensitivity(
    pretreatment_case: str = "press_mill_only",
    feed_price: float = 0.02,
):
    print("\n" + "=" * 65)
    print(f"BIOMETHANE NPV — pretreatment: {pretreatment_case} | feed=${feed_price:.3f}/kg")
    print("=" * 65)

    tea, msp, sys, npv_results = run_case(
        pretreatment_case=pretreatment_case,
        feed_price_per_kg_wet=feed_price,
        silent=True,
        target_biomethane_prices_mmbtu=BIOMETHANE_MARKET_MMBTU,
    )
    msp_mmbtu = msp.get("usd_per_mmbtu", float("nan"))
    print(f"  MSP: ${msp_mmbtu:.2f}/MMBtu\n")
    print(f"  {'Market ($/MMBtu)':<22} {'NPV ($M)':<20} {'Viable?'}")
    print("  " + "-" * 48)
    for price, npv in npv_results.items():
        viable = "✓ profitable" if npv > 0 else "✗ loss"
        print(f"  ${price:<21.1f} ${npv/1e6:<19.2f} {viable}")


def run_biostimulant_sensitivity(
    pretreatment_case: str = "press_mill_only",
    feed_price: float = 0.02,
):
    """
    Show impact of biostimulant co-product price on biomethane MSP.
    This is a key sensitivity because biostimulant revenue directly offsets
    the AD capital burden, potentially making the pathway more competitive.
    """
    print("\n" + "=" * 65)
    print(
        f"BIOSTIMULANT PRICE SENSITIVITY — "
        f"pretreatment: {pretreatment_case} | feed=${feed_price:.3f}/kg"
    )
    print("=" * 65)
    print(
        f"  {'Biostimulant ($/kg)':<24} {'MSP ($/MMBtu)':<20} "
        f"{'Biostimulant rev ($M/yr)'}"
    )
    print("  " + "-" * 65)

    for bs_price in BIOSTIMULANT_SENSITIVITY_USD_PER_KG:
        tea, msp, sys, _ = run_case(
            pretreatment_case=pretreatment_case,
            feed_price_per_kg_wet=feed_price,
            biostimulant_price=bs_price,
            silent=True,
        )
        msp_mmbtu = msp.get("usd_per_mmbtu", float("nan"))
        # Get biostimulant stream flow
        try:
            conc = sys.flowsheet.stream["biostimulant_product"]
            bs_revenue = bs_price * conc.F_mass * 330 * 24
        except Exception:
            bs_revenue = 0.0
        label = "no revenue" if bs_price == 0 else f"${bs_price:.2f}/kg"
        print(
            f"  {label:<24} ${msp_mmbtu:<19.3f} ${bs_revenue/1e6:.2f}M"
        )


# =============================================================
# Main
# =============================================================

if __name__ == "__main__":

    # 1. Base case with full diagnostics
    print("\n>>> BASE CASE: press_mill_only | base-case | biostimulant=$0.00/kg")
    run_case(
        pretreatment_case="press_mill_only",
        feed_price_per_kg_wet=0.02,
        case_label="near_zero",
        run_diagnostics=True,
        target_biomethane_prices_mmbtu=BIOMETHANE_MARKET_MMBTU,
    )

    print("\n>>> BASE CASE: press_mill_only | tipping_fee | biostimulant=$0.00/kg")
    run_case(
        pretreatment_case="press_mill_only",
        feed_price_per_kg_wet=-0.02,
        case_label="tipping_fee",
        run_diagnostics=False,
        target_biomethane_prices_mmbtu=BIOMETHANE_MARKET_MMBTU,
    )

    # 2. Pretreatment comparison (base + optimistic biostimulant)
    run_pretreatment_comparison(feed_price=0.02, biostimulant_price=0.50)
    run_pretreatment_comparison(feed_price=0.02, biostimulant_price=1.00)

    # 3. Feed price sensitivity
    run_feed_price_sensitivity(pretreatment_case="press_mill_only")

    # 4. Biomethane NPV sensitivity for each pretreatment case
    for case in PRETREATMENT_CASES:
        run_biomethane_npv_sensitivity(pretreatment_case=case, feed_price=0.02)

    # 5. Biostimulant price sensitivity
    run_biostimulant_sensitivity(pretreatment_case="press_mill_only", feed_price=0.02)
    run_biostimulant_sensitivity(pretreatment_case="combined_PE", feed_price=0.02)

    # 6. Unit CAPEX and OPEX by pretreatment case
    run_unit_capex_opex_table()

    print("\n>>> COMBINED_PE | near_zero | biostimulant=$0.00/kg — FULL DIAGNOSTICS")
    run_case(
        pretreatment_case="combined_PE",
        feed_price_per_kg_wet=0.00,
        case_label="near_zero",
        run_diagnostics=True,
        target_biomethane_prices_mmbtu=BIOMETHANE_MARKET_MMBTU,
    )

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


PRETREATMENT_LABELS = {
    "press_mill_only": "Press + mill only",
    "enzymatic": "Enzymatic",
    "peroxide": "Peroxide",
    "combined_PE": "Combined PE",
    "combined_PTE": "Combined PTE",
}

FEED_LABELS = {
    "tipping_fee": "Tipping fee\n(-$0.02/kg)",
    "near_zero": "Near-zero\n($0.00/kg)",
    "low_cost_collect": "Low-cost collect\n($0.02/kg)",
    "beach_midpoint": "Beach midpoint\n($0.05/kg)",
}


def collect_feed_pretreatment_results(
    biostimulant_price: float = 0.0,
):
    rows = []

    for pretreatment_case in PRETREATMENT_CASES:
        for feed_label, feed_price in FEED_PRICE_CASES:
            tea, msp, sys, _ = run_case(
                pretreatment_case=pretreatment_case,
                feed_price_per_kg_wet=feed_price,
                biostimulant_price=biostimulant_price,
                silent=True,
            )

            # Use AD methane production if available
            try:
                ch4_prod_kgph = float(
                    sys.flowsheet.unit.AD.design_results["Methane production (kg/hr)"]
                )
            except Exception:
                # fallback to biomethane methane flow
                ch4_prod_kgph = float(sys.flowsheet.stream.biomethane.imass["Methane"])

            rows.append({
                "pretreatment_case": pretreatment_case,
                "pretreatment_label": PRETREATMENT_LABELS.get(pretreatment_case, pretreatment_case),
                "feed_label": feed_label,
                "feed_label_pretty": FEED_LABELS.get(feed_label, feed_label),
                "feed_price_usd_per_kg_wet": feed_price,
                "msp_usd_per_mmbtu": float(msp.get("usd_per_mmbtu", np.nan)),
                "msp_usd_per_kg_ch4": float(msp.get("usd_per_kg_ch4", np.nan)),
                "methane_kgph": ch4_prod_kgph,
                "TCI_MUSD": tea.TCI / 1e6,
                "VOC_MUSD_per_yr": tea.VOC / 1e6,
            })

    return pd.DataFrame(rows)


def _plot_annotated_heatmap(
    data: pd.DataFrame,
    value_col: str,
    title: str,
    cbar_label: str,
    outfile: str,
    fmt: str = "{:.2f}",
):
    pivot = data.pivot(
        index="pretreatment_label",
        columns="feed_label_pretty",
        values=value_col,
    )

    # Keep desired order
    pretreatment_order = [PRETREATMENT_LABELS[c] for c in PRETREATMENT_CASES]
    feed_order = [FEED_LABELS[c[0]] for c in FEED_PRICE_CASES]
    pivot = pivot.reindex(index=pretreatment_order, columns=feed_order)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    im = ax.imshow(pivot.values, aspect="auto")

    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=0, ha="center")
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            ax.text(j, i, fmt.format(val), ha="center", va="center", fontsize=9)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    ax.set_title(title, pad=12)
    ax.set_xlabel("Feed price scenario")
    ax.set_ylabel("Pretreatment case")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close()


def make_feed_pretreatment_figures(
    biostimulant_price: float = 0.0,
):
    df = collect_feed_pretreatment_results(
        biostimulant_price=biostimulant_price
    )

    _plot_annotated_heatmap(
        data=df,
        value_col="msp_usd_per_mmbtu",
        title="Biomethane MSP across feed price and pretreatment scenarios",
        cbar_label="MSP ($/MMBtu)",
        outfile=str(OUT / "fig_feed_pretreatment_msp_heatmap.png"),
        fmt="{:.2f}",
    )

    _plot_annotated_heatmap(
        data=df,
        value_col="methane_kgph",
        title="Methane production across feed price and pretreatment scenarios",
        cbar_label="Methane production (kg/h)",
        outfile=str(OUT / "fig_feed_pretreatment_ch4_heatmap.png"),
        fmt="{:.0f}",
    )

    # Also save the underlying table
    df.to_csv(OUT / "feed_pretreatment_results.csv", index=False)

    print("\nSaved:")
    print("  fig_feed_pretreatment_msp_heatmap.png")
    print("  fig_feed_pretreatment_ch4_heatmap.png")
    print("  feed_pretreatment_results.csv")

if __name__ == "__main__":
    make_feed_pretreatment_figures(biostimulant_price=0.0)
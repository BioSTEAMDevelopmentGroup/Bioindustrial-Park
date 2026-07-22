# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
plot_vfa_pathway_figures.py
---------------------------
Figures:
1) VFA composition in acidogenic broth (kg/hr) --> bar chart of mass flow of each VFA in acidogenic broth
2) Crude microbial oil MSP vs feed price ($/kg oil) --> MSP vs benchmark ($5)
3) Product scenario comparison at near-zero feed price ($/kg product)

"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import biosteam as bst

# Make sure imports work regardless of the invocation cwd or sys.path state
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
for _p in (REPO_ROOT, SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from vfa_fermentation_tea import (
    VFA_IDS,
    FEED_PRICE_CASES,
    PRODUCT_SCENARIOS,
    build_and_simulate,
    build_and_simulate_scenario,
    run_case,
    _patch_ev607,
    _apply_disposal_costs,
)
from biorefineries.sabre._tea import create_tea, solve_product_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# -------------------------
# Plot style
# -------------------------
plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        10,
    "axes.titlesize":   10,
    "axes.labelsize":   10,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "figure.dpi":       150,
    "axes.linewidth":   0.8,
    "axes.edgecolor":   "black",
    "xtick.direction":  "in",
    "ytick.direction":  "in",
    "xtick.top":        True,
    "ytick.right":      True,
})


# ============================================================
# Figure 1 — VFA composition in acidogenic broth
# ============================================================
def make_vfa_composition_figure(feed_price: float = 0.02):
    streams, units, full_sys = build_and_simulate(feed_price)
    vfa_broth = streams["vfa_broth"]

    vfa_names = {
        "AceticAcid": "Acetate",
        "PropionicAcid": "Propionate",
        "ButyricAcid": "Butyrate",
        "ValericAcid": "Valerate",
        "HexanoicAcid": "Hexanoate",
    }

    labels = []
    flows = []
    for cid in VFA_IDS:
        if cid in vfa_broth.chemicals:
            labels.append(vfa_names.get(cid, cid))
            flows.append(float(vfa_broth.imass[cid]))

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    bars = ax.bar(labels, flows, edgecolor="black", linewidth=0.8, zorder=3)

    for bar, val in zip(bars, flows):
        ax.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(flows) * 0.02,
            f"{val:,.0f}",
            ha="center", va="bottom", fontsize=8
        )

    ax.set_ylabel("VFA flow in acidogenic broth (kg/hr)")
    ax.set_title("VFA composition of acidogenic broth")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)

    fig.tight_layout()
    fig.savefig(OUT / "fig_vfa_composition.png", bbox_inches="tight")

    print("\nVFA composition (kg/hr):")
    for lab, val in zip(labels, flows):
        print(f"  {lab:<12} {val:,.2f}")


# ============================================================
# Figure 2 — Crude microbial oil MSP vs feed price
# ============================================================
def make_oil_msp_vs_feed_price_figure():
    feed_prices = [price for _, price in FEED_PRICE_CASES]
    feed_labels = {
        "tipping_fee": "Tipping fee",
        "near_zero": "Near-zero",
        "low_cost_collect": "Low-cost collect",
        "beach_midpoint": "Beach midpoint",
    }

    msp_usd_per_kg_oil = []

    for label, price in FEED_PRICE_CASES:
        tea, msp, streams, units, full_sys = run_case(
            feed_price_per_kg_wet=price,
            case_label=label,
            run_diagnostics=False,
            silent=True,
        )
        msp_usd_per_kg_oil.append(float(msp["usd_per_kg"]))

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.plot(
        feed_prices,
        msp_usd_per_kg_oil,
        marker="o",
        linewidth=1.8,
        markersize=6,
        markeredgecolor="black",
        zorder=3,
    )

    for x, y, (label, _) in zip(feed_prices, msp_usd_per_kg_oil, FEED_PRICE_CASES):
        ax.text(
            x,
            y + max(msp_usd_per_kg_oil) * 0.015,
            f"{y:.2f}",
            ha="center", va="bottom", fontsize=8
        )

    # $5/kg as the microbial oil market reference
    ax.axhline(5.0, color="black", linewidth=1.0, linestyle="--", zorder=1, label="$5/kg market reference")

    ax.set_xlabel("Feed price ($/kg wet Sargassum)")
    ax.set_ylabel("Crude microbial oil MSP ($/kg)")
    ax.set_title("Effect of feed price on crude microbial oil MSP")
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(OUT / "fig_crude_oil_msp_vs_feed_price.png", bbox_inches="tight")

    print("\nCrude microbial oil MSP by feed price:")
    for (label, price), msp in zip(FEED_PRICE_CASES, msp_usd_per_kg_oil):
        print(f"  {label:<18} feed={price:+.2f} $/kg -> MSP={msp:.3f} $/kg")


# ============================================================
# Figure 3 — Product scenario comparison at near-zero feed price
# ============================================================
def make_product_scenario_comparison_figure(feed_price: float = 0.00):
    labels = []
    msp_vals = []
    annual_prod_t = []

    for sc in PRODUCT_SCENARIOS:
        streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=feed_price,
            product_yield=sc["yield"],
            residence_time_h=sc["residence_h"],
        )

        _patch_ev607(full_sys, silent=True)
        _apply_disposal_costs(streams)

        # Putting purification cost into OE.add_OPEX exactly as in vfa_fermentation_tea.py
        oil_stream = streams["microbial_oil"]
        product_kg_hr = float(oil_stream.imass["MicrobialOil"])
        extraction_usd_per_hr = product_kg_hr * sc["extraction_usd_per_kg"]

        try:
            oe_unit = bst.main_flowsheet.unit["OE"]
            oe_unit.add_OPEX = {"Product extraction/purification": extraction_usd_per_hr}
        except Exception:
            pass

        tea = create_tea(full_sys)

        msp_dict = solve_product_msp(
            tea=tea,
            product_stream=oil_stream,
        )

        labels.append(sc["label"])
        msp_vals.append(float(msp_dict["usd_per_kg"]))
        annual_prod_t.append(float(msp_dict["annual_product_kg"]) / 1000.0)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    bars = ax.bar(labels, msp_vals, edgecolor="black", linewidth=0.8, zorder=3)

    for bar, val in zip(bars, msp_vals):
        ax.text(
            bar.get_x() + bar.get_width()/2,
            bar.get_height() + max(msp_vals) * 0.015,
            f"{val:,.1f}",
            ha="center", va="bottom", fontsize=8
        )

    ax.set_ylabel("MSP ($/kg product)")
    ax.set_title("Product scenario comparison at near-zero feed price")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)

    fig.tight_layout()
    fig.savefig(OUT / "fig_product_scenario_comparison.png", bbox_inches="tight")

    print("\nProduct scenario comparison at near-zero feed:")
    for lab, msp, prod in zip(labels, msp_vals, annual_prod_t):
        print(f"  {lab:<18} MSP={msp:,.2f} $/kg | annual production={prod:,.1f} t/yr")


if __name__ == "__main__":
    make_vfa_composition_figure(feed_price=0.02)
    make_oil_msp_vs_feed_price_figure()
    make_product_scenario_comparison_figure(feed_price=0.02)

    print("\nSaved:")
    print(OUT / "fig_vfa_composition.png")
    print(OUT / "fig_crude_oil_msp_vs_feed_price.png")
    print(OUT / "fig_product_scenario_comparison.png")
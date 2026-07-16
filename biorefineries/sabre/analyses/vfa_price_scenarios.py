# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
vfa_price_scenarios.py
---------------------------

Figures:
1) fig_vfa_feed_price.png
   - crude microbial oil MSP vs feed price

2) fig_vfa_biostimulant_price.png
   - crude microbial oil MSP vs assumed biostimulant price

3) fig_vfa_product_scenarios.png
   - MSP by product scenario, with market-price markers
   - initial scoping assumptions

"""

from __future__ import annotations

from pathlib import Path
import sys

import biosteam as bst
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
for _p in (REPO_ROOT, SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from vfa_fermentation_tea import (
    FEED_PRICE_CASES,
    PRODUCT_SCENARIOS,
    OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL,
    SOLIDS_DISPOSAL_USD_PER_KG,
    build_and_simulate,
    build_and_simulate_scenario,
    run_case,
    _patch_ev607,
    _apply_disposal_costs,
)
from biorefineries.sabre._tea import make_baseline_tea, solve_product_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Base assumptions ──────────────────────────────────────────────────────────
BIOSTIMULANT_PRICE_CASES = [0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00]
FEED_PRICE_BASE = 0.02

# Soybean oil market price range ($/kg) — used as benchmark band
SOYBEAN_OIL_MARKET_LOW_USD_PER_KG  = 0.62
SOYBEAN_OIL_MARKET_HIGH_USD_PER_KG = 1.50

# Product scenario colors — consistent with matplotlib default cycle
PRODUCT_COLORS = [
    "#1f77b4",   # blue   — Microbial oil
    "#ff7f0e",   # orange — Omega-3 EPA oil
    "#2ca02c",   # green  — Astaxanthin
]

# ── Shared plot style ─────────────────────────────────────────────────────────
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


# ── Helpers ───────────────────────────────────────────────────────────────────

def _solve_msp_from_system(
    full_sys,
    streams,
    extraction_usd_per_kg_oil: float = OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL,
    solids_disposal_usd_per_kg: float = SOLIDS_DISPOSAL_USD_PER_KG,
) -> float:
    """Solve crude microbial oil MSP on the same basis as run_case()."""
    _patch_ev607(full_sys, silent=True)
    _apply_disposal_costs(streams, solids_disposal_usd_per_kg=solids_disposal_usd_per_kg)

    oil_stream = streams["backend_oil"]
    oil_kg_hr = float(oil_stream.imass["MicrobialOil"])

    extraction_usd_per_hr = oil_kg_hr * extraction_usd_per_kg_oil
    try:
        oe_unit = bst.main_flowsheet.unit["OE"]
        oe_unit.add_OPEX = {"Oil extraction reagent": extraction_usd_per_hr}
    except Exception:
        pass

    tea = make_baseline_tea(full_sys)
    msp = solve_product_msp(
        tea=tea,
        product_stream=oil_stream,
        product_ID="MicrobialOil",
    )
    return float(msp["usd_per_kg_product"])


def _apply_biostimulant_price(streams, price_per_kg: float) -> bool:
    candidate_ids = [
        "biostimulant_membrane_concentrate",
        "biostimulant_concentrate",
        "pressate_concentrate",
    ]
    for sid in candidate_ids:
        try:
            bst.main_flowsheet.stream[sid].price = price_per_kg
            return True
        except Exception:
            pass
    if isinstance(streams, dict):
        for sid in candidate_ids:
            try:
                streams[sid].price = price_per_kg
                return True
            except Exception:
                pass
    return False


def _add_soybean_band(ax):
    """Add soybean oil market price band to an axes."""
    ax.axhline(
        SOYBEAN_OIL_MARKET_HIGH_USD_PER_KG,
        color="black", linewidth=1.0, linestyle="--", zorder=2,
        label="Soybean oil market price ($0.62–1.50/kg)",
    )
    ax.axhline(
        SOYBEAN_OIL_MARKET_LOW_USD_PER_KG,
        color="black", linewidth=1.0, linestyle="--", zorder=2,
    )
    ax.axhspan(
        SOYBEAN_OIL_MARKET_LOW_USD_PER_KG,
        SOYBEAN_OIL_MARKET_HIGH_USD_PER_KG,
        alpha=0.10, color="black", zorder=1,
    )


# ── Figure 1 — Feed Price ─────────────────────────────────────────────────────
def make_feed_price_figure():
    feed_prices = [price for _, price in FEED_PRICE_CASES]
    msp_vals = []

    for label, price in FEED_PRICE_CASES:
        tea, msp, streams, units, full_sys = run_case(
            feed_price_per_kg_wet=price,
            case_label=label,
            run_diagnostics=False,
            silent=True,
        )
        msp_vals.append(float(msp["usd_per_kg_product"]))

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.plot(
        feed_prices,
        msp_vals,
        marker="o",
        linewidth=1.8,
        markersize=6,
        markeredgecolor="black",
        zorder=3,
    )

    for x, y in zip(feed_prices, msp_vals):
        ax.text(
            x,
            y + max(msp_vals) * 0.015,
            f"{y:.2f}",
            ha="center", va="bottom",
            fontsize=8,
        )

    # Base case star at $0.02/kg
    base_feed_price = 0.02
    base_idx = feed_prices.index(base_feed_price)
    ax.scatter(
        [base_feed_price],
        [msp_vals[base_idx]],
        marker="*", s=220,
        color="gold", edgecolors="black", linewidths=0.8,
        zorder=5, label="Base case",
    )

    _add_soybean_band(ax)

    ax.set_xlabel("Feed price ($/kg wet Sargassum)")
    ax.set_ylabel("Crude microbial oil MSP ($/kg)")
    ax.set_title("Effect of feed price on crude microbial oil MSP")
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_ylim(0, max(msp_vals) * 1.15)

    fig.tight_layout()
    fig.savefig(OUT / "fig_vfa_feed_price.png", bbox_inches="tight")
    print("Saved: fig_vfa_feed_price.png")


# ── Figure 2 — Biostimulant Price ────────────────────────────────────────────
def make_biostimulant_price_figure():
    msp_vals = []

    for bs_price in BIOSTIMULANT_PRICE_CASES:
        streams, units, full_sys = build_and_simulate(FEED_PRICE_BASE)
        _apply_biostimulant_price(streams, bs_price)
        msp_val = _solve_msp_from_system(full_sys, streams)
        msp_vals.append(msp_val)

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    ax.plot(
        BIOSTIMULANT_PRICE_CASES,
        msp_vals,
        marker="o",
        linewidth=1.8,
        markersize=6,
        markeredgecolor="black",
        zorder=3,
    )

    for x, y in zip(BIOSTIMULANT_PRICE_CASES, msp_vals):
        ax.text(
            x,
            y + max(msp_vals) * 0.015,
            f"{y:.2f}",
            ha="center", va="bottom",
            fontsize=8,
        )

    # Base case star at $0.00/kg biostimulant
    base_idx = BIOSTIMULANT_PRICE_CASES.index(0.00)
    ax.scatter(
        [0.00],
        [msp_vals[base_idx]],
        marker="*", s=220,
        color="gold", edgecolors="black", linewidths=0.8,
        zorder=5, label="Base case",
    )

    _add_soybean_band(ax)

    ax.set_xlabel("Biostimulant price ($/kg)")
    ax.set_ylabel("Crude microbial oil MSP ($/kg)")
    ax.set_title("Effect of biostimulant price on crude microbial oil MSP")
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_ylim(0, max(msp_vals) * 1.15)

    fig.tight_layout()
    fig.savefig(OUT / "fig_vfa_biostimulant_price.png", bbox_inches="tight")
    print("Saved: fig_vfa_biostimulant_price.png")


# ── Figure 3 — Product Scenarios ─────────────────────────────────────────────
def make_product_scenarios_figure():
    labels = []
    msp_vals = []
    market_vals = []

    for sc in PRODUCT_SCENARIOS:
        streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=FEED_PRICE_BASE,
            product_yield=sc["yield"],
            residence_time_h=sc["residence_h"],
        )

        _patch_ev607(full_sys, silent=True)
        _apply_disposal_costs(streams)

        oil_stream = streams["backend_oil"]
        product_kg_hr = float(oil_stream.imass["MicrobialOil"])
        extraction_usd_per_hr = product_kg_hr * sc["extraction_usd_per_kg"]

        try:
            bst.main_flowsheet.unit["OE"].add_OPEX = {
                "Product extraction/purification": extraction_usd_per_hr
            }
        except Exception:
            pass

        tea = make_baseline_tea(full_sys)
        msp_dict = solve_product_msp(
            tea=tea,
            product_stream=oil_stream,
            product_ID="MicrobialOil",
        )

        labels.append(sc["label"])
        msp_vals.append(float(msp_dict["usd_per_kg_product"]))
        market_vals.append(float(sc["market_price"]))

    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(7.0, 4.8))

    # Colored bars — one color per product
    bars = ax.bar(
        x,
        msp_vals,
        color=PRODUCT_COLORS[:len(labels)],
        edgecolor="black",
        linewidth=0.8,
        width=0.55,
        zorder=2,
    )

    # MSP label on each bar, offset above to clear the diamond
    for bar, val in zip(bars, msp_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            val * 1.6,
            f"${val:,.2f}/kg",
            ha="center", va="bottom",
            fontsize=8,
        )

    # Market reference diamonds
    ax.scatter(
        x,
        market_vals,
        marker="D",
        s=65,
        color="gold",
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Price ($/kg product)")
    ax.set_title("Product scenario comparison at near-zero feed price")

    # Log scale to show all three products on one axis
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"${v:,.0f}"))
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0, which="both")
    ax.set_ylim(bottom=0, top=max(msp_vals) * 5)

    # Legend: colored patch per product MSP + market reference diamonds
    legend_handles = [
        Patch(facecolor=PRODUCT_COLORS[i], edgecolor="black",
              label=f"{labels[i]} MSP")
        for i in range(len(labels))
    ] + [
        Line2D([0], [0], marker="D", color="w",
               markerfacecolor="gold", markeredgecolor="black",
               markersize=8,
               label=f"{lab}: ${mv:,.0f}/kg market")
        for lab, mv in zip(labels, market_vals)
    ]
    ax.legend(handles=legend_handles, frameon=False, fontsize=8, loc="upper left")

    fig.tight_layout()
    fig.savefig(OUT / "fig_vfa_product_scenarios.png", bbox_inches="tight")
    print("Saved: fig_vfa_product_scenarios.png")

    print("\nProduct scenario results:")
    for lab, msp, mkt in zip(labels, msp_vals, market_vals):
        print(f"  {lab:<18} MSP={msp:,.2f} $/kg | market ref={mkt:,.2f} $/kg")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    make_feed_price_figure()
    make_biostimulant_price_figure()
    make_product_scenarios_figure()

    print("\nSaved:")
    print(OUT / "fig_vfa_feed_price.png")
    print(OUT / "fig_vfa_biostimulant_price.png")
    print(OUT / "fig_vfa_product_scenarios.png")
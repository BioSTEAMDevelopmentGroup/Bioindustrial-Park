# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
plot_vfa_yield_sensitivities.py --> see how yield effects the VFA and oil production
-------------------------------

1) Crude microbial oil MSP vs acidogenic VFA yield
   (kg VFA / kg VS destroyed)

2) Crude microbial oil MSP vs fermentation oil yield
   (kg oil / kg VFA consumed)

"""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import biosteam as bst

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
for _p in (REPO_ROOT, SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from vfa_fermentation_tea import (
    build_and_simulate,
    build_and_simulate_scenario,
    _patch_ev607,
    _apply_disposal_costs,
)
from biorefineries.sabre._tea import make_baseline_tea, solve_product_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# -------------------------
# Base assumptions / ranges
# -------------------------
FEED_PRICE_USD_PER_KG_WET = 0.02

# From assumptions.yaml fit_note: sensitivity range 0.35–0.65
VFA_YIELD_RANGE = np.array([0.35, 0.40, 0.47, 0.55, 0.65])  # kg VFA / kg VS destroyed

# Base fermentation yield is 0.15 kg/kg VFA consumed.
# Use a conservative-to-optimistic scoping range.
FERM_YIELD_RANGE = np.array([0.10, 0.144, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])  # kg oil / kg VFA consumed

BASE_FERM_YIELD = 0.144
BASE_RESIDENCE_H = 48.0

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


def solve_oil_msp_from_full_system(
    full_sys,
    streams,
    extraction_usd_per_kg_oil: float = 1.00,
):
    _patch_ev607(full_sys, silent=True)
    _apply_disposal_costs(streams)

    oil_stream = streams["backend_oil"]
    oil_kg_hr = float(oil_stream.imass["MicrobialOil"])

    # match run_case() logic in vfa_fermentation_tea.py
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


def override_acidogenic_vfa_yield(sys, target_vfa_yield):
    """
    Find the acidogenic digester unit generically and override
    vfa_kg_per_kg_vs, then re-simulate the (single, unified) system.
    """
    found = False
    for u in sys.units:
        if hasattr(u, "vfa_kg_per_kg_vs"):
            u.vfa_kg_per_kg_vs = float(target_vfa_yield)
            found = True
    if not found:
        raise RuntimeError("Could not find acidogenic digester unit with attribute 'vfa_kg_per_kg_vs'.")
    sys.simulate()


# ============================================================
if __name__ == "__main__":
    # Figure 1 — MSP vs acidogenic VFA yield
    # ============================================================
    msp_vs_vfa_yield = []

    for vfa_yield in VFA_YIELD_RANGE:
        # Build normal base case first
        streams, units, full_sys = build_and_simulate(FEED_PRICE_USD_PER_KG_WET)

        # Override acidogenic yield and re-simulate the (single, unified) system
        override_acidogenic_vfa_yield(full_sys, vfa_yield)

        msp_val = solve_oil_msp_from_full_system(full_sys, streams, extraction_usd_per_kg_oil=0.50)
        msp_vs_vfa_yield.append(msp_val)

    fig1, ax1 = plt.subplots(figsize=(6.6, 4.4))
    ax1.plot(
        VFA_YIELD_RANGE,
        msp_vs_vfa_yield,
        marker="o",
        linewidth=1.8,
        markersize=6,
        markeredgecolor="black",
        zorder=3,
    )

    ax1.axhline(
        y=1.50,
        linestyle="--",
        linewidth=1.2,
        color="black",
        zorder=2,
        label="Soybean oil market price ($0.62–1.50/kg)",
    )
    ax1.axhline(
        y=0.62,
        linestyle="--",
        linewidth=1.2,
        color="black",
        zorder=2,
    )
    ax1.axhspan(0.62, 1.50, alpha=0.10, color="black", zorder=1)

    for x, y in zip(VFA_YIELD_RANGE, msp_vs_vfa_yield):
        ax1.text(
            x,
            y + max(msp_vs_vfa_yield) * 0.015,
            f"{y:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    # base-case star
    base_vfa_yield = 0.55
    base_idx = list(VFA_YIELD_RANGE).index(base_vfa_yield)
    base_msp_vfa = msp_vs_vfa_yield[base_idx]

    ax1.scatter(
        [base_vfa_yield],
        [base_msp_vfa],
        marker="*",
        s=220,
        color="gold",
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
        label="Base case",
    )

    ax1.set_xlabel("Acidogenic VFA yield (kg VFA/kg VS destroyed)")
    ax1.set_ylabel("Crude microbial oil MSP ($/kg)")
    ax1.set_title("Effect of acidogenic VFA yield on crude microbial oil MSP")
    ax1.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax1.set_ylim(0, max(msp_vs_vfa_yield) * 1.10)
    ax1.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax1.legend(frameon=False, fontsize=8, loc="upper right")

    fig1.tight_layout()
    fig1.savefig(OUT / "fig_msp_vs_acidogenic_vfa_yield.png", bbox_inches="tight")


    # ============================================================
    # Figure 2 — MSP vs fermentation oil yield
    # ============================================================
    msp_vs_ferm_yield = []

    for ferm_yield in FERM_YIELD_RANGE:
        streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=FEED_PRICE_USD_PER_KG_WET,
            product_yield=float(ferm_yield),
            residence_time_h=BASE_RESIDENCE_H,
        )

        msp_val = solve_oil_msp_from_full_system(
            full_sys,
            streams,
            extraction_usd_per_kg_oil=0.50,
            )
        msp_vs_ferm_yield.append(msp_val)

    fig2, ax2 = plt.subplots(figsize=(6.6, 4.4))
    ax2.plot(
        FERM_YIELD_RANGE,
        msp_vs_ferm_yield,
        marker="o",
        linewidth=1.8,
        markersize=6,
        markeredgecolor="black",
        zorder=3,
    )

    ax2.axhline(
        y=1.50,
        linestyle="--",
        linewidth=1.2,
        color="black",
        zorder=2,
        label="Soybean oil market price ($0.62–1.50/kg)",
    )
    ax2.axhline(
        y=0.62,
        linestyle="--",
        linewidth=1.2,
        color="black",
        zorder=2,
    )
    ax2.axhspan(0.62, 1.50, alpha=0.10, color="black", zorder=1)

    for x, y in zip(FERM_YIELD_RANGE, msp_vs_ferm_yield):
        ax2.text(
            x,
            y + max(msp_vs_ferm_yield) * 0.015,
            f"{y:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    # base-case star
    base_ferm_yield = 0.144
    base_idx = list(FERM_YIELD_RANGE).index(base_ferm_yield)
    base_msp_ferm = msp_vs_ferm_yield[base_idx]

    ax2.scatter(
        [base_ferm_yield],
        [base_msp_ferm],
        marker="*",
        s=220,
        color="gold",
        edgecolors="black",
        linewidths=0.8,
        zorder=5,
        label="Base case",
    )

    ax2.set_xlabel("Fermentation oil yield (kg oil/kg VFA consumed)")
    ax2.set_ylabel("Crude microbial oil MSP ($/kg)")
    ax2.set_title("Effect of fermentation oil yield on crude microbial oil MSP")
    ax2.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax2.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax2.set_ylim(0, max(msp_vs_ferm_yield) * 1.10)
    ax2.legend(frameon=False, fontsize=8, loc="upper right")

    fig2.tight_layout()
    fig2.savefig(OUT / "fig_msp_vs_fermentation_yield.png", bbox_inches="tight")

    print("\nMSP vs acidogenic VFA yield:")
    for x, y in zip(VFA_YIELD_RANGE, msp_vs_vfa_yield):
        print(f"  {x:.2f} kg VFA/kg VS destroyed -> ${y:.3f}/kg oil")

    print("\nMSP vs fermentation oil yield:")
    for x, y in zip(FERM_YIELD_RANGE, msp_vs_ferm_yield):
        print(f"  {x:.2f} kg oil/kg VFA consumed -> ${y:.3f}/kg oil")

    print("\nSaved:")
    print(OUT / "fig_msp_vs_acidogenic_vfa_yield.png")
    print(OUT / "fig_msp_vs_fermentation_yield.png")
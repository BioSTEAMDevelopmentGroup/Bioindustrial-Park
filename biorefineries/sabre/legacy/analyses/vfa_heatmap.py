# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
vfa_heatmap.py
-----------------------
2D heatmap of crude microbial oil MSP ($/kg) as a function of
feedstock price (x-axis) and biostimulant price (y-axis)
for the acidogenic AD–fermentation pathway.

Outputs:
  - fig_vfa_msp_heatmap.png   (single panel, main text)

Run from sabre_project root:
    python scripts/vfa_heatmap.py
"""

from pathlib import Path
import sys

import biosteam as bst
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
for _p in (REPO_ROOT, SCRIPT_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from vfa_fermentation_tea import (
    OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL,
    SOLIDS_DISPOSAL_USD_PER_KG,
    build_and_simulate,
    _patch_ev607,
    _apply_disposal_costs,
)
from biorefineries.sabre._tea import create_tea, solve_product_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Grid definition ───────────────────────────────────────────────────────────

FEED_PRICES    = [-0.02, 0.00, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10]  # $/kg wet
BIOSTIM_PRICES = [0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00]  # $/kg

# Soybean oil market price range for reference ($/kg)
SOYBEAN_OIL_LOW  = 0.62
SOYBEAN_OIL_HIGH = 1.50


# ── Helpers ───────────────────────────────────────────────────────────────────

def _apply_biostimulant_price(price_per_kg: float):
    # biostimulant_product (BiostimulantEvaporator outlet) is the actual
    # system-boundary product stream; biostimulant_membrane_concentrate is
    # now internal (consumed by BiostimulantEvaporator), so pricing it has
    # no effect on TEA revenue.
    try:
        bst.main_flowsheet.stream["biostimulant_product"].price = price_per_kg
    except Exception:
        pass


def build_case(feed_price: float, biostim_price: float) -> float:
    streams, units, full_sys = build_and_simulate(feed_price)

    _apply_biostimulant_price(biostim_price)
    _patch_ev607(full_sys, silent=True)
    _apply_disposal_costs(streams, solids_disposal_usd_per_kg=SOLIDS_DISPOSAL_USD_PER_KG)

    oil_stream = streams["microbial_oil"]
    oil_kg_hr = float(oil_stream.imass["MicrobialOil"])

    extraction_usd_per_hr = oil_kg_hr * OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL
    try:
        bst.main_flowsheet.unit["OE"].add_OPEX = {
            "Oil extraction reagent": extraction_usd_per_hr
        }
    except Exception:
        pass

    tea = create_tea(full_sys)
    msp = solve_product_msp(
        tea=tea,
        product_stream=oil_stream,
    )
    return float(msp["usd_per_kg_stream"])


# ── Run grid ──────────────────────────────────────────────────────────────────

# grid shape: (len(BIOSTIM_PRICES), len(FEED_PRICES))
# rows = biostimulant price (y-axis), cols = feed price (x-axis)
if __name__ == "__main__":
    print("Running VFA MSP heatmap grid...")
    grid = np.zeros((len(BIOSTIM_PRICES), len(FEED_PRICES)))

    for i, bs in enumerate(BIOSTIM_PRICES):
        for j, fp in enumerate(FEED_PRICES):
            msp = build_case(fp, bs)
            grid[i, j] = msp
            print(f"  feed={fp:+.2f}  bs={bs:.2f}  -> MSP=${msp:.2f}/kg")


    # ── Plot ──────────────────────────────────────────────────────────────────────

    plt.rcParams.update({
        "font.family":     "DejaVu Sans",
        "font.size":       10,
        "axes.titlesize":  11,
        "axes.labelsize":  10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.dpi":      150,
        "axes.linewidth":  0.8,
    })

    CMAP = "RdYlGn_r"
    VMIN = np.percentile(grid, 2)
    VMAX = np.percentile(grid, 98)

    n_rows, n_cols = grid.shape

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    im = ax.imshow(
        grid,
        aspect="auto",
        origin="lower",
        cmap=CMAP,
        vmin=VMIN,
        vmax=VMAX,
        extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5],
    )

    # Cell annotations
    for i in range(n_rows):
        for j in range(n_cols):
            val = grid[i, j]
            normed = (val - VMIN) / max(VMAX - VMIN, 1)
            color = "white" if normed < 0.35 or normed > 0.72 else "black"
            ax.text(
                j, i, f"${val:.0f}",
                ha="center", va="center",
                fontsize=7, color=color, fontweight="bold",
            )

    # Highlight cells within soybean oil market price range
    for i in range(n_rows):
        for j in range(n_cols):
            val = grid[i, j]
            if val <= SOYBEAN_OIL_HIGH:
                edgecolor = "#1a6faf" if val <= SOYBEAN_OIL_LOW else "black"
                rect = plt.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    linewidth=2.5, edgecolor=edgecolor,
                    facecolor="none", zorder=6,
                )
                ax.add_patch(rect)

    # Tick labels
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([f"${fp:.2f}" for fp in FEED_PRICES], fontsize=8)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f"${bs:.2f}" for bs in BIOSTIM_PRICES], fontsize=8)

    ax.set_xlabel("Feedstock price ($/kg wet)")
    ax.set_ylabel("Biostimulant price ($/kg)")
    ax.set_title("Crude microbial oil MSP ($/kg)\nAcidogenic AD–fermentation pathway", pad=6)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Crude microbial oil MSP ($/kg)", fontsize=9)

    # Legend for borders
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="none", edgecolor="black",   linewidth=2.0, label=f"MSP ≤ ${SOYBEAN_OIL_HIGH}/kg (soybean oil high)"),
        Patch(facecolor="none", edgecolor="#1a6faf", linewidth=2.0, label=f"MSP ≤ ${SOYBEAN_OIL_LOW}/kg (soybean oil low)"),
    ]
    fig.legend(handles=legend_handles, fontsize=8, frameon=False,
               loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.04))

    fig.tight_layout()
    fig.savefig(OUT / "fig_vfa_msp_heatmap.png", bbox_inches="tight")
    print("\nSaved: fig_vfa_msp_heatmap.png")

    # ── Console summary ───────────────────────────────────────────────────────────

    print("\nMSP grid — Crude microbial oil ($/kg):")
    print(f"{'':>14}", end="")
    for fp in FEED_PRICES:
        print(f"  feed={fp:+.2f}", end="")
    print()
    for i, bs in enumerate(BIOSTIM_PRICES):
        print(f"bs=${bs:.2f}/kg    ", end="")
        for j in range(len(FEED_PRICES)):
            print(f"  {grid[i, j]:9.2f}", end="")
        print()
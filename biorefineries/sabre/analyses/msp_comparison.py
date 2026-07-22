# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_comparison.py
------------------
Compare minimum selling price (MSP) across the four standalone SaBRe
flowsheets. The 'integrated' flowsheet (alpha-split between the
methanogenic and VFA-to-oil pathways) is intentionally excluded.

    - Biostimulant     -> biostimulant liquid product, $/kg product
    - AD-biomethane    -> biomethane, $/kg CH4
    - AD-VFA           -> VFA broth, $/kg total VFA
    - AD-fermentation  -> crude microbial oil, $/kg MicrobialOil

The product of each pathway is different, MSP of each pathway is the price
per kg of the pathway's own value-carrying product needed to hit the
target IRR, assuming no revenue from any other product.

All four cases share one baseline economic basis as set in tea.yaml

Run from the repo root:
    python biorefineries/sabre/analyses/msp_comparison.py
"""

import sys
from pathlib import Path

import biosteam as bst
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre.systems import (
    price_biostimulant_system,
    price_ad_biomethane_system,
    price_ad_vfa_system,
    price_ad_fermentation_system,
)

OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


# System builders -- each returns a common result dict
BUILDERS = (
    price_biostimulant_system,
    price_ad_biomethane_system,
    price_ad_vfa_system,
    price_ad_fermentation_system,
)


# =================================================================
# Main
# =================================================================

if __name__ == "__main__":
    results = [builder() for builder in BUILDERS]

    # ── Plot ──────────────────────────────────────────────────────────────
    plt.rcParams.update({
        "font.family":      "DejaVu Sans",
        "font.size":        10,
        "axes.titlesize":   11,
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

    UNIT_SUBLABELS = {
        "Biostimulant":     "$/kg product",
        "AD-biomethane":    "$/kg CH4",
        "AD-VFA":           "$/kg total VFA",
        "AD-fermentation":  "$/kg crude oil",
    }
    tick_labels = [
        f"{r['label']}\n({UNIT_SUBLABELS[r['label']]})" for r in results
    ]
    msp_vals = [r["msp_usd_per_kg"] for r in results]

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars = ax.bar(tick_labels, msp_vals, edgecolor="black", linewidth=0.8, zorder=3)

    for bar, val in zip(bars, msp_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(msp_vals) * 0.015,
            f"${val:,.2f}",
            ha="center", va="bottom", fontsize=9,
        )

    ax.set_ylabel("Minimum selling price ($/kg product)")
    ax.set_title(
        "MSP by SaBRe pathway (near-zero feed price)\n"
        "Each bar is priced on its own pathway's value-carrying product",
        fontsize=10,
    )
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.set_ylim(0, max(msp_vals) * 1.15)

    fig.tight_layout()
    fig.savefig(OUT / "fig_msp_comparison_four_systems.png", bbox_inches="tight")

    print("\nSaved:")
    print(f"  {OUT / 'fig_msp_comparison_four_systems.png'}")

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
    create_biostimulant_system,
    create_ad_biomethane_system,
    create_ad_vfa_system,
    create_ad_fermentation_system,
)
from biorefineries.sabre._tea import solve_product_msp, CH4_MMBTU_PER_KG

OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# =================================================================
# Shared baseline economics
# Values match the defaults already used across
# legacy_analyses/ad_feed_tea.py and legacy_analyses/vfa_fermentation_tea.py
# so all four systems are compared on a common basis.
# =================================================================

# Not used now, kept as a reference
# def _patch_ev607():
#     """
#     Replace Ev607's cost + utility with a low-duty placeholder when V < 0.02.
#     BioSTEAM's MultiEffectEvaporator cost correlation produces nonsensical
#     vessel geometry at near-zero evaporation duty (see
#     legacy_analyses/vfa_fermentation_tea.py::_patch_ev607).
#     """
#     try:
#         ev607 = bst.main_flowsheet.unit["Ev607"]
#     except Exception:
#         return
#     v = getattr(ev607, "V", None)
#     if v is None or v >= 0.02:
#         return
#     feed = ev607.ins[0]
#     feed_m3h = max(feed.F_mass / 1000.0, 1.0)
#     placeholder_usd = 50000.0 * (feed_m3h ** 0.6)
#     for k in list(ev607.baseline_purchase_costs.keys()):
#         ev607.baseline_purchase_costs[k] = 0.0
#     ev607.baseline_purchase_costs["Evaporator (low-duty placeholder)"] = placeholder_usd
#     ev607.heat_utilities.clear()


# =================================================================
# System builders -- each returns a common result dict
# =================================================================

def price_biostimulant_system() -> dict:
    bst.main_flowsheet.clear()

    sys = create_biostimulant_system()
    sys.simulate()

    product = sys.flowsheet.stream.biostimulant_product
    msp = solve_product_msp(tea=sys.TEA, product_stream=product)

    return {
        "label": "Biostimulant",
        "product_desc": "biostimulant liquid product",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
    }


def price_ad_biomethane_system() -> dict:
    bst.main_flowsheet.clear()

    sys = create_ad_biomethane_system(feedstock="pelagic")
    sys.simulate()

    product = sys.flowsheet.stream.biomethane
    msp = solve_product_msp(
        tea=sys.TEA, product_stream=product,
        energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
    )

    return {
        "label": "AD-biomethane",
        "product_desc": "biomethane (whole-stream basis)",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "msp_usd_per_mmbtu": msp["usd_per_mmbtu"],
        "annual_product_mmbtu": msp["annual_product_mmbtu"],
    }


def price_ad_vfa_system() -> dict:
    bst.main_flowsheet.clear()

    sys = create_ad_vfa_system(feedstock="pelagic")
    sys.simulate()
    
    product = sys.flowsheet.stream.pure_vfa
    msp = solve_product_msp(tea=sys.TEA, product_stream=product)

    return {
        "label": "AD-VFA",
        "product_desc": "VFA broth (total-VFA basis)",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
    }


def price_ad_fermentation_system() -> dict:
    bst.main_flowsheet.clear()

    sys = create_ad_fermentation_system(feedstock="pelagic")
    sys.simulate()
    # _patch_ev607() not used now, kept as a reference

    product = sys.flowsheet.stream.microbial_oil
    msp = solve_product_msp(tea=sys.TEA, product_stream=product)

    return {
        "label": "AD-fermentation",
        "product_desc": "crude microbial oil",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
    }


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

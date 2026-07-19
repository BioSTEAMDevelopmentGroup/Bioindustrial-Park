# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
plot_feed_price_sensitivity.py
------------------------------
Plot biomethane MSP as a function of feedstock price for each AD pretreatment case.

Feed price scenarios ($/kg wet):
    -0.02  gate fee / tipping fee (operator is paid to take feedstock)
     0.00  free feedstock
     0.02  baseline purchase price
     0.05  high harvest cost

Run from sabre_project root:
    python scripts/plot_feed_price_sensitivity.py
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

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.systems import create_ad_biomethane_system
from biorefineries.sabre._tea import create_tea, solve_biomethane_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PRETREATMENT_CASES = [
    "press_mill_only",
    "enzymatic",
    "peroxide",
    "combined_PE",
    "combined_PTE",
]

CASE_LABELS = {
    "press_mill_only": "Press + mill only",
    "enzymatic": "Enzymatic",
    "peroxide": "Peroxide",
    "combined_PE": "Combined PE",
    "combined_PTE": "Combined PTE",
}

FEED_PRICES_USD_PER_KG_WET = [-0.02, 0.00, 0.02, 0.05]

FEED_PRICE_LABELS = {
    -0.02: "−$0.02/kg (gate fee)",
     0.00: "$0.00/kg (free)",
     0.02: "$0.02/kg (baseline)",
     0.05: "$0.05/kg (premium)",
}

BIOSTIMULANT_PRICE_USD_PER_KG  = 0.00
LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG  = -0.002
SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG  = -0.02


def apply_stream_economics(sys, biostimulant_price: float = BIOSTIMULANT_PRICE_USD_PER_KG):
    for sid in ("pressate_permeate",):
        try:
            sys.flowsheet.stream[sid].price = 0.0
        except Exception:
            pass
    for sid in ("biostimulant_product",):
        try:
            sys.flowsheet.stream[sid].price = biostimulant_price
        except Exception:
            pass
    for sid in ("soil_amendment",):
        try:
            sys.flowsheet.stream[sid].price = SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG
        except Exception:
            pass
    for sid in ("liquid_digestate",):
        try:
            sys.flowsheet.stream[sid].price = LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG
        except Exception:
            pass


def build_case(pretreatment_case: str, feed_price: float):
    bst.main_flowsheet.clear()
    set_thermo()

    sys = create_ad_biomethane_system(
        feedstock_type="pelagic",
        pretreatment_case=pretreatment_case,
    )
    sys.feeds[0].price = feed_price
    sys.simulate()

    apply_stream_economics(sys)

    tea = create_tea(sys)
    msp = solve_biomethane_msp(tea, sys.flowsheet.stream.biomethane)

    # Annual feedstock cost for reporting
    feed_kgph = float(sys.feeds[0].F_mass)
    annual_feed_cost = feed_price * feed_kgph * 330 * 24  # USD/yr

    return float(msp["usd_per_mmbtu"]), annual_feed_cost / 1e6


# ── Run simulations ──────────────────────────────────────────────────────────
if __name__ == "__main__":

    results = {case: [] for case in PRETREATMENT_CASES}
    feed_costs = {case: [] for case in PRETREATMENT_CASES}

    for case in PRETREATMENT_CASES:
        for fp in FEED_PRICES_USD_PER_KG_WET:
            msp_val, fc_musd = build_case(case, fp)
            results[case].append(msp_val)
            feed_costs[case].append(fc_musd)
            print(f"{case:20s}  feed=${fp:+.2f}/kg  MSP=${msp_val:.2f}/MMBtu  feed cost=${fc_musd:.1f}M/yr")

    # ── Plot ─────────────────────────────────────────────────────────────────────

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

    markers = {
        "press_mill_only": "o",
        "enzymatic":       "s",
        "peroxide":        "^",
        "combined_PE":     "D",
        "combined_PTE":    "v",
    }

    fig, ax = plt.subplots(figsize=(7.0, 4.8))

    for case in PRETREATMENT_CASES:
        ax.plot(
            FEED_PRICES_USD_PER_KG_WET,
            results[case],
            marker=markers[case],
            linewidth=1.8,
            markersize=6,
            markeredgecolor="black",
            label=CASE_LABELS[case],
            zorder=3,
        )

    # Market price reference lines
    ax.axhline(3.0,  color="black", linewidth=1.3, linestyle="--",  zorder=2, label="$3/MMBtu")
    ax.axhline(10.0, color="black", linewidth=1.3, linestyle=":",   zorder=2, label="$10/MMBtu")
    ax.axhline(14.0, color="black", linewidth=1.3, linestyle="-.",  zorder=2, label="$14/MMBtu")

    ax.set_xlabel("Feedstock price ($/kg wet)")
    ax.set_ylabel("Biomethane MSP ($/MMBtu)")
    ax.set_title("Effect of feedstock price on biomethane MSP")

    ax.axvline(
        0.02,
        color="black",
        linewidth=1.0,
        linestyle="--",
        zorder=2,
        label="Baseline ($0.02/kg)",
    )

    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("$%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")

    fig.tight_layout()
    fig.savefig(OUT / "fig_feed_price_sensitivity.png", bbox_inches="tight")

    # ── Console summary ───────────────────────────────────────────────────────────

    print("\nMSP results ($/MMBtu):")
    for case in PRETREATMENT_CASES:
        print(f"\n{CASE_LABELS[case]}")
        for fp, msp, fc in zip(FEED_PRICES_USD_PER_KG_WET, results[case], feed_costs[case]):
            sign = "income" if fc < 0 else "cost"
            print(f"  feed = {fp:+.2f} $/kg  ->  MSP = ${msp:.2f}/MMBtu  |  annual feed {sign} = ${abs(fc):.1f}M/yr")

    print("\nSaved:")
    print(OUT / "fig_feed_price_sensitivity.png")
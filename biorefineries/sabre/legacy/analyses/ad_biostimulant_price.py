# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
plot_biostimulant_sensitivity.py --> feed price and biostimulant price sensitivity for AD pretreatment case
--------------------------------
Plot biomethane MSP as a function of biostimulant selling price for the
methanogenic AD pretreatment cases

Run from sabre_project root:
    python scripts/plot_biostimulant_sensitivity.py
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

# Same sensitivity points used in ad_tea_final.py
BIOSTIMULANT_PRICES = [0.00, 0.50, 1.00, 1.50, 2.00, 2.50, 3.00]

FEED_PRICE_USD_PER_KG_WET = 0.02
LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG = -0.002
SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG = -0.02


def apply_stream_economics(sys, biostimulant_price: float):
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


def build_case(pretreatment_case: str, biostimulant_price: float):
    bst.main_flowsheet.clear()
    set_thermo()

    sys = create_ad_biomethane_system(
        feedstock="pelagic",
        pretreatment_case=pretreatment_case,
    )
    sys.feeds[0].price = FEED_PRICE_USD_PER_KG_WET
    sys.simulate()

    apply_stream_economics(sys, biostimulant_price=biostimulant_price)

    tea = create_tea(sys)
    msp = solve_biomethane_msp(tea, sys.flowsheet.stream.biomethane)

    annual_revenue = 0.0
    try:
        conc = sys.flowsheet.stream["biostimulant_product"]
        annual_revenue = biostimulant_price * conc.F_mass * 330 * 24
    except Exception:
        pass

    return float(msp["usd_per_mmbtu"]), annual_revenue / 1e6


results = {case: [] for case in PRETREATMENT_CASES}
if __name__ == "__main__":
    revenues = {case: [] for case in PRETREATMENT_CASES}

    for case in PRETREATMENT_CASES:
        for bs_price in BIOSTIMULANT_PRICES:
            msp_mmbtu, annual_rev_musd = build_case(case, bs_price)
            results[case].append(msp_mmbtu)
            revenues[case].append(annual_rev_musd)

    # plotting styles
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

    fig, ax = plt.subplots(figsize=(7.0, 4.8))

    markers = {
        "press_mill_only": "o",
        "enzymatic": "s",
        "peroxide": "^",
        "combined_PE": "D",
        "combined_PTE": "v",
    }

    for case in PRETREATMENT_CASES:
        ax.plot(
            BIOSTIMULANT_PRICES,
            results[case],
            marker=markers[case],
            linewidth=1.8,
            markersize=6,
            markeredgecolor="black",
            label=CASE_LABELS[case],
            zorder=3,
        )

    ax.axhline(3.0, color="black", linewidth=1.3, linestyle="--", zorder=2, label="$3/MMBtu")
    ax.axhline(10.0, color="black", linewidth=1.3, linestyle=":", zorder=2, label="$10/MMBtu")
    ax.axhline(14.0, color="black", linewidth=1.3, linestyle="-.", zorder=2, label="$14/MMBtu")

    ax.set_xlabel("Biostimulant price ($/kg)")
    ax.set_ylabel("Biomethane MSP ($/MMBtu)")
    ax.set_title("Effect of biostimulant price on biomethane MSP")
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="both", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
    ax.axvline(
        0.00,
        color="black",
        linewidth=1.0,
        linestyle="--",
        zorder=2,
        label="Baseline ($0.00/kg)",
    )

    fig.tight_layout()
    fig.savefig(OUT / "fig_biostimulant_sensitivity.png", bbox_inches="tight")

    print("\nMSP results ($/MMBtu):")
    for case in PRETREATMENT_CASES:
        print(f"\n{CASE_LABELS[case]}")
        for bs_price, msp, rev in zip(BIOSTIMULANT_PRICES, results[case], revenues[case]):
            print(f"  biostimulant = ${bs_price:.2f}/kg -> MSP = ${msp:.3f}/MMBtu | annual revenue = ${rev:.2f}M/yr")

    print("\nSaved:")
    print(OUT / "fig_biostimulant_sensitivity.png")
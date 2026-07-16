# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
plot_two_pretreatment_figures_from_tea_with_errorbars.py
--------------------------------------------------------
Recalculate methane production and MSP directly from the current SABRE model and TEA,
including error bars based on the low/high methane-yield estimates in assumptions.yaml.

Outputs:
1) fig_methane_production_by_pretreatment_errorbars.png
2) fig_msp_by_pretreatment_errorbars.png

BASECASE: biostimulant: $0/kg (no revenue), feedstock: $0.02/kg wet
"""

from __future__ import annotations

import sys
from pathlib import Path

import biosteam as bst
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.systems import create_ad_biogas_system
from biorefineries.sabre._tea import make_baseline_tea, solve_biomethane_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Reads the methane yield effects for each pretreatment case from assumptions.yaml
PRETREATMENT_CASES = [
    "press_mill_only",
    "enzymatic",
    "peroxide",
    "combined_PE",
    "combined_PTE",
]

CASE_LABELS = {
    "press_mill_only": "Press+mill\nonly",
    "enzymatic": "Enzymatic",
    "peroxide": "Peroxide",
    "combined_PE": "Combined\nPE",
    "combined_PTE": "Combined\nPTE",
}

LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG = -0.002
SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG = -0.02
BIOSTIMULANT_PRICE_BASE_USD_PER_KG = 0.00
FEED_PRICE_USD_PER_KG_WET = 0.02


def apply_stream_economics(
    sys,
    biostimulant_price: float = BIOSTIMULANT_PRICE_BASE_USD_PER_KG,
):
    # Pressate permeate
    for sid in ("pressate_permeate",):
        try:
            sys.flowsheet.stream[sid].price = 0.0
        except Exception:
            pass

    # Biostimulant concentrate
    for sid in ("biostimulant_membrane_concentrate",):
        try:
            sys.flowsheet.stream[sid].price = biostimulant_price
        except Exception:
            pass

    # Solid digestate disposal
    for sid in ("soil_amendment",):
        try:
            sys.flowsheet.stream[sid].price = SOLIDS_DIGESTATE_DISPOSAL_USD_PER_KG
        except Exception:
            pass

    # Liquid digestate disposal
    for sid in ("liquid_digestate",):
        try:
            sys.flowsheet.stream[sid].price = LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG
        except Exception:
            pass


def build_case(
    pretreatment_case: str,
    feed_price: float = FEED_PRICE_USD_PER_KG_WET,
    yield_override: float | None = None,
):
    bst.main_flowsheet.clear()
    set_thermo()

    sys = create_ad_biogas_system(
        quality="pelagic_high_quality",
        pretreatment_case=pretreatment_case,
        ch4_override=yield_override,
    )

    if yield_override is not None:
        sys.flowsheet.unit.AD.ch4_kg_per_kg_vs_fed = float(yield_override)

    sys.feeds[0].price = feed_price
    sys.simulate()

    apply_stream_economics(sys, biostimulant_price=BIOSTIMULANT_PRICE_BASE_USD_PER_KG)

    tea = make_baseline_tea(sys)
    msp = solve_biomethane_msp(tea, sys.flowsheet.stream.biomethane)
    ch4_kgph = float(sys.flowsheet.stream.biomethane.imass["Methane"])

    return sys, tea, msp, ch4_kgph


A = load_assumptions()
CASE_DEFS = A["ad_pretreatment_cases"]


def get_case_yields(case: str):
    ad_effects = CASE_DEFS[case]["ad_effects"]
    central = float(ad_effects["ch4_kg_per_kg_vs_fed"])
    low = float(ad_effects.get("ch4_kg_per_kg_vs_fed_low", central))
    high = float(ad_effects.get("ch4_kg_per_kg_vs_fed_high", central))
    return central, low, high


# Plotting style
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

methane_kgph = []
if __name__ == "__main__":
    methane_err_low = []
    methane_err_high = []

    msp_mmbtu = []
    msp_err_low = []
    msp_err_high = []

    labels = [CASE_LABELS[c] for c in PRETREATMENT_CASES]

    for case in PRETREATMENT_CASES:
        central_y, low_y, high_y = get_case_yields(case)

        sys_c, tea_c, msp_c, ch4_c = build_case(case, yield_override=central_y)
        sys_l, tea_l, msp_l, ch4_l = build_case(case, yield_override=low_y)
        sys_h, tea_h, msp_h, ch4_h = build_case(case, yield_override=high_y)

        methane_kgph.append(ch4_c)
        methane_err_low.append(max(ch4_c - ch4_l, 0.0))
        methane_err_high.append(max(ch4_h - ch4_c, 0.0))

        msp_c_val = float(msp_c["usd_per_mmbtu"])
        msp_l_val = float(msp_l["usd_per_mmbtu"])
        msp_h_val = float(msp_h["usd_per_mmbtu"])

        msp_mmbtu.append(msp_c_val)
        msp_err_low.append(max(msp_c_val - msp_h_val, 0.0))
        msp_err_high.append(max(msp_l_val - msp_c_val, 0.0))

    x = np.arange(len(labels))
    best_idx = int(np.argmin(msp_mmbtu))


    def bar_colors(n, best_idx, default, best):
        return [best if i == best_idx else default for i in range(n)]


    # ============================================================
    # Figure 1 — Methane production by pretreatment case
    # ============================================================
    fig1, ax1 = plt.subplots(figsize=(6.5, 4))

    bars1 = ax1.bar(
        x,
        methane_kgph,
        color=bar_colors(len(labels), best_idx, "#74c476", "#238b45"),
        edgecolor="black",
        linewidth=0.8,
        width=0.55,
        yerr=np.array([methane_err_low, methane_err_high]),
        capsize=6,
        error_kw={
            "elinewidth": 1.0,
            "ecolor": "#222222",
            "capthick": 1.0,
            "zorder": 5,
        },
        zorder=3,
    )

    for bar, val, ehi in zip(bars1, methane_kgph, methane_err_high):
        ax1.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ehi + max(methane_kgph) * 0.008,
            f"{val:,.0f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.set_ylabel("CH₄ production (kg/hr)")
    ax1.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax1.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))

    # force room above bars so caps are not cramped
    ax1.set_ylim(0, max(np.array(methane_kgph) + np.array(methane_err_high)) * 1.10)

    fig1.tight_layout()
    fig1.savefig(OUT / "fig_methane_production_by_pretreatment_errorbars.png", bbox_inches="tight")


    # ============================================================
    # Figure 2 — MSP by pretreatment case
    # ============================================================
    fig2, ax2 = plt.subplots(figsize=(6.5, 4))

    bars2 = ax2.bar(
        x,
        msp_mmbtu,
        color=bar_colors(len(labels), best_idx, "#6baed6", "#2171b5"),
        edgecolor="black",
        linewidth=0.8,
        width=0.55,
        yerr=np.array([msp_err_low, msp_err_high]),
        capsize=6,
        error_kw={
            "elinewidth": 1.0,
            "ecolor": "#222222",
            "capthick": 1.0,
            "zorder": 5,
        },
        zorder=3,
    )

    for bar, val, ehi in zip(bars2, msp_mmbtu, msp_err_high):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + ehi + max(msp_mmbtu) * 0.015,
            f"${val:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax2.axhline(
        3.0,
        color="black",
        linewidth=1.3,
        linestyle="--",
        zorder=2,
        label="Market price ~$3/MMBtu",
    )

    ax2.axhline(
        10.0,
        color="black",
        linewidth=1.3,
        linestyle=":",
        zorder=2,
        label="Market price ~$10/MMBtu",
    )

    ax2.axhline(
        14.0,
        color="black",
        linewidth=1.3,
        linestyle="-.",
        zorder=2,
        label="Market price ~$14/MMBtu",
    )

    ax2.set_xticks(x)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Biomethane MSP ($/MMBtu)")
    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax2.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax2.legend(fontsize=8, frameon=False)

    fig2.tight_layout()
    fig2.savefig(OUT / "fig_msp_by_pretreatment_errorbars.png", bbox_inches="tight")


    print("Methane production (kg/hr) with low/high bounds:")
    for label, val, elo, ehi in zip(labels, methane_kgph, methane_err_low, methane_err_high):
        print(f"  {label:<18} {val:,.2f} (-{elo:,.2f}, +{ehi:,.2f})")

    print("\nMSP ($/MMBtu) with low/high bounds:")
    for label, val, elo, ehi in zip(labels, msp_mmbtu, msp_err_low, msp_err_high):
        print(f"  {label:<18} {val:.3f} (-{elo:.3f}, +{ehi:.3f})")

    print("\nSaved:")
    print(OUT / "fig_methane_production_by_pretreatment_errorbars.png")
    print(OUT / "fig_msp_by_pretreatment_errorbars.png")
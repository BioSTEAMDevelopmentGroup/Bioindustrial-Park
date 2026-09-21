# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_hydrolysis_3hp.py
---------------------
Sensitivity of the 3-HP (dried Ca(3HP)2) MSP to the EnzymaticPress
hydrolysis conversion (data/preprocessing.yaml `enzymatic_press.
hydrolysis_conversion`, currently 0.50: the fraction of Glucan and of
Alginate hydrolyzed to Glucose / AlginateMonomer), swept from 0.20 to 0.90
at the 15,000 wet ton/day baseline and the current 20 wt% fermentation-feed
solids target.

Setting `PR.hydrolysis_conversion` after the build is NOT enough:
EnzymaticPress builds its `hydrolysis_rxns` ParallelReaction in __init__
with the conversion baked in, so the reactions' own `.X` must be set too
(done here, and checked below against the measured Glucan conversion).

The enzyme dose is set by EnzymaticPress._run from the Glucan + Alginate
feed mass alone (20 mg/g + 10 % excess), independent of the conversion, and
the unit has no saccharification-reactor cost of its own. So a higher
conversion here costs no extra enzyme, time, or equipment; the sweep shows
the value of the conversion, not the cost of achieving it.

The fermentation feed is held at the same total-solids target, so sugar
concentration and titer move with the conversion (the non-sugar solubles
are unchanged, the sugars are not).

Outputs (relative to this script's directory):
    ../results/tables/msp_hydrolysis_3hp.csv
    ../results/figures/fig_msp_hydrolysis_3hp.png

Run from the repo root:
    python biorefineries/sabre/analyses/msp_hydrolysis_3hp.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import biosteam as bst

from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.systems import _3hp_system
from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.analyses.msp_scale_scenarios_3hp import (
    SURFACE, INK, INK_2, GRID, AXIS, SERIES_WITH, TIGHT_TOLERANCE,
)
from biorefineries.sabre.analyses.msp_fermentation_solids_3hp import SERIES_TITER, SUGAR_IDS
from biorefineries.sabre.analyses.msp_yield_productivity_3hp import _style_axes

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

CURRENT_CONVERSION = load_assumptions("preprocessing.yaml")["enzymatic_press"]["hydrolysis_conversion"]
CONVERSIONS = np.round(np.arange(0.20, 0.9001, 0.05), 4)


def run_case(hydrolysis_conversion: float | None = None) -> dict:
    """Build/simulate the baseline system at one hydrolysis conversion (None = as built)."""
    bst.main_flowsheet.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sys = _3hp_system.create_3hp_system()
        sys.set_tolerance(subsystems=True, **TIGHT_TOLERANCE)
        PR = sys.flowsheet.unit.PR
        if hydrolysis_conversion is not None:
            PR.hydrolysis_conversion = hydrolysis_conversion
            PR.hydrolysis_rxns.X = np.full(len(PR.hydrolysis_rxns.X), hydrolysis_conversion)
        sys.simulate()
        msp = solve_product_msp(sys.TEA, sys.flowsheet.stream.ca3hp2_product)

    # Check the override actually took effect: measured Glucan conversion.
    chems = sys.flowsheet.unit.PR.chemicals
    feed, pressate = PR.ins[0], PR.outs[1]
    measured = float(
        pressate.imass["Glucose"] * chems.Glucan.MW / chems.Glucose.MW / feed.imass["Glucan"]
    )
    if hydrolysis_conversion is not None and abs(measured - hydrolysis_conversion) > 1e-6:
        raise RuntimeError(
            f"hydrolysis conversion override not applied: wanted {hydrolysis_conversion}, "
            f"measured {measured:.6f}"
        )

    ff = sys.flowsheet.stream.fermentation_feed
    F401 = sys.flowsheet.unit.F401
    F_vol = float(ff.F_vol)
    sugar_g_per_L = {i: float(ff.imass[i]) / F_vol for i in SUGAR_IDS}
    return {
        "hydrolysis_conversion": PR.hydrolysis_conversion,
        "measured_glucan_conversion": measured,
        "sugar_g_per_L": sum(sugar_g_per_L.values()),
        **{f"{i}_g_per_L": v for i, v in sugar_g_per_L.items()},
        "titer_g_per_L": float(F401.design_results["Titer (g HP-equivalent/L)"]),
        "fermentation_batch_h": float(F401.tau),
        "F401_installed_usd": float(F401.installed_cost),
        "FCI_usd": sys.TEA.FCI,
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
        "msp_usd_per_kg": msp["usd_per_kg"],
    }


def plot_sweep(current: dict, df: pd.DataFrame, out_path: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "figure.dpi": 100, "savefig.dpi": 200,
    })
    fig = plt.figure(figsize=(9.6, 5.6), facecolor=SURFACE)
    left, right = 0.08, 0.66
    ax_msp = fig.add_axes([left, 0.55, right - left, 0.29], facecolor=SURFACE)
    ax_sug = fig.add_axes([left, 0.14, right - left, 0.29], facecolor=SURFACE, sharex=ax_msp)

    x = df["hydrolysis_conversion"] * 100.0
    x_cur = current["hydrolysis_conversion"] * 100.0

    msp_top = float(np.ceil(df["msp_usd_per_kg"].max() * 1.12 * 2) / 2)
    sug_top = float(np.ceil(df["sugar_g_per_L"].max() * 1.12 / 100) * 100)
    panels = (
        (ax_msp, r"MSP (USD/kg Ca(3HP)$_2$)", (0, msp_top), list(np.arange(0, msp_top + 0.01, 0.5)), "{:.2f}",
         (("msp_usd_per_kg", SERIES_WITH, None),)),
        (ax_sug, "Concentration (g/L)", (0, sug_top), list(np.arange(0, sug_top + 1, 100)), "{:.0f}",
         (("sugar_g_per_L", SERIES_WITH, "Sugars in feed"),
          ("titer_g_per_L", SERIES_TITER, "Titer (3-HP equivalent)"))),
    )
    for ax, ylabel, ylim, yticks, fmt, series in panels:
        ax.set_xlim(17, 93)
        _style_axes(ax, ylim, yticks)
        ax.set_ylabel(ylabel, color=INK_2, fontsize=10)
        ax.axvline(x_cur, color=AXIS, linewidth=1.0, zorder=1)
        y_span = ylim[1] - ylim[0]
        for col, color, legend_label in series:
            ax.plot(x, df[col], color=color, linewidth=2.0, solid_capstyle="round",
                    solid_joinstyle="round", zorder=3, label=legend_label)
            ax.plot([x_cur], [current[col]], marker="o", markersize=9, color=color,
                    markeredgecolor=SURFACE, markeredgewidth=2, linestyle="none", zorder=4)
            ends = [(x.iloc[-1], df[col].iloc[-1], "right")]
            if len(series) == 1:
                ends.insert(0, (x.iloc[0], df[col].iloc[0], "left"))
            for xi, yi, ha in ends:
                ax.text(xi, yi + 0.06 * y_span, fmt.format(yi), ha=ha, va="bottom",
                        color=INK, fontsize=10, zorder=5)
        if len(series) > 1:
            leg = ax.legend(loc="upper left", frameon=False, fontsize=9.5, handlelength=1.6,
                            handletextpad=0.6, borderaxespad=0.2)
            for t in leg.get_texts():
                t.set_color(INK_2)

    plt.setp(ax_msp.get_xticklabels(), visible=False)
    ticks = list(range(20, 91, 10))
    ax_sug.set_xticks(ticks)
    ax_sug.set_xticklabels([f"{v}%" for v in ticks])
    ax_sug.tick_params(axis="x", pad=6)
    ax_sug.set_xlabel("EnzymaticPress hydrolysis conversion of Glucan and Alginate (%)",
                      color=INK_2, fontsize=10, labelpad=8)

    fig.text(left, 0.945, "Ca(3HP)$_2$ MSP vs. enzymatic hydrolysis conversion",
             color=INK, fontsize=12.5, fontweight="semibold", ha="left", va="center")
    fig.text(left, 0.895, "15,000 t/day baseline, 20 wt% fermentation-feed solids",
             color=INK_2, fontsize=10, ha="left", va="center")

    nx = 0.71
    fig.text(nx, 0.78, "Current assumption", color=INK, fontsize=10.5,
             fontweight="semibold", ha="left", va="center")
    for txt, y in (
        (f"{x_cur:.0f}% hydrolysis conversion", 0.735),
        (f"{current['sugar_g_per_L']:.1f} g/L sugars", 0.695),
        (f"{current['titer_g_per_L']:.1f} g/L titer", 0.655),
        (f"MSP ${current['msp_usd_per_kg']:.2f}/kg", 0.615),
    ):
        fig.text(nx, y, txt, color=INK_2, fontsize=10, ha="left", va="center")
    fig.text(nx, 0.54,
             "Enzyme dose is fixed and the\nunit has no reactor cost, so\nhigher conversion is free here.\n"
             "No titer limit is modeled.",
             color=INK_2, fontsize=8.5, ha="left", va="top", linespacing=1.35)

    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    current = {"case": "current", **run_case()}
    df = pd.DataFrame([{"case": "sweep", **run_case(float(c))} for c in CONVERSIONS])
    csv_path = TABLES_OUT / "msp_hydrolysis_3hp.csv"
    fig_path = FIGURES_OUT / "fig_msp_hydrolysis_3hp.png"

    pd.concat([pd.DataFrame([current]), df], ignore_index=True).to_csv(csv_path, index=False)
    plot_sweep(current, df, fig_path)

    print(f"Saved table:  {csv_path}")
    print(f"Saved figure: {fig_path}\n")
    cols = ["hydrolysis_conversion", "sugar_g_per_L", "titer_g_per_L", "fermentation_batch_h",
            "annual_product_t", "msp_usd_per_kg"]
    with pd.option_context("display.width", 200, "display.max_columns", None):
        fmt = lambda v: f"{v:,.4g}"
        print("current:\n", pd.DataFrame([current])[cols].to_string(index=False, float_format=fmt), "\n")
        print(df[cols].to_string(index=False, float_format=fmt))

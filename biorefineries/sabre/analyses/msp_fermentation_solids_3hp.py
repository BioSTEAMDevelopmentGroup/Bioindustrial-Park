# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_fermentation_solids_3hp.py
------------------------------
Sensitivity of the 3-HP (dried Ca(3HP)2) MSP to the total-solids target
the pressate is concentrated to ahead of fermentation, at the 15,000
wet ton/day baseline (data/3hp.yaml `fermentation_feed_evaporator.
target_solids_wt_frac`, currently 0.20), swept from 5 to 25 wt%.

The target is set on `EV.target_solids_wt_frac` after the build. It drives
both the evaporator itself and the PFS split that decides how much pressate
goes through the membrane concentrator (systems._3hp_system's
`adjust_pressate_split` reads it at run time), so the membrane/evaporator
duty split moves with it.

Also reports the concentration of the fermentable sugars in the
fermentation feed (Glucose + Mannitol + AlginateMonomer, the
`pressate_concentrator.retained_solute_IDs`), in g/L, and the fermentation
titer (HPFermentation's simulated "Titer (g HP-equivalent/L)", an output of
the feed sugar level and the fixed conversions in data/3hp.yaml, not an
input) on the same g/L axis.

Outputs (relative to this script's directory):
    ../results/tables/msp_fermentation_solids_3hp.csv
    ../results/figures/fig_msp_fermentation_solids_3hp.png

Run from the repo root:
    python biorefineries/sabre/analyses/msp_fermentation_solids_3hp.py
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

SERIES_TITER = "#eb6834"  # categorical slot 2 (orange)

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

CURRENT_TARGET = load_assumptions("3hp.yaml")["fermentation_feed_evaporator"]["target_solids_wt_frac"]
SUGAR_IDS = tuple(load_assumptions("3hp.yaml")["pressate_concentrator"]["retained_solute_IDs"])
TARGETS = np.round(np.arange(0.05, 0.2501, 0.01), 4)


def run_target(target_solids_wt_frac: float) -> dict:
    """Build/simulate the baseline system at one fermentation-feed solids target."""
    bst.main_flowsheet.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sys = _3hp_system.create_3hp_system()
        sys.set_tolerance(subsystems=True, **TIGHT_TOLERANCE)
        sys.flowsheet.unit.EV.target_solids_wt_frac = target_solids_wt_frac
        sys.simulate()
        msp = solve_product_msp(sys.TEA, sys.flowsheet.stream.ca3hp2_product)

    feed = sys.flowsheet.stream.fermentation_feed
    F_vol = float(feed.F_vol)
    sugar_g_per_L = {i: float(feed.imass[i]) / F_vol for i in SUGAR_IDS}  # kg/hr / (m3/hr) = g/L
    return {
        "target_solids_wt_frac": target_solids_wt_frac,
        "feed_solids_wt_frac": 1.0 - float(feed.imass["Water"]) / float(feed.F_mass),
        "sugar_g_per_L": sum(sugar_g_per_L.values()),
        **{f"{i}_g_per_L": v for i, v in sugar_g_per_L.items()},
        "frac_pressate_to_membrane": (
            float(sys.flowsheet.unit.PFS.outs[0].F_mass)
            / float(sys.flowsheet.stream.pressate.F_mass)
        ),
        "titer_g_per_L": float(sys.flowsheet.unit.F401.design_results["Titer (g HP-equivalent/L)"]),
        "fermentation_batch_h": float(sys.flowsheet.unit.F401.tau),
        "FCI_usd": sys.TEA.FCI,
        "msp_usd_per_kg": msp["usd_per_kg"],
    }


def plot_sweep(df: pd.DataFrame, out_path: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "figure.dpi": 100, "savefig.dpi": 200,
    })
    fig = plt.figure(figsize=(9.6, 5.6), facecolor=SURFACE)
    left, right = 0.08, 0.66
    ax_msp = fig.add_axes([left, 0.55, right - left, 0.29], facecolor=SURFACE)
    ax_sug = fig.add_axes([left, 0.14, right - left, 0.29], facecolor=SURFACE, sharex=ax_msp)

    x = df["target_solids_wt_frac"] * 100.0
    cur = df.loc[np.isclose(df["target_solids_wt_frac"], CURRENT_TARGET)].iloc[0]
    x_cur = CURRENT_TARGET * 100.0

    # (axes, ylabel, ylim, yticks, label format, series); each series is
    # (column, color, legend label, end-label sides).
    panels = (
        (ax_msp, r"MSP (USD/kg Ca(3HP)$_2$)", (0, 2.5), [0, 0.5, 1.0, 1.5, 2.0, 2.5], "{:.2f}",
         (("msp_usd_per_kg", SERIES_WITH, None, "both"),)),
        # Sugars and titer share one g/L axis. Only the right ends are
        # labeled: at the left end the two lines are closer than a label is tall.
        (ax_sug, "Concentration (g/L)", (0, 300), [0, 100, 200, 300], "{:.0f}",
         (("sugar_g_per_L", SERIES_WITH, "Sugars in feed", "right"),
          ("titer_g_per_L", SERIES_TITER, "Titer (3-HP equivalent)", "right"))),
    )
    for ax, ylabel, ylim, yticks, fmt, series in panels:
        ax.set_xlim(4, 26)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8, linestyle="-")
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(AXIS)
        ax.tick_params(axis="both", length=0, colors=INK_2, labelsize=9.5)
        ax.set_ylabel(ylabel, color=INK_2, fontsize=10)

        # Current-assumption guide: solid hairline, recessive.
        ax.axvline(x_cur, color=AXIS, linewidth=1.0, zorder=1)
        y_span = ylim[1] - ylim[0]
        for col, color, legend_label, ends in series:
            ax.plot(x, df[col], color=color, linewidth=2.0, solid_capstyle="round",
                    solid_joinstyle="round", zorder=3, label=legend_label)
            # Current point: >= 8 px marker with a 2 px surface ring.
            ax.plot([x_cur], [cur[col]], marker="o", markersize=9, color=color,
                    markeredgecolor=SURFACE, markeredgewidth=2, linestyle="none", zorder=4)

            # Direct labels at the ends of the sweep, above the line.
            end_points = [(x.iloc[-1], df[col].iloc[-1], "right")]
            if ends == "both":
                end_points.insert(0, (x.iloc[0], df[col].iloc[0], "left"))
            for xi, yi, ha in end_points:
                ax.text(xi, yi + 0.06 * y_span, fmt.format(yi), ha=ha, va="bottom",
                        color=INK, fontsize=10, zorder=5)

        if len(series) > 1:
            leg = ax.legend(loc="upper left", frameon=False, fontsize=9.5, handlelength=1.6,
                            handletextpad=0.6, borderaxespad=0.2)
            for t in leg.get_texts():
                t.set_color(INK_2)

    plt.setp(ax_msp.get_xticklabels(), visible=False)
    ax_sug.set_xticks(range(5, 26, 5))
    ax_sug.set_xticklabels([f"{v}%" for v in range(5, 26, 5)])
    ax_sug.tick_params(axis="x", pad=6)
    ax_sug.set_xlabel("Target total solids of fermentation feed after concentration (wt%)",
                      color=INK_2, fontsize=10, labelpad=8)

    fig.text(left, 0.945, "Ca(3HP)$_2$ MSP vs. fermentation-feed solids target",
             color=INK, fontsize=12.5, fontweight="semibold", ha="left", va="center")
    fig.text(left, 0.895, "15,000 t/day baseline; pressate concentrated by membrane + evaporator",
             color=INK_2, fontsize=10, ha="left", va="center")

    # Side note: the current assumption, with the sugar composition at that point.
    nx = 0.71
    fig.text(nx, 0.78, "Current assumption", color=INK, fontsize=10.5,
             fontweight="semibold", ha="left", va="center")
    fig.text(nx, 0.735, f"{CURRENT_TARGET * 100:.0f} wt% total solids", color=INK_2,
             fontsize=10, ha="left", va="center")
    fig.text(nx, 0.695, f"{cur['sugar_g_per_L']:.1f} g/L sugars:", color=INK_2,
             fontsize=10, ha="left", va="center")
    names = {"Glucose": "Glucose", "Mannitol": "Mannitol", "AlginateMonomer": "Alginate monomer"}
    for k, sid in enumerate(SUGAR_IDS):
        y = 0.657 - 0.036 * k
        conc = cur[f"{sid}_g_per_L"]
        fig.text(nx + 0.012, y, names.get(sid, sid), color=INK_2, fontsize=9,
                 ha="left", va="center")
        fig.text(nx + 0.205, y, f"{conc:.1f} g/L", color=INK_2, fontsize=9,
                 ha="right", va="center")
        fig.text(nx + 0.270, y, f"{100 * conc / cur['sugar_g_per_L']:.1f}%", color=INK_2,
                 fontsize=9, ha="right", va="center")
    y_titer = 0.657 - 0.036 * len(SUGAR_IDS) - 0.012
    fig.text(nx, y_titer, f"{cur['titer_g_per_L']:.1f} g/L titer", color=INK_2,
             fontsize=10, ha="left", va="center")
    y_msp = y_titer - 0.04
    fig.text(nx, y_msp, f"MSP ${cur['msp_usd_per_kg']:.2f}/kg", color=INK_2,
             fontsize=10, ha="left", va="center")
    fig.text(nx, y_msp - 0.075,
             "Shares are % of sugar mass.\nNo inhibition or titer limit\nis modeled, so high-solids\npoints may not be fermentable.",
             color=INK_2, fontsize=8.5, ha="left", va="top", linespacing=1.35)

    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    df = pd.DataFrame([run_target(float(t)) for t in TARGETS])
    csv_path = TABLES_OUT / "msp_fermentation_solids_3hp.csv"
    fig_path = FIGURES_OUT / "fig_msp_fermentation_solids_3hp.png"

    df.to_csv(csv_path, index=False)
    plot_sweep(df, fig_path)

    print(f"Saved table:  {csv_path}")
    print(f"Saved figure: {fig_path}\n")
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(df.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))

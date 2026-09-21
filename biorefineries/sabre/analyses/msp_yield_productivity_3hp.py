# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_yield_productivity_3hp.py
-----------------------------
One-at-a-time sensitivity of the 3-HP (dried Ca(3HP)2) MSP to the
fermentation's yield and to its productivity, at the 15,000 wet ton/day
baseline and the current 20 wt% fermentation-feed solids target.

Yield sweep (0.30-0.80 g 3-HP / g sugar):
    "Yield" is the overall 3-HP mass formed per mass of fermentable sugar
    supplied (Glucose + Mannitol + AlginateMonomer), the same for all three
    sugars. HPFermentation forms 3-HP as conversion * mass yield per sugar
    (data/3hp.yaml fermentation.substrates), so each sugar's `hp_conversion`
    is set to `target / hp_yield_kg_per_kg_consumed` on the F401 instance
    after the build. The current model's per-sugar split (0.49 glucose /
    0.396 mannitol / 0.371 alginate monomer) is a weighted 0.427 overall at
    the baseline feed composition; that unmodified run is the "current" point.
    The titer is not an input: it follows from yield and the feed sugar
    concentration, so it is reported alongside.

Productivity sweep (0.3-1.5 g/L/h):
    Set on `F401.productivity_g_per_L_per_h`. It only enters through the
    batch time (tau = titer / productivity), i.e. fermenter size and cost;
    the yield is left at the current per-sugar values.

Neither the system builder nor the yaml files are edited.

Outputs (relative to this script's directory):
    ../results/tables/msp_yield_productivity_3hp.csv
    ../results/figures/fig_msp_yield_productivity_3hp.png

Run from the repo root:
    python biorefineries/sabre/analyses/msp_yield_productivity_3hp.py
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

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

_3HP_YAML = load_assumptions("3hp.yaml")
CURRENT_PRODUCTIVITY = _3HP_YAML["fermentation"]["productivity_g_per_L_per_h"]

YIELDS = np.round(np.arange(0.30, 0.8001, 0.05), 4)
# 0.55 = Nature Communications 2026 (I. orientalis, fed-batch); 0.76 = current model.
PRODUCTIVITIES = sorted({*np.round(np.arange(0.3, 1.5001, 0.1), 4), 0.55, CURRENT_PRODUCTIVITY})


def run_case(yield_g_per_g: float | None = None, productivity_g_per_L_per_h: float | None = None) -> dict:
    """Build/simulate the baseline system with an optional yield / productivity override."""
    bst.main_flowsheet.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sys = _3hp_system.create_3hp_system()
        sys.set_tolerance(subsystems=True, **TIGHT_TOLERANCE)
        F401 = sys.flowsheet.unit.F401
        if yield_g_per_g is not None:
            for params in F401.substrates.values():
                X = yield_g_per_g / params["hp_yield_kg_per_kg_consumed"]
                if X + params["biomass_conversion"] > 1.0:
                    raise ValueError(f"yield {yield_g_per_g} needs conversion {X:.3f}, > 1 with biomass")
                params["hp_conversion"] = X
        if productivity_g_per_L_per_h is not None:
            F401.productivity_g_per_L_per_h = productivity_g_per_L_per_h
        sys.simulate()
        msp = solve_product_msp(sys.TEA, sys.flowsheet.stream.ca3hp2_product)

    dr = F401.design_results
    return {
        "yield_g_per_g": dr["HP formed (kg/h)"] / dr["Substrate available (kg/h)"],
        "productivity_g_per_L_per_h": F401.productivity_g_per_L_per_h,
        "titer_g_per_L": float(dr["Titer (g HP-equivalent/L)"]),
        "fermentation_batch_h": float(F401.tau),
        "F401_installed_usd": float(F401.installed_cost),
        "FCI_usd": sys.TEA.FCI,
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
        "msp_usd_per_kg": msp["usd_per_kg"],
    }


def collect_results():
    current = {"case": "current", **run_case()}
    yield_rows = [{"case": "yield", **run_case(yield_g_per_g=float(y))} for y in YIELDS]
    prod_rows = [{"case": "productivity", **run_case(productivity_g_per_L_per_h=float(p))} for p in PRODUCTIVITIES]
    return current, pd.DataFrame(yield_rows), pd.DataFrame(prod_rows)


def _style_axes(ax, ylim, yticks):
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, linestyle="-")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(axis="both", length=0, colors=INK_2, labelsize=9.5)


def _plot_line(ax, x, y, x_cur, y_cur, y_span):
    ax.plot(x, y, color=SERIES_WITH, linewidth=2.0, solid_capstyle="round",
            solid_joinstyle="round", zorder=3)
    ax.axvline(x_cur, color=AXIS, linewidth=1.0, zorder=1)
    ax.plot([x_cur], [y_cur], marker="o", markersize=9, color=SERIES_WITH,
            markeredgecolor=SURFACE, markeredgewidth=2, linestyle="none", zorder=4)
    ax.text(x_cur, y_cur + 0.05 * y_span, "current", ha="center", va="bottom",
            color=INK_2, fontsize=9, zorder=5)
    for xi, yi, ha in ((x.iloc[0], y.iloc[0], "left"), (x.iloc[-1], y.iloc[-1], "right")):
        ax.text(xi, yi + 0.05 * y_span, f"{yi:.2f}", ha=ha, va="bottom",
                color=INK, fontsize=10, zorder=5)


def plot_sensitivity(current: dict, ydf: pd.DataFrame, pdf: pd.DataFrame, out_path: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "figure.dpi": 100, "savefig.dpi": 200,
    })
    fig = plt.figure(figsize=(10.4, 5.4), facecolor=SURFACE)
    ax_y = fig.add_axes([0.075, 0.24, 0.41, 0.56], facecolor=SURFACE)
    ax_p = fig.add_axes([0.565, 0.24, 0.41, 0.56], facecolor=SURFACE, sharey=ax_y)

    y_top = float(np.ceil(max(ydf["msp_usd_per_kg"].max(), pdf["msp_usd_per_kg"].max()) * 1.12 * 2) / 2)
    ylim = (0, y_top)
    yticks = np.arange(0, y_top + 0.01, 0.5)
    _style_axes(ax_y, ylim, yticks)
    _style_axes(ax_p, ylim, yticks)
    plt.setp(ax_p.get_yticklabels(), visible=False)
    ax_y.set_ylabel(r"MSP (USD/kg Ca(3HP)$_2$)", color=INK_2, fontsize=10)
    y_span = ylim[1] - ylim[0]

    # -- yield panel; second tick row = implied titer -------------------------
    _plot_line(ax_y, ydf["yield_g_per_g"], ydf["msp_usd_per_kg"],
               current["yield_g_per_g"], current["msp_usd_per_kg"], y_span)
    ax_y.set_xlim(0.27, 0.83)
    tenths = ydf["yield_g_per_g"] * 10
    tick_rows = ydf[np.isclose(tenths, np.round(tenths), atol=1e-6)]
    ax_y.set_xticks(tick_rows["yield_g_per_g"])
    ax_y.set_xticklabels([f"{a:.1f}\n{b:.0f} g/L" for a, b in
                          zip(tick_rows["yield_g_per_g"], tick_rows["titer_g_per_L"])], color=INK_2)
    ax_y.tick_params(axis="x", pad=6)
    ax_y.set_xlabel("Overall 3-HP yield (g / g sugar)", color=INK_2, fontsize=10, labelpad=8)

    # -- productivity panel; second tick row = batch time ---------------------
    _plot_line(ax_p, pdf["productivity_g_per_L_per_h"], pdf["msp_usd_per_kg"],
               current["productivity_g_per_L_per_h"], current["msp_usd_per_kg"], y_span)
    ax_p.set_xlim(0.2, 1.6)
    ticks = [0.3, 0.55, CURRENT_PRODUCTIVITY, 1.0, 1.5]
    sel = pdf.set_index("productivity_g_per_L_per_h").loc[ticks]
    ax_p.set_xticks(ticks)
    ax_p.set_xticklabels([f"{t:g}\n{h:.0f} h" for t, h in zip(ticks, sel["fermentation_batch_h"])], color=INK_2)
    ax_p.tick_params(axis="x", pad=6)
    ax_p.set_xlabel("Fermentation productivity (g/L/h)", color=INK_2, fontsize=10, labelpad=8)

    fig.text(0.075, 0.945, "Ca(3HP)$_2$ MSP sensitivity to fermentation yield and productivity",
             color=INK, fontsize=12.5, fontweight="semibold", ha="left", va="center")
    fig.text(0.075, 0.895, "15,000 t/day baseline, 20 wt% fermentation-feed solids; one parameter varied at a time",
             color=INK_2, fontsize=10, ha="left", va="center")
    fig.text(0.075, 0.835, f"Yield, productivity fixed at {CURRENT_PRODUCTIVITY:g} g/L/h", color=INK_2,
             fontsize=10, ha="left", va="center")
    fig.text(0.565, 0.835, f"Productivity, yield fixed at current {current['yield_g_per_g']:.2f} g/g",
             color=INK_2, fontsize=10, ha="left", va="center")
    fig.text(0.075, 0.045,
             "Lower tick rows: implied titer (yield panel) and batch time (productivity panel). "
             "No titer limit is modeled.",
             color=INK_2, fontsize=8.5, ha="left", va="center")

    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    current, ydf, pdf = collect_results()
    csv_path = TABLES_OUT / "msp_yield_productivity_3hp.csv"
    fig_path = FIGURES_OUT / "fig_msp_yield_productivity_3hp.png"

    pd.concat([pd.DataFrame([current]), ydf, pdf], ignore_index=True).to_csv(csv_path, index=False)
    plot_sensitivity(current, ydf, pdf, fig_path)

    print(f"Saved table:  {csv_path}")
    print(f"Saved figure: {fig_path}\n")
    cols = ["yield_g_per_g", "productivity_g_per_L_per_h", "titer_g_per_L", "fermentation_batch_h",
            "F401_installed_usd", "annual_product_t", "msp_usd_per_kg"]
    with pd.option_context("display.width", 220, "display.max_columns", None):
        fmt = lambda v: f"{v:,.4g}"
        print("current:\n", pd.DataFrame([current])[cols].to_string(index=False, float_format=fmt), "\n")
        print("yield sweep:\n", ydf[cols].to_string(index=False, float_format=fmt), "\n")
        print("productivity sweep:\n", pdf[cols].to_string(index=False, float_format=fmt))

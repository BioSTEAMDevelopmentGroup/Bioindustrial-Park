# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_scale_scenarios_3hp.py
--------------------------
3-HP (dried Ca(3HP)2) MSP under three plant-scale / campaign-length
scenarios, each with and without the solid-disposal cost:

    baseline    15,000 wet ton/day, 330 operating days/yr (data/feedstock.yaml
                + data/tea.yaml as-is)
    330 days    10,000 wet ton/yr spread over 330 operating days (~30 ton/day)
    100 days    10,000 wet ton/yr processed in 100 operating days (100 ton/day)

"Solid disposal cost" is the one solid-waste cost in this system: the
BoilerTurbogenerator's ash disposal (data/tea.yaml `price.disposal_solid`,
set on `BT.ash_disposal_price` by systems._3hp_system). pressed_cake,
milling_losses and cell_mass are burned in that boiler rather than priced
as disposal, and the liquid purges are priced as wastewater, not solids.
"Without" sets the ash disposal price to 0; nothing else changes.

Scaling is applied from here, without editing the system builder:
  * feed rate: systems._3hp_system.get_scale_feed_kgph is swapped for a
    fixed value while the system is built.
  * operating days: `tea.operating_days` is set after the build, and
    units._biostimulant.OPERATING_HOURS_PER_YEAR (which
    PressateConcentrator uses to turn its annual maintenance cost into an
    hourly add_OPEX) is swapped to match, so that annual maintenance is
    still `maintenance_frac_of_capex_per_yr * capex` rather than being
    silently rescaled by operating_days / 330.
  * flowsheet tolerances are tightened (set_tolerance) for every scenario.
    BioSTEAM's defaults (1 kmol/hr absolute or 1 % relative) are loose
    enough at ~30 ton/day for the mother-liquor recycle loop to stop
    short of steady state, giving a ~4 % different product output for the
    same annual feed.

Outputs (relative to this script's directory):
    ../results/tables/msp_scale_scenarios_3hp.csv
    ../results/figures/fig_msp_scale_scenarios_3hp.png

Run from the repo root:
    python biorefineries/sabre/analyses/msp_scale_scenarios_3hp.py
"""

import re
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import biosteam as bst

from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.systems import _3hp_system
from biorefineries.sabre.units import _biostimulant
from biorefineries.sabre.utils import load_assumptions

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

_TEA_YAML = load_assumptions("tea.yaml")
_FEEDSTOCK_YAML = load_assumptions("feedstock.yaml")
BASELINE_DAYS = _TEA_YAML["operating_days"]
BASELINE_TPD = _FEEDSTOCK_YAML["scale"]["wet_feed_ton_per_day"]
ASH_DISPOSAL_PRICE = _TEA_YAML["price"]["disposal_solid"]["baseline"]  # USD/kg, negative = cost
ANNUAL_FEED_SMALL = 10_000  # wet ton/yr

# (key, feed [wet ton/day], operating days/yr). Feed is metric, matching
# data/feedstock.yaml's `ton_definition`.
SCENARIOS = (
    ("baseline", BASELINE_TPD, BASELINE_DAYS),
    ("330 days", ANNUAL_FEED_SMALL / BASELINE_DAYS, BASELINE_DAYS),
    ("100 days", ANNUAL_FEED_SMALL / 100, 100),
)

# Relative + absolute molar tolerance for every flowsheet loop -- see module docstring.
TIGHT_TOLERANCE = dict(mol=1e-6, rmol=1e-6)


@contextmanager
def _scaled_build(feed_tpd: float, operating_days: float):
    """Swap the module-level scale inputs for the duration of a system build."""
    saved_feed = _3hp_system.get_scale_feed_kgph
    saved_hours = _biostimulant.OPERATING_HOURS_PER_YEAR
    saved_ash_price = bst.stream_utility_prices["Ash disposal"]
    _3hp_system.get_scale_feed_kgph = lambda A: feed_tpd * 1000.0 / 24.0
    _biostimulant.OPERATING_HOURS_PER_YEAR = operating_days * 24.0
    try:
        yield
    finally:
        _3hp_system.get_scale_feed_kgph = saved_feed
        _biostimulant.OPERATING_HOURS_PER_YEAR = saved_hours
        bst.stream_utility_prices["Ash disposal"] = saved_ash_price


def run_scenario(feed_tpd: float, operating_days: float, ash_disposal_price: float) -> dict:
    """Build, simulate, and solve the Ca(3HP)2 MSP at one scale/campaign/ash-price."""
    bst.main_flowsheet.clear()
    with _scaled_build(feed_tpd, operating_days), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        sys = _3hp_system.create_3hp_system()
        sys.flowsheet.unit.BT.ash_disposal_price = ash_disposal_price
        sys.TEA.operating_days = operating_days
        sys.set_tolerance(subsystems=True, **TIGHT_TOLERANCE)
        sys.simulate()
        msp = solve_product_msp(sys.TEA, sys.flowsheet.stream.ca3hp2_product)

    # BioSTEAM's cost/design correlations warn when a unit is sized outside
    # the range they were fitted on -- collect which units, for the caveat.
    out_of_range = sorted({
        m.group(1)
        for w in caught
        for m in [re.search(r"<\w+: ([^>\x1b]+?)(?:\x1b\[[0-9;]*m)?>", str(w.message))]
        if m and "out of bounds" in str(w.message)
    })
    return {
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
        "ash_kg_per_h": float(sys.flowsheet.unit.BT.outs[2].F_mass),
        "FCI_usd": sys.TEA.FCI,
        "out_of_range_units": out_of_range,
    }


def collect_results() -> pd.DataFrame:
    rows = []
    for key, feed_tpd, days in SCENARIOS:
        with_cost = run_scenario(feed_tpd, days, ASH_DISPOSAL_PRICE)
        without_cost = run_scenario(feed_tpd, days, 0.0)
        rows.append({
            "scenario": key,
            "feed_wet_ton_per_day": feed_tpd,
            "operating_days": days,
            "annual_feed_wet_ton": feed_tpd * days,
            "annual_product_t": with_cost["annual_product_t"],
            "ash_kg_per_h": with_cost["ash_kg_per_h"],
            "FCI_usd": with_cost["FCI_usd"],
            "msp_with_solid_disposal_usd_per_kg": with_cost["msp_usd_per_kg"],
            "msp_without_solid_disposal_usd_per_kg": without_cost["msp_usd_per_kg"],
            "units_outside_cost_correlation_range": " ".join(with_cost["out_of_range_units"]),
        })
    return pd.DataFrame(rows)


# -- plotting ---------------------------------------------------------------

# Chart tokens (light mode), from the dataviz reference palette.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES_WITH = "#2a78d6"     # categorical slot 1 (blue)
SERIES_WITHOUT = "#eb6834"  # categorical slot 2 (orange)

FIG_W_IN, FIG_H_IN = 7.4, 4.7
BAR_W_IN = 0.40
GAP_IN = 2 / 96      # 2 px surface gap between adjacent bars
RADIUS_IN = 4 / 96   # 4 px rounded data-end


def _rounded_column(ax, x_center, height, width, color, rx, ry):
    """
    Column with a rounded top and a square base at y=0. The rounded box is
    drawn extending below the baseline and clipped to y >= 0, so only the
    top corners keep their radius. rx/ry are the same radius in data units
    on each axis (they differ because the axes aren't square).
    """
    x0 = x_center - width / 2
    bottom = -3 * ry
    patch = FancyBboxPatch(
        (x0, bottom), width, height - bottom,
        boxstyle=f"round,pad=0,rounding_size={rx}",
        mutation_aspect=ry / rx,
        facecolor=color, edgecolor="none", linewidth=0, zorder=3,
    )
    ax.add_patch(patch)
    patch.set_clip_path(Rectangle((x0, 0), width, height + ry, transform=ax.transData))


def _scenario_label(row) -> str:
    tpd = row["feed_wet_ton_per_day"]
    days = int(row["operating_days"])
    annual = row["annual_feed_wet_ton"]
    if row["scenario"] == "baseline":
        return f"{tpd:,.0f} t/day (baseline)\n({days} d/yr, {annual / 1e6:.2f} Mt/yr)"
    return f"{annual:,.0f} t/yr over {days} d\n({tpd:,.0f} t/day)"


def plot_msp_columns(df: pd.DataFrame, out_path: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "figure.dpi": 100,
        "savefig.dpi": 200,
    })
    fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN), facecolor=SURFACE)
    left, right, bottom, top = 0.10, 0.97, 0.16, 0.77
    ax = fig.add_axes([left, bottom, right - left, top - bottom], facecolor=SURFACE)

    n = len(df)
    y_max = 12.0
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, y_max)

    ax_w_in = FIG_W_IN * (right - left)
    ax_h_in = FIG_H_IN * (top - bottom)
    x_per_in = (n) / ax_w_in
    y_per_in = y_max / ax_h_in
    bar_w = BAR_W_IN * x_per_in
    offset = (BAR_W_IN + GAP_IN) / 2 * x_per_in
    rx, ry = RADIUS_IN * x_per_in, RADIUS_IN * y_per_in

    # Hairline grid + baseline, recessive.
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_yticks(range(0, int(y_max) + 1, 2))
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="both", length=0, colors=INK_2, labelsize=9.5)
    ax.tick_params(axis="x", pad=8)
    ax.set_xticks(range(n))
    ax.set_xticklabels([_scenario_label(r) for _, r in df.iterrows()], color=INK_2)
    ax.set_ylabel(r"MSP (USD/kg Ca(3HP)$_2$)", color=INK_2, fontsize=10)

    col_with = "msp_with_solid_disposal_usd_per_kg"
    col_without = "msp_without_solid_disposal_usd_per_kg"
    for i, (_, r) in enumerate(df.iterrows()):
        for dx, col, color in ((-offset, col_with, SERIES_WITH), (+offset, col_without, SERIES_WITHOUT)):
            v = r[col]
            _rounded_column(ax, i + dx, v, bar_w, color, rx, ry)
            ax.text(i + dx, v + 0.18, f"{v:.2f}", ha="center", va="bottom",
                    color=INK, fontsize=10, zorder=4)

    fig.text(left, 0.94, "Ca(3HP)$_2$ MSP by plant scale and campaign length",
             color=INK, fontsize=12.5, fontweight="semibold", ha="left", va="center")
    fig.text(left, 0.885, "Sargassum to dried Ca(3HP)$_2$; solid disposal = boiler ash",
             color=INK_2, fontsize=10, ha="left", va="center")

    handles = [
        Patch(facecolor=SERIES_WITH, edgecolor="none",
              label=f"With solid disposal cost (${abs(ASH_DISPOSAL_PRICE):.2f}/kg ash)"),
        Patch(facecolor=SERIES_WITHOUT, edgecolor="none", label="Without solid disposal cost"),
    ]
    leg = fig.legend(handles=handles, loc="center left", bbox_to_anchor=(left - 0.005, 0.83),
                     ncol=2, frameon=False, fontsize=9.5, handlelength=1.0, handleheight=1.0,
                     handletextpad=0.5, columnspacing=1.8)
    for t in leg.get_texts():
        t.set_color(INK_2)

    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    df = collect_results()
    csv_path = TABLES_OUT / "msp_scale_scenarios_3hp.csv"
    fig_path = FIGURES_OUT / "fig_msp_scale_scenarios_3hp.png"

    df.to_csv(csv_path, index=False)
    plot_msp_columns(df, fig_path)

    print(f"Saved table:  {csv_path}")
    print(f"Saved figure: {fig_path}\n")
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(df.drop(columns="units_outside_cost_correlation_range").to_string(index=False))
    print()
    for _, r in df.iterrows():
        print(f"{r['scenario']}: units outside cost-correlation range: "
              f"{r['units_outside_cost_correlation_range'] or 'none'}")

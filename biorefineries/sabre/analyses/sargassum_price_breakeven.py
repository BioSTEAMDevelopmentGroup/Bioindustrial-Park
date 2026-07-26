# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
sargassum_price_breakeven.py
-----------------------------
For the four standalone SaBRe flowsheets (the 'integrated' flowsheet is
excluded, same as `msp_comparison.py`/`disposal_price_breakeven.py`), solve
the sargassum feedstock price ($/kg wet) that hits the target IRR (data/
tea.yaml `tea.IRR`), holding every product at its data/tea.yaml baseline
selling price:

    - Biostimulant     -> biostimulant_product
    - AD-biomethane    -> biomethane, data/pretreatment.yaml
                          `pretreatment_ad` case 'combined_PE' only (same
                          case cost_breakdown.py/disposal_price_breakeven.py
                          use as their single AD-biomethane representative)
    - AD-VFA           -> pure_vfa
    - AD-fermentation  -> microbial_oil

This is the feedstock-side mirror of `msp_comparison.py`: there, product
price is solved holding sargassum_feed's price fixed at its tea.yaml
baseline; here, sargassum_feed's price is solved holding the pathway's own
product price fixed at its own tea.yaml baseline (biostimulant priced as a
co-product on each AD pathway's combined TEA, same as `msp_comparison.py`'s
`credit_tipping_fee=False` basis -- i.e. no tipping-fee credit from
biostimulant's pressed_cake disposal cost).

A solved sargassum price above the tea.yaml baseline (0.0) means the
pathway can afford to pay for feedstock and still hit target IRR at that
product price; a solved price below baseline (i.e. the pathway would need
to be paid, a more negative "tipping fee") means the pathway can't hit
target IRR at the baseline product price without a feedstock subsidy.

A second figure (run_sargassum_vs_product_price_four_systems()) sweeps this
same break-even sargassum price across each pathway's own product over its
own data/tea.yaml price range (`price.biostimulant.range`, `price.vfa.
range`, `price.microbial_oil.range`, `price.biomethane_mmbtu.range`), all
four pathways as one figure/one set of axes. Since each pathway's product
has a different unit and price scale, the x-axis (product price) is
normalized to % of that pathway's own tea.yaml baseline -- so 100% always
means "this pathway's own product at its own tea.yaml baseline," for every
line. The y-axis (break-even sargassum price, $/kg) is already common to
every pathway, so it isn't normalized.

Run from the repo root:
    python biorefineries/sabre/analyses/sargassum_price_breakeven.py
"""

import sys
from pathlib import Path

import numpy as np
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
from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre._tea import usd_per_mmbtu_to_usd_per_kg, CH4_MMBTU_PER_KG

OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

_TEA_PRICE = load_assumptions("tea.yaml")["price"]
SARGASSUM_BASELINE = _TEA_PRICE["sargassum"]["baseline"]
SARGASSUM_RANGE = _TEA_PRICE["sargassum"]["range"]

# AD-biomethane pretreatment case (data/pretreatment.yaml `pretreatment_ad`),
# same single representative case cost_breakdown.py/disposal_price_breakeven.py use.
BIOMETHANE_PRETREATMENT_CASE = "combined_PE"

PRODUCT_STREAM_BY_LABEL = {
    "Biostimulant":     "biostimulant_product",
    "AD-VFA":           "pure_vfa",
    "AD-fermentation":  "microbial_oil",
    "AD-biomethane":    "biomethane",
}


def solve_sargassum_breakeven(sys: bst.System, label: str) -> dict:
    """
    Solve the sargassum_feed price ($/kg wet) that zeroes this system's NPV
    at its TEA's target IRR, holding every already-priced product/reagent/
    waste stream fixed at its current price (i.e. the data/tea.yaml
    baseline each create_*_system() sets by default).
    """
    fw = sys.flowsheet.stream
    feed = fw.sargassum_feed

    breakeven_usd_per_kg = sys.TEA.solve_price(feed)

    annual_hours = sys.TEA.operating_days * 24
    annual_feed_kg = float(feed.F_mass) * annual_hours

    product = getattr(fw, PRODUCT_STREAM_BY_LABEL[label])

    return {
        "label": label,
        "breakeven_usd_per_kg": breakeven_usd_per_kg,
        "annual_feed_kg": annual_feed_kg,
        "product_stream_name": product.ID,
        "product_price_usd_per_kg": float(product.price),
        "sys": sys,
    }


# =================================================================
# Main
# =================================================================

def run_sargassum_price_breakeven() -> None:
    results = []

    bst.main_flowsheet.clear()
    sys = create_biostimulant_system()
    sys.simulate()
    results.append(solve_sargassum_breakeven(sys, "Biostimulant"))

    bst.main_flowsheet.clear()
    sys = create_ad_biomethane_system(pretreatment_case=BIOMETHANE_PRETREATMENT_CASE)
    sys.simulate()
    results.append(solve_sargassum_breakeven(sys, "AD-biomethane"))

    bst.main_flowsheet.clear()
    sys = create_ad_vfa_system()
    sys.simulate()
    results.append(solve_sargassum_breakeven(sys, "AD-VFA"))

    bst.main_flowsheet.clear()
    sys = create_ad_fermentation_system()
    sys.simulate()
    results.append(solve_sargassum_breakeven(sys, "AD-fermentation"))

    for r in results:
        r["bar_label"] = r["label"]

    all_results = results

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
        "ytick.minor.visible": True,
        "xtick.minor.visible": True,
    })

    fig, ax = plt.subplots(figsize=(10.5, 5.5))

    tick_labels = [r["bar_label"] for r in all_results]
    vals = [r["breakeven_usd_per_kg"] for r in all_results]
    colors = [
        "#4C72B0" if r["label"] != "AD-biomethane" else "#C86E5A" for r in all_results
    ]

    bars = ax.bar(tick_labels, vals, edgecolor="black", linewidth=0.8, zorder=3, color=colors)

    span = max(vals + [0]) - min(vals + [0]) or 1
    for bar, val in zip(bars, vals):
        above = val >= 0
        x = bar.get_x() + bar.get_width() / 2
        y = bar.get_height() + (span * 0.02 if above else -span * 0.02)
        ax.text(x, y, f"${val:,.3f}", ha="center", va="bottom" if above else "top", fontsize=9)

    ax.axhline(
        SARGASSUM_BASELINE, color="black", linewidth=1.0, linestyle="--", zorder=2,
    )
    ax.text(
        -0.5, SARGASSUM_BASELINE,
        f"tea.yaml baseline (${SARGASSUM_BASELINE:.2f}/kg) ",
        ha="left", va="bottom", fontsize=8, color="#444444",
    )

    ax.set_ylabel("Break-even sargassum feedstock price ($/kg wet)")
    ax.set_title(
        "Sargassum feedstock price for break-even (target IRR), products at tea.yaml baseline",
        fontsize=10,
    )
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.3f"))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))
    ax.xaxis.set_minor_locator(mticker.NullLocator())  # x is categorical, no minor ticks
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.axhline(0, color="black", linewidth=0.8, zorder=3)
    lo, hi = min(0, min(vals), SARGASSUM_BASELINE), max(0, max(vals), SARGASSUM_BASELINE)
    pad = (hi - lo) * 0.22 or 1
    ax.set_ylim(lo - pad, hi + pad)
    ax.tick_params(axis="x", labelrotation=15)

    fig.tight_layout()
    out_path = OUT / "fig_sargassum_price_breakeven_four_systems.png"
    fig.savefig(out_path, bbox_inches="tight")

    print("\nSaved:")
    print(f"  {out_path}")

    # ── Summary table ────────────────────────────────────────────────────
    print(
        f"\nSargassum feedstock break-even price by pathway "
        f"(products at tea.yaml baseline; tea.yaml sargassum baseline = "
        f"${SARGASSUM_BASELINE:.4f}/kg, range {SARGASSUM_RANGE}):"
    )
    header = (
        f"  {'Pathway':<38}{'Product @ price':<26}"
        f"{'Break-even ($/kg)':>20}{'Headroom vs baseline':>22}{'Annual feed (t/yr)':>22}"
    )
    print(header)
    for r in all_results:
        label = (
            r["label"] if r["label"] != "AD-biomethane"
            else f"AD-biomethane ({BIOMETHANE_PRETREATMENT_CASE})"
        )
        product_str = f"{r['product_stream_name']} @ ${r['product_price_usd_per_kg']:.3f}"
        headroom = r["breakeven_usd_per_kg"] - SARGASSUM_BASELINE
        print(
            f"  {label:<38}{product_str:<26}"
            f"{r['breakeven_usd_per_kg']:>19.4f}{headroom:>+21.4f}"
            f"{r['annual_feed_kg'] / 1000:>21,.1f}"
        )
    print(
        "  (Headroom > 0: pathway can afford a more expensive sargassum feed "
        "than tea.yaml baseline and still hit target IRR at the given product "
        "price -- i.e. it clears target IRR at baseline with room to spare. "
        "Headroom < 0: pathway needs a cheaper (or subsidized/tipping-fee) "
        "feedstock than tea.yaml baseline to hit target IRR.)"
    )


# Per-pathway sweep configuration: which product to vary, over which
# data/tea.yaml range (in the product's own natural unit), and how to
# convert that unit to the $/kg biosteam actually prices streams in
# (None means the natural unit already is $/kg).
SWEEP_CONFIGS = (
    {
        "label": "Biostimulant",
        "build_fn": create_biostimulant_system,
        "product_stream_name": "biostimulant_product",
        "x_range": tuple(_TEA_PRICE["biostimulant"]["range"]),
        "x_baseline": _TEA_PRICE["biostimulant"]["baseline"],
        "x_label": "Biostimulant selling price ($/kg product)",
        "to_usd_per_kg": None,
    },
    {
        "label": "AD-VFA",
        "build_fn": create_ad_vfa_system,
        "product_stream_name": "pure_vfa",
        "x_range": tuple(_TEA_PRICE["vfa"]["range"]),
        "x_baseline": _TEA_PRICE["vfa"]["baseline"],
        "x_label": "VFA selling price ($/kg total VFA)",
        "to_usd_per_kg": None,
    },
    {
        "label": "AD-fermentation",
        "build_fn": create_ad_fermentation_system,
        "product_stream_name": "microbial_oil",
        "x_range": tuple(_TEA_PRICE["microbial_oil"]["range"]),
        "x_baseline": _TEA_PRICE["microbial_oil"]["baseline"],
        "x_label": "Microbial oil selling price ($/kg crude oil)",
        "to_usd_per_kg": None,
    },
    {
        "label": "AD-biomethane",
        "build_fn": lambda: create_ad_biomethane_system(pretreatment_case=BIOMETHANE_PRETREATMENT_CASE),
        "product_stream_name": "biomethane",
        "x_range": tuple(_TEA_PRICE["biomethane_mmbtu"]["range"]),
        "x_baseline": _TEA_PRICE["biomethane_mmbtu"]["baseline"],
        "x_label": "Biomethane selling price ($/mmbtu)",
        "to_usd_per_kg": lambda x: usd_per_mmbtu_to_usd_per_kg(x, CH4_MMBTU_PER_KG),
    },
)


def sweep_sargassum_vs_product_price(cfg: dict, n_points: int) -> dict:
    """
    Build cfg['label']'s system once, simulate once, then sweep
    cfg['product_stream_name']'s price over cfg['x_range'] (converted to
    $/kg via cfg['to_usd_per_kg'] if given) and solve the break-even
    sargassum_feed price at each point via solve_sargassum_breakeven()'s
    same TEA.solve_price() call. Re-simulation isn't needed between price
    points -- price only feeds the TEA's cost/sales bookkeeping, not the
    mass/energy balances solved by sys.simulate().
    """
    label = cfg["label"]
    bst.main_flowsheet.clear()
    sys = cfg["build_fn"]()
    sys.simulate()

    fw = sys.flowsheet.stream
    product = getattr(fw, cfg["product_stream_name"])
    feed = fw.sargassum_feed

    x_values = np.linspace(cfg["x_range"][0], cfg["x_range"][1], n_points)
    breakevens = []
    for x in x_values:
        product.price = float(x if cfg["to_usd_per_kg"] is None else cfg["to_usd_per_kg"](x))
        breakevens.append(sys.TEA.solve_price(feed))

    return {**cfg, "x_values": x_values, "breakevens": breakevens}


PATHWAY_COLORS = {
    "Biostimulant":     "#4C72B0",
    "AD-VFA":           "#55A868",
    "AD-fermentation":  "#8172B2",
    "AD-biomethane":    "#C86E5A",
}


def run_sargassum_vs_product_price_four_systems(n_points: int = 16) -> None:
    """
    For each of the four standalone pathways, sweep its own product's
    selling price over its own data/tea.yaml price range and solve the
    break-even sargassum feedstock price at each point. The x-axis (product
    price) is normalized to a % of that pathway's own tea.yaml baseline --
    despite having different products/units/price scales, all four
    pathways then share one x-axis, with 100% always meaning "this
    pathway's own product at its own tea.yaml baseline" (labeled with each
    pathway's actual baseline price in the legend). The y-axis is the
    actual break-even sargassum price ($/kg wet) -- common to every
    pathway already, so no normalization needed there.
    """
    sweep_results = [sweep_sargassum_vs_product_price(cfg, n_points) for cfg in SWEEP_CONFIGS]
    for r in sweep_results:
        r["x_pct"] = [x / r["x_baseline"] * 100 for x in r["x_values"]]

    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    for r in sweep_results:
        unit = "/mmbtu" if r["label"] == "AD-biomethane" else "/kg"
        legend_label = f"{r['label']} (baseline ${r['x_baseline']:.2f}{unit})"
        ax.plot(
            r["x_pct"], r["breakevens"], marker="o", markersize=4, linewidth=1.5,
            color=PATHWAY_COLORS[r["label"]], label=legend_label, zorder=3,
        )

    ax.axhline(SARGASSUM_BASELINE, color="black", linewidth=1.0, linestyle="--", zorder=2)
    ax.text(
        ax.get_xlim()[0], SARGASSUM_BASELINE,
        f" tea.yaml sargassum baseline (${SARGASSUM_BASELINE:.2f}/kg) ",
        ha="left", va="bottom", fontsize=8, color="#444444",
    )
    # Vertical line marks 100% (each pathway's own product at its own
    # tea.yaml baseline) -- the actual baseline price is in the legend
    # instead of a text label here, since it differs by pathway/unit.
    ax.axvline(100, color="black", linewidth=1.0, linestyle=":", zorder=2)

    ax.set_xlabel("Product selling price (% of its own tea.yaml baseline)")
    ax.set_ylabel("Break-even sargassum feedstock price ($/kg wet)")
    ax.set_title(
        "Break-even sargassum feedstock price vs. product selling price\n"
        "(product price normalized to each pathway's own tea.yaml baseline)",
        fontsize=11,
    )
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.0f%%"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.3f"))
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))
    ax.grid(linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.legend(fontsize=9, loc="best")

    fig.tight_layout()
    out_path = OUT / "fig_sargassum_price_vs_selling_price_pct_four_systems.png"
    fig.savefig(out_path, bbox_inches="tight")

    print("\nSaved:")
    print(f"  {out_path}")

    for r in sweep_results:
        print(
            f"\n{r['label']} break-even sargassum price vs. {r['product_stream_name']} "
            f"selling price (range {r['x_range']}, % of baseline ${r['x_baseline']:.3f}):"
        )
        print(
            f"  {'Price (% of baseline)':<24}{'Price':<14}"
            f"{'Break-even sargassum ($/kg)':>28}"
        )
        for x, x_pct, breakeven in zip(r["x_values"], r["x_pct"], r["breakevens"]):
            print(f"  {x_pct:<24.1f}{x:<14.4f}{breakeven:>28.4f}")


if __name__ == "__main__":
    run_sargassum_price_breakeven()
    run_sargassum_vs_product_price_four_systems()

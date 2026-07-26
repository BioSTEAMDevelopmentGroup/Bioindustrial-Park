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
    - AD-biomethane    -> biomethane, $/mmbtu, one bar per data/pretreatment.yaml
                          `pretreatment_ad` case (press_mill_only, enzymatic,
                          peroxide, combined_PE, combined_PTE)
    - AD-VFA           -> VFA broth, $/kg total VFA
    - AD-fermentation  -> crude microbial oil, $/kg MicrobialOil

The product of each pathway is different, MSP of each pathway is the price
per kg (or, for biomethane, per mmbtu) of the pathway's own value-carrying
product needed to hit the target IRR, assuming no revenue from any other
product.

All cases share one baseline economic basis as set in tea.yaml

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


# System builders for the single-case, $/kg pathways
KG_BUILDERS = (
    price_biostimulant_system,
    price_ad_vfa_system,
    price_ad_fermentation_system,
)

# AD-biomethane pretreatment cases (data/pretreatment.yaml `pretreatment_ad`),
# priced on a $/mmbtu basis
PRETREATMENT_CASES = (
    "press_mill_only",
    "enzymatic",
    "peroxide",
    "combined_PE",
    "combined_PTE",
)
PRETREATMENT_LABELS = {
    "press_mill_only": "Press + mill only",
    "enzymatic":       "Enzymatic",
    "peroxide":        "Peroxide",
    "combined_PE":     "Combined PE",
    "combined_PTE":    "Combined PTE",
}

# Final value-carrying product stream (sys.flowsheet.stream) for each AD
# pathway, used for the dry-feed mass-yield summary below. Biostimulant is
# excluded -- its own product IS the biostimulant_product stream, already
# covered by the nonwater-fraction/wet-yield summary.
PRODUCT_STREAM_BY_LABEL = {
    "AD-biomethane":   "biomethane",
    "AD-VFA":          "pure_vfa",
    "AD-fermentation": "microbial_oil",
}


def compute_product_purity(sys: bst.System, label: str) -> float:
    """
    Mass frac of the pathway's value-carrying component(s) in its own final
    product stream:
      - Biostimulant:     non-water mass frac of biostimulant_product
        (i.e. its solids content)
      - AD-VFA:           total-VFA mass frac of pure_vfa (VFA_IDs per the
        MF (VFAMicrofilter) unit)
      - AD-biomethane:    Methane mass frac of the biomethane stream
      - AD-fermentation:  microbial-oil-chemical mass frac of the
        microbial_oil stream (product_ID per the OE (OilExtraction) unit)
    """
    fw = sys.flowsheet.stream
    un = sys.flowsheet.unit

    if label == "Biostimulant":
        product = fw.biostimulant_product
        water = float(product.imass["Water"]) if "Water" in product.chemicals.IDs else 0.0
        value_mass = float(product.F_mass) - water
    elif label == "AD-VFA":
        product = fw.pure_vfa
        vfa_ids = un.MF.vfa_IDs
        value_mass = sum(
            float(product.imass[cid]) for cid in vfa_ids if cid in product.chemicals.IDs
        )
    elif label == "AD-biomethane":
        product = fw.biomethane
        value_mass = float(product.imass["Methane"]) if "Methane" in product.chemicals.IDs else 0.0
    elif label == "AD-fermentation":
        product = fw.microbial_oil
        product_id = un.OE.product_ID
        value_mass = float(product.imass[product_id]) if product_id in product.chemicals.IDs else 0.0
    else:
        raise ValueError(f"Unknown label {label!r}")

    total_mass = float(product.F_mass)
    return value_mass / total_mass if total_mass else float("nan")


def summarize_mass_fractions(sys: bst.System, label: str) -> dict:
    """
    For a simulated pathway system, return:
      - product_yield_dry_pct: the pathway's final product mass as a % of
        the initial dry (moisture-free) sargassum_feed mass, for AD
        pathways (see PRODUCT_STREAM_BY_LABEL); for the Biostimulant
        pathway, this is instead biostimulant_product mass as a % of the
        initial *wet* sargassum_feed mass (product_yield_is_wet is set to
        flag the differing basis)
      - product_purity_frac: see compute_product_purity()
    """
    fw = sys.flowsheet.stream
    feed = fw.sargassum_feed

    feed_mass = float(feed.F_mass)
    feed_water = float(feed.imass["Water"]) if "Water" in feed.chemicals.IDs else 0.0
    feed_dry_mass = feed_mass - feed_water

    result = {"label": label}

    if label == "Biostimulant":
        biostim_mass = float(fw.biostimulant_product.F_mass)
        result["product_yield_dry_pct"] = (
            biostim_mass / feed_mass * 100 if feed_mass else float("nan")
        )
        result["product_yield_is_wet"] = True
    else:
        product_stream_name = PRODUCT_STREAM_BY_LABEL.get(label)
        if product_stream_name is not None:
            product = getattr(fw, product_stream_name)
            result["product_yield_dry_pct"] = (
                float(product.F_mass) / feed_dry_mass * 100 if feed_dry_mass else float("nan")
            )

    result["product_purity_frac"] = compute_product_purity(sys, label)

    return result


# =================================================================
# Main
# =================================================================

def run_msp_comparison(credit_tipping_fee: bool) -> None:
    """
    Build both the $/kg (Biostimulant/AD-VFA/AD-fermentation) and $/mmbtu
    (AD-biomethane by pretreatment case) comparisons, save one figure, and
    print the mass-yield summary table.

    Parameters
    ----------
    credit_tipping_fee : bool
        Forwarded to price_ad_vfa_system()/price_ad_fermentation_system()/
        price_ad_biomethane_system(). If False, each AD pathway's product
        price is solved on the fixed data/tea.yaml assumption basis
        (pressed_cake arrives free, no disposal-cost interaction with
        biostimulant). If True, each AD pathway is additionally credited a
        tipping fee for taking on pressed_cake that biostimulant would
        otherwise have paid data/tea.yaml `price.disposal_solid.baseline`
        to dispose of.
    """
    # Each builder's bst.main_flowsheet.clear() call wipes the *shared*
    # stream registry of every previously-built system, so mass fractions
    # must be captured right after each system is built/simulated -- not
    # after the fact from a stored 'sys' reference, once a later builder
    # has already cleared it out from under us.
    kg_builders = (
        price_biostimulant_system,
        lambda: price_ad_vfa_system(credit_tipping_fee=credit_tipping_fee),
        lambda: price_ad_fermentation_system(credit_tipping_fee=credit_tipping_fee),
    )
    kg_results = []
    summary_rows = []
    for builder in kg_builders:
        r = builder()
        kg_results.append(r)
        summary_rows.append(summarize_mass_fractions(r["sys"], r["label"]))

    biomethane_results = []
    for case in PRETREATMENT_CASES:
        r = price_ad_biomethane_system(pretreatment_case=case, credit_tipping_fee=credit_tipping_fee)
        biomethane_results.append(r)
        row = summarize_mass_fractions(r["sys"], r["label"])
        row["label"] = f"AD-biomethane ({PRETREATMENT_LABELS[case]})"
        summary_rows.append(row)

    scenario_tag = "tipping_fee" if credit_tipping_fee else "fixed"
    scenario_desc = (
        "AD pathways credited a tipping fee for taking pressed_cake off biostimulant's hands"
        if credit_tipping_fee else
        "fixed data/tea.yaml assumption basis (biostimulant priced at its flat baseline, "
        "combined capital across biostimulant + AD units)"
    )

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
        "AD-VFA":           "$/kg total VFA",
        "AD-fermentation":  "$/kg crude oil",
    }

    fig, (ax_kg, ax_biomethane) = plt.subplots(1, 2, figsize=(11.5, 5.2))

    def _fmt_usd_per_yr(v: float) -> str:
        """Format an annual USD figure with a K/M/B suffix, e.g. '$4.20M'."""
        av = abs(v)
        if av >= 1e9:
            return f"${v / 1e9:,.2f}B"
        if av >= 1e6:
            return f"${v / 1e6:,.2f}M"
        if av >= 1e3:
            return f"${v / 1e3:,.1f}K"
        return f"${v:,.0f}"

    def _label_bars(ax, bars, vals, sales_vals):
        span = max(vals + [0]) - min(vals + [0]) or 1
        for bar, val, sales in zip(bars, vals, sales_vals):
            above = val >= 0
            x = bar.get_x() + bar.get_width() / 2
            msp_y = bar.get_height() + (span * 0.015 if above else -span * 0.015)
            sales_y = bar.get_height() + (span * 0.075 if above else -span * 0.075)
            ax.text(
                x, msp_y, f"${val:,.2f}",
                ha="center", va="bottom" if above else "top", fontsize=9,
            )
            ax.text(
                x, sales_y, f"{_fmt_usd_per_yr(sales)}/yr sales",
                ha="center", va="bottom" if above else "top", fontsize=7,
                color="#444444",
            )

    # -- Left panel: single-case, $/kg pathways ------------------------------
    kg_tick_labels = [
        f"{r['label']}\n({UNIT_SUBLABELS[r['label']]})" for r in kg_results
    ]
    kg_vals = [r["msp_usd_per_kg"] for r in kg_results]
    kg_sales = [r["msp_usd_per_kg"] * r["annual_product_kg"] for r in kg_results]

    bars = ax_kg.bar(kg_tick_labels, kg_vals, edgecolor="black", linewidth=0.8, zorder=3)
    _label_bars(ax_kg, bars, kg_vals, kg_sales)

    ax_kg.set_ylabel("Minimum selling price ($/kg product)")
    ax_kg.set_title("Biostimulant, AD-VFA, AD-fermentation", fontsize=10)
    ax_kg.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.2f"))
    ax_kg.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax_kg.axhline(0, color="black", linewidth=0.8, zorder=3)
    kg_lo, kg_hi = min(0, min(kg_vals)), max(0, max(kg_vals))
    pad = (kg_hi - kg_lo) * 0.22 or 1
    ax_kg.set_ylim(kg_lo - pad, kg_hi + pad)

    # -- Right panel: AD-biomethane by pretreatment case, $/mmbtu ------------
    pt_tick_labels = [
        PRETREATMENT_LABELS[r["pretreatment_case"]] for r in biomethane_results
    ]
    mmbtu_vals = [r["msp_usd_per_mmbtu"] for r in biomethane_results]
    mmbtu_sales = [
        r["msp_usd_per_mmbtu"] * r["annual_product_mmbtu"] for r in biomethane_results
    ]

    bars = ax_biomethane.bar(
        pt_tick_labels, mmbtu_vals, edgecolor="black", linewidth=0.8, zorder=3,
        color="#C86E5A",
    )
    _label_bars(ax_biomethane, bars, mmbtu_vals, mmbtu_sales)

    ax_biomethane.set_ylabel("Minimum selling price ($/mmbtu biomethane)")
    ax_biomethane.set_title("AD-biomethane, by pretreatment case", fontsize=10)
    ax_biomethane.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.2f"))
    ax_biomethane.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax_biomethane.axhline(0, color="black", linewidth=0.8, zorder=3)
    mmbtu_lo, mmbtu_hi = min(0, min(mmbtu_vals)), max(0, max(mmbtu_vals))
    pad = (mmbtu_hi - mmbtu_lo) * 0.22 or 1
    ax_biomethane.set_ylim(mmbtu_lo - pad, mmbtu_hi + pad)
    ax_biomethane.tick_params(axis="x", labelrotation=15)

    fig.suptitle(
        "MSP by SaBRe pathway (near-zero feed price)\n"
        f"{scenario_desc}",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out_path = OUT / f"fig_msp_comparison_four_systems_{scenario_tag}.png"
    fig.savefig(out_path, bbox_inches="tight")

    print("\nSaved:")
    print(f"  {out_path}")

    # ── Mass-yield summary ──────────────────────────────────────────────────
    print(f"\nMass-yield summary ({scenario_desc}):")
    header = (
        f"  {'Pathway':<38}"
        f"{'product yield (wt% of dry feed)':>34}{'product purity':>18}"
    )
    print(header)
    for row in summary_rows:
        if "product_yield_dry_pct" in row:
            product_str = f"{row['product_yield_dry_pct']:.2f}%"
            if row.get("product_yield_is_wet"):
                product_str += " (wet)"
        else:
            product_str = "n/a"
        print(
            f"  {row['label']:<38}"
            f"{product_str:>34}"
            f"{row['product_purity_frac'] * 100:>17.2f}%"
        )


if __name__ == "__main__":
    run_msp_comparison(credit_tipping_fee=False)
    run_msp_comparison(credit_tipping_fee=True)

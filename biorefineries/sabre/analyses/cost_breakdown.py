# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
cost_breakdown.py
------------------
Generate a capital- and operating-cost breakdown for the four standalone
SaBRe flowsheets (the 'integrated' flowsheet is excluded, same as
`msp_comparison.py`):

    - Biostimulant
    - AD-biomethane   (data/pretreatment.yaml `pretreatment_ad` case
                        'combined_PE')
    - AD-VFA
    - AD-fermentation

Each pathway is built via its `price_*_system()` function (not `create_*_
system()` directly) with `credit_tipping_fee=False` -- the default economic
basis, also used by `msp_comparison.py`'s first pass -- where AD-biomethane/
AD-VFA/AD-fermentation each embed their own copy of the biostimulant
preprocessing subsystem (Press, PressateConcentrator, Evaporator, etc.) and
price biostimulant_product as a co-product at its flat tea.yaml baseline, on
one combined TEA. So FCI/AOC for those three pathways include a share of
biostimulant's own capital and OPEX -- this doesn't change with
credit_tipping_fee (that flag only changes which TEA scope biosteam.solve_price
solves against, not which units the `sys` returned here contains).

For each system this script builds a per-unit breakdown of:
    - FCI  (fixed capital investment): each unit's installed cost, scaled up
      by FCI / installed_equipment_cost so the per-unit values sum to FCI
      exactly (this scaling is uniform across units, so it preserves each
      unit's share of the total).
    - AOC  (annual operating cost = FOC + VOC): each unit's utility cost +
      add_OPEX, plus its allocated share of FOC (using the same FCI-share
      weighting as above, since FOC is itself defined as a fraction of FCI).
      Costs that aren't tied to any one unit -- e.g. feedstock/reagent
      material cost, which is a property of priced boundary streams, not
      of equipment -- get their own non-unit category instead of being
      forced onto a unit (currently $0 for most SaBRe systems, but this
      keeps the breakdown correct if that changes).

Outputs (relative to this script's directory):
    ../results/tables/cost_breakdown_<system>.csv      per-unit + category
                                                         breakdown table
    ../results/tables/cost_breakdown_summary.csv        cross-system summary
    ../results/figures/fig_cost_breakdown.png            FCI/AOC totals by
                                                         system
    ../results/figures/fig_cost_breakdown_<system>.png   FCI/AOC breakdown
                                                         by unit (small
                                                         contributors grouped
                                                         into "Others")

Run from the repo root:
    python biorefineries/sabre/analyses/cost_breakdown.py
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

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

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

MATERIAL_COST_LABEL = "Feedstock & reagents (material cost)"

# Small unit/category contributors get folded into "Others" in the by-unit
# pie charts so the legend doesn't get overwhelmed with slivers.
GROUP_THRESHOLD_FRAC = 0.04
GROUP_MAX_SLICES = 7


SYSTEM_BUILDERS = {
    "Biostimulant": price_biostimulant_system,
    "AD-biomethane": lambda: price_ad_biomethane_system(pretreatment_case="combined_PE"),
    "AD-VFA": price_ad_vfa_system,
    "AD-fermentation": price_ad_fermentation_system,
}


def _unit_annual_add_opex(unit, operating_days: float) -> float:
    """Annualize a unit's add_OPEX ($/hr dict or scalar), mirroring
    `SaBReTEA._annual_unit_add_opex` but for a single unit."""
    add_opex = getattr(unit, "add_OPEX", None)
    if not add_opex:
        return 0.0
    if isinstance(add_opex, dict):
        hourly = sum(float(v or 0.0) for v in add_opex.values())
    else:
        hourly = float(add_opex or 0.0)
    return hourly * operating_days * 24.0


def unit_cost_breakdown(tea) -> pd.DataFrame:
    """
    Per-unit + non-unit-category breakdown of FCI and AOC for one TEA.

    Columns sum exactly to `tea.FCI` and `tea.AOC`, respectively (modulo
    floating-point error), so the table can be used both for reporting and
    for the grouped-pie-chart breakdown.
    """
    units = sorted(
        [u for u in tea.system.units if u._design or u._cost], key=lambda x: x.line
    )
    operating_hours = tea.operating_days * 24.0
    installed_equipment_cost = tea.installed_equipment_cost
    fci_scale = tea.FCI / installed_equipment_cost if installed_equipment_cost else 0.0

    rows = []
    for u in units:
        installed_cost = u.installed_cost
        utility_cost_annual = u.utility_cost * operating_hours
        add_opex_annual = _unit_annual_add_opex(u, tea.operating_days)
        fci_allocated = installed_cost * fci_scale
        foc_allocated = fci_allocated * tea.foc_frac_of_fci
        rows.append({
            "Unit": u.ID,
            "Unit operation": u.line,
            "Purchase cost (USD)": u.purchase_cost,
            "Installed cost (USD)": installed_cost,
            "Utility cost (USD/yr)": utility_cost_annual,
            "Other add_OPEX (USD/yr)": add_opex_annual,
            "FCI allocated (USD)": fci_allocated,
            "FOC allocated (USD/yr)": foc_allocated,
            "AOC allocated (USD/yr)": utility_cost_annual + add_opex_annual + foc_allocated,
        })

    # Non-unit category: costs tied to priced boundary streams (feedstock,
    # reagents), not to any one piece of equipment.
    if tea.material_cost:
        rows.append({
            "Unit": MATERIAL_COST_LABEL,
            "Unit operation": MATERIAL_COST_LABEL,
            "Purchase cost (USD)": 0.0,
            "Installed cost (USD)": 0.0,
            "Utility cost (USD/yr)": 0.0,
            "Other add_OPEX (USD/yr)": 0.0,
            "FCI allocated (USD)": 0.0,
            "FOC allocated (USD/yr)": 0.0,
            "AOC allocated (USD/yr)": tea.material_cost,
        })

    return pd.DataFrame(rows)


def cost_breakdown_for_system(label: str, builder) -> tuple[pd.DataFrame, dict]:
    """
    Run one SaBRe pathway's `price_*_system()` function (which builds,
    simulates, and solves its product's MSP internally) and return its
    per-unit/category cost table and a system-level summary dict.
    """
    result = builder()
    sys = result["sys"]
    tea = sys.TEA

    unit_table = unit_cost_breakdown(tea)
    unit_table.insert(0, "System", label)

    summary = {
        "system": label,
        "installed_equipment_cost_usd": tea.installed_equipment_cost,
        "purchase_cost_usd": tea.purchase_cost,
        "FCI_usd": tea.FCI,
        "TCI_usd": tea.TCI,
        "material_cost_usd_per_yr": tea.material_cost,
        "utility_cost_usd_per_yr": tea.utility_cost,
        "other_opex_usd_per_yr": tea._annual_unit_add_opex(),
        "FOC_usd_per_yr": tea.FOC,
        "VOC_usd_per_yr": tea.VOC,
        "AOC_usd_per_yr": tea.AOC,
        "sales_usd_per_yr": tea.sales,
        "msp_usd_per_kg": result.get("msp_usd_per_kg"),
        "msp_usd_per_mmbtu": result.get("msp_usd_per_mmbtu"),
    }
    return unit_table, summary


def group_small_contributors(
    df: pd.DataFrame,
    value_col: str,
    label_col: str = "Unit",
    threshold_frac: float = GROUP_THRESHOLD_FRAC,
    max_slices: int = GROUP_MAX_SLICES,
) -> pd.DataFrame:
    """
    Collapse small contributors (below `threshold_frac` of the total, or
    beyond the top `max_slices`) into a single "Others" row, so per-unit
    pie/bar breakdowns stay readable regardless of how many units a system
    has.
    """
    df = df.loc[df[value_col] > 0, [label_col, value_col]].copy()
    total = df[value_col].sum()
    if df.empty or total <= 0:
        return df

    df = df.sort_values(value_col, ascending=False).reset_index(drop=True)
    keep = df[value_col] >= threshold_frac * total
    if keep.sum() > max_slices:
        keep = pd.Series(False, index=df.index)
        keep.iloc[:max_slices] = True

    dropped = ~keep
    if dropped.sum() <= 1:
        # Nothing to gain from an "Others" bucket of size <=1 -- keep as is.
        return df
    grouped = pd.concat([
        df.loc[keep],
        pd.DataFrame([{label_col: "Others", value_col: df.loc[dropped, value_col].sum()}]),
    ], ignore_index=True)
    return grouped.sort_values(value_col, ascending=False).reset_index(drop=True)


def plot_breakdown_pie(ax, df: pd.DataFrame, value_col: str, label_col: str, title: str):
    if df.empty:
        ax.axis("off")
        ax.set_title(title, fontsize=10)
        return

    values = df[value_col].to_numpy()
    labels = df[label_col].to_numpy()
    wedges, _, _ = ax.pie(
        values,
        autopct=lambda p: f"{p:.0f}%" if p >= 4 else "",
        pctdistance=0.75,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(edgecolor="black", linewidth=0.6),
        textprops=dict(fontsize=8),
    )
    legend_labels = [f"{lab}  (${val / 1e6:,.2f}M)" for lab, val in zip(labels, values)]
    ax.legend(
        wedges, legend_labels, loc="center left", bbox_to_anchor=(1.02, 0.5),
        fontsize=7, frameon=False,
    )
    ax.set_title(title, fontsize=10)


# =================================================================
# Main
# =================================================================

if __name__ == "__main__":
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

    unit_tables = []
    summaries = []

    for label, builder in SYSTEM_BUILDERS.items():
        unit_table, summary = cost_breakdown_for_system(label, builder)
        unit_tables.append(unit_table)
        summaries.append(summary)

        csv_name = f"cost_breakdown_{label.lower().replace('-', '_')}.csv"
        unit_table.to_csv(TABLES_OUT / csv_name, index=False)

        # -- Per-system breakdown by unit (small contributors -> "Others") --
        fci_group = group_small_contributors(unit_table, "FCI allocated (USD)")
        aoc_group = group_small_contributors(unit_table, "AOC allocated (USD/yr)")

        fig, (ax_fci, ax_aoc) = plt.subplots(1, 2, figsize=(12.5, 5.2))
        plot_breakdown_pie(
            ax_fci, fci_group, "FCI allocated (USD)", "Unit",
            f"FCI by unit (${summary['FCI_usd'] / 1e6:,.1f}M total)",
        )
        plot_breakdown_pie(
            ax_aoc, aoc_group, "AOC allocated (USD/yr)", "Unit",
            f"AOC by unit (${summary['AOC_usd_per_yr'] / 1e6:,.1f}M/yr total)",
        )
        fig.suptitle(f"{label}: cost breakdown by unit", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        fig_name = f"fig_cost_breakdown_{label.lower().replace('-', '_')}.png"
        fig.savefig(FIGURES_OUT / fig_name, bbox_inches="tight")
        plt.close(fig)

    summary_df = pd.DataFrame(summaries).set_index("system")
    summary_df.to_csv(TABLES_OUT / "cost_breakdown_summary.csv")

    print("\nSaved per-unit cost tables and summary to:")
    print(f"  {TABLES_OUT}")

    # ── Cross-system overview: FCI and AOC totals ───────────────────────────
    labels = list(summary_df.index)
    fci_musd = summary_df["FCI_usd"] / 1e6
    opex_categories = [
        ("material_cost_usd_per_yr", "Material"),
        ("utility_cost_usd_per_yr", "Utility"),
        ("other_opex_usd_per_yr", "Other add_OPEX"),
        ("FOC_usd_per_yr", "Fixed OPEX"),
    ]

    fig, (ax_capex, ax_opex) = plt.subplots(1, 2, figsize=(11.5, 5.2))

    bars = ax_capex.bar(labels, fci_musd, edgecolor="black", linewidth=0.8, zorder=3)
    for bar, val in zip(bars, fci_musd):
        ax_capex.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + fci_musd.max() * 0.015,
            f"${val:,.2f}M",
            ha="center", va="bottom", fontsize=9,
        )
    ax_capex.set_ylabel("Fixed capital investment (10$^6$ USD)")
    ax_capex.set_title("FCI by system", fontsize=10)
    ax_capex.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax_capex.set_ylim(0, fci_musd.max() * 1.15)
    ax_capex.tick_params(axis="x", labelrotation=15)

    bottom = pd.Series(0.0, index=labels)
    for key, cat_label in opex_categories:
        vals = summary_df[key] / 1e6
        ax_opex.bar(
            labels, vals, bottom=bottom, edgecolor="black", linewidth=0.8,
            zorder=3, label=cat_label,
        )
        bottom = bottom + vals

    ax_opex.set_ylabel("Annual operating cost (10$^6$ USD/yr)")
    ax_opex.set_title("AOC breakdown by category, by system", fontsize=10)
    ax_opex.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax_opex.legend(fontsize=8, loc="upper right")
    ax_opex.tick_params(axis="x", labelrotation=15)

    fig.suptitle(
        "SaBRe cost breakdown by pathway\n"
        "See fig_cost_breakdown_<system>.png for the by-unit breakdown of each bar",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIGURES_OUT / "fig_cost_breakdown.png", bbox_inches="tight")
    plt.close(fig)

    print(f"  {FIGURES_OUT}")

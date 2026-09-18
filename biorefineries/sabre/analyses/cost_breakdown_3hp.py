# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
cost_breakdown_3hp.py
----------------------
Baseline capital- and operating-cost breakdown for the standalone 3hp
flowsheet (systems._3hp_system.create_3hp_system()), at the fixed
data/tea.yaml assumption basis (systems._3hp_system.price_3hp_system()'s
only mode -- no credit_tipping_fee variant, unlike the AD pathways').

Kept as its own script rather than folded into cost_breakdown.py's
SYSTEM_BUILDERS dict: the 3hp system needs
sabre._chemicals.create_chemicals(include_hp3=True)'s chemical superset,
which must not share a Python process with any of the other four systems
(see systems/_3hp_system.py's own module docstring) -- cost_breakdown.py
builds all four of those together in one process, so 3hp can't safely
join that loop. Reuses cost_breakdown.py's per-unit breakdown/grouping/
plotting helpers directly (pure functions, no side effects at import
time), rather than duplicating them.

Outputs (relative to this script's directory):
    ../results/tables/cost_breakdown_3hp.csv     per-unit + category breakdown
    ../results/tables/cost_breakdown_summary_3hp.csv   one-row system summary
    ../results/figures/fig_cost_breakdown_3hp.png       FCI/AOC breakdown by unit

Run from the repo root:
    python biorefineries/sabre/analyses/cost_breakdown_3hp.py
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre.systems._3hp_system import price_3hp_system
from biorefineries.sabre.analyses.cost_breakdown import (
    unit_cost_breakdown, group_small_contributors, plot_breakdown_pie,
)

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

LABEL = "3-HP"


def cost_breakdown_for_3hp() -> tuple[pd.DataFrame, dict]:
    """
    Run price_3hp_system() (builds, simulates, and solves the Ca(3HP)2
    product's MSP internally) and return its per-unit/category cost table
    and a system-level summary dict -- same shape as
    cost_breakdown.py::cost_breakdown_for_system(), so the two scripts'
    outputs stay directly comparable.
    """
    result = price_3hp_system()
    sys = result["sys"]
    tea = sys.TEA

    unit_table = unit_cost_breakdown(tea)
    unit_table.insert(0, "System", LABEL)

    summary = {
        "system": LABEL,
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
    }
    return unit_table, summary


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

    unit_table, summary = cost_breakdown_for_3hp()

    unit_table.to_csv(TABLES_OUT / "cost_breakdown_3hp.csv", index=False)
    pd.DataFrame([summary]).set_index("system").to_csv(
        TABLES_OUT / "cost_breakdown_summary_3hp.csv"
    )

    print("\nSaved per-unit cost table and summary to:")
    print(f"  {TABLES_OUT}")

    # -- Breakdown by unit (small contributors -> "Others") -----------------
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
    fig.suptitle(f"{LABEL}: cost breakdown by unit", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIGURES_OUT / "fig_cost_breakdown_3hp.png", bbox_inches="tight")
    plt.close(fig)

    print(f"  {FIGURES_OUT}")

    print(f"\n{LABEL} summary:")
    for k, v in summary.items():
        if k == "system":
            continue
        print(f"  {k:<32}{v:,.4f}" if isinstance(v, float) else f"  {k:<32}{v}")

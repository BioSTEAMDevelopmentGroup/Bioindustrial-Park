# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
ad_heatmap.py
-------------------
2D heatmap of biomethane MSP ($/MMBtu) as a function of
feedstock price (x-axis) and biostimulant price (y-axis).

Viable cells are highlighted with bold borders:
  - black border : MSP <= $14/MMBtu
  - blue border  : MSP <= $3/MMBtu

Outputs one heatmap per pretreatment case (main text: combined_PE only,
supplementary: all five cases in a panel).

Run from sabre_project root:
    python scripts/ad_heatmap.py
"""

import sys
from pathlib import Path

import biosteam as bst
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.systems import create_ad_biomethane_system
from biorefineries.sabre._tea import create_tea, solve_biomethane_msp


OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ── Grid definition ───────────────────────────────────────────────────────────

FEED_PRICES    = [-0.02, 0.00, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10]  # $/kg wet
BIOSTIM_PRICES = [0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00]  # $/kg

PRETREATMENT_CASES = [
    "press_mill_only",
    "enzymatic",
    "peroxide",
    "combined_PE",
    "combined_PTE",
]

CASE_LABELS = {
    "press_mill_only": "Press + mill only",
    "enzymatic":       "Enzymatic",
    "peroxide":        "Peroxide",
    "combined_PE":     "Combined PE",
    "combined_PTE":    "Combined PTE",
}

LIQUID_DIGESTATE_USD_PER_KG = -0.002
SOLIDS_DIGESTATE_USD_PER_KG = -0.02


# ── Helpers ───────────────────────────────────────────────────────────────────

def apply_stream_economics(sys, biostimulant_price: float):
    for sid, price in [
        ("pressate_permeate",              0.0),
        ("biostimulant_product",           biostimulant_price),
        ("soil_amendment",                 SOLIDS_DIGESTATE_USD_PER_KG),
        ("liquid_digestate",               LIQUID_DIGESTATE_USD_PER_KG),
    ]:
        try:
            sys.flowsheet.stream[sid].price = price
        except Exception:
            pass


def build_case(pretreatment_case: str, feed_price: float, biostim_price: float) -> float:
    bst.main_flowsheet.clear()
    create_chemicals()

    sys = create_ad_biomethane_system(
        feedstock="pelagic",
        pretreatment_case=pretreatment_case,
    )
    sys.feeds[0].price = feed_price
    sys.simulate()

    apply_stream_economics(sys, biostimulant_price=biostim_price)

    tea = create_tea(sys)
    msp = solve_biomethane_msp(tea, sys.flowsheet.stream.biomethane)
    return float(msp["usd_per_mmbtu"])


# ── Run grid ──────────────────────────────────────────────────────────────────
if __name__ == "__main__":

    # msp_grids[case] = 2-D array, shape (len(BIOSTIM_PRICES), len(FEED_PRICES))
    # rows = biostimulant price (y-axis), cols = feed price (x-axis)
    msp_grids: dict[str, np.ndarray] = {}

    for case in PRETREATMENT_CASES:
        print(f"\n=== {CASE_LABELS[case]} ===")
        grid = np.zeros((len(BIOSTIM_PRICES), len(FEED_PRICES)))

        for i, bs in enumerate(BIOSTIM_PRICES):
            for j, fp in enumerate(FEED_PRICES):
                msp = build_case(case, fp, bs)
                grid[i, j] = msp
                print(f"  feed={fp:+.2f}  bs={bs:.2f}  -> MSP=${msp:.2f}")

        msp_grids[case] = grid


    # ── Shared colormap setup ─────────────────────────────────────────────────────

    plt.rcParams.update({
        "font.family":     "DejaVu Sans",
        "font.size":       10,
        "axes.titlesize":  11,
        "axes.labelsize":  10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "figure.dpi":      150,
        "axes.linewidth":  0.8,
    })

    # Colormap: green = cheap/viable, red = expensive
    CMAP = "RdYlGn_r"

    # Global vmin/vmax across all cases for panel comparability
    all_vals = np.concatenate([g.ravel() for g in msp_grids.values()])
    VMIN, VMAX = np.percentile(all_vals, 2), np.percentile(all_vals, 98)


    def draw_heatmap(ax, grid, case, show_xlabel=True, show_ylabel=True, show_cbar=False, fig=None):
        """Draw one heatmap panel onto ax."""

        n_rows, n_cols = grid.shape

        # Plot on integer indices so every cell is exactly the same pixel size
        im = ax.imshow(
            grid,
            aspect="auto",
            origin="lower",
            cmap=CMAP,
            vmin=VMIN,
            vmax=VMAX,
            extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5],
        )

        # Cell annotations
        for i in range(n_rows):
            for j in range(n_cols):
                val = grid[i, j]
                normed = (val - VMIN) / max(VMAX - VMIN, 1)
                color = "white" if normed < 0.35 or normed > 0.72 else "black"
                ax.text(
                    j, i, f"${val:.0f}",
                    ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold",
                )

        # Tick positions = integer indices; labels = actual price values
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels([f"${fp:.2f}" for fp in FEED_PRICES], fontsize=8)
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels([f"${bs:.2f}" for bs in BIOSTIM_PRICES], fontsize=8)

        if show_xlabel:
            ax.set_xlabel("Feedstock price ($/kg wet)")
        if show_ylabel:
            ax.set_ylabel("Biostimulant price ($/kg)")

        ax.set_title(CASE_LABELS[case], fontsize=10, pad=4)

        if show_cbar and fig is not None:
            cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("MSP ($/MMBtu)", fontsize=9)

        return im


    # ── Figure 1: Combined PE only (main text) ────────────────────────────────────

    fig1, ax1 = plt.subplots(figsize=(6.5, 4.5))
    im = draw_heatmap(ax1, msp_grids["combined_PE"], "combined_PE", fig=fig1, show_cbar=False)

    cbar = fig1.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label("Biomethane MSP ($/MMBtu)", fontsize=9)

    fig1.tight_layout()
    fig1.savefig(OUT / "fig_msp_heatmap_combined_PE.png", bbox_inches="tight")
    print("\nSaved: fig_msp_heatmap_combined_PE.png")


    # ── Figure 2: All five cases panel (supplementary) ────────────────────────────
    # 3 rows × 2 cols; last slot used for shared colorbar

    THRESHOLD_TTF = 14.0   # $/MMBtu — black border
    THRESHOLD_HH  =  3.0   # $/MMBtu — blue border

    fig2, axes = plt.subplots(3, 2, figsize=(11, 14), constrained_layout=False)
    fig2.subplots_adjust(hspace=0.20, wspace=0.25, bottom=0.10)
    axes_flat = axes.flatten()

    last_im = None
    for idx, case in enumerate(PRETREATMENT_CASES):
        ax = axes_flat[idx]
        show_x = idx >= 3          # bottom row only
        show_y = idx % 2 == 0      # left column only
        im = draw_heatmap(
            ax, msp_grids[case], case,
            show_xlabel=show_x,
            show_ylabel=show_y,
            fig=fig2,
        )
        last_im = im

        # Cell border highlighting
        grid = msp_grids[case]
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                val = grid[i, j]
                if val <= THRESHOLD_HH:
                    edgecolor = "#1a6faf"
                elif val <= THRESHOLD_TTF:
                    edgecolor = "black"
                else:
                    continue
                rect = plt.Rectangle(
                    (j - 0.5, i - 0.5), 1.0, 1.0,
                    linewidth=2.0, edgecolor=edgecolor,
                    facecolor="none", zorder=6,
                )
                ax.add_patch(rect)

    # Hide the unused 6th panel
    axes_flat[-1].set_visible(False)

    # Shared colorbar in the empty 6th slot
    cbar_ax = axes_flat[-1].get_position()
    cb_ax = fig2.add_axes([cbar_ax.x0, cbar_ax.y0 + cbar_ax.height * 0.15,
                            0.025, cbar_ax.height * 0.70])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(vmin=VMIN, vmax=VMAX))
    sm.set_array([])
    cbar = fig2.colorbar(sm, cax=cb_ax)
    cbar.set_label("Biomethane MSP ($/MMBtu)", fontsize=10)

    # Bottom legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="none", edgecolor="black",   linewidth=2.0,
              label="MSP ≤ $14/MMBtu"),
        Patch(facecolor="none", edgecolor="#1a6faf", linewidth=2.0,
              label="MSP ≤ $3/MMBtu"),
    ]
    fig2.legend(
        handles=legend_handles,
        fontsize=9, frameon=False, ncol=2,
        loc="lower center",
        bbox_to_anchor=(0.45, 0.01),
    )

    fig2.suptitle(
        "Biomethane MSP ($/MMBtu) — feedstock price × biostimulant price sensitivity",
        fontsize=12, y=1.01,
    )

    fig2.savefig(OUT / "fig_msp_heatmap_all_cases.png", bbox_inches="tight")
    print("Saved: fig_msp_heatmap_all_cases.png")


    # ── Console summary ───────────────────────────────────────────────────────────

    print("\nMSP grid — Combined PE ($/MMBtu):")
    print(f"{'':>12}", end="")
    for fp in FEED_PRICES:
        print(f"  feed={fp:+.2f}", end="")
    print()
    for i, bs in enumerate(BIOSTIM_PRICES):
        print(f"bs={bs:.2f}/kg  ", end="")
        for j in range(len(FEED_PRICES)):
            print(f"  {msp_grids['combined_PE'][i, j]:9.2f}", end="")
        print()
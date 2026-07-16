# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Integrated-system TEA figures
===========================================
Figures:
  1) Alpha sweep with market price uncertainty bands — four scenarios
  2) NPV bar chart at alpha=0 — four scenarios × three biomethane prices
  3) 2D heatmap — biomethane × oil price → NPV at alpha=0, base scenario
  4) 2×3 panel heatmaps — alpha=0/0.25/0.50/0.75/1.0, three feed prices
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import biosteam as bst
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parents[2]
for _p in (REPO_ROOT, THIS_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from integrated_tea import (
    BIOMETHANE_MARKET_MMBTU,
    BIOMETHANE_PRICES,
    OIL_MARKET_USD_PER_KG,
    _apply_stream_prices,
    _compute_npv_at_market,
    _patch_ev607,
    _wire_oil_reagent,
    run_alpha_sweep,
)
from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.systems import create_integrated_biorefinery
from biorefineries.sabre._tea import create_tea

# ── Output directories ────────────────────────────────────────────────────────
OUTPUT_DIR = THIS_DIR.parent / "results" / "integrated_figures"
CSV_DIR    = OUTPUT_DIR / "csv"
FIG_DIR    = OUTPUT_DIR / "figures"

# ── Scenario definitions ──────────────────────────────────────────────────────
@dataclass(frozen=True)
class Scenario:
    name:               str
    label:              str
    short_label:        str
    feed_price:         float
    biostimulant_price: float
    pretreatment_case:  str   = "press_mill_only"
    market_mmbtu:       float = BIOMETHANE_MARKET_MMBTU
    market_oil:         float = OIL_MARKET_USD_PER_KG


SCENARIOS: tuple[Scenario, ...] = (
    Scenario(
        name="base",
        label="Base | feed $0.02/kg | biostim $0.50/kg",
        short_label="Base\n($0.02/kg feed)",
        feed_price=0.02,
        biostimulant_price=0.50,
    ),
    Scenario(
        name="tipping_fee",
        label="Tipping fee | feed -$0.02/kg | biostim $0.50/kg",
        short_label="Tipping fee\n(-$0.02/kg)",
        feed_price=-0.02,
        biostimulant_price=0.50,
    ),
    Scenario(
        name="biostim",
        label="Base | feed \$0.02/kg | biostim \$3.00/kg",
        short_label="Biostimulant\n(\$3.00/kg)",
        feed_price=0.02,
        biostimulant_price=3.00,
    ),
    Scenario(
        name="best_case",
        label="Best case | tipping fee | biostim \$3.00/kg",
        short_label="Best case\n(tip + biostim)",
        feed_price=-0.02,
        biostimulant_price=3.00,
    ),
)

BASE_SCENARIO_NAME = "base"

SCENARIO_COLORS = {
    "base":        "#1f77b4",
    "tipping_fee": "#ff7f0e",
    "biostim":     "#2ca02c",
    "best_case":   "#d62728",
}

BIOMETHANE_PRICE_LABELS = {
    3.0:  "$3/MMBtu (Henry Hub)",
    10.0: "$10/MMBtu (TTF)",
    14.0: "$14/MMBtu (JKM)",
}

TEXT = "#2C2C2A"


# ── Helpers ───────────────────────────────────────────────────────────────────
def _save_dataframe(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"Saved CSV: {path}")


def _build_npv_at_alpha(
    alpha: float,
    feed_price: float,
    biostimulant_price: float,
    pretreatment_case: str,
    market_mmbtu: float,
    market_oil: float,
) -> float:
    """Run one simulation at given alpha and return NPV ($M)."""
    bst.main_flowsheet.clear()
    set_thermo()

    sys, streams, units, _ = create_integrated_biorefinery(
        alpha=alpha,
        pretreatment_case=pretreatment_case,
    )
    streams["feed"].price = feed_price
    sys.simulate()

    _patch_ev607()
    _apply_stream_prices(streams, biostimulant_price)
    _wire_oil_reagent(streams, units)

    tea = create_tea(sys)
    npv = _compute_npv_at_market(tea, streams, market_mmbtu, market_oil)
    return npv / 1e6


def _style_ax(ax):
    ax.grid(False)
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("#1a1a1a")
    ax.tick_params(labelsize=9)


def _add_npv0_contour(ax, grid, n_rows, n_cols):
    """Add NPV=0 dashed contour to a heatmap panel."""
    try:
        from scipy.interpolate import RegularGridInterpolator
        bm_f    = np.linspace(0, n_cols - 1, 300)
        oil_f   = np.linspace(0, n_rows - 1, 300)
        BM, OIL = np.meshgrid(bm_f, oil_f)
        interp  = RegularGridInterpolator(
            (np.arange(n_rows), np.arange(n_cols)), grid,
            method="linear", bounds_error=False)
        Z  = interp(np.stack([OIL.ravel(), BM.ravel()], axis=1)).reshape(BM.shape)
        cs = ax.contour(BM, OIL, Z, levels=[0], colors=["black"],
                        linewidths=1.8, linestyles=["--"], zorder=5)
        ax.clabel(cs, fmt={0: "NPV=0"}, fontsize=6.5, inline=True)
    except Exception:
        pass


# ── Figure 1 — Alpha sweep with uncertainty bands ─────────────────────────────
def plot_alpha_sweep_with_bands(outpath: Path) -> None:
    """
    Four scenarios, each with:
      - Central line : biomethane $10/MMBtu · oil $1.00/kg
      - Shaded band  : [bm=$3, oil=$0.62] to [bm=$14, oil=$1.50]

    Scenarios:
      1. feed=$0.02,  biostim=$0.50
      2. feed=-$0.02, biostim=$0.50
      3. feed=$0.02,  biostim=$3.00
      4. feed=-$0.02, biostim=$3.00
    """
    BAND_SCENARIOS = [
        {"feed":  0.02, "biostim": 0.50,
        "label": r"Feed \$0.02/kg $\cdot$ Biostim \$0.50/kg",
        "color": "#1f77b4", "ls": "-",  "marker": "o"},
        {"feed": -0.02, "biostim": 0.50,
        "label": r"Tipping fee $\cdot$ Biostim \$0.50/kg",
        "color": "#ff7f0e", "ls": "--", "marker": "^"},
        {"feed":  0.02, "biostim": 3.00,
        "label": r"Feed \$0.02/kg $\cdot$ Biostim \$3.00/kg",
        "color": "#2ca02c", "ls": "-",  "marker": "s"},
        {"feed": -0.02, "biostim": 3.00,
        "label": r"Tipping fee $\cdot$ Biostim \$3.00/kg",
        "color": "#d62728", "ls": "--", "marker": "D"},
    ]

    BM_LOW,  BM_MID,  BM_HIGH  =  3.0, 10.0, 14.0
    OIL_LOW, OIL_MID, OIL_HIGH = 0.62,  1.00,  1.50

    fig, ax = plt.subplots(figsize=(9.0, 5.5))
    fig.subplots_adjust(left=0.11, right=0.68, top=0.88, bottom=0.22)

    for sc in BAND_SCENARIOS:
        print(f"\nRunning band sweep: {sc['label']}")

        r_mid  = run_alpha_sweep(sc["feed"], "band", "press_mill_only",
                                 sc["biostim"], BM_MID,  OIL_MID,  False)
        r_low  = run_alpha_sweep(sc["feed"], "band", "press_mill_only",
                                 sc["biostim"], BM_LOW,  OIL_LOW,  False)
        r_high = run_alpha_sweep(sc["feed"], "band", "press_mill_only",
                                 sc["biostim"], BM_HIGH, OIL_HIGH, False)

        alphas   = [r["alpha"]          for r in r_mid  if r["ok"]]
        npv_mid  = [r["combined_npv_M"] for r in r_mid  if r["ok"]]
        npv_low  = [r["combined_npv_M"] for r in r_low  if r["ok"]]
        npv_high = [r["combined_npv_M"] for r in r_high if r["ok"]]

        ax.fill_between(alphas, npv_low, npv_high,
                        color=sc["color"], alpha=0.15, zorder=2)
        ax.plot(alphas, npv_mid,
                color=sc["color"], ls=sc["ls"], lw=1.8,
                marker=sc["marker"], markersize=5,
                markerfacecolor="white", markeredgewidth=1.4,
                markeredgecolor=sc["color"],
                label=sc["label"], zorder=3)

    ax.axhline(0, color="#555555", lw=0.9, ls="-", zorder=1)
    ax.text(1.02, 0, "NPV = 0", transform=ax.get_yaxis_transform(),
            va="center", fontsize=7.5, color="#555555")

    ax.set_xlabel("α  (fraction of milled biomass to biomethane pathway)",
                  fontsize=10, color=TEXT)
    ax.set_ylabel("Combined NPV  ($M)", fontsize=10, color=TEXT)
    ax.set_title(
        r"Integrated biorefinery NPV vs. pathway split ($\alpha$)" + "\n"
        r"Central line: \$10/MMBtu, \$1.00/kg oil  |  "
        r"Band: [\$3$-$\$14/MMBtu] $\times$ [\$0.62$-$\$1.50/kg oil]",
        fontsize=10, color=TEXT, pad=8,
    )
    ax.set_xlim(-0.03, 1.03)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"${v:,.0f}M"))
    _style_ax(ax)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24),
              ncol=2, frameon=True, fontsize=8.5)

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Figure 2 — NPV bar chart at alpha=1 ──────────────────────────────────────
def plot_npv_bar_chart(npv_table: dict, outpath: Path) -> None:
    sc_names  = [sc.name        for sc in SCENARIOS]
    sc_labels = [sc.short_label for sc in SCENARIOS]

    width   = 0.22
    x       = np.arange(len(sc_names))
    HATCHES = ["", "//", "xx"]
    COLORS  = ["#aec7e8", "#1f77b4", "#17517d"]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for k, (bm_price, hatch, color) in enumerate(
            zip(BIOMETHANE_PRICES, HATCHES, COLORS)):
        offsets  = x + (k - 1) * width
        npv_vals = [npv_table[sc][bm_price] for sc in sc_names]
        ax.bar(offsets, npv_vals, width=width, color=color,
               edgecolor="white", hatch=hatch, linewidth=0.5, zorder=3,
               label=BIOMETHANE_PRICE_LABELS[bm_price])

    ax.axhline(0, color="#555555", lw=0.9, zorder=2)
    ax.set_ylabel("Combined NPV at alpha = 1  (\$M)", fontsize=10, color=TEXT)
    ax.set_title(
        f"Integrated biorefinery NPV at alpha = 1\n"
        f"Oil \${OIL_MARKET_USD_PER_KG:.2f}/kg  ·  Press-mill  ·  Biostimulant \$0.50/kg",
        fontsize=10.5, color=TEXT, pad=8,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(sc_labels, fontsize=9)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"\${v:,.0f}M"))
    _style_ax(ax)
    ax.legend(fontsize=8.5, frameon=True, loc="upper left")

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Figure 3 — Single heatmap at alpha=0 ─────────────────────────────────────
def plot_npv_heatmap(outpath: Path) -> None:
    BM_PRICES_GRID  = [3.0, 5.0, 8.0, 10.0, 12.0, 14.0]
    OIL_PRICES_GRID = [0.62, 0.80, 1.00, 1.20, 1.50]

    base = next(s for s in SCENARIOS if s.name == BASE_SCENARIO_NAME)

    print("\nBuilding NPV heatmap (alpha=0, base scenario)...")
    grid = np.zeros((len(OIL_PRICES_GRID), len(BM_PRICES_GRID)))

    for i, oil_p in enumerate(OIL_PRICES_GRID):
        for j, bm_p in enumerate(BM_PRICES_GRID):
            npv = _build_npv_at_alpha(
                0.0, base.feed_price, base.biostimulant_price,
                base.pretreatment_case, bm_p, oil_p)
            grid[i, j] = npv
            print(f"  bm=${bm_p:.0f}  oil=${oil_p:.2f}  NPV=${npv:.1f}M")

    n_rows, n_cols = grid.shape
    vmax = max(abs(grid.min()), abs(grid.max()))
    vmin = -vmax

    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    im = ax.imshow(grid, aspect="auto", origin="lower", cmap="RdYlGn",
                   vmin=vmin, vmax=vmax,
                   extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5])

    for i in range(n_rows):
        for j in range(n_cols):
            val    = grid[i, j]
            normed = (val - vmin) / max(vmax - vmin, 1)
            color  = "white" if normed < 0.25 or normed > 0.75 else "black"
            ax.text(j, i, f"${val:.0f}M", ha="center", va="center",
                    fontsize=7.5, color=color)

    _add_npv0_contour(ax, grid, n_rows, n_cols)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([f"${p:.0f}" for p in BM_PRICES_GRID], fontsize=8.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f"${p:.2f}" for p in OIL_PRICES_GRID], fontsize=8.5)
    ax.set_xlabel("Biomethane price ($/MMBtu)", fontsize=10, color=TEXT)
    ax.set_ylabel("Microbial oil price ($/kg)", fontsize=10, color=TEXT)
    ax.set_title(
        "Integrated biorefinery NPV ($M) at alpha = 0\n"
        "Base scenario: feed $0.02/kg · biostimulant $0.50/kg · press-mill",
        fontsize=10.5, color=TEXT, pad=8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("NPV ($M)", fontsize=9)

    bm_ttf_idx  = BM_PRICES_GRID.index(10.0)
    bm_jkm_idx  = BM_PRICES_GRID.index(14.0)
    oil_mid_idx = OIL_PRICES_GRID.index(1.00)
    ax.axvline(bm_ttf_idx,  color="#1f77b4", lw=1.2, ls=":", zorder=4,
               label="$10/MMBtu (TTF)")
    ax.axvline(bm_jkm_idx,  color="#ff7f0e", lw=1.2, ls=":", zorder=4,
               label="$14/MMBtu (JKM)")
    ax.axhline(oil_mid_idx, color="#2ca02c", lw=1.2, ls=":", zorder=4,
               label="$1.00/kg oil")
    ax.legend(fontsize=7.5, frameon=True, loc="upper left")

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Figure 4 — Panel heatmaps varying alpha, one file per feed price ──────────
def plot_npv_heatmap_panel(
    outpath: Path,
    feed_price: float,
    feed_label: str,
    biostimulant_price: float = 0.50,
) -> None:
    """2x3 panel of NPV heatmaps at alpha=0/0.25/0.50/0.75/1.0."""
    BM_PRICES_GRID  = [3.0, 5.0, 8.0, 10.0, 12.0, 14.0]
    OIL_PRICES_GRID = [0.62, 0.80, 1.00, 1.20, 1.50]
    ALPHA_VALUES    = [0.0, 0.25, 0.50, 0.75, 1.0]

    grids: dict[float, np.ndarray] = {}

    for alpha in ALPHA_VALUES:
        print(f"\nHeatmap panel [{feed_label}] alpha={alpha:.2f}...")
        grid = np.zeros((len(OIL_PRICES_GRID), len(BM_PRICES_GRID)))
        for i, oil_p in enumerate(OIL_PRICES_GRID):
            for j, bm_p in enumerate(BM_PRICES_GRID):
                npv = _build_npv_at_alpha(
                    alpha, feed_price, biostimulant_price,
                    "press_mill_only", bm_p, oil_p)
                grid[i, j] = npv
                print(f"  a={alpha:.2f} bm=${bm_p:.0f} oil=${oil_p:.2f} NPV=${npv:.1f}M")
        grids[alpha] = grid

    all_vals = np.concatenate([g.ravel() for g in grids.values()])
    vmax = np.percentile(np.abs(all_vals), 98)
    vmin = -vmax
    n_rows, n_cols = len(OIL_PRICES_GRID), len(BM_PRICES_GRID)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=False)
    fig.subplots_adjust(hspace=0.32, wspace=0.25, bottom=0.08)
    axes_flat = axes.flatten()

    for idx, alpha in enumerate(ALPHA_VALUES):
        ax   = axes_flat[idx]
        grid = grids[alpha]

        ax.imshow(grid, aspect="auto", origin="lower", cmap="RdYlGn",
                  vmin=vmin, vmax=vmax,
                  extent=[-0.5, n_cols - 0.5, -0.5, n_rows - 0.5])

        for i in range(n_rows):
            for j in range(n_cols):
                val    = grid[i, j]
                normed = (val - vmin) / max(vmax - vmin, 1)
                color  = "white" if normed < 0.25 or normed > 0.75 else "black"
                ax.text(j, i, f"${val:.0f}M", ha="center", va="center",
                        fontsize=6.5, color=color)

        _add_npv0_contour(ax, grid, n_rows, n_cols)

        ax.set_xticks(range(n_cols))
        ax.set_yticks(range(n_rows))
        if idx >= 3:
            ax.set_xticklabels([f"${p:.0f}" for p in BM_PRICES_GRID], fontsize=7.5)
            ax.set_xlabel("Biomethane price ($/MMBtu)", fontsize=8.5, color=TEXT)
        else:
            ax.set_xticklabels([])
        if idx % 3 == 0:
            ax.set_yticklabels([f"${p:.2f}" for p in OIL_PRICES_GRID], fontsize=7.5)
            ax.set_ylabel("Microbial oil price ($/kg)", fontsize=8.5, color=TEXT)
        else:
            ax.set_yticklabels([])
        ax.set_title(f"alpha = {alpha:.2f}", fontsize=10, color=TEXT, pad=4)

    axes_flat[-1].set_visible(False)
    cbar_ax = fig.add_axes([0.69, 0.10, 0.02, 0.35])
    sm = plt.cm.ScalarMappable(cmap="RdYlGn",
                                norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label("NPV ($M)", fontsize=9)

    fig.suptitle(
        f"Integrated biorefinery NPV ($M) — biomethane x oil price sensitivity\n"
        f"Feed: ${feed_price:+.2f}/kg ({feed_label}) · "
        f"biostimulant ${biostimulant_price:.2f}/kg · press-mill",
        fontsize=11, color=TEXT, y=1.01,
    )

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {outpath}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Figure 1 — Alpha sweep with uncertainty bands
    plot_alpha_sweep_with_bands(FIG_DIR / "figure1_alpha_sweep_bands.png")

    # Figure 2 — NPV bar chart at alpha=0
    print("\nBuilding NPV bar chart data...")
    npv_table: dict[str, dict[float, float]] = {sc.name: {} for sc in SCENARIOS}
    for sc in SCENARIOS:
        for bm_price in BIOMETHANE_PRICES:
            npv = _build_npv_at_alpha(
                1.0, sc.feed_price, sc.biostimulant_price,
                sc.pretreatment_case, bm_price, OIL_MARKET_USD_PER_KG)
            npv_table[sc.name][bm_price] = npv
            print(f"  {sc.name}  bm=${bm_price:.0f}/MMBtu  NPV=${npv:.1f}M")
    plot_npv_bar_chart(npv_table, FIG_DIR / "figure2_npv_bar_chart.png")

    # Figure 3 — Single heatmap at alpha=0
    plot_npv_heatmap(FIG_DIR / "figure3_npv_heatmap.png")

    # Figure 4 — Panel heatmaps at five alpha values, three feed prices
    for feed_price, feed_label in [
        (-0.02, "tipping_fee"),
        ( 0.00, "near_zero"),
        ( 0.02, "base"),
    ]:
        plot_npv_heatmap_panel(
            outpath=FIG_DIR / f"figure4_alpha_panel_{feed_label}.png",
            feed_price=feed_price,
            feed_label=feed_label,
        )

    print(f"\nDone.\nFigures: {FIG_DIR.resolve()}\nCSVs:    {CSV_DIR.resolve()}")


if __name__ == "__main__":
    main()
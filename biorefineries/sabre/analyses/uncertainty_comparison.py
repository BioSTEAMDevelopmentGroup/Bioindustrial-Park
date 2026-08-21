# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
uncertainty_comparison.py
--------------------------
Compare MSP uncertainty (Monte Carlo, N=1000) and the top Spearman-correlated
parameters across the four standalone SaBRe flowsheets, using the
biorefineries/sabre/results/tables/model_*_N1000.xlsx files produced by
analyses/model_biostimulant.py, model_ad_vfa.py, model_ad_biomethane.py, and
model_ad_fermentation.py. The 'integrated' flowsheet is not covered (no
model file exists for it, per prior scope decision).

Each system's 'MSP' metric is USD/kg of that system's own value-carrying
product (biostimulant liquid, VFA broth, biomethane, crude microbial oil,
respectively) -- same $/kg-of-product basis, but not the same product, same
caveat as analyses/msp_comparison.py's single-point comparison.

Run from the repo root:
    python biorefineries/sabre/analyses/uncertainty_comparison.py
"""
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
TABLES_DIR = SCRIPT_DIR.parent / "results" / "tables"
OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

SYSTEMS = [
    ("Biostimulant", "model_biostimulant_N1000.xlsx"),
    ("AD-VFA", "model_ad_vfa_N1000.xlsx"),
    ("AD-biomethane", "model_ad_biomethane_N1000.xlsx"),
    ("AD-fermentation", "model_ad_fermentation_N1000.xlsx"),
]

# Fixed-order categorical palette, dataviz skill reference palette (validated
# CVD-safe in this order): blue, orange, aqua, yellow.
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
# Diverging pair (positive/negative Spearman rho), same reference palette.
DIVERGING_POS = "#e34948"  # red
DIVERGING_NEG = "#2a78d6"  # blue
GRIDLINE = "#e1e0d9"

PLT_RCPARAMS = {
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
}


def load_data(filename: str) -> pd.DataFrame:
    return pd.read_excel(TABLES_DIR / filename, sheet_name="Data", header=[0, 1], index_col=0)


def get_msp_series(df: pd.DataFrame) -> pd.Series:
    # Metric columns are tagged Element == '-'; ad_biomethane also has an
    # 'MSP (energy basis) [USD/mmbtu]' metric, excluded by requiring the
    # feature name to start with 'MSP [' (the $/kg one).
    col = next(c for c in df.columns if c[0] == '-' and c[1].startswith('MSP ['))
    return df[col].astype(float)


def get_parameter_columns(df: pd.DataFrame) -> list:
    return [c for c in df.columns if c[0] != '-']


def plot_msp_uncertainty(all_msp: dict) -> None:
    labels = [label for label, _ in SYSTEMS]
    data = [all_msp[label].values for label in labels]

    fig, ax = plt.subplots(figsize=(8, 5.5))

    bp = ax.boxplot(
        data, whis=(5, 95), showfliers=True,
        patch_artist=True, widths=0.55,
        flierprops=dict(
            marker='o', markersize=2.5, markerfacecolor='#898781',
            markeredgecolor='none', alpha=0.4,
        ),
        medianprops=dict(color='black', linewidth=1.3),
        whiskerprops=dict(color='black', linewidth=0.8),
        capprops=dict(color='black', linewidth=0.8),
        boxprops=dict(linewidth=0.8, edgecolor='black'),
    )
    for patch, color in zip(bp['boxes'], CATEGORICAL):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_yscale('symlog', linthresh=1)
    ax.set_ylabel("MSP (USD/kg product)")
    ax.axhline(0, color='black', linewidth=0.8, zorder=1)
    ax.grid(axis='y', linewidth=0.4, color=GRIDLINE, zorder=0)
    ax.set_title(
        "MSP uncertainty across SaBRe pathways (N=1000 Monte Carlo, ±50% broad sweep)\n"
        "box = IQR, whiskers = 5th–95th pct, dots = points beyond; symlog y-axis (linear within ±1)",
        fontsize=9,
    )

    out_path = OUT / "fig_msp_uncertainty_four_systems.png"
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


def plot_sensitivity(all_data: dict, top_n: int = 6) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.ravel()

    for ax, (label, _) in zip(axes, SYSTEMS):
        df = all_data[label]
        msp = get_msp_series(df)
        param_cols = get_parameter_columns(df)

        rows = []
        for col in param_cols:
            x = df[col].astype(float).values
            rho, p = spearmanr(x, msp.values)
            name = col[1].split(' [')[0]  # strip trailing " [units]"
            rows.append((name, rho, p))
        rows.sort(key=lambda r: abs(r[1]), reverse=True)
        top = rows[:top_n][::-1]  # reversed so the largest |rho| plots at the top

        names = [r[0] for r in top]
        rhos = [r[1] for r in top]
        colors = [DIVERGING_POS if r >= 0 else DIVERGING_NEG for r in rhos]

        ax.barh(names, rhos, color=colors, edgecolor='black', linewidth=0.6, zorder=3)
        ax.axvline(0, color='black', linewidth=0.8, zorder=3)
        ax.set_xlim(-1, 1)
        ax.set_title(label, fontsize=10)
        ax.grid(axis='x', linewidth=0.4, color=GRIDLINE, zorder=0)
        ax.set_xlabel("Spearman ρ vs MSP")
        ax.tick_params(axis='y', labelsize=8)

    fig.suptitle(
        f"Top {top_n} parameters by |Spearman ρ| vs MSP, per pathway (N=1000)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_path = OUT / "fig_sensitivity_comparison_four_systems.png"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


def main():
    plt.rcParams.update(PLT_RCPARAMS)

    all_data = {label: load_data(fn) for label, fn in SYSTEMS}
    all_msp = {label: get_msp_series(df) for label, df in all_data.items()}

    print("\nMSP summary (N=1000, USD/kg product):")
    for label, _ in SYSTEMS:
        s = all_msp[label]
        print(
            f"  {label:<18} median={s.median():>10.2f}  "
            f"IQR=[{s.quantile(.25):>9.2f}, {s.quantile(.75):>9.2f}]  "
            f"5-95%=[{s.quantile(.05):>10.2f}, {s.quantile(.95):>10.2f}]"
        )

    plot_msp_uncertainty(all_msp)
    plot_sensitivity(all_data)


if __name__ == "__main__":
    main()

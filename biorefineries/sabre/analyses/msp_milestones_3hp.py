# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
msp_milestones_3hp.py
---------------------
Ca(3HP)2 minimum selling price (MSP) if the preprocessing, hydrolysis and
fermentation targets of the ASPECT proposal's budget-period milestones
(BP1/GNG1, BP2, BP3/GNG2, BP5) are just met, at the 15,000 wet ton/day
baseline, for the 'enzymatic_press' configuration (the only one whose
hydrolysis step this script maps; see below). Each milestone's metrics are
imposed on the 3-HP system as follows.

PREPROCESSING
Carbohydrate retention (>= 85% BP2, >= 90% BP3/GNG2/BP5):
    The target is for preprocessing as a whole, of which pressing is one part
    (alongside stabilization, handling and milling). The model's press loses
    no carbohydrate (solids go to the cake, solubles to the pressate), so the
    total loss, 1 - retention, is lumped at the mill: `MI.loss_frac = 1 -
    retention`. The loss stream has the feed's composition, so it is a uniform
    loss of carbohydrate (and everything else) ahead of the press, and goes to
    the boiler. This replaces the model's 3% milling loss (97% retention).
    BP1 has no retention target and keeps the model's 3%. The milestone's
    0.5 dry t/day design scale is not applied.

HYDROLYSIS (EnzymaticPress, `PR`)
Monomer release (>= 60% BP2, >= 80% BP3/GNG2/BP5):
    Sets `PR.hydrolysis_conversion` (Glucan -> Glucose and Alginate ->
    AlginateMonomer). Solubilization and monomer release targets are equal
    in every milestone, so the one value is used. BP1 keeps the model's 50%.
Enzyme loading (<= 25 mg/g dry weight biomass from BP2 on):
    Sets `PR.enzyme_loading_mg_per_g_substrate` so that the total enzyme
    dose (including the model's 10% excess) is 25 mg per g of dry biomass fed
    to the press. The model's own dose is only 0.5 mg/g glucan+alginate, far
    below that cap, so the cap is the worst case. Each milestone is also run
    at the model's own dose ("msp_model_enzyme_usd_per_kg").
Fermentable carbohydrate concentration (>= 40 g/L BP2, 60 BP3/GNG2, 120 BP5):
    Taken as the sugar (Glucose + Mannitol + AlginateMonomer) concentration of
    the fermentation feed, solved through `EV.target_solids_wt_frac` together
    with the titer (below); whichever of the two is binding sets the feed
    concentration.
Hydrolysis time (<= 96 h, fed batch per the project team):
    Not represented: EnzymaticPress is an instantaneous step with no tank,
    so there is no time cost to constrain.

FERMENTATION
The three fermentation metrics of each milestone are imposed as follows.

Yield (g 3-HP / g sugar):
    Sets each sugar's `hp_conversion` on the HPFermentation (F401) instance,
    as in msp_yield_productivity_3hp.py: conversion = overall yield /
    `hp_yield_kg_per_kg_consumed`. The model has no galactose, so Glucose gets
    the milestone's glucose yield and Mannitol and AlginateMonomer get its
    mannitol/galactose yield (the same value for all three from BP3 on, where
    the milestone is stated per consumed carbohydrate).
    Milestones stated on a "consumed carbohydrate" basis (BP3, GNG2, BP5) are
    multiplied by the milestone's carbohydrate-utilization target (BP3 and
    GNG2: 85%, BP5: 90%) to get the overall yield on sugar supplied, which is
    what `hp_conversion` represents (see data/3hp.yaml, glucose_conversion).
    BP1 and BP2 yields are taken as overall g/g with no extra utilization
    factor.

Titer (g/L):
    Not a model input: it follows from the sugar concentration of the
    fermentation feed and the yield. It is hit by solving
    `EV.target_solids_wt_frac` (the fermentation-feed solids target, see
    msp_fermentation_solids_3hp.py) with a secant iteration until the binding
    one of {simulated "Titer (g HP-equivalent/L)" vs the titer target, feed
    sugar concentration vs the fermentable-carbohydrate target} is within 1%
    of its milestone. The milestones are minimums, so this is the least
    concentrated (most costly) feed that still meets them.

Productivity (g/L/h):
    Set on `F401.productivity_g_per_L_per_h` (only enters through batch time,
    i.e. fermenter size).

BP2 gives different titer and productivity for glucose (40 g/L, 0.24 g/L/h)
and for mannitol/galactose (15 g/L, 0.10 g/L/h), but the model has one
fermentation broth, so both are run: "BP2 (glucose titer/productivity)" and
"BP2 (mannitol/galactose titer/productivity)". Yields are per sugar in both.
BP4 has no fermentation-relevant targets (scale-up only) and is not run. The
0.5 dry t/day pre-pilot scale targets are not applied (the model is at the
commercial baseline scale). The proposal's alginate/AHG
targets (qualitative, or 5 g/L) are not represented.

The unmodified model ("current model") is included for reference.

LCA
Each case also reports cradle-to-gate GWP100 (kg CO2e per kg Ca(3HP)2 and per
kg 3-HP) from sabre/_lca.py, with and without the credit for exported
electricity, and the enzyme, lime, natural gas and electricity contributions.
Characterization factors are placeholders (data/lca.yaml); the Sargassum
feedstock carries no burden, and wastewater, mother-liquor and ash treatment
are not included. See sabre/_lca.py for the full scope.

Outputs (relative to this script's directory):
    ../results/tables/msp_milestones_3hp_<configuration>.csv
    ../results/figures/fig_msp_milestones_3hp_<configuration>.png

Run from the repo root (optional argument: configuration, default
enzymatic_press):
    python biorefineries/sabre/analyses/msp_milestones_3hp.py
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

from biorefineries.sabre import _lca
from biorefineries.sabre._tea import solve_product_msp
from biorefineries.sabre.systems import _3hp_system
from biorefineries.sabre.analyses.msp_scale_scenarios_3hp import (
    SURFACE, INK, INK_2, GRID, AXIS, SERIES_WITH, TIGHT_TOLERANCE,
)

TABLES_OUT = SCRIPT_DIR.parent / "results" / "tables"
FIGURES_OUT = SCRIPT_DIR.parent / "results" / "figures"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
FIGURES_OUT.mkdir(parents=True, exist_ok=True)

TITER_RTOL = 0.01
MAX_TITER_ITER = 6
SOLIDS_BOUNDS = (0.02, 0.30)  # fermentation-feed solids wt frac search range

SUGAR_IDS = ("Glucose", "Mannitol", "AlginateMonomer")

# label, titer [g/L], yield for Glucose, yield for Mannitol/AlginateMonomer
# [g/g], carbohydrate utilization (1 = yields are already overall),
# productivity [g/L/h]; preprocessing/hydrolysis: carbohydrate retention,
# monomer release (hydrolysis conversion), fermentable carbohydrate [g/L],
# enzyme cap [mg/g dry biomass]. None = keep the model's value (BP1).
_BP2 = dict(retention=0.85, hydrolysis=0.60, carbs=40, enzyme=25)
_BP3 = dict(retention=0.90, hydrolysis=0.80, carbs=60, enzyme=25)
MILESTONES = (
    dict(label="BP1 / GNG1", titer=25, y_glucose=0.20, y_other=0.20, utilization=1.00, productivity=0.14,
         retention=None, hydrolysis=None, carbs=None, enzyme=None),
    dict(label="BP2 (glucose titer/productivity)", titer=40, y_glucose=0.45, y_other=0.20, utilization=1.00, productivity=0.24, **_BP2),
    dict(label="BP2 (mannitol/galactose titer/productivity)", titer=15, y_glucose=0.45, y_other=0.20, utilization=1.00, productivity=0.10, **_BP2),
    dict(label="BP3", titer=23, y_glucose=0.45, y_other=0.45, utilization=0.85, productivity=0.20, **_BP3),
    dict(label="GNG2", titer=40, y_glucose=0.45, y_other=0.45, utilization=0.85, productivity=0.25, **_BP3),
    dict(label="BP5", titer=60, y_glucose=0.55, y_other=0.55, utilization=0.90, productivity=0.40,
         **{**_BP3, "carbs": 120}),
)


def _simulate(sys):
    """Simulate; return (titer [g/L], fermentation-feed sugar [g/L])."""
    sys.simulate()
    fu, fs = sys.flowsheet.unit, sys.flowsheet.stream
    feed = fs.fermentation_feed
    sugar = sum(float(feed.imass[i]) for i in SUGAR_IDS) / float(feed.F_vol)  # kg/h / (m3/h) = g/L
    return float(fu.F401.design_results["Titer (g HP-equivalent/L)"]), sugar


def run_case(configuration: str, milestone: dict | None = None, enzyme_cap: bool = True) -> dict:
    """
    Build the system, impose a milestone's preprocessing/hydrolysis/
    fermentation targets (all unmodified if `milestone` is None), and solve
    the MSP. `enzyme_cap=False` keeps the model's own enzyme dose instead of
    the milestone's mg/g dry biomass cap.
    """
    if configuration != "enzymatic_press":
        raise ValueError("Only the 'enzymatic_press' configuration's hydrolysis step is mapped.")
    bst.main_flowsheet.clear()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sys = _3hp_system.create_3hp_system(configuration=configuration)
        sys.set_tolerance(subsystems=True, **TIGHT_TOLERANCE)
        F401 = sys.flowsheet.unit.F401
        EV = sys.flowsheet.unit.EV
        MI = sys.flowsheet.unit.MI
        PR = sys.flowsheet.unit.PR

        if milestone is not None:
            # preprocessing / hydrolysis
            if milestone["retention"] is not None:
                MI.loss_frac = 1.0 - milestone["retention"]
            if milestone["hydrolysis"] is not None:
                PR.hydrolysis_conversion = milestone["hydrolysis"]
                for rxn in PR.hydrolysis_rxns:
                    rxn.X = milestone["hydrolysis"]
            # fermentation
            for sugar, params in F401.substrates.items():
                y = milestone["y_glucose"] if sugar == "Glucose" else milestone["y_other"]
                X = y * milestone["utilization"] / params["hp_yield_kg_per_kg_consumed"]
                if X + params["biomass_conversion"] > 1.0:
                    raise ValueError(f"{milestone['label']}: {sugar} needs conversion {X:.3f}, > 1 with biomass")
                params["hp_conversion"] = X
            F401.productivity_g_per_L_per_h = milestone["productivity"]

        titer, sugar = _simulate(sys)
        if milestone is not None:
            if enzyme_cap and milestone["enzyme"] is not None:
                # Total dose (incl. excess) = cap mg per g dry biomass fed to the press.
                feed = PR.ins[0]
                dry = float(feed.F_mass) - float(feed.imass["Water"])
                substrate = float(feed.imass["Glucan"] + feed.imass["Alginate"])
                PR.enzyme_loading_mg_per_g_substrate = (
                    milestone["enzyme"] * dry / substrate / (1.0 + PR.enzyme_excess_frac)
                )
                titer, sugar = _simulate(sys)

            # Secant iteration on the feed-solids target. Both titer and sugar rise
            # with it and both milestones must be met, so the binding (smaller)
            # ratio of {titer, sugar} to its milestone is driven to 1.
            carbs = milestone["carbs"] or 0.0

            def ratio(t, s):
                return min(
                    t / milestone["titer"] if milestone["titer"] else np.inf,
                    s / carbs if carbs else np.inf,
                )

            r0, x0 = ratio(titer, sugar), EV.target_solids_wt_frac
            x1 = float(np.clip(x0 / r0, *SOLIDS_BOUNDS))
            EV.target_solids_wt_frac = x1
            titer, sugar = _simulate(sys)
            r1 = ratio(titer, sugar)
            for _ in range(MAX_TITER_ITER):
                if abs(r1 - 1.0) <= TITER_RTOL or r1 == r0 or x1 == x0:
                    break
                x2 = float(np.clip(x1 + (1.0 - r1) * (x1 - x0) / (r1 - r0), *SOLIDS_BOUNDS))
                if x2 == x1:
                    break
                x0, r0, x1 = x1, r1, x2
                EV.target_solids_wt_frac = x1
                titer, sugar = _simulate(sys)
                r1 = ratio(titer, sugar)
            if abs(r1 - 1.0) > TITER_RTOL:
                raise RuntimeError(
                    f"{milestone['label']}: titer {titer:.1f} g/L, sugar {sugar:.1f} g/L did not converge "
                    f"(feed solids target {x1:.3f})"
                )

        product = sys.flowsheet.stream.ca3hp2_product
        msp = solve_product_msp(sys.TEA, product)
        _lca.set_characterization_factors(sys)
        gwp = _lca.get_gwp(sys, product)
        chem = product.chemicals
        hp_basis = 2 * chem.HP.MW / chem.Ca3HP2.MW  # kg 3-HP per kg pure Ca(3HP)2

    dr = F401.design_results
    return {
        "milestone": milestone["label"] if milestone else "current model",
        "retention": 1.0 - MI.loss_frac,
        "monomer_release": PR.hydrolysis_conversion,
        "enzyme_mg_per_g_dwb": (
            PR.enzyme_loading_mg_per_g_substrate * (1.0 + PR.enzyme_excess_frac)
            * (float(PR.ins[0].imass["Glucan"] + PR.ins[0].imass["Alginate"]))
            / (float(PR.ins[0].F_mass) - float(PR.ins[0].imass["Water"]))
        ),
        "target_sugar_g_per_L": milestone["carbs"] if milestone else np.nan,
        "sugar_g_per_L": sugar,
        "target_titer_g_per_L": (milestone["titer"] or np.nan) if milestone else np.nan,
        "titer_g_per_L": titer,
        "yield_g_per_g": dr["HP formed (kg/h)"] / dr["Substrate available (kg/h)"],
        "productivity_g_per_L_per_h": F401.productivity_g_per_L_per_h,
        "feed_solids_target_wt_frac": EV.target_solids_wt_frac,
        "fermentation_batch_h": float(F401.tau),
        "FCI_usd": sys.TEA.FCI,
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
        "msp_usd_per_kg": msp["usd_per_kg"],
        "msp_usd_per_kg_3hp_basis": msp["usd_per_kg"] / hp_basis,
        **gwp,
    }


def collect_results(configuration: str) -> pd.DataFrame:
    rows = [run_case(configuration)]
    for m in MILESTONES:
        row = run_case(configuration, m)
        # same milestone at the model's own (much lower) enzyme dose
        alt = run_case(configuration, m, enzyme_cap=False) if m["enzyme"] is not None else row
        row["msp_model_enzyme_usd_per_kg"] = alt["msp_usd_per_kg"]
        row["gwp_model_enzyme_kg_per_kg"] = alt["gwp_kg_per_kg"]
        rows.append(row)
    df = pd.DataFrame(rows)
    cur = df["milestone"] == "current model"
    df.loc[cur, "msp_model_enzyme_usd_per_kg"] = df["msp_usd_per_kg"]
    df.loc[cur, "gwp_model_enzyme_kg_per_kg"] = df["gwp_kg_per_kg"]
    return df


def collect_sugar_sensitivity(configuration: str, sugar_targets=(60, 80, 100, 120)) -> pd.DataFrame:
    """
    BP5 with the fermentable-carbohydrate target relaxed and the titer target
    waived (titer then follows from the sugar level and the yield), since the
    BP5 titer of 60 g/L at 0.495 g/g already implies ~120 g/L of sugar.
    """
    bp5 = next(m for m in MILESTONES if m["label"] == "BP5")
    rows = []
    for c in sugar_targets:
        m = {**bp5, "label": f"BP5, sugar >= {c} g/L (titer free)", "carbs": c, "titer": None}
        rows.append(run_case(configuration, m))
    return pd.DataFrame(rows)


def plot_results(df: pd.DataFrame, configuration: str, out_path: Path):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "figure.dpi": 100, "savefig.dpi": 200,
    })
    fig, ax = plt.subplots(figsize=(9.6, 5.4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    labels = [m.replace(" (", "\n(") for m in df["milestone"]]
    vals = df["msp_usd_per_kg"]
    colors = [AXIS if m == "current model" else SERIES_WITH for m in df["milestone"]]
    bars = ax.bar(range(len(df)), vals, color=colors, zorder=3)
    span = max(vals.max(), 0) - min(vals.min(), 0)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.015 * span, f"{v:.2f}",
                ha="center", va="bottom", color=INK, fontsize=10)
    ax.plot(range(len(df)), df["msp_model_enzyme_usd_per_kg"], linestyle="none", marker="D",
            markersize=6, color=INK, zorder=4, label="at the model's own enzyme dose")
    ax.legend(frameon=False, loc="upper right", labelcolor=INK_2, fontsize=9)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(labels, color=INK_2, fontsize=8)
    ax.set_ylabel(r"MSP (USD/kg Ca(3HP)$_2$)", color=INK_2)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)
    ax.tick_params(axis="y", length=0, colors=INK_2)
    ax.tick_params(axis="x", length=0)
    ax.set_title(f"Ca(3HP)$_2$ MSP at the ASPECT milestone targets ({configuration})",
                 color=INK, fontsize=12, loc="left")
    fig.tight_layout()
    fig.savefig(out_path, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    configuration = sys.argv[1] if len(sys.argv) > 1 else "enzymatic_press"
    df = collect_results(configuration)
    csv_path = TABLES_OUT / f"msp_milestones_3hp_{configuration}.csv"
    fig_path = FIGURES_OUT / f"fig_msp_milestones_3hp_{configuration}.png"
    df.to_csv(csv_path, index=False)
    sens = collect_sugar_sensitivity(configuration)
    sens_path = TABLES_OUT / f"msp_milestones_3hp_bp5_sugar_{configuration}.csv"
    sens.to_csv(sens_path, index=False)
    plot_results(df, configuration, fig_path)
    print(f"Saved table:  {csv_path}")
    print(f"Saved figure: {fig_path}\n")
    with pd.option_context("display.width", 220, "display.max_columns", None):
        print(df.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))
        print("\nBP5 fermentable-sugar sensitivity:")
        print(sens.to_string(index=False, float_format=lambda v: f"{v:,.4g}"))

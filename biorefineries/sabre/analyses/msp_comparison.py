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
    - AD-biomethane    -> biomethane, $/kg CH4
    - AD-VFA           -> VFA broth, $/kg total VFA
    - AD-fermentation  -> crude microbial oil, $/kg MicrobialOil

Each pathway's product is chemically different (an aqueous seaweed-extract
concentrate, a purified gas, a mixed organic-acid broth, and a crude
lipid), so "MSP" here means the same thing it does in
legacy_analyses/ad_feed_tea.py's product-scenario comparison: the price
per kg of the pathway's own value-carrying product needed to hit the
target IRR, assuming no revenue from any other product. It is a
production-cost comparison, not a claim that the four products compete
in the same market.

All four cases share one baseline economic basis (feed price, disposal
costs, reagent/nutrient prices) already established in
legacy_analyses/ad_feed_tea.py and legacy_analyses/vfa_fermentation_tea.py,
so the comparison sits on a common footing rather than four independently
tuned base cases.

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

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.systems import (
    create_biostimulant_system,
    create_ad_biomethane_system,
    create_ad_vfa_system,
    create_ad_fermentation_system,
)
from biorefineries.sabre._tea import create_tea, solve_product_msp, CH4_MMBTU_PER_KG

OUT = SCRIPT_DIR.parent / "results" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# =================================================================
# Shared baseline economics
# Values match the defaults already used across
# legacy_analyses/ad_feed_tea.py and legacy_analyses/vfa_fermentation_tea.py
# so all four systems are compared on a common basis.
# =================================================================
FEED_PRICE_USD_PER_KG_WET = 0.00  # near-zero feedstock price (base case)

SOLIDS_DISPOSAL_USD_PER_KG           = -0.04   # wet solids cake, "standard biosolids" basis
SOIL_AMENDMENT_DISPOSAL_USD_PER_KG   = -0.02   # AD digestate solids; heavy metals rule out land application
LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG = -0.002  # ~$2/m3 wastewater basis
FERM_WASTEWATER_DISPOSAL_USD_PER_KG  = -0.005  # fermentation wastewater basis; also covers
                                                # evaporator condensate + microfilter retentate,
                                                # which are now mixed into this stream upstream
BIOSTIMULANT_CREDIT_USD_PER_KG       = 0.00    # no byproduct revenue in the base case

OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL = 0.50

PRICE_MGSO4_USD_PER_KG     = 0.40
PRICE_AMMONIA_USD_PER_KG   = 0.50
PRICE_PHOSPHATE_USD_PER_KG = 0.80
PRICE_NAOH_USD_PER_KG      = 0.35

VFA_IDS = ["AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"]


def _set_price(flowsheet, stream_id: str, price: float):
    """Best-effort price assignment; a stream may not exist in every case."""
    try:
        flowsheet.stream[stream_id].price = price
    except Exception:
        pass


def _patch_ev607():
    """
    Replace Ev607's cost + utility with a low-duty placeholder when V < 0.02.
    BioSTEAM's MultiEffectEvaporator cost correlation produces nonsensical
    vessel geometry at near-zero evaporation duty (see
    legacy_analyses/vfa_fermentation_tea.py::_patch_ev607).
    """
    try:
        ev607 = bst.main_flowsheet.unit["Ev607"]
    except Exception:
        return
    v = getattr(ev607, "V", None)
    if v is None or v >= 0.02:
        return
    feed = ev607.ins[0]
    feed_m3h = max(feed.F_mass / 1000.0, 1.0)
    placeholder_usd = 50000.0 * (feed_m3h ** 0.6)
    for k in list(ev607.baseline_purchase_costs.keys()):
        ev607.baseline_purchase_costs[k] = 0.0
    ev607.baseline_purchase_costs["Evaporator (low-duty placeholder)"] = placeholder_usd
    ev607.heat_utilities.clear()


# =================================================================
# System builders -- each returns a common result dict
# =================================================================

def price_biostimulant_system() -> dict:
    bst.main_flowsheet.clear()
    create_chemicals()

    sys, streams, units = create_biostimulant_system()
    streams["feed"].price = FEED_PRICE_USD_PER_KG_WET
    streams["pressed_cake"].price = SOLIDS_DISPOSAL_USD_PER_KG
    sys.simulate()

    tea = create_tea(sys)
    product = streams["biostimulant_product"]
    msp = solve_product_msp(tea=tea, product_stream=product)

    annual_hours = tea.operating_days * 24
    return {
        "label": "Biostimulant",
        "product_desc": "biostimulant liquid product",
        "msp_usd_per_kg": msp["usd_per_kg_stream"],
        "tea": tea,
        "annual_product_t": product.F_mass * annual_hours / 1000.0,
    }


def price_ad_biomethane_system() -> dict:
    bst.main_flowsheet.clear()
    create_chemicals()

    sys = create_ad_biomethane_system(feedstock="pelagic")
    sys.flowsheet.stream.sargassum_feed.price = FEED_PRICE_USD_PER_KG_WET

    fs = sys.flowsheet
    _set_price(fs, "residual_permeate", 0.0)
    _set_price(fs, "biostimulant_product", BIOSTIMULANT_CREDIT_USD_PER_KG)
    _set_price(fs, "soil_amendment", SOIL_AMENDMENT_DISPOSAL_USD_PER_KG)
    _set_price(fs, "liquid_digestate", LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG)

    sys.simulate()

    tea = create_tea(sys)
    biomethane = fs.stream.biomethane
    msp = solve_product_msp(
        tea=tea, product_stream=biomethane, product_ID="Methane",
        energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
    )

    return {
        "label": "AD-biomethane",
        "product_desc": "biomethane (CH4 basis)",
        "msp_usd_per_kg": msp["usd_per_kg_product"],
        "tea": tea,
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
    }


def price_ad_vfa_system() -> dict:
    bst.main_flowsheet.clear()
    create_chemicals()

    sys = create_ad_vfa_system(feedstock="pelagic")
    sys.flowsheet.stream.sargassum_feed.price = FEED_PRICE_USD_PER_KG_WET

    fs = sys.flowsheet
    _set_price(fs, "residual_permeate", 0.0)
    _set_price(fs, "biostimulant_product", BIOSTIMULANT_CREDIT_USD_PER_KG)
    _set_price(fs, "acidogenic_residual_solids", SOLIDS_DISPOSAL_USD_PER_KG)

    sys.simulate()

    tea = create_tea(sys)
    vfa_broth = fs.stream.vfa_broth
    msp_stream = tea.solve_price(vfa_broth)

    vfa_mass_frac = sum(
        float(vfa_broth.imass[cid]) for cid in VFA_IDS if cid in vfa_broth.chemicals
    ) / float(vfa_broth.F_mass)
    msp_per_kg_vfa = msp_stream / vfa_mass_frac if vfa_mass_frac > 0 else float("nan")

    annual_hours = tea.operating_days * 24
    annual_vfa_t = vfa_broth.F_mass * vfa_mass_frac * annual_hours / 1000.0

    return {
        "label": "AD-VFA",
        "product_desc": "VFA broth (total-VFA basis)",
        "msp_usd_per_kg": msp_per_kg_vfa,
        "tea": tea,
        "annual_product_t": annual_vfa_t,
    }


def price_ad_fermentation_system() -> dict:
    bst.main_flowsheet.clear()
    create_chemicals()

    sys, streams, units = create_ad_fermentation_system(feedstock="pelagic")
    streams["feed"].price = FEED_PRICE_USD_PER_KG_WET

    fs = sys.flowsheet
    for sid, price in {
        "fermentation_mgso4":     PRICE_MGSO4_USD_PER_KG,
        "fermentation_ammonia":   PRICE_AMMONIA_USD_PER_KG,
        "fermentation_phosphate": PRICE_PHOSPHATE_USD_PER_KG,
        "fermentation_base":      PRICE_NAOH_USD_PER_KG,
    }.items():
        _set_price(fs, sid, price)

    sys.simulate()
    _patch_ev607()

    _set_price(fs, "residual_permeate", 0.0)
    _set_price(fs, "biostimulant_product", BIOSTIMULANT_CREDIT_USD_PER_KG)
    _set_price(fs, "acidogenic_residual_solids", SOLIDS_DISPOSAL_USD_PER_KG)
    _set_price(fs, "wastewater", FERM_WASTEWATER_DISPOSAL_USD_PER_KG)

    oil_stream = streams["microbial_oil"]
    oil_kg_hr = float(oil_stream.imass["MicrobialOil"])
    try:
        fs.unit["OE"].add_OPEX = {
            "Oil extraction reagent": oil_kg_hr * OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL
        }
    except Exception:
        pass

    tea = create_tea(sys)
    msp = solve_product_msp(tea=tea, product_stream=oil_stream, product_ID="MicrobialOil")

    return {
        "label": "AD-fermentation",
        "product_desc": "crude microbial oil",
        "msp_usd_per_kg": msp["usd_per_kg_product"],
        "tea": tea,
        "annual_product_t": msp["annual_product_kg"] / 1000.0,
    }


BUILDERS = (
    price_biostimulant_system,
    price_ad_biomethane_system,
    price_ad_vfa_system,
    price_ad_fermentation_system,
)


# =================================================================
# Main
# =================================================================

if __name__ == "__main__":
    results = [builder() for builder in BUILDERS]

    print("\n" + "=" * 78)
    print(f"SaBRe FOUR-SYSTEM MSP COMPARISON  (feed = ${FEED_PRICE_USD_PER_KG_WET:.2f}/kg wet)")
    print("=" * 78)
    print(f"  {'Pathway':<18} {'Product basis':<28} {'MSP ($/kg)':>12} {'TCI ($M)':>10} {'Output (t/yr)':>14}")
    print("  " + "-" * 84)
    for r in results:
        print(
            f"  {r['label']:<18} {r['product_desc']:<28} "
            f"${r['msp_usd_per_kg']:>10.2f} ${r['tea'].TCI/1e6:>8.1f} "
            f"{r['annual_product_t']:>13,.0f}"
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
        "AD-biomethane":    "$/kg CH4",
        "AD-VFA":           "$/kg total VFA",
        "AD-fermentation":  "$/kg crude oil",
    }
    tick_labels = [
        f"{r['label']}\n({UNIT_SUBLABELS[r['label']]})" for r in results
    ]
    msp_vals = [r["msp_usd_per_kg"] for r in results]

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars = ax.bar(tick_labels, msp_vals, edgecolor="black", linewidth=0.8, zorder=3)

    for bar, val in zip(bars, msp_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(msp_vals) * 0.015,
            f"${val:,.2f}",
            ha="center", va="bottom", fontsize=9,
        )

    ax.set_ylabel("Minimum selling price ($/kg product)")
    ax.set_title(
        "MSP by SaBRe pathway (near-zero feed price)\n"
        "Each bar is priced on its own pathway's value-carrying product",
        fontsize=10,
    )
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("$%.0f"))
    ax.grid(axis="y", linewidth=0.4, color="#D3D1C7", zorder=0)
    ax.set_ylim(0, max(msp_vals) * 1.15)

    fig.tight_layout()
    fig.savefig(OUT / "fig_msp_comparison_four_systems.png", bbox_inches="tight")

    print("\nSaved:")
    print(f"  {OUT / 'fig_msp_comparison_four_systems.png'}")

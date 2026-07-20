# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Integrated Biorefinery TEA
===========================================
Sweeps the split fraction alpha (0 → 1) between the two pathways:
  alpha=0   → 100% VFA-to-oil
  alpha=0.5 → 50/50 split
  alpha=1   → 100% biomethane

For each alpha, reports:
  - Biomethane MSP ($/MMBtu)
  - Microbial oil MSP ($/kg)
  - Combined NPV at assumed market prices
  - TCI and VOC breakdown

Products:
  - Biomethane at assumed $10/MMBtu (European TTF benchmark)
  - Microbial oil at assumed $1.00/kg (midpoint of $0.62–1.50/kg soybean oil range)
  - Biostimulant at $0.50/kg (standard assumption)
"""

import sys
from pathlib import Path

import biosteam as bst
import math

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.systems import create_ad_integrated_system
from biorefineries.sabre._tea import create_tea, solve_product_msp, solve_biomethane_msp

# -------------------------
# Market price assumptions
# -------------------------
BIOMETHANE_MARKET_MMBTU    = 10.00   # $/MMBtu — European TTF benchmark
OIL_MARKET_USD_PER_KG      =  1.00   # $/kg    — midpoint of soybean oil $0.62–1.50/kg range
BIOSTIMULANT_USD_PER_KG    =  0.50   # $/kg    — standard assumption

# Biomethane price sensitivity range
BIOMETHANE_PRICES = [3.0, 10.0, 14.0]   # Henry Hub, TTF, JKM

# -------------------------
# Reagent cost (oil extraction)
# -------------------------
OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL = 0.50

# -------------------------
# Disposal costs
# -------------------------
# Covers fermentation wastewater plus evaporator condensate and microfilter
# retentate, mixed upstream into the same stream
FERM_WASTEWATER_DISPOSAL_USD_PER_KG = -0.005
SOLIDS_DISPOSAL_USD_PER_KG          = -0.04
LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG = -0.002
SOLID_DIGESTATE_DISPOSAL_USD_PER_KG  = -0.02

# -------------------------
# Feed price base case
# -------------------------
FEED_PRICE_BASE = 0.02   # $/kg wet — standard collection cost assumption

# -------------------------
# Alpha sweep points
# -------------------------
ALPHA_SWEEP = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

CH4_MMBTU_PER_KG = 0.0526


def _get_integrated_stream(stream_id: str):
    try:
        return bst.main_flowsheet.stream[stream_id]
    except Exception:
        return None


def _apply_stream_prices(streams, biostimulant_price=BIOSTIMULANT_USD_PER_KG):
    if streams.get("biostimulant_product") is not None:
        streams["biostimulant_product"].price = biostimulant_price

    permeate = _get_integrated_stream("pressate_permeate")
    if permeate is not None:
        permeate.price = 0.0

    for sid, price in [
        ("wastewater",                 FERM_WASTEWATER_DISPOSAL_USD_PER_KG),
        ("acidogenic_residual_solids", SOLIDS_DISPOSAL_USD_PER_KG),
    ]:
        s = streams.get(sid)
        if s is not None:
            s.price = price

    for sid, price in [
        ("soil_amendment",   SOLID_DIGESTATE_DISPOSAL_USD_PER_KG),
        ("liquid_digestate", LIQUID_DIGESTATE_DISPOSAL_USD_PER_KG),
    ]:
        s = streams.get(sid)
        if s is not None:
            s.price = price


def _wire_oil_reagent(streams, units):
    oil_stream = streams.get("backend_oil")
    oe = units.get("OE")
    if oil_stream is not None and oe is not None:
        oil_kg_hr = float(oil_stream.imass["MicrobialOil"])
        reagent_usd_per_hr = oil_kg_hr * OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL
        oe.add_OPEX = {"Oil extraction reagent": reagent_usd_per_hr}


def _patch_ev607():
    try:
        ev607 = bst.main_flowsheet.unit["Ev607"]
        v = getattr(ev607, "V", None)
        if v is not None and v < 0.02:
            feed = ev607.ins[0]
            feed_m3h = max(feed.F_mass / 1000.0, 1.0)
            placeholder_usd = 50000.0 * (feed_m3h ** 0.6)
            for k in list(ev607.baseline_purchase_costs.keys()):
                ev607.baseline_purchase_costs[k] = 0.0
            ev607.baseline_purchase_costs["Evaporator (low-duty placeholder)"] = placeholder_usd
            ev607.power_utility.consumption = 0.0
            ev607.heat_utilities.clear()
    except Exception:
        pass


def _compute_npv_at_market(tea, streams, market_mmbtu, market_oil_usd_per_kg):
    biomethane = streams.get("biomethane")
    oil_stream = streams.get("backend_oil")

    old_bm_price  = biomethane.price  if biomethane  is not None else None
    old_oil_price = oil_stream.price  if oil_stream  is not None else None

    try:
        if biomethane is not None and float(biomethane.F_mass) > 0:
            ch4_mass   = float(biomethane.imass["Methane"])
            total_mass = float(biomethane.F_mass)
            ch4_frac   = ch4_mass / total_mass if total_mass > 0 else 0.0
            biomethane.price = market_mmbtu * CH4_MMBTU_PER_KG * ch4_frac

        if oil_stream is not None and float(oil_stream.F_mass) > 0:
            oil_mass_frac = float(oil_stream.imass["MicrobialOil"]) / float(oil_stream.F_mass)
            oil_stream.price = market_oil_usd_per_kg * oil_mass_frac

        npv = tea.NPV

    finally:
        if biomethane  is not None and old_bm_price  is not None:
            biomethane.price  = old_bm_price
        if oil_stream  is not None and old_oil_price is not None:
            oil_stream.price  = old_oil_price

    return npv


def run_alpha_sweep(
    feed_price: float = FEED_PRICE_BASE,
    case_label: str = "base",
    pretreatment_case: str = "press_mill_only",
    biostimulant_price: float = BIOSTIMULANT_USD_PER_KG,
    market_mmbtu: float = BIOMETHANE_MARKET_MMBTU,
    market_oil: float = OIL_MARKET_USD_PER_KG,
    print_summary: bool = True,
):
    results = []

    for alpha in ALPHA_SWEEP:
        bst.main_flowsheet.clear()
        set_thermo()

        try:
            sys, streams, units, _ = create_ad_integrated_system(
                alpha=alpha,
                pretreatment_case=pretreatment_case,
            )
            streams["feed"].price = feed_price
            sys.simulate()

            _patch_ev607()
            _apply_stream_prices(streams, biostimulant_price)
            _wire_oil_reagent(streams, units)

            tea = create_tea(sys)

            biomethane = streams.get("biomethane")
            oil_stream = streams.get("backend_oil")
            msp_mmbtu = float("nan")
            msp_ch4   = float("nan")
            msp_oil   = float("nan")
            oil_kg_yr = 0.0

            # Biomethane MSP (oil price = 0)
            if biomethane is not None and alpha > 0 and float(biomethane.F_mass) > 0:
                old_oil_price = oil_stream.price if oil_stream is not None else None
                try:
                    if oil_stream is not None:
                        oil_stream.price = 0.0
                    bm_msp    = solve_biomethane_msp(tea, biomethane)
                    msp_mmbtu = bm_msp.get("usd_per_mmbtu", float("nan"))
                    msp_ch4   = bm_msp.get("usd_per_kg_ch4", float("nan"))
                finally:
                    if oil_stream is not None and old_oil_price is not None:
                        oil_stream.price = old_oil_price

            # Oil MSP (biomethane price = 0)
            if oil_stream is not None and alpha < 1.0 and float(oil_stream.F_mass) > 0:
                old_bm_price = biomethane.price if biomethane is not None else None
                try:
                    if biomethane is not None:
                        biomethane.price = 0.0
                    oil_msp  = solve_product_msp(tea, oil_stream, product_ID="MicrobialOil")
                    msp_oil  = oil_msp.get("usd_per_kg_product", float("nan"))
                    oil_kg_yr = oil_msp.get("annual_product_kg", 0.0)
                finally:
                    if biomethane is not None and old_bm_price is not None:
                        biomethane.price = old_bm_price

            npv = _compute_npv_at_market(tea, streams, market_mmbtu, market_oil)

            row = {
                "alpha": alpha,
                "msp_biomethane_mmbtu": msp_mmbtu,
                "msp_biomethane_ch4":   msp_ch4,
                "msp_oil_usd_per_kg":   msp_oil,
                "combined_npv_M":       npv / 1e6,
                "tci_M":                tea.TCI / 1e6,
                "voc_M":                tea.VOC / 1e6,
                "foc_M":                tea.FOC / 1e6,
                "oil_kg_yr":            oil_kg_yr,
                "ok": True,
            }

        except Exception as e:
            row = {
                "alpha": alpha,
                "msp_biomethane_mmbtu": float("nan"),
                "msp_biomethane_ch4":   float("nan"),
                "msp_oil_usd_per_kg":   float("nan"),
                "combined_npv_M":       float("nan"),
                "tci_M":                float("nan"),
                "voc_M":                float("nan"),
                "foc_M":                float("nan"),
                "oil_kg_yr":            0.0,
                "ok":    False,
                "error": str(e),
            }
            print(f"  [alpha={alpha:.1f}] ERROR: {e}")

        results.append(row)

    if print_summary:
        _print_sweep_table(results, case_label, pretreatment_case,
                           feed_price, biostimulant_price, market_mmbtu, market_oil)

    return results


def _print_sweep_table(results, case_label, pretreatment_case,
                       feed_price, biostimulant_price, market_mmbtu, market_oil):
    print("\n" + "=" * 95)
    print(
        f"INTEGRATED BIOREFINERY — ALPHA SWEEP\n"
        f"  Feed: ${feed_price:.3f}/kg [{case_label}]  |  "
        f"Pretreatment: {pretreatment_case}  |  "
        f"Biostimulant: ${biostimulant_price:.2f}/kg\n"
        f"  Market: biomethane ${market_mmbtu:.1f}/MMBtu  |  "
        f"microbial oil ${market_oil:.2f}/kg"
    )
    print("=" * 95)
    print(
        f"  {'Alpha':>6}  {'→CH4':>6}  {'→Oil':>6}  "
        f"{'MSP CH4 ($/MMBtu)':>20}  {'MSP Oil ($/kg)':>16}  "
        f"{'NPV @ mkt ($M)':>15}  {'TCI ($M)':>10}"
    )
    print("  " + "-" * 93)

    best_npv = max(
        (r["combined_npv_M"] for r in results if r["ok"]),
        default=float("nan"),
    )

    for r in results:
        if not r["ok"]:
            print(f"  {r['alpha']:>6.1f}  ERROR")
            continue

        alpha = r["alpha"]
        msp_ch4_str = f"${r['msp_biomethane_mmbtu']:.2f}" if not math.isnan(r["msp_biomethane_mmbtu"]) else "  n/a  "
        msp_oil_str = f"${r['msp_oil_usd_per_kg']:.3f}"   if not math.isnan(r["msp_oil_usd_per_kg"])   else "  n/a  "
        star = " ◄ best" if abs(r["combined_npv_M"] - best_npv) < 0.01 else ""

        print(
            f"  {alpha:>6.1f}  {alpha*100:>5.0f}%  {(1-alpha)*100:>5.0f}%  "
            f"{msp_ch4_str:>20}  {msp_oil_str:>16}  "
            f"${r['combined_npv_M']:>13.1f}M  ${r['tci_M']:>8.1f}M{star}"
        )

    valid = [r for r in results if r["ok"] and not math.isnan(r["combined_npv_M"])]
    if valid:
        best = max(valid, key=lambda r: r["combined_npv_M"])
        print(f"\n  Optimal alpha = {best['alpha']:.1f} (NPV = ${best['combined_npv_M']:.1f}M)")


# =============================================================
# Main — runs updated four scenarios at three biomethane prices
# =============================================================
if __name__ == "__main__":

    SCENARIOS = [
        ("base",       0.02,  0.50),
        ("tipping",   -0.02,  0.50),
        ("biostim",    0.02,  1.00),
        ("best_case", -0.02,  1.00),
    ]

    for name, feed, biostim in SCENARIOS:
        for bm_price in BIOMETHANE_PRICES:
            print(f"\n>>> {name} | feed=${feed:.2f}/kg | biostim=${biostim:.2f}/kg | "
                  f"biomethane=${bm_price:.0f}/MMBtu")
            run_alpha_sweep(
                feed_price=feed,
                case_label=name,
                pretreatment_case="press_mill_only",
                biostimulant_price=biostim,
                market_mmbtu=bm_price,
                market_oil=OIL_MARKET_USD_PER_KG,
            )
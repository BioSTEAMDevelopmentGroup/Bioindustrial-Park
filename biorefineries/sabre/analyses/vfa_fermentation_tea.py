# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import sys
from pathlib import Path

import biosteam as bst

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from biorefineries.sabre._chemicals import set_thermo
from biorefineries.sabre.systems import create_vfa_ad_system
from biorefineries.sabre.systems import create_vfa_fermentation_system
from biorefineries.sabre._tea import make_baseline_tea, solve_product_msp

# -------------------------
# Nutrient / reagent prices
# -------------------------
PRICE_MGSO4_USD_PER_KG     = 0.40  # MgSO4 (industrial grade)
PRICE_AMMONIA_USD_PER_KG   = 0.50  # liquid ammonia (N source)
PRICE_PHOSPHATE_USD_PER_KG = 0.80  # KH2PO4 (P source)
PRICE_NAOH_USD_PER_KG      = 0.35  # NaOH (pH base)

# -------------------------
# Disposal costs (negative stream prices = cost in BioSTEAM VOC)
# -------------------------
# Retentate: microfilter reject → modeled as wastewater disposal
#   Basis: industrial wastewater $3-8/m3 ≈ $0.003-0.008/kg
RETENTATE_DISPOSAL_USD_PER_KG = -0.005

# Fermentation wastewater: residual organics + salts + cell debris
#   Basis: biological wastewater treatment $3-8/m3 ≈ $0.003-0.008/kg
FERM_WASTEWATER_DISPOSAL_USD_PER_KG = -0.005

# Acidogenic residual solids: screw press cake
#   Basis: standard biosolids ~$30-80/dry ton ≈ $0.03-0.08/kg wet
#   Base case $0.04/kg — sensitivity tested below
SOLIDS_DISPOSAL_USD_PER_KG = -0.04

# -------------------------
# Oil extraction reagent cost
OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL = 0.50

VFA_IDS = ["AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"]

# -------------------------
# Feed price scenarios
# -------------------------
FEED_PRICE_CASES = [
    ("tipping_fee",      -0.02),
    ("near_zero",         0.00),
    ("low_cost_collect",  0.02),
    ("beach_midpoint",    0.05),
]

# -------------------------
# Sensitivity ranges
# -------------------------
REAGENT_SENSITIVITY   = [0.50, 1.00, 1.50]   # $/kg oil
SOLIDS_SENSITIVITY    = [-0.02, -0.04, -0.10] # $/kg wet cake
BIOSTIMULANT_PRICE_BASE_USD_PER_KG = 0.00
BIOSTIMULANT_SENSITIVITY_USD_PER_KG = [0.00, 0.50, 1.00, 2.00]

# =============================================================
# Helper functions
# =============================================================

def _print_stream_vfa(label: str, stream):
    print(f"\n  [{label}]")
    print(f"    Total F_mass:  {stream.F_mass:>12.2f} kg/hr")
    if "Water" in stream.chemicals:
        print(f"    Water:         {stream.imass['Water']:>12.2f} kg/hr")
    total_vfa = 0.0
    for vfa in VFA_IDS:
        if vfa in stream.chemicals:
            kg = float(stream.imass[vfa])
            total_vfa += kg
            print(f"    {vfa:<22} {kg:>10.2f} kg/hr")
    print(f"    TOTAL VFAs:    {total_vfa:>12.2f} kg/hr")
    return total_vfa


def _patch_ev607(full_sys=None, silent: bool = False):
    """
    Replace Ev607 cost + utility with a low-duty placeholder when V < 0.02.
    MultiEffectEvaporator cost correlation produces nonsensical vessel geometry
    at near-zero evaporation duty.
    """
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
            if not silent:
                print(f"  [Ev607 patch] V={v:.4f} → low-duty placeholder: ${placeholder_usd:,.0f}")
        else:
            if not silent:
                print(f"  [Ev607] V={v:.4f} → cost correlation in bounds, no patch needed")
    except Exception as e:
        if not silent:
            print(f"  [WARNING] Could not patch Ev607: {e}")


def _apply_disposal_costs(
    streams,
    solids_disposal_usd_per_kg: float = SOLIDS_DISPOSAL_USD_PER_KG,
):
    """
    Assign disposal costs to waste outlet streams as negative stream prices.
    Returns a summary dict of {stream_name: annual_cost_usd}.
    """
    annual_hours = 330.0 * 24.0
    summary = {}

    try:
        s = streams["vfa_retentate"]
        s.price = RETENTATE_DISPOSAL_USD_PER_KG
        summary["vfa_retentate"] = abs(s.price * s.F_mass * annual_hours)
    except Exception as e:
        print(f"  [WARNING] retentate disposal: {e}")

    try:
        s = streams["fermentation_wastewater"]
        s.price = FERM_WASTEWATER_DISPOSAL_USD_PER_KG
        summary["fermentation_wastewater"] = abs(s.price * s.F_mass * annual_hours)
    except Exception as e:
        print(f"  [WARNING] ferm wastewater disposal: {e}")

    try:
        s = bst.main_flowsheet.stream["acidogenic_residual_solids"]
        s.price = solids_disposal_usd_per_kg
        summary["acidogenic_residual_solids"] = abs(s.price * s.F_mass * annual_hours)
    except Exception as e:
        print(f"  [WARNING] solids disposal: {e}")

    return summary


# =============================================================
# Core simulation + TEA
# =============================================================

def build_and_simulate(feed_price_per_kg_wet: float):
    bst.main_flowsheet.clear()
    set_thermo()

    vfa_sys = create_vfa_ad_system()
    vfa_sys.feeds[0].price = feed_price_per_kg_wet

    fer_sys, streams, units = create_vfa_fermentation_system(
        vfa_broth=vfa_sys.flowsheet.stream.vfa_broth,
    )

    for sid, price in {
        "fermentation_mgso4":     PRICE_MGSO4_USD_PER_KG,
        "fermentation_ammonia":   PRICE_AMMONIA_USD_PER_KG,
        "fermentation_phosphate": PRICE_PHOSPHATE_USD_PER_KG,
        "fermentation_base":      PRICE_NAOH_USD_PER_KG,
    }.items():
        try:
            bst.main_flowsheet.stream[sid].price = price
        except Exception:
            pass

    full_sys = bst.System.from_units(
        "full_vfa_to_oil_sys",
        units=list(vfa_sys.units) + list(fer_sys.units),
    )
    full_sys.simulate()

    return vfa_sys, fer_sys, streams, units, full_sys

def run_case(
    feed_price_per_kg_wet: float,
    case_label: str = "",
    run_diagnostics: bool = True,
    reagent_usd_per_kg_oil: float = OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL,
    solids_disposal_usd_per_kg: float = SOLIDS_DISPOSAL_USD_PER_KG,
    silent: bool = False,
):

    vfa_sys, fer_sys, streams, units, full_sys = build_and_simulate(feed_price_per_kg_wet)

    _patch_ev607(full_sys)
    disposal_costs = _apply_disposal_costs(
        streams,
        solids_disposal_usd_per_kg=solids_disposal_usd_per_kg,
    )

    biostim_summary = _apply_biostimulant_credit(BIOSTIMULANT_PRICE_BASE_USD_PER_KG)

    oil_stream = streams["backend_oil"]
    oil_kg_hr = float(oil_stream.imass["MicrobialOil"])
    annual_hours = 330.0 * 24.0

    # -------------------------
    # Reagent cost into OE.add_OPEX ($/hr)
    # SABREBaselineTEA.VOC reads _annual_unit_add_opex() which sums
    # add_OPEX across all units
    # -------------------------
    reagent_usd_per_hr = oil_kg_hr * reagent_usd_per_kg_oil
    oil_extraction_reagent_annual = reagent_usd_per_hr * annual_hours  # for printing
    try:
        oe_unit = bst.main_flowsheet.unit["OE"]
        oe_unit.add_OPEX = {"Oil extraction reagent": reagent_usd_per_hr}
    except Exception as e:
        if not silent:
            print(f"  [WARNING] Could not set OE.add_OPEX: {e}")

    tea = make_baseline_tea(full_sys)

    msp = solve_product_msp(
        tea=tea,
        product_stream=oil_stream,
        product_ID="MicrobialOil",
    )

    if silent:
        return tea, msp, streams, units, full_sys

    # ======================================================
    # PRINT RESULTS
    # ======================================================
    label_str = f" [{case_label}]" if case_label else ""
    print("\n" + "=" * 60)
    print(f"FEED PRICE = {feed_price_per_kg_wet:.4f} $/kg wet Sargassum{label_str}")
    print(f"Reagent = ${reagent_usd_per_kg_oil:.2f}/kg oil  |  "
          f"Solids disposal = ${solids_disposal_usd_per_kg:.3f}/kg wet")
    print("=" * 60)

    # -------------------------
    # DIAGNOSTIC: Mass balance trace
    # -------------------------
    if run_diagnostics:
        print("\n" + "-" * 60)
        print("DIAGNOSTIC: VFA MASS BALANCE TRACE")
        print("-" * 60)

        sargassum_feed = vfa_sys.feeds[0]
        print(f"\n  [sargassum_feed]")
        print(f"    Total F_mass:  {sargassum_feed.F_mass:>12.2f} kg/hr")
        try:
            dry = sargassum_feed.F_mass - float(sargassum_feed.imass["Water"])
            print(f"    Dry mass:      {dry:>12.2f} kg/hr")
        except Exception:
            pass

        try:
            vfa_broth = vfa_sys.flowsheet.stream.vfa_broth
            total_vfa_in = _print_stream_vfa(
                "vfa_broth (VFA AD outlet → fermenter inlet)", vfa_broth
            )
        except Exception as e:
            print(f"  Could not read vfa_broth: {e}")
            total_vfa_in = 0.0

        try:
            _print_stream_vfa("vfa_permeate (after microfilter)", streams["vfa_permeate"])
        except Exception as e:
            print(f"  Could not read vfa_permeate: {e}")

        try:
            _print_stream_vfa(
                "conditioned_vfa_broth (fermenter inlet)", streams["conditioned_vfa_broth"]
            )
        except Exception as e:
            print(f"  Could not read conditioned_vfa_broth: {e}")

        try:
            fbroth = streams["fermentation_broth"]
            print(f"\n  [fermentation_broth (R601 outlet)]")
            print(f"    Total F_mass:   {fbroth.F_mass:>12.2f} kg/hr")
            for cid in ["MicrobialOil", "CellMass", "Water"]:
                if cid in fbroth.chemicals:
                    print(f"    {cid:<22} {float(fbroth.imass[cid]):>10.2f} kg/hr")
            residual = sum(
                float(fbroth.imass[v]) for v in VFA_IDS if v in fbroth.chemicals
            )
            print(f"    Residual VFAs:  {residual:>12.2f} kg/hr")
        except Exception as e:
            print(f"  Could not read fermentation_broth: {e}")

        try:
            oil_kg_hr_diag = float(oil_stream.imass["MicrobialOil"])
            print(f"\n  [backend_oil (final product)]")
            print(f"    Total F_mass:   {oil_stream.F_mass:>12.2f} kg/hr")
            print(f"    MicrobialOil:   {oil_kg_hr_diag:>12.2f} kg/hr")
            if total_vfa_in > 0:
                try:
                    r601 = bst.main_flowsheet.unit["R601"]
                    conv = float(r601.conversion)
                    oil_yield = float(r601.product_yield_kg_per_kg_vfa_consumed)
                except Exception:
                    conv = 0.85
                    oil_yield = 0.15

                exp_consumed = total_vfa_in * conv
                exp_oil = exp_consumed * oil_yield

                print(f"\n  [Yield check]")
                print(f"    VFA in (vfa_broth):           {total_vfa_in:>10.2f} kg/hr")
                print(f"    Expected VFA consumed:        {exp_consumed:>10.2f} kg/hr")
                print(f"    Expected oil ({oil_yield:.3f} yield): {exp_oil:>10.2f} kg/hr")
                print(f"    Actual oil produced:          {oil_kg_hr_diag:>10.2f} kg/hr")
                if exp_oil > 0:
                    print(f"    Actual vs expected:           {oil_kg_hr_diag / exp_oil * 100:>10.1f}%")
        except Exception as e:
            print(f"  Could not read backend_oil: {e}")

        print("\n" + "-" * 60)

    # -------------------------
    # Disposal & waste cost summary
    # -------------------------
    print("\n=== DISPOSAL & WASTE COSTS ===")
    for name, cost in disposal_costs.items():
        print(f"  {name:<35} ${cost:>12,.0f}/yr")
    print(f"  {'Oil extraction reagent':<35} ${oil_extraction_reagent_annual:>12,.0f}/yr")
    print(f"  {'TOTAL':<35} ${sum(disposal_costs.values()) + oil_extraction_reagent_annual:>12,.0f}/yr")

    # -------------------------
    # Feeds
    # -------------------------
    print("\n=== FEEDS SEEN BY FULL SYSTEM ===")
    for s in full_sys.feeds:
        if s.F_mass == 0.0 and s.source is None and s.ID.startswith("s"):
            continue
        print(f"  {s.ID:<30} {round(s.F_mass, 3):>12.3f} kg/hr   price={s.price:.4f} $/kg")

    # -------------------------
    # TEA summary
    # -------------------------
    print("\n=== TEA SUMMARY ===")
    print(f"  TCI ($):                     {tea.TCI:>18,.0f}")
    print(f"  FCI ($):                     {tea.FCI:>18,.0f}")
    print(f"  VOC ($/yr):                  {tea.VOC:>18,.0f}")
    print(f"  FOC ($/yr):                  {tea.FOC:>18,.0f}")
    reagent_annual = reagent_usd_per_hr * annual_hours
    print(f"  → Oil extraction reagent:    {reagent_annual:>18,.0f} $/yr (in VOC via OE.add_OPEX)")
    print(f"  Material cost ($/yr):        {tea.material_cost:>18,.0f}")
    print(f"  Utility cost ($/yr):         {tea.utility_cost:>18,.0f}")

    # -------------------------
    # MSP — includes reagent cost automatically
    # -------------------------
    print("\n=== MSP ===")
    for k, v in msp.items():
        print(f"  {k:<40} {v}")

    # -------------------------
    # Product rate
    # -------------------------
    print("\n=== PRODUCT RATE ===")
    op_days = getattr(tea, "operating_days", 330.0)
    print(f"  backend_oil total (kg/hr):      {oil_stream.F_mass:>10.2f}")
    print(f"  MicrobialOil in stream (kg/hr): {float(oil_stream.imass['MicrobialOil']):>10.2f}")
    print(f"  MicrobialOil annual (kg/yr):    {float(oil_stream.imass['MicrobialOil']) * op_days * 24:>10.0f}")

    # -------------------------
    # Unit purchase costs
    # -------------------------
    print("\n=== TOP UNIT PURCHASE COSTS ===")
    purchase = [
        (u.ID, sum(u.baseline_purchase_costs.values()) if u.baseline_purchase_costs else 0.0)
        for u in full_sys.units
    ]
    purchase.sort(key=lambda x: x[1], reverse=True)
    for uid, cost in purchase[:15]:
        print(f"  {uid:<20} ${cost:>15,.0f}")

    return tea, msp, streams, units, full_sys


# =============================================================
# Sensitivity analysis helpers
# =============================================================

def run_reagent_sensitivity(feed_price: float = 0.00):
    """
    Sensitivity on oil extraction reagent cost at fixed feed price.
    Tests: $0.50, $1.00, $1.50 per kg oil.
    """
    print("\n" + "=" * 60)
    print(f"SENSITIVITY: Oil Extraction Reagent Cost")
    print(f"(feed price = ${feed_price:.3f}/kg wet, solids = ${abs(SOLIDS_DISPOSAL_USD_PER_KG):.2f}/kg)")
    print("=" * 60)
    print(f"  {'Reagent ($/kg oil)':<22} {'MSP ($/kg oil)':<20} {'TCI ($M)':<15} {'VOC ($/yr M)'}")
    print("  " + "-" * 70)

    for reagent in REAGENT_SENSITIVITY:
        tea, msp, *_ = run_case(
            feed_price_per_kg_wet=feed_price,
            reagent_usd_per_kg_oil=reagent,
            silent=True,
        )
        msp_val = msp.get("usd_per_kg_product", float("nan"))
        print(
            f"  ${reagent:<21.2f} ${msp_val:<19.4f} "
            f"${tea.TCI/1e6:<14.1f} ${tea.VOC/1e6:.1f}M"
        )


def run_solids_sensitivity(feed_price: float = 0.00):
    """
    Sensitivity on acidogenic solids disposal cost at fixed feed price.
    Tests: $0.02, $0.04, $0.10 per kg wet cake --> disposal costs

    """
    print("\n" + "=" * 60)
    print(f"SENSITIVITY: Acidogenic Solids Disposal Cost")
    print(f"(feed price = ${feed_price:.3f}/kg wet, reagent = ${OIL_EXTRACTION_REAGENT_USD_PER_KG_OIL:.2f}/kg oil)")
    print("=" * 60)
    print(f"  {'Disposal ($/kg wet)':<24} {'Label':<20} {'MSP ($/kg oil)':<20} {'Disposal cost ($/yr)'}")
    print("  " + "-" * 80)

    labels = {-0.02: "standard biosolids", -0.04: "base case", -0.10: "hazardous waste"}
    for rate in SOLIDS_SENSITIVITY:
        tea, msp, streams, *_ = run_case(
            feed_price_per_kg_wet=feed_price,
            solids_disposal_usd_per_kg=rate,
            silent=True,
        )
        msp_val = msp.get("usd_per_kg_product", float("nan"))
        # compute disposal cost for solids only
        annual_hours = 330.0 * 24.0
        try:
            solids_stream = bst.main_flowsheet.stream["acidogenic_residual_solids"]
            disposal_annual = abs(rate * solids_stream.F_mass * annual_hours)
        except Exception:
            disposal_annual = float("nan")
        label = labels.get(rate, "")
        print(
            f"  ${abs(rate):<23.3f} {label:<20} ${msp_val:<19.4f} ${disposal_annual/1e6:.2f}M"
        )
    
def _apply_biostimulant_credit(
    price_usd_per_kg: float = BIOSTIMULANT_PRICE_BASE_USD_PER_KG,
):
    annual_hours = 330.0 * 24.0
    summary = {}

    for sid in ("biostimulant_membrane_concentrate", "biostimulant_product"):
        try:
            s = bst.main_flowsheet.stream[sid]
            s.price = price_usd_per_kg
            summary[sid] = s.price * s.F_mass * annual_hours
        except Exception:
            pass

    return summary


# =============================================================
# High-value product scenario comparison
# =============================================================

# Product scenarios — same VFA fermentation infrastructure,
# different yield, residence time, extraction cost, and market price.
#
# Modeling scope:
#   - Yield and residence time are passed directly into create_vfa_fermentation_system()
#     so the fermenter (R601) auto-scales its volume and CAPEX correctly.
#   - Extraction cost replaces the oil reagent in OE.add_OPEX ($/kg product).
#   - Downstream purification (HPLC for astaxanthin, winterization for EPA)
#     is lumped into the extraction cost rather than modeled as explicit units.
#   - All three scenarios use the same VFA feedstock and preprocessing.
#   - Engineered Yarrowia strains are assumed available — metabolic engineering
#     costs are not capitalized
#
# Sources:
#   Microbial oil:   Cortés-Peña et al. 2024 GCB Bioenergy; Knoshaug et al. 2018 NREL
#   EPA omega-3:     Xue et al. 2013 Nat Biotechnol (15% DCW); Xie et al. 2017 Biotechnol Bioeng
#   Astaxanthin:     Abdullah et al. 2025 Synth Syst Biotechnol (2.75 g/L, 156h fed-batch)
#   Market prices:   Oil $2-5/kg; EPA bulk $50-150/kg; Astaxanthin natural $1500-2000/kg

PRODUCT_SCENARIOS = [
    {
        "label":          "Microbial oil",
        "product_ID":     "MicrobialOil",
        "yield":           0.144,    # Gao et al. 2020, Y_L/S at pH 8 #0.144
        "residence_h":    48.0,
        "extraction_usd_per_kg": 0.50,
        "market_price":   1.50,    # $/kg — vegetable oil commodity proxy
        "market_label":   "$1.50/kg",
        "global_market_tonne_yr": None,
        "source": "Gao et al. 2020 Biotechnol Biofuels; World Bank commodity price",
    },
    {
        "label":          "Omega-3 EPA oil",
        "product_ID":     "MicrobialOil",  # same chemical class, re-use stream
        "yield":           0.10,    # ~15% DCW on VFA basis — Xue et al. 2013
        "residence_h":   168.0,     # 7-day fed-batch — Xie et al. 2017
        "extraction_usd_per_kg": 50.00,  # winterization + polishing chromatography
        "market_price":   80.00,    # $/kg bulk EPA oil — DuPont Newharvest basis
        "market_label":  "$80/kg",
        "global_market_tonne_yr": 50000,  # fish oil EPA equivalent — no constraint at this scale
        "source": "Xue et al. 2013 Nat Biotechnol; Xie et al. 2017 Biotechnol Bioeng",
    },
    {
        "label":          "Astaxanthin",
        "product_ID":     "MicrobialOil",  # same extraction stream, re-use
        "yield":           0.003,   # g/g VFA — back-calc from 2.75 g/L at ~50 g/L DCW
        "residence_h":   156.0,     # fed-batch — Abdullah et al. 2025
        "extraction_usd_per_kg": 500.00, # cell disruption + HPLC + saponification
        "market_price":  1800.00,   # $/kg natural astaxanthin wholesale
        "market_label":  "$1,800/kg",
        "global_market_tonne_yr": 250,   # global natural astaxanthin — BINDING constraint
        "source": "Abdullah et al. 2025 Synth Syst Biotechnol; market $447M/250t",
    },
]

def build_and_simulate_scenario(
    feed_price_per_kg_wet: float,
    product_yield: float,
    residence_time_h: float,
):
    """
    Build and simulate the VFA fermentation system with custom
    yield and residence time parameters.
    Returns (vfa_sys, fer_sys, streams, units, full_sys).
    """
    bst.main_flowsheet.clear()
    set_thermo()

    vfa_sys = create_vfa_ad_system()
    vfa_sys.feeds[0].price = feed_price_per_kg_wet

    fer_sys, streams, units = create_vfa_fermentation_system(
        vfa_broth=vfa_sys.flowsheet.stream.vfa_broth,
        product_yield_kg_per_kg_vfa_consumed=product_yield,   # scenario override --> from scenario input
        residence_time_h=residence_time_h,                    # scenario override --> from scenario input
    )

    for sid, price in {
        "fermentation_mgso4":     PRICE_MGSO4_USD_PER_KG,
        "fermentation_ammonia":   PRICE_AMMONIA_USD_PER_KG,
        "fermentation_phosphate": PRICE_PHOSPHATE_USD_PER_KG,
        "fermentation_base":      PRICE_NAOH_USD_PER_KG,
    }.items():
        try:
            bst.main_flowsheet.stream[sid].price = price
        except Exception:
            pass

    full_sys = bst.System.from_units(
        "full_vfa_to_product_sys",
        units=list(vfa_sys.units) + list(fer_sys.units),
    )
    full_sys.simulate()

    return vfa_sys, fer_sys, streams, units, full_sys


def run_product_scenario_comparison(feed_price: float = 0.00):
    """
    Compare microbial oil vs EPA vs astaxanthin at the same VFA feedstock
    Reports MSP, TCI, annual production, and market absorption for each scenario
    """
    annual_hours = 330.0 * 24.0

    print("\n" + "=" * 100)
    print(f"PRODUCT SCENARIO COMPARISON — feed = ${feed_price:.3f}/kg wet Sargassum")
    print(
        "Downstream purification lumped into extraction cost. "
        "Engineered strain costs not capitalized."
    )
    print("=" * 100)

    # Header
    print(
        f"  {'Product':<22} {'Yield':>8} {'t_res':>7} {'Extr.':>8} "
        f"{'MSP ($/kg)':>12} {'TCI ($M)':>10} {'Prod. (t/yr)':>13} "
        f"{'Market (t/yr)':>14} {'Mkt absorb':>11} {'NPV @ mkt ($M)':>15}"
    )
    print("  " + "-" * 116)

    for sc in PRODUCT_SCENARIOS:
        try:
            vfa_sys, fer_sys, streams, units, full_sys = build_and_simulate_scenario(
                feed_price_per_kg_wet=feed_price,
                product_yield=sc["yield"],
                residence_time_h=sc["residence_h"],
            )

            # Patch Ev607
            _patch_ev607(full_sys, silent=True)

            # Disposal costs
            _apply_disposal_costs(streams)

            # Wire extraction cost as OE.add_OPEX
            oil_stream = streams["backend_oil"]
            product_kg_hr = float(oil_stream.imass["MicrobialOil"])
            extraction_usd_per_hr = product_kg_hr * sc["extraction_usd_per_kg"]
            try:
                oe_unit = bst.main_flowsheet.unit["OE"]
                oe_unit.add_OPEX = {
                    "Product extraction/purification": extraction_usd_per_hr
                }
            except Exception:
                pass

            tea = make_baseline_tea(full_sys)

            msp_dict = solve_product_msp(
                tea=tea,
                product_stream=oil_stream,
                product_ID="MicrobialOil",
            )
            msp_val = msp_dict.get("usd_per_kg_product", float("nan"))
            annual_product_kg = msp_dict.get("annual_product_kg", 0.0)
            annual_product_t = annual_product_kg / 1000.0

            # NPV at market price
            old_price = oil_stream.price
            oil_mass_frac = float(oil_stream.imass["MicrobialOil"]) / float(oil_stream.F_mass)
            oil_stream.price = sc["market_price"] * oil_mass_frac
            npv_at_market = tea.NPV / 1e6
            oil_stream.price = old_price

            # Market absorption (how much percentage of the global market)
            mkt = sc["global_market_tonne_yr"]
            if mkt is not None:
                absorb_pct = annual_product_t / mkt * 100
                absorb_str = f"{absorb_pct:.0f}% of mkt"
                mkt_str = f"{mkt:,.0f}"
            else:
                absorb_str = "no constraint"
                mkt_str = "n/a"

            print(
                f"  {sc['label']:<22} {sc['yield']:>8.3f} {sc['residence_h']:>6.0f}h "
                f"${sc['extraction_usd_per_kg']:>6.0f}/kg "
                f"${msp_val:>11.2f} ${tea.TCI/1e6:>9.1f}M "
                f"{annual_product_t:>12,.1f}t {mkt_str:>14} "
                f"{absorb_str:>11} ${npv_at_market:>13.1f}M"
            )

        except Exception as e:
            print(f"  {sc['label']:<22} ERROR: {e}")

    print()
    print("  Notes:")
    print("  - MSP = minimum selling price assuming no revenue from other products")
    print("  - NPV @ market = project NPV if product sells at assumed market price")
    print("  - Astaxanthin: global natural market ~250 t/yr (market absorption is binding)")
    print("  - EPA: global market ~50,000 t/yr fish-oil-equivalent (not binding at this scale)")
    print("  - All scenarios use same VFA feedstock, preprocessing, and acidogenic AD")
    print("  - Engineered Yarrowia strains assumed available; strain dev. costs not included")


# =============================================================
# Main
# =============================================================

if __name__ == "__main__":
    # -------------------------
    # Base cases: near_zero and tipping_fee
    # -------------------------
    print("\n>>> BASE CASE: near_zero (feed = $0.00/kg)")
    run_case(0.00, case_label="near_zero", run_diagnostics=True)

    print("\n>>> BASE CASE: tipping_fee (feed = -$0.02/kg)")
    run_case(-0.02, case_label="tipping_fee", run_diagnostics=False)

    # -------------------------
    # Additional feed price scenarios (no diagnostics)
    # -------------------------
    print("\n>>> ADDITIONAL FEED PRICE SCENARIOS")
    for label, price in FEED_PRICE_CASES:
        if price in (0.00, -0.02):
            continue
        print(f"\n  --- {label} (${price:.3f}/kg wet) ---")
        run_case(price, case_label=label, run_diagnostics=False)

    # -------------------------
    # Sensitivity 1: Oil extraction reagent cost
    # -------------------------
    run_reagent_sensitivity(feed_price=0.00)

    # -------------------------
    # Sensitivity 2: Acidogenic solids disposal cost
    # -------------------------
    run_solids_sensitivity(feed_price=0.00)

    # -------------------------
    # Product scenario comparison
    # Microbial oil vs EPA omega-3 vs Astaxanthin
    # -------------------------
    print("\n>>> PRODUCT SCENARIO COMPARISON")
    run_product_scenario_comparison(feed_price=0.00)

    print("\n>>> PRODUCT SCENARIO COMPARISON — tipping_fee")
    run_product_scenario_comparison(feed_price=-0.02)
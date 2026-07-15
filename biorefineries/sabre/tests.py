# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression tests for the SaBRe biorefineries.sabre.analyses scripts.

Baseline values were captured by running the current, pre-migration code
(assumptions.yaml loaded and threaded through as kwargs at the analyses-script
layer). As that data is relocated into unit __init__ defaults and cost
decorators, these tests should keep passing unchanged, since the same
underlying numbers are being moved, not altered. A failure here means a
migration step silently changed a value instead of just relocating it.

Most analyses/*.py scripts originally ran their full simulation sweep
unconditionally at import time (no `if __name__ == "__main__":` guard),
which made them expensive and side-effecting (writing figures/spreadsheets)
just to import their reusable build_case()/run_case() functions. Those
scripts were given a guard as a prerequisite for this test module -- a pure
indentation change, no logic or values touched. Only a handful of spot
points are checked per script (not the full sweeps those scripts print),
since the goal here is drift detection, not re-running every analysis.
"""
import math
import sys
from pathlib import Path

import numpy as np

_ANALYSES_DIR = Path(__file__).resolve().parent / "analyses"
if str(_ANALYSES_DIR) not in sys.path:
    sys.path.insert(0, str(_ANALYSES_DIR))

from vfa_fermentation_tea import (  # noqa: E402
    run_case,
    build_and_simulate_scenario,
    PRODUCT_SCENARIOS,
)
from biorefineries.sabre._tea import make_baseline_tea, solve_product_msp  # noqa: E402

_RTOL = 1e-6
_ATOL = 1e-6


def _assert_close(actual, expected, label):
    assert np.allclose(actual, expected, rtol=_RTOL, atol=_ATOL), (
        f"{label}: expected {expected!r}, got {actual!r}"
    )


def _assert_nan(actual, label):
    assert math.isnan(actual), f"{label}: expected NaN, got {actual!r}"


def _check_base_case(feed_price, expected):
    tea, msp, streams, units, full_sys = run_case(
        feed_price, run_diagnostics=False, silent=True,
    )
    oil_stream = streams["backend_oil"]

    _assert_close(tea.TCI, expected["TCI"], "TCI")
    _assert_close(tea.FCI, expected["FCI"], "FCI")
    _assert_close(tea.VOC, expected["VOC"], "VOC")
    _assert_close(tea.FOC, expected["FOC"], "FOC")
    _assert_close(tea.material_cost, expected["material_cost"], "material_cost")
    _assert_close(tea.utility_cost, expected["utility_cost"], "utility_cost")
    _assert_close(tea.sales, expected["sales"], "sales")
    _assert_close(tea.NPV, expected["NPV"], "NPV")
    _assert_close(msp["usd_per_kg_product"], expected["MSP_usd_per_kg"], "MSP")
    _assert_close(
        float(oil_stream.imass["MicrobialOil"]), expected["oil_kg_hr"], "oil_kg_hr"
    )


def test_vfa_fermentation_near_zero_feed_price():
    _check_base_case(0.00, {
        "TCI": 283650633.1905655,
        "FCI": 270143460.18149096,
        "VOC": 25319059.33262662,
        "FOC": 10805738.40725964,
        "material_cost": 0.0,
        "utility_cost": 18602208.87289611,
        "sales": -21669346.949496437,
        "NPV": -772484936.7882472,
        "MSP_usd_per_kg": 6.9514255681728745,
        "oil_kg_hr": 1685.646383051919,
    })


def test_vfa_fermentation_tipping_fee_feed_price():
    _check_base_case(-0.02, {
        "TCI": 283650633.1905656,
        "FCI": 270143460.181491,
        "VOC": -73680940.6673734,
        "FOC": 10805738.407259641,
        "material_cost": -99000000.00000001,
        "utility_cost": 18602208.87289611,
        "sales": -21669346.949496437,
        "NPV": 45842016.47741586,
        "MSP_usd_per_kg": -0.46412745216541507,
        "oil_kg_hr": 1685.646383051919,
    })


def test_vfa_fermentation_product_scenarios():
    """
    Microbial oil / omega-3 EPA / astaxanthin scenarios from
    vfa_fermentation_tea.PRODUCT_SCENARIOS, at feed price = $0.00/kg wet.
    """
    expected = {
        "Microbial oil": {
            "TCI": 283650633.1905656,
            "VOC": 18643899.655741017,
            "MSP_usd_per_kg": 4.828292327488782,
            "annual_product_kg": 13350319.353771199,
        },
        "Omega-3 EPA oil": {
            "TCI": 577359449.7001549,
            "VOC": 20609437.1094456,
            "MSP_usd_per_kg": 12.28174425161931,
            "annual_product_kg": 9271055.106785554,
        },
        "Astaxanthin": {
            "TCI": 547534273.5079175,
            "VOC": 20264448.373959795,
            "MSP_usd_per_kg": 390.8306206438654,
            "annual_product_kg": 278131.6532035666,
        },
    }

    for sc in PRODUCT_SCENARIOS:
        vfa_sys, fer_sys, streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=0.00,
            product_yield=sc["yield"],
            residence_time_h=sc["residence_h"],
        )
        tea = make_baseline_tea(full_sys)
        oil_stream = streams["backend_oil"]
        msp_dict = solve_product_msp(
            tea=tea, product_stream=oil_stream, product_ID="MicrobialOil",
        )

        exp = expected[sc["label"]]
        label = sc["label"]
        _assert_close(tea.TCI, exp["TCI"], f"{label} TCI")
        _assert_close(tea.VOC, exp["VOC"], f"{label} VOC")
        _assert_close(
            msp_dict["usd_per_kg_product"], exp["MSP_usd_per_kg"], f"{label} MSP"
        )
        _assert_close(
            msp_dict["annual_product_kg"],
            exp["annual_product_kg"],
            f"{label} annual_product_kg",
        )


# =================================================================
# ad_tea_final.py -- methanogenic AD TEA
# =================================================================
import ad_tea_final as _adf  # noqa: E402


def _check_ad_tea_final_case(pretreatment_case, feed_price, expected):
    tea, msp, sys_, npv_results = _adf.run_case(
        pretreatment_case=pretreatment_case,
        feed_price_per_kg_wet=feed_price,
        silent=True,
        target_biomethane_prices_mmbtu=_adf.BIOMETHANE_MARKET_MMBTU,
    )
    _assert_close(tea.TCI, expected["TCI"], "TCI")
    _assert_close(tea.VOC, expected["VOC"], "VOC")
    _assert_close(tea.FOC, expected["FOC"], "FOC")
    _assert_close(tea.material_cost, expected["material_cost"], "material_cost")
    _assert_close(tea.utility_cost, expected["utility_cost"], "utility_cost")
    _assert_close(msp["usd_per_mmbtu"], expected["msp_usd_per_mmbtu"], "MSP $/MMBtu")
    _assert_close(msp["usd_per_kg_ch4"], expected["msp_usd_per_kg_ch4"], "MSP $/kg CH4")
    _assert_close(
        float(sys_.flowsheet.stream.biomethane.F_mass),
        expected["biomethane_kg_hr"],
        "biomethane_kg_hr",
    )
    _assert_close(npv_results[3.0], expected["npv_at_3_mmbtu"], "NPV @ $3/MMBtu")


def test_ad_tea_final_press_mill_only_base_feed_price():
    _check_ad_tea_final_case("press_mill_only", 0.02, {
        "TCI": 234329909.5434592,
        "VOC": 111097636.11626092,
        "FOC": 8926853.696893685,
        "material_cost": 99000000.00000001,
        "utility_cost": 6444913.296881207,
        "msp_usd_per_mmbtu": 77.32617374310271,
        "msp_usd_per_kg_ch4": 4.067356738887203,
        "biomethane_kg_hr": 5880.880607140298,
        "npv_at_3_mmbtu": -1323258108.6818225,
    })


def test_ad_tea_final_press_mill_only_tipping_fee():
    _check_ad_tea_final_case("press_mill_only", -0.02, {
        "TCI": 234329909.5434592,
        "VOC": -86902363.88373911,
        "FOC": 8926853.696893685,
        "material_cost": -99000000.00000001,
        "utility_cost": 6444913.296881208,
        "msp_usd_per_mmbtu": -15.683469674919225,
        "msp_usd_per_kg_ch4": -0.8249505049007513,
        "biomethane_kg_hr": 5880.880607140298,
        "npv_at_3_mmbtu": 278430085.57203114,
    })


def test_ad_tea_final_combined_pe_near_zero():
    _check_ad_tea_final_case("combined_PE", 0.00, {
        "TCI": 260748060.61126694,
        "VOC": 22433166.024744563,
        "FOC": 9933259.451857788,
        "material_cost": 0.0,
        "utility_cost": 7585557.70201066,
        "msp_usd_per_mmbtu": 13.543093261271574,
        "msp_usd_per_kg_ch4": 0.7123667055428848,
        "biomethane_kg_hr": 14395.227832748194,
        "npv_at_3_mmbtu": -515738681.5999446,
    })


# =================================================================
# integrated_tea.py -- alpha sweep across the split biorefinery
# =================================================================
import integrated_tea as _itea  # noqa: E402


def test_integrated_tea_alpha_sweep():
    rows = _itea.run_alpha_sweep(print_summary=False)
    by_alpha = {r["alpha"]: r for r in rows}

    expected = {
        0.0: {
            "msp_oil_usd_per_kg": 11.858154113399467,
            "combined_npv_M": -1210.5891785004417,
            "tci_M": 283.7706558745471,
            "voc_M": 117.56345963121451,
            "oil_kg_yr": 13350319.353771199,
        },
        0.5: {
            "msp_biomethane_mmbtu": 139.65967888996116,
            "msp_oil_usd_per_kg": 22.2698143184909,
            "combined_npv_M": -1096.662612438269,
            "tci_M": 262.9567856842591,
            "voc_M": 114.3342005534255,
            "oil_kg_yr": 6675159.676885599,
        },
        1.0: {
            "msp_biomethane_mmbtu": 64.69794194966727,
            "combined_npv_M": -972.2687745173178,
            "tci_M": 234.32990954345922,
            "voc_M": 111.09763611626092,
            "oil_kg_yr": 0.0,
        },
    }

    _assert_nan(by_alpha[0.0]["msp_biomethane_mmbtu"], "alpha=0 msp_biomethane_mmbtu")
    _assert_nan(by_alpha[1.0]["msp_oil_usd_per_kg"], "alpha=1 msp_oil_usd_per_kg")

    for alpha, exp in expected.items():
        r = by_alpha[alpha]
        for key, val in exp.items():
            _assert_close(r[key], val, f"alpha={alpha} {key}")


# =================================================================
# ad_biostimulant_price.py
# =================================================================
import ad_biostimulant_price as _abp  # noqa: E402


def test_ad_biostimulant_price_build_case():
    msp_mmbtu, annual_rev_musd = _abp.build_case("press_mill_only", 0.00)
    _assert_close(msp_mmbtu, 77.32617374310271, "press_mill_only/$0.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 0.0, "press_mill_only/$0.00 annual_rev_musd")

    msp_mmbtu, annual_rev_musd = _abp.build_case("combined_PE", 1.00)
    _assert_close(msp_mmbtu, 21.213980682219965, "combined_PE/$1.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 53.7662505435227, "combined_PE/$1.00 annual_rev_musd")


# =================================================================
# ad_feed_tea.py
# =================================================================
import ad_feed_tea as _aft  # noqa: E402


def test_ad_feed_tea_build_case():
    msp_mmbtu, fc_musd = _aft.build_case("press_mill_only", 0.02)
    _assert_close(msp_mmbtu, 77.32617374310271, "press_mill_only/0.02 msp_mmbtu")
    _assert_close(fc_musd, 99.00000000000001, "press_mill_only/0.02 annual_feed_cost_musd")

    msp_mmbtu, fc_musd = _aft.build_case("enzymatic", -0.02)
    _assert_close(msp_mmbtu, -12.302964057653245, "enzymatic/-0.02 msp_mmbtu")
    _assert_close(fc_musd, -99.00000000000001, "enzymatic/-0.02 annual_feed_cost_musd")


# =================================================================
# ad_heatmap.py
# =================================================================
import ad_heatmap as _adh  # noqa: E402


def test_ad_heatmap_build_case():
    msp_mmbtu = _adh.build_case("combined_PE", 0.02, 0.50)
    _assert_close(msp_mmbtu, 25.77290912984433, "combined_PE/0.02/0.50 msp_mmbtu")


# =================================================================
# plot_two_pretreatment_figures_from_tea.py
# =================================================================
import plot_two_pretreatment_figures_from_tea as _p2p  # noqa: E402


def test_plot_two_pretreatment_figures_press_mill_only():
    central, low, high = _p2p.get_case_yields("press_mill_only")
    _assert_close(central, 0.1, "press_mill_only central_yield")
    _assert_close(low, 0.1, "press_mill_only low_yield")
    _assert_close(high, 0.1, "press_mill_only high_yield")

    sys_, tea, msp, ch4_kgph = _p2p.build_case("press_mill_only", yield_override=central)
    _assert_close(msp["usd_per_mmbtu"], 77.32617374310271, "msp_usd_per_mmbtu")
    _assert_close(ch4_kgph, 5110.063367928744, "ch4_kgph")
    _assert_close(tea.TCI, 234329909.54345918, "TCI")
    _assert_close(tea.VOC, 111097636.11626092, "VOC")


# =================================================================
# vfa_heatmap.py
# =================================================================
import vfa_heatmap as _vhm  # noqa: E402


def test_vfa_heatmap_build_case():
    msp = _vhm.build_case(0.02, 0.50)
    _assert_close(msp, 12.353309447221273, "feed=0.02/biostim=0.50 msp_usd_per_kg")


# =================================================================
# hauke_stream_tables.py
# =================================================================
import hauke_stream_tables as _hst  # noqa: E402


def test_hauke_stream_tables_press_mill_only():
    rows = _hst.run_case("press_mill_only")
    by_stream = {r["stream"]: r for r in rows}

    expected = {
        "sargassum_feed": {
            "F_mass_kg_hr": 625000.0000000001,
            "TS_kg_hr": 82937.50000000012,
            "VS_kg_hr": 54299.18125000011,
            "water_wt_frac": 0.8672999999999998,
        },
        "biomethane": {
            "F_mass_kg_hr": 5880.880607140298,
            "TS_kg_hr": 5880.880607140298,
            "VS_kg_hr": 5880.880607140298,
            "water_wt_frac": 0.0,
        },
        "soil_amendment": {
            "F_mass_kg_hr": 53856.6173953014,
            "TS_kg_hr": 29621.139567415772,
            "VS_kg_hr": 18731.705245915764,
            "water_wt_frac": 0.44999999999999996,
        },
        "liquid_digestate": {
            "F_mass_kg_hr": 450459.5712948056,
            "TS_kg_hr": 44431.7093511236,
            "VS_kg_hr": 28097.5578688736,
            "water_wt_frac": 0.9013636024573556,
        },
    }

    for stream_id, exp in expected.items():
        row = by_stream[stream_id]
        for key, val in exp.items():
            _assert_close(row[key], val, f"{stream_id} {key}")


# =================================================================
# export_ad_report_fixed.py
# =================================================================
import export_ad_report_fixed as _exr  # noqa: E402


def test_export_ad_report_fixed_press_mill_only():
    sys_, unit_map = _exr.build_case("pelagic_high_quality", "press_mill_only")
    ad = unit_map["AD"]
    up = unit_map["UP"]

    _assert_close(_exr.installed_cost(ad), 140889529.11585945, "AD installed_cost")
    _assert_close(_exr.annualize_add_opex(ad), 4204820.666666668, "AD add_opex_annual")
    _assert_close(_exr.installed_cost(up), 33149424.62270938, "UP installed_cost")
    _assert_close(_exr.stream_total_mass(up.outs[0]), 5880.880607140298, "biomethane F_mass")
    _assert_close(
        _exr.gas_molpct(up.outs[0], "Methane"), 94.78801573491452, "biomethane CH4 mol%"
    )


# =================================================================
# vfa_pathway_figures.py -- VFA composition of the acidogenic broth
# =================================================================
import vfa_pathway_figures as _vpf  # noqa: E402


def test_vfa_pathway_figures_vfa_composition():
    vfa_sys, fer_sys, streams, units, full_sys = _vpf.build_and_simulate(0.02)
    vfa_broth = vfa_sys.flowsheet.stream.vfa_broth

    expected = {
        "AceticAcid": 9198.11406227174,
        "PropionicAcid": 2640.1994067631845,
        "ButyricAcid": 1192.3481191833737,
        "ValericAcid": 1163.9588782504363,
        "HexanoicAcid": 0.0,
    }
    for cid, val in expected.items():
        _assert_close(float(vfa_broth.imass[cid]), val, f"vfa_broth {cid}")


# =================================================================
# vfa_price_scenarios.py -- biostimulant price -> MSP transform
# =================================================================
import vfa_price_scenarios as _vps  # noqa: E402


def test_vfa_price_scenarios_biostimulant_price():
    for bs_price, expected_msp in (
        (0.00, 14.366978528098292),
        (1.00, 10.33964036545701),
    ):
        vfa_sys, fer_sys, streams, units, full_sys = _vps.build_and_simulate(
            _vps.FEED_PRICE_BASE
        )
        _vps._apply_biostimulant_price(streams, bs_price)
        msp_val = _vps._solve_msp_from_system(full_sys, streams)
        _assert_close(msp_val, expected_msp, f"biostim=${bs_price:.2f} msp_usd_per_kg")


# =================================================================
# vfa_yield_sensitivities.py
# =================================================================
import vfa_yield_sensitivities as _vys  # noqa: E402


def test_vfa_yield_sensitivities_acidogenic_yield():
    vfa_sys, fer_sys, streams, units, full_sys = _vys.build_and_simulate(
        _vys.FEED_PRICE_USD_PER_KG_WET
    )
    _vys.override_acidogenic_vfa_yield(vfa_sys, 0.47)
    full_sys.simulate()
    msp_val = _vys.solve_oil_msp_from_full_system(
        full_sys, streams, extraction_usd_per_kg_oil=0.50
    )
    _assert_close(msp_val, 16.65218743436552, "acidogenic_vfa_yield=0.47 msp_usd_per_kg")


def test_vfa_yield_sensitivities_fermentation_yield():
    vfa_sys, fer_sys, streams, units, full_sys = _vys.build_and_simulate_scenario(
        feed_price_per_kg_wet=_vys.FEED_PRICE_USD_PER_KG_WET,
        product_yield=0.20,
        residence_time_h=_vys.BASE_RESIDENCE_H,
    )
    msp_val = _vys.solve_oil_msp_from_full_system(
        full_sys, streams, extraction_usd_per_kg_oil=0.50
    )
    _assert_close(msp_val, 10.490020104025023, "fermentation_yield=0.20 msp_usd_per_kg")


# =================================================================
# integrated_system_sensitivity.py
# =================================================================
import integrated_system_sensitivity as _iss  # noqa: E402


def test_integrated_system_sensitivity_npv_at_alpha():
    npv_M = _iss._build_npv_at_alpha(
        1.0, 0.02, 0.50, "press_mill_only", 3.0, _iss.OIL_MARKET_USD_PER_KG
    )
    _assert_close(npv_M, -1097.4418092934945, "alpha=1.0 base scenario @ $3/MMBtu NPV ($M)")


# =================================================================
# feed_test.py -- Sargassum feed construction from yaml quality params
# (tested directly rather than importing feed_test.py, since that script
# has no reusable functions of its own -- it's a linear diagnostic printout)
# =================================================================
import biosteam as _bst  # noqa: E402
from biorefineries.sabre._chemicals import set_thermo as _set_thermo  # noqa: E402
from biorefineries.sabre._process_settings import (  # noqa: E402
    load_assumptions as _load_assumptions,
    get_quality_params as _get_quality_params,
    get_scale_feed_kgph as _get_scale_feed_kgph,
)
from biorefineries.sabre.streams import make_sargassum_feed as _make_sargassum_feed  # noqa: E402


def test_feed_construction_pelagic_high_quality():
    _bst.main_flowsheet.clear()
    _set_thermo()
    A = _load_assumptions()
    q = _get_quality_params(A, "pelagic_high_quality")
    fresh_feed_kgph = _get_scale_feed_kgph(A)
    feed = _make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph,
        moisture_frac=q["moisture_frac"],
        quality="pelagic_high_quality",
    )
    _assert_close(float(feed.F_mass), 625000.0000000001, "F_mass")
    _assert_close(float(feed.imass["Water"]), 542062.5, "Water")
    _assert_close(float(feed.imass["Ash"]), 28638.31875000001, "Ash")
    _assert_close(q["moisture_frac"], 0.8673, "moisture_frac_yaml")

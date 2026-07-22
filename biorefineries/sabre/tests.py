# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression tests for the SaBRe biorefineries.sabre.legacy_analyses scripts.

Baseline values were captured by running the current, pre-migration code
(assumptions.yaml loaded and threaded through as kwargs at the analyses-script
layer). As that data is relocated into unit __init__ defaults and cost
decorators, these tests should keep passing unchanged, since the same
underlying numbers are being moved, not altered. A failure here means a
migration step silently changed a value instead of just relocating it.

Most legacy_analyses/*.py scripts originally ran their full simulation sweep
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

_LEGACY_ANALYSES_DIR = Path(__file__).resolve().parent / "legacy" / "analyses"
if str(_LEGACY_ANALYSES_DIR) not in sys.path:
    sys.path.insert(0, str(_LEGACY_ANALYSES_DIR))

from vfa_fermentation_tea import (  # noqa: E402
    run_case,
    build_and_simulate_scenario,
    PRODUCT_SCENARIOS,
)
from biorefineries.sabre._tea import create_tea, solve_product_msp  # noqa: E402

_RTOL = 1e-4
_ATOL = 1e-4


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
    oil_stream = streams["microbial_oil"]

    _assert_close(tea.TCI, expected["TCI"], "TCI")
    _assert_close(tea.FCI, expected["FCI"], "FCI")
    _assert_close(tea.VOC, expected["VOC"], "VOC")
    _assert_close(tea.FOC, expected["FOC"], "FOC")
    _assert_close(tea.material_cost, expected["material_cost"], "material_cost")
    _assert_close(tea.utility_cost, expected["utility_cost"], "utility_cost")
    _assert_close(tea.sales, expected["sales"], "sales")
    _assert_close(tea.NPV, expected["NPV"], "NPV")
    _assert_close(msp["usd_per_kg"], expected["MSP_usd_per_kg"], "MSP")
    _assert_close(
        float(oil_stream.imass["MicrobialOil"]), expected["oil_kg_hr"], "oil_kg_hr"
    )


def test_vfa_fermentation_near_zero_feed_price():
    _check_base_case(0.00, {
        "TCI": 321913115.0113817,
        "FCI": 306583919.05845875,
        "VOC": 28510808.06602466,
        "FOC": 12263356.76233835,
        "material_cost": 657159.5976143774,
        "utility_cost": 21335880.031990513,
        "sales": -19927307.16929145,
        "NPV": -835623034.1111149,
        "MSP_usd_per_kg": 8.665822837953183,
        "oil_kg_hr": 1634.1526795831005,
    })


def test_vfa_fermentation_tipping_fee_feed_price():
    _check_base_case(-0.02, {
        "TCI": 321913115.0113817,
        "FCI": 306583919.05845875,
        "VOC": -70489191.93397532,
        "FOC": 12263356.76233835,
        "material_cost": -98342840.4023856,
        "utility_cost": 21335880.031990513,
        "sales": -19927307.16929145,
        "NPV": -10754525.061558746,
        "MSP_usd_per_kg": 1.109424881023777,
        "oil_kg_hr": 1634.1526795831005,
    })


def test_vfa_fermentation_product_scenarios():
    """
    Microbial oil / omega-3 EPA / astaxanthin scenarios from
    vfa_fermentation_tea.PRODUCT_SCENARIOS, at feed price = $0.00/kg wet.
    """
    expected = {
        "Microbial oil": {
            "TCI": 321913115.0113816,
            "VOC": 28576174.173207987,
            "MSP_usd_per_kg": 3.892999084169563,
            "annual_product_kg": 13101480.584426785,
        },
        "Omega-3 EPA oil": {
            "TCI": 500596243.5973739,
            "VOC": 36278993.85905999,
            "MSP_usd_per_kg": 9.575235479205613,
            "annual_product_kg": 9146892.39595412,
        },
        "Astaxanthin": {
            "TCI": 482052584.6274711,
            "VOC": 30765243.58078747,
            "MSP_usd_per_kg": 184.48076779287175,
            "annual_product_kg": 428822.98045753106,
        },
    }

    for sc in PRODUCT_SCENARIOS:
        streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=0.00,
            product_yield=sc["yield"],
            residence_time_h=sc["residence_h"],
        )
        tea = create_tea(full_sys)
        oil_stream = streams["microbial_oil"]
        msp_dict = solve_product_msp(tea=tea, product_stream=oil_stream)

        exp = expected[sc["label"]]
        label = sc["label"]
        _assert_close(tea.TCI, exp["TCI"], f"{label} TCI")
        _assert_close(tea.VOC, exp["VOC"], f"{label} VOC")
        _assert_close(
            msp_dict["usd_per_kg"], exp["MSP_usd_per_kg"], f"{label} MSP"
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
        "TCI": 314254303.7698967,
        "VOC": 110356909.8710158,
        "FOC": 11971592.524567492,
        "material_cost": 98999999.99999999,
        "utility_cost": 6109050.465734033,
        "msp_usd_per_mmbtu": 85.51925191724928,
        "msp_usd_per_kg_ch4": 4.498312650847312,
        "biomethane_kg_hr": 5055.060612937506,
        "npv_at_3_mmbtu": -1422699279.9872646,
    })


def test_ad_tea_final_press_mill_only_tipping_fee():
    _check_ad_tea_final_case("press_mill_only", -0.02, {
        "TCI": 314254303.76989675,
        "VOC": -87643090.12898417,
        "FOC": 11971592.524567494,
        "material_cost": -98999999.99999999,
        "utility_cost": 6109050.465734033,
        "msp_usd_per_mmbtu": -10.421214236518416,
        "msp_usd_per_kg_ch4": -0.5481558688408688,
        "biomethane_kg_hr": 5055.060612937506,
        "npv_at_3_mmbtu": 199805212.71141493,
    })


def test_ad_tea_final_combined_pe_near_zero():
    _check_ad_tea_final_case("combined_PE", 0.00, {
        "TCI": 340278613.448074,
        "VOC": 21769866.467900738,
        "FOC": 12962994.798021866,
        "material_cost": 0.0,
        "utility_cost": 7217818.726045765,
        "msp_usd_per_mmbtu": 16.054388865578613,
        "msp_usd_per_kg_ch4": 0.8444608543294351,
        "biomethane_kg_hr": 14002.517897836899,
        "npv_at_3_mmbtu": -618397463.4113625,
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
            "msp_oil_usd_per_kg": 13.303083384165571,
            "combined_npv_M": -1347.29595876633,
            "tci_M": 321.9131150113816,
            "voc_M": 127.51080806602465,
            "oil_kg_yr": 13101480.584426785,
        },
        0.5: {
            "msp_biomethane_mmbtu": 153.26724546305803,
            "msp_oil_usd_per_kg": 24.14306146548551,
            "combined_npv_M": -1179.4581479923638,
            "tci_M": 321.25433966291683,
            "voc_M": 118.96126016913509,
            "oil_kg_yr": 6550740.292213392,
        },
        1.0: {
            "msp_biomethane_mmbtu": 68.29251657424578,
            "combined_npv_M": -1002.7142670493796,
            "tci_M": 314.2543037698967,
            "voc_M": 110.3569098710158,
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
    _assert_close(msp_mmbtu, 85.51925191724928, "press_mill_only/$0.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 0.0, "press_mill_only/$0.00 annual_rev_musd")

    msp_mmbtu, annual_rev_musd = _abp.build_case("combined_PE", 1.00)
    _assert_close(msp_mmbtu, 20.934086433423587, "combined_PE/$1.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 71.10437847149987, "combined_PE/$1.00 annual_rev_musd")


# =================================================================
# ad_feed_tea.py
# =================================================================
import ad_feed_tea as _aft  # noqa: E402


def test_ad_feed_tea_build_case():
    msp_mmbtu, fc_musd = _aft.build_case("press_mill_only", 0.02)
    _assert_close(msp_mmbtu, 85.51925191724928, "press_mill_only/0.02 msp_mmbtu")
    _assert_close(fc_musd, 98.99999999999999, "press_mill_only/0.02 annual_feed_cost_musd")

    msp_mmbtu, fc_musd = _aft.build_case("enzymatic", -0.02)
    _assert_close(msp_mmbtu, -7.885754437292184, "enzymatic/-0.02 msp_mmbtu")
    _assert_close(fc_musd, -98.99999999999999, "enzymatic/-0.02 annual_feed_cost_musd")


# =================================================================
# ad_heatmap.py
# =================================================================
import ad_heatmap as _adh  # noqa: E402


def test_ad_heatmap_build_case():
    msp_mmbtu = _adh.build_case("combined_PE", 0.02, 0.50)
    _assert_close(msp_mmbtu, 27.153124458522804, "combined_PE/0.02/0.50 msp_mmbtu")


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
    _assert_close(msp["usd_per_mmbtu"], 85.51925191724928, "msp_usd_per_mmbtu")
    _assert_close(ch4_kgph, 4953.959400678757, "ch4_kgph")
    _assert_close(tea.TCI, 314254303.76989675, "TCI")
    _assert_close(tea.VOC, 110356909.8710158, "VOC")


# =================================================================
# vfa_heatmap.py
# =================================================================
import vfa_heatmap as _vhm  # noqa: E402


def test_vfa_heatmap_build_case():
    msp = _vhm.build_case(0.02, 0.50)
    _assert_close(msp, 13.508619843460652, "feed=0.02/biostim=0.50 msp_usd_per_kg")


# =================================================================
# hauke_stream_tables.py
# =================================================================
import hauke_stream_tables as _hst  # noqa: E402


def test_hauke_stream_tables_press_mill_only():
    rows = _hst.run_case("press_mill_only")
    by_stream = {r["stream"]: r for r in rows}

    expected = {
        "sargassum_feed": {
            "F_mass_kg_hr": 624999.9999999999,
            "TS_kg_hr": 82937.49999999988,
            "VS_kg_hr": 54299.18124999988,
            "water_wt_frac": 0.8673000000000002,
        },
        "biomethane": {
            "F_mass_kg_hr": 5055.060612937506,
            "TS_kg_hr": 5055.060612937506,
            "VS_kg_hr": 5055.060612937506,
            "water_wt_frac": 0.0,
        },
        "methanogenic_solid_digestate": {
            "F_mass_kg_hr": 52927.759659600015,
            "TS_kg_hr": 29110.26781278001,
            "VS_kg_hr": 18220.833491280006,
            "water_wt_frac": 0.44999999999999996,
        },
        "liquid_digestate": {
            "F_mass_kg_hr": 441526.6414433647,
            "TS_kg_hr": 43665.40171916998,
            "VS_kg_hr": 27331.250236919976,
            "water_wt_frac": 0.9011035855584469,
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
    sys_, unit_map = _exr.build_case("pelagic", "press_mill_only")
    ad = unit_map["AD"]
    up = unit_map["UP"]

    _assert_close(_exr.installed_cost(ad), 127825506.10408078, "AD installed_cost")
    _assert_close(_exr.annualize_add_opex(ad), 3838081.637052034, "AD add_opex_annual")
    _assert_close(_exr.installed_cost(up), 32136764.637289904, "UP installed_cost")
    _assert_close(_exr.stream_total_mass(up.outs[0]), 5055.060612937506, "biomethane F_mass")
    _assert_close(
        _exr.gas_molpct(up.outs[0], "Methane"), 99.26147149151141, "biomethane CH4 mol%"
    )


# =================================================================
# vfa_pathway_figures.py -- VFA composition of the vfa_broth stream
# (post-microfiltration permeate, ready for fermentation)
# =================================================================
import vfa_pathway_figures as _vpf  # noqa: E402


def test_vfa_pathway_figures_vfa_composition():
    streams, units, full_sys = _vpf.build_and_simulate(0.02)
    vfa_broth = streams["vfa_broth"]

    expected = {
        "AceticAcid": 8738.78438279733,
        "PropionicAcid": 2508.354776543678,
        "ButyricAcid": 1132.8053829552096,
        "ValericAcid": 1105.8338262181805,
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
        (0.00, 16.2222207267543),
        (1.00, 10.795018958056952),
    ):
        streams, units, full_sys = _vps.build_and_simulate(
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
    streams, units, full_sys = _vys.build_and_simulate(
        _vys.FEED_PRICE_USD_PER_KG_WET
    )
    _vys.override_acidogenic_vfa_yield(full_sys, 0.47)
    msp_val = _vys.solve_oil_msp_from_full_system(
        full_sys, streams, extraction_usd_per_kg_oil=0.50
    )
    _assert_close(msp_val, 12.554239982667346, "acidogenic_vfa_yield=0.47 msp_usd_per_kg")


def test_vfa_yield_sensitivities_fermentation_yield():
    streams, units, full_sys = _vys.build_and_simulate_scenario(
        feed_price_per_kg_wet=_vys.FEED_PRICE_USD_PER_KG_WET,
        product_yield=0.20,
        residence_time_h=_vys.BASE_RESIDENCE_H,
    )
    msp_val = _vys.solve_oil_msp_from_full_system(
        full_sys, streams, extraction_usd_per_kg_oil=0.50
    )
    _assert_close(msp_val, 7.943572290878948, "fermentation_yield=0.20 msp_usd_per_kg")


# =================================================================
# integrated_system_sensitivity.py
# =================================================================
import integrated_system_sensitivity as _iss  # noqa: E402


def test_integrated_system_sensitivity_npv_at_alpha():
    npv_M = _iss._build_npv_at_alpha(
        1.0, 0.02, 0.50, "press_mill_only", 3.0, _iss.OIL_MARKET_USD_PER_KG
    )
    _assert_close(npv_M, -1124.063473055855, "alpha=1.0 base scenario @ $3/MMBtu NPV ($M)")


# =================================================================
# feed_test.py -- Sargassum feed construction from yaml quality params
# (tested directly rather than importing feed_test.py, since that script
# has no reusable functions of its own -- it's a linear diagnostic printout)
# =================================================================
import biosteam as _bst  # noqa: E402
from biorefineries.sabre._chemicals import create_chemicals as _create_chemicals  # noqa: E402
from biorefineries.sabre.utils import (  # noqa: E402
    load_assumptions as _load_assumptions,
    get_feedstock_type_params as _get_feedstock_type_params,
    get_scale_feed_kgph as _get_scale_feed_kgph,
    make_sargassum_feed as _make_sargassum_feed,
)


def test_feed_construction_pelagic():
    _bst.main_flowsheet.clear()
    _create_chemicals()
    feedstock_assumptions = _load_assumptions("feedstock.yaml")
    params = _get_feedstock_type_params(feedstock_assumptions, "pelagic")
    fresh_feed_kgph = _get_scale_feed_kgph(feedstock_assumptions)
    feed = _make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph,
        moisture_frac=params["moisture_frac"],
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )
    _assert_close(float(feed.F_mass), 625000.0000000001, "F_mass")
    _assert_close(float(feed.imass["Water"]), 542062.5, "Water")
    _assert_close(float(feed.imass["Ash"]), 28638.31875000001, "Ash")
    _assert_close(params["moisture_frac"], 0.8673, "moisture_frac_yaml")

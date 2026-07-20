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
        "TCI": 284703875.28279954,
        "FCI": 271146547.8883805,
        "VOC": 26141205.5724501,
        "FOC": 10845861.91553522,
        "material_cost": 672247.8455353046,
        "utility_cost": 18744087.558802776,
        "sales": -30892391.26741339,
        "NPV": -858266589.1805022,
        "MSP_usd_per_kg": 7.71659906991003,
        "oil_kg_hr": 1685.646383051923,
    })


def test_vfa_fermentation_tipping_fee_feed_price():
    _check_base_case(-0.02, {
        "TCI": 284703875.28279954,
        "FCI": 271146547.8883805,
        "VOC": -72858794.42754988,
        "FOC": 10845861.91553522,
        "material_cost": -98327752.15446469,
        "utility_cost": 18744087.558802776,
        "sales": -30892391.26741339,
        "NPV": -30313167.11144813,
        "MSP_usd_per_kg": 0.30104604652366995,
        "oil_kg_hr": 1685.646383051923,
    })


def test_vfa_fermentation_product_scenarios():
    """
    Microbial oil / omega-3 EPA / astaxanthin scenarios from
    vfa_fermentation_tea.PRODUCT_SCENARIOS, at feed price = $0.00/kg wet.
    """
    expected = {
        "Microbial oil": {
            "TCI": 284703875.2827994,
            "VOC": 19466045.89556448,
            "MSP_usd_per_kg": 4.902617610182285,
            "annual_product_kg": 13350319.35377123,
        },
        "Omega-3 EPA oil": {
            "TCI": 578894201.1955233,
            "VOC": 21434131.76771388,
            "MSP_usd_per_kg": 12.39743639530954,
            "annual_product_kg": 9271055.106785577,
        },
        "Astaxanthin": {
            "TCI": 549048352.6077971,
            "VOC": 21088994.01350955,
            "MSP_usd_per_kg": 394.6744844983792,
            "annual_product_kg": 278131.6532035673,
        },
    }

    for sc in PRODUCT_SCENARIOS:
        streams, units, full_sys = build_and_simulate_scenario(
            feed_price_per_kg_wet=0.00,
            product_yield=sc["yield"],
            residence_time_h=sc["residence_h"],
        )
        tea = create_tea(full_sys)
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
        "TCI": 234607314.6014704,
        "VOC": 111337839.81391096,
        "FOC": 8937421.508627443,
        "material_cost": 98999999.99999999,
        "utility_cost": 6666574.558810216,
        "msp_usd_per_mmbtu": 77.46005777678312,
        "msp_usd_per_kg_ch4": 4.074399039058792,
        "biomethane_kg_hr": 5880.880607140312,
        "npv_at_3_mmbtu": -1325645300.12332,
    })


def test_ad_tea_final_press_mill_only_tipping_fee():
    _check_ad_tea_final_case("press_mill_only", -0.02, {
        "TCI": 234607314.6014704,
        "VOC": -86662160.18608901,
        "FOC": 8937421.508627443,
        "material_cost": -98999999.99999999,
        "utility_cost": 6666574.558810216,
        "msp_usd_per_mmbtu": -15.549585661015861,
        "msp_usd_per_kg_ch4": -0.8179082057694343,
        "biomethane_kg_hr": 5880.880607140312,
        "npv_at_3_mmbtu": 276550669.48272127,
    })


def test_ad_tea_final_combined_pe_near_zero():
    _check_ad_tea_final_case("combined_PE", 0.00, {
        "TCI": 261063483.58637446,
        "VOC": 22976566.46277774,
        "FOC": 9945275.565195218,
        "material_cost": 0.0,
        "utility_cost": 7809748.798318315,
        "msp_usd_per_mmbtu": 13.643885758613232,
        "msp_usd_per_kg_ch4": 0.717668390903056,
        "biomethane_kg_hr": 14395.227832748227,
        "npv_at_3_mmbtu": -520723361.3442967,
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
            "msp_oil_usd_per_kg": 12.980795778504042,
            "combined_npv_M": -1336.4609103990551,
            "tci_M": 284.70387528279946,
            "voc_M": 117.79379805002915,
            "oil_kg_yr": 13350319.35377123,
        },
        0.5: {
            "msp_biomethane_mmbtu": 149.42187802216412,
            "msp_oil_usd_per_kg": 23.826472358536844,
            "combined_npv_M": -1183.930695721619,
            "tci_M": 263.56210861162157,
            "voc_M": 114.5694716231035,
            "oil_kg_yr": 6675159.676885615,
        },
        1.0: {
            "msp_biomethane_mmbtu": 67.41976757012725,
            "combined_npv_M": -1020.9331805064713,
            "tci_M": 234.60731460147042,
            "voc_M": 111.33783981391096,
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
    _assert_close(msp_mmbtu, 77.46005777678312, "press_mill_only/$0.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 0.0, "press_mill_only/$0.00 annual_rev_musd")

    msp_mmbtu, annual_rev_musd = _abp.build_case("combined_PE", 1.00)
    _assert_close(msp_mmbtu, 23.183323062463288, "combined_PE/$1.00 msp_mmbtu")
    _assert_close(annual_rev_musd, 42.74777083871248, "combined_PE/$1.00 annual_rev_musd")


# =================================================================
# ad_feed_tea.py
# =================================================================
import ad_feed_tea as _aft  # noqa: E402


def test_ad_feed_tea_build_case():
    msp_mmbtu, fc_musd = _aft.build_case("press_mill_only", 0.02)
    _assert_close(msp_mmbtu, 77.46005777678312, "press_mill_only/0.02 msp_mmbtu")
    _assert_close(fc_musd, 99.00000000000001, "press_mill_only/0.02 annual_feed_cost_musd")

    msp_mmbtu, fc_musd = _aft.build_case("enzymatic", -0.02)
    _assert_close(msp_mmbtu, -12.089341350460483, "enzymatic/-0.02 msp_mmbtu")
    _assert_close(fc_musd, -99.00000000000001, "enzymatic/-0.02 annual_feed_cost_musd")


# =================================================================
# ad_heatmap.py
# =================================================================
import ad_heatmap as _adh  # noqa: E402


def test_ad_heatmap_build_case():
    msp_mmbtu = _adh.build_case("combined_PE", 0.02, 0.50)
    _assert_close(msp_mmbtu, 26.807976568008115, "combined_PE/0.02/0.50 msp_mmbtu")


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
    _assert_close(msp["usd_per_mmbtu"], 77.46005777678312, "msp_usd_per_mmbtu")
    _assert_close(ch4_kgph, 5110.063367928756, "ch4_kgph")
    _assert_close(tea.TCI, 234607314.6014704, "TCI")
    _assert_close(tea.VOC, 111337839.81391096, "VOC")


# =================================================================
# vfa_heatmap.py
# =================================================================
import vfa_heatmap as _vhm  # noqa: E402


def test_vfa_heatmap_build_case():
    msp = _vhm.build_case(0.02, 0.50)
    _assert_close(msp, 13.53115021844506, "feed=0.02/biostim=0.50 msp_usd_per_kg")


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
    sys_, unit_map = _exr.build_case("pelagic", "press_mill_only")
    ad = unit_map["AD"]
    up = unit_map["UP"]

    _assert_close(_exr.installed_cost(ad), 140924589.88401794, "AD installed_cost")
    _assert_close(_exr.annualize_add_opex(ad), 4215343.394006237, "AD add_opex_annual")
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
    streams, units, full_sys = _vpf.build_and_simulate(0.02)
    vfa_broth = streams["vfa_broth"]

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
        (0.00, 15.132152027791374),
        (1.00, 11.93014840867312),
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
    _assert_close(msp_val, 17.628029213265773, "acidogenic_vfa_yield=0.47 msp_usd_per_kg")


def test_vfa_yield_sensitivities_fermentation_yield():
    streams, units, full_sys = _vys.build_and_simulate_scenario(
        feed_price_per_kg_wet=_vys.FEED_PRICE_USD_PER_KG_WET,
        product_yield=0.20,
        residence_time_h=_vys.BASE_RESIDENCE_H,
    )
    msp_val = _vys.solve_oil_msp_from_full_system(
        full_sys, streams, extraction_usd_per_kg_oil=0.50
    )
    _assert_close(msp_val, 11.041010251377745, "fermentation_yield=0.20 msp_usd_per_kg")


# =================================================================
# integrated_system_sensitivity.py
# =================================================================
import integrated_system_sensitivity as _iss  # noqa: E402


def test_integrated_system_sensitivity_npv_at_alpha():
    npv_M = _iss._build_npv_at_alpha(
        1.0, 0.02, 0.50, "press_mill_only", 3.0, _iss.OIL_MARKET_USD_PER_KG
    )
    _assert_close(npv_M, -1146.1062152826482, "alpha=1.0 base scenario @ $3/MMBtu NPV ($M)")


# =================================================================
# feed_test.py -- Sargassum feed construction from yaml quality params
# (tested directly rather than importing feed_test.py, since that script
# has no reusable functions of its own -- it's a linear diagnostic printout)
# =================================================================
import biosteam as _bst  # noqa: E402
from biorefineries.sabre._chemicals import set_thermo as _set_thermo  # noqa: E402
from biorefineries.sabre.utils import (  # noqa: E402
    load_assumptions as _load_assumptions,
    get_feedstock_type_params as _get_feedstock_type_params,
    get_scale_feed_kgph as _get_scale_feed_kgph,
    make_sargassum_feed as _make_sargassum_feed,
)


def test_feed_construction_pelagic():
    _bst.main_flowsheet.clear()
    _set_thermo()
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

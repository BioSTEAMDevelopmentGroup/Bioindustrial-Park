# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Regression tests for the current SaBRe systems (biostimulant, biomethane, VFA,
microbial oil, and 3-HP), via each system's `price_*_system()` function, which
builds, simulates, and solves the minimum selling price (MSP) of the product.

Golden values were recorded from the current code. The tolerance is loose
because results move by several percent with upstream biosteam/thermosteam/
hensmith updates (e.g., heat exchanger network results), so these tests guard
against gross regressions rather than exact reproduction.

The 3-HP system rebuilds the process-global chemical set (see
`systems/_3hp_system.py`), so it is run in a separate Python process.
"""
import json
import subprocess
import sys

import numpy as np
import pytest

from biorefineries.sabre.systems import (
    price_biostimulant_system,
    price_biomethane_system,
    price_vfa_system,
    price_microbial_oil_system,
)

_RTOL = 1e-1


def _assert_close(actual, expected, label):
    assert np.allclose(actual, expected, rtol=_RTOL), (
        f"{label}: expected {expected!r}, got {actual!r}"
    )


def test_biostimulant():
    r = price_biostimulant_system()
    _assert_close(r["msp_usd_per_kg"], 2.7226, "MSP")
    _assert_close(r["annual_product_kg"], 71104378.5, "annual product")


def test_biomethane_press_mill_only():
    r = price_biomethane_system(pretreatment_case="press_mill_only")
    _assert_close(r["msp_usd_per_kg"], 0.6109, "MSP $/kg")
    _assert_close(r["msp_usd_per_mmbtu"], 11.6145, "MSP $/MMBtu")
    _assert_close(r["annual_product_kg"], 40036080.1, "annual product")


def test_vfa():
    r = price_vfa_system()
    _assert_close(r["msp_usd_per_kg"], 0.10329, "MSP")
    _assert_close(r["annual_product_kg"], 104671217.4, "annual product")


def test_microbial_oil():
    r = price_microbial_oil_system()
    _assert_close(r["msp_usd_per_kg"], 3.8144, "MSP")
    _assert_close(r["annual_product_kg"], 13101480.6, "annual product")


_3HP_SCRIPT = """
import json, sys
from biorefineries.sabre.systems import price_3hp_system
r = price_3hp_system(scenario=sys.argv[1])
print("RESULT" + json.dumps({k: r[k] for k in ("msp_usd_per_kg", "annual_product_kg")}))
"""


@pytest.mark.parametrize("scenario, msp, annual", [
    ("baseline", 1.3769, 75215630.6),
    ("acid_press", 1.149, 75497210.2),
])
def test_3hp(scenario, msp, annual):
    proc = subprocess.run(
        [sys.executable, "-c", _3HP_SCRIPT, scenario],
        capture_output=True, text=True, check=True,
    )
    line = next(l for l in proc.stdout.splitlines() if l.startswith("RESULT"))
    r = json.loads(line[len("RESULT"):])
    _assert_close(r["msp_usd_per_kg"], msp, f"3-HP MSP ({scenario})")
    _assert_close(r["annual_product_kg"], annual, f"3-HP annual product ({scenario})")


_3HP_INVENTORY_SCRIPT = """
import json, sys
from biorefineries.sabre import _lca
from biorefineries.sabre.systems import price_3hp_system
sys = price_3hp_system(scenario=sys.argv[1])["sys"]
_lca.set_characterization_factors(sys)
sys.simulate()
rows = _lca.inventory_check(sys, sys.flowsheet.stream.ca3hp2_product)
print("RESULT" + json.dumps(dict(
    missing=[r["stream"] for r in rows if r["status"] != "ok"], n=len(rows),
    gwp=_lca.get_gwp(sys, sys.flowsheet.stream.ca3hp2_product)["gwp_no_credit_per_kg_3hp"],
)))
"""


@pytest.mark.parametrize("scenario, gwp", [
    ("baseline", 1.3561),
    ("acid_press", 1.8987),
])
def test_3hp_inventory_complete(scenario, gwp):
    """Every inlet and outlet of each 3-HP scenario has a cost and an impact basis."""
    proc = subprocess.run(
        [sys.executable, "-c", _3HP_INVENTORY_SCRIPT, scenario],
        capture_output=True, text=True, check=True,
    )
    line = next(l for l in proc.stdout.splitlines() if l.startswith("RESULT"))
    r = json.loads(line[len("RESULT"):])
    assert r["missing"] == [], r["missing"]
    assert r["n"] > 20
    _assert_close(r["gwp"], gwp, "GWP100 without electricity credit, per kg 3-HP")

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
r = price_3hp_system(configuration=sys.argv[1])
print("RESULT" + json.dumps({k: r[k] for k in ("msp_usd_per_kg", "annual_product_kg")}))
"""


@pytest.mark.parametrize("configuration, msp, annual_kg", [
    ("enzymatic_press", 0.8823, 49136479.8),
    ("biostimulant_cake", -1.6780, 67836637.4),
    ("juice_cake", 1.1562, 72930492.9),
])
def test_3hp(configuration, msp, annual_kg):
    proc = subprocess.run(
        [sys.executable, "-c", _3HP_SCRIPT, configuration],
        capture_output=True, text=True, check=True,
    )
    line = next(l for l in proc.stdout.splitlines() if l.startswith("RESULT"))
    r = json.loads(line[len("RESULT"):])
    _assert_close(r["msp_usd_per_kg"], msp, f"{configuration} MSP")
    _assert_close(r["annual_product_kg"], annual_kg, f"{configuration} annual product")

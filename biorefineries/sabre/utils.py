# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Process settings (scenario assumptions) for the SaBRe flowsheets.
"""
from pathlib import Path
import yaml

__all__ = ('load_assumptions', 'wet_tpd_to_kgph', 'get_scale_feed_kgph', 'get_quality_params')


def load_assumptions():
    path = Path(__file__).resolve().parent / "data" / "assumptions.yaml"
    with open(path, "r") as f:
        return yaml.safe_load(f)


def wet_tpd_to_kgph(wet_tpd: float, ton_definition: str) -> float:
    if ton_definition == "metric":
        kg_per_ton = 1000.0
    elif ton_definition == "short":
        kg_per_ton = 907.185
    else:
        raise ValueError("ton_definition must be 'metric' or 'short'")
    return wet_tpd * kg_per_ton / 24.0


def get_scale_feed_kgph(A: dict) -> float:
    wet_tpd = A["scale"]["wet_feed_ton_per_day"]
    ton_def = A["scale"]["ton_definition"]
    return wet_tpd_to_kgph(wet_tpd, ton_def)


def get_quality_params(A: dict, quality: str) -> dict:
    qb = A["quality_bins"]
    if quality not in qb:
        raise KeyError(f"Unknown quality '{quality}'. Options: {list(qb.keys())}")
    return qb[quality]

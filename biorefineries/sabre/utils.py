# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Process settings (scenario assumptions) for the SaBRe flowsheets, and
feedstock stream construction from those assumptions.
"""
from pathlib import Path
import yaml
import biosteam as bst

__all__ = (
    'load_assumptions', 'wet_tpd_to_kgph', 'get_scale_feed_kgph',
    'get_feedstock_type_params', 'get_ad_temperature_K', 'make_sargassum_feed',
)


def load_assumptions(filename: str = "assumptions.yaml"):
    # TODO: drop the default once all assumptions have migrated out of
    # assumptions.yaml into topic-specific files (e.g. feedstock.yaml)
    # and every caller passes filename explicitly.
    path = Path(__file__).resolve().parent / "data" / filename
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


def get_feedstock_type_params(A: dict, feedstock_type: str) -> dict:
    ft = A["feedstock_type"]
    if feedstock_type not in ft:
        raise KeyError(f"Unknown feedstock_type '{feedstock_type}'. Options: {list(ft.keys())}")
    return ft[feedstock_type]


def get_ad_temperature_K(ad_assumptions: dict, temperature_regime: str) -> float:
    """
    Resolve the AD operating temperature (K) for a given regime
    ('mesophilic' or 'thermophilic') from data/ad.yaml's
    `ad.temperature_regimes` section.
    """
    regimes = ad_assumptions["temperature_regimes"]
    if temperature_regime not in regimes:
        raise KeyError(
            f"Unknown temperature_regime '{temperature_regime}'. "
            f"Options: {list(regimes.keys())}"
        )
    return float(regimes[temperature_regime]["temperature_K"])


# Dry-basis composition (data/feedstock.yaml `dry_composition` section).
# Structural ratios among non-water components; per-feedstock-type moisture
# and ash targets (`feedstock_type` section) scale/override this baseline
# in make_sargassum_feed() below.
DRY_COMPOSITION = load_assumptions("feedstock.yaml")["dry_composition"]


def make_sargassum_feed(fresh_feed_kgph: float, moisture_frac: float, ash_wt_frac_dry: float):
    water_kgph = fresh_feed_kgph * moisture_frac
    dry_kgph   = fresh_feed_kgph * (1 - moisture_frac)

    base_ash = DRY_COMPOSITION["Ash"]
    base_nonash_sum = 1.0 - base_ash
    target_nonash_sum = 1.0 - ash_wt_frac_dry

    if base_nonash_sum <= 1e-12:
        raise RuntimeError("Base non-ash fraction is ~0; cannot scale.")
    scale = target_nonash_sum / base_nonash_sum

    # Build final dry-basis fractions with ash overridden, other components
    # rescaled to fill the remaining dry mass proportionally.
    dry_fracs = {"Ash": ash_wt_frac_dry}
    for k, v in DRY_COMPOSITION.items():
        if k == "Ash":
            continue
        dry_fracs[k] = v * scale

    # Normalize for numerical safety
    s = sum(dry_fracs.values())
    dry_fracs = {k: v / s for k, v in dry_fracs.items()}

    # Allocate component mass flows
    kwargs = {k: dry_kgph * frac for k, frac in dry_fracs.items()}

    return bst.Stream(
        "sargassum_feed",
        Water=water_kgph,
        units="kg/hr",
        phase="l",
        **kwargs,
    )

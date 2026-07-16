# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Stream functions (feed/product streams)

Purpose:
- Create BioSTEAM Stream objects with correct mass flowrates and composition
- Convert scenario inputs (fresh feed kg/hr, moisture, quality bin) into component flows
"""

import biosteam as bst

__all__ = ('make_sargassum_feed',)

# Wet-basis composition
WET_COMPOSITION = dict(
    Water=0.8673,
    Ash=0.04580,        # 34.53% dry × 0.1327
    OtherSolids=0.03902,  # 29.41% dry × 0.1327
    Glucan=0.01460,     # 11.00% dry × 0.1327
    Alginate=0.01393,   # 10.50% dry × 0.1327
    Protein=0.01003,    # 7.56% dry × 0.1327
    Mannitol=0.00265,   # 2.00% dry × 0.1327
    Fucoidan=0.00663,   # 5.00% dry × 0.1327
)

# Feedstock quality bins (assumptions.yaml `quality_bins` section), inlined
# here as defaults so this module doesn't need to read the yaml at runtime.
# Source: Milledge et al. 2020, Energies, Table 1 -- Turks and Caicos
# pelagic species averages for `pelagic_high_quality`; mixed Sargassum
# inundation for `beach_wrack_low_quality`.
QUALITY_BINS = {
    "pelagic_high_quality": dict(
        moisture_frac=0.8673,
        ash_wt_frac_dry=0.3453,
    ),
    "beach_wrack_low_quality": dict(
        moisture_frac=0.8198,
        ash_wt_frac_dry=0.4694,
    ),
}

def make_sargassum_feed(fresh_feed_kgph: float, moisture_frac: float, quality: str):
    if quality not in QUALITY_BINS:
        raise KeyError(f"Unknown quality '{quality}'. Options: {list(QUALITY_BINS.keys())}")
    ash_target = QUALITY_BINS[quality]["ash_wt_frac_dry"]  # dry-basis ash

    water_kgph = fresh_feed_kgph * moisture_frac
    dry_kgph   = fresh_feed_kgph * (1 - moisture_frac)

    # Convert wet composition -> base dry-basis fractions (excluding water)
    solids = {k: v for k, v in WET_COMPOSITION.items() if k != "Water"}
    solids_sum = sum(solids.values())
    base_dry = {k: v / solids_sum for k, v in solids.items()}  # sums to 1

    base_ash = base_dry["Ash"]
    base_nonash_sum = 1.0 - base_ash
    target_nonash_sum = 1.0 - ash_target

    if base_nonash_sum <= 1e-12:
        raise RuntimeError("Base non-ash fraction is ~0; cannot scale.")
    scale = target_nonash_sum / base_nonash_sum

    # Build final dry-basis fractions with ash overridden
    dry_fracs = {"Ash": ash_target}
    for k, v in base_dry.items():
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

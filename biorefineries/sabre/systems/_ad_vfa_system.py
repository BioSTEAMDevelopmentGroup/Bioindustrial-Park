# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
VFA-producing acidogenic AD system builder for the SaBRe flowsheets.
"""
import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
    get_ad_temperature_K, non_none,
)
from biorefineries.sabre.units import (
    AcidogenicDigester, DigestateScrewPress, PressateConcentrator,
    BiostimulantEvaporator, Press, Mill,
)

__all__ = ('create_ad_vfa_system',)


# data/ad.yaml, loaded once -- AD sizing/operating parameters, shared
# with systems._ad_biomethane_system and/or systems._integrated_system.
# (AcidogenicDigester/DigestateScrewPress source their own defaults
# directly from this same file, so their sub-blocks aren't re-extracted
# here.)
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]


def create_ad_vfa_system(
    feedstock_type: str = "pelagic",
    milled_biomass_stream=None,
    enable_heat_shock: bool = None,
    hs_target_temperature_K: float = None,
    hs_events_per_day: float = None,
    hs_heated_fraction_of_liquid: float = None,
    hs_duration_min: float = None,
    temperature_regime: str = "mesophilic",
):
    """
    Build the VFA acidogenic AD subsystem.

    Pathway:
        Standalone:
            raw Sargassum
              -> Press
              -> optional PressateConcentrator / optional BiostimulantEvaporator
              -> Mill
              -> VFA_AD
              -> ScrewPress

        Integrated:
            milled_biomass_stream
              -> VFA_AD
              -> ScrewPress

    Key outputs (accessible via flowsheet):
        - offgas
        - vfa_broth
        - acidogenic_residual_solids
        - biostimulant_membrane_concentrate
        - pressate_permeate
        - biostimulant_product
    """
    vfaS = {**_AD_SHARED, **_AD_SHARED.get("acidogenic", {})}
    vfaS["temperature_K"] = get_ad_temperature_K(_AD_SHARED, temperature_regime)

    path = []

    if milled_biomass_stream is not None:
        # -------------------------
        # INTEGRATED MODE
        # Preprocessing is handled upstream (shared train).
        # -------------------------
        ad_inlet = milled_biomass_stream

    else:
        # -------------------------
        # STANDALONE MODE
        # Build full preprocessing from raw feed.
        # -------------------------
        feedstock_assumptions = load_assumptions("feedstock.yaml")
        params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
        fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
        moisture_frac = params["moisture_frac"]

        feed = make_sargassum_feed(
            fresh_feed_kgph=fresh_feed_kgph,
            moisture_frac=moisture_frac,
            ash_wt_frac_dry=params["ash_wt_frac_dry"],
        )

        # Press/PC/EV/Mill all have their own yaml-sourced defaults
        # (data/preprocessing.yaml, data/biostimulant.yaml) -- not
        # re-declared here.
        PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))
        path.append(PR)

        PC = PressateConcentrator(
            "PC_VFA",
            ins=PR - 1,
            outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
        )
        path.append(PC)

        biostimulant_fresh_water_vfa = bst.Stream(
            "biostimulant_fresh_water_vfa", Water=0.0, units="kg/hr"
        )
        EV = BiostimulantEvaporator(
            "EV_VFA",
            ins=(PC - 0, PC - 1, biostimulant_fresh_water_vfa),
            outs=("biostimulant_product", "biostimulant_vapor", "residual_permeate"),
        )
        path.append(EV)

        ML = Mill("ML", ins=PR - 0, outs=("milled_biomass", "milling_losses"))
        path.append(ML)
        ad_inlet = ML - 0

    # -------------------------
    # Acidogenic digester
    # -------------------------
    AD = AcidogenicDigester(
        "VFA_AD",
        ins=ad_inlet,
        outs=("offgas", "acidogenic_broth"),
        # vs_destruction/vfa_kg_per_kg_vs/vfa_split/digestible_IDs/hrt_days/
        # headspace_frac/max_single_digester_volume_MG/mixing_W_per_m3 all
        # match AcidogenicDigester's own data/ad.yaml-sourced defaults
        # exactly (verified) -- not re-declared here.
        target_temperature_K=float(vfaS["temperature_K"]),
        **non_none(
            enable_heat_shock=enable_heat_shock,
            hs_target_temperature_K=hs_target_temperature_K,
            hs_events_per_day=hs_events_per_day,
            hs_heated_fraction_of_liquid=hs_heated_fraction_of_liquid,
            hs_duration_min=hs_duration_min,
        ),
    )

    # DigestateScrewPress sources these from data/ad.yaml as its own
    # class-level defaults -- not re-declared here.
    SP = DigestateScrewPress(ID="SP_VFA", ins=AD - 1, outs=("acidogenic_residual_solids", "vfa_broth"))

    path.extend([AD, SP])

    sys_id = "AD_VFA_sys" if milled_biomass_stream is None else "AD_VFA_subsys"
    return bst.System(sys_id, path=tuple(path))

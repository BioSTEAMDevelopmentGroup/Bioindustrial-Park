# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
AD/biomethane system builder for the SaBRe (Sargassum Biorefinery) flowsheets.

Purpose:
- Build a flowsheet block for the methanogenic AD system
- Return a BioSTEAM System for simulation and diagramming

Key entry points:
- create_ad_biomethane_system(...)
- _build_methanogenic_pathway(...) — shared with systems._integrated_system,
  which routes only a fraction of milled biomass through this same pathway.

Notes:
- feed -> press -> pressate concentrator -> mill
  -> optional pretreatment -> AD -> H2S removal -> biogas upgrading -> digestate separation
- Uses plant-scale throughput from YAML (e.g. 15,000 ton/day wet feed)
- Final product is biomethane (post-upgrading), not the raw AD biogas.
"""

import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
    get_ad_temperature_K, non_none,
)
from biorefineries.sabre.units import (
    AnaerobicDigester, BiogasUpgrading, H2SRemoval, DigestateScrewPress,
    Press, Mill, HeatingPretreatment, EnzymaticPretreatment,
    PeroxidePretreatment, PressateConcentrator, BiostimulantEvaporator,
)

__all__ = ('create_ad_biomethane_system',)


# Pretreatment case definitions, shared by _build_methanogenic_pathway()
# below (data/pretreatment.yaml `pretreatment_ad` section).
_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]

# data/ad.yaml, loaded once -- AD sizing/operating parameters and
# performance parameters, shared with systems._ad_vfa_system and/or
# systems._integrated_system. (H2SRemoval/BiogasUpgrading/
# DigestateScrewPress source their own defaults directly from this same
# file, so their sub-blocks aren't re-extracted here.)
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]


def _apply_biodegradability_overrides(base_dict, override_factors):
    """
    Multiply selected biodegradability values by override factors, clamped to [0, 1].
    """
    out = dict(base_dict)
    for cid, factor in (override_factors or {}).items():
        if cid not in out:
            continue
        out[cid] = min(max(out[cid] * float(factor), 0.0), 1.0)
    return out


def _build_methanogenic_pathway(
    feed_stream, pretreatment_case, temperature_regime="mesophilic",
    biogas_ids=("biogas", "digestate"),
):
    """
    Build [optional pretreatment] -> AD -> H2S removal -> biogas upgrading
    -> digestate screw press, starting from an already-milled feed stream.

    Shared by create_ad_biomethane_system() (which supplies a Press+Mill-derived
    feed) and systems._integrated_system.create_integrated_biorefinery()
    (which supplies a splitter-derived feed). Returns (path_units, streams,
    units) so callers can assemble their own bst.System.
    """
    adS = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
    adS["temperature_K"] = get_ad_temperature_K(_AD_SHARED, temperature_regime)
    adp = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
    pretreatments = _PRETREATMENT_AD

    enable_feed_dilution = bool(adS.get("enable_feed_dilution", True))
    target_feed_moisture_frac = float(adS.get("target_feed_moisture_frac", 0.92))

    selected_pretreatment = pretreatment_case or adS.get("pretreatment_case", "press_mill_only")
    pt_case = pretreatments.get(selected_pretreatment, {})
    pt_kind = pt_case.get("kind", "none")

    ad_feed = feed_stream
    pt_units = []
    HT = EZ = PX = None

    if pt_kind == "none":
        pass

    # EnzymaticPretreatment/PeroxidePretreatment/HeatingPretreatment all
    # source these same parameters from the same pretreatment.yaml paths
    # as their own class-level defaults (enzymatic/enzymatic,
    # peroxide/peroxide, combined_PTE/heating respectively) -- not
    # re-declared here.
    elif pt_kind == "enzymatic":
        EZ = EnzymaticPretreatment("EZ", ins=ad_feed, outs=("enzyme_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.append(EZ)

    elif pt_kind == "peroxide":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        ad_feed = PX - 0
        pt_units.append(PX)

    elif pt_kind == "combined_PE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=PX - 0, outs=("combined_PE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, EZ])

    elif pt_kind == "combined_PTE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        HT = HeatingPretreatment("HT", ins=PX - 0, outs=("heated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=HT - 0, outs=("combined_PTE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, HT, EZ])

    else:
        raise ValueError(f"Unknown pretreatment_case '{selected_pretreatment}'")

    ad_effects = pt_case["ad_effects"]
    vs_destruction = float(ad_effects["vs_destruction"])
    ch4_kg_per_kg_vs_fed = float(ad_effects["ch4_kg_per_kg_vs_fed"])
    raw_biogas_molfrac = dict(ad_effects["raw_biogas_molfrac"])
    biodegradability = _apply_biodegradability_overrides(
        dict(adp["biodegradability"]), ad_effects.get("biodegradability_factor_overrides", {})
    )

    MX = None
    if enable_feed_dilution:
        if not (0.0 < target_feed_moisture_frac < 1.0):
            raise ValueError(
                f"target_feed_moisture_frac must be between 0 and 1, got {target_feed_moisture_frac}"
            )
        dilution_water = bst.Stream("dilution_water", Water=0.0, units="kg/hr")
        MX = bst.Mixer("MX", ins=(ad_feed, dilution_water), outs=("diluted_ad_feed",))

        @MX.add_specification(run=False)
        def _set_dilution_water():
            feed = MX.ins[0]
            total_mass = float(feed.F_mass)
            water_mass = float(feed.imass["Water"]) if "Water" in feed.chemicals else 0.0
            dry_mass = max(total_mass - water_mass, 0.0)
            if dry_mass <= 0.0:
                water_to_add = 0.0
            else:
                target_water_mass = dry_mass * target_feed_moisture_frac / (1.0 - target_feed_moisture_frac)
                water_to_add = max(target_water_mass - water_mass, 0.0)
            dilution_water.empty()
            dilution_water.imass["Water"] = water_to_add

        ad_feed = MX - 0
        pt_units.append(MX)

    # vs_destruction/ch4_kg_per_kg_vs_fed/raw_biogas_molfrac/biodegradability
    # are pretreatment-case-specific (data/pretreatment.yaml ad_effects) and
    # target_temperature_K depends on temperature_regime -- both must stay
    # explicit here. The rest (hrt_days, headspace_frac,
    # max_single_digester_volume_MG, maintenance_usd_per_m3_yr,
    # mixing_W_per_m3, digestible_IDs) match AnaerobicDigester's own
    # data/ad.yaml-sourced defaults exactly, so aren't re-declared.
    AD = AnaerobicDigester(
        "AD", ins=ad_feed, outs=biogas_ids,
        vs_destruction=vs_destruction,
        ch4_kg_per_kg_vs_fed=ch4_kg_per_kg_vs_fed,
        raw_biogas_molfrac=raw_biogas_molfrac,
        biodegradability=biodegradability,
        target_temperature_K=adS["temperature_K"],
    )

    # H2SRemoval/BiogasUpgrading/DigestateScrewPress all source these same
    # parameters from data/ad.yaml as their own class-level defaults -- not
    # re-declared here.
    H2SR = H2SRemoval("H2SR", ins=AD - 0, outs=("treated_biogas", "spent_h2s_media"))

    UP = BiogasUpgrading("UP", ins=H2SR - 0, outs=("biomethane", "offgas"))

    SP = DigestateScrewPress(ID="SP", ins=AD - 1, outs=("soil_amendment", "liquid_digestate"))

    path_units = [*pt_units, AD, H2SR, UP, SP]
    streams = {
        "biomethane": UP.outs[0], "offgas": UP.outs[1],
        "soil_amendment": SP.outs[0], "liquid_digestate": SP.outs[1],
    }
    units = {"AD": AD, "H2SR": H2SR, "UP": UP, "SP": SP}
    return path_units, streams, units


def create_ad_biomethane_system(
    feedstock_type: str = "pelagic",
    pretreatment_case: str | None = None,
    press_cake_solids_wt_frac: float | None = None,
    ch4_override=None,
    temperature_regime: str = "mesophilic",
):
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)

    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    # Press/PC/EV/Mill all have their own yaml-sourced defaults
    # (data/preprocessing.yaml, data/biostimulant.yaml) -- not
    # re-declared here.
    PR = Press(
        "PR", ins=feed, outs=("pressed_cake", "pressate"),
        **non_none(cake_solids_wt_frac=press_cake_solids_wt_frac),
    )

    PC = PressateConcentrator(
        "PC", ins=PR - 1,
        outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
    )

    biostimulant_fresh_water = bst.Stream("biostimulant_fresh_water", Water=0.0, units="kg/hr")
    EV = BiostimulantEvaporator(
        "EV", ins=(PC - 0, PC - 1, biostimulant_fresh_water),
        outs=("biostimulant_product", "biostimulant_vapor", "residual_permeate"),
    )

    ML = Mill("ML", ins=PR - 0, outs=("milled_biomass", "milling_losses"))

    resolved_pretreatment_case = pretreatment_case or _AD_SHARED.get("methanogenic", {}).get("pretreatment_case", "press_mill_only")
    path_units, streams, units = _build_methanogenic_pathway(
        ML - 0, resolved_pretreatment_case, temperature_regime=temperature_regime,
    )

    if ch4_override is not None:
        units["AD"].ch4_kg_per_kg_vs_fed = float(ch4_override)

    path = [PR, PC, EV, ML]
    path.extend(path_units)

    sys = bst.System("AD_Biomethane_sys", path=path)
    return sys

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
AD/biomethane system builder for the SaBRe (Sargassum Biorefinery) flowsheets.
"""

import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
    get_ad_temperature_K,
)
from biorefineries.sabre.units import (
    AnaerobicDigester, BiogasUpgrading, H2SRemoval, DigestateScrewPress,
    Press, Mill, HeatingPretreatment, EnzymaticPretreatment,
    PeroxidePretreatment, PressateConcentrator, BiostimulantEvaporator,
)

__all__ = ('create_ad_biomethane_system',)


# Load assumptions
_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]


def _build_methanogenic_pathway(
    feed_stream, pretreatment_case, temperature_regime="mesophilic",
):
    """
    Build [optional pretreatment] -> AD -> H2S removal -> biogas upgrading
    -> digestate screw press, starting from an already-milled feed stream.
    """
    adS = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
    adS["temperature_K"] = get_ad_temperature_K(_AD_SHARED, temperature_regime)
    adp = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
    pretreatments = _PRETREATMENT_AD

    ad_feed = feed_stream
    pt_units = []
    HT = EZ = PX = None

    if pretreatment_case == "press_mill_only":
        pass

    elif pretreatment_case == "enzymatic":
        EZ = EnzymaticPretreatment("EZ", ins=ad_feed, outs=("enzyme_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.append(EZ)

    elif pretreatment_case == "peroxide":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        ad_feed = PX - 0
        pt_units.append(PX)

    elif pretreatment_case == "combined_PE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=PX - 0, outs=("combined_PE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, EZ])

    elif pretreatment_case == "combined_PTE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        HT = HeatingPretreatment("HT", ins=PX - 0, outs=("heated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=HT - 0, outs=("combined_PTE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, HT, EZ])

    else:
        raise ValueError(f"Unknown pretreatment_case '{pretreatment_case}'")

    pt_case = pretreatments[pretreatment_case]
    ad_effects = pt_case["ad_effects"]

    AD = AnaerobicDigester(
        "AD", ins=ad_feed, outs=("biogas", "digestate"),
        vs_destruction=float(ad_effects["vs_destruction"]),
        ch4_kg_per_kg_vs_fed=float(ad_effects["ch4_kg_per_kg_vs_fed"]),
        raw_biogas_molfrac=dict(ad_effects["raw_biogas_molfrac"]),
        biodegradability=dict(adp["biodegradability"]),
        target_temperature_K=adS["temperature_K"],
    )


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
    pretreatment_case: str = 'press_mill_only',
):
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    feed = make_sargassum_feed(
        fresh_feed_kgph=get_scale_feed_kgph(feedstock_assumptions),
        moisture_frac=params["moisture_frac"],
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))

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

    path_units, streams, units = _build_methanogenic_pathway(ML - 0, pretreatment_case)

    path = [PR, PC, EV, ML]
    path.extend(path_units)

    sys = bst.System("AD_Biomethane_sys", path=path)
    return sys

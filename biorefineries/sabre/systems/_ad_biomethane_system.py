# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BIOSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
AD/biomethane system builder for the SaBRe (Sargassum Biorefinery) flowsheets.
"""

import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import load_assumptions, get_ad_temperature_K
from biorefineries.sabre.units import (
    AnaerobicDigester, BiogasUpgrading, H2SRemoval, DigestateScrewPress,
    Mill, HeatingPretreatment, EnzymaticPretreatment, PeroxidePretreatment,
)
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system

__all__ = ('create_ad_biomethane_system',)


# Load assumptions
_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]


def create_ad_biomethane_system(
    feedstock: str | bst.Stream = "pelagic",
    pretreatment_case: str = 'press_mill_only',
):
    """
    Build the AD/biomethane system: [optional pretreatment] -> AD -> H2S
    removal -> biogas upgrading -> digestate screw press.

    Parameters
    ----------
    feedstock : str or stream
        If a str, it's a feedstock type (data/feedstock.yaml
        `feedstock_type`) and create_biostimulant_system() is called to
        build Press -> PressateConcentrator -> BiostimulantEvaporator,
        followed by a Mill on the pressed cake, to get a milled feed.
        If a stream, it's used directly as the already-milled feed (e.g.
        a splitter-derived stream from systems._ad_integrated_system,
        which builds its own shared preprocessing once).
    pretreatment_case : str
        data/pretreatment.yaml `pretreatment_ad` case name.
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    path = []

    if isinstance(feedstock, str):
        bio_sys, bio_streams, bio_units = create_biostimulant_system(feedstock_type=feedstock)
        path.extend(bio_sys.units)

        ML = Mill("ML", ins=bio_streams["pressed_cake"], outs=("milled_biomass", "milling_losses"))
        path.append(ML)
        ad_feed = ML - 0
    else:
        ad_feed = feedstock

    adS = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
    adS["temperature_K"] = get_ad_temperature_K(_AD_SHARED, "mesophilic")
    adp = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
    pretreatments = _PRETREATMENT_AD

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

    path.extend([*pt_units, AD, H2SR, UP, SP])

    sys = bst.System("ad_biomethane_sys", path=path)
    return sys

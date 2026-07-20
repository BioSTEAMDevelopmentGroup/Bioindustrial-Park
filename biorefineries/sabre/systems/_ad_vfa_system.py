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
)
from biorefineries.sabre.units import (
    AcidogenicDigester, DigestateScrewPress, PressateConcentrator,
    BiostimulantEvaporator, Press, Mill,
)

__all__ = ('create_ad_vfa_system',)


def _build_acidogenic_pathway(feed_stream):
    """
    Build VFA_AD -> SP_VFA, starting from an already-milled feed stream.
    It does not include downstream usage of the VFA sream.
    """
    AD = AcidogenicDigester("VFA_AD", ins=feed_stream, outs=("offgas", "acidogenic_broth"))
    SP = DigestateScrewPress(ID="SP_VFA", ins=AD - 1, outs=("acidogenic_residual_solids", "vfa_broth"))

    path_units = [AD, SP]
    streams = {
        "offgas": AD.outs[0],
        "acidogenic_residual_solids": SP.outs[0],
        "vfa_broth": SP.outs[1],
    }
    units = {"VFA_AD": AD, "SP_VFA": SP}
    return path_units, streams, units


def create_ad_vfa_system(feedstock_type: str = "pelagic"):
    """
    Build the standalone VFA acidogenic AD subsystem, from raw feedstock:

        raw Sargassum
          -> Press
          -> optional PressateConcentrator / optional BiostimulantEvaporator
          -> Mill
          -> VFA_AD
          -> ScrewPress

    For an already-milled feed stream (e.g. shared preprocessing in
    systems._ad_integrated_system), use _build_acidogenic_pathway() directly
    instead.

    Key outputs (accessible via flowsheet):
        - offgas
        - vfa_broth
        - acidogenic_residual_solids
        - biostimulant_membrane_concentrate
        - pressate_permeate
        - biostimulant_product
    """
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph,
        moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))

    PC = PressateConcentrator(
        "PC_VFA",
        ins=PR - 1,
        outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
    )

    biostimulant_fresh_water_vfa = bst.Stream(
        "biostimulant_fresh_water_vfa", Water=0.0, units="kg/hr"
    )
    EV = BiostimulantEvaporator(
        "EV_VFA",
        ins=(PC - 0, PC - 1, biostimulant_fresh_water_vfa),
        outs=("biostimulant_product", "biostimulant_vapor", "residual_permeate"),
    )

    ML = Mill("ML", ins=PR - 0, outs=("milled_biomass", "milling_losses"))

    path_units, streams, units = _build_acidogenic_pathway(ML - 0)
    path = [PR, PC, EV, ML, *path_units]

    return bst.System("AD_VFA_sys", path=tuple(path))

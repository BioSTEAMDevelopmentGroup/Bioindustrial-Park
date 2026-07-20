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

from biorefineries.sabre.units import AcidogenicDigester, DigestateScrewPress, Mill
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system

__all__ = ('create_ad_vfa_system',)


def create_ad_vfa_system(feedstock: str | bst.Stream = "pelagic"):
    """
    Build the VFA acidogenic AD system: VFA_AD -> SP_VFA.

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

    Key outputs (accessible via flowsheet):
        - offgas
        - vfa_broth
        - acidogenic_residual_solids
        - biostimulant_membrane_concentrate (standalone mode only)
        - pressate_permeate (standalone mode only)
        - biostimulant_product (standalone mode only)
    """
    path = []

    if isinstance(feedstock, str):
        bio_sys, bio_streams, bio_units = create_biostimulant_system(feedstock_type=feedstock)
        path.extend(bio_sys.units)

        ML = Mill("ML", ins=bio_streams["pressed_cake"], outs=("milled_biomass", "milling_losses"))
        path.append(ML)
        ad_inlet = ML - 0
    else:
        ad_inlet = feedstock

    AD = AcidogenicDigester("VFA_AD", ins=ad_inlet, outs=("offgas", "acidogenic_broth"))
    SP = DigestateScrewPress(ID="SP_VFA", ins=AD - 1, outs=("acidogenic_residual_solids", "vfa_broth"))
    path.extend([AD, SP])

    return bst.System("AD_VFA_sys", path=tuple(path))

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

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.units import (
    AcidogenicDigester,
    DigestateScrewPress,
    Mill,
    VFAMicrofilter,
)
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system

__all__ = ('create_ad_vfa_system',)


def create_ad_vfa_system(feedstock: str | bst.Stream = "pelagic"):
    """
    Build the VFA acidogenic AD system: VFA_AD -> SP_VFA -> MF.

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
        - vfa_broth (post-microfiltration permeate; ready for fermentation)
        - vfa_retentate
        - acidogenic_residual_solids
        - biostimulant_membrane_concentrate (standalone mode only)
        - pressate_permeate (standalone mode only)
        - biostimulant_product (standalone mode only)
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    path = []

    if isinstance(feedstock, str):
        bio_sys, bio_streams, bio_units = create_biostimulant_system(feedstock_type=feedstock)
        # Fold in the biostimulant subsystem's units, but not its own HXN facility --
        # this system gets its own HXN below, scoped to all units visible here, so
        # nesting the subsystem's narrower one would double-count already-optimized
        # utilities.
        path.extend(u for u in bio_sys.units if not isinstance(u, bst.HeatExchangerNetwork))

        ML = Mill("ML", ins=bio_streams["pressed_cake"], outs=("milled_biomass", "milling_losses"))
        path.append(ML)
        ad_inlet = ML - 0
    else:
        ad_inlet = feedstock

    AD = AcidogenicDigester("VFA_AD", ins=ad_inlet, outs=("offgas", "acidogenic_broth"))
    SP = DigestateScrewPress(ID="SP_VFA", ins=AD - 1, outs=("acidogenic_residual_solids", "raw_vfa_broth"))
    MF = VFAMicrofilter("MF", ins=SP - 1, outs=("vfa_broth", "vfa_retentate"))
    path.extend([AD, SP, MF])

    HXN = bst.HeatExchangerNetwork("HXN", units=tuple(path))
    path.append(HXN)

    return bst.System("ad_vfa_sys", path=tuple(path))

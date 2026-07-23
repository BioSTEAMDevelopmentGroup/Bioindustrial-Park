# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Integrated Sargassum Biorefinery System Builder
================================================
Shared preprocessing (create_biostimulant_system's Press -> PC -> PSP/DIL ->
EV, plus a Mill on the pressed cake) then splits milled biomass:
  alpha     -> Methanogenic AD pathway -> biomethane
  (1-alpha) -> VFA fermentation pathway -> microbial oil

Edge cases handled cleanly:
  alpha=0.0 -> pure VFA-to-oil  (methanogenic pathway not built)
  alpha=1.0 -> pure biomethane  (VFA/fermentation pathway not built)

Each pathway is built by calling create_ad_biomethane_system() /
create_ad_fermentation_system() directly, passing the splitter-derived
stream as `feedstock` (so each skips its own create_biostimulant_system
call and uses the shared preprocessing built here instead) -- so this
integrated system and the standalone AD systems never drift apart.
"""

from pathlib import Path

import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.units import Mill
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system
from biorefineries.sabre.systems._ad_biomethane_system import create_ad_biomethane_system
from biorefineries.sabre.systems._ad_fermentation_system import create_ad_fermentation_system

__all__ = ('create_ad_integrated_system', 'MassSplitter')


class MassSplitter(bst.Unit):
    """
    Splits all components by the same mass fraction.
    outs[0] gets fraction alpha (-> methanogenic AD).
    outs[1] gets fraction (1-alpha) (-> VFA fermentation).
    """
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID="", ins=None, outs=(), alpha=0.5, **kwargs):
        super().__init__(ID, ins, outs, **kwargs)
        self.alpha = float(alpha)

    def _run(self):
        feed = self.ins[0]
        to_methane, to_vfa = self.outs
        to_methane.empty()
        to_vfa.empty()
        to_methane.phase = feed.phase
        to_vfa.phase = feed.phase
        alpha = min(max(self.alpha, 0.0), 1.0)
        for cid in feed.chemicals.IDs:
            m = float(feed.imass[cid])
            to_methane.imass[cid] = alpha * m
            to_vfa.imass[cid] = (1.0 - alpha) * m

    def _design(self):
        self.design_results["Alpha (to methane AD)"] = self.alpha
        self.design_results["Alpha (to VFA oil)"] = 1.0 - self.alpha

    def _cost(self):
        pass


def create_ad_integrated_system(
    alpha: float = 0.5,
    feedstock: str = "pelagic",
    pretreatment_case: str = "press_mill_only",
):
    """
    Build the full integrated Sargassum biorefinery.

    Parameters
    ----------
    alpha : float in [0, 1]
        Fraction of milled biomass routed to methanogenic AD.
        alpha=0 -> pure VFA-to-oil (methane path not built).
        alpha=1 -> pure biomethane (VFA path not built).
    feedstock : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    pretreatment_case : str
        data/pretreatment.yaml `pretreatment_ad` case name.

    Returns
    -------
    sys : bst.System
        Key streams and units are accessible via `sys.flowsheet.stream`
        and `sys.flowsheet.unit` (e.g. 'feed', 'biomethane',
        'microbial_oil', 'PR', 'SPL'). Streams belonging to a pathway
        that wasn't built (alpha=0 or alpha=1) are simply absent.
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    alpha = float(alpha)
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha must be in [0, 1], got {alpha}")

    build_methane = alpha > 0.0
    build_vfa     = alpha < 1.0

    # =========================================================
    # SHARED PREPROCESSING: create_biostimulant_system's Press -> PC ->
    # PSP/DIL -> EV, plus a Mill on the pressed cake.
    # =========================================================
    bio_sys = create_biostimulant_system(feedstock=feedstock)

    # milling_losses: no price -- pure mass loss, not a disposed waste stream.
    ML = Mill("ML", ins=bio_sys.flowsheet.stream.pressed_cake, outs=("milled_biomass", "milling_losses"))

    # =========================================================
    # SPLITTER
    # =========================================================
    to_methane_ad = bst.Stream("to_methane_ad")
    to_vfa_ad     = bst.Stream("to_vfa_ad")

    SPL = MassSplitter(
        "SPL",
        ins=ML - 0,
        outs=(to_methane_ad, to_vfa_ad),
        alpha=alpha,
    )

    # =========================================================
    # BUILD PATHWAYS CONDITIONALLY
    # =========================================================
    methane_units = []
    if build_methane:
        methane_sys = create_ad_biomethane_system(
            feedstock=SPL - 0, pretreatment_case=pretreatment_case,
        )
        methane_units = list(methane_sys.units)

    vfa_units = []
    if build_vfa:
        vfa_fer_sys = create_ad_fermentation_system(feedstock=SPL - 1)
        vfa_units = list(vfa_fer_sys.units)

    # =========================================================
    # ASSEMBLE FULL SYSTEM
    # =========================================================
    # Fold in every subsystem's units, but not their own HXN facilities -- this
    # system gets its own HXN below, scoped to all units visible here, so nesting
    # any subsystem's narrower one would double-count already-optimized utilities.
    preprocessing = [u for u in bio_sys.units if not isinstance(u, bst.HeatExchangerNetwork)] + [ML]
    methane_units = [u for u in methane_units if not isinstance(u, bst.HeatExchangerNetwork)]
    vfa_units = [u for u in vfa_units if not isinstance(u, bst.HeatExchangerNetwork)]
    all_units = preprocessing + [SPL] + methane_units + vfa_units

    HXN = bst.HeatExchangerNetwork("HXN", units=tuple(all_units))
    all_units = all_units + [HXN]

    sys = bst.System.from_units("ad_integrated_sys", units=all_units)

    return sys


if __name__ == '__main__':
    sys = create_ad_integrated_system(feedstock="pelagic")
    sys.simulate()

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")
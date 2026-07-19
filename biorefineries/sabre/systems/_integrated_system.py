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
Shared preprocessing (Press -> PC -> Mill) then splits milled biomass:
  alpha     -> Methanogenic AD pathway -> biomethane
  (1-alpha) -> VFA fermentation pathway -> microbial oil

Edge cases handled cleanly:
  alpha=0.0 -> pure VFA-to-oil  (methanogenic pathway not built)
  alpha=1.0 -> pure biomethane  (VFA/fermentation pathway not built)

The methanogenic pathway (pretreatment -> AD -> H2S removal -> biogas
upgrading -> digestate screw press) is built via the shared
_build_methanogenic_pathway() helper in systems._ad_biomethane_system, so this
integrated system and the standalone AD/biomethane system never drift apart.

Returns (sys, streams, units, alpha).
"""

import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
    non_none,
)
from biorefineries.sabre.units import Press, Mill, PressateConcentrator, BiostimulantEvaporator
from biorefineries.sabre.systems._ad_biomethane_system import _build_methanogenic_pathway
from biorefineries.sabre.systems._ad_vfa_system import create_ad_vfa_system
from biorefineries.sabre.systems._ad_fermentation_system import (
    _create_vfa_fermentation_system, _VFA_DOWNSTREAM,
)

__all__ = ('create_integrated_biorefinery', 'MassSplitter')


def _get_stream(stream_id: str):
    """Safe stream lookup — returns None if not in registry."""
    try:
        return bst.main_flowsheet.stream[stream_id]
    except Exception:
        return None


def _get_unit(unit_id: str):
    """Safe unit lookup — returns None if not in registry."""
    try:
        return bst.main_flowsheet.unit[unit_id]
    except Exception:
        return None


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


def _build_methane_pathway(ad_feed_in, pretreatment_case, temperature_regime="mesophilic"):
    """
    Thin wrapper: builds AD -> H2SR -> UP -> SP via the shared
    methanogenic-pathway builder in systems._ad_biomethane_system, so this
    integrated system and the standalone AD/biomethane system never drift
    apart. Returns (units_list, streams_dict, units_dict).
    """
    path_units, streams, units = _build_methanogenic_pathway(
        ad_feed_in, pretreatment_case, temperature_regime=temperature_regime,
    )
    return path_units, streams, units


def _build_vfa_pathway(vfa_stream, ferm_kwargs, temperature_regime="mesophilic"):
    """Build VFA_AD -> SP_VFA -> fermentation chain. Returns (units_list, streams_dict, units_dict)."""
    vfa_subsys = create_ad_vfa_system(
        milled_biomass_stream=vfa_stream, temperature_regime=temperature_regime,
    )

    vfa_broth = _get_stream("vfa_broth")
    if vfa_broth is None:
        raise RuntimeError(
            "Could not find 'vfa_broth' stream. Check SP_VFA output IDs in _ad_vfa_system.py."
        )

    fer_sys, fer_streams, fer_units = _create_vfa_fermentation_system(
        vfa_broth=vfa_broth, **ferm_kwargs
    )

    path_units = list(vfa_subsys.units) + list(fer_sys.units)

    streams = {
        "vfa_broth":               vfa_broth,
        "backend_oil":             fer_streams.get("backend_oil"),
        "fermentation_wastewater": fer_streams.get("fermentation_wastewater"),
        "vfa_retentate":           fer_streams.get("vfa_retentate"),
        "acidogenic_residual_solids": _get_stream("acidogenic_residual_solids"),
    }
    units = {
        "VFA_AD": _get_unit("VFA_AD"),
        "SP_VFA": _get_unit("SP_VFA"),
        **fer_units,
    }

    return path_units, streams, units


def create_integrated_biorefinery(
    alpha: float = 0.5,
    feedstock_type: str = "pelagic",
    pretreatment_case: str = "press_mill_only",
    vfa_conversion: float = None,
    vfa_product_yield: float = None,
    vfa_biomass_yield: float = None,
    vfa_co2_yield: float = None,
    vfa_o2_demand: float = None,
    ferm_residence_time_h: float = None,
    ferm_target_pH: float = None,
    ferm_mgso4_dose: float = None,
    target_oil_and_solids_content: float = _VFA_DOWNSTREAM["target_oil_and_solids_content_g_per_L"],
    temperature_regime: str = "mesophilic",
):
    """
    Build the full integrated Sargassum biorefinery.

    Parameters
    ----------
    alpha : float in [0, 1]
        Fraction of milled biomass routed to methanogenic AD.
        alpha=0 -> pure VFA-to-oil (methane path not built).
        alpha=1 -> pure biomethane (VFA path not built).

    Returns
    -------
    sys : bst.System
    streams : dict of key streams (None for streams not built at edge cases)
    units : dict of key units
    alpha : float
    """
    alpha = float(alpha)
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha must be in [0, 1], got {alpha}")

    build_methane = alpha > 0.0
    build_vfa     = alpha < 1.0

    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)

    # =========================================================
    # SHARED PREPROCESSING: Press -> PC -> Mill
    # Press/PC/EV/Mill all have their own yaml-sourced defaults
    # (data/preprocessing.yaml, data/biostimulant.yaml) -- not
    # re-declared here since none are exposed as overridable
    # parameters of create_integrated_biorefinery() itself.
    # =========================================================
    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph,
        moisture_frac=params["moisture_frac"],
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))

    PC = PressateConcentrator(
        "PC",
        ins=PR - 1,
        outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
    )

    biostimulant_fresh_water = bst.Stream(
        "biostimulant_fresh_water", Water=0.0, units="kg/hr"
    )
    EV = BiostimulantEvaporator(
        "EV",
        ins=(PC - 0, PC - 1, biostimulant_fresh_water),
        outs=("biostimulant_product", "biostimulant_vapor", "residual_permeate"),
    )

    ML = Mill("ML", ins=PR - 0, outs=("milled_biomass", "milling_losses"))

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
    ferm_kwargs = non_none(
        conversion=vfa_conversion,
        product_yield_kg_per_kg_vfa_consumed=vfa_product_yield,
        biomass_yield_kg_per_kg_vfa_consumed=vfa_biomass_yield,
        co2_yield_kg_per_kg_vfa_consumed=vfa_co2_yield,
        oxygen_kg_per_kg_vfa_consumed=vfa_o2_demand,
        residence_time_h=ferm_residence_time_h,
        target_pH=ferm_target_pH,
        magnesium_sulfate_dose_kg_per_m3=ferm_mgso4_dose,
        target_oil_and_solids_content=target_oil_and_solids_content,
    )

    methane_units   = []
    methane_streams = {}
    methane_units_d = {}
    if build_methane:
        methane_units, methane_streams, methane_units_d = _build_methane_pathway(
            SPL - 0, pretreatment_case, temperature_regime=temperature_regime,
        )

    vfa_units   = []
    vfa_streams = {}
    vfa_units_d = {}
    if build_vfa:
        vfa_units, vfa_streams, vfa_units_d = _build_vfa_pathway(
            SPL - 1, ferm_kwargs, temperature_regime=temperature_regime,
        )

    # =========================================================
    # ASSEMBLE FULL SYSTEM
    # =========================================================
    preprocessing = [PR, PC, EV, ML]
    all_units = preprocessing + [SPL] + methane_units + vfa_units

    sys = bst.System.from_units("integrated_biorefinery", units=all_units)

    # =========================================================
    # STREAMS DICT
    # =========================================================
    streams = {
        "feed":           feed,
        "milled_biomass": ML.outs[0],
        "to_methane_ad":  to_methane_ad,
        "to_vfa_ad":      to_vfa_ad,
        "biostimulant_membrane_concentrate": _get_stream("biostimulant_membrane_concentrate"),
        "biostimulant_product": _get_stream("biostimulant_product"),
        # Methane pathway (None if alpha=0)
        **{k: methane_streams.get(k) for k in [
            "biomethane", "offgas", "soil_amendment", "liquid_digestate"
        ]},
        # VFA pathway (None if alpha=1)
        **{k: vfa_streams.get(k) for k in [
            "vfa_broth", "backend_oil", "fermentation_wastewater",
            "vfa_retentate", "acidogenic_residual_solids",
        ]},
    }

    units = {
        "PR": PR, "PC": PC, "EV": EV, "ML": ML, "SPL": SPL,
        **methane_units_d,
        **vfa_units_d,
    }

    return sys, streams, units, alpha

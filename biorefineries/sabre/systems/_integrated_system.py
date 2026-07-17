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
)
from biorefineries.sabre.units import Press, Mill, PressateConcentrator
from biorefineries.sabre.systems._ad_biomethane_system import _build_methanogenic_pathway
from biorefineries.sabre.systems._ad_vfa_system import create_ad_vfa_system
from biorefineries.sabre.systems._ad_fermentation_system import _create_vfa_fermentation_system

__all__ = ('create_integrated_biorefinery', 'MassSplitter')


SOLIDS_IDS = [
    "Ash", "Protein", "Lignin", "Glucan", "Xylan", "Mannan", "Galactan", "Arabinan",
    "Alginate", "Fucoidan", "Mannitol", "OtherSolids",
]


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


def _build_methane_pathway(A, ad_feed_in, pretreatment_case):
    """
    Thin wrapper: builds AD -> H2SR -> UP -> SP via the shared
    methanogenic-pathway builder in systems._ad_biomethane_system, so this
    integrated system and the standalone AD/biomethane system never drift
    apart. Returns (units_list, streams_dict, units_dict).
    """
    path_units, streams, units = _build_methanogenic_pathway(A, ad_feed_in, pretreatment_case)
    return path_units, streams, units


def _build_vfa_pathway(vfa_stream, ferm_kwargs):
    """Build VFA_AD -> SP_VFA -> fermentation chain. Returns (units_list, streams_dict, units_dict)."""
    vfa_subsys = create_ad_vfa_system(milled_biomass_stream=vfa_stream)

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
    vfa_conversion: float = 0.85,
    vfa_product_yield: float = 0.144,
    vfa_biomass_yield: float = 0.40,
    vfa_co2_yield: float = 0.20,
    vfa_o2_demand: float = 0.80,
    ferm_residence_time_h: float = 48.0,
    ferm_target_pH: float = 8.0,
    ferm_mgso4_dose: float = 0.49,
    target_oil_and_solids_content: float = 70.0,
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

    A = load_assumptions()
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    preprocessing_assumptions = load_assumptions("preprocessing.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)

    # =========================================================
    # SHARED PREPROCESSING: Press -> PC -> Mill
    # =========================================================
    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph,
        moisture_frac=params["moisture_frac"],
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    pp  = preprocessing_assumptions.get("preprocessing", {})
    prA = pp.get("press", {})
    mlA = pp.get("mill", {})

    PR = Press(
        "PR",
        ins=feed,
        outs=("pressed_cake", "pressate"),
        solids_IDs=SOLIDS_IDS,
        solids_capture_frac=prA.get("solids_capture_frac", 0.98),
        cake_solids_wt_frac=prA.get("cake_solids_wt_frac", 0.35),
        solubles_to_pressate_frac=prA.get("solubles_to_pressate_frac", 1.0),
        power_kWh_per_dry_ton_TS=prA.get("power_kWh_per_dry_ton_TS"),
        ref_capacity_tph_wet=(prA.get("ref_capacity_tph_wet") or 50.0),
        capex_installed_ref_usd=(prA.get("capex_installed_ref_usd") or 5e6),
        scale_exponent=(prA.get("scale_exponent") or 0.6),
        F_BM=(prA.get("F_BM") or 1.0),
    )

    pb  = preprocessing_assumptions.get("pressate_biostimulant", {})
    pcA = pb.get("concentrator", {})
    PC  = None
    if pb.get("enabled", False) and pb.get("concentrate_pressate", False):
        PC = PressateConcentrator(
            "PC",
            ins=PR - 1,
            outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
            retained_solute_IDs=tuple(pcA.get(
                "retained_solute_IDs",
                ["Alginate", "Fucoidan", "Mannitol", "Protein", "OtherSolids"],
            )),
            water_recovery_to_permeate=pcA.get("water_recovery_to_permeate", 0.70),
            retained_solute_recovery_to_concentrate=pcA.get(
                "retained_solute_recovery_to_concentrate", 0.95
            ),
            design_flux_L_m2_h=pcA.get("design_flux_L_m2_h", 35.0),
            operating_pressure_bar=pcA.get("operating_pressure_bar", 5.0),
            electricity_kWh_per_m3_feed=pcA.get("electricity_kWh_per_m3_feed", 0.8),
            capex_usd_per_m2=pcA.get("capex_usd_per_m2", 120.0),
            maintenance_frac_of_capex_per_yr=pcA.get(
                "maintenance_frac_of_capex_per_yr", 0.035
            ),
        )

    ML = Mill(
        "ML",
        ins=PR - 0,
        outs=("milled_biomass", "milling_losses"),
        loss_frac=mlA.get("loss_frac", 0.15),
        power_kWh_per_dry_ton_dry=mlA.get("power_kWh_per_dry_ton_dry"),
        ref_capacity_dry_ton_per_hr=mlA.get("ref_capacity_dry_ton_per_hr", 10.0),
        purchase_cost_ref_usd=mlA.get("purchase_cost_ref_usd", 206400.0),
        install_factor=mlA.get("install_factor", 1.8),
        scale_exponent=mlA.get("scale_exponent", 0.6),
        F_BM=mlA.get("F_BM", 1.0),
    )

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
    ferm_kwargs = dict(
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
            A, SPL - 0, pretreatment_case
        )

    vfa_units   = []
    vfa_streams = {}
    vfa_units_d = {}
    if build_vfa:
        vfa_units, vfa_streams, vfa_units_d = _build_vfa_pathway(SPL - 1, ferm_kwargs)

    # =========================================================
    # ASSEMBLE FULL SYSTEM
    # =========================================================
    preprocessing = [PR] + ([PC] if PC else []) + [ML]
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
        "biostimulant_membrane_concentrate": (
            _get_stream("biostimulant_membrane_concentrate") if PC else None
        ),
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
        "PR": PR, "PC": PC, "ML": ML, "SPL": SPL,
        **methane_units_d,
        **vfa_units_d,
    }

    return sys, streams, units, alpha

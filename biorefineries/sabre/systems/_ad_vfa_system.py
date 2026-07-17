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


SOLIDS_IDS = [
    "Ash", "Protein", "Lignin", "Glucan", "Xylan", "Mannan", "Galactan", "Arabinan",
    "Alginate", "Fucoidan", "Mannitol", "OtherSolids",
]


def _get_ad_vfa_case(assumptions: dict) -> dict:
    """
    Merge top-level vfa_ad_performance with the selected case.
    Returns a flat dict of scalar parameters only.
    NOTE: vfa_split is a nested dict and must be read separately
    via _get_ad_vfa_split() to avoid being silently dropped on merge.
    """
    perf = assumptions["vfa_ad_performance"]
    case_name = perf.get("case")
    cases = perf.get("cases")
    if case_name and isinstance(cases, dict):
        if case_name not in cases:
            raise KeyError(f"Unknown VFA AD case: {case_name}")
        merged = dict(perf)
        merged.update(cases[case_name])
        return merged
    return perf


def _get_ad_vfa_split(assumptions: dict) -> dict | None:
    """
    Read vfa_split directly from the selected case dict.
    Bypasses the shallow-merge issue in _get_ad_vfa_case() that silently
    drops nested dicts. Returns a plain {str: float} dict or None.
    """
    perf = assumptions["vfa_ad_performance"]
    case_name = perf.get("case")
    cases = perf.get("cases", {})

    if not case_name or case_name not in cases:
        return None

    raw_split = cases[case_name].get("vfa_split")
    if not raw_split:
        return None

    return {str(k): float(v) for k, v in raw_split.items()}


def create_ad_vfa_system(
    feedstock_type: str = "pelagic",
    milled_biomass_stream=None,
    enable_heat_shock: bool = False,
    hs_target_temperature_K: float = 338.15,
    hs_events_per_day: float = 1.0 / 7.0,
    hs_heated_fraction_of_liquid: float = 0.10,
    hs_duration_min: float = 15.0,
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
    A = load_assumptions()

    vfaS = A["vfa_ad"]
    vfaP = _get_ad_vfa_case(A)
    vfa_split = _get_ad_vfa_split(A)

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
        preprocessing_assumptions = load_assumptions("preprocessing.yaml")
        params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
        fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
        moisture_frac = params["moisture_frac"]

        feed = make_sargassum_feed(
            fresh_feed_kgph=fresh_feed_kgph,
            moisture_frac=moisture_frac,
            ash_wt_frac_dry=params["ash_wt_frac_dry"],
        )

        pp = preprocessing_assumptions.get("preprocessing", {})
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

        path.append(PR)

        # -------------------------
        # Pressate biostimulant side-stream
        # -------------------------
        pb = preprocessing_assumptions.get("pressate_biostimulant", {})
        pcA = pb.get("concentrator", {})
        evA = pb.get("evaporator", {})

        PC = None
        EV = None

        if pb.get("enabled", False) and pb.get("concentrate_pressate", False):
            PC = PressateConcentrator(
                "PC_VFA",
                ins=PR - 1,
                outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
                retained_solute_IDs=tuple(
                    pcA.get(
                        "retained_solute_IDs",
                        ["Alginate", "Fucoidan", "Mannitol", "Protein", "OtherSolids"],
                    )
                ),
                water_recovery_to_permeate=pcA.get("water_recovery_to_permeate", 0.93),
                retained_solute_recovery_to_concentrate=pcA.get(
                    "retained_solute_recovery_to_concentrate", 0.95
                ),
                design_flux_L_m2_h=pcA.get("design_flux_L_m2_h", 35.0),
                operating_pressure_bar=pcA.get("operating_pressure_bar", 5.0),
                electricity_kWh_per_m3_feed=pcA.get("electricity_kWh_per_m3_feed", 2.5),
                capex_usd_per_m2=pcA.get("capex_usd_per_m2", 500.0),
                maintenance_frac_of_capex_per_yr=pcA.get(
                    "maintenance_frac_of_capex_per_yr", 0.035
                ),
            )
            path.append(PC)

            if evA.get("enabled", False):
                EV = BiostimulantEvaporator(
                    "EV_VFA",
                    ins=PC - 0,
                    outs=("biostimulant_product", "biostimulant_vapor"),
                    target_solids_wt_frac=evA.get("target_solids_wt_frac", 0.20),
                    boiling_temperature_K=evA.get("boiling_temperature_K", 333.15),
                    latent_heat_kJ_per_kg=evA.get("latent_heat_kJ_per_kg", 2350.0),
                    sensible_cp_kJ_per_kgK=evA.get("sensible_cp_kJ_per_kgK", 4.18),
                    electricity_kWh_per_kg_water_evap=evA.get(
                        "electricity_kWh_per_kg_water_evap", 0.0
                    ),
                    nonwater_recovery_to_product=evA.get(
                        "nonwater_recovery_to_product", 0.995
                    ),
                    capex_ref_usd=evA.get("capex_ref_usd", 750000.0),
                    ref_evaporation_kgph=evA.get("ref_evaporation_kgph", 10000.0),
                    scale_exponent=evA.get("scale_exponent", 0.60),
                    maintenance_frac_of_capex_per_yr=evA.get(
                        "maintenance_frac_of_capex_per_yr", 0.035
                    ),
                    F_BM=evA.get("F_BM", 1.0),
                )
                path.append(EV)

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

        path.append(ML)
        ad_inlet = ML - 0

    # -------------------------
    # Acidogenic digester
    # -------------------------
    AD = AcidogenicDigester(
        "VFA_AD",
        ins=ad_inlet,
        outs=("offgas", "acidogenic_broth"),
        vs_destruction=float(vfaP.get("vs_destruction", 0.50)),
        vfa_kg_per_kg_vs=float(vfaP.get("vfa_kg_per_kg_vs", 0.47)),
        vfa_split=vfa_split,
        digestible_IDs=vfaP.get(
            "digestible_IDs", A.get("ad_performance", {}).get("digestible_IDs")
        ),
        produce_offgas_co2=bool(vfaP.get("produce_offgas_co2", True)),
        hrt_days=float(vfaS.get("hrt_days", 15.0)),
        slurry_density_kg_per_m3=float(vfaS.get("slurry_density_kg_per_m3", 1000.0)),
        headspace_frac=float(vfaS.get("gas_storage_frac_of_total_volume", 0.2)),
        max_single_digester_volume_MG=float(vfaS.get("max_single_digester_volume_MG", 1.5)),
        mixing_W_per_m3=float(vfaS.get("mixing_W_per_m3", 5.0)),
        influent_temperature_K=float(vfaS.get("influent_temperature_K", 298.15)),
        target_temperature_K=float(vfaS.get("target_temperature_K", 308.15)),
        cp_kJ_per_kgK=float(vfaS.get("cp_kJ_per_kgK", 4.18)),
        enable_heat_shock=enable_heat_shock,
        hs_target_temperature_K=hs_target_temperature_K,
        hs_events_per_day=hs_events_per_day,
        hs_heated_fraction_of_liquid=hs_heated_fraction_of_liquid,
        hs_duration_min=hs_duration_min,
    )

    SP = DigestateScrewPress(
        ID="SP_VFA",
        ins=AD - 1,
        outs=("acidogenic_residual_solids", "vfa_broth"),
        solids_IDs=SOLIDS_IDS,
        ts_capture_frac=A.get("digestate_screw_press", {}).get("ts_capture_frac", 0.40),
        cake_moisture_frac=A.get("digestate_screw_press", {}).get("cake_moisture_frac", 0.50),
        capacity_tph_each=A.get("digestate_screw_press", {}).get("capacity_tph_each", 6.0),
        kWh_per_m3=A.get("digestate_screw_press", {}).get("kWh_per_m3", 0.67),
        eur_to_usd=A.get("digestate_screw_press", {}).get("eur_to_usd", 1.19),
        capex_eur_table=A.get("digestate_screw_press", {}).get("capex_eur_table"),
        include_polymer_dosing=A.get("digestate_screw_press", {}).get("include_polymer_dosing", False),
        polymer_dosing_cost_eur_each=A.get("digestate_screw_press", {}).get("polymer_dosing_cost_eur_each", 0.0),
        F_BM=A.get("digestate_screw_press", {}).get("F_BM", 1.0),
    )

    path.extend([AD, SP])

    sys_id = "AD_VFA_sys" if milled_biomass_stream is None else "AD_VFA_subsys"
    return bst.System(sys_id, path=tuple(path))

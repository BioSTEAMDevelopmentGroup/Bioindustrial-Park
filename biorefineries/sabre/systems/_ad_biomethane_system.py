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
    get_ad_temperature_K,
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

# data/ad.yaml, loaded once -- AD sizing/operating parameters, performance
# parameters, and the H2S removal / biogas upgrading / digestate screw
# press sections, all shared with systems._ad_vfa_system and/or
# systems._integrated_system.
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]
_H2S_REMOVAL = _AD_YAML["h2s_removal"]
_BIOGAS_UPGRADING = _AD_YAML["biogas_upgrading"]
_DIGESTATE_SCREW_PRESS = _AD_YAML["digestate_screw_press"]


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
    adC = _AD_SHARED["cost"]
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

    elif pt_kind == "heating":
        hA = pt_case.get("heating", {})
        HT = HeatingPretreatment(
            "HT", ins=ad_feed, outs=("heated_biomass",),
            target_temperature_K=hA.get("target_temperature_K", 338.15),
            residence_time_hr=hA.get("residence_time_hr", 0.25),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            cp_kJ_per_kgK=adS.get("cp_kJ_per_kgK", 4.18),
            capex_usd=hA.get("capex_usd", 0.0),
            maintenance_frac_of_capex_per_yr=hA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        ad_feed = HT - 0
        pt_units.append(HT)

    elif pt_kind == "enzymatic":
        eA = pt_case.get("enzymatic", {})
        EZ = EnzymaticPretreatment(
            "EZ", ins=ad_feed, outs=("enzyme_treated_biomass",),
            temperature_K=eA.get("temperature_K", 308.15),
            residence_time_hr=eA.get("residence_time_hr", 24.0),
            enzyme_dose_kg_per_kg_dry_feed=eA.get("enzyme_dose_kg_per_kg_dry_feed", 0.02),
            treated_fraction=eA.get("treated_fraction", 1.0),
            enzyme_recycle_factor=eA.get("enzyme_recycle_factor", 1.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=eA.get("capex_usd", 0.0),
            enzyme_price_usd_per_kg=eA.get("enzyme_price_usd_per_kg", 7.0),
            maintenance_frac_of_capex_per_yr=eA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        ad_feed = EZ - 0
        pt_units.append(EZ)

    elif pt_kind == "peroxide":
        pA = pt_case.get("peroxide", {})
        PX = PeroxidePretreatment(
            "PX", ins=ad_feed, outs=("peroxide_treated_biomass",),
            h2o2_wt_frac_on_dry_feed=pA.get("h2o2_wt_frac_on_dry_feed", 0.025),
            temperature_K=pA.get("temperature_K", 298.15),
            residence_time_hr=pA.get("residence_time_hr", 2.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=pA.get("capex_usd", 0.0),
            h2o2_price_usd_per_kg=pA.get("h2o2_price_usd_per_kg", 0.37),
            maintenance_frac_of_capex_per_yr=pA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        ad_feed = PX - 0
        pt_units.append(PX)

    elif pt_kind == "combined_PE":
        pA = pretreatments.get("peroxide", {}).get("peroxide", {})
        eA = pretreatments.get("enzymatic", {}).get("enzymatic", {})
        PX = PeroxidePretreatment(
            "PX", ins=ad_feed, outs=("peroxide_treated_biomass",),
            h2o2_wt_frac_on_dry_feed=pA.get("h2o2_wt_frac_on_dry_feed", 0.025),
            temperature_K=pA.get("temperature_K", 298.15),
            residence_time_hr=pA.get("residence_time_hr", 2.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=pA.get("capex_usd", 0.0),
            h2o2_price_usd_per_kg=pA.get("h2o2_price_usd_per_kg", 0.37),
            maintenance_frac_of_capex_per_yr=pA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        EZ = EnzymaticPretreatment(
            "EZ", ins=PX - 0, outs=("combined_PE_treated_biomass",),
            temperature_K=eA.get("temperature_K", 308.15),
            residence_time_hr=eA.get("residence_time_hr", 24.0),
            enzyme_dose_kg_per_kg_dry_feed=eA.get("enzyme_dose_kg_per_kg_dry_feed", 0.02),
            treated_fraction=eA.get("treated_fraction", 1.0),
            enzyme_recycle_factor=eA.get("enzyme_recycle_factor", 1.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=eA.get("capex_usd", 0.0),
            enzyme_price_usd_per_kg=eA.get("enzyme_price_usd_per_kg", 7.0),
            maintenance_frac_of_capex_per_yr=eA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        ad_feed = EZ - 0
        pt_units.extend([PX, EZ])

    elif pt_kind == "combined_PTE":
        pA = pretreatments.get("peroxide", {}).get("peroxide", {})
        hA = pt_case.get("heating", {})
        eA = pretreatments.get("enzymatic", {}).get("enzymatic", {})
        PX = PeroxidePretreatment(
            "PX", ins=ad_feed, outs=("peroxide_treated_biomass",),
            h2o2_wt_frac_on_dry_feed=pA.get("h2o2_wt_frac_on_dry_feed", 0.025),
            temperature_K=pA.get("temperature_K", 298.15),
            residence_time_hr=pA.get("residence_time_hr", 2.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=pA.get("capex_usd", 0.0),
            h2o2_price_usd_per_kg=pA.get("h2o2_price_usd_per_kg", 0.37),
            maintenance_frac_of_capex_per_yr=pA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        HT = HeatingPretreatment(
            "HT", ins=PX - 0, outs=("heated_biomass",),
            target_temperature_K=hA.get("target_temperature_K", 393.15),
            residence_time_hr=hA.get("residence_time_hr", 0.25),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            cp_kJ_per_kgK=adS.get("cp_kJ_per_kgK", 4.18),
            capex_usd=hA.get("capex_usd", 0.0),
            maintenance_frac_of_capex_per_yr=hA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        EZ = EnzymaticPretreatment(
            "EZ", ins=HT - 0, outs=("combined_PTE_treated_biomass",),
            temperature_K=eA.get("temperature_K", 308.15),
            residence_time_hr=eA.get("residence_time_hr", 24.0),
            enzyme_dose_kg_per_kg_dry_feed=eA.get("enzyme_dose_kg_per_kg_dry_feed", 0.02),
            treated_fraction=eA.get("treated_fraction", 1.0),
            enzyme_recycle_factor=eA.get("enzyme_recycle_factor", 1.0),
            slurry_density_kg_per_m3=adS.get("slurry_density_kg_per_m3", 1000.0),
            capex_usd=eA.get("capex_usd", 0.0),
            enzyme_price_usd_per_kg=eA.get("enzyme_price_usd_per_kg", 7.0),
            maintenance_frac_of_capex_per_yr=eA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
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

    AD = AnaerobicDigester(
        "AD", ins=ad_feed, outs=biogas_ids,
        vs_destruction=vs_destruction,
        ch4_kg_per_kg_vs_fed=ch4_kg_per_kg_vs_fed,
        raw_biogas_molfrac=raw_biogas_molfrac,
        digestible_IDs=tuple(adp["digestible_IDs"]),
        biodegradability=biodegradability,
        hrt_days=adS["hrt_days"],
        slurry_density_kg_per_m3=adS["slurry_density_kg_per_m3"],
        headspace_frac=adS["gas_storage_frac_of_total_volume"],
        max_single_digester_volume_MG=adS["max_single_digester_volume_MG"],
        maintenance_usd_per_m3_yr=adC.get("maintenance_usd_per_m3_yr", None),
        mixing_W_per_m3=adS["mixing_W_per_m3"],
        influent_temperature_K=adS["influent_temperature_K"],
        target_temperature_K=adS["temperature_K"],
        cp_kJ_per_kgK=adS["cp_kJ_per_kgK"],
    )

    h2sA = _H2S_REMOVAL
    H2SR = H2SRemoval(
        "H2SR", ins=AD - 0, outs=("treated_biogas", "spent_h2s_media"),
        h2s_removal_efficiency=h2sA.get("h2s_removal_efficiency", 0.99),
        reagent_cost_usd_per_Nm3_raw=h2sA.get("reagent_cost_usd_per_Nm3_raw", 0.002),
    )

    upA = _BIOGAS_UPGRADING
    UP = BiogasUpgrading(
        "UP", ins=H2SR - 0, outs=("biomethane", "offgas"),
        ch4_recovery=upA["ch4_recovery"],
        co2_removal=upA["co2_removal"],
        electricity_kwh_per_Nm3_raw=upA["electricity_kWh_per_Nm3_raw"],
        capex_usd_per_Nm3ph_raw=upA["capex_usd_per_Nm3ph_raw"],
        maintenance_frac_of_capex_per_yr=upA.get("maintenance_frac_of_capex_per_yr", 0.035),
    )

    sp = _DIGESTATE_SCREW_PRESS
    SP = DigestateScrewPress(
        ID="SP", ins=AD - 1, outs=("soil_amendment", "liquid_digestate"),
        ts_capture_frac=sp.get("ts_capture_frac", 0.40),
        cake_moisture_frac=sp.get("cake_moisture_frac", 0.50),
        capacity_tph_each=sp.get("capacity_tph_each", 6.0),
        kWh_per_m3=sp.get("kWh_per_m3", 0.67),
        include_polymer_dosing=sp.get("include_polymer_dosing", False),
        polymer_dosing_cost_usd_each=sp.get("polymer_dosing_cost_usd_each", 0.0),
        F_BM=sp.get("F_BM", 1.0),
    )

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
    preprocessing_assumptions = load_assumptions("preprocessing.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)

    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    pp = preprocessing_assumptions.get("preprocessing", {})
    prA = pp.get("press", {})
    mlA = pp.get("mill", {})

    PR = Press(
        "PR", ins=feed, outs=("pressed_cake", "pressate"),
        solids_IDs=tuple(prA.get("solids_IDs", [
            "Glucan", "Xylan", "Mannan", "Galactan", "Arabinan", "Alginate",
            "Fucoidan", "Mannitol", "Protein", "OtherSolids", "Lignin", "Ash",
        ])),
        solids_capture_frac=prA.get("solids_capture_frac", 0.98),
        cake_solids_wt_frac=(
            float(press_cake_solids_wt_frac) if press_cake_solids_wt_frac is not None
            else prA.get("cake_solids_wt_frac", 0.35)
        ),
        solubles_to_pressate_frac=prA.get("solubles_to_pressate_frac", 1.0),
        power_kWh_per_dry_ton_TS=prA.get("power_kWh_per_dry_ton_TS", None),
        ref_capacity_tph_wet=prA.get("ref_capacity_tph_wet", 50.0),
        capex_installed_ref_usd=prA.get("capex_installed_ref_usd", 5e6),
        scale_exponent=prA.get("scale_exponent", 0.6),
        F_BM=prA.get("F_BM", 1.0),
    )

    pb = preprocessing_assumptions.get("pressate_biostimulant", {})
    pcA = pb.get("concentrator", {})
    evA = pb.get("evaporator", {})
    PC = EV = None

    if pb.get("enabled", False) and pb.get("concentrate_pressate", False):
        PC = PressateConcentrator(
            "PC", ins=PR - 1,
            outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
            retained_solute_IDs=tuple(pcA.get(
                "retained_solute_IDs", ["Alginate", "Fucoidan", "Mannitol", "Protein", "OtherSolids"],
            )),
            water_recovery_to_permeate=pcA.get("water_recovery_to_permeate", 0.70),
            retained_solute_recovery_to_concentrate=pcA.get("retained_solute_recovery_to_concentrate", 0.95),
            design_flux_L_m2_h=pcA.get("design_flux_L_m2_h", 35.0),
            operating_pressure_bar=pcA.get("operating_pressure_bar", 5.0),
            electricity_kWh_per_m3_feed=pcA.get("electricity_kWh_per_m3_feed", 0.8),
            capex_usd_per_m2=pcA.get("capex_usd_per_m2", 120.0),
            maintenance_frac_of_capex_per_yr=pcA.get("maintenance_frac_of_capex_per_yr", 0.035),
        )
        if evA.get("enabled", False):
            EV = BiostimulantEvaporator(
                "EV", ins=PC - 0, outs=("biostimulant_product", "biostimulant_vapor"),
                target_solids_wt_frac=evA.get("target_solids_wt_frac", 0.20),
                boiling_temperature_K=evA.get("boiling_temperature_K", 333.15),
                latent_heat_kJ_per_kg=evA.get("latent_heat_kJ_per_kg", 2350.0),
                sensible_cp_kJ_per_kgK=evA.get("sensible_cp_kJ_per_kgK", 4.18),
                electricity_kWh_per_kg_water_evap=evA.get("electricity_kWh_per_kg_water_evap", 0.0),
                nonwater_recovery_to_product=evA.get("nonwater_recovery_to_product", 0.995),
                capex_ref_usd=evA.get("capex_ref_usd", 750000.0),
                ref_evaporation_kgph=evA.get("ref_evaporation_kgph", 10000.0),
                scale_exponent=evA.get("scale_exponent", 0.60),
                maintenance_frac_of_capex_per_yr=evA.get("maintenance_frac_of_capex_per_yr", 0.035),
                F_BM=evA.get("F_BM", 1.0),
            )

    ML = Mill(
        "ML", ins=PR - 0, outs=("milled_biomass", "milling_losses"),
        loss_frac=mlA.get("loss_frac", 0.15),
        power_kWh_per_dry_ton_dry=mlA.get("power_kWh_per_dry_ton_dry", None),
        ref_capacity_dry_ton_per_hr=mlA.get("ref_capacity_dry_ton_per_hr", 10.0),
        purchase_cost_ref_usd=mlA.get("purchase_cost_ref_usd", 206400.0),
        install_factor=mlA.get("install_factor", 1.8),
        scale_exponent=mlA.get("scale_exponent", 0.6),
        F_BM=mlA.get("F_BM", 1.0),
    )

    resolved_pretreatment_case = pretreatment_case or _AD_SHARED.get("methanogenic", {}).get("pretreatment_case", "press_mill_only")
    path_units, streams, units = _build_methanogenic_pathway(
        ML - 0, resolved_pretreatment_case, temperature_regime=temperature_regime,
    )

    if ch4_override is not None:
        units["AD"].ch4_kg_per_kg_vs_fed = float(ch4_override)

    path = [PR]
    if PC is not None:
        path.append(PC)
    if EV is not None:
        path.append(EV)
    path.append(ML)
    path.extend(path_units)

    sys = bst.System("AD_Biomethane_sys", path=path)
    return sys

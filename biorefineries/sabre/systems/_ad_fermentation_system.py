# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Acidogenic-AD-to-microbial-oil system builder for the SaBRe flowsheets.

Public entry point:
- create_ad_fermentation_system(...): the full feedstock-to-product
  system. Wraps create_ad_vfa_system() (raw Sargassum -> acidogenic AD ->
  vfa_broth) followed by _create_vfa_fermentation_system(), and returns
  one assembled bst.System so the whole feedstock-to-oil chain simulates
  (and propagates changes) as a single graph, mirroring how
  create_ad_biomethane_system() is a single self-contained entry point.

Internal helper:
- _create_vfa_fermentation_system(vfa_broth, ...): the fermentation train
  alone (microfilter -> medium conditioning -> Yarrowia lipolytica
  fermenter -> oil recovery). Requires an already-built vfa_broth stream
  as input, so it cannot run standalone -- not part of this module's
  public API. Used by create_ad_fermentation_system() above, and directly
  by systems._integrated_system, which already has a vfa_broth stream
  from a shared/partial preprocessing train and just needs the
  fermentation half.
"""
import biosteam as bst
import flexsolve as flx

from biorefineries.sabre.utils import load_assumptions, non_none
from biorefineries.sabre.units import (
    YarrowiaLipidFermenter,
    OilExtraction,
    VFAMicrofilter,
    FermentationMediumTank,
)
from biorefineries.sabre.systems._ad_vfa_system import create_ad_vfa_system

__all__ = ('create_ad_fermentation_system',)


# Default parameters for _create_vfa_fermentation_system() below
# (data/fermentation.yaml and data/downstream_processing.yaml).
_FERMENTATION_YAML = load_assumptions("fermentation.yaml")
_VFA_FERM = _FERMENTATION_YAML["vfa"]
_VFA_CASE = _VFA_FERM["cases"][_VFA_FERM["case"]]

_DOWNSTREAM_PROCESSING_YAML = load_assumptions("downstream_processing.yaml")
_VFA_DOWNSTREAM = _DOWNSTREAM_PROCESSING_YAML["oil_extraction"]


def _create_vfa_fermentation_system(
    vfa_broth,
    *,
    product_ID: str = _VFA_CASE["product_ID"],
    vfa_IDs: list[str] | None = None,
    conversion: float = None,
    product_yield_kg_per_kg_vfa_consumed: float = None,
    biomass_yield_kg_per_kg_vfa_consumed: float = None,
    co2_yield_kg_per_kg_vfa_consumed: float = None,
    oxygen_kg_per_kg_vfa_consumed: float = None,
    residence_time_h: float = None,
    target_pH: float = None,
    ammonia_dose_kg_per_m3: float = None,
    phosphate_dose_kg_per_m3: float = None,
    base_dose_kg_per_m3: float = None,
    magnesium_sulfate_dose_kg_per_m3: float = None,
    seed_water_kgph: float = _VFA_CASE["seed_water_kgph"],
    seed_cellmass_kgph: float = _VFA_CASE["seed_cellmass_kgph"],

    vfa_to_permeate_frac: float = None,
    water_to_permeate_frac: float = None,
    solids_to_permeate_frac: float = None,
    dissolved_other_to_permeate_frac: float = None,
    microfilter_SEC_kWh_per_m3_feed: float = None,
    microfilter_design_flux_L_m2_h: float = None,

    medium_tank_residence_time_h: float = None,
    medium_tank_mixing_kW_per_m3: float = None,

    target_oil_and_solids_content: float = _VFA_DOWNSTREAM["target_oil_and_solids_content_g_per_L"],
    target_wastewater_concentration: float = 60.0,  # no yaml counterpart
    backend_oil_recovery: float = _VFA_DOWNSTREAM["backend_oil_recovery"],
    backend_oil_water_split: float = _VFA_DOWNSTREAM["backend_oil_water_split"],

    recycle_total_fraction: float = _VFA_DOWNSTREAM["recycle_total_fraction"],
    recycle_cellmass_wt_frac: float = _VFA_DOWNSTREAM["recycle_cellmass_wt_frac"],

    homogenization_kWh_per_kg_dry_biomass: float = None,
):

    chem_ids = set(bst.settings.thermo.chemicals.IDs)
    for req in (product_ID, "CarbonDioxide", "CellMass", "Water"):
        if req not in chem_ids:
            raise RuntimeError(f"Missing required chemical in thermo: {req}")

    if vfa_IDs is None:
        vfa_IDs = list(_VFA_FERM["vfa_IDs"])

    ammonia = bst.Stream("fermentation_ammonia")
    phosphate = bst.Stream("fermentation_phosphate")
    base = bst.Stream("fermentation_base")
    mgso4 = bst.Stream("fermentation_mgso4")
    fresh_seed = bst.Stream(
        "fresh_seed_makeup",
        Water=seed_water_kgph,
        CellMass=seed_cellmass_kgph,
        units="kg/hr",
    )
    recycle_biomass = bst.Stream("recycle_biomass")
    M602 = bst.Mixer("M602", ins=(fresh_seed, recycle_biomass), outs=("seed_to_fermenter",))

    # -------------------------------------------------
    # Upstream conditioning
    # -------------------------------------------------
    MF = VFAMicrofilter(
        "MF",
        ins=vfa_broth,
        outs=("vfa_permeate", "vfa_retentate"),
        vfa_IDs=vfa_IDs,
        **non_none(
            vfa_to_permeate_frac=vfa_to_permeate_frac,
            water_to_permeate_frac=water_to_permeate_frac,
            solids_to_permeate_frac=solids_to_permeate_frac,
            dissolved_other_to_permeate_frac=dissolved_other_to_permeate_frac,
            SEC_kWh_per_m3_feed=microfilter_SEC_kWh_per_m3_feed,
            design_flux_L_m2_h=microfilter_design_flux_L_m2_h,
        ),
    )

    T601 = FermentationMediumTank(
        "T601",
        ins=(MF - 0, ammonia, phosphate, base, mgso4),
        outs=("conditioned_vfa_broth",),
        **non_none(
            ammonia_dose_kg_per_m3=ammonia_dose_kg_per_m3,
            phosphate_dose_kg_per_m3=phosphate_dose_kg_per_m3,
            base_dose_kg_per_m3=base_dose_kg_per_m3,
            magnesium_sulfate_dose_kg_per_m3=magnesium_sulfate_dose_kg_per_m3,
            target_pH=target_pH,
            tau=medium_tank_residence_time_h,
            mixing_kW_per_m3=medium_tank_mixing_kW_per_m3,
        ),
    )

    # -------------------------------------------------
    # Aerated fermenter (Yarrowia lipolytica)
    # -------------------------------------------------
    R601 = YarrowiaLipidFermenter(
        "R601",
        ins=(T601 - 0, M602 - 0),
        outs=("fermentation_vent", "fermentation_broth"),
        vfa_IDs=vfa_IDs,
        product_ID=product_ID,
        V_wf=0.8,
        V_max=150.0,
        kW_per_m3=0.06,
        tau_0=3.0,
        **non_none(
            conversion=conversion,
            product_yield_kg_per_kg_vfa_consumed=product_yield_kg_per_kg_vfa_consumed,
            biomass_yield_kg_per_kg_vfa_consumed=biomass_yield_kg_per_kg_vfa_consumed,
            co2_yield_kg_per_kg_vfa_consumed=co2_yield_kg_per_kg_vfa_consumed,
            oxygen_kg_per_kg_vfa_consumed=oxygen_kg_per_kg_vfa_consumed,
            tau=residence_time_h,
            target_pH=target_pH,
        ),
    )

    # -------------------------------------------------
    # Post-fermentation oil separation
    # -------------------------------------------------
    V605 = bst.MixTank("V605", ins=R601 - 1, outs=("mixed_fermentation_broth",))
    P606 = bst.Pump("P606", ins=V605 - 0, outs=("pumped_fermentation_broth",))

    Ev607 = bst.MultiEffectEvaporator(
        "Ev607",
        ins=P606 - 0,
        outs=("fermentation_concentrate", "evaporator_vapor"),
        P=(101325, 69682, 47057, 30953),
        V=0.90,
        V_definition="First-effect",
        thermo=(R601.outs[1].thermo.ideal()),
        flash=False,
    )
    Ev607.target_oil_and_solids_content = target_oil_and_solids_content
    Ev607.remove_evaporators = True

    P_original = tuple(Ev607.P)
    Pstart = P_original[0]
    Plast = P_original[-1]
    N = len(P_original)

    def concentration_objective(V):
        Ev607.V = V
        Ev607.run()
        effluent = Ev607.outs[0]
        total = effluent.F_mass
        if total <= 0:
            return 0.0
        water = effluent.imass["Water"]
        nonwater_conc = 1000.0 * (1.0 - water / total)
        return Ev607.target_oil_and_solids_content - nonwater_conc

    @Ev607.add_specification(run=False)
    def adjust_evaporation():
        V_last = Ev607.V
        x0 = 0.0
        x1 = 0.5
        Ev607.P = P_original
        Ev607._reload_components = True

        y0 = concentration_objective(x0)
        if y0 <= 0.0:
            Ev607.V = x0
            return
        else:
            Ev607._load_components()
            for i in range(1, N):
                if concentration_objective(1e-6) < 0.0:
                    Ev607.P = tuple(__import__("numpy").linspace(Pstart, Plast, N - 1))
                    Ev607._reload_components = True
                else:
                    break
            y1 = concentration_objective(x1)
            Ev607.V = flx.IQ_interpolation(
                concentration_objective,
                x0, x1, y0, y1,
                x=V_last,
                ytol=1e-5,
                xtol=1e-6,
            )

    P607 = bst.Pump("P607", ins=Ev607 - 0, outs=("pumped_concentrate",), P=101325.0)

    # -------------------------------------------------
    # Cell disruption + oil extraction placeholder
    # Slots between P607 and C603_2.
    # Pass-through unit: carries capital + electricity costs for
    # high-pressure homogenization and lipid extraction.
    # Capital anchor: NREL/TP-5100-55431 (2012), Davis et al.
    # -------------------------------------------------
    OE = OilExtraction(
        "OE",
        ins=P607 - 0,
        outs=("extracted_broth",),
        product_ID=product_ID,
        cellmass_ID="CellMass",
        **non_none(homogenization_kWh_per_kg_dry_biomass=homogenization_kWh_per_kg_dry_biomass),
    )

    C603_2 = bst.LiquidsSplitCentrifuge(
        "C603_2",
        ins=OE - 0,   # now takes from OE, not P607
        outs=("backend_oil", "cellmass_plus_aqueous"),
        split={product_ID: backend_oil_recovery, "Water": backend_oil_water_split},
    )

    S602 = bst.MockSplitter(
        "S602",
        ins=C603_2 - 1,
        outs=(recycle_biomass, "residual_purge"),
    )
    S602.recycle_total_fraction = recycle_total_fraction
    S602.recycle_cellmass_wt_frac = recycle_cellmass_wt_frac

    @S602.add_specification(run=True)
    def adjust_biomass_recycle():
        feed = S602.ins[0]
        recycle, purge = S602.outs

        recycle.empty()
        purge.copy_like(feed)

        cellmass_available = float(feed.imass["CellMass"]) if "CellMass" in feed.chemicals else 0.0
        water_available = float(feed.imass["Water"]) if "Water" in feed.chemicals else 0.0

        if feed.F_mass <= 0 or cellmass_available <= 0:
            recycle.empty()
            purge.copy_like(feed)
            return

        target_recycle_mass = S602.recycle_total_fraction * feed.F_mass
        target_cellmass = target_recycle_mass * S602.recycle_cellmass_wt_frac
        cellmass_recycle = min(cellmass_available, target_cellmass)
        water_recycle = cellmass_recycle * (1.0 / S602.recycle_cellmass_wt_frac - 1.0)
        water_recycle = min(water_available, water_recycle)

        recycle.imass["CellMass"] = cellmass_recycle
        recycle.imass["Water"] = water_recycle
        purge.imass["CellMass"] -= cellmass_recycle
        purge.imass["Water"] -= water_recycle
        purge.mol.remove_negatives()
        recycle.T = purge.T = feed.T

    S601 = bst.Splitter(
        "S601",
        ins=Ev607 - 1,
        outs=("condensate_to_wastewater", "evaporator_condensate"),
        split=0.5,
    )

    M601 = bst.Mixer(
        "M601",
        ins=(S601 - 0, S602 - 1),
        outs=("fermentation_wastewater",),
    )

    M601.target_wastewater_concentration = target_wastewater_concentration

    @M601.add_specification(run=True, impacted_units=[S601])
    def adjust_wastewater_concentration():
        concentrated_wastewater = M601.ins[1]
        waste = concentrated_wastewater.F_mass - concentrated_wastewater.imass["Water"]
        if concentrated_wastewater.F_vol <= 0 or waste <= 0:
            S601.split[:] = 0.0
            return

        current_concentration = waste / concentrated_wastewater.F_vol
        required_water = (
            (1.0 / M601.target_wastewater_concentration) - (1.0 / current_concentration)
        ) * waste * 1000.0

        F_mass = S601.ins[0].F_mass
        if F_mass > 0:
            split = required_water / F_mass
            split = min(max(split, 0.0), 1.0)
            S601.split[:] = split

    # OE added to system path between P607 and C603_2
    sys = bst.System(
        "VFA_FER_sys",
        path=(MF, T601, M602, R601, V605, P606, Ev607, P607, OE, C603_2, S602, S601, M601)
    )

    key_streams = {
        "vfa_broth": vfa_broth,
        "vfa_permeate": MF.outs[0],
        "vfa_retentate": MF.outs[1],
        "conditioned_vfa_broth": T601.outs[0],
        "fresh_seed_makeup": fresh_seed,
        "recycle_biomass": recycle_biomass,
        "seed": M602.outs[0],
        "vent": R601.outs[0],
        "fermentation_broth": R601.outs[1],
        "extracted_broth": OE.outs[0],
        "backend_oil": C603_2.outs[0],
        "residual_slurry": C603_2.outs[1],
        "residual_purge": S602.outs[1],
        "evaporator_condensate": S601.outs[1],
        "fermentation_wastewater": M601.outs[0],
        "ammonia": ammonia,
        "phosphate": phosphate,
        "base": base,
        "mgso4": mgso4,
    }

    units = {
        "microfilter": MF,
        "medium_tank": T601,
        "seed_mixer": M602,
        "fermenter": R601,
        "post_mix_tank": V605,
        "post_feed_pump": P606,
        "evaporator": Ev607,
        "backend_oil_pump": P607,
        "oil_extraction": OE,
        "oil_centrifuge": C603_2,
        "biomass_recycle_splitter": S602,
        "condensate_splitter": S601,
        "wastewater_mixer": M601,
        # legacy aliases
        "centrifuge": C603_2,
        "lipid_recovery": C603_2,
    }

    return sys, key_streams, units


def create_ad_fermentation_system(
    feedstock_type: str = "pelagic",
    milled_biomass_stream=None,
    enable_heat_shock: bool = False,
    hs_target_temperature_K: float = 338.15,
    hs_events_per_day: float = 1.0 / 7.0,
    hs_heated_fraction_of_liquid: float = 0.10,
    hs_duration_min: float = 15.0,
    temperature_regime: str = "mesophilic",
    **fermentation_kwargs,
):
    """
    Build the full feedstock-to-product system: raw Sargassum (or an
    already-milled biomass stream) -> acidogenic AD -> VFA broth -> VFA
    fermentation -> microbial oil.

    Parameters
    ----------
    feedstock_type, milled_biomass_stream, enable_heat_shock,
    hs_target_temperature_K, hs_events_per_day,
    hs_heated_fraction_of_liquid, hs_duration_min, temperature_regime
        Forwarded to create_ad_vfa_system().
    **fermentation_kwargs
        Forwarded to _create_vfa_fermentation_system() (everything except
        vfa_broth, which is supplied internally from the AD-VFA
        subsystem's output).

    Returns
    -------
    sys : bst.System
        The full feedstock -> fermentation-product system.
    streams : dict
        Key streams, including the raw feedstock ('feed'), the AD-VFA
        subsystem's 'vfa_broth', and all of
        _create_vfa_fermentation_system()'s streams (final product in
        'backend_oil').
    units : dict
        Key units from both subsystems.
    """
    ad_vfa_sys = create_ad_vfa_system(
        feedstock_type=feedstock_type,
        milled_biomass_stream=milled_biomass_stream,
        enable_heat_shock=enable_heat_shock,
        hs_target_temperature_K=hs_target_temperature_K,
        hs_events_per_day=hs_events_per_day,
        hs_heated_fraction_of_liquid=hs_heated_fraction_of_liquid,
        hs_duration_min=hs_duration_min,
        temperature_regime=temperature_regime,
    )
    feed = ad_vfa_sys.feeds[0]
    vfa_broth = ad_vfa_sys.flowsheet.stream.vfa_broth

    fer_sys, fer_streams, fer_units = _create_vfa_fermentation_system(
        vfa_broth=vfa_broth, **fermentation_kwargs,
    )

    sys = bst.System.from_units(
        "AD_Fermentation_sys",
        units=list(ad_vfa_sys.units) + list(fer_sys.units),
    )

    streams = {"feed": feed, "vfa_broth": vfa_broth, **fer_streams}
    units = {u.ID: u for u in ad_vfa_sys.units}
    units.update(fer_units)

    return sys, streams, units

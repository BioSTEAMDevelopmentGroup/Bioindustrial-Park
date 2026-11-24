# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Sargassum-to-dried-Ca(3HP)2 system builder for the SaBRe flowsheets.

Pathway, following the collaborators' process:
    raw Sargassum -> Mill (MI) -> Press (PR: juice + pulp)
      juice -> PFS/PC/DIL/EV (membrane + evaporator concentration)
      pulp  -> [baseline scenario only: mild thermochemical pretreatment
               (1% H2SO4, 1 M NaCl, 121 degC, 20 min; Chundawat) in
               ACID_ST/AM/M_PT/R_PT/F_PT/H_PT, and ammonia neutralization
               (NH3_ST/NT)]
            -> M_SA (enzyme) -> H_SA -> SA (96 h hydrolysis tank) -> SLS
               (residue to the boiler) -> HPFS/HPC/HDIL/HEV (concentration)
      concentrated juice + concentrated hydrolysate -> FM -> HPFermentation
      -> the shared downstream below (cell mass to the boiler as before).
    Two scenarios: `baseline` (with the mild thermochemical pretreatment) and
    `pretreatment_press` (no pretreatment reactor; the same 1 wt% H2SO4 and 1 M NaCl
    dosed into the pre-press biomass, assumed to give the same hydrolysis as the baseline
    pretreatment, then enzymatic hydrolysis with ammonia neutralization).
    The hydrolysate is concentrated as well as the juice, because at 20 wt%
    insoluble solids it carries only about 4 wt% sugar.

Downstream of the fermentation feed:
      -> HPFermentation (Ca(OH)2 dosed in situ; Ca3HP2 already formed by
         the time broth leaves this unit)
      -> S401 (SolidsCentrifuge: removes CellMass)
      -> M_ML (mixes in recycled mother liquor)
      -> broth_evaporator (MultiEffectEvaporator: concentrates to the
         patent's ~300-600 g/L crystallization threshold)
      -> CaHPCrystallizer -> S402 (CrystalCentrifuge: separates crystal
         cake, which entrains some mother liquor and so its impurities,
         from the mother liquor)
      -> SP_ML (splits mother liquor: `mother_liquor_recycle` back to M_ML,
         data/3hp.yaml mother_liquor_recycle.recycle_frac; remainder
         `mother_liquor_purge` -> disposal_wastewater)
      -> DrumDryer -> dried_product (Ca3HP2 + entrained impurities)
      -> PS (MockSplitter, accounting only: `ca3hp2_product` = the pure
         Ca3HP2 content, priced and used for MSP; `product_impurities` =
         the rest, unpriced)
    (M_ML ... SP_ML form a recycle loop, built as the `ml_loop` subsystem.)

    milling_losses, pressed_cake, cell_mass (all combustible biomass) --
    mixed together (M_BT) -> BoilerTurbogenerator (BT, a facility, not in
    the main path above) for steam/electricity credit -- NOT priced as
    disposal. This is the one place this system departs from every other
    sabre system's own convention of pricing solid/liquid waste directly
    rather than modeling a real BT/WWT facility (no other sabre system has
    a BT), per explicit user instruction to route both press cake and
    cell mass to the boiler.

Chemical-set note: this pathway needs sabre._chemicals.create_chemicals's
include_hp3=True chemical superset (Glucose, Uronics, Enzyme, HP,
CalciumDihydroxide, Ca3HP2), which every other sabre system deliberately
excludes to avoid a small (~0.08%) NumPy-pairwise-summation-driven drift on
their own golden-value regression tests (see create_chemicals's own
docstring). Because thermosteam's chemical set is process-global, this
system must not be built in the same Python process as any of the other
sabre systems -- create_3hp_system() unconditionally (re)builds the
superset every time it's called, which would silently corrupt the other
five if they ran afterward in that same process. This is also why
create_3hp_system() is not wired into sabre.__init__.load()'s dispatcher.
"""
from pathlib import Path

import biosteam as bst
import flexsolve as flx
import numpy as np
import thermosteam as tmo

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
    get_solids_group_IDs,
)
from biorefineries.sabre.units import (
    Mill, Press, PressateConcentrator, BiostimulantEvaporator,
    HPFermentation, CaHPCrystallizer, CrystalCentrifuge, HydrolysisTank,
    DiluteAcidReactor,
)
from biorefineries.sabre._tea import create_tea

__all__ = ('create_3hp_system', 'price_3hp_system', 'SCENARIOS')

_TEA_PRICE = load_assumptions("tea.yaml")["price"]
_3HP_YAML = load_assumptions("3hp.yaml")
_HYDROLYSATE_SEPARATOR = _3HP_YAML["hydrolysate_separator"]
_PULP_HYDROLYSIS = _3HP_YAML["pulp_enzymatic_hydrolysis"]
_MILD_PRETREATMENT = _3HP_YAML["mild_pretreatment"]
SCENARIOS = ("baseline", "pretreatment_press")
_PRETREATMENT_PRESS = _3HP_YAML["pretreatment_press"]

_PRESS_CAKE_SOLIDS_WT_FRAC_OVERRIDE = _3HP_YAML["press_cake_solids_wt_frac_override"]
_PRESSATE_CONCENTRATOR = _3HP_YAML["pressate_concentrator"]
_FEED_EVAPORATOR = _3HP_YAML["fermentation_feed_evaporator"]
_BROTH_EVAPORATOR = _3HP_YAML["broth_evaporator"]
_CELL_MASS_CENTRIFUGE = _3HP_YAML["cell_mass_centrifuge"]
_CRYSTALLIZATION = _3HP_YAML["crystallization"]
_CRYSTAL_SEPARATOR = _3HP_YAML["crystal_separator"]
_ML_RECYCLE = _3HP_YAML["mother_liquor_recycle"]
_DRYER = _3HP_YAML["dryer"]

_PRETREATMENT_PRESS_DISSOLVED_IDS = ("Glucose", "Uronics", "H2SO4", "NaCl")
# Database solids that behave as dissolved salts (not insoluble solids).
_DISSOLVED_SALT_IDS = ("AmmoniumSulfate", "NaCl")
_NOT_HEAT_INTEGRATED = frozenset({"M_PR", "NJ", "H_SA", "M_PT", "R_PT", "F_PT", "H_PT", "AM", "NT"})

FRESH_WATER_PRICE_USD_PER_KG = bst.stream_utility_prices['Process water']

# PressateConcentrator/BiostimulantEvaporator settings of the 3-HP sugar
# concentration chain (data/3hp.yaml), applied to both the juice and
# the hydrolysate.
_SUGAR_PC_KWARGS = dict(
    retained_solute_IDs=tuple(_PRESSATE_CONCENTRATOR["retained_solute_IDs"]),
    water_recovery_to_permeate=_PRESSATE_CONCENTRATOR["water_recovery_to_permeate"],
    retained_solute_recovery_to_concentrate=_PRESSATE_CONCENTRATOR["retained_solute_recovery_to_concentrate"],
    nontarget_solute_recovery_to_permeate=_PRESSATE_CONCENTRATOR["nontarget_solute_recovery_to_permeate"],
    design_flux_L_m2_h=_PRESSATE_CONCENTRATOR["design_flux_L_m2_h"],
    operating_pressure_bar=_PRESSATE_CONCENTRATOR["operating_pressure_bar"],
    electricity_kWh_per_m3_feed=_PRESSATE_CONCENTRATOR["electricity_kWh_per_m3_feed"],
    capex_usd_per_m2=_PRESSATE_CONCENTRATOR["capex_usd_per_m2"],
    maintenance_frac_of_capex_per_yr=_PRESSATE_CONCENTRATOR["maintenance_frac_of_capex_per_yr"],
)
_SUGAR_EV_KWARGS = dict(
    target_solids_wt_frac=_FEED_EVAPORATOR["target_solids_wt_frac"],
    boiling_temperature_K=_FEED_EVAPORATOR["boiling_temperature_K"],
    condensing_temperature_K=_FEED_EVAPORATOR["condensing_temperature_K"],
    nonwater_recovery_to_product=_FEED_EVAPORATOR["nonwater_recovery_to_product"],
    tau=_FEED_EVAPORATOR["tau"],
    vacuum_system_preference=_FEED_EVAPORATOR["vacuum_system_preference"],
)


def _make_pressate_concentration_chain(
    feed, *, unit_prefix, stream_prefix, fresh_water_ID, product_ID, condensate_ID,
    pc_kwargs, ev_kwargs,
):
    """
    PFS -> PressateConcentrator -> DIL -> BiostimulantEvaporator: the same
    split-and-bypass concentration chain
    systems._biostimulant_system.create_biostimulant_system builds, factored
    out here so it can be built twice ('' prefix for the juice, 'H' for the
    pulp hydrolysate). PFS's
    specification splits the feed into only the fraction that needs
    concentrating and a raw bypass, drawing fresh water if even bypassing
    everything is still above the EV target -- see
    create_biostimulant_system for the derivation.

    Returns (PFS, PC, DIL, EV).
    """
    fresh_water = bst.Stream(fresh_water_ID, Water=0.0, units="kg/hr")
    fresh_water.price = FRESH_WATER_PRICE_USD_PER_KG

    PFS = bst.MockSplitter(
        f"{unit_prefix}PFS", ins=feed,
        outs=(f"{stream_prefix}to_concentrator", f"{stream_prefix}bypass_pressate"),
    )
    PC = PressateConcentrator(
        f"{unit_prefix}PC", ins=PFS - 0,
        outs=(f"{stream_prefix}membrane_concentrate", f"{stream_prefix}permeate"),
        **pc_kwargs,
    )
    PC.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]
    DIL = bst.Mixer(
        f"{unit_prefix}DIL", ins=(PC - 0, PFS - 1, fresh_water),
        outs=(f"{stream_prefix}diluted_pressate",),
    )
    EV = BiostimulantEvaporator(
        f"{unit_prefix}EV", ins=DIL - 0, outs=(product_ID, condensate_ID), **ev_kwargs,
    )
    EV.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

    @PFS.add_specification(run=True)
    def adjust_pressate_split():
        pressate = PFS.ins[0]
        to_conc, bypass = PFS.outs

        fresh_water.empty()

        F = float(pressate.F_mass)
        ids = set(pressate.chemicals.IDs)
        water_raw = float(pressate.imass["Water"]) if "Water" in ids else 0.0
        nonwater_raw = F - water_raw

        target = EV.target_solids_wt_frac

        if (
            target is None or not (0.0 < target < 1.0)
            or nonwater_raw <= 0 or not PC.enabled
        ):
            to_conc.copy_like(pressate)
            bypass.empty()
            return

        retained = set(PC.retained_solute_IDs)
        conc_water = water_raw * (1.0 - PC.water_recovery_to_permeate)
        conc_nonwater = 0.0
        for cid in pressate.chemicals.IDs:
            if cid == "Water":
                continue
            m = float(pressate.imass[cid])
            if m <= 0:
                continue
            rec = (
                PC.retained_solute_recovery_to_concentrate if cid in retained
                else (1.0 - PC.nontarget_solute_recovery_to_permeate)
            )
            conc_nonwater += m * rec
        conc_mass_full = conc_water + conc_nonwater

        if conc_mass_full <= 0:
            to_conc.copy_like(pressate)
            bypass.empty()
            return

        x_raw = nonwater_raw / F
        x_conc = conc_nonwater / conc_mass_full

        if target <= x_raw:
            f = 0.0
            water_needed_total = nonwater_raw * (1.0 - target) / target
            fresh_water.imass["Water"] = max(water_needed_total - water_raw, 0.0)
        elif target >= x_conc:
            f = 1.0
        else:
            A = conc_mass_full * (x_conc - target)
            B = F * (x_raw - target)
            f = B / (B - A)

        to_conc.copy_like(pressate)
        to_conc.mol *= f
        bypass.copy_like(pressate)
        bypass.mol *= (1.0 - f)

    return PFS, PC, DIL, EV


def _make_pulp_pretreatment_train(cake):
    """
    Mild thermochemical pretreatment (1% H2SO4, 1 M NaCl, 121 degC, 20 min)
    and neutralization of the press pulp (the `baseline` scenario), placed
    ahead of the enzymatic hydrolysis:
    acid dosing -> steam mixer -> pressure reactor -> flash -> condensate,
    flashed slurry -> ammonia neutralization. Mirrors biorefineries.cellulosic's
    corn stover train (systems/pretreatment/dilute_acid.py) with sabre-local
    units; see data/3hp.yaml (`mild_pretreatment`).

    Acid and ammonia are dosed by specifications on their storage tanks,
    so the fresh-feed streams (the ones priced) carry the actual demand.
    Steam for the steam mixer is booked as a heat utility by
    bst.SteamMixer itself, so there is no priced steam stream.

    Returns (units, neutralized_slurry): the units in path order and the
    neutralized slurry outlet, which feeds the enzymatic hydrolysis.
    """
    A = _MILD_PRETREATMENT
    price_wastewater = _TEA_PRICE["disposal_wastewater"]["baseline"]

    def dry_pulp_kgph():
        return cake.F_mass - cake.imass["Water"]

    # --- Acid dosing ------------------------------------------------------
    sulfuric_acid = bst.Stream(
        "sulfuric_acid",
        H2SO4=A["acid_solution_H2SO4_kg"], Water=A["acid_solution_water_kg"], units="kg/hr",
    )
    sulfuric_acid.price = _TEA_PRICE["sulfuric_acid"]["baseline"]
    ACID_ST = bst.StorageTank("ACID_ST", ins=sulfuric_acid, outs=("acid_from_storage",))

    warm_water_1 = bst.Stream(
        "warm_water_1", Water=0.0, units="kg/hr",
        T=A["warm_water_T_K"], P=A["warm_water_P_atm"] * 101325,
    )
    warm_water_1.price = FRESH_WATER_PRICE_USD_PER_KG
    # NaCl of the mild thermochemical pretreatment (1 M in the slurry liquid),
    # dosed on dry pulp like the acid; bought fresh, no recycle.
    sodium_chloride = bst.Stream("sodium_chloride", NaCl=0.0, units="kg/hr", phase="s")
    sodium_chloride.price = _TEA_PRICE["sodium_chloride"]["baseline"]
    AM = bst.Mixer("AM", ins=(warm_water_1, ACID_ST - 0, sodium_chloride), outs=("acid_liquor",))

    @ACID_ST.add_specification(run=True)
    def dose_acid():
        solution_kgph = A["acid_solution_loading_kg_per_kg_dry"] * dry_pulp_kgph()
        total_ref = A["acid_solution_H2SO4_kg"] + A["acid_solution_water_kg"]
        sulfuric_acid.imass["H2SO4"] = solution_kgph * A["acid_solution_H2SO4_kg"] / total_ref
        sulfuric_acid.imass["Water"] = solution_kgph * A["acid_solution_water_kg"] / total_ref

    @AM.add_specification(run=True)
    def dose_warm_water():
        warm_water_1.imass["Water"] = A["warm_water_kg_per_kg_dry"] * dry_pulp_kgph()
        sodium_chloride.imass["NaCl"] = A["nacl_kg_per_kg_dry"] * dry_pulp_kgph()

    # --- Steam mixer -> pressure reactor -> flash -> condenser --------------
    warm_water_2 = bst.Stream(
        "warm_water_2", Water=0.0, units="kg/hr",
        T=A["warm_water_T_K"], P=A["warm_water_P_atm"] * 101325,
    )
    warm_water_2.price = FRESH_WATER_PRICE_USD_PER_KG
    pretreatment_steam = bst.Stream("pretreatment_steam", phase="g")
    M_PT = bst.SteamMixer(
        "M_PT", ins=(cake, pretreatment_steam, warm_water_2, AM - 0),
        outs=("steam_mixed_slurry",),
        P=A["reactor_P_atm"] * 101325, T=A["reactor_T_K"],
        solids_loading=A["solids_loading"],
    )
    R_PT = DiluteAcidReactor("R_PT", ins=M_PT - 0, outs=("pretreated_slurry",))
    F_PT = bst.Flash(
        "F_PT", ins=R_PT - 0, outs=("pretreatment_flash_vapor", "flashed_slurry"),
        P=A["flash_P_Pa"], Q=0,
    )
    # Saturated liquid out (V=0) rather than a fixed T: the 1 atm flash vapor
    # leaves at water's boiling point, so a hard-coded 373.15 K would ask for
    # an outlet hotter than the inlet.
    H_PT = bst.HXutility(
        "H_PT", ins=F_PT - 0, outs=("pretreatment_condensate",), V=0,
    )
    H_PT.outs[0].price = price_wastewater

    # --- Neutralization: 2 NH3 + H2SO4 -> (NH4)2SO4 -------------------------
    ammonia = bst.Stream("ammonia", Ammonia=1.0, units="kg/hr", phase="l")
    ammonia.price = _TEA_PRICE["ammonia"]["baseline"]
    NH3_ST = bst.StorageTank("NH3_ST", ins=ammonia, outs=("ammonia_from_storage",))
    NT = bst.MixTank("NT", ins=(F_PT - 1, NH3_ST - 0), outs=("neutralized_slurry",))

    @NH3_ST.add_specification(run=True)
    def dose_ammonia():
        ammonia.imol["Ammonia"] = A["ammonia_mol_per_mol_H2SO4"] * F_PT.outs[1].imol["H2SO4"]

    neutralization_rxn = tmo.Reaction(
        "2 Ammonia + H2SO4 -> AmmoniumSulfate", reactant="H2SO4", X=1.0,
    )

    @NT.add_specification
    def neutralize():
        NT._run()
        neutralization_rxn.adiabatic_reaction(NT.outs[0])

    units = [ACID_ST, AM, M_PT, R_PT, F_PT, H_PT, NH3_ST, NT]
    return units, NT - 0


def _build_juice_pulp_upstream(MI, press_params=None, scenario="baseline"):
    """
    Upstream section, following the
    collaborators' process description: mill (MI) -> plain screw press (PR,
    juice = `pressate` and `pressed_cake` = pulp) -> juice concentration
    (PFS/PC/DIL/EV); in the `baseline` scenario the pulp is first
    given the mild thermochemical pretreatment and neutralized (`_make_pulp_pretreatment_train`),
    in `pretreatment_press` the same acid and NaCl are dosed ahead of the screw press instead; the pulp goes through an enzymatic hydrolysis tank (SA,
    96 h) and a pressure filter (SLS, residual solids to the boiler); the
    concentrated juice and the pulp hydrolysate are blended (FM) into the
    single fermentation feed.

    Scenario parameters live on the units: `M_SA.enzyme_loading_mg_per_g_dwb`
    (total enzyme per g dry biomass fed to the press), `SA.hydrolysis_conversion` and (with
    `SA.hydrolysis_rxns.X`) the monomer release, and `FM.target_solids_wt_frac`
    (dissolved-solids wt frac of the fermentation feed).

    Returns a dict with 'units' (path order, MI included),
    'fermentation_feed' and 'boiler_solids'.
    """
    # A feedstock type with a measured press split (feedstock.yaml `press`)
    # sets the cake solids and the per-component cake fractions.
    press_params = press_params or {}
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    A = _MILD_PRETREATMENT
    cake_fraction = dict(press_params.get("cake_fraction") or {})
    press_feed = MI - 0
    press_units = []
    # Neutralization of dissolved H2SO4 with ammonia: 2 NH3 + H2SO4 -> (NH4)2SO4.
    neutralization_rxn = tmo.Reaction(
        "2 Ammonia + H2SO4 -> AmmoniumSulfate", reactant="H2SO4", X=1.0)
    if scenario == "pretreatment_press":
        # The mild thermochemical pretreatment of the baseline (1 wt% H2SO4 and
        # 1 M NaCl in the liquid, same basis) applied to the milled biomass
        # ahead of the screw press: concentrated H2SO4 (the same ~93.3 wt%
        # solution as the baseline acid) so that H2SO4 is `H2SO4_wt_frac_of_liquid`
        # of the liquid of the pre-press mixture, and NaCl at
        # `NaCl_kg_per_kg_liquid_water`. No heating or hold time is modeled in
        # the press. The press is assumed to give the baseline pretreatment's
        # hydrolysis (same Glucan and Alginate conversions), applied to the
        # mixture. The monomers, acid and NaCl are dissolved, so they follow the
        # water between cake and pressate (set below); Mannitol keeps the
        # measured split.
        for cid in _PRETREATMENT_PRESS_DISSOLVED_IDS:
            cake_fraction[cid] = 0.0
        press_acid = bst.Stream("press_sulfuric_acid", units="kg/hr")
        press_acid.price = _TEA_PRICE["sulfuric_acid"]["baseline"]
        press_salt = bst.Stream("press_sodium_chloride", NaCl=0.0, units="kg/hr", phase="s")
        press_salt.price = _TEA_PRICE["sodium_chloride"]["baseline"]
        M_PR = bst.Mixer(
            "M_PR", ins=(MI - 0, press_acid, press_salt), outs=("biomass_with_acid",))
        M_PR.H2SO4_wt_frac_of_liquid = _PRETREATMENT_PRESS["H2SO4_wt_frac_of_liquid"]
        M_PR.NaCl_kg_per_kg_liquid_water = _PRETREATMENT_PRESS["NaCl_kg_per_kg_liquid_water"]
        press_hydrolysis_rxns = tmo.ParallelReaction([
            tmo.Reaction("Glucan + Water -> Glucose", reactant="Glucan",
                         X=A["glucan_to_glucose_conversion"]),
            tmo.Reaction("Alginate + Water -> Uronics", reactant="Alginate",
                         X=A["alginate_to_uronics_conversion"]),
        ])

        @M_PR.add_specification
        def dose_acid_and_hydrolyze():
            biomass = MI.outs[0]
            solution_H2SO4 = A["acid_solution_H2SO4_kg"] / (
                A["acid_solution_H2SO4_kg"] + A["acid_solution_water_kg"])
            # solution mass a: H2SO4 = x * liquid, with liquid = biomass water
            # + the water of the acid solution: solution_H2SO4 * a =
            # x * (W + (1 - solution_H2SO4) * a)
            x = M_PR.H2SO4_wt_frac_of_liquid
            W = float(biomass.imass["Water"])
            a = x * W / (solution_H2SO4 - x * (1.0 - solution_H2SO4))
            press_acid.empty()
            press_acid.imass["H2SO4"] = a * solution_H2SO4
            press_acid.imass["Water"] = a * (1.0 - solution_H2SO4)
            press_salt.empty()
            press_salt.phase = "s"
            press_salt.imass["NaCl"] = M_PR.NaCl_kg_per_kg_liquid_water * (W + a * (1.0 - solution_H2SO4))
            M_PR._run()
            press_hydrolysis_rxns(M_PR.outs[0])

        press_feed = M_PR - 0
        press_units = [M_PR]
    PR = Press(
        "PR", ins=press_feed, outs=("pressed_cake", "pressate"),
        cake_solids_wt_frac=press_params.get(
            "cake_solids_wt_frac", _PRESS_CAKE_SOLIDS_WT_FRAC_OVERRIDE),
        component_cake_frac=cake_fraction or None,
    )
    if scenario == "pretreatment_press":
        @PR.add_specification
        def dissolved_follow_water():
            PR._run()
            cake_out, pressate_out = PR.outs
            water = float(cake_out.imass["Water"] + pressate_out.imass["Water"])
            to_cake_water = float(cake_out.imass["Water"]) / water if water > 0 else 0.0
            for cid in _PRETREATMENT_PRESS_DISSOLVED_IDS:
                total = float(cake_out.imass[cid] + pressate_out.imass[cid])
                cake_out.imass[cid] = total * to_cake_water
                pressate_out.imass[cid] = total * (1.0 - to_cake_water)
    S = _PULP_HYDROLYSIS
    H = _HYDROLYSATE_SEPARATOR
    cake = PR - 0
    pretreatment_units = []
    if scenario == "baseline":
        pretreatment_units, cake = _make_pulp_pretreatment_train(PR - 0)

    # pretreatment_press: the juice carries part of the acid and is neutralized before
    # its concentration (the cake's acid is neutralized in the hydrolysis, M_SA).
    juice = PR - 1
    juice_units = []
    if scenario == "pretreatment_press":
        juice_ammonia = bst.Stream("juice_ammonia", Ammonia=1.0, units="kg/hr", phase="l")
        juice_ammonia.price = _TEA_PRICE["ammonia"]["baseline"]
        NJ = bst.MixTank("NJ", ins=(PR - 1, juice_ammonia), outs=("neutralized_pressate",))

        @NJ.add_specification
        def neutralize_juice():
            juice_ammonia.empty()
            juice_ammonia.phase = "l"
            juice_ammonia.imol["Ammonia"] = (
                A["ammonia_mol_per_mol_H2SO4"] * float(NJ.ins[0].imol["H2SO4"]))
            NJ._run()
            neutralization_rxn.adiabatic_reaction(NJ.outs[0])

        juice = NJ - 0
        juice_units = [NJ]

    enzyme = bst.Stream("hydrolysis_enzyme", units="kg/hr")
    enzyme.price = _TEA_PRICE["enzyme"]["baseline"]
    hydrolysis_water = bst.Stream("hydrolysis_water", Water=0.0, units="kg/hr")
    hydrolysis_water.price = FRESH_WATER_PRICE_USD_PER_KG
    # pretreatment_press: ammonia for neutralizing the cake's acid joins the hydrolysis mixer.
    hydrolysis_ammonia = bst.Stream("hydrolysis_ammonia", Ammonia=0.0, units="kg/hr", phase="l")
    hydrolysis_ammonia.price = _TEA_PRICE["ammonia"]["baseline"]
    M_SA = bst.Mixer(
        "M_SA",
        ins=(cake, enzyme, hydrolysis_water, *((hydrolysis_ammonia,) if scenario == "pretreatment_press" else ())),
        outs=("pulp_with_enzyme",),
    )
    M_SA.enzyme_loading_mg_per_g_dwb = _PULP_HYDROLYSIS["enzyme_loading_mg_per_g_dwb"]
    insoluble_IDs = [i for i in get_solids_group_IDs() if i not in _DISSOLVED_SALT_IDS]

    @M_SA.add_specification
    def dose_enzyme_and_water():
        pulp = M_SA.ins[0]
        enzyme.empty()
        enzyme.phase = "l"
        hydrolysis_water.empty()
        hydrolysis_ammonia.empty()
        hydrolysis_ammonia.phase = "l"
        if scenario == "pretreatment_press":
            hydrolysis_ammonia.imol["Ammonia"] = (
                A["ammonia_mol_per_mol_H2SO4"] * float(pulp.imol["H2SO4"]))
        biomass = MI.outs[0]  # milled biomass fed to the press
        dry_kgph = float(biomass.F_mass - biomass.imass["Water"])
        enzyme.imass["Enzyme"] = M_SA.enzyme_loading_mg_per_g_dwb / 1000.0 * dry_kgph
        insoluble_kgph = float(pulp.imass[insoluble_IDs].sum())
        total_kgph = pulp.F_mass + enzyme.F_mass + hydrolysis_ammonia.F_mass
        hydrolysis_water.imass["Water"] = max(insoluble_kgph / S["solids_loading"] - total_kgph, 0.0)
        M_SA._run()
        if scenario == "pretreatment_press":
            neutralization_rxn.adiabatic_reaction(M_SA.outs[0])

    H_SA = bst.HXutility("H_SA", ins=M_SA - 0, outs=("pulp_to_hydrolysis",), T=S["T_K"])
    SA = HydrolysisTank(
        "SA", ins=H_SA - 0, outs=("hydrolyzed_pulp",), tau=S["tau_hr"], V_max=S["V_max_m3"],
    )
    SA.hydrolysis_conversion = _PULP_HYDROLYSIS["hydrolysis_conversion"]
    SA.hydrolysis_rxns = tmo.ParallelReaction([
        tmo.Reaction("Glucan + Water -> Glucose", reactant="Glucan", X=SA.hydrolysis_conversion),
        tmo.Reaction("Alginate + Water -> Uronics", reactant="Alginate", X=SA.hydrolysis_conversion),
    ])

    @SA.add_specification
    def hydrolyze():
        SA._run()
        SA.hydrolysis_rxns.X[:] = SA.hydrolysis_conversion
        SA.hydrolysis_rxns(SA.outs[0])

    split = {
        cid: (H["solids_split_to_cake"] if cid in insoluble_IDs else H["solutes_split_to_cake"])
        for cid in SA.chemicals.IDs if cid != "Water"
    }
    SLS = bst.PressureFilter(
        "SLS", ins=SA - 0, outs=("hydrolysis_residue", "pulp_hydrolysate"),
        split=split, moisture_content=H["moisture_content"],
    )

    # Juice concentration chain (PFS/PC/DIL/EV) on the pressate.
    PFS, PC, DIL, EV = _make_pressate_concentration_chain(
        juice, unit_prefix="", stream_prefix="",
        fresh_water_ID="fermentation_fresh_water",
        product_ID="concentrated_juice", condensate_ID="EV_condensed_vapor",
        pc_kwargs=_SUGAR_PC_KWARGS, ev_kwargs=_SUGAR_EV_KWARGS,
    )
    # The pulp hydrolysate carries only a few wt% sugar (the hydrolysis slurry is
    # 20 wt% insoluble solids), and boiling off that water would dominate the
    # plant's heat demand, so it goes through the same membrane/evaporator
    # concentration chain ('H' prefix) before blending. Its evaporator target
    # follows the juice evaporator's, the one knob the scenario scripts set.
    HPFS, HPC, HDIL, HEV = _make_pressate_concentration_chain(
        SLS - 1, unit_prefix="H", stream_prefix="hydrolysate_",
        fresh_water_ID="hydrolysate_fresh_water",
        product_ID="concentrated_hydrolysate",
        condensate_ID="HEV_condensed_vapor",
        pc_kwargs=_SUGAR_PC_KWARGS, ev_kwargs=_SUGAR_EV_KWARGS,
    )

    @SLS.add_specification(run=True)
    def sync_hydrolysate_target():
        HEV.target_solids_wt_frac = EV.target_solids_wt_frac

    FM = bst.Mixer("FM", ins=(EV - 0, HEV - 0), outs=("fermentation_feed",))
    units = (MI, *press_units, PR, *juice_units, PFS, PC, DIL, EV, *pretreatment_units, M_SA, H_SA, SA, SLS, HPFS, HPC, HDIL, HEV, FM)
    return dict(units=units, fermentation_feed=FM - 0, boiler_solids=(SLS - 0,))


def create_3hp_system(
    feedstock: str = "wiatrowski2026",
    ca3hp2_price: float | None = None,
    scenario: str = "baseline",
):
    """
    Build the full feedstock-to-product system: raw Sargassum -> Mill ->
    pulp/juice upstream -> HPFermentation -> cell-mass removal ->
    broth concentration -> CaHPCrystallizer -> crystal separation -> drying
    -> dried Ca(3HP)2, split (accounting only) into a pure-Ca3HP2 stream and
    an impurities stream. Upstream (see the module docstring): mill, plain
    press; juice concentrated; pulp enzymatically hydrolyzed in a 96 h tank,
    filtered and concentrated; the two concentrates blended into one
    fermentation feed. Scenario knobs: `M_SA.enzyme_loading_mg_per_g_dwb`,
    `SA.hydrolysis_conversion`, `EV.target_solids_wt_frac`.

    Parameters
    ----------
    feedstock : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    ca3hp2_price : float, optional
        Price (USD/kg of pure Ca3HP2) to set on the `ca3hp2_product`
        stream, i.e. on the Ca3HP2 content of the dried product only.
        Defaults to data/tea.yaml `price.3hp.baseline` when not given.
    scenario : str
        'baseline' (milled biomass -> press -> mild thermochemical pretreatment and
        neutralization of the pulp -> enzymatic hydrolysis) or
        'pretreatment_press' (1 wt% H2SO4 and 1 M NaCl dosed into the pre-press biomass, same
        hydrolysis as the baseline pretreatment assumed in the press; ammonia
        neutralization in the hydrolysis step and for the juice).
    Returns
    -------
    sys : bst.System
        Key streams and units are accessible via `sys.flowsheet.stream`
        and `sys.flowsheet.unit` (e.g. 'sargassum_feed', 'pressed_cake',
        'ca3hp2_product', 'MI', 'PR', 'SA', 'SLS', 'F401', 'C401').
    """
    # Always rebuilds the include_hp3 chemical superset -- see module
    # docstring for why this can't reuse an already-set thermo the way
    # every other sabre system's `try: get_chemicals() except: create_...`
    # guard does.
    if scenario not in SCENARIOS:
        raise ValueError(f"scenario must be one of {SCENARIOS}, got {scenario!r}")
    create_chemicals(
        set_thermo=True, include_hp3=True,
        include_dilute_acid=True,
    )

    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
        dry_composition=params.get("dry_composition"),
    )
    feed.price = _TEA_PRICE["sargassum"]["baseline"]

    # -------------------------------------------------
    # Preprocessing: Mill -> press and pulp/juice upstream
    # -------------------------------------------------
    MI = Mill("MI", ins=feed, outs=("milled_biomass", "milling_losses"))
    # milling_losses is combustible biomass, same as pressed_cake -- routed
    # to BT below instead of priced as disposal.

    # What the upstream section hands to the fermentation-to-product section:
    #   upstream_units     -- units in path order, MI included
    #   fermentation_feed  -- stream feeding HPFermentation
    #   boiler_solids      -- solid streams (besides milling_losses and
    #                         cell_mass) that go to BT
    upstream = _build_juice_pulp_upstream(MI, params.get("press"), scenario)
    upstream_units = upstream["units"]
    fermentation_feed = upstream["fermentation_feed"]
    boiler_solids = upstream["boiler_solids"]

    # -------------------------------------------------
    # Fermentation
    # -------------------------------------------------
    lime = bst.Stream("lime")
    lime.price = _TEA_PRICE["lime"]["baseline"]

    F401 = HPFermentation("F401", ins=(fermentation_feed, lime), outs=("fermentation_vent", "fermentation_broth"))

    # -------------------------------------------------
    # Cell-mass removal
    # -------------------------------------------------
    S401 = bst.units.SolidsCentrifuge(
        "S401",
        ins=F401 - 1,
        outs=("cell_mass", "clarified_broth"),
        split=dict(CellMass=_CELL_MASS_CENTRIFUGE["cellmass_split_to_cake"], Ca3HP2=0.0),
        moisture_content=_CELL_MASS_CENTRIFUGE["moisture_content"],
    )
    # cell_mass is combustible biomass, same as pressed_cake/milling_losses
    # -- routed to BT below instead of priced as disposal, per the user's
    # explicit instruction to route cell mass to the boiler.

    # -------------------------------------------------
    # Broth concentration ahead of crystallization
    # -------------------------------------------------
    # Mother-liquor recycle enters here, ahead of F402, so the evaporator's
    # concentration target (adjust_evaporation below) sees the recycled
    # solutes too. `ml_recycle` is an empty placeholder until SP_ML (below)
    # takes it over as its first outlet.
    ml_recycle = bst.Stream("mother_liquor_recycle")
    M_ML = bst.Mixer("M_ML", ins=(S401 - 1, ml_recycle), outs=("broth_to_evaporator",))

    F402 = bst.MultiEffectEvaporator(
        "F402",
        ins=M_ML - 0,
        outs=("broth_concentrate", "F402_evaporator_vapor"),
        P=tuple(_BROTH_EVAPORATOR["P_Pa"]),
        V=_BROTH_EVAPORATOR["V"],
        V_definition="First-effect",
        thermo=(F401.outs[1].thermo.ideal()),
        flash=False,
    )
    F402.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]
    F402.target_product_concentration_g_per_L = _BROTH_EVAPORATOR["target_product_concentration_g_per_L"]

    P_original = tuple(F402.P)
    Pstart, Plast, N = P_original[0], P_original[-1], len(P_original)

    def _concentration_objective(V):
        F402.V = V
        F402.run()
        effluent = F402.outs[0]
        total = effluent.F_mass
        if total <= 0:
            return 0.0
        water = effluent.imass["Water"]
        conc_g_per_L = 1000.0 * (1.0 - water / total)
        return F402.target_product_concentration_g_per_L - conc_g_per_L

    @F402.add_specification(run=False)
    def adjust_evaporation():
        V_last = F402.V
        x0, x1 = 0.0, 0.5
        F402.P = P_original
        F402._reload_components = True

        y0 = _concentration_objective(x0)
        if y0 <= 0.0:
            F402.V = x0
            return
        F402._load_components()
        for i in range(1, N):
            if _concentration_objective(1e-6) < 0.0:
                F402.P = tuple(np.linspace(Pstart, Plast, N - 1))
                F402._reload_components = True
            else:
                break
        y1 = _concentration_objective(x1)
        F402.V = flx.IQ_interpolation(
            _concentration_objective,
            x0, x1, y0, y1,
            x=V_last, ytol=1e-5, xtol=1e-6,
        )

    F402_P = bst.Pump("F402_P", ins=F402 - 0, outs=("pumped_broth_concentrate",), P=101325.0)

    # -------------------------------------------------
    # Crystallization + crystal separation
    # -------------------------------------------------
    C401 = CaHPCrystallizer(
        "C401",
        ins=F402_P - 0,
        outs=("crystal_slurry",),
        target_recovery=_CRYSTALLIZATION["target_recovery"],
        T=_CRYSTALLIZATION["T_K"],
        tau=_CRYSTALLIZATION["tau"],
        V=_CRYSTALLIZATION["V_m3"],
    )

    # The recovery fraction here -- not C401's own phase tags -- is what
    # actually enforces target_recovery downstream. The cake also entrains
    # mother liquor at `moisture_content` water, so the dried product carries
    # the liquor's dissolved impurities (see data/3hp.yaml crystal_separator
    # and CrystalCentrifuge's own docstring).
    S402 = CrystalCentrifuge(
        "S402",
        ins=C401 - 0,
        outs=("crystal_cake", "mother_liquor"),
        product_recovery=_CRYSTALLIZATION["target_recovery"],
        moisture_content=_CRYSTAL_SEPARATOR["moisture_content"],
    )

    # Mother-liquor split: `recycle_frac` back to M_ML (crystallization
    # only, NOT to the fermenter -- see data/3hp.yaml mother_liquor_recycle),
    # remainder purged as liquid waste. Nothing in this loop consumes the
    # unconverted sugars, so the purge carries all of them out at steady
    # state whatever recycle_frac is.
    SP_ML = bst.Splitter(
        "SP_ML", ins=S402 - 1,
        outs=(ml_recycle, "mother_liquor_purge"),
        split=_ML_RECYCLE["recycle_frac"],
    )
    SP_ML.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

    ml_loop = bst.System(
        "ml_loop", path=[M_ML, F402, F402_P, C401, S402, SP_ML], recycle=ml_recycle,
    )

    # -------------------------------------------------
    # Drying
    # -------------------------------------------------
    F403 = bst.units.DrumDryer(
        "F403",
        ins=S402 - 0,
        outs=("dried_product", "dryer_vent"),
        split={"Water": 1},
        moisture_content=_DRYER["target_moisture_content"],
    )

    # Accounting split only (MockSplitter: no cost, no physical separation).
    # The dried product ships with the impurities entrained from the mother
    # liquor (see CrystalCentrifuge), but only its Ca3HP2 content is priced,
    # so `ca3hp2_product` (the pure-Ca3HP2 stream, and the stream MSP is
    # solved on) is credited at the product price and `product_impurities`
    # earns nothing.
    PS = bst.MockSplitter(
        "PS", ins=F403 - 0, outs=("ca3hp2_product", "product_impurities"),
    )
    PS.outs[0].price = (
        _TEA_PRICE["3hp"]["baseline"] if ca3hp2_price is None
        else float(ca3hp2_price)
    )

    @PS.add_specification(run=True)
    def split_pure_product():
        dried = PS.ins[0]
        pure, impurities = PS.outs
        pure.empty()
        impurities.empty()
        pure.phase = impurities.phase = dried.phase
        pure.imol["Ca3HP2"] = dried.imol["Ca3HP2"]
        impurities.mol[:] = dried.mol - pure.mol
        pure.T = impurities.T = dried.T
        pure.P = impurities.P = dried.P

    path = [
        *upstream_units,
        F401, S401, ml_loop, F403, PS,
    ]
    HXN = bst.HeatExchangerNetwork(
        "HXN",
        units=(
            # H_SA heats a solids-laden slurry (Glucan, Alginate, ...) whose
            # components have no vapor-phase thermo, which breaks the pinch
            # analysis; it and the acid-pretreatment units (same slurry, plus
            # the flash) are left out of heat integration (their steam is
            # charged in full).
            *(u for u in upstream_units if u.ID not in _NOT_HEAT_INTEGRATED),
            F401, S401, M_ML, F402, F402_P, C401, S402, SP_ML, F403,
        ),
    )
    path.append(HXN)

    # -------------------------------------------------
    # Boiler/turbogenerator: combusts the three combustible solid streams
    # (pressed_cake, milling_losses, cell_mass) for steam/electricity
    # credit, per the design spec sec. 5 and the user's explicit
    # instruction to route both press cake and cell mass to the boiler --
    # not priced as disposal (see the comments at PR/MI/S401 above).
    # Standalone bst.BoilerTurbogenerator, plus a CoolingTower and a
    # ChilledWaterPackage (below) so the cooling and chilled water the units
    # draw have their capital cost, power and makeup water in the TEA and the
    # LCA, and a ProcessWaterCenter for their makeup water. Not
    # biorefineries.cellulosic's create_facilities wrapper, which also builds
    # other facilities (CIPpackage, AirDistributionPackage, FireWaterTank)
    # sabre has no use for.
    # bst defaults (boiler_efficiency=0.80, turbogenerator_efficiency=0.85,
    # fuel_price=0.218) are used as-is; these are standard literature
    # baseline values already, not sabre-specific assumptions, so they're
    # not threaded through data/3hp.yaml. The one exception is
    # ash_disposal_price, which is set to sabre's own
    # data/tea.yaml price.disposal_solid baseline (not bst's default
    # -0.0318) so BT's ash is charged on the same basis as every other
    # sabre solid waste stream.
    M_BT = bst.Mixer(
        "M_BT",
        ins=(*boiler_solids, MI - 1, S401 - 0),
        outs=("solids_to_boiler",),
    )
    path.append(M_BT)

    BT = bst.BoilerTurbogenerator(
        "BT", ins=(M_BT - 0,),
        ash_disposal_price=_TEA_PRICE["disposal_solid"]["baseline"],
    )
    # BT.outs: [0] emissions (gas, vents to atmosphere -- unpriced, same
    # convention as every other sabre vent/off-gas stream); [1] blowdown
    # water (genuine liquid waste, priced as disposal_wastewater, same as
    # every other sabre liquid waste stream); [2] ash_disposal (its cost
    # is already handled internally via BT's own `ash_disposal_price` ->
    # define_utility('Ash disposal', ...) mechanism, so it is NOT also
    # priced here -- doing so would double-count the same cost).
    BT.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

    CT = bst.CoolingTower("CT")
    # Cooling-tower blowdown is liquid waste like the other purged streams.
    CT.outs[0].price = _TEA_PRICE["disposal_wastewater"]["baseline"]
    CWP = bst.ChilledWaterPackage("CWP")
    # Process water center, as in biosteam's create_facilities and
    # biorefineries.cellulosic: balances the boiler and cooling-tower makeup
    # water (reverse-osmosis grade) and any fresh process water in the system
    # (both defaulted at run time), charging their prices and capital.
    PWC = bst.ProcessWaterCenter("PWC")
    # Facility order matters: they run in the order given, and the cooling
    # tower and PWC size themselves from the duties and makeup flows of the
    # facilities before them. The cooling tower must come after both the
    # boiler (turbogenerator condenser, ~6 MW) and the chilled water package
    # (~5 MW), which put their duty on cooling water, so this follows
    # biorefineries.cellulosic's order (BT, CWP, CT, PWC). Neither
    # bst.create_facilities (CT before CWP) nor bst.Facility.ordered_facilities
    # (CWP, CT, BT, PWC) does, and either sizes the cooling tower 25-40% too
    # small here; bst.System(..., facilities=) uses the order as given.
    facilities = (BT, CWP, CT, PWC)

    sys = bst.System(
        "hp3_sys",
        path=path, facilities=facilities,
    )
    create_tea(sys)

    return sys


def price_3hp_system(scenario: str = "baseline") -> dict:
    """
    Build, simulate and price the 3-HP system (`create_3hp_system()`): solves
    the minimum selling price of the pure Ca(3HP)2 stream. Returns a dict with
    the MSP (`msp_usd_per_kg`), annual product (`annual_product_kg`) and the
    simulated system (`sys`).
    """
    from biorefineries.sabre._tea import solve_product_msp

    bst.main_flowsheet.clear()
    sys = create_3hp_system(scenario=scenario)
    sys.simulate()
    msp = solve_product_msp(tea=sys.TEA, product_stream=sys.flowsheet.stream.ca3hp2_product)
    return {
        "label": "Calcium 3-hydroxypropionate",
        "product_desc": "dried crystalline Ca(3HP)2",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "sys": sys,
    }


if __name__ == '__main__':
    results = price_3hp_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Sargassum-to-dried-Ca(3HP)2 system builder for the SaBRe flowsheets.

Pathway:
    raw Sargassum
      -> Mill (milled_biomass, milling_losses)
      -> EnzymaticPress (pressed_cake, pressate)
      -> PFS (splits pressate into only the fraction that needs
         concentrating) -> PressateConcentrator -> DIL (blends concentrate
         + bypass pressate + fresh water) -> fermentation_feed_evaporator
         (targets ~20 wt% dissolved solids, analogous to
         biorefineries.HP's EnzymeHydrolysateMixer.solid_loading)
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
include_hp3=True chemical superset (Glucose, AlginateMonomer, Enzyme, HP,
CalciumDihydroxide, Ca3HP2), which every other sabre system deliberately
excludes to avoid a small (~0.08%) NumPy-pairwise-summation-driven drift on
their own golden-value regression tests (see create_chemicals's own
docstring). Because thermosteam's chemical set is process-global, this
system must not be built in the same Python process as any of the other
five sabre systems -- create_3hp_system() unconditionally (re)builds the
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
    Mill, Press, EnzymaticPress, PressateConcentrator, BiostimulantEvaporator,
    HPFermentation, CaHPCrystallizer, CrystalCentrifuge, DiluteAcidReactor,
)
from biorefineries.sabre._tea import create_tea

__all__ = ('create_3hp_system', 'price_3hp_system', 'CONFIGURATIONS')

# Upstream configurations (see create_3hp_system's docstring).
CONFIGURATIONS = ("enzymatic_press", "biostimulant_cake", "juice_cake")

_TEA_PRICE = load_assumptions("tea.yaml")["price"]
_3HP_YAML = load_assumptions("3hp.yaml")
_ACID_PRETREATMENT = _3HP_YAML["cake_dilute_acid_pretreatment"]
_SACCHARIFICATION = _3HP_YAML["cake_saccharification"]
_HYDROLYSATE_SEPARATOR = _3HP_YAML["hydrolysate_separator"]
_ENZYMATIC_PRESS = load_assumptions("preprocessing.yaml")["enzymatic_press"]

# Unit IDs of the "biostimulant_cake" configuration's biostimulant subsystem
# (its Press plus the prefixed PFS/PC/DIL/EV chain), excluded from the 3-HP-
# only TEA scope when price_3hp_system(credit_tipping_fee=True) -- same role
# as systems._biostimulant_system.BIOSTIMULANT_UNIT_IDS for the AD pathways.
_BIOSTIMULANT_SCOPE_UNIT_IDS = frozenset({"PR", "BPFS", "BPC", "BDIL", "BEV"})
_PRESS_CAKE_SOLIDS_WT_FRAC_OVERRIDE = _3HP_YAML["press_cake_solids_wt_frac_override"]
_PRESSATE_CONCENTRATOR = _3HP_YAML["pressate_concentrator"]
_FEED_EVAPORATOR = _3HP_YAML["fermentation_feed_evaporator"]
_BROTH_EVAPORATOR = _3HP_YAML["broth_evaporator"]
_CELL_MASS_CENTRIFUGE = _3HP_YAML["cell_mass_centrifuge"]
_CRYSTALLIZATION = _3HP_YAML["crystallization"]
_CRYSTAL_SEPARATOR = _3HP_YAML["crystal_separator"]
_ML_RECYCLE = _3HP_YAML["mother_liquor_recycle"]
_DRYER = _3HP_YAML["dryer"]

FRESH_WATER_PRICE_USD_PER_KG = bst.stream_utility_prices['Process water']

# PressateConcentrator/BiostimulantEvaporator settings of the 3-HP sugar
# concentration chain (data/3hp.yaml), same values the enzymatic_press
# configuration passes inline.
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
    systems._biostimulant_system.create_biostimulant_system builds (and the
    `enzymatic_press` configuration builds inline above), factored out here
    so the cake configurations can build it twice with different settings
    ('B'-prefixed biostimulant chain and the 3-HP sugar chain). PFS's
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


def _make_cake_hydrolysis_train(cake):
    """
    Dilute-acid pretreatment -> neutralization -> saccharification ->
    solid-liquid separation of the press cake, mirroring
    biorefineries.cellulosic's corn stover train
    (systems/pretreatment/dilute_acid.py, then its saccharification system)
    with sabre-local units. See data/3hp.yaml (`cake_dilute_acid_pretreatment`,
    `cake_saccharification`, `hydrolysate_separator`).

    Acid and ammonia are dosed by specifications on their storage tanks,
    so the fresh-feed streams (the ones priced) carry the actual demand.
    Steam for the steam mixer is booked as a heat utility by
    bst.SteamMixer itself, so there is no priced steam stream.

    Returns (units, residue, hydrolysate): the units in path order, the
    filter's residual-solids outlet, and its hydrolysate outlet.
    """
    A = _ACID_PRETREATMENT
    S = _SACCHARIFICATION
    H = _HYDROLYSATE_SEPARATOR
    price_wastewater = _TEA_PRICE["disposal_wastewater"]["baseline"]

    def dry_cake_kgph():
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
    AM = bst.Mixer("AM", ins=(warm_water_1, ACID_ST - 0), outs=("acid_liquor",))

    @ACID_ST.add_specification(run=True)
    def dose_acid():
        solution_kgph = A["acid_solution_loading_kg_per_kg_dry"] * dry_cake_kgph()
        total_ref = A["acid_solution_H2SO4_kg"] + A["acid_solution_water_kg"]
        sulfuric_acid.imass["H2SO4"] = solution_kgph * A["acid_solution_H2SO4_kg"] / total_ref
        sulfuric_acid.imass["Water"] = solution_kgph * A["acid_solution_water_kg"] / total_ref

    @AM.add_specification(run=True)
    def dose_warm_water():
        warm_water_1.imass["Water"] = A["warm_water_kg_per_kg_dry"] * dry_cake_kgph()

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

    # --- Saccharification -------------------------------------------------
    enzyme = bst.Stream("saccharification_enzyme", units="kg/hr")
    enzyme.price = _TEA_PRICE["enzyme"]["baseline"]
    saccharification_water = bst.Stream("saccharification_water", Water=0.0, units="kg/hr")
    saccharification_water.price = FRESH_WATER_PRICE_USD_PER_KG
    M_SA = bst.Mixer(
        "M_SA", ins=(NT - 0, enzyme, saccharification_water),
        outs=("slurry_with_enzyme",),
    )
    # Dissolved ammonium sulfate is locked to the solid phase (like corn
    # stover's) but is not an insoluble solid for loading/filtration purposes.
    solids_IDs = get_solids_group_IDs()
    insoluble_IDs = [i for i in solids_IDs if i != "AmmoniumSulfate"]
    enzyme_dose_frac = (
        _ENZYMATIC_PRESS["enzyme_loading_mg_per_g_substrate"] / 1000.0
        * (1.0 + _ENZYMATIC_PRESS["enzyme_excess_frac"])
    )

    @M_SA.add_specification(run=True)
    def dose_enzyme_and_water():
        slurry = M_SA.ins[0]
        enzyme.empty()
        enzyme.phase = "l"
        saccharification_water.empty()
        enzyme.imass["Enzyme"] = enzyme_dose_frac * (
            slurry.imass["Glucan"] + slurry.imass["Alginate"]
        )
        insoluble_kgph = float(slurry.imass[insoluble_IDs].sum())
        total_kgph = slurry.F_mass + enzyme.F_mass
        saccharification_water.imass["Water"] = max(
            insoluble_kgph / S["solids_loading"] - total_kgph, 0.0
        )

    H_SA = bst.HXutility(
        "H_SA", ins=M_SA - 0, outs=("slurry_to_saccharification",), T=S["T_K"],
    )
    SA = bst.MixTank("SA", ins=H_SA - 0, outs=("saccharified_slurry",), tau=S["tau_hr"])
    hydrolysis_rxns = tmo.ParallelReaction([
        tmo.Reaction("Glucan + Water -> Glucose", reactant="Glucan",
                     X=S["hydrolysis_conversion"]),
        tmo.Reaction("Alginate + Water -> AlginateMonomer", reactant="Alginate",
                     X=S["hydrolysis_conversion"]),
    ])

    @SA.add_specification
    def saccharify():
        SA._run()
        hydrolysis_rxns(SA.outs[0])

    # --- Solid-liquid separation -------------------------------------------
    split = {
        cid: (H["solids_split_to_cake"] if cid in insoluble_IDs else H["solutes_split_to_cake"])
        for cid in SA.chemicals.IDs if cid != "Water"
    }
    SLS = bst.PressureFilter(
        "SLS", ins=SA - 0, outs=("hydrolysis_residue", "hydrolysate"),
        split=split, moisture_content=H["moisture_content"],
    )

    units = [ACID_ST, AM, M_PT, R_PT, F_PT, H_PT, NH3_ST, NT, M_SA, H_SA, SA, SLS]
    return units, SLS - 0, SLS - 1


def _build_cake_saccharification_upstream(configuration, MI, biostimulant_price):
    """
    Upstream section of the 'biostimulant_cake' and 'juice_cake'
    configurations: plain Press on the milled biomass, the cake through
    `_make_cake_hydrolysis_train`, and the 3-HP sugar concentration chain.
    In 'biostimulant_cake' the pressate goes to a separate biostimulant
    chain (unit IDs prefixed 'B') instead; in 'juice_cake' it is mixed with
    the hydrolysate (M_JH), like biorefineries.cane's hydrolysate-and-juice
    mixer.

    Returns a dict with 'units' (path order, MI included),
    'fermentation_feed' and 'boiler_solids'.
    """
    PR = Press(
        "PR", ins=MI - 0, outs=("pressed_cake", "pressate"),
        cake_solids_wt_frac=_PRESS_CAKE_SOLIDS_WT_FRAC_OVERRIDE,
    )
    units = [MI, PR]

    if configuration == "biostimulant_cake":
        BPFS, BPC, BDIL, BEV = _make_pressate_concentration_chain(
            PR - 1, unit_prefix="B", stream_prefix="biostimulant_",
            fresh_water_ID="biostimulant_fresh_water",
            product_ID="biostimulant_product",
            condensate_ID="biostimulant_condensed_vapor",
            pc_kwargs={}, ev_kwargs={},  # biostimulant.yaml defaults
        )
        BEV.outs[0].price = (
            _TEA_PRICE["biostimulant"]["baseline"] if biostimulant_price is None
            else float(biostimulant_price)
        )
        units += [BPFS, BPC, BDIL, BEV]

    train_units, residue, hydrolysate = _make_cake_hydrolysis_train(PR.outs[0])
    units += train_units

    if configuration == "juice_cake":
        M_JH = bst.Mixer("M_JH", ins=(PR - 1, hydrolysate), outs=("juice_and_hydrolysate",))
        units.append(M_JH)
        sugar_feed = M_JH - 0
    else:
        sugar_feed = hydrolysate

    PFS, PC, DIL, EV = _make_pressate_concentration_chain(
        sugar_feed, unit_prefix="", stream_prefix="",
        fresh_water_ID="fermentation_fresh_water",
        product_ID="fermentation_feed", condensate_ID="EV_condensed_vapor",
        pc_kwargs=_SUGAR_PC_KWARGS, ev_kwargs=_SUGAR_EV_KWARGS,
    )
    units += [PFS, PC, DIL, EV]

    return dict(units=tuple(units), fermentation_feed=EV - 0, boiler_solids=(residue,))


def create_3hp_system(
    feedstock: str = "pelagic",
    ca3hp2_price: float | None = None,
    configuration: str = "enzymatic_press",
    biostimulant_price: float | None = None,
):
    """
    Build the full feedstock-to-product system: raw Sargassum -> Mill ->
    [upstream, per `configuration`] -> HPFermentation -> cell-mass removal ->
    broth concentration -> CaHPCrystallizer -> crystal separation -> drying
    -> dried Ca(3HP)2, split (accounting only) into a pure-Ca3HP2 stream and
    an impurities stream. Everything from HPFermentation on is identical
    across configurations.

    Parameters
    ----------
    feedstock : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    ca3hp2_price : float, optional
        Price (USD/kg of pure Ca3HP2) to set on the `ca3hp2_product`
        stream, i.e. on the Ca3HP2 content of the dried product only.
        Defaults to data/tea.yaml `price.3hp.baseline` when not given.
    configuration : str
        One of `CONFIGURATIONS`:

        - 'enzymatic_press' (default): EnzymaticPress hydrolyzes part of the
          glucan/alginate in the press; the pressate is concentrated and
          fermented and the pressed cake goes to the boiler.
        - 'biostimulant_cake': plain Press; the pressate becomes a
          biostimulant product and the cake is dilute-acid pretreated,
          saccharified and filtered into the fermentation feed (residual
          solids to the boiler), analogous to sugarcane bagasse.
        - 'juice_cake': plain Press; the pressate ("juice") and the cake's
          hydrolysate (same train as above) are combined into the
          fermentation feed, analogous to biorefineries.cane's combined
          1G+2G design.
    biostimulant_price : float, optional
        Price (USD/kg) set on the `biostimulant_product` stream
        ('biostimulant_cake' only). Defaults to data/tea.yaml
        `price.biostimulant.baseline`.

    Returns
    -------
    sys : bst.System
        Key streams and units are accessible via `sys.flowsheet.stream`
        and `sys.flowsheet.unit` (e.g. 'sargassum_feed', 'pressed_cake',
        'ca3hp2_product', 'MI', 'PR', 'F401', 'C401'; plus 'PT'-area units
        such as 'R_PT', 'SA', 'SLS' in the cake configurations).
    """
    if configuration not in CONFIGURATIONS:
        raise ValueError(
            f"Unknown configuration {configuration!r}. Choose from {CONFIGURATIONS}."
        )

    # Always rebuilds the include_hp3 chemical superset -- see module
    # docstring for why this can't reuse an already-set thermo the way
    # every other sabre system's `try: get_chemicals() except: create_...`
    # guard does. The cake configurations additionally need the dilute-acid
    # chemicals (H2SO4, AmmoniumSulfate); the baseline's own set is left
    # exactly as it was.
    create_chemicals(
        set_thermo=True, include_hp3=True,
        include_dilute_acid=(configuration != "enzymatic_press"),
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
    # Preprocessing: Mill -> [Press, per configuration]
    # -------------------------------------------------
    MI = Mill("MI", ins=feed, outs=("milled_biomass", "milling_losses"))
    # milling_losses is combustible biomass, same as pressed_cake -- routed
    # to BT below instead of priced as disposal.

    # What each configuration's upstream section hands to the shared
    # fermentation-to-product section below:
    #   upstream_units     -- units in path order, MI included
    #   fermentation_feed  -- stream feeding HPFermentation
    #   boiler_solids      -- solid streams (besides milling_losses and
    #                         cell_mass) that go to BT
    if configuration == "enzymatic_press":
        enzyme = bst.Stream("enzyme_cocktail")
        enzyme.price = _TEA_PRICE["enzyme"]["baseline"]

        PR = EnzymaticPress(
            "PR", ins=(MI - 0, enzyme), outs=("pressed_cake", "pressate"),
            # Overrides preprocessing.yaml's own cake_solids_wt_frac (0.15) --
            # see data/3hp.yaml press_cake_solids_wt_frac_override's own
            # sources note for why.
            cake_solids_wt_frac=_PRESS_CAKE_SOLIDS_WT_FRAC_OVERRIDE,
        )
        # pressed_cake (unconverted Glucan/Alginate, Ash, Protein, Fucoidan,
        # OtherSolids) is NOT priced as disposal -- it's routed to BT (the
        # BoilerTurbogenerator facility, wired near the end of this function)
        # for steam/electricity credit, per the design spec sec. 5.

        # -------------------------------------------------
        # Pressate concentration: PFS (splits raw pressate) -> PC -> DIL ->
        # fermentation_feed_evaporator (~20 wt% dissolved solids target).
        # Mirrors systems._biostimulant_system.create_biostimulant_system's own
        # PFS/PC/DIL/EV chain exactly (same adjust_pressate_split logic),
        # retargeted to this pathway's own solute set/yaml section.
        # -------------------------------------------------
        fresh_water = bst.Stream("fermentation_fresh_water", Water=0.0, units="kg/hr")
        fresh_water.price = FRESH_WATER_PRICE_USD_PER_KG

        PFS = bst.MockSplitter(
            "PFS",
            ins=PR - 1,
            outs=("to_concentrator", "bypass_pressate"),
        )

        PC = PressateConcentrator(
            "PC", ins=PFS - 0,
            outs=("membrane_concentrate", "permeate"),
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
        PC.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

        DIL = bst.Mixer(
            "DIL",
            ins=(PC - 0, PFS - 1, fresh_water),
            outs=("diluted_pressate",),
        )

        EV = BiostimulantEvaporator(
            "EV",
            ins=DIL - 0,
            outs=("fermentation_feed", "EV_condensed_vapor"),
            target_solids_wt_frac=_FEED_EVAPORATOR["target_solids_wt_frac"],
            boiling_temperature_K=_FEED_EVAPORATOR["boiling_temperature_K"],
            condensing_temperature_K=_FEED_EVAPORATOR["condensing_temperature_K"],
            nonwater_recovery_to_product=_FEED_EVAPORATOR["nonwater_recovery_to_product"],
            tau=_FEED_EVAPORATOR["tau"],
            vacuum_system_preference=_FEED_EVAPORATOR["vacuum_system_preference"],
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

        upstream_units = (MI, PR, PFS, PC, DIL, EV)
        fermentation_feed = EV - 0
        boiler_solids = (PR - 0,)

    else:
        upstream = _build_cake_saccharification_upstream(configuration, MI, biostimulant_price)
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
            *upstream_units,
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
    # Minimal standalone bst.BoilerTurbogenerator usage (no CoolingTower/
    # ProcessWaterCenter/etc.), matching the class's own docstring example
    # rather than biorefineries.cellulosic's create_facilities wrapper
    # (which also builds four *other* facilities sabre has no use for --
    # ChilledWaterPackage, CIPpackage, AirDistributionPackage,
    # FireWaterTank -- none of which exist anywhere else in sabre either).
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

    sys = bst.System(
        "hp3_sys" if configuration == "enzymatic_press" else f"hp3_{configuration}_sys",
        path=path, facilities=(BT,),
    )
    create_tea(sys)

    return sys


def price_3hp_system(
    configuration: str = "enzymatic_press",
    credit_tipping_fee: bool = False,
) -> dict:
    """
    Parameters
    ----------
    configuration : str
        One of `CONFIGURATIONS` (see create_3hp_system).
    credit_tipping_fee : bool
        Only for configuration='biostimulant_cake', mirroring
        systems._biomethane_system.price_biomethane_system's switch of the
        same name.

        If False (default): biostimulant is priced at its flat data/tea.yaml
        `price.biostimulant.baseline`, and Ca(3HP)2's price is solved on a
        single TEA covering every unit, biostimulant's included.

        If True: biostimulant is instead priced at its own standalone MSP
        (price_biostimulant_system()), and Ca(3HP)2's price is solved
        against a 3-HP-only TEA scope that excludes the biostimulant
        subsystem's capital (its Press and 'B'-prefixed chain) -- which
        biostimulant's own price already recovers -- plus a tipping-fee
        credit equal to data/tea.yaml `price.disposal_solid.baseline` on the
        pressed_cake mass taken in, since biostimulant would otherwise have
        paid that same fee to dispose of it.
    """
    from biorefineries.sabre._tea import (
        solve_product_msp, disposal_avoided_credit_usd_per_yr, apply_revenue_credit,
    )
    if credit_tipping_fee and configuration != "biostimulant_cake":
        raise ValueError("credit_tipping_fee only applies to configuration='biostimulant_cake'.")

    biostimulant_price = None
    if configuration == "biostimulant_cake":
        if credit_tipping_fee:
            from biorefineries.sabre.systems._biostimulant_system import price_biostimulant_system
            # Base chemical set first, so the standalone biostimulant price
            # doesn't depend on whichever set an earlier call left behind.
            create_chemicals(set_thermo=True)
            biostimulant_price = price_biostimulant_system()["msp_usd_per_kg"]
        else:
            biostimulant_price = _TEA_PRICE["biostimulant"]["baseline"]

    bst.main_flowsheet.clear()

    sys = create_3hp_system(configuration=configuration, biostimulant_price=biostimulant_price)
    sys.simulate()

    product = sys.flowsheet.stream.ca3hp2_product

    tipping_fee_usd_per_yr = 0.0
    if credit_tipping_fee:
        specific_units = [u for u in sys.units if u.ID not in _BIOSTIMULANT_SCOPE_UNIT_IDS]
        specific_sys = bst.System.from_units("hp3_specific_sys", units=specific_units)
        specific_tea = create_tea(specific_sys)
        msp = solve_product_msp(tea=specific_tea, product_stream=product)
        msp_before_tipping_credit = msp["usd_per_kg"]
        tipping_fee_usd_per_yr = disposal_avoided_credit_usd_per_yr(
            sys.flowsheet.stream.pressed_cake, _TEA_PRICE["disposal_solid"]["baseline"], specific_tea,
        )
        msp = apply_revenue_credit(msp, tipping_fee_usd_per_yr)
    else:
        msp = solve_product_msp(tea=sys.TEA, product_stream=product)

    result = {
        "label": "Calcium 3-hydroxypropionate",
        "product_desc": "dried crystalline Ca(3HP)2",
        "configuration": configuration,
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "sys": sys,
    }
    if configuration == "biostimulant_cake":
        result.update(
            biostimulant_price_usd_per_kg=biostimulant_price,
            credit_tipping_fee=credit_tipping_fee,
            tipping_fee_usd_per_yr=tipping_fee_usd_per_yr,
        )
        if credit_tipping_fee:
            # The 3-HP-only-scope MSP before the tipping-fee credit, and that
            # scope's TCI -- the credit is large enough to dominate the result.
            result.update(
                msp_before_tipping_credit_usd_per_kg=msp_before_tipping_credit,
                specific_scope_TCI=specific_tea.TCI,
            )
    return result


if __name__ == '__main__':
    results = price_3hp_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")

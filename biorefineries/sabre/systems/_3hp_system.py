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
      -> CaHPCrystallizer -> S402 (SolidsCentrifuge: separates crystal
         cake from mother liquor)
      -> SP_ML (splits mother liquor: `mother_liquor_recycle` back to M_ML,
         data/3hp.yaml mother_liquor_recycle.recycle_frac; remainder
         `mother_liquor_purge` -> disposal_wastewater)
      -> DrumDryer -> dried Ca(3HP)2 product
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

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
)
from biorefineries.sabre.units import (
    Mill, EnzymaticPress, PressateConcentrator, BiostimulantEvaporator,
    HPFermentation, CaHPCrystallizer,
)
from biorefineries.sabre._tea import create_tea

__all__ = ('create_3hp_system', 'price_3hp_system')

_TEA_PRICE = load_assumptions("tea.yaml")["price"]
_3HP_YAML = load_assumptions("3hp.yaml")
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


def create_3hp_system(
    feedstock: str = "pelagic",
    ca3hp2_price: float | None = None,
):
    """
    Build the full feedstock-to-product system: raw Sargassum -> Mill ->
    EnzymaticPress -> pressate concentration -> HPFermentation ->
    cell-mass removal -> broth concentration -> CaHPCrystallizer ->
    crystal separation -> drying -> dried Ca(3HP)2.

    Parameters
    ----------
    feedstock : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    ca3hp2_price : float, optional
        Price (USD/kg) to set on the dried Ca(3HP)2 product stream.
        Defaults to data/tea.yaml `price.3hp.baseline` when not given.

    Returns
    -------
    sys : bst.System
        Key streams and units are accessible via `sys.flowsheet.stream`
        and `sys.flowsheet.unit` (e.g. 'sargassum_feed', 'pressed_cake',
        'ca3hp2_product', 'MI', 'PR', 'F401', 'C401').
    """
    # Always rebuilds the include_hp3 chemical superset -- see module
    # docstring for why this can't reuse an already-set thermo the way
    # every other sabre system's `try: get_chemicals() except: create_...`
    # guard does.
    create_chemicals(set_thermo=True, include_hp3=True)

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
    # Preprocessing: Mill -> EnzymaticPress
    # -------------------------------------------------
    MI = Mill("MI", ins=feed, outs=("milled_biomass", "milling_losses"))
    # milling_losses is combustible biomass, same as pressed_cake -- routed
    # to BT below instead of priced as disposal.

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

    # -------------------------------------------------
    # Fermentation
    # -------------------------------------------------
    lime = bst.Stream("lime")
    lime.price = _TEA_PRICE["lime"]["baseline"]

    F401 = HPFermentation("F401", ins=(EV - 0, lime), outs=("fermentation_vent", "fermentation_broth"))

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

    # See data/3hp.yaml crystal_separator's own comment: the split fraction
    # here -- not C401's own phase tags -- is what actually enforces
    # target_recovery downstream (confirmed by direct testing that
    # bst.units.SolidsCentrifuge collapses a multi-phase feed's per-chemical
    # mol across phases before splitting).
    S402 = bst.units.SolidsCentrifuge(
        "S402",
        ins=C401 - 0,
        outs=("crystal_cake", "mother_liquor"),
        split=dict(Ca3HP2=_CRYSTALLIZATION["target_recovery"]),
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
        outs=("ca3hp2_product", "dryer_vent"),
        split={"Water": 1},
        moisture_content=_DRYER["target_moisture_content"],
    )
    F403.outs[0].price = (
        _TEA_PRICE["3hp"]["baseline"] if ca3hp2_price is None
        else float(ca3hp2_price)
    )

    path = [
        MI, PR, PFS, PC, DIL, EV,
        F401, S401, ml_loop, F403,
    ]
    HXN = bst.HeatExchangerNetwork(
        "HXN",
        units=(
            MI, PR, PFS, PC, DIL, EV,
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
        ins=(PR - 0, MI - 1, S401 - 0),
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

    sys = bst.System("hp3_sys", path=path, facilities=(BT,))
    create_tea(sys)

    return sys


def price_3hp_system() -> dict:
    from biorefineries.sabre._tea import solve_product_msp
    bst.main_flowsheet.clear()

    sys = create_3hp_system()
    sys.simulate()

    product = sys.flowsheet.stream.ca3hp2_product
    msp = solve_product_msp(tea=sys.TEA, product_stream=product)

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

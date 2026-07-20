# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Standalone biostimulant system builder for the SaBRe flowsheets.

Pathway:
    raw Sargassum
      -> Press (pressed_cake, pressate)
      -> PressateConcentrator (concentrate, permeate)
      -> PSP (MockSplitter: draws only as much of the permeate as needed
         for dilution, per DIL's specification; the rest passes through
         as residual_permeate)
      -> DIL (Mixer: concentrate + drawn permeate + drawn fresh water)
      -> BiostimulantEvaporator (evaporates any remaining excess water to
         hit target_solids_wt_frac; a no-op when DIL already diluted the
         concentrate down to the target)
      -> biostimulant_product
"""
import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
)
from biorefineries.sabre.units import Press, PressateConcentrator, BiostimulantEvaporator

__all__ = ('create_biostimulant_system',)


def create_biostimulant_system(
    feedstock_type: str = "pelagic",
    target_solids_wt_frac: float | None = 0.20,
    pressed_cake_disposal_usd_per_kg: float = 0.0,
    fresh_water_price_usd_per_kg: float = 0.0,
):
    """
    Build the standalone biostimulant system: Press -> PressateConcentrator
    -> PSP/DIL (dilute, if needed) -> BiostimulantEvaporator (concentrate,
    if needed), producing a biostimulant product at `target_solids_wt_frac`
    regardless of whether that's above or below the concentrator's natural
    output concentration.

    Parameters
    ----------
    feedstock_type : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    target_solids_wt_frac : float or None
        Target dry-solids weight fraction for the final biostimulant
        product. If lower than PressateConcentrator's natural output, the
        PSP/DIL stage dilutes with permeate water first, then fresh
        water. If higher, BiostimulantEvaporator evaporates water to
        reach it. `None` disables both stages (direct pass-through, no
        CAPEX/OPEX, no mass flow change).
    pressed_cake_disposal_usd_per_kg : float
        Price assigned to the pressed cake stream, since this pathway
        has no use for it (default 0.0; pass a negative value for a
        disposal cost).
    fresh_water_price_usd_per_kg : float
        Price assigned to the fresh water make-up stream (default 0.0).

    Returns
    -------
    sys : bst.System
    streams : dict
        'feed', 'pressed_cake', 'pressate_permeate' (PC's permeate,
        before any is drawn off), 'residual_permeate' (after PSP draws
        off whatever DIL needed for dilution, if any), 'fresh_water',
        'biostimulant_product', 'vapor'.
    units : dict
        'PR', 'PC', 'PSP', 'DIL', 'EV'.
    """
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    # Press/PC both have their own yaml-sourced defaults
    # (data/preprocessing.yaml, data/biostimulant.yaml) -- not
    # re-declared here.
    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))
    PR.outs[0].price = float(pressed_cake_disposal_usd_per_kg)

    PC = PressateConcentrator(
        "PC", ins=PR - 1,
        outs=("biostimulant_membrane_concentrate", "pressate_permeate"),
    )

    fresh_water = bst.Stream("biostimulant_fresh_water", Water=0.0, units="kg/hr")
    fresh_water.price = float(fresh_water_price_usd_per_kg)

    # Splits PC's permeate into just the amount DIL needs for dilution
    # (outs[0]) and whatever's left over (outs[1], residual_permeate).
    PSP = bst.MockSplitter(
        "PSP",
        ins=PC - 1,
        outs=("permeate_to_dilution", "residual_permeate"),
    )

    # Mixes the concentrate with whatever water PSP drew off (permeate
    # + fresh water make-up).
    DIL = bst.Mixer(
        "DIL",
        ins=(PC - 0, PSP - 0, fresh_water),
        outs=("diluted_concentrate",),
    )

    EV = BiostimulantEvaporator(
        "EV",
        ins=DIL - 0,
        outs=("biostimulant_product", "biostimulant_vapor"),
        target_solids_wt_frac=target_solids_wt_frac,
    )

    @PSP.add_specification(run=True)
    def adjust_permeate_draw():
        concentrate_in = PC.outs[0]
        permeate_in = PSP.ins[0]
        draw, residual = PSP.outs

        draw.empty()
        residual.copy_like(permeate_in)
        fresh_water.empty()

        target = EV.target_solids_wt_frac
        if target is None or not (0.0 < target < 1.0):
            return

        water_in = (
            float(concentrate_in.imass["Water"])
            if "Water" in concentrate_in.chemicals.IDs else 0.0
        )
        nonwater = concentrate_in.F_mass - water_in
        if nonwater <= 0:
            return

        actual_solids_wt_frac = nonwater / concentrate_in.F_mass
        if target >= actual_solids_wt_frac:
            # Already at or above target -- EV concentrates the rest of
            # the way by evaporating; no dilution water needed here.
            return

        # Diluting: draw the shortfall from permeate first, then fresh water
        water_to_product = nonwater * (1.0 - target) / target
        shortfall = max(water_to_product - water_in, 0.0)
        if shortfall <= 0:
            return

        permeate_water_avail = (
            float(permeate_in.imass["Water"])
            if "Water" in permeate_in.chemicals.IDs else 0.0
        )
        permeate_water_used = min(shortfall, permeate_water_avail)
        fresh_water_used = shortfall - permeate_water_used

        if permeate_water_used > 0:
            draw.imass["Water"] = permeate_water_used
            residual.imass["Water"] = permeate_water_avail - permeate_water_used
        if fresh_water_used > 0:
            fresh_water.imass["Water"] = fresh_water_used

    sys = bst.System("Biostimulant_sys", path=(PR, PC, PSP, DIL, EV))

    streams = {
        "feed": feed,
        "pressed_cake": PR.outs[0],
        "pressate_permeate": PC.outs[1],
        "residual_permeate": PSP.outs[1],
        "fresh_water": fresh_water,
        "biostimulant_product": EV.outs[0],
        "vapor": EV.outs[1],
    }
    units = {"PR": PR, "PC": PC, "PSP": PSP, "DIL": DIL, "EV": EV}

    return sys, streams, units

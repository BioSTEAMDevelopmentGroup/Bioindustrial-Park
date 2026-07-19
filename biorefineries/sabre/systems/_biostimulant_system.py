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
      -> BiostimulantEvaporator (adjusts concentrate to target_solids_wt_frac,
         drawing make-up water from the permeate first, then fresh water)
      -> biostimulant_product

The pressed cake is not used by this pathway -- it's priced as a
disposal liability (`pressed_cake_disposal_usd_per_kg`, default 0.0) via
this builder's own kwarg, not read from assumptions.yaml (this codebase
doesn't store stream prices in the yaml; every other price in these
system builders -- feed price, disposal costs -- is likewise set as a
Python-level default that callers can override, e.g.
`legacy_analyses/vfa_fermentation_tea.py`'s `*_DISPOSAL_USD_PER_KG` constants).
Fresh water make-up is priced the same way
(`fresh_water_price_usd_per_kg`, default 0.0).
"""
import biosteam as bst

from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
)
from biorefineries.sabre.units import Press, PressateConcentrator, BiostimulantEvaporator

__all__ = ('create_biostimulant_system',)


def create_biostimulant_system(
    feedstock_type: str = "pelagic",
    target_solids_wt_frac: float = 0.20,
    pressed_cake_disposal_usd_per_kg: float = 0.0,
    fresh_water_price_usd_per_kg: float = 0.0,
):
    """
    Build the standalone biostimulant system: Press -> PressateConcentrator
    -> BiostimulantEvaporator, producing a biostimulant product at
    `target_solids_wt_frac` regardless of whether that's above or below
    the concentrator's natural output concentration.

    Parameters
    ----------
    feedstock_type : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    target_solids_wt_frac : float
        Target dry-solids weight fraction for the final biostimulant
        product. If higher than what PressateConcentrator alone
        produces, BiostimulantEvaporator evaporates water to reach it.
        If lower, it dilutes with permeate water first, then fresh water.
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
        'feed', 'pressed_cake', 'pressate_permeate' (PC's permeate, before
        any is drawn off), 'residual_permeate' (after dilution draw-off,
        if any), 'fresh_water', 'biostimulant_product', 'vapor'.
    units : dict
        'PR', 'PC', 'EV'.
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

    EV = BiostimulantEvaporator(
        "EV",
        ins=(PC - 0, PC - 1, fresh_water),
        outs=("biostimulant_product", "biostimulant_vapor", "residual_permeate"),
        target_solids_wt_frac=target_solids_wt_frac,
    )

    sys = bst.System("Biostimulant_sys", path=(PR, PC, EV))

    streams = {
        "feed": feed,
        "pressed_cake": PR.outs[0],
        "pressate_permeate": PC.outs[1],
        "residual_permeate": EV.outs[2],
        "fresh_water": fresh_water,
        "biostimulant_product": EV.outs[0],
        "vapor": EV.outs[1],
    }
    units = {"PR": PR, "PC": PC, "EV": EV}

    return sys, streams, units

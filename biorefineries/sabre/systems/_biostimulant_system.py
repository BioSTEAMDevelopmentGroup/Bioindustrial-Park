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
      -> PFS (MockSplitter: splits pressate into only the fraction that
         needs concentrating and a raw, unconcentrated bypass, per DIL's
         specification -- see PFS's specification for the split-fraction
         derivation)
      -> PressateConcentrator (concentrate, permeate; only ever sees the
         to-concentrator fraction -- permeate is a pure byproduct, never
         drawn from for dilution)
      -> DIL (Mixer: concentrate + bypass pressate + drawn fresh water)
      -> BiostimulantEvaporator (evaporates any remaining excess water to
         hit target_solids_wt_frac; a no-op when DIL already diluted the
         concentrate down to the target)
      -> biostimulant_product
"""
import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
)
from biorefineries.sabre.units import Press, PressateConcentrator, BiostimulantEvaporator

__all__ = ('create_biostimulant_system',)


def create_biostimulant_system(
    feedstock_type: str = "pelagic",
):
    """
    Build the standalone biostimulant system: Press -> PFS (splits off only
    the fraction of pressate that needs concentrating) -> PressateConcentrator
    -> DIL (blends concentrate + bypass pressate + fresh water, if needed)
    -> BiostimulantEvaporator (concentrate further, if needed), producing a
    biostimulant product at its target solids content regardless of whether
    that's above or below the concentrator's natural output concentration.

    Parameters
    ----------
    feedstock_type : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).

    Returns
    -------
    sys : bst.System
    streams : dict
        'feed', 'pressed_cake', 'permeate' (PC's permeate -- a pure
        byproduct stream, never drawn from), 'fresh_water',
        'biostimulant_product', 'vapor'.
    units : dict
        'PR', 'PFS', 'PC', 'DIL', 'EV'.
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock_type)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )

    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))

    fresh_water = bst.Stream("biostimulant_fresh_water", Water=0.0, units="kg/hr")

    # Splits raw pressate into only the fraction that needs concentrating
    # (outs[0], to PC) and a raw, unconcentrated bypass (outs[1], to DIL) --
    # see adjust_pressate_split below for the split-fraction derivation.
    PFS = bst.MockSplitter(
        "PFS",
        ins=PR - 1,
        outs=("to_concentrator", "bypass_pressate"),
    )

    PC = PressateConcentrator(
        "PC", ins=PFS - 0,
        outs=("biostimulant_membrane_concentrate", "permeate"),
    )

    # Mixes PC's concentrate with PFS's bypass pressate and fresh water
    # make-up.
    DIL = bst.Mixer(
        "DIL",
        ins=(PC - 0, PFS - 1, fresh_water),
        outs=("diluted_concentrate",),
    )

    EV = BiostimulantEvaporator(
        "EV",
        ins=DIL - 0,
        outs=("biostimulant_product", "biostimulant_vapor"),
    )

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

        # Default/no-op: send everything to the concentrator. Covers no
        # target set, a degenerate (all-water or empty) feed, and a
        # disabled PC (which passes its feed through unchanged regardless
        # of how much it receives).
        if (
            target is None or not (0.0 < target < 1.0)
            or nonwater_raw <= 0 or not PC.enabled
        ):
            to_conc.copy_like(pressate)
            bypass.empty()
            return

        # PC's concentrate composition is flow-independent (its recovery
        # fractions are ratios), so what it WOULD yield if fed 100% of the
        # current pressate can be computed directly, without PC having run.
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
            # Bypassing alone can't dilute enough -- bypass everything and
            # make up the rest with fresh water.
            f = 0.0
            water_needed_total = nonwater_raw * (1.0 - target) / target
            fresh_water.imass["Water"] = max(water_needed_total - water_raw, 0.0)
        elif target >= x_conc:
            # Even full concentration isn't enough -- EV finishes by
            # evaporating; no bypass, no fresh water.
            f = 1.0
        else:
            # Blend f*concentrate + (1-f)*raw pressate to hit target exactly.
            A = conc_mass_full * (x_conc - target)
            B = F * (x_raw - target)
            f = B / (B - A)

        to_conc.copy_like(pressate)
        to_conc.mol *= f
        bypass.copy_like(pressate)
        bypass.mol *= (1.0 - f)

    path = [PR, PFS, PC, DIL, EV]

    # EV is the only HXN-eligible unit in this scope.
    # HXN would crash if there is no HXN-eligible units in the sysem.
    if EV.target_solids_wt_frac is not None:
        HXN = bst.HeatExchangerNetwork("HXN", units=tuple(path))
        path.append(HXN)

    sys = bst.System("biostimulant_sys", path=path)

    streams = {
        "feed": feed,
        "pressed_cake": PR.outs[0],
        "permeate": PC.outs[1],
        "fresh_water": fresh_water,
        "biostimulant_product": EV.outs[0],
        "vapor": EV.outs[1],
    }
    units = {"PR": PR, "PFS": PFS, "PC": PC, "DIL": DIL, "EV": EV}

    return sys, streams, units

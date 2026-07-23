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
from pathlib import Path

import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import (
    load_assumptions, get_feedstock_type_params, get_scale_feed_kgph, make_sargassum_feed,
)
from biorefineries.sabre.units import Press, PressateConcentrator, BiostimulantEvaporator
from biorefineries.sabre._tea import create_tea

__all__ = ('create_biostimulant_system', 'price_biostimulant_system', 'BIOSTIMULANT_UNIT_IDS')

_TEA_PRICE = load_assumptions("tea.yaml")["price"]

# Unit IDs built by create_biostimulant_system(). Used by the AD pathway
# systems (systems._ad_biomethane_system, systems._ad_vfa_system,
# systems._ad_fermentation_system) to carve an AD-only TEA scope out of
# their combined system -- excluding these units' capital -- so their own
# product's price recovers only that pathway's own capital, not
# biostimulant's (which biostimulant's own standalone price already
# recovers on its own; bst.System has no free attribute slot to tag this
# per-instance, so it's a fixed constant instead).
BIOSTIMULANT_UNIT_IDS = frozenset({'PR', 'PFS', 'PC', 'DIL', 'EV'})

FRESH_WATER_PRICE_USD_PER_KG = bst.stream_utility_prices['Process water']


def create_biostimulant_system(
    feedstock: str = "pelagic",
    biostimulant_price: float | None = None,
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
    feedstock : str
        Feedstock type (data/feedstock.yaml `feedstock_type`).
    biostimulant_price : float, optional
        Price (USD/kg) to set on the biostimulant_product stream. Defaults
        to data/tea.yaml `price.biostimulant.baseline` when not given -- pass
        this in (e.g. an IRR=0 cost-basis price solved from the standalone
        biostimulant system) when this subsystem is embedded as a byproduct
        credit inside a larger pathway (see systems._ad_biomethane_system,
        systems._ad_vfa_system).

    Returns
    -------
    sys : bst.System
        Key streams and units are accessible via `sys.flowsheet.stream`
        and `sys.flowsheet.unit` (e.g. 'feed', 'pressed_cake',
        'biostimulant_product', 'PR', 'PC', 'EV').
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    feedstock_assumptions = load_assumptions("feedstock.yaml")
    params = get_feedstock_type_params(feedstock_assumptions, feedstock)
    fresh_feed_kgph = get_scale_feed_kgph(feedstock_assumptions)
    moisture_frac = params["moisture_frac"]

    feed = make_sargassum_feed(
        fresh_feed_kgph=fresh_feed_kgph, moisture_frac=moisture_frac,
        ash_wt_frac_dry=params["ash_wt_frac_dry"],
    )
    feed.price = _TEA_PRICE["sargassum"]["baseline"]

    PR = Press("PR", ins=feed, outs=("pressed_cake", "pressate"))
    # Solid waste when this system runs standalone (pressed_cake is a leaf
    # stream of *this* bst.System). When embedded in a caller that feeds
    # pressed_cake into its own Mill instead, this price becomes inert --
    # BioSTEAM's TEA only counts price on that caller's own boundary streams.
    PR.outs[0].price = _TEA_PRICE["disposal_solid"]["baseline"]

    fresh_water = bst.Stream("biostimulant_fresh_water", Water=0.0, units="kg/hr")
    fresh_water.price = FRESH_WATER_PRICE_USD_PER_KG

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
    # Liquid waste: "a pure byproduct stream, never drawn from" (see docstring).
    PC.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

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
        outs=("biostimulant_product", "condensed_vapor"),
    )
    EV.outs[0].price = (
        _TEA_PRICE["biostimulant"]["baseline"] if biostimulant_price is None
        else float(biostimulant_price)
    )
    # Condensed by EV's own auxiliary condenser (V=0) -- leaves as a liquid,
    # not an off-gas -- so it's liquid waste, not exempt from pricing.
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
    create_tea(sys)

    return sys

def price_biostimulant_system(IRR: float | None = None) -> dict:
    """
    Parameters
    ----------
    IRR : float, optional
        Overrides the TEA's target IRR (data/tea.yaml `tea.IRR`, normally
        10%) before solving biostimulant_product's price. Pass IRR=0.0 to
        get biostimulant's zero-return, cost-basis price -- e.g. for use as
        a byproduct credit in a larger pathway (see
        systems._ad_biomethane_system.price_ad_biomethane_system(), which
        uses this to avoid crediting the other pathways with a return on
        the biostimulant subsystem's own capital on top of their own
        target-IRR product price).
    """
    from biorefineries.sabre._tea import solve_product_msp
    bst.main_flowsheet.clear()

    sys = create_biostimulant_system()
    sys.simulate()
    if IRR is not None:
        sys.TEA.IRR = float(IRR)

    product = sys.flowsheet.stream.biostimulant_product
    msp = solve_product_msp(tea=sys.TEA, product_stream=product)

    return {
        "label": "Biostimulant",
        "product_desc": "biostimulant liquid product",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        'sys': sys,
    }

if __name__ == '__main__':
    results = price_biostimulant_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")
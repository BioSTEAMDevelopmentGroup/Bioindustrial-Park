#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Standalone solvent-free isobutanol-ethanol-water separation system.

Recovers high-purity ethanol and isobutanol product streams from a corn
fermentation effluent (composition of stream ``P301-0`` of the isobutanol
biorefinery, scenario B) by heterogeneous-azeotropic distillation with a
decanter -- no extraction solvent and no isobutanol molecular sieve.

Design rationale (Dortmund-UNIFAC equilibria via thermosteam, 1 atm):

* Isobutanol-water forms a minimum-boiling heteroazeotrope at 89.8 C,
  y_IBO = 0.326 mol (66.6 wt%). Condensed and decanted at 40 C it splits
  into an organic phase at 87.5 wt% IBO (x_IBO = 0.543) and an aqueous
  phase at 10.7 wt% IBO (x_IBO = 0.029).
* Dilute isobutanol in water is extremely volatile (alpha ~ 26 vs. water;
  ethanol ~ 13), so the beer column carries both alcohols overhead.
* On the IBO-rich side of the azeotrope, water is the light key: distilling
  the decanted organic phase gives dry isobutanol as the *bottoms* product,
  with only the azeotropic overhead recycled to the decanter.

Flowsheet::

    broth ---> D101 beer column ----------------------> stillage (solids)
                 |  (LHK = (Isobutanol, Water): EtOH is a light non-key,
                 |   distillate held below the azeotrope, y_IBO = 0.25)
                 v
    M201 <--- MS201 water-rich retentate  <---- MS201 molecular sieve
      |  <--- S301 aqueous phase                 ^ (final EtOH polish)
      v                                          |
    D102 EtOH rectifier --- distillate vapor -- H202 superheater
      |   (LHK = (Ethanol, Water), y_top = 0.80805;    MS201 EtOH-rich
      |    IBO is a heavy non-key -> bottoms)          --> H201 --> ethanol
      v
    D103 IBO stripper ------ bottoms -------------> D103_bottoms (~IBO-free)
      |  (LHK = (Isobutanol, Water), y_IBO = 0.25)
      v distillate
    M301 --> H301 (40 C) --> S301 decanter
                               |  organic (87.5 wt% IBO)
                               v
                             D104 drying column (LHK = (Water, Isobutanol))
                               |    distillate (azeotrope) --> back to M301
                               v bottoms
                             H302 --> isobutanol product (>99.9 wt%)

Every column is a McCabe-Thiele ``BinaryDistillation`` operated on a
homogeneous, monotonic branch of the VLE (distillate specs kept clear of
the azeotropes); the azeotrope is crossed only by the decanter (rigorous
LLE). The molecular sieve is used only as the final ethanol polishing step
(stock-factory split, ~99.2 wt% product).

Feed adaptivity (0 -> ~200 g/L of either alcohol): the column design
specifications are re-derived from the current feed composition by process
specifications on every simulation pass, so one flowsheet covers the full
titer range without infeasible specs:

* D101 sizes its overhead water to serve BOTH downstream demands --
  ``water = max(IBO*(1-y1)/y1, EtOH*(1-xt)/xt)`` with ``y1`` the keys-basis
  distillate-IBO ceiling (0.25, azeotrope margin) and ``xt`` the maximum
  EtOH keys-fraction fed to the rectifier (0.55, well below its 0.808
  distillate spec) -- so the rectifier is never starved of water at low
  IBO titers. At the baseline the IBO term governs and the historical
  y_top = 0.25 design is recovered exactly.
* Bottoms specs scale with the feed (``x_bot = min(default, 0.1*z_feed)``),
  bounding light-key losses at ~10% instead of letting a fixed absolute
  spec cross the feed composition at low titers.
* D103's distillate spec rises with its feed (up to 0.31) so the spec
  always stays on the feasible side of the feed composition while keeping
  a margin to the 0.326 azeotrope.
* D101 additionally enforces an ethanol-water minimum-reflux floor: the
  IBO-water McCabe-Thiele alone yields a near-zero reflux, pricing the
  column's ~10x ethanol (light non-key) enrichment at nothing; the floor
  uses the same q-line construction on the ethanol-water pair, making the
  duty continuous with the zero-IBO ethanol-water design.
* A light key below ``min_key_flow`` (1e-2 kmol/hr ~ 0.5 kg/hr, i.e.
  absent at plant scale) has no column duty: D101 switches to a plain
  ethanol-water beer column
  (LHK = (Ethanol, Water), same overhead water -- the continuous limit of
  the IBO-mode design), while D102/D103 pass their feed to the bottoms
  outlet with design/cost skipped. Zero-throughput auxiliaries (decanter
  loop, ethanol polish train) are likewise skipped when their feeds
  vanish; recoveries at very low titers plateau near 1 - 2*(0.1) rather
  than collapsing.

Known model idealizations (inherited from ``BinaryDistillation``'s
boiling-point-based non-key routing): ethanol (light non-key in D101/D103)
reports 100% to distillates, and isobutanol (heavy non-key in D102)
reports 100% to the rectifier bottoms. Both routings match the physical
direction of travel but are sharper than reality.

Usage (standalone, without importing the biorefinery package)::

    import importlib.util
    spec = importlib.util.spec_from_file_location('separations', path)
    separations = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(separations)

    import biosteam as bst, thermosteam as tmo
    tmo.settings.set_thermo(separations.create_separation_chemicals())
    bst.main_flowsheet.set_flowsheet('IBO_EtOH_sep')
    feed = separations.create_scenario_B_feed()
    sys = separations.create_IBO_EtOH_separation_system(ins=feed)
    sys.simulate()

(``from biorefineries.isobutanol import separations`` also works; the package
import is build-free — the biorefinery build + baseline simulation run only at
an explicit ``isobutanol.load()`` call.)
"""

import numpy as np
import biosteam as bst
import thermosteam as tmo
from scipy.optimize import brentq

__all__ = (
    'create_separation_chemicals',
    'create_scenario_B_feed',
    'scenario_B_P301_flows',
    'create_IBO_EtOH_separation_system',
)

#%% Chemicals and reference feed

def create_separation_chemicals():
    """Return a compiled chemicals object: corn set + Isobutanol +
    AceticAcid (mirrors system.py, without the extraction solvent)."""
    from biorefineries import corn
    chems = [c for c in corn.chemicals.create_chemicals()]
    chems.append(tmo.Chemical('Isobutanol'))
    chems.append(tmo.Chemical('AceticAcid'))
    return chems

#: Fermentation effluent P301-0 of the isobutanol biorefinery at the
#: scenario B baseline (kg/hr; snapshot of the IBO_2026 stack, 2026-08-20,
#: taken under the former solvent-extraction recovery train, whose
#: purity-adjusted ethanol MPSP was 0.943 $/kg -- NOT the current baseline,
#: which is ~0.396 $/kg with this train integrated; commit 2797aa80).
#: EtOH titer 60.0 g/L, IBO titer 47.7 g/L.
scenario_B_P301_flows = {
    'Water':            116979.09,
    'Ethanol':            8946.50,
    'Glucose':              57.40,
    'Ash':                 795.22,
    'Yeast':              1867.72,
    'CaO':                   5.81,
    'TriOlein':           1931.24,
    'H2SO4':                48.43,
    'Fiber':              6060.68,
    'SolubleProtein':     1931.30,
    'InsolubleProtein':   2800.29,
    'Isobutanol':         7110.37,
    'AceticAcid':          135.12,
}

def create_scenario_B_feed(ID='broth', EtOH_mult=1.0, IBO_mult=1.0):
    """Return a Stream with the scenario B ``P301-0`` composition, with
    optional multipliers on the ethanol and isobutanol flows."""
    feed = tmo.Stream(ID, T=305.15, P=607950.0)
    for chem_ID, flow in scenario_B_P301_flows.items():
        feed.imass[chem_ID] = flow
    feed.imass['Ethanol'] *= EtOH_mult
    feed.imass['Isobutanol'] *= IBO_mult
    return feed

#%% Separation system factory

@bst.SystemFactory(
    ID='IBO_EtOH_separation_sys',
    ins=[dict(ID='broth')],
    outs=[dict(ID='ethanol_product'),
          dict(ID='isobutanol_product'),
          dict(ID='stillage'),
          dict(ID='D103_bottoms')],
)
def create_IBO_EtOH_separation_system(
        ins, outs,
        P=101325.0,
        beer_column_y_top_IBO=0.25,   # keys-basis CEILING; IBO-water azeotrope is 0.326
        beer_column_x_bot_IBO=1e-4,   # ceiling; scaled down with the feed at low titers
        rectifier_y_top_EtOH=0.80805, # stock ethanol-factory spec
        rectifier_x_bot_EtOH=3.9106e-6, # ceiling; scaled down with the feed at low titers
        rectifier_Rmin=1.2,           # EtOH-water tangent-pinch floor
        stripper_y_top_IBO=0.25,      # floor; raised toward the ceiling with richer feeds
        stripper_x_bot_IBO=1e-5,      # ceiling; scaled down with the feed at low titers
        drying_column_y_top_water=0.60, # IBO-rich branch; azeotrope is 0.674
        drying_column_x_bot_water=1e-4,
        decanter_T=313.15,
        product_T=305.15,
        rectifier_feed_x_EtOH_max=0.55, # max EtOH keys-fraction fed to D102 (D101
                                        # takes enough water overhead to ensure this)
        stripper_y_top_IBO_max=0.31,  # adaptive D103 distillate ceiling (< azeotrope)
        bottoms_key_loss_frac=0.10,   # x_bot cap as a fraction of feed keys-composition
        min_key_flow=1e-2,            # kmol/hr; a light key below this is absent, and a
                                      # unit fed less than this in total has no duty
                                      # (~0.5 kg/hr vs ~1e5 kg/hr of feed). Kept ABOVE
                                      # the system's molar flow tolerance so shut-off
                                      # recycle loops drain to exactly zero instead of
                                      # parking dust at the tolerance floor.
    ):
    broth, = ins
    ethanol_product, isobutanol_product, stillage, D103_bottoms = outs

    # Beer column: both alcohols overhead, solids/heavies to stillage.
    # LHK = (Isobutanol, Water) makes ethanol (Tb below both keys) a light
    # non-key -> 100% distillate, consistent with its dilute-aqueous
    # volatility (alpha ~ 13-26); the keys-basis distillate composition is
    # held below the IBO-water azeotrope so the McCabe-Thiele staircase
    # stays on the homogeneous water-rich branch.
    D101 = bst.BinaryDistillation(
        'D101', ins=broth, outs=('D101_distillate', stillage),
        LHK=('Isobutanol', 'Water'),
        y_top=beer_column_y_top_IBO, x_bot=beer_column_x_bot_IBO,
        k=1.2, P=P,
        partial_condenser=False,
    )

    # Ethanol rectifier: EtOH overhead as near-azeotropic vapor; IBO (heavy
    # non-key) to bottoms. Distillate/bottoms specs follow the stock
    # biorefineries.ethanol purification factory. Rmin enforces the
    # EtOH-water tangent-pinch minimum reflux (~1.0 for a near-bubble-point
    # feed), which biosteam's feed-line construction underestimates.
    M201 = bst.Mixer('M201', ins=(D101-0, '', ''))
    D102 = bst.BinaryDistillation(
        'D102', ins=M201-0, outs=('D102_distillate', 'D102_bottoms'),
        LHK=('Ethanol', 'Water'),
        y_top=rectifier_y_top_EtOH, x_bot=rectifier_x_bot_EtOH,
        k=1.2, Rmin=rectifier_Rmin, P=P,
        partial_condenser=True,
    )

    # Superheat vapor ahead of the molecular sieve (stock-factory practice).
    H202 = bst.HXutility('H202', ins=D102-0, T=115 + 273.15, V=1,
                         heat_only=True)

    # Molecular sieve, final ethanol polish only. Stock-factory splits and
    # stream convention: outs[0] is the water-rich retentate (recycled),
    # outs[1] the ethanol-rich product (~99.2 wt%). approx_duty is disabled
    # as in biorefineries.corn (the built-in duty approximation assumes
    # outs[1] is the small water-rich stream).
    MS201 = bst.MolecularSieve(
        'MS201', ins=H202-0, outs=('MS201_water_rich', 'MS201_ethanol_rich'),
        split=(2165.14/13356.04, 1280.06/1383.85),
        order=('Ethanol', 'Water'),
    )
    MS201.approx_duty = False
    MS201-0-1-M201  # water-rich retentate back to the rectifier

    H201 = bst.HXutility('H201', ins=MS201-1, outs=ethanol_product,
                         T=product_T, rigorous=True)

    # IBO stripper: re-concentrates isobutanol from the rectifier bottoms
    # toward the azeotrope; any residual ethanol is a light non-key ->
    # overhead, so the bottoms leave essentially alcohol-free.
    D103 = bst.BinaryDistillation(
        'D103', ins=D102-1, outs=('D103_distillate', D103_bottoms),
        LHK=('Isobutanol', 'Water'),
        y_top=stripper_y_top_IBO, x_bot=stripper_x_bot_IBO,
        k=1.2, P=P,
        partial_condenser=False,
    )

    # Decanter loop: condensed overheads cooled to 40 C split across the
    # IBO-water miscibility gap (rigorous LLE). The aqueous phase returns
    # to the rectifier feed -- this is also the escape path that prevents
    # ethanol from accumulating in the loop.
    M301 = bst.Mixer('M301', ins=(D103-0, ''))
    H301 = bst.HXutility('H301', ins=M301-0, T=decanter_T, rigorous=False)
    S301 = bst.LLESettler('S301', ins=H301-0,
                          outs=('S301_organic', 'S301_aqueous'),
                          top_chemical='Isobutanol')
    S301-1-2-M201  # aqueous phase back to the rectifier feed

    # IBO drying column: the organic phase lies on the IBO-rich branch of
    # the azeotrope, where water is the light key. Bottoms is dry
    # isobutanol; the azeotropic overhead returns to the decanter.
    D104 = bst.BinaryDistillation(
        'D104', ins=S301-0, outs=('D104_distillate', 'D104_bottoms'),
        LHK=('Water', 'Isobutanol'),
        y_top=drying_column_y_top_water, x_bot=drying_column_x_bot_water,
        k=1.2, P=P,
        partial_condenser=False,
    )
    D104-0-1-M301  # azeotropic overhead back to the decanter

    H302 = bst.HXutility('H302', ins=D104-1, outs=isobutanol_product,
                         T=product_T, rigorous=False)

    # ------------------------------------------------------------------
    # Feed-adaptive design specifications (0 -> ~200 g/L of either alcohol)
    # ------------------------------------------------------------------
    # Design specifications are re-derived from the current feed on every
    # simulation pass so that column specs always stay on the feasible side
    # of the feed composition. Where a stream genuinely vanishes (a light
    # key absent, or a unit fed nothing) the unit is run as nonexistent:
    # outlets emptied / feed passed through, design and cost skipped.
    # `Unit._setup` clears utilities and cost dictionaries at the start of
    # every simulation, so a skipped unit reports all-zero results with
    # nothing stale.

    def _design_and_cost_toggle(unit, _no_op=lambda: None):
        """Return a setter that enables/disables `unit`'s design & costing."""
        design, cost = unit._design, unit._cost
        def set_active(active):
            if active:
                unit._design, unit._cost = design, cost
            else:
                unit._design = unit._cost = _no_op
        return set_active

    D101_set_active = _design_and_cost_toggle(D101)
    D102_set_active = _design_and_cost_toggle(D102)
    D103_set_active = _design_and_cost_toggle(D103)

    def _empty_outs(unit):
        # BinaryDistillation's mass balance only assigns the chemical
        # indices of its CURRENT light/heavy(-non)key classification, so
        # outlet entries from a previous classification (an LHK switch, or
        # a bypass pass-through) would survive as stale flows -- observed
        # as the distillate "creating" the previous operating point's IBO
        # after D101 switched to ethanol-water mode. Start every active
        # run from empty outlets.
        for i in unit.outs: i.empty()

    # Overhead-water demands per kmol of each alcohol taken overhead in D101:
    # enough to hold the distillate at/below the keys-basis IBO ceiling, AND
    # enough to hold the rectifier feed at/below its target EtOH keys-fraction.
    water_per_IBO_ovhd = (1.0 - beer_column_y_top_IBO)/beer_column_y_top_IBO
    water_per_EtOH_ovhd = ((1.0 - rectifier_feed_x_EtOH_max)
                           /rectifier_feed_x_EtOH_max)

    def _EtOH_water_Rmin(column, E, W, W_ovhd):
        """Minimum reflux for the ethanol-water enrichment `column`
        performs on its ethanol (a light non-key when LHK is IBO-water):
        the IBO-water McCabe-Thiele staircase alone computes a near-zero
        reflux, because the ~10x ethanol enrichment the column also
        delivers costs nothing in the non-key idealization. This is the
        same q-line/equilibrium-intersection construction
        ``BinaryDistillation`` uses for its keys, applied to the
        ethanol-water pair (feed ``z = E/(E+W)``, distillate
        ``x_D = E/(E+W_ovhd)``), so the enforced-``Rmin`` design is
        CONTINUOUS with the LHK = (Ethanol, Water) design that takes over
        when isobutanol vanishes. Call only after ``column._run()`` (the
        feed-quality estimate uses the freshly mixed feed)."""
        q = column.get_feed_quality()
        if abs(q - 1.0) < 1e-4: q = 1.0 - 1e-4
        z = E/(E + W)
        x_D = E/(E + W_ovhd)
        bp = column.outs[1].get_bubble_point(('Ethanol', 'Water'))
        P = column.P
        q_line = lambda x: (q*x - z)/(q - 1.0)
        x_Rmin = brentq(lambda x: q_line(x)
                        - bp.solve_Ty(np.array((x, 1.0 - x)), P)[1][0],
                        0.0, 1.0)
        y_Rmin = q_line(x_Rmin)
        m = (y_Rmin - x_D)/(x_Rmin - x_D)
        if m >= 1.0: return 10.0  # numerical-pinch guard; not reached for x_D <= ~0.55
        return max(0.0, m/(1.0 - m))

    @D101.add_specification(run=False)
    def adapt_beer_column_to_feed():
        broth_feed, = D101.ins
        E = broth_feed.imol['Ethanol']
        I = broth_feed.imol['Isobutanol']
        W = broth_feed.imol['Water']
        if I > min_key_flow:
            # IBO present: IBO/water McCabe-Thiele with EtOH as light non-key
            # (dilute IBO is more volatile than EtOH). The distillate carries
            # whichever overhead-water demand governs; at the baseline the
            # IBO term does and y_top = beer_column_y_top_IBO is recovered.
            if D101.LHK != ('Isobutanol', 'Water'):
                D101.LHK = ('Isobutanol', 'Water')
            W_ovhd = max(I*water_per_IBO_ovhd, E*water_per_EtOH_ovhd)
            W_ovhd = min(W_ovhd, 0.8*W)  # y_top > z_feed guaranteed
            D101.y_top = I/(I + W_ovhd)
            D101.x_bot = min(beer_column_x_bot_IBO,
                             bottoms_key_loss_frac*I/(I + W))
            D101_set_active(True)
            _empty_outs(D101)
            D101._run()
            # Reflux floor for the ethanol enrichment (see _EtOH_water_Rmin);
            # only the mass balance above depends on the specs set so far,
            # so setting Rmin here still precedes the design calculations.
            # Never drop below biosteam's 0.3 default enforced Rmin: the
            # computed IBO-water Rmin can be NEGATIVE (dilute IBO is so
            # volatile that y*(z_feed) > y_top), and a zero enforced floor
            # would let R = 0 crash the operating-line construction.
            D101.Rmin = max(0.3, _EtOH_water_Rmin(D101, E, W, W_ovhd)
                                 if E > min_key_flow else 0.0)
        elif E > min_key_flow:
            # No IBO: a plain ethanol-water beer column. Same overhead water
            # as the IBO-mode limit I -> 0, so the design is continuous.
            if D101.LHK != ('Ethanol', 'Water'):
                D101.LHK = ('Ethanol', 'Water')
            W_ovhd = min(E*water_per_EtOH_ovhd, 0.8*W)
            D101.y_top = E/(E + W_ovhd)
            D101.x_bot = min(beer_column_x_bot_IBO,
                             bottoms_key_loss_frac*E/(E + W))
            D101.Rmin = 0.3  # the EtOH-water design computes its own Rmin
            D101_set_active(True)
            _empty_outs(D101)
            D101._run()
        else:
            # No alcohols at all: nothing to boil overhead.
            D101_set_active(False)
            D101.outs[0].empty()
            D101.outs[1].copy_like(broth_feed)

    @D102.add_specification(run=False)
    def adapt_rectifier_to_feed():
        rectifier_feed, = D102.ins
        E = rectifier_feed.imol['Ethanol']
        W = rectifier_feed.imol['Water']
        if E > min_key_flow:
            D102.x_bot = min(rectifier_x_bot_EtOH,
                             bottoms_key_loss_frac*E/(E + W))
            D102_set_active(True)
            _empty_outs(D102)
            D102._run()
        else:
            # No ethanol: no rectifier; feed continues to the IBO stripper.
            D102_set_active(False)
            D102.outs[0].empty()
            D102.outs[1].copy_like(rectifier_feed)

    @D103.add_specification(run=False)
    def adapt_stripper_to_feed():
        stripper_feed, = D103.ins
        I = stripper_feed.imol['Isobutanol']
        W = stripper_feed.imol['Water']
        if I > min_key_flow:
            z = I/(I + W)
            # Keep the distillate spec above the feed composition (feasible
            # enrichment) but below the ceiling's azeotrope margin.
            D103.y_top = max(stripper_y_top_IBO,
                             min(1.11*z, stripper_y_top_IBO_max))
            D103.x_bot = min(stripper_x_bot_IBO, bottoms_key_loss_frac*z)
            D103_set_active(True)
            _empty_outs(D103)
            D103._run()
        else:
            # No isobutanol: no stripper; feed leaves with the bottoms.
            D103_set_active(False)
            D103.outs[0].empty()
            D103.outs[1].copy_like(stripper_feed)

    def _add_low_flow_guard(unit):
        """Run `unit` as nonexistent (outlets emptied, design/cost skipped)
        when fed less than `min_key_flow` in total: zero-throughput
        equipment has no duty, and draining recycle loops (e.g. the
        decanter loop after D103 shuts off) reach exactly zero instead of
        parking dust at the flow tolerance."""
        set_active = _design_and_cost_toggle(unit)
        @unit.add_specification(run=False)
        def low_flow_guard():
            if sum([i.F_mol for i in unit.ins]) > min_key_flow:
                set_active(True)
                unit._run()
            else:
                set_active(False)
                for i in unit.outs: i.empty()

    for unit in (H202, MS201, H201, H301, S301, D104, H302):
        _add_low_flow_guard(unit)

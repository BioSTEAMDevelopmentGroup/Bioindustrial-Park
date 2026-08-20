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
    D103 IBO stripper ------ bottoms -------------> wastewater (~IBO-free)
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

(``from biorefineries.isobutanol import separations`` also works, but
triggers the package's import-time baseline simulation.)
"""

import biosteam as bst
import thermosteam as tmo

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
#: scenario B baseline (kg/hr; snapshot of the IBO_2026 stack, 2026-08-20;
#: purity-adjusted ethanol MPSP = 0.943 $/kg). EtOH titer 60.0 g/L,
#: IBO titer 47.7 g/L.
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
          dict(ID='wastewater')],
)
def create_IBO_EtOH_separation_system(
        ins, outs,
        P=101325.0,
        beer_column_y_top_IBO=0.25,   # keys-basis; IBO-water azeotrope is 0.326
        beer_column_x_bot_IBO=1e-4,
        rectifier_y_top_EtOH=0.80805, # stock ethanol-factory spec
        rectifier_x_bot_EtOH=3.9106e-6,
        rectifier_Rmin=1.2,           # EtOH-water tangent-pinch floor
        stripper_y_top_IBO=0.25,
        stripper_x_bot_IBO=1e-5,
        drying_column_y_top_water=0.60, # IBO-rich branch; azeotrope is 0.674
        drying_column_x_bot_water=1e-4,
        decanter_T=313.15,
        product_T=305.15,
    ):
    broth, = ins
    ethanol_product, isobutanol_product, stillage, wastewater = outs

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
        'D103', ins=D102-1, outs=('D103_distillate', wastewater),
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

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import math
import biosteam as bst
import thermosteam as tmo

from biorefineries.sabre.utils import get_solids_group_IDs, load_assumptions

__all__ = ('Press', 'EnzymaticPress', 'Mill')

# --- constants ---
KG_PER_METRIC_TON = 1000.0
KG_PER_DRY_TON = 907.18474  # US short ton (2000 lb)
HR_PER_DAY = 24.0

# Loaded assumptions
_PREPROCESSING = load_assumptions("preprocessing.yaml")
_PRESS = _PREPROCESSING["press"]
_ENZYMATIC_PRESS = _PREPROCESSING["enzymatic_press"]
_MILL = _PREPROCESSING["mill"]


class Press(bst.Unit):
    """
    Multi-train mechanical dewatering press. Splits wet biomass into:
      - pressed_cake: retains most solids + enough water to hit target solids wt%
      - pressate: remaining water + uncaptured solids + (optional) solubles

    Economics:
    - CAPEX can be set using an installed system correlation
    - Electricity set using kWh per dry ton TS

    Parameters
    ----------
    ins : stream
        Wet biomass feed.
    outs : tuple[stream, stream]
        Pressed cake and pressate.
    solids_IDs : Iterable[str], optional
        Chemical IDs captured (at `solids_capture_frac`) to the cake;
        everything else (except Water) is split by
        `solubles_to_pressate_frac`. Defaults to the "solids" chemical
        group when not given (see Sources below).
    solids_capture_frac : float
        Fraction of solids-ID mass captured to the cake (0-1).
    cake_solids_wt_frac : float
        Target cake solids weight fraction; entrained water is added to
        meet it.
    solubles_to_pressate_frac : float
        Fraction of non-solids, non-water mass routed to pressate (0-1);
        the remainder stays in the cake.
    power_kWh_per_dry_ton_TS : float
        Electricity intensity per dry ton of TS processed (preferred
        basis).
    ref_capacity_tph_wet : float
        Reference wet-throughput capacity (metric ton/h) of a single
        press train for CAPEX scaling; larger throughputs are split
        across parallel trains.
    capex_installed_ref_usd : float
        Installed cost at `ref_capacity_tph_wet`, scaled by
        `scale_exponent` for other train sizes.
    scale_exponent : float
        Power-law scaling exponent for installed CAPEX vs. train
        capacity.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/preprocessing.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2  # pressed_cake, pressate
    _F_BM_default = {"Press system": _PRESS["F_BM"]}

    def __init__(
        self, ID="", ins=None, outs=(),
        solids_IDs=None,  # if not given, defaults to the "solids" chemical group (see utils.get_solids_group_IDs)
        solids_capture_frac=_PRESS["solids_capture_frac"],
        cake_solids_wt_frac=_PRESS["cake_solids_wt_frac"],
        solubles_to_pressate_frac=_PRESS["solubles_to_pressate_frac"],

        # --- utilities ---
        power_kWh_per_dry_ton_TS=_PRESS["power_kWh_per_dry_ton_TS"],

        # --- costing ---
        ref_capacity_tph_wet=_PRESS["ref_capacity_tph_wet"],
        capex_installed_ref_usd=_PRESS["capex_installed_ref_usd"],
        scale_exponent=_PRESS["scale_exponent"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        if solids_IDs is None:
            solids_IDs = get_solids_group_IDs(self.chemicals)
        self.solids_IDs = tuple(solids_IDs)
        self.solids_capture_frac = float(solids_capture_frac)
        self.cake_solids_wt_frac = float(cake_solids_wt_frac)
        self.solubles_to_pressate_frac = float(solubles_to_pressate_frac)
        self.power_kWh_per_dry_ton_TS = power_kWh_per_dry_ton_TS
        self.ref_capacity_tph_wet = float(ref_capacity_tph_wet)
        self.capex_installed_ref_usd = float(capex_installed_ref_usd)
        self.scale_exponent = float(scale_exponent)

    def _get_mass(self, stream, chem_id):
        if chem_id not in stream.chemicals:
            return 0.0
        return float(stream.imass[chem_id])

    def _available_solids(self, stream):
        return [sid for sid in self.solids_IDs if sid in stream.chemicals]

    def _run(self):
        self._split_to_cake_and_pressate(self.ins[0], *self.outs)

    def _split_to_cake_and_pressate(self, feed, cake, pressate):
        """
        Core cake/pressate mass-split logic, factored out of `_run` so
        that `EnzymaticPress` can apply it to a hydrolyzed feed stream
        instead of `self.ins[0]` directly, without duplicating the split
        rules.

        The by-component split (everything except water) delegates to
        `Stream.split_to` -- the same native primitive
        `bst.units.Splitter`/`SolidsSeparator`/`ScrewPress` themselves use
        internally -- instead of a hand-rolled per-chemical loop.
        `energy_balance=False` keeps this an exact behavioral match for
        the original hand-rolled version (which built `cake`/`pressate`
        via bare `.imass[...]` assignment after `.empty()`, so it never
        copied `feed`'s T/P/phase either); using the default
        `energy_balance=True` would additionally copy `feed`'s T/P/phase
        onto both outlets, a real behavior change out of scope here.

        Water is NOT delegated to `bst.separations.adjust_moisture_content`
        (the native utility `SolidsSeparator` itself uses for this exact
        step): confirmed by direct testing that its `strict=False`
        fallback can manufacture water out of nothing when the feed
        doesn't have enough to hit the target moisture content (e.g. 5 kg
        available -> 35 kg assigned) -- a real mass-balance bug in that
        utility for the under-supply case. This keeps the original
        clamp-to-what's-available arithmetic instead.
        """
        cake.empty()
        pressate.empty()

        cake.phase = "l"
        pressate.phase = "l"

        solids = self._available_solids(feed)
        solids_set = set(solids)

        cap = min(max(self.solids_capture_frac, 0.0), 1.0)
        sol_to_p = min(max(self.solubles_to_pressate_frac, 0.0), 1.0)

        # Fraction of each chemical routed to `cake` (the 1st split_to
        # argument); everything not listed here (i.e. Water) defaults to
        # 0 -- routed to cake below via the moisture-content target
        # instead.
        split = {sid: cap for sid in solids}
        for chem_id in feed.chemicals.IDs:
            if chem_id in solids_set or chem_id == "Water":
                continue
            split[chem_id] = 1.0 - sol_to_p

        isplit = feed.chemicals.isplit(split)
        feed.split_to(cake, pressate, isplit.data, energy_balance=False)

        # Allocate water to hit cake solids wt% target
        TS_cake = sum(self._get_mass(cake, sid) for sid in solids)
        other_nonwater_cake = sum(
            self._get_mass(cake, i) for i in feed.chemicals.IDs
            if i not in solids_set and i != "Water"
        )

        f = self.cake_solids_wt_frac
        if TS_cake > 0 and 0 < f < 1:
            water_needed = TS_cake * (1 - f) / f - other_nonwater_cake
            water_needed = max(water_needed, 0.0)
        else:
            water_needed = 0.0

        water_avail = self._get_mass(feed, "Water")
        water_to_cake = min(water_needed, water_avail)

        cake.imass["Water"] = water_to_cake
        pressate.imass["Water"] = water_avail - water_to_cake

    def _design(self):
        feed = self.ins[0]
        solids = self._available_solids(feed)

        # TS through the press (kg/h) based on solids_IDs present in the stream
        TS_kgph = sum(self._get_mass(feed, sid) for sid in solids)
        dry_ton_per_hr_TS = TS_kgph / KG_PER_DRY_TON
        dtpd = dry_ton_per_hr_TS * HR_PER_DAY

        self.design_results["TS (kg/h)"] = TS_kgph
        self.design_results["TS (dry ton/h)"] = dry_ton_per_hr_TS
        self.design_results["Capacity (dry ton/day)"] = dtpd

        kW = float(self.power_kWh_per_dry_ton_TS) * dry_ton_per_hr_TS
        self.power_utility(kW)

    def _cost(self):
        feed = self.ins[0]

        wet_tph = feed.F_mass / 1000.0  # metric ton/h
        Q0 = self.ref_capacity_tph_wet
        C0 = self.capex_installed_ref_usd
        n = self.scale_exponent

        if wet_tph <= 0:
            capex = 0.0
            N = 0
            Q_each = 0.0
        else:
            N = max(1, math.ceil(wet_tph / Q0))
            Q_each = wet_tph / N
            capex = N * C0 * (Q_each / Q0) ** n

        self.design_results["Wet throughput (tph)"] = wet_tph
        self.design_results["Number of press trains"] = int(N)
        self.design_results["Train throughput (tph)"] = Q_each
        self.design_results["Installed CAPEX ($)"] = capex

        self.baseline_purchase_costs["Press system"] = capex


class EnzymaticPress(Press):
    """
    `Press` subclass that hydrolyzes a fraction of Glucan and Alginate to
    their soluble monomers (Glucose, AlginateMonomer) using an added
    enzyme cocktail (cellulase + alginate lyase), before the inherited
    cake/pressate mass split runs. Once hydrolyzed, Glucose/AlginateMonomer
    are not in the "solids" chemical group (see
    `biorefineries.sabre._chemicals`), so they follow the existing
    "solubles" path into the pressate, while unconverted Glucan/Alginate
    stay solid and are captured to the cake at `solids_capture_frac`, same
    as Ash/Protein/Fucoidan/OtherSolids.

    This unit replaces both a dilute-acid pretreatment reactor and a
    separate enzymatic saccharification step -- there is no intermediate
    oligomer chemistry; hydrolysis goes directly to the monomer form.

    Parameters
    ----------
    ins : tuple[stream, stream]
        Wet biomass feed and enzyme cocktail feed. The enzyme stream's
        flow is set internally by `_run` (dose scales with Glucan +
        Alginate mass in the feed); any pre-existing content is
        overwritten.
    outs : tuple[stream, stream]
        Pressed cake and pressate (same as `Press`).
    hydrolysis_conversion : float
        Fraction of Glucan and Alginate hydrolyzed to Glucose/
        AlginateMonomer (0-1), applied identically to both substrates.
    enzyme_loading_mg_per_g_substrate : float
        Enzyme dose per g of Glucan + Alginate in the feed, before the
        excess allowance.
    enzyme_excess_frac : float
        Fractional excess added on top of `enzyme_loading_mg_per_g_substrate`
        (e.g. 0.10 for 10% excess).
    **kwargs
        Forwarded to `Press.__init__` (solids_capture_frac,
        cake_solids_wt_frac, power_kWh_per_dry_ton_TS, etc.) and, from
        there, to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/preprocessing.yaml (`enzymatic_press`) for the default
    values and references.
    """

    _N_ins = 2  # wet biomass feed, enzyme cocktail

    def __init__(
        self, ID="", ins=None, outs=(),
        hydrolysis_conversion=_ENZYMATIC_PRESS["hydrolysis_conversion"],
        enzyme_loading_mg_per_g_substrate=_ENZYMATIC_PRESS["enzyme_loading_mg_per_g_substrate"],
        enzyme_excess_frac=_ENZYMATIC_PRESS["enzyme_excess_frac"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.hydrolysis_conversion = float(hydrolysis_conversion)
        self.enzyme_loading_mg_per_g_substrate = float(enzyme_loading_mg_per_g_substrate)
        self.enzyme_excess_frac = float(enzyme_excess_frac)

        # Internal working stream: feed + enzyme, post-hydrolysis. Kept
        # separate from self.ins[0] so the upstream feed stream is never
        # mutated in place.
        self._hydrolyzed_feed = bst.Stream(None, thermo=self.thermo)

        self.hydrolysis_rxns = tmo.ParallelReaction([
            tmo.Reaction('Glucan + Water -> Glucose', reactant='Glucan',
                         X=self.hydrolysis_conversion),
            tmo.Reaction('Alginate + Water -> AlginateMonomer', reactant='Alginate',
                         X=self.hydrolysis_conversion),
        ])

    def _run(self):
        feed, enzyme = self.ins
        substrate_kgph = feed.imass['Glucan'] + feed.imass['Alginate']
        dose_frac = (self.enzyme_loading_mg_per_g_substrate / 1000.0
                     * (1.0 + self.enzyme_excess_frac))
        enzyme.empty()
        enzyme.phase = "l"
        enzyme.imass['Enzyme'] = dose_frac * substrate_kgph

        combined = self._hydrolyzed_feed
        combined.mix_from([feed, enzyme])
        self.hydrolysis_rxns(combined)

        self._split_to_cake_and_pressate(combined, *self.outs)


class Mill(bst.Unit):
    """
    Multi-train hammer mill. Applies explicit mass loss during milling/shredding.
    Sends lost material to a 'losses' stream (same composition as feed).

    Economics:
    - Electricity: kWh per dry ton of dry material
    - CAPEX: anchor scaling for hammer mill

    Parameters
    ----------
    ins : stream
        Feed to be milled.
    outs : tuple[stream, stream]
        Milled biomass and milling losses (same composition as feed).
    loss_frac : float
        Fraction of feed mass lost during milling/shredding (0-1), sent
        to the losses stream.
    power_kWh_per_dry_ton_dry : float
        Electricity intensity per dry ton of dry material processed.
    ref_capacity_dry_ton_per_hr : float
        Reference dry-throughput capacity (dry ton/h) of a single mill
        for CAPEX scaling; larger throughputs are split across parallel
        mills.
    purchase_cost_ref_usd : float
        Purchase cost at `ref_capacity_dry_ton_per_hr`, scaled by
        `scale_exponent` for other mill sizes.
    scale_exponent : float
        Power-law scaling exponent for purchase cost vs. mill capacity.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/preprocessing.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2  # milled_biomass, milling_losses
    _F_BM_default = {"Hammer mill": _MILL["F_BM"]}

    def __init__(
        self, ID="", ins=None, outs=(),
        loss_frac=_MILL["loss_frac"],

        # utilities
        power_kWh_per_dry_ton_dry=_MILL["power_kWh_per_dry_ton_dry"],

        # costing
        ref_capacity_dry_ton_per_hr=_MILL["ref_capacity_dry_ton_per_hr"],
        purchase_cost_ref_usd=_MILL["purchase_cost_ref_usd"],
        scale_exponent=_MILL["scale_exponent"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.loss_frac = float(loss_frac)

        self.power_kWh_per_dry_ton_dry = power_kWh_per_dry_ton_dry

        self.ref_capacity_dry_ton_per_hr = float(ref_capacity_dry_ton_per_hr)
        self.purchase_cost_ref_usd = float(purchase_cost_ref_usd)
        self.scale_exponent = float(scale_exponent)

    def _run(self):
        feed = self.ins[0]
        milled, losses = self.outs
        milled.empty()
        losses.empty()

        milled.phase = feed.phase
        losses.phase = feed.phase

        # Same fraction lost across every chemical, i.e. a plain scalar
        # split -- delegated to the native Stream.split_to primitive
        # (same one bst.units.Splitter uses internally) instead of a
        # hand-rolled per-chemical loop. energy_balance=False matches the
        # original's behavior exactly (T/P were never copied from feed
        # here, only phase, set explicitly above).
        lf = min(max(self.loss_frac, 0.0), 1.0)
        feed.split_to(milled, losses, 1.0 - lf, energy_balance=False)

    def _design(self):
        feed = self.ins[0]
        water_kgph = float(feed.imass["Water"]) if "Water" in feed.chemicals.IDs else 0.0
        dry_kgph = max(feed.F_mass - water_kgph, 0.0)

        dry_ton_per_hr = dry_kgph / KG_PER_DRY_TON
        self.design_results["Dry mass (kg/h)"] = dry_kgph
        self.design_results["Dry throughput (dry ton/h)"] = dry_ton_per_hr

        kW = float(self.power_kWh_per_dry_ton_dry) * dry_ton_per_hr
        self.power_utility(kW)

    def _cost(self):
        dry_ton_per_hr = float(self.design_results.get("Dry throughput (dry ton/h)", 0.0))

        Q0 = max(self.ref_capacity_dry_ton_per_hr, 1e-9)
        C0 = self.purchase_cost_ref_usd
        n = self.scale_exponent

        # Size with parallel mills if needed
        N = max(1, math.ceil(dry_ton_per_hr / Q0)) if dry_ton_per_hr > 0 else 1
        Q_each = dry_ton_per_hr / N if N else dry_ton_per_hr

        purchase_each = C0 * (max(Q_each, 1e-9) / Q0) ** n
        purchase_total = N * purchase_each

        self.design_results["Number of mills"] = N
        self.design_results["Mill capacity each (dry ton/h)"] = Q_each
        self.design_results["Purchase cost each ($)"] = purchase_each

        self.baseline_purchase_costs["Hammer mill"] = purchase_total

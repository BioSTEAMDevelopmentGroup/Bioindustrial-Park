# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Preprocessing: mechanical dewatering (press) and milling / size reduction.
"""

import math
import biosteam as bst

__all__ = ('Press', 'Mill')

# --- constants ---
KG_PER_METRIC_TON = 1000.0
KG_PER_DRY_TON = 907.18474  # US short ton (2000 lb)
HR_PER_DAY = 24.0


class Press(bst.Unit):
    """
    Multi-train mechanical dewatering press. Constructed in three places in
    this codebase (`systems/_ad_biogas_system.py`, `systems/_ad_vfa_system.py`,
    `systems/_integrated_system.py`), all with the same
    `assumptions.yaml`-sourced parameters.

    Splits wet biomass into:
      - pressed_cake: retains most solids + enough water to hit target solids wt%
      - pressate: remaining water + uncaptured solids + (optional) solubles

    Economics:
    - CAPEX can be set using an installed system correlation
    - Electricity set using kWh per dry ton TS (preferred) or legacy kWh per wet ton (optional)

    Sources
    -------
    solids_capture_frac (0.98), solubles_to_pressate_frac (1.0),
    capex_installed_ref_usd ($3,100,000), scale_exponent (0.60),
    power_kWh_per_dry_ton_TS (5.0):
        assumptions.yaml `preprocessing.press` section.
        `capex_installed_ref_usd` sourced to `sources.capex`: "City of
        Bellingham Biosolids Facility Planning Report, TM 12: Solids
        Screening, Thickening, and Dewatering Technology Evaluation (2021)
        -- installed screw press system anchor based on municipal screw
        press dewatering cost tables; adapted to this model using
        dry-solids throughput scaling." `power_kWh_per_dry_ton_TS` sourced
        to `sources.power`: "HUBER screw press case study (2023) -- use ~5
        kWh/t dry solids for press-only electricity."
    cake_solids_wt_frac (0.15):
        assumptions.yaml `preprocessing.press.cake_solids_wt_frac`. NOTE:
        the yaml's own `sources.performance` citation (Vincent Corporation
        seaweed pressing notes; brown macroalgae dewatering literature)
        states "Sargassum-specific pressing data indicate cake solids
        around 26-29 wt% TS; 27 wt% TS selected as central baseline" -- but
        the actual modeled value is 0.15 (15%), not 0.27. This mismatch was
        found during this conversion and deliberately NOT changed (fixing
        it would alter simulation output, out of scope for a conversion
        pass) -- flagged in CONVERSION_NOTES.md.
    ref_capacity_tph_wet (50.0) -- basis mismatch, flagged not fixed:
        No yaml counterpart. assumptions.yaml actually declares
        `preprocessing.press.basis: dry_tph` with
        `ref_capacity_tph_dry: 0.953`, but this class only ever implements
        a wet-tph reference capacity -- the yaml's dry-basis fields are
        never read by any code path in this codebase. This 50.0 default
        (matching what every real call site already falls back to, since
        no call site supplies a dry-basis equivalent either) was not
        changed; see CONVERSION_NOTES.md.
    solids_IDs:
        assumptions.yaml `preprocessing.press.solids_IDs`. The
        pre-conversion default additionally listed `Xylan`, `Mannan`,
        `Galactan` -- real thermo chemicals, just absent from this yaml
        list, dropped to match (same pattern seen in `AnaerobicDigester`
        and `AcidogenicDigester`). Unlike `DigestateScrewPress`, this
        parameter *is* actively used by `_run()`/`_design()`.
    capex_model / "pca_screwpress_curve":
        REMOVED. The pre-conversion `_cost()` dispatched on a
        `capex_model` string with two branches ("scaled_anchor" and
        "pca_screwpress_curve"). Grep across this codebase and
        assumptions.yaml confirmed "pca_screwpress_curve" is never set
        anywhere -- yaml's `preprocessing.press.capex_model` is always
        "scaled_anchor", the only branch ever exercised at any of the 3
        real call sites. Collapsed to a single hardcoded formula (the old
        "scaled_anchor" branch), matching the "hardcoded anchor values
        instead of a capex_model string-dispatch pattern" convention used
        for `AnaerobicDigester`/`AcidogenicDigester`.
    """

    _N_ins = 1
    _N_outs = 2  # pressed_cake, pressate

    def __init__(
        self, ID="", ins=None, outs=(),
        solids_IDs=(
            "Glucan",
            "Arabinan",
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
            "Lignin",
            "Ash",
        ),
        solids_capture_frac=0.98,
        cake_solids_wt_frac=0.15,
        solubles_to_pressate_frac=1.0,

        # --- utilities ---
        power_kWh_per_dry_ton_TS=5.0,
        power_kWh_per_ton_wet=None,

        # --- costing ---
        F_BM=1.0,
        ref_capacity_tph_wet=50.0,
        capex_installed_ref_usd=3_100_000.0,
        scale_exponent=0.6,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.solids_IDs = tuple(solids_IDs)
        self.solids_capture_frac = float(solids_capture_frac)
        self.cake_solids_wt_frac = float(cake_solids_wt_frac)
        self.solubles_to_pressate_frac = float(solubles_to_pressate_frac)

        self.power_kWh_per_dry_ton_TS = power_kWh_per_dry_ton_TS
        self.power_kWh_per_ton_wet = power_kWh_per_ton_wet

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Press system"] = float(F_BM)
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
        feed = self.ins[0]
        cake, pressate = self.outs
        cake.empty()
        pressate.empty()

        cake.phase = "l"
        pressate.phase = "l"

        solids = self._available_solids(feed)

        # Split defined solids by capture
        cap = min(max(self.solids_capture_frac, 0.0), 1.0)
        for sid in solids:
            m = self._get_mass(feed, sid)
            m_cake = cap * m
            cake.imass[sid] = m_cake
            pressate.imass[sid] = m - m_cake

        # Partition everything else except water
        sol_to_p = min(max(self.solubles_to_pressate_frac, 0.0), 1.0)
        for chem_id in feed.chemicals.IDs:
            if chem_id in solids or chem_id == "Water":
                continue
            m = self._get_mass(feed, chem_id)
            m_p = sol_to_p * m
            pressate.imass[chem_id] += m_p
            cake.imass[chem_id] += (m - m_p)

        # Allocate water to hit cake solids wt% target
        TS_cake = sum(self._get_mass(cake, sid) for sid in solids)
        other_nonwater_cake = sum(
            self._get_mass(cake, i) for i in feed.chemicals.IDs
            if i not in solids and i != "Water"
        )

        f = self.cake_solids_wt_frac
        if TS_cake > 0 and 0 < f < 1:
            water_needed = TS_cake * (1 - f) / f - other_nonwater_cake
            water_needed = max(water_needed, 0.0)
        else:
            water_needed = 0.0

        water_avail = self._get_mass(feed, "Water")
        water_to_cake = min(water_needed, water_avail)

        cake.imass["Water"] += water_to_cake
        pressate.imass["Water"] += (water_avail - water_to_cake)

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

        # Power (preferred TS basis; fallback wet basis)
        if self.power_kWh_per_dry_ton_TS is not None:
            kW = float(self.power_kWh_per_dry_ton_TS) * dry_ton_per_hr_TS
            self.power_utility(kW)
        elif self.power_kWh_per_ton_wet is not None:
            wet_ton_per_hr = feed.F_mass / KG_PER_METRIC_TON
            kW = float(self.power_kWh_per_ton_wet) * wet_ton_per_hr
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


class Mill(bst.Unit):
    """
    Multi-train hammer mill. Constructed in three places in this codebase
    (`systems/_ad_biogas_system.py`, `systems/_ad_vfa_system.py`,
    `systems/_integrated_system.py`), all with the same
    `assumptions.yaml`-sourced parameters.

    Applies explicit mass loss during milling/shredding.
    Sends lost material to a 'losses' stream (same composition as feed).

    Economics:
    - Electricity: kWh per dry ton of dry material (preferred) or legacy wet ton basis
    - CAPEX: anchor scaling for hammer mill (purchase cost), then apply install factor

    Sources
    -------
    loss_frac (0.03):
        assumptions.yaml `preprocessing.mill.loss_frac`. Cited in yaml's
        `preprocessing.sources.loss` (a sibling of `press`/`mill` under
        `preprocessing`, but content-wise specific to Mill -- see note
        below): "INL preprocessing/logistics accounting -- dry matter loss
        across preprocessing is typically on the order of ~1-5%; 3%
        selected as central baseline."
    power_kWh_per_dry_ton_dry (25.0):
        assumptions.yaml `preprocessing.mill.power_kWh_per_dry_ton_dry`.
        Cited in `preprocessing.sources.power`: "Oyedeji et al. (2020);
        NREL FY11 grinding basis -- herbaceous biomass coarse grinding
        typically falls in the ~11-28 kWh/t range; 25 kWh/t selected as a
        defensible central proxy for dried Sargassum."
    ref_capacity_dry_ton_per_hr (10.0), purchase_cost_ref_usd ($206,400),
    install_factor (1.8), scale_exponent (0.6):
        assumptions.yaml `preprocessing.mill` section.
        `purchase_cost_ref_usd` cited in `preprocessing.sources.capex`:
        "INL hammer mill cost anchor
        (https://inldigitallibrary.inl.gov/sites/sti/sti/7323596.pdf) --
        hammer mill purchase cost anchor $206,400 used for capacity
        scaling with exponent 0.6 and install factor."
    Note on yaml structure: `preprocessing.sources` is indented as a
    sibling of `press:`/`mill:` (not nested inside `mill:`), but every
    citation under it is specific to Mill's own parameters (grinding
    power, hammer mill capex, preprocessing dry-matter loss) -- verified
    by reading the citation text itself, not inferred from yaml
    indentation, to avoid misattributing a Press citation here.
    capex_model / "inl_hammermill_anchor":
        REMOVED. The pre-conversion `_cost()` gated its only cost formula
        behind a `capex_model == "inl_hammermill_anchor"` string check.
        yaml's `preprocessing.mill.capex_model` is always exactly that
        string at every real call site (confirmed by grep across this
        codebase and assumptions.yaml) -- there was never a second
        branch, so the check was always true in practice. Collapsed to
        run unconditionally, same "hardcoded anchor, no string dispatch"
        convention used for `AnaerobicDigester`/`AcidogenicDigester`/`Press`.
    """

    _N_ins = 1
    _N_outs = 2  # milled_biomass, milling_losses

    def __init__(
        self, ID="", ins=None, outs=(),
        loss_frac=0.03,

        # utilities
        power_kWh_per_dry_ton_dry=25.0,  # preferred
        power_kWh_per_ton_wet=None,      # legacy optional

        # costing
        ref_capacity_dry_ton_per_hr=10.0,
        purchase_cost_ref_usd=206400.0,
        install_factor=1.8,
        scale_exponent=0.6,
        F_BM=1.0,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.loss_frac = float(loss_frac)

        self.power_kWh_per_dry_ton_dry = power_kWh_per_dry_ton_dry
        self.power_kWh_per_ton_wet = power_kWh_per_ton_wet

        self.ref_capacity_dry_ton_per_hr = float(ref_capacity_dry_ton_per_hr)
        self.purchase_cost_ref_usd = float(purchase_cost_ref_usd)
        self.install_factor = float(install_factor)
        self.scale_exponent = float(scale_exponent)

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Hammer mill"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        milled, losses = self.outs
        milled.empty()
        losses.empty()

        milled.phase = feed.phase
        losses.phase = feed.phase

        lf = min(max(self.loss_frac, 0.0), 1.0)
        for chem_id in feed.chemicals.IDs:
            m = float(feed.imass[chem_id])
            m_loss = lf * m
            losses.imass[chem_id] = m_loss
            milled.imass[chem_id] = m - m_loss

    def _design(self):
        feed = self.ins[0]
        water_kgph = float(feed.imass["Water"]) if "Water" in feed.chemicals.IDs else 0.0
        dry_kgph = max(feed.F_mass - water_kgph, 0.0)

        dry_ton_per_hr = dry_kgph / KG_PER_DRY_TON
        self.design_results["Dry mass (kg/h)"] = dry_kgph
        self.design_results["Dry throughput (dry ton/h)"] = dry_ton_per_hr

        # Power: prefer dry-basis, fallback wet-basis
        if self.power_kWh_per_dry_ton_dry is not None:
            kW = float(self.power_kWh_per_dry_ton_dry) * dry_ton_per_hr
            self.power_utility(kW)
        elif self.power_kWh_per_ton_wet is not None:
            wet_ton_per_hr = feed.F_mass / KG_PER_METRIC_TON
            kW = float(self.power_kWh_per_ton_wet) * wet_ton_per_hr
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
        installed_total = N * purchase_each * self.install_factor

        self.design_results["Number of mills"] = N
        self.design_results["Mill capacity each (dry ton/h)"] = Q_each
        self.design_results["Purchase cost each ($)"] = purchase_each
        self.design_results["Installed cost total ($)"] = installed_total

        self.baseline_purchase_costs["Hammer mill"] = installed_total

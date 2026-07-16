# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Enzymatic pretreatment before methanogenic AD.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('EnzymaticPretreatment',)


class EnzymaticPretreatment(bst.Unit):
    """
    Enzymatic pretreatment reactor, one of the optional pre-AD pretreatment
    stages in `systems._ad_biogas_system._build_methanogenic_pathway`.
    Constructed in exactly one place in this codebase (that function),
    reached only when the selected pretreatment case has `kind: enzymatic`
    or `kind: combined_PE`/`combined_PTE` (which include an enzymatic
    stage in sequence).

    Not converted to an `@cost` decorator: `capex_usd` is a flat,
    non-scaling placeholder cost (not a function of any design basis), so
    the decorator's power-law machinery wouldn't apply cleanly -- kept as
    the simplest possible custom `_cost()` instead.

    Sources
    -------
    capex_usd ($7,280,000):
        assumptions.yaml `ad_pretreatment_cases.enzymatic.sources.capex`,
        key `nrel_2015_saccharification_tank`, report NREL/TP-5100-62498.
    enzyme_price_usd_per_kg ($7.00):
        assumptions.yaml `ad_pretreatment_cases.enzymatic.sources.
        enzyme_price`, key `nrel_2015_enzyme_cost`, report
        NREL/TP-5100-62498 (same report as capex, different section).
    temperature_K, residence_time_hr, enzyme_dose_kg_per_kg_dry_feed
    (0.005), treated_fraction (0.25), enzyme_recycle_factor (2.0):
        assumptions.yaml `ad_pretreatment_cases.enzymatic.enzymatic`
        section. No named literature source for these five in yaml --
        only `capex` and `enzyme_price` have `sources:` entries. The
        yaml's `sources.methane_yield`/`biogas_composition`/
        `vs_destruction` entries (Chikani-Cabrera et al. 2022) describe the
        *downstream AD effect* of this pretreatment case, not this
        reactor's own operating parameters -- not cited here to avoid
        misattributing an AD-outcome citation to a pretreatment-reactor
        parameter.

    KNOWN BUG, VERIFIED, NOT FIXED (see CONVERSION_NOTES.md): `_cost()`
    computes `annual_maintenance` into `design_results` but never assigns
    it to `self.add_OPEX`, unlike `HeatingPretreatment` (which does this
    correctly). Since `capex_usd` is nonzero at this unit's real call site,
    this silently omits real maintenance cost from every pretreatment case
    that includes an enzymatic stage (enzymatic, combined_PE, combined_PTE)
    -- understating VOC/MSP for those cases. NOT fixed here: those exact
    cases are pinned by `tests.py` assertions
    (`test_ad_tea_final_combined_pe_near_zero`,
    `test_ad_feed_tea_build_case`'s "enzymatic" case) captured from the
    current, buggy behavior -- fixing this would change simulated output,
    out of scope for a conversion pass. Needs explicit sign-off and a
    `tests.py` baseline update before ever fixing.
    """

    _N_ins = 1
    _N_outs = 1  # enzyme_treated_biomass

    def __init__(
        self, ID="", ins=None, outs=(),
        temperature_K=308.15,
        residence_time_hr=24.0,
        enzyme_dose_kg_per_kg_dry_feed=0.005,
        treated_fraction=0.25,
        enzyme_recycle_factor=2.0,
        slurry_density_kg_per_m3=1000.0,

        # economics
        capex_usd=7_280_000.0,
        enzyme_price_usd_per_kg=7.0,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.temperature_K = float(temperature_K)
        self.residence_time_hr = float(residence_time_hr)
        self.enzyme_dose_kg_per_kg_dry_feed = float(enzyme_dose_kg_per_kg_dry_feed)
        self.treated_fraction = float(treated_fraction)
        self.enzyme_recycle_factor = float(enzyme_recycle_factor)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)

        self.capex_usd = float(capex_usd)
        self.enzyme_price_usd_per_kg = float(enzyme_price_usd_per_kg)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Enzymatic pretreatment"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(feed)
        treated.phase = feed.phase

        water = float(feed.imass["Water"]) if "Water" in feed.chemicals.IDs else 0.0
        dry_mass = max(feed.F_mass - water, 0.0)

        treated_fraction = min(max(self.treated_fraction, 0.0), 1.0)
        recycle = max(self.enzyme_recycle_factor, 1.0)

        treated_dry_mass = dry_mass * treated_fraction
        enzyme_kgph = self.enzyme_dose_kg_per_kg_dry_feed * treated_dry_mass / recycle

        self._dry_mass_kgph = dry_mass
        self._treated_dry_mass_kgph = treated_dry_mass
        self._enzyme_kgph = enzyme_kgph

        if "Enzyme" in treated.chemicals.IDs:
            treated.imass["Enzyme"] += enzyme_kgph

    def _design(self):
        feed = self.ins[0]

        m_kgph = float(feed.F_mass)
        vol_m3ph = m_kgph / self.slurry_density_kg_per_m3
        V = vol_m3ph * self.residence_time_hr

        dry_mass = getattr(self, "_dry_mass_kgph", 0.0)
        treated_dry_mass = getattr(self, "_treated_dry_mass_kgph", 0.0)
        enzyme_kgph = getattr(self, "_enzyme_kgph", 0.0)
        enzyme_cost_usd_per_hr = enzyme_kgph * self.enzyme_price_usd_per_kg

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Temperature (K)"] = self.temperature_K
        self.design_results["Dry feed basis (kg/h)"] = dry_mass
        self.design_results["Treated dry feed basis (kg/h)"] = treated_dry_mass
        self.design_results["Treated fraction"] = self.treated_fraction
        self.design_results["Enzyme recycle factor"] = self.enzyme_recycle_factor
        self.design_results["Enzyme addition (kg/h)"] = enzyme_kgph
        self.design_results["Enzyme cost ($/h)"] = enzyme_cost_usd_per_hr

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Enzymatic pretreatment"] = capex

        enzyme_kgph = getattr(self, "_enzyme_kgph", 0.0)
        enzyme_cost_usd_per_hr = enzyme_kgph * self.enzyme_price_usd_per_kg

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        add_opex = {}
        if enzyme_cost_usd_per_hr > 0:
            add_opex["Enzyme reagent cost"] = enzyme_cost_usd_per_hr
        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            add_opex["Enzymatic pretreatment maintenance"] = (
                annual_maintenance / operating_hours_per_yr
            )
        if add_opex:
            self.add_OPEX = add_opex

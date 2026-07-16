# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Hydrogen peroxide pretreatment before methanogenic AD.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('PeroxidePretreatment',)


class PeroxidePretreatment(bst.Unit):
    """
    Hydrogen peroxide pretreatment reactor, one of the optional pre-AD
    pretreatment stages in
    `systems._ad_biogas_system._build_methanogenic_pathway`. Constructed in
    exactly one place in this codebase (that function), reached only when
    the selected pretreatment case has `kind: peroxide` or
    `kind: combined_PE`/`combined_PTE` (which include a peroxide stage in
    sequence).

    Not converted to an `@cost` decorator: same reasoning as
    `EnzymaticPretreatment` -- `capex_usd` is a flat, non-scaling
    placeholder cost, not a function of any design basis.

    Sources
    -------
    capex_usd ($1,000,000):
        assumptions.yaml `ad_pretreatment_cases.peroxide.sources.capex`,
        key `nrel_2011_conditioning_tank_anchor`, report
        NREL/TP-5100-47764.
    h2o2_wt_frac_on_dry_feed, temperature_K, residence_time_hr,
    h2o2_price_usd_per_kg ($0.37/kg):
        assumptions.yaml `ad_pretreatment_cases.peroxide.peroxide` section.
        No named literature source for these four in yaml -- only `capex`
        has a `sources:` entry. Same caveat as `EnzymaticPretreatment`:
        yaml's `sources.methane_yield`/`biogas_composition`/
        `vs_destruction` entries describe the downstream AD effect of this
        pretreatment case, not this reactor's own operating parameters --
        not cited here.

    KNOWN BUG, VERIFIED, NOT FIXED (see CONVERSION_NOTES.md): same
    missing `add_OPEX` wiring for maintenance as `EnzymaticPretreatment` --
    `_cost()` computes `annual_maintenance` but never assigns it to
    `self.add_OPEX`. Understates VOC/MSP for every pretreatment case that
    includes a peroxide stage (peroxide, combined_PE, combined_PTE). NOT
    fixed here for the same reason: pinned by `tests.py` assertions
    captured from the current, buggy behavior.
    """

    _N_ins = 1
    _N_outs = 1  # peroxide_treated_biomass

    def __init__(
        self, ID="", ins=None, outs=(),
        h2o2_wt_frac_on_dry_feed=0.025,
        temperature_K=298.15,
        residence_time_hr=2.0,
        slurry_density_kg_per_m3=1000.0,

        # economics
        capex_usd=1_000_000.0,
        h2o2_price_usd_per_kg=0.37,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.h2o2_wt_frac_on_dry_feed = float(h2o2_wt_frac_on_dry_feed)
        self.temperature_K = float(temperature_K)
        self.residence_time_hr = float(residence_time_hr)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)

        self.capex_usd = float(capex_usd)
        self.h2o2_price_usd_per_kg = float(h2o2_price_usd_per_kg)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Peroxide pretreatment"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(feed)
        treated.phase = feed.phase

        water = float(feed.imass["Water"]) if "Water" in feed.chemicals.IDs else 0.0
        dry_mass = max(feed.F_mass - water, 0.0)

        h2o2_kgph = self.h2o2_wt_frac_on_dry_feed * dry_mass
        self._dry_mass_kgph = dry_mass
        self._h2o2_kgph = h2o2_kgph

        if "HydrogenPeroxide" in treated.chemicals.IDs:
            treated.imass["HydrogenPeroxide"] += h2o2_kgph

    def _design(self):
        feed = self.ins[0]

        m_kgph = float(feed.F_mass)
        vol_m3ph = m_kgph / self.slurry_density_kg_per_m3
        V = vol_m3ph * self.residence_time_hr

        dry_mass = getattr(self, "_dry_mass_kgph", 0.0)
        h2o2_kgph = getattr(self, "_h2o2_kgph", 0.0)
        h2o2_cost_usd_per_hr = h2o2_kgph * self.h2o2_price_usd_per_kg

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Temperature (K)"] = self.temperature_K
        self.design_results["Dry feed basis (kg/h)"] = dry_mass
        self.design_results["H2O2 addition (kg/h)"] = h2o2_kgph
        self.design_results["H2O2 cost ($/h)"] = h2o2_cost_usd_per_hr

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Peroxide pretreatment"] = capex

        h2o2_kgph = getattr(self, "_h2o2_kgph", 0.0)
        h2o2_cost_usd_per_hr = h2o2_kgph * self.h2o2_price_usd_per_kg

        if h2o2_cost_usd_per_hr > 0:
            self.add_OPEX = {
                "H2O2 reagent cost": h2o2_cost_usd_per_hr
            }

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

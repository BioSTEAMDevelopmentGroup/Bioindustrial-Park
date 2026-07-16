# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Milling / size reduction

Applies explicit mass loss during milling/shredding.
Sends lost material to a 'losses' stream (same composition as feed).

Economics:
- Electricity: kWh per dry ton of dry material (preferred) or legacy wet ton basis
- CAPEX: anchor scaling for hammer mill (purchase cost), then apply install factor
"""

import math
import biosteam as bst

__all__ = ('Mill',)

KG_PER_METRIC_TON = 1000.0
KG_PER_DRY_TON = 907.18474  # US short ton
HR_PER_DAY = 24.0


class Mill(bst.Unit):
    """
    Multi-train hammer mill. Constructed in three places in this codebase
    (`systems/_ad_biogas_system.py`, `systems/_vfa_ad_system.py`,
    `systems/_integrated_system.py`), all with the same
    `assumptions.yaml`-sourced parameters.

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

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
    _N_ins = 1
    _N_outs = 2  # milled_biomass, milling_losses

    def __init__(
        self, ID="", ins=None, outs=(),
        loss_frac=0.15,

        # utilities
        power_kWh_per_dry_ton_dry=None,  # preferred
        power_kWh_per_ton_wet=None,      # legacy optional

        # costing
        capex_model=None,               # e.g., "inl_hammermill_anchor"
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

        self.capex_model = capex_model
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

        installed_total = 0.0
        if (self.capex_model or "").lower() == "inl_hammermill_anchor":
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

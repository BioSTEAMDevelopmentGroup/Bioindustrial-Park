# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Membrane biogas upgrading to pipeline-quality biomethane.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('BiogasUpgrading',)


class BiogasUpgrading(bst.Unit):
    _N_ins = 1
    _N_outs = 2  # biomethane, offgas

    COST_ITEM = "Membrane upgrading skid"
    F_BM = {COST_ITEM: 1.0}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        ch4_recovery=0.98,
        co2_removal=0.95,
        electricity_kwh_per_Nm3_raw=0.25,
        capex_usd_per_Nm3ph_raw=2200.0,
        maintenance_frac_of_capex_per_yr=0.035,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.ch4_recovery = float(ch4_recovery)
        self.co2_removal = float(co2_removal)
        self.electricity_kwh_per_Nm3_raw = float(electricity_kwh_per_Nm3_raw)
        self.capex_usd_per_Nm3ph_raw = float(capex_usd_per_Nm3ph_raw)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        # make sure instance has the same F_BM mapping
        self.F_BM = dict(type(self).F_BM)

        if not (0.0 <= self.ch4_recovery <= 1.0):
            raise ValueError("ch4_recovery must be between 0 and 1.")
        if not (0.0 <= self.co2_removal <= 1.0):
            raise ValueError("co2_removal must be between 0 and 1.")
        if self.electricity_kwh_per_Nm3_raw < 0:
            raise ValueError("electricity_kwh_per_Nm3_raw must be >= 0.")
        if self.capex_usd_per_Nm3ph_raw < 0:
            raise ValueError("capex_usd_per_Nm3ph_raw must be >= 0.")
        if not (0.0 <= self.maintenance_frac_of_capex_per_yr <= 1.0):
            raise ValueError("maintenance_frac_of_capex_per_yr must be between 0 and 1.")

    def _run(self):
        raw = self.ins[0]
        biomethane, offgas = self.outs

        biomethane.empty()
        offgas.empty()

        biomethane.phase = "g"
        offgas.phase = "g"

        ch4_in = float(raw.imol["Methane"])
        co2_in = float(raw.imol["CarbonDioxide"])

        ch4_to_bm = self.ch4_recovery * ch4_in
        biomethane.imol["Methane"] = ch4_to_bm
        offgas.imol["Methane"] = ch4_in - ch4_to_bm

        co2_to_off = self.co2_removal * co2_in
        offgas.imol["CarbonDioxide"] = co2_to_off
        biomethane.imol["CarbonDioxide"] = co2_in - co2_to_off

        for cid in raw.chemicals.IDs:
            if cid in ("Methane", "CarbonDioxide"):
                continue
            n = float(raw.imol[cid])
            if n > 0:
                offgas.imol[cid] = n

    def _design(self):
        raw = self.ins[0]
        biomethane, offgas = self.outs

        # Use dry raw biogas basis for sizing and power
        dry_gas_IDs = ["Methane", "CarbonDioxide", "HydrogenSulfide", "Nitrogen", "Oxygen"]
        n_kmolph_dry = sum(float(raw.imol[i]) for i in dry_gas_IDs if i in raw.chemicals.IDs)

        Q_Nm3ph_dry = 22.414 * n_kmolph_dry
        self.design_results["Raw biogas flow (Nm3/h, dry)"] = Q_Nm3ph_dry

        kW = self.electricity_kwh_per_Nm3_raw * Q_Nm3ph_dry
        self.design_results["Electricity demand (kW)"] = kW

        # Explicit overwrite
        self.power_utility.consumption = kW
        self.power_utility.production = 0.0

        if biomethane.F_mol > 0:
            self.design_results["Biomethane CH4 mol%"] = (
                100.0 * float(biomethane.imol["Methane"]) / float(biomethane.F_mol)
            )
        else:
            self.design_results["Biomethane CH4 mol%"] = 0.0

        self.design_results["Methane slip (kmol/h)"] = float(offgas.imol["Methane"])

        ch4_in = float(raw.imol["Methane"])
        if ch4_in > 0:
            self.design_results["Methane recovery actual"] = (
                float(biomethane.imol["Methane"]) / ch4_in
            )
        else:
            self.design_results["Methane recovery actual"] = 0.0

    def _cost(self):
        Q_Nm3ph = self.design_results.get("Raw biogas flow (Nm3/h, dry)")
        if Q_Nm3ph is None:
            raw = self.ins[0]
            dry_gas_IDs = ["Methane", "CarbonDioxide", "HydrogenSulfide", "Nitrogen", "Oxygen"]
            n_kmolph_dry = sum(float(raw.imol[i]) for i in dry_gas_IDs if i in raw.chemicals.IDs)
            Q_Nm3ph = 22.414 * n_kmolph_dry
            self.design_results["Raw biogas flow (Nm3/h, dry)"] = Q_Nm3ph

        capex = self.capex_usd_per_Nm3ph_raw * float(Q_Nm3ph)
        self.baseline_purchase_costs[self.COST_ITEM] = capex

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Membrane upgrading maintenance": annual_maintenance / operating_hours_per_yr
            }

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Thermal (heating) pretreatment before methanogenic AD.
"""

from __future__ import annotations

import biosteam as bst

__all__ = ('HeatingPretreatment',)


class HeatingPretreatment(bst.Unit):
    _N_ins = 1
    _N_outs = 1  # heated_biomass

    def __init__(
        self, ID="", ins=None, outs=(),
        target_temperature_K=338.15,
        residence_time_hr=0.25,
        slurry_density_kg_per_m3=1000.0,
        cp_kJ_per_kgK=4.18,

        # costing
        capex_usd=0.0,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.target_temperature_K = float(target_temperature_K)
        self.residence_time_hr = float(residence_time_hr)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)
        self.cp_kJ_per_kgK = float(cp_kJ_per_kgK)

        self.capex_usd = float(capex_usd)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Heating pretreatment"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        heated = self.outs[0]
        heated.copy_like(feed)
        heated.phase = feed.phase

    def _design(self):
        feed = self.ins[0]

        m_kgph = float(feed.F_mass)
        vol_m3ph = m_kgph / self.slurry_density_kg_per_m3

        T_in = float(getattr(feed, "T", self.target_temperature_K))
        T_out = self.target_temperature_K
        dT = max(T_out - T_in, 0.0)

        V = vol_m3ph * self.residence_time_hr
        Q_kJph = m_kgph * self.cp_kJ_per_kgK * dT

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Inlet T (K)"] = T_in
        self.design_results["Target T (K)"] = T_out
        self.design_results["Heating duty (kJ/h)"] = Q_kJph

        if Q_kJph > 0:
            try:
                self.add_heat_utility(Q_kJph, T_in, T_out)
            except TypeError:
                self.add_heat_utility(Q_kJph, T_in)

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Heating pretreatment"] = capex

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Heating pretreatment maintenance": annual_maintenance / operating_hours_per_yr
            }

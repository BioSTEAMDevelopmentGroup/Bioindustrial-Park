# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Pressate membrane concentrator for the biostimulant side-stream.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('PressateConcentrator',)


class PressateConcentrator(bst.Unit):
    _N_ins = 1
    _N_outs = 2  # concentrated_biostimulant, permeate

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        retained_solute_IDs=(
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
        ),
        water_recovery_to_permeate=0.70,
        retained_solute_recovery_to_concentrate=0.95,
        design_flux_L_m2_h=35.0,
        operating_pressure_bar=5.0,
        electricity_kWh_per_m3_feed=0.8,
        capex_usd_per_m2=120.0,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.retained_solute_IDs = tuple(retained_solute_IDs)
        self.water_recovery_to_permeate = float(water_recovery_to_permeate)
        self.retained_solute_recovery_to_concentrate = float(
            retained_solute_recovery_to_concentrate
        )
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)
        self.operating_pressure_bar = float(operating_pressure_bar)
        self.electricity_kWh_per_m3_feed = float(electricity_kWh_per_m3_feed)
        self.capex_usd_per_m2 = float(capex_usd_per_m2)
        self.maintenance_frac_of_capex_per_yr = float(
            maintenance_frac_of_capex_per_yr
        )

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Pressate concentrator"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        conc, perm = self.outs

        conc.empty()
        perm.empty()

        conc.phase = "l"
        perm.phase = "l"

        ids = set(feed.chemicals.IDs)

        water = float(feed.imass["Water"]) if "Water" in ids else 0.0
        water_to_perm = self.water_recovery_to_permeate * water
        water_to_conc = water - water_to_perm

        if "Water" in ids:
            conc.imass["Water"] = water_to_conc
            perm.imass["Water"] = water_to_perm

        retained = set(self.retained_solute_IDs)

        for cid in feed.chemicals.IDs:
            if cid == "Water":
                continue

            m = float(feed.imass[cid])
            if m <= 0:
                continue

            if cid in retained:
                m_to_conc = self.retained_solute_recovery_to_concentrate * m
                m_to_perm = m - m_to_conc
            else:
                # base case: non-target solubles mostly follow permeate/water
                m_to_perm = 0.70 * m
                m_to_conc = m - m_to_perm

            conc.imass[cid] = m_to_conc
            perm.imass[cid] = m_to_perm

        self._feed_m3ph = float(feed.F_vol)
        self._conc_m3ph = float(conc.F_vol)
        self._perm_m3ph = float(perm.F_vol)

    def _design(self):
        feed = self.ins[0]

        feed_m3ph = float(feed.F_vol)
        flux_m3_m2_h = self.design_flux_L_m2_h / 1000.0  # LMH -> m3/m2/h
        membrane_area_m2 = feed_m3ph / flux_m3_m2_h if flux_m3_m2_h > 0 else 0.0

        power_kW = self.electricity_kWh_per_m3_feed * feed_m3ph

        self.design_results["Feed flow (m3/h)"] = feed_m3ph
        self.design_results["Design flux (L/m2-h)"] = self.design_flux_L_m2_h
        self.design_results["Operating pressure (bar)"] = self.operating_pressure_bar
        self.design_results["Membrane area (m2)"] = membrane_area_m2
        self.design_results["Electricity demand (kW)"] = power_kW

        self.power_utility(power_kW)

    def _cost(self):
        membrane_area_m2 = self.design_results.get("Membrane area (m2)", 0.0)
        capex = membrane_area_m2 * self.capex_usd_per_m2

        self.baseline_purchase_costs["Pressate concentrator"] = capex

        # Annual maintenance -> add_OPEX so it flows through
        # SABREBaselineTEA.VOC -> solve_product_msp correctly
        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Pressate concentrator maintenance": annual_maintenance / operating_hours_per_yr
            }

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Finishing evaporator for the pressate-derived biostimulant product.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('BiostimulantEvaporator',)


class BiostimulantEvaporator(bst.Unit):
    """
    Finishing evaporator for the pressate-derived biostimulant product.

    Not constructed anywhere in this codebase currently: both real call
    sites (`systems/_ad_biogas_system.py`, `systems/_vfa_ad_system.py`)
    gate construction behind `evA.get("enabled", False)`, and
    `assumptions.yaml` has no `pressate_biostimulant.evaporator` section
    at all (confirmed by grep across the whole yaml file) -- so `evA` is
    always `{}` and this unit is never actually instantiated. Same
    situation as `DigestateDecanterCentrifuge`.

    Sources
    -------
    All parameters below (target_solids_wt_frac, boiling_temperature_K,
    latent_heat_kJ_per_kg, sensible_cp_kJ_per_kgK,
    electricity_kWh_per_kg_water_evap, nonwater_recovery_to_product,
    capex_ref_usd, ref_evaporation_kgph, scale_exponent,
    maintenance_frac_of_capex_per_yr):
        No assumptions.yaml section exists for this unit at all -- unlike
        every other unwired-but-yaml-documented unit in this conversion
        (e.g. `DigestateDecanterCentrifuge`), there isn't even an orphaned
        yaml block to check these against. Every default here is a pure
        code-level scoping assumption with no literature or yaml grounding
        whatsoever. This module's own pre-existing docstring had no
        citations either (just "mild vacuum evaporation ~60 C" as an
        inline comment on `boiling_temperature_K`).

    Maintenance OPEX wiring: `_cost()` computed annual maintenance into
    `design_results` but never assigned it to `self.add_OPEX`, unlike
    every sibling unit with capex-dependent maintenance in this codebase
    (`PressateConcentrator`, `BiogasUpgrading`) -- meaning it would
    silently not affect TEA/MSP if this unit were ever wired in. Fixed
    during this conversion (zero regression risk, since the unit carries
    no test coverage either way -- it's never simulated) to match the
    add_OPEX pattern used everywhere else in this file family, rather than
    left as a latent gap for whoever eventually wires this unit in.
    """

    _N_ins = 1
    _N_outs = 2  # concentrated_product, evaporated_water

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        target_solids_wt_frac=0.20,
        boiling_temperature_K=333.15,  # mild vacuum evaporation ~60 C
        latent_heat_kJ_per_kg=2350.0,
        sensible_cp_kJ_per_kgK=4.18,
        electricity_kWh_per_kg_water_evap=0.0,
        nonwater_recovery_to_product=0.995,
        capex_ref_usd=750000.0,
        ref_evaporation_kgph=10000.0,
        scale_exponent=0.60,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.target_solids_wt_frac = float(target_solids_wt_frac)
        self.boiling_temperature_K = float(boiling_temperature_K)
        self.latent_heat_kJ_per_kg = float(latent_heat_kJ_per_kg)
        self.sensible_cp_kJ_per_kgK = float(sensible_cp_kJ_per_kgK)
        self.electricity_kWh_per_kg_water_evap = float(electricity_kWh_per_kg_water_evap)
        self.nonwater_recovery_to_product = float(nonwater_recovery_to_product)

        self.capex_ref_usd = float(capex_ref_usd)
        self.ref_evaporation_kgph = float(ref_evaporation_kgph)
        self.scale_exponent = float(scale_exponent)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Biostimulant evaporator"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        product, vapor = self.outs

        product.empty()
        vapor.empty()

        product.phase = "l"
        vapor.phase = "g"

        ids = set(feed.chemicals.IDs)

        water_in = float(feed.imass["Water"]) if "Water" in ids else 0.0
        nonwater_in = float(feed.F_mass) - water_in

        # Keep nearly all non-water in product
        rec = min(max(self.nonwater_recovery_to_product, 0.0), 1.0)

        for cid in feed.chemicals.IDs:
            if cid == "Water":
                continue
            m = float(feed.imass[cid])
            if m <= 0:
                continue
            m_prod = rec * m
            m_vap = m - m_prod
            product.imass[cid] = m_prod
            vapor.imass[cid] = m_vap

        nonwater_prod = sum(
            float(product.imass[cid]) for cid in feed.chemicals.IDs if cid != "Water"
        )

        x_target = self.target_solids_wt_frac
        if not (0.0 < x_target < 1.0):
            raise ValueError("target_solids_wt_frac must be between 0 and 1")

        # water needed in final liquid to hit target solids
        if nonwater_prod <= 0:
            water_to_product = 0.0
        else:
            water_to_product = nonwater_prod * (1.0 - x_target) / x_target

        water_to_product = min(max(water_to_product, 0.0), water_in)
        water_to_vapor = water_in - water_to_product

        if "Water" in ids:
            product.imass["Water"] = water_to_product
            vapor.imass["Water"] = water_to_vapor

        self._water_evaporated_kgph = water_to_vapor
        self._product_solids_wt_frac = (
            nonwater_prod / product.F_mass if product.F_mass > 0 else 0.0
        )

    def _design(self):
        feed = self.ins[0]

        water_evap = getattr(self, "_water_evaporated_kgph", 0.0)

        sensible = 0.0
        feed_T = float(feed.T)
        if water_evap > 0 and feed_T < self.boiling_temperature_K:
            sensible = float(feed.F_mass) * self.sensible_cp_kJ_per_kgK * (
                self.boiling_temperature_K - feed_T
            )

        latent = water_evap * self.latent_heat_kJ_per_kg
        total_duty_kJph = sensible + latent
        power_kW = water_evap * self.electricity_kWh_per_kg_water_evap

        self.design_results["Feed flow (kg/h)"] = float(feed.F_mass)
        self.design_results["Water evaporated (kg/h)"] = water_evap
        self.design_results["Target solids wt frac"] = self.target_solids_wt_frac
        self.design_results["Actual product solids wt frac"] = getattr(
            self, "_product_solids_wt_frac", 0.0
        )
        self.design_results["Boiling temperature (K)"] = self.boiling_temperature_K
        self.design_results["Sensible duty (kJ/h)"] = sensible
        self.design_results["Latent duty (kJ/h)"] = latent
        self.design_results["Total duty (kJ/h)"] = total_duty_kJph
        self.design_results["Electricity demand (kW)"] = power_kW

        if total_duty_kJph > 0:
            self.add_heat_utility(total_duty_kJph, T_in=feed.T, T_out=self.boiling_temperature_K)
        if power_kW > 0:
            self.power_utility(power_kW)

    def _cost(self):
        evap = max(float(self.design_results.get("Water evaporated (kg/h)", 0.0)), 0.0)
        if self.ref_evaporation_kgph <= 0:
            raise ValueError("ref_evaporation_kgph must be > 0")

        if evap <= 0:
            capex = 0.0
        else:
            capex = self.capex_ref_usd * (evap / self.ref_evaporation_kgph) ** self.scale_exponent

        self.baseline_purchase_costs["Biostimulant evaporator"] = capex

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Biostimulant evaporator maintenance": annual_maintenance / operating_hours_per_yr
            }

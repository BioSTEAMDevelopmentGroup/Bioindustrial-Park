# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst

from biorefineries.sabre.utils import OPERATING_HOURS_PER_YEAR, load_assumptions

__all__ = ('PressateConcentrator', 'BiostimulantEvaporator')

# data/biostimulant.yaml, loaded once
_BIOSTIMULANT = load_assumptions("biostimulant.yaml")
_CONCENTRATOR = _BIOSTIMULANT["pressate_concentrator"]
_EVAPORATOR = _BIOSTIMULANT["biostimulant_evaporator"]


class PressateConcentrator(bst.Unit):
    """
    Membrane concentrator for the pressate-derived biostimulant side-stream.

    Parameters
    ----------
    ins : stream
        Pressate feed (e.g. `Press`'s pressate outlet).
    outs : tuple[stream, stream]
        Concentrate (retained solutes) and permeate.
    enabled : bool
        If False, this unit is a no-op bypass: the feed passes straight
        through to the concentrate outlet unchanged (permeate stays
        empty), and `_design`/`_cost` report zero membrane area,
        electricity, installed cost, and maintenance OPEX. Lets this
        unit sit permanently wired into a flowsheet's `path` while being
        switched on/off without restructuring stream connectivity.
    retained_solute_IDs : Iterable[str]
        Chemical IDs retained (at `retained_solute_recovery_to_concentrate`)
        to the concentrate; everything else (except Water) is split by
        `nontarget_solute_recovery_to_permeate`.
    water_recovery_to_permeate : float
        Fraction of feed water routed to the permeate (0-1).
    retained_solute_recovery_to_concentrate : float
        Fraction of `retained_solute_IDs` mass captured to the
        concentrate (0-1).
    nontarget_solute_recovery_to_permeate : float
        Fraction of non-`retained_solute_IDs`, non-Water mass (0-1)
        routed to the permeate; the remainder stays in the concentrate.
    design_flux_L_m2_h : float
        Design membrane flux, used to size membrane area.
    operating_pressure_bar : float
        Operating pressure (design-basis only; reported in
        `design_results`, not used in a cost/duty calculation).
    electricity_kWh_per_m3_feed : float
        Electricity intensity per m3 of feed processed.
    capex_usd_per_m2 : float
        Installed cost per m2 of membrane area.
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of installed cost, added to
        `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/biostimulant.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2  # concentrated_biostimulant, permeate
    _F_BM_default = {"Pressate concentrator": _CONCENTRATOR["F_BM"]}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        enabled=_CONCENTRATOR["enabled"],
        retained_solute_IDs=tuple(_CONCENTRATOR["retained_solute_IDs"]),
        water_recovery_to_permeate=_CONCENTRATOR["water_recovery_to_permeate"],
        retained_solute_recovery_to_concentrate=_CONCENTRATOR["retained_solute_recovery_to_concentrate"],
        nontarget_solute_recovery_to_permeate=_CONCENTRATOR["nontarget_solute_recovery_to_permeate"],
        design_flux_L_m2_h=_CONCENTRATOR["design_flux_L_m2_h"],
        operating_pressure_bar=_CONCENTRATOR["operating_pressure_bar"],
        electricity_kWh_per_m3_feed=_CONCENTRATOR["electricity_kWh_per_m3_feed"],
        capex_usd_per_m2=_CONCENTRATOR["capex_usd_per_m2"],
        maintenance_frac_of_capex_per_yr=_CONCENTRATOR["maintenance_frac_of_capex_per_yr"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.enabled = bool(enabled)
        self.retained_solute_IDs = tuple(retained_solute_IDs)
        self.water_recovery_to_permeate = float(water_recovery_to_permeate)
        self.retained_solute_recovery_to_concentrate = float(
            retained_solute_recovery_to_concentrate
        )
        self.nontarget_solute_recovery_to_permeate = float(
            nontarget_solute_recovery_to_permeate
        )
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)
        self.operating_pressure_bar = float(operating_pressure_bar)
        self.electricity_kWh_per_m3_feed = float(electricity_kWh_per_m3_feed)
        self.capex_usd_per_m2 = float(capex_usd_per_m2)
        self.maintenance_frac_of_capex_per_yr = float(
            maintenance_frac_of_capex_per_yr
        )

    def _run(self):
        feed = self.ins[0]
        conc, perm = self.outs

        conc.empty()
        perm.empty()

        conc.phase = "l"
        perm.phase = "l"

        if not self.enabled:
            conc.copy_like(feed)
            self._feed_m3ph = float(feed.F_vol)
            self._conc_m3ph = float(conc.F_vol)
            self._perm_m3ph = 0.0
            return

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
                m_to_perm = self.nontarget_solute_recovery_to_permeate * m
                m_to_conc = m - m_to_perm

            conc.imass[cid] = m_to_conc
            perm.imass[cid] = m_to_perm

        self._feed_m3ph = float(feed.F_vol)
        self._conc_m3ph = float(conc.F_vol)
        self._perm_m3ph = float(perm.F_vol)

    def _design(self):
        feed = self.ins[0]
        feed_m3ph = float(feed.F_vol)

        if not self.enabled:
            self.design_results["Feed flow (m3/h)"] = feed_m3ph
            self.design_results["Membrane area (m2)"] = 0.0
            self.design_results["Electricity demand (kW)"] = 0.0
            self.power_utility(0.0)
            return

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
        # SaBReTEA.VOC -> solve_product_msp correctly
        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            self.add_OPEX = {
                "Pressate concentrator maintenance": annual_maintenance / OPERATING_HOURS_PER_YEAR
            }


class BiostimulantEvaporator(bst.Unit):
    """
    Finishing step for the pressate-derived biostimulant product: adjusts
    a concentrate stream to a target solids content, in either direction.

    Direction logic (in `_run`):
    - If the concentrate is below target_solids_wt_frac (needs
      concentrating): evaporates the excess water, same as this unit's
      original press-pretreatment-side design -- utility duty and CAPEX
      scale with the amount of water evaporated (`_design`/`_cost`,
      unchanged from before).
    - If the concentrate is above target_solids_wt_frac (needs diluting):
      adds water to hit the target. Water is drawn from ins[1] (permeate)
      first, up to its availability; any remaining shortfall is drawn
      from ins[2] (fresh water). No evaporation duty or CAPEX applies to
      this direction -- diluting is just mixing, not a heat/mass-transfer
      operation, so `_design`/`_cost` naturally come out to zero when
      `_water_evaporated_kgph` is zero.

    Parameters
    ----------
    ins : tuple[stream, stream, stream]
        Concentrate to adjust, permeate (water source when diluting), and
        fresh water make-up (settable/phantom stream, filled in by
        `_run()` with exactly the shortfall permeate doesn't cover).
    outs : tuple[stream, stream, stream]
        Biostimulant product (adjusted to `target_solids_wt_frac`), vapor
        (water removed by evaporation; zero when diluting), and residual
        permeate (permeate minus whatever water it contributed).
    enabled : bool
        If False, this unit is a no-op bypass: the concentrate passes
        straight through to the product outlet unchanged (vapor stays
        empty, permeate passes straight through to residual permeate
        unchanged, no fresh water is drawn), and `_design`/`_cost` report
        zero duty, electricity, installed cost, and maintenance OPEX.
        Lets this unit sit permanently wired into a flowsheet's `path`
        while being switched on/off without restructuring stream
        connectivity.
    target_solids_wt_frac : float
        Target dry-solids weight fraction of the product (0-1); the unit
        evaporates or dilutes to reach it, whichever direction applies.
    boiling_temperature_K : float
        Evaporation temperature; also the reference temperature for the
        latent-heat lookup (see Notes below).
    electricity_kWh_per_kg_water_evap : float
        Electricity intensity per kg of water evaporated.
    nonwater_recovery_to_product : float
        Fraction of non-water solute mass recovered to the product (0-1);
        the remainder is carried out with the vapor.
    capex_ref_usd : float
        Installed cost at `ref_evaporation_kgph`, scaled by
        `scale_exponent` for other evaporation rates.
    ref_evaporation_kgph : float
        Reference water-evaporation rate (kg/h) for CAPEX scaling.
    scale_exponent : float
        Power-law scaling exponent for installed CAPEX vs. evaporation
        rate.
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of installed cost, added to
        `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/biostimulant.yaml for the default values and references.
    """

    _N_ins = 3
    _N_outs = 3  # biostimulant_product, vapor, residual_permeate
    _F_BM_default = {"Biostimulant evaporator": _EVAPORATOR["F_BM"]}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        enabled=_EVAPORATOR["enabled"],
        target_solids_wt_frac=_EVAPORATOR["target_solids_wt_frac"],
        boiling_temperature_K=_EVAPORATOR["boiling_temperature_K"],
        electricity_kWh_per_kg_water_evap=_EVAPORATOR["electricity_kWh_per_kg_water_evap"],
        nonwater_recovery_to_product=_EVAPORATOR["nonwater_recovery_to_product"],
        capex_ref_usd=_EVAPORATOR["capex_ref_usd"],
        ref_evaporation_kgph=_EVAPORATOR["ref_evaporation_kgph"],
        scale_exponent=_EVAPORATOR["scale_exponent"],
        maintenance_frac_of_capex_per_yr=_EVAPORATOR["maintenance_frac_of_capex_per_yr"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.enabled = bool(enabled)
        self.target_solids_wt_frac = float(target_solids_wt_frac)
        self.boiling_temperature_K = float(boiling_temperature_K)
        self.electricity_kWh_per_kg_water_evap = float(electricity_kWh_per_kg_water_evap)
        self.nonwater_recovery_to_product = float(nonwater_recovery_to_product)

        self.capex_ref_usd = float(capex_ref_usd)
        self.ref_evaporation_kgph = float(ref_evaporation_kgph)
        self.scale_exponent = float(scale_exponent)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

    def _run(self):
        conc_in, permeate_in, fresh_water_in = self.ins
        product, vapor, residual_permeate = self.outs

        product.empty()
        vapor.empty()
        fresh_water_in.empty()
        residual_permeate.copy_like(permeate_in)

        product.phase = "l"
        vapor.phase = "g"
        residual_permeate.phase = "l"
        fresh_water_in.phase = "l"

        if not self.enabled:
            product.copy_like(conc_in)
            self._water_evaporated_kgph = 0.0
            self._permeate_water_used_kgph = 0.0
            self._fresh_water_used_kgph = 0.0
            water = float(conc_in.imass["Water"]) if "Water" in conc_in.chemicals.IDs else 0.0
            nonwater = product.F_mass - water
            self._product_solids_wt_frac = nonwater / product.F_mass if product.F_mass > 0 else 0.0
            return

        ids = set(conc_in.chemicals.IDs)

        water_in = float(conc_in.imass["Water"]) if "Water" in ids else 0.0

        # Keep nearly all non-water in product
        rec = min(max(self.nonwater_recovery_to_product, 0.0), 1.0)

        for cid in conc_in.chemicals.IDs:
            if cid == "Water":
                continue
            m = float(conc_in.imass[cid])
            if m <= 0:
                continue
            m_prod = rec * m
            m_vap = m - m_prod
            product.imass[cid] = m_prod
            vapor.imass[cid] = m_vap

        nonwater_prod = sum(
            float(product.imass[cid]) for cid in conc_in.chemicals.IDs if cid != "Water"
        )

        x_target = self.target_solids_wt_frac
        if not (0.0 < x_target < 1.0):
            raise ValueError("target_solids_wt_frac must be between 0 and 1")

        # water needed in final liquid to hit target solids
        if nonwater_prod <= 0:
            water_to_product = 0.0
        else:
            water_to_product = nonwater_prod * (1.0 - x_target) / x_target
        water_to_product = max(water_to_product, 0.0)

        permeate_water_used = 0.0
        fresh_water_used = 0.0

        if water_to_product <= water_in:
            # Concentrating: evaporate the excess water
            water_to_vapor = water_in - water_to_product
        else:
            # Diluting: draw the shortfall from permeate first, then fresh water
            water_to_vapor = 0.0
            shortfall = water_to_product - water_in

            permeate_water_avail = (
                float(permeate_in.imass["Water"])
                if "Water" in permeate_in.chemicals.IDs else 0.0
            )
            permeate_water_used = min(shortfall, permeate_water_avail)
            shortfall -= permeate_water_used
            fresh_water_used = max(shortfall, 0.0)

        product.imass["Water"] = water_to_product
        vapor.imass["Water"] = water_to_vapor

        if permeate_water_used > 0:
            residual_permeate.imass["Water"] = max(
                float(permeate_in.imass["Water"]) - permeate_water_used, 0.0
            )
        if fresh_water_used > 0:
            fresh_water_in.imass["Water"] = fresh_water_used

        self._water_evaporated_kgph = water_to_vapor
        self._permeate_water_used_kgph = permeate_water_used
        self._fresh_water_used_kgph = fresh_water_used
        self._product_solids_wt_frac = (
            nonwater_prod / product.F_mass if product.F_mass > 0 else 0.0
        )

    def _design(self):
        conc_in = self.ins[0]

        if not self.enabled:
            self.design_results["Feed flow (kg/h)"] = float(conc_in.F_mass)
            self.design_results["Water evaporated (kg/h)"] = 0.0
            self.design_results["Total duty (kJ/h)"] = 0.0
            self.design_results["Electricity demand (kW)"] = 0.0
            self.power_utility(0.0)
            return

        water_evap = getattr(self, "_water_evaporated_kgph", 0.0)

        sensible = 0.0
        feed_T = float(conc_in.T)
        if water_evap > 0 and feed_T < self.boiling_temperature_K:
            sensible = float(conc_in.F_mass) * conc_in.Cp * (
                self.boiling_temperature_K - feed_T
            )

        water_chem = self.chemicals.Water
        latent_heat_kJ_per_kg = water_chem.Hvap(self.boiling_temperature_K) / water_chem.MW
        latent = water_evap * latent_heat_kJ_per_kg
        total_duty_kJph = sensible + latent
        power_kW = water_evap * self.electricity_kWh_per_kg_water_evap

        self.design_results["Feed flow (kg/h)"] = float(conc_in.F_mass)
        self.design_results["Water evaporated (kg/h)"] = water_evap
        self.design_results["Permeate water used (kg/h)"] = getattr(
            self, "_permeate_water_used_kgph", 0.0
        )
        self.design_results["Fresh water used (kg/h)"] = getattr(
            self, "_fresh_water_used_kgph", 0.0
        )
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
            self.add_heat_utility(total_duty_kJph, T_in=conc_in.T, T_out=self.boiling_temperature_K)
        
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
            self.add_OPEX = {
                "Biostimulant evaporator maintenance": annual_maintenance / OPERATING_HOURS_PER_YEAR
            }

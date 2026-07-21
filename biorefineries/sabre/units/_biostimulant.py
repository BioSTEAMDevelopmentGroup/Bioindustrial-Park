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
    Finishing step for the pressate-derived biostimulant product: evaporates
    water to reach a target solids content.

    Parameters
    ----------
    ins : stream
        Concentrate to evaporate.
    outs : tuple[stream, stream]
        Biostimulant product (adjusted to `target_solids_wt_frac`) and
        condensed vapor (water removed by evaporation, plus the small
        non-recovered-solute fraction, cooled to `condensing_temperature_K`
        by the auxiliary condenser -- see Notes -- so it leaves this unit
        as a liquid, e.g. ready for a wastewater mixer downstream).
    target_solids_wt_frac : float or None
        Target dry-solids weight fraction of the product (0-1). `None` is
        a no-op bypass: the concentrate passes straight through to the
        product outlet unchanged.
    boiling_temperature_K : float
        Evaporation temperature; also the reference temperature for the
        latent-heat lookup (see Notes below).
    condensing_temperature_K : float
        Temperature the auxiliary condenser cools the vapor down to.
    nonwater_recovery_to_product : float
        Fraction of non-water solute mass recovered to the product (0-1);
        the remainder is carried out with the vapor.
    tau : float
        Nominal vessel holdup time, in hr. No literature/prior-yaml grounding.
    vacuum_system_preference : str
        Forwarded to `bst.VacuumSystem`. Defaults to 'Liquid-ring pump'.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    Notes
    -----
    No literature or prior-yaml cost basis exists for a flat or
    evaporation-rate-scaled evaporator vessel cost, so this unit has no
    cost item of its own -- installed cost comes entirely from the
    `evaporation_hx`, `condenser`, and `vacuum_system` auxiliaries' own
    BioSTEAM-computed design and costing.

    `boiling_temperature_K` (333.15 K / 60 C) is below water's atmospheric
    boiling point, i.e. this represents a vacuum evaporator. `vacuum_system`
    (a `bst.VacuumSystem`, sized from the raw vapor flow and the water
    saturation pressure at `boiling_temperature_K`) costs the vacuum
    pump/ejector itself; tried 'Steam-jet ejector' at this unit's scale and
    got a steam duty larger than the evaporator's own heating duty --
    implausible, looks like the correlation extrapolating outside its valid
    range (same failure mode seen sizing HeatingPretreatment's auxiliary HX
    at combined_PTE scale) -- so 'Liquid-ring pump' is the default instead.

    See Also
    --------
    Refer to data/biostimulant.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2  # biostimulant_product, condensed_vapor
    auxiliary_unit_names = ('evaporation_hx', 'condenser', 'vacuum_system')

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        target_solids_wt_frac=_EVAPORATOR["target_solids_wt_frac"],
        boiling_temperature_K=_EVAPORATOR["boiling_temperature_K"],
        condensing_temperature_K=_EVAPORATOR["condensing_temperature_K"],
        nonwater_recovery_to_product=_EVAPORATOR["nonwater_recovery_to_product"],
        tau=_EVAPORATOR["tau"],
        vacuum_system_preference=_EVAPORATOR["vacuum_system_preference"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.target_solids_wt_frac = (
            None if target_solids_wt_frac is None else float(target_solids_wt_frac)
        )
        self.boiling_temperature_K = float(boiling_temperature_K)
        self.condensing_temperature_K = float(condensing_temperature_K)
        self.nonwater_recovery_to_product = float(nonwater_recovery_to_product)
        self.tau = float(tau)  # hr
        self.vacuum_system_preference = vacuum_system_preference

        # Auxiliary evaporation_hx: heats the incoming concentrate up to the
        # combined (product + raw_vapor) exit state -- not a simple T-target,
        # since this is a partial, recovery-fraction-driven evaporation (not
        # an equilibrium flash), so ins[0]/outs[0] are hand-populated in _run
        # to match the real mass balance exactly, and its own _run is never
        # called (see .simulate(run=False) in _design).
        self.auxiliary(
            "evaporation_hx", bst.HXutility,
            ins=[None], outs=[None],
            T=self.boiling_temperature_K,
        )

        # Auxiliary condenser: cools/condenses the raw vapor (a pure-Water
        # sub-stream, set up in _run) down to condensing_temperature_K
        # before it leaves as outs[1]. Its ins[0] is an internal auxlet
        # stream (populated each _run); its outs[0] writes directly into
        # this unit's own real outs[1], so no extra system-level unit or
        # stream is needed downstream to do the condensing.
        self.auxiliary(
            "condenser", bst.HXutility,
            ins=[None], outs=[self.outs[1]],
            T=self.condensing_temperature_K, V=0,
        )

    def _run(self):
        conc_in, = self.ins
        product, condensate = self.outs
        evap_hx = self.evaporation_hx
        condenser = self.condenser
        raw_vapor = condenser.ins[0]
        heated_in, heated_out = evap_hx.ins[0], evap_hx.outs[0]

        product.empty()
        product.phase = "l"
        raw_vapor.empty()
        raw_vapor.phase = "g"

        if self.target_solids_wt_frac is None:
            product.copy_like(conc_in)
            condenser._run()
            heated_in.empty()
            heated_out.empty()
            self._water_evaporated_kgph = 0.0
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
            raw_vapor.imass[cid] = m_vap

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

        # Evaporate the excess water.
        water_to_vapor = max(water_in - water_to_product, 0.0)
        product.imass["Water"] = water_in - water_to_vapor
        raw_vapor.imass["Water"] = water_to_vapor

        # product and raw_vapor both actually leave the evaporator at its
        # operating (boiling) temperature.
        product.T = raw_vapor.T = self.boiling_temperature_K

        # evaporation_hx: ins[0] = incoming concentrate; outs[0] = the combined
        # (product + raw_vapor) exit state, hand-built to match the real mass
        # balance above exactly (a partial, recovery-fraction-driven
        # evaporation, not an equilibrium flash that HXutility's own T/V/H
        # targeting could reproduce).
        heated_in.copy_like(conc_in)
        heated_out.copy_like(conc_in)
        heated_out.phases = ("g", "l")
        heated_out["g"].copy_like(raw_vapor)
        heated_out["l"].copy_like(product)

        # Condense the raw vapor down to condensing_temperature_K; writes
        # directly into `condensate` (this unit's real outs[1]).
        condenser._run()

        self._water_evaporated_kgph = water_to_vapor
        self._product_solids_wt_frac = (
            nonwater_prod / product.F_mass if product.F_mass > 0 else 0.0
        )

    def _design(self):
        conc_in = self.ins[0]

        if self.target_solids_wt_frac is None:
            self.design_results["Feed flow (kg/h)"] = float(conc_in.F_mass)
            self.design_results["Water evaporated (kg/h)"] = 0.0
            self.design_results["Total duty (kJ/h)"] = 0.0
            self.evaporation_hx.simulate(run=False)
            self.condenser.simulate(run=False)
            # No vacuum_system here: unlike the HXutility auxiliaries (which
            # correctly go to exactly zero cost/utility on empty streams),
            # bst.VacuumSystem has a nonzero baseline floor (a small minimum
            # pump size) even with all-zero inputs -- so for this true no-op
            # bypass, skip creating it entirely rather than report a
            # nonzero cost/utility for a vacuum system that doesn't exist.
            return

        water_evap = getattr(self, "_water_evaporated_kgph", 0.0)

        self.design_results["Feed flow (kg/h)"] = float(conc_in.F_mass)
        self.design_results["Water evaporated (kg/h)"] = water_evap
        self.design_results["Target solids wt frac"] = self.target_solids_wt_frac
        self.design_results["Actual product solids wt frac"] = getattr(
            self, "_product_solids_wt_frac", 0.0
        )
        self.design_results["Boiling temperature (K)"] = self.boiling_temperature_K

        # Auxiliary evaporation_hx and condenser: their own duty/area/cost,
        # from BioSTEAM's own HXutility design+costing algorithms (run=False
        # since _run already populated their ins/outs above). No cost item of
        # our own is needed -- BioSTEAM rolls up all three auxiliaries'
        # installed costs automatically via auxiliary_unit_names.
        self.evaporation_hx.simulate(run=False)
        self.condenser.simulate(run=False)
        self.design_results["Total duty (kJ/h)"] = self.evaporation_hx.net_duty

        # Auxiliary vacuum_system: costs the vacuum pump/ejector itself, since
        # boiling_temperature_K is below water's atmospheric boiling point
        # (P_suction = water's saturation pressure there). vessel_volume is a
        # nominal holdup-time-based estimate (see `tau`), not a real vessel
        # design -- BioSTEAM's own air-inleakage sizing scales with it.
        raw_vapor = self.condenser.ins[0]
        vessel_volume = conc_in.F_vol * self.tau
        Psat = self.chemicals.Water.Psat(self.boiling_temperature_K)
        self.vacuum_system = bst.VacuumSystem(
            None, self.vacuum_system_preference,
            F_mass=raw_vapor.F_mass, F_vol=raw_vapor.F_vol,
            P_suction=Psat, vessel_volume=vessel_volume,
        )

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst

from biorefineries.sabre.utils import OPERATING_HOURS_PER_YEAR, load_assumptions

__all__ = ('EnzymaticPretreatment', 'HeatingPretreatment', 'PeroxidePretreatment')

# Load assumptions
_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]
_ENZYMATIC = _PRETREATMENT_AD["enzymatic"]["enzymatic"]
_PEROXIDE = _PRETREATMENT_AD["peroxide"]["peroxide"]
_HEATING = _PRETREATMENT_AD["combined_PTE"]["heating"]


class EnzymaticPretreatment(bst.Unit):
    """
    Enzymatic pretreatment reactor, one of the optional pre-AD pretreatment
    stages. Used when the selected pretreatment case has `kind: enzymatic`
    or `kind: combined_PE`/`combined_PTE` (which include an enzymatic
    stage in sequence).

    Parameters
    ----------
    ins : stream
        Biomass feed to be enzymatically treated.
    outs : stream
        Enzyme-treated biomass.
    temperature_K : float
        Reactor operating temperature. Feed is assumed to enter at its own
        stream temperature; the duty to raise it to `temperature_K` is
        charged as a heat utility.
    residence_time_hr : float
        Reactor residence time, used to size reactor volume.
    enzyme_dose_kg_per_kg_dry_feed : float
        Enzyme dose per kg of dry (non-water) feed processed.
    treated_fraction : float
        Fraction of dry feed mass actually treated (0-1); only this
        fraction's dry mass is used as the enzyme-dose basis.
    enzyme_recycle_factor : float
        Factor (>=1) by which enzyme is recycled; fresh enzyme addition is
        divided by this factor.
    capex_usd : float
        Flat installed capital cost (not a function of throughput).
    enzyme_price_usd_per_kg : float
        Enzyme reagent unit price, added to `add_OPEX`.
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of `capex_usd`, added to
        `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/pretreatment.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 1  # enzyme_treated_biomass
    _F_BM_default = {"Enzymatic pretreatment": _ENZYMATIC["F_BM"]}
    auxiliary_unit_names = ("heat_exchanger",)

    def __init__(
        self, ID="", ins=None, outs=(),
        temperature_K=_ENZYMATIC["temperature_K"],
        residence_time_hr=_ENZYMATIC["residence_time_hr"],
        enzyme_dose_kg_per_kg_dry_feed=_ENZYMATIC["enzyme_dose_kg_per_kg_dry_feed"],
        treated_fraction=_ENZYMATIC["treated_fraction"],
        enzyme_recycle_factor=_ENZYMATIC["enzyme_recycle_factor"],

        # economics
        capex_usd=_ENZYMATIC["capex_usd"],
        enzyme_price_usd_per_kg=_ENZYMATIC["enzyme_price_usd_per_kg"],
        maintenance_frac_of_capex_per_yr=_ENZYMATIC["maintenance_frac_of_capex_per_yr"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.temperature_K = float(temperature_K)
        self.residence_time_hr = float(residence_time_hr)
        self.enzyme_dose_kg_per_kg_dry_feed = float(enzyme_dose_kg_per_kg_dry_feed)
        self.treated_fraction = float(treated_fraction)
        self.enzyme_recycle_factor = float(enzyme_recycle_factor)

        self.capex_usd = float(capex_usd)
        self.enzyme_price_usd_per_kg = float(enzyme_price_usd_per_kg)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

    def _run(self):
        feed = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(feed)
        treated.phase = feed.phase
        treated.T = self.temperature_K

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
        vol_m3ph = feed.F_vol
        V = vol_m3ph * self.residence_time_hr

        dry_mass = getattr(self, "_dry_mass_kgph", 0.0)
        treated_dry_mass = getattr(self, "_treated_dry_mass_kgph", 0.0)
        enzyme_kgph = getattr(self, "_enzyme_kgph", 0.0)
        enzyme_cost_usd_per_hr = enzyme_kgph * self.enzyme_price_usd_per_kg
        self._enzyme_cost_usd_per_hr = enzyme_cost_usd_per_hr

        T_in = float(feed.T)
        T_out = self.temperature_K
        hx = self.auxiliary("heat_exchanger", bst.HXutility, ins=self.ins[0], T=T_out)
        hx.simulate()

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Inlet T (K)"] = T_in
        self.design_results["Temperature (K)"] = self.temperature_K
        self.design_results["Heating duty (kJ/h)"] = hx.net_duty
        self.design_results["Dry feed basis (kg/h)"] = dry_mass
        self.design_results["Treated dry feed basis (kg/h)"] = treated_dry_mass
        self.design_results["Treated fraction"] = self.treated_fraction
        self.design_results["Enzyme recycle factor"] = self.enzyme_recycle_factor
        self.design_results["Enzyme addition (kg/h)"] = enzyme_kgph
        self.design_results["Enzyme cost ($/h)"] = enzyme_cost_usd_per_hr

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Enzymatic pretreatment"] = capex
        # NOTE: `capex_usd` is sourced from an NREL saccharification-tank cost
        # anchor and does not include heat-exchange equipment.
        enzyme_cost_usd_per_hr = getattr(self, "_enzyme_cost_usd_per_hr", 0.0)

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        add_opex = {}
        if enzyme_cost_usd_per_hr > 0:
            add_opex["Enzyme reagent cost"] = enzyme_cost_usd_per_hr
        if annual_maintenance > 0:
            add_opex["Enzymatic pretreatment maintenance"] = (
                annual_maintenance / OPERATING_HOURS_PER_YEAR
            )
        if add_opex:
            self.add_OPEX = add_opex


class HeatingPretreatment(bst.Unit):
    """
    Thermal (heating) pretreatment reactor, one of the optional pre-AD
    pretreatment stages. Used in the `combined_PTE` case's peroxide->heating->
    enzymatic sequence.

    Parameters
    ----------
    ins : stream
        Biomass feed to be heated.
    outs : stream
        Heated biomass.
    target_temperature_K : float
        Target reactor outlet temperature; the duty to raise the feed
        (at its own stream temperature) to this value is charged as a
        heat utility.
    residence_time_hr : float
        Reactor residence time, used to size reactor volume.
    capex_usd : float
        Flat installed capital cost (not a function of throughput).
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of `capex_usd`, added to
        `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/pretreatment.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 1  # heated_biomass
    _F_BM_default = {"Heating pretreatment": _HEATING["F_BM"]}
    auxiliary_unit_names = ("heat_exchanger",)

    def __init__(
        self, ID="", ins=None, outs=(),
        target_temperature_K=_HEATING["target_temperature_K"],
        residence_time_hr=_HEATING["residence_time_hr"],

        # costing
        capex_usd=_HEATING["capex_usd"],
        maintenance_frac_of_capex_per_yr=_HEATING["maintenance_frac_of_capex_per_yr"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.target_temperature_K = float(target_temperature_K)
        self.residence_time_hr = float(residence_time_hr)

        self.capex_usd = float(capex_usd)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

    def _run(self):
        feed = self.ins[0]
        heated = self.outs[0]
        heated.copy_like(feed)
        heated.phase = feed.phase
        heated.T = self.target_temperature_K

    def _design(self):
        feed = self.ins[0]

        m_kgph = float(feed.F_mass)
        vol_m3ph = feed.F_vol
        V = vol_m3ph * self.residence_time_hr

        T_in = float(feed.T)
        T_out = self.target_temperature_K
        # Auxiliary HXutility: duty comes from BioSTEAM's own enthalpy balance
        # instead of a fixed single-point feed.Cp, and is HXN-eligible by default.
        hx = self.auxiliary("heat_exchanger", bst.HXutility, ins=self.ins[0], T=T_out)
        hx.simulate()

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Inlet T (K)"] = T_in
        self.design_results["Target T (K)"] = T_out
        self.design_results["Heating duty (kJ/h)"] = hx.net_duty

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Heating pretreatment"] = capex
        # NOTE: `capex_usd` has no cited cost source in data/pretreatment.yaml;
        # assumed here to NOT include heat-exchange
        # equipment.

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            self.add_OPEX = {
                "Heating pretreatment maintenance": annual_maintenance / OPERATING_HOURS_PER_YEAR
            }


class PeroxidePretreatment(bst.Unit):
    """
    Hydrogen peroxide pretreatment reactor, one of the optional pre-AD
    pretreatment stages. Used when the selected pretreatment case has `kind: peroxide` or
    `kind: combined_PE`/`combined_PTE` (which include a peroxide stage in
    sequence).

    Parameters
    ----------
    ins : stream
        Biomass feed to be peroxide-treated.
    outs : stream
        Peroxide-treated biomass.
    h2o2_wt_frac_on_dry_feed : float
        Hydrogen peroxide dose as a weight fraction of dry (non-water)
        feed processed.
    residence_time_hr : float
        Reactor residence time, used to size reactor volume.
    capex_usd : float
        Flat installed capital cost (not a function of throughput).
    h2o2_price_usd_per_kg : float
        Hydrogen peroxide reagent unit price, added to `add_OPEX`.
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of `capex_usd`, added to
        `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    Notes
    -----
    Feed slurry volumetric flow is read directly from the inlet stream
    (`feed.F_vol`), not from a fixed assumption.

    See Also
    --------
    Refer to data/pretreatment.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 1  # peroxide_treated_biomass
    _F_BM_default = {"Peroxide pretreatment": _PEROXIDE["F_BM"]}

    def __init__(
        self, ID="", ins=None, outs=(),
        h2o2_wt_frac_on_dry_feed=_PEROXIDE["h2o2_wt_frac_on_dry_feed"],
        residence_time_hr=_PEROXIDE["residence_time_hr"],

        # economics
        capex_usd=_PEROXIDE["capex_usd"],
        h2o2_price_usd_per_kg=_PEROXIDE["h2o2_price_usd_per_kg"],
        maintenance_frac_of_capex_per_yr=_PEROXIDE["maintenance_frac_of_capex_per_yr"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.h2o2_wt_frac_on_dry_feed = float(h2o2_wt_frac_on_dry_feed)
        self.residence_time_hr = float(residence_time_hr)

        self.capex_usd = float(capex_usd)
        self.h2o2_price_usd_per_kg = float(h2o2_price_usd_per_kg)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

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
        vol_m3ph = feed.F_vol
        V = vol_m3ph * self.residence_time_hr

        dry_mass = getattr(self, "_dry_mass_kgph", 0.0)
        h2o2_kgph = getattr(self, "_h2o2_kgph", 0.0)
        h2o2_cost_usd_per_hr = h2o2_kgph * self.h2o2_price_usd_per_kg
        self._h2o2_cost_usd_per_hr = h2o2_cost_usd_per_hr

        self.design_results["Mass flow (kg/h)"] = m_kgph
        self.design_results["Slurry flow (m3/h)"] = vol_m3ph
        self.design_results["Residence time (h)"] = self.residence_time_hr
        self.design_results["Reactor volume (m3)"] = V
        self.design_results["Dry feed basis (kg/h)"] = dry_mass
        self.design_results["H2O2 addition (kg/h)"] = h2o2_kgph
        self.design_results["H2O2 cost ($/h)"] = h2o2_cost_usd_per_hr

    def _cost(self):
        capex = self.capex_usd
        self.baseline_purchase_costs["Peroxide pretreatment"] = capex

        h2o2_cost_usd_per_hr = getattr(self, "_h2o2_cost_usd_per_hr", 0.0)

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        add_opex = {}
        if h2o2_cost_usd_per_hr > 0:
            add_opex["H2O2 reagent cost"] = h2o2_cost_usd_per_hr
        if annual_maintenance > 0:
            add_opex["Peroxide pretreatment maintenance"] = (
                annual_maintenance / OPERATING_HOURS_PER_YEAR
            )
        if add_opex:
            self.add_OPEX = add_opex

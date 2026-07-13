# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Oil Extraction Placeholder Unit
================================
Represents cell disruption + lipid extraction for microbial oil recovery
from Yarrowia lipolytica fermentation broth

All separation is handled downstream by C603_2 (LiquidsSplitCentrifuge).
Its purpose is to carry realistic capital and operating costs for:
    1. Cell disruption (high-pressure homogenization)
    2. Solvent or aqueous extraction of intracellular lipids

Capital cost anchor:
    NREL/TP-5100-55431 (2012), Davis et al., "Techno-Economic Analysis of
    Autotrophic Microalgae for Fuel Production". Cell disruption + extraction
    section reported at ~$8.5M installed for ~10 dry ton/hr biomass throughput.
    Scaled here using a 0.6 power law on dry biomass flow.

Operating cost:
    - Electricity: high-pressure homogenization ~1.5 kWh/kg dry biomass
      (Doucha & Livansky 2008; Postma et al. 2017)
    - Reagent/solvent cost included as $/kg oil surrogate operating cost,
      passed through as a fixed annual cost adder in the TEA script
      Literature range: $0.50-1.50/kg oil for solvent + recovery
      (Knoshaug et al. 2018, NREL; Laurens et al. 2017, Green Chem.)
"""

from __future__ import annotations
import biosteam as bst

__all__ = ('OilExtractionPlaceholder',)


class OilExtractionPlaceholder(bst.Unit):
    """
    Pass-through unit representing cell disruption and lipid extraction

    Inputs:
        ins[0]: concentrated fermentation broth (from upstream pump/evaporator)

    Outputs:
        outs[0]: extracted broth (same composition — split handled by C603_2)

    """

    _N_ins = 1
    _N_outs = 1
    _F_BM_default = {"Oil extraction": 1.0}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        product_ID: str = "MicrobialOil",
        cellmass_ID: str = "CellMass",
        homogenization_kWh_per_kg_dry_biomass: float = 1.5,
        ref_dry_biomass_tph: float = 10.0,
        ref_installed_cost_usd: float = 8_500_000.0,
        scale_exponent: float = 0.6,
        F_BM: float = 1.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.product_ID = product_ID
        self.cellmass_ID = cellmass_ID
        self.homogenization_kWh_per_kg_dry_biomass = float(
            homogenization_kWh_per_kg_dry_biomass
        )
        self.ref_dry_biomass_tph = float(ref_dry_biomass_tph)
        self.ref_installed_cost_usd = float(ref_installed_cost_usd)
        self.scale_exponent = float(scale_exponent)
        self.F_BM = {"Oil extraction": float(F_BM)}

    def _run(self):
        # Pass-through: composition unchanged.
        # Separation is handled downstream by C603_2.
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)

    def _design(self):
        feed = self.ins[0]
        chem_ids = set(feed.chemicals.IDs)

        # Dry biomass = cell mass + any residual oil in feed
        dry_biomass_kgph = 0.0
        if self.cellmass_ID in chem_ids:
            dry_biomass_kgph += float(feed.imass[self.cellmass_ID])
        if self.product_ID in chem_ids:
            dry_biomass_kgph += float(feed.imass[self.product_ID])

        dry_biomass_tph = dry_biomass_kgph / 1000.0

        # Homogenization power (kW)
        homogenization_kW = (
            self.homogenization_kWh_per_kg_dry_biomass * dry_biomass_kgph
        )
        self.power_utility.consumption = homogenization_kW

        # Oil produced
        oil_kgph = float(feed.imass[self.product_ID]) if self.product_ID in chem_ids else 0.0

        self.design_results["Feed flow (kg/h)"] = feed.F_mass
        self.design_results["Dry biomass feed (kg/h)"] = dry_biomass_kgph
        self.design_results["Dry biomass feed (dry ton/h)"] = dry_biomass_tph
        self.design_results["Oil in feed (kg/h)"] = oil_kgph
        self.design_results[
            "Homogenization power (kW)"
        ] = homogenization_kW
        self.design_results[
            "Electricity intensity (kWh/kg dry biomass)"
        ] = self.homogenization_kWh_per_kg_dry_biomass

    def _cost(self):
        dry_biomass_tph = (
            self.design_results.get("Dry biomass feed (dry ton/h)", 0.0)
        )
        ref = self.ref_dry_biomass_tph
        ref_cost = self.ref_installed_cost_usd

        if dry_biomass_tph <= 0 or ref <= 0:
            installed_cost = ref_cost
        else:
            installed_cost = ref_cost * (dry_biomass_tph / ref) ** self.scale_exponent

        self.baseline_purchase_costs["Oil extraction"] = installed_cost
        self.design_results["Installed cost ($)"] = installed_cost

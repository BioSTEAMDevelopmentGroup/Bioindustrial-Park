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

Capital cost anchor (UNVERIFIED -- see CONVERSION_NOTES.md, "Citation-accuracy
issues" -- nobody involved in this conversion has read the primary source):
    assumptions.yaml states the $7.848M installed cost (at 10 dry ton/hr
    biomass throughput) is "back-calculated from NREL 55431 cell-disruption
    plus extraction/separation cost bases," without showing the calculation.
    A separate, pre-existing docstring in this codebase attributed the figure
    to NREL/TP-5100-55431 (2012), Davis et al., "Techno-Economic Analysis of
    Autotrophic Microalgae for Fuel Production," and claimed that report
    itself gives ~$8.5M for the same throughput. Neither claim has been
    checked against the actual report. Scaled here using a 0.6 power law on
    dry biomass flow.

Operating cost:
    - Electricity: high-pressure homogenization 0.203 kWh/kg dry biomass.
      assumptions.yaml describes this only as "harmonized to the NREL
      algal-lipid basis," no specific paper named. Doucha & Livansky (2008)
      and Postma et al. (2017) are cited elsewhere in this codebase for the
      superseded 1.5 kWh/kg figure -- NOT verified to support 0.203.
      UNVERIFIED, see CONVERSION_NOTES.md ("Citation-accuracy issues").
    - Reagent/solvent cost included as $/kg oil surrogate operating cost,
      passed through as a fixed annual cost adder in the TEA script
      Literature range: $0.50-1.50/kg oil for solvent + recovery
      (Knoshaug et al. 2018, NREL; Laurens et al. 2017, Green Chem.) --
      also not independently verified against the primary sources.
"""

import biosteam as bst
from biosteam.units.decorators import cost

from biorefineries.sabre.utils import load_assumptions

__all__ = ('OilExtractionPlaceholder',)

# Loaded yaml assumptions
_VFA_FERMENTATION_YAML = load_assumptions("vfa_fermentation.yaml")
_DOWNSTREAM_RECOVERY = _VFA_FERMENTATION_YAML["vfa_fermentation"]["downstream_recovery"]


@cost('Dry biomass feed (dry ton/h)', 'Oil extraction', units='dry ton/h',
      CE=567.5, cost=_DOWNSTREAM_RECOVERY["oil_extraction_ref_installed_cost_usd"],
      S=_DOWNSTREAM_RECOVERY["oil_extraction_ref_dry_biomass_tph"],
      n=_DOWNSTREAM_RECOVERY["oil_extraction_scale_exponent"], BM=1.0)
class OilExtractionPlaceholder(bst.Unit):
    """
    Pass-through unit representing cell disruption and lipid extraction

    Inputs:
        ins[0]: concentrated fermentation broth (from upstream pump/evaporator)

    Outputs:
        outs[0]: extracted broth (same composition — split handled by C603_2)

    Sources
    -------
    installed cost anchor ($7.848M at 10 dry ton/h, n=0.6):
        assumptions.yaml states this is back-calculated from NREL/TP-5100-55431
        (2012) cell-disruption plus extraction/separation cost bases; the
        calculation itself is not shown/verified.
    homogenization_kWh_per_kg_dry_biomass:
        UNVERIFIED. assumptions.yaml says only "harmonized to the NREL
        algal-lipid basis." A pre-existing docstring in this codebase
        attributed Doucha & Livansky (2008) and Postma et al. (2017) to the
        superseded 1.5 kWh/kg figure this replaced -- that attribution itself
        was never checked against the papers either. See CONVERSION_NOTES.md
        ("Citation-accuracy issues") before citing either paper for anything
        here.
    """

    _N_ins = 1
    _N_outs = 1
    _units = {'Dry biomass feed (dry ton/h)': 'dry ton/h'}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        product_ID: str = "MicrobialOil",
        cellmass_ID: str = "CellMass",
        homogenization_kWh_per_kg_dry_biomass: float = 0.203,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.product_ID = product_ID
        self.cellmass_ID = cellmass_ID
        self.homogenization_kWh_per_kg_dry_biomass = float(
            homogenization_kWh_per_kg_dry_biomass
        )

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

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam.units.decorators import cost

from biorefineries.sabre.utils import load_assumptions

__all__ = ('OilExtraction',)

# Loaded yaml assumptions (downstream_processing.yaml was split into
# vfa.yaml and microbial_oil.yaml -- oil_extraction now lives here).
_MICROBIAL_OIL_YAML = load_assumptions("microbial_oil.yaml")
_OIL_EXTRACTION = _MICROBIAL_OIL_YAML["oil_extraction"]


@cost('Dry biomass feed (dry ton/h)', 'Oil extraction', units='dry ton/h',
      CE=567.5, cost=_OIL_EXTRACTION["ref_installed_cost_usd"],
      S=_OIL_EXTRACTION["ref_dry_biomass_tph"],
      n=_OIL_EXTRACTION["scale_exponent"], BM=_OIL_EXTRACTION["F_BM"])
class OilExtraction(bst.Unit):
    """
    Cell disruption and lipid extraction for microbial oil recovery
    from Yarrowia lipolytica fermentation broth.
    All separation should be handled downstream.

    Parameters
    ----------
    ins : stream
        Concentrated fermentation broth (from upstream pump/evaporator).
    outs : stream
        Extracted broth (same composition; split should be handled downstream).
    product_ID : str
        Chemical ID of the fermentation product (oil) in the feed.
    cellmass_ID : str
        Chemical ID of cell mass in the feed.
    homogenization_kWh_per_kg_dry_biomass : float
        Electricity intensity for high-pressure homogenization, per kg
        of dry biomass (cell mass + product) in the feed.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/downstream_processing.yaml for the default values and references.
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
        product_ID: str = _OIL_EXTRACTION["product_ID"],
        cellmass_ID: str = _OIL_EXTRACTION["cellmass_ID"],
        homogenization_kWh_per_kg_dry_biomass: float = _OIL_EXTRACTION["homogenization_kWh_per_kg_dry_biomass"],
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
        self.power_utility(homogenization_kW)

        # Oil produced
        oil_kgph = float(feed.imass[self.product_ID]) if self.product_ID in chem_ids else 0.0

        self.design_results["Feed flow (kg/h)"] = feed.F_mass
        self.design_results["Dry biomass feed (kg/h)"] = dry_biomass_kgph
        self.design_results["Dry biomass feed (dry ton/h)"] = dry_biomass_tph
        self.design_results["Oil in feed (kg/h)"] = oil_kgph
        self.design_results["Homogenization power (kW)"] = homogenization_kW
        self.design_results[
            "Electricity intensity (kWh/kg dry biomass)"
        ] = self.homogenization_kWh_per_kg_dry_biomass

        self.add_OPEX = {
            "Oil extraction reagent": oil_kgph * _OIL_EXTRACTION["oil_extraction_reagent_usd_per_kg_oil"]
        }
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam.units.decorators import cost

from biorefineries.sabre.utils import load_assumptions, get_solids_group_IDs

__all__ = ('VFAMicrofilter', 'OilExtractionPlaceholder')

# Loaded yaml assumptions
_DOWNSTREAM_PROCESSING_YAML = load_assumptions("downstream_processing.yaml")
_VFA_MICROFILTER = _DOWNSTREAM_PROCESSING_YAML["vfa_microfilter"]
_OIL_EXTRACTION = _DOWNSTREAM_PROCESSING_YAML["oil_extraction"]


@cost('Membrane area (m2)', 'Microfilter', units='m2',
      CE=567.5, cost=_VFA_MICROFILTER["membrane_cost_usd_per_m2"], S=1., n=1., BM=_VFA_MICROFILTER["F_BM"])
class VFAMicrofilter(bst.Unit):
    """
    Split-based representation of a VFA-rich permeate step.
    Includes first-pass power draw and area-based membrane cost.

    Sources
    -------
    design_flux_L_m2_h, membrane_cost_usd_per_m2:
        assumptions.yaml gives no citation for these two values (20 LMH,
        $200/m2) beyond descriptive notes -- "20 LMH selected as a
        conservative design flux for low-pressure clarification" and "$200/m2
        selected as an upper-end polymeric membrane cost proxy." No named
        source. (Do not confuse with the *different* PressateConcentrator
        unit's 35 LMH flux, which yaml does cite to seaweed/algal membrane
        concentration literature -- Sievers et al. 2017; Diaz-Reinoso and
        Dominguez 2020 -- that citation does not apply here.)
    """
    _N_ins = 1
    _N_outs = 2
    _units = {'Membrane area (m2)': 'm2'}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        vfa_IDs=None,
        solids_IDs=None,
        vfa_to_permeate_frac: float = 0.98,
        water_to_permeate_frac: float = 0.97,
        solids_to_permeate_frac: float = 0.05,
        dissolved_other_to_permeate_frac: float = 0.90,
        broth_density_kg_per_m3: float = 1000.0,
        SEC_kWh_per_m3_feed: float = 0.08,
        design_flux_L_m2_h: float = 20.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.vfa_IDs = tuple(vfa_IDs or [
            "AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"
        ])
        if solids_IDs is None:
            solids_IDs = get_solids_group_IDs(self.chemicals)
        self.solids_IDs = tuple(solids_IDs)
        self.vfa_to_permeate_frac = float(vfa_to_permeate_frac)
        self.water_to_permeate_frac = float(water_to_permeate_frac)
        self.solids_to_permeate_frac = float(solids_to_permeate_frac)
        self.dissolved_other_to_permeate_frac = float(dissolved_other_to_permeate_frac)
        self.broth_density_kg_per_m3 = float(broth_density_kg_per_m3)
        self.SEC_kWh_per_m3_feed = float(SEC_kWh_per_m3_feed)
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)

    def _run(self):
        feed = self.ins[0]
        permeate, retentate = self.outs

        permeate.empty()
        retentate.empty()
        permeate.phase = "l"
        retentate.phase = "l"

        for cid in feed.chemicals.IDs:
            m = float(feed.imass[cid])
            if m <= 0:
                continue

            if cid in self.vfa_IDs:
                frac = self.vfa_to_permeate_frac
            elif cid == "Water":
                frac = self.water_to_permeate_frac
            elif cid in self.solids_IDs:
                frac = self.solids_to_permeate_frac
            else:
                frac = self.dissolved_other_to_permeate_frac

            frac = min(max(frac, 0.0), 1.0)
            permeate.imass[cid] = m * frac
            retentate.imass[cid] = m * (1.0 - frac)

    def _design(self):
        feed = self.ins[0]
        feed_m3h = (
            feed.F_mass / self.broth_density_kg_per_m3
            if self.broth_density_kg_per_m3 > 0 else 0.0
        )
        membrane_area_m2 = 0.0
        if self.design_flux_L_m2_h > 0:
            membrane_area_m2 = feed_m3h * 1000.0 / self.design_flux_L_m2_h

        self.design_results["Feed flow (kg/h)"] = feed.F_mass
        self.design_results["Feed flow (m3/h)"] = feed_m3h
        self.design_results["Permeate flow (kg/h)"] = self.outs[0].F_mass
        self.design_results["Retentate flow (kg/h)"] = self.outs[1].F_mass
        self.design_results["Membrane area (m2)"] = membrane_area_m2
        self.power_utility.consumption = self.SEC_kWh_per_m3_feed * feed_m3h


@cost('Dry biomass feed (dry ton/h)', 'Oil extraction', units='dry ton/h',
      CE=567.5, cost=_OIL_EXTRACTION["oil_extraction_ref_installed_cost_usd"],
      S=_OIL_EXTRACTION["oil_extraction_ref_dry_biomass_tph"],
      n=_OIL_EXTRACTION["oil_extraction_scale_exponent"], BM=_OIL_EXTRACTION["oil_extraction_F_BM"])
class OilExtractionPlaceholder(bst.Unit):
    """
    Pass-through unit representing cell disruption and lipid extraction

    Inputs:
        ins[0]: concentrated fermentation broth (from upstream pump/evaporator)

    Outputs:
        outs[0]: extracted broth (same composition — split handled by C603_2)

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

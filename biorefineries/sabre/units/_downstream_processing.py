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

__all__ = ('VFAMicrofilter', 'OilExtraction')

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

    Parameters
    ----------
    ins : stream
        VFA-rich broth feed.
    outs : tuple[stream, stream]
        Permeate and retentate.
    vfa_IDs : Iterable[str]
        Chemical IDs treated as VFA.
    solids_IDs : Iterable[str]
        Chemical IDs treated as solids for the split. Defaults to the
        "solids" chemical group when not given (see
        `utils.get_solids_group_IDs`).
    vfa_to_permeate_frac : float
        Fraction of VFA mass routed to permeate (0-1).
    water_to_permeate_frac : float
        Fraction of water mass routed to permeate (0-1).
    solids_to_permeate_frac : float
        Fraction of solids-ID mass routed to permeate (0-1).
    dissolved_other_to_permeate_frac : float
        Fraction of remaining dissolved mass routed to permeate (0-1).
    SEC_kWh_per_m3_feed : float
        Specific electricity consumption per m3 of feed.
    design_flux_L_m2_h : float
        Design membrane flux, used to size membrane area.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/downstream_processing.yaml for the default values and references.
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
        vfa_IDs: list[str] = _VFA_MICROFILTER["vfa_IDs"],
        solids_IDs=None,  # if not given, defaults to the "solids" chemical group (see utils.get_solids_group_IDs)
        vfa_to_permeate_frac: float = _VFA_MICROFILTER["vfa_to_permeate_frac"],
        water_to_permeate_frac: float = _VFA_MICROFILTER["water_to_permeate_frac"],
        solids_to_permeate_frac: float = _VFA_MICROFILTER["solids_to_permeate_frac"],
        dissolved_other_to_permeate_frac: float = _VFA_MICROFILTER["dissolved_other_to_permeate_frac"],
        SEC_kWh_per_m3_feed: float = _VFA_MICROFILTER["SEC_kWh_per_m3_feed"],
        design_flux_L_m2_h: float = _VFA_MICROFILTER["design_flux_L_m2_h"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.vfa_IDs = tuple(vfa_IDs)
        if solids_IDs is None:
            solids_IDs = get_solids_group_IDs(self.chemicals)
        self.solids_IDs = tuple(solids_IDs)
        self.vfa_to_permeate_frac = float(vfa_to_permeate_frac)
        self.water_to_permeate_frac = float(water_to_permeate_frac)
        self.solids_to_permeate_frac = float(solids_to_permeate_frac)
        self.dissolved_other_to_permeate_frac = float(dissolved_other_to_permeate_frac)
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
        feed_m3h = feed.F_vol
        membrane_area_m2 = 0.0
        if self.design_flux_L_m2_h > 0:
            membrane_area_m2 = feed_m3h * 1000.0 / self.design_flux_L_m2_h

        self.design_results["Feed flow (kg/h)"] = feed.F_mass
        self.design_results["Feed flow (m3/h)"] = feed_m3h
        self.design_results["Permeate flow (kg/h)"] = self.outs[0].F_mass
        self.design_results["Retentate flow (kg/h)"] = self.outs[1].F_mass
        self.design_results["Membrane area (m2)"] = membrane_area_m2
        self.power_utility(self.SEC_kWh_per_m3_feed * feed_m3h)


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

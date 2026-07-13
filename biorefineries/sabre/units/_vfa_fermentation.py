# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
VFA-fed microbial-oil fermenter.
"""
from __future__ import annotations

from typing import Optional
import biosteam as bst

__all__ = ('YarrowiaLipidFermenter',)


class YarrowiaLipidFermenter(bst.AeratedBioreactor):
    """
    VFA-fed microbial-oil fermenter structured after the oilcane
    AeratedFermentation class

    This is now using the *reactor architecture* of the oilcane code
    (AeratedBioreactor with automatic aeration / power optimization),
    but it still uses yield-based VFA conversion because MicrobialOil
    is a pseudo-component in the current Sabre thermo

    Inputs:
    - conditioned VFA broth
    - optional seed stream

    Outputs:
    - vent
    - fermentation broth
    """

    V_max_default = 150

    def _init(
        self,
        vfa_IDs: Optional[list[str]] = None,
        product_ID: str = "MicrobialOil",
        conversion: float = 0.85,
        product_yield_kg_per_kg_vfa_consumed: float = 0.45,
        biomass_yield_kg_per_kg_vfa_consumed: float = 0.10,
        co2_yield_kg_per_kg_vfa_consumed: float = 0.20,
        oxygen_kg_per_kg_vfa_consumed: float = 0.80,
        residence_time_h: float = 48.0,
        dT_hx_loop: float = 8.0,
        Q_O2_consumption: float = -460240.0,
        batch: bool = True,
        broth_density_kg_per_m3: float = 1000.0,
        target_pH: float = 8.0,

        # real stirred-tank sizing / power args
        V_wf: float = 0.8,
        V_max: float = 300.0,
        kW_per_m3: float = 0.06,
        tau_0: float = 3.4,
        N_reactors: int | None = None,

        **kwargs,
    ):
        bst.AeratedBioreactor._init(
        self,
        reactions=None,
        batch=batch,
        dT_hx_loop=dT_hx_loop,
        Q_O2_consumption=Q_O2_consumption,
        optimize_power=True,
        tau=residence_time_h,
        V_wf=V_wf,
        V_max=V_max,
        kW_per_m3=kW_per_m3,
        tau_0=tau_0,
        **kwargs,
    )

        self.vfa_IDs = vfa_IDs or [
            "AceticAcid",
            "PropionicAcid",
            "ButyricAcid",
            "ValericAcid",
            "HexanoicAcid",
        ]
        self.product_ID = product_ID
        self.conversion = float(conversion)
        self.product_yield_kg_per_kg_vfa_consumed = float(product_yield_kg_per_kg_vfa_consumed)
        self.biomass_yield_kg_per_kg_vfa_consumed = float(biomass_yield_kg_per_kg_vfa_consumed)
        self.co2_yield_kg_per_kg_vfa_consumed = float(co2_yield_kg_per_kg_vfa_consumed)
        self.oxygen_kg_per_kg_vfa_consumed = float(oxygen_kg_per_kg_vfa_consumed)
        self.broth_density_kg_per_m3 = float(broth_density_kg_per_m3)
        self.target_pH = float(target_pH)

    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ("CarbonDioxide", "Oxygen", "Nitrogen"), remove=True)

    def _run_reactions(self, effluent):
        ids = set(self.chemicals.IDs)

        required = (
            self.product_ID,
            "CellMass",
            "CarbonDioxide",
            "Water",
            "Oxygen",
            "Nitrogen",
        )
        missing = [i for i in required if i not in ids]
        if missing:
            raise RuntimeError(f"Missing required chemicals in thermo: {missing}")

        # Total available VFAs (mass basis)
        vfa_available = 0.0
        for cid in self.vfa_IDs:
            if cid in ids:
                vfa_available += float(effluent.imass[cid])

        if vfa_available <= 1e-12:
            self.design_results["Target pH"] = self.target_pH
            self.design_results["VFA available (kg/h)"] = 0.0
            self.design_results["VFA consumed (kg/h)"] = 0.0
            self.design_results["Microbial oil formed (kg/h)"] = 0.0
            self.design_results["Biomass formed (kg/h)"] = 0.0
            self.design_results["CO2 formed (kg/h)"] = 0.0
            self.design_results["O2 demand (kg/h)"] = 0.0
            return

        vfa_consumed = self.conversion * vfa_available
        product_formed = self.product_yield_kg_per_kg_vfa_consumed * vfa_consumed
        biomass_formed = self.biomass_yield_kg_per_kg_vfa_consumed * vfa_consumed
        co2_formed = self.co2_yield_kg_per_kg_vfa_consumed * vfa_consumed
        o2_demand = self.oxygen_kg_per_kg_vfa_consumed * vfa_consumed

        accounted = product_formed + biomass_formed + co2_formed
        unaccounted = max(0.0, vfa_consumed - accounted)

        # Remove VFAs in proportion to composition
        for cid in self.vfa_IDs:
            if cid in ids:
                m = float(effluent.imass[cid])
                take = vfa_consumed * (m / vfa_available) if vfa_available > 0 else 0.0
                effluent.imass[cid] -= take

        # Form product, biomass, CO2
        effluent.imass[self.product_ID] += product_formed
        effluent.imass["CellMass"] += biomass_formed
        effluent.imass["CarbonDioxide"] += co2_formed

        # Consume oxygen supplied by AeratedBioreactor
        available_o2 = float(effluent.imass["Oxygen"])
        effluent.imass["Oxygen"] = max(0.0, available_o2 - o2_demand)

        # Put leftover converted mass into soluble organics if available
        if "SolubleOrganics" in ids:
            effluent.imass["SolubleOrganics"] += unaccounted
        else:
            effluent.imass["Water"] += unaccounted

        if effluent.imass["Water"] < 0.0:
            effluent.imass["Water"] = 0.0

        self.design_results["Target pH"] = self.target_pH
        self.design_results["VFA available (kg/h)"] = vfa_available
        self.design_results["VFA consumed (kg/h)"] = vfa_consumed
        self.design_results["Microbial oil formed (kg/h)"] = product_formed
        self.design_results["Biomass formed (kg/h)"] = biomass_formed
        self.design_results["CO2 formed (kg/h)"] = co2_formed
        self.design_results["O2 demand (kg/h)"] = o2_demand

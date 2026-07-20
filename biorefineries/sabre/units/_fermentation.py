# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from typing import Optional
import biosteam as bst
from biosteam.units.design_tools.tank_design import mix_tank_purchase_cost_algorithms

from biorefineries.sabre.utils import load_assumptions

__all__ = ('YarrowiaLipidFermenter', 'FermentationMediumTank')

# Loaded yaml assumptions
_FERMENTATION_YAML = load_assumptions("fermentation.yaml")
_VFA = _FERMENTATION_YAML["vfa"]
_VFA_CASE = _VFA["cases"][_VFA["case"]]
_FERMENTATION_MEDIUM_TANK = _FERMENTATION_YAML["fermentation_medium_tank"]


class YarrowiaLipidFermenter(bst.AeratedBioreactor):
    """
    VFA-fed microbial-oil fermenter structured after the oilcane
    AeratedFermentation class. It uses AeratedBioreactor
    with automatic aeration / power optimization),
    but it still uses yield-based VFA conversion because MicrobialOil
    is a pseudo-component in the current Sabre thermo

    Parameters
    ----------
    ins : stream
        Conditioned VFA broth, recycle biomass, and seed makeup.
    outs : tuple[stream, stream]
        Vent and fermentation broth.
    vfa_IDs : list[str]
        Chemical IDs treated as VFA substrate.
    product_ID : str
        Chemical ID of the fermentation product.
    conversion : float
        Fraction of available VFA mass consumed (0-1).
    product_yield_kg_per_kg_vfa_consumed : float
        Product yield per kg of VFA consumed.
    biomass_yield_kg_per_kg_vfa_consumed : float
        Cell mass yield per kg of VFA consumed.
    co2_yield_kg_per_kg_vfa_consumed : float
        CO2 yield per kg of VFA consumed.
    oxygen_kg_per_kg_vfa_consumed : float
        O2 demand per kg of VFA consumed.
    tau : float
        Fermentation residence time [hr] (bst.AeratedBioreactor's own
        attribute name).
    Q_O2_consumption : float
        Heat of reaction per kmol O2 consumed; forwarded to
        `bst.AeratedBioreactor`. Not yaml-sourced.
        Cannot be calculated from the current sabre thermo because MicrobialOil is a pseudo-component.
    target_pH : float
        Target fermentation pH; recorded in `design_results`, not
        enforced.
    seed_water_kgph : float
        Water flow (kg/hr) for the `seed` influent stream.
    seed_cellmass_kgph : float
        Cell mass flow (kg/hr) for the `seed` influent stream.
    **kwargs
        Forwarded to `bst.AeratedBioreactor.__init__`.

    See Also
    --------
    Refer to data/fermentation.yaml for the default values and references.
    """

    def _init(
        self,
        vfa_IDs: list[str] = _VFA["vfa_IDs"],
        product_ID: str = _VFA_CASE["product_ID"],
        conversion: float = _VFA_CASE["conversion"],
        product_yield_kg_per_kg_vfa_consumed: float = _VFA_CASE["product_yield_kg_per_kg_vfa_consumed"],
        biomass_yield_kg_per_kg_vfa_consumed: float = _VFA_CASE["biomass_yield_kg_per_kg_vfa_consumed"],
        co2_yield_kg_per_kg_vfa_consumed: float = _VFA_CASE["co2_yield_kg_per_kg_vfa_consumed"],
        oxygen_kg_per_kg_vfa_consumed: float = _VFA_CASE["oxygen_kg_per_kg_vfa_consumed"],
        tau: float = _VFA_CASE["tau"],
        Q_O2_consumption: float = -460240.0,
        target_pH: float = _VFA_CASE["target_pH"],
        seed_water_kgph: float = _VFA_CASE["seed_water_kgph"],
        seed_cellmass_kgph: float = _VFA_CASE["seed_cellmass_kgph"],
        **kwargs,
    ):
        kwargs.setdefault('optimize_power', True)
        bst.AeratedBioreactor._init(
        self,
        reactions=None,
        Q_O2_consumption=Q_O2_consumption,
        tau=tau,
        **kwargs,
    )

        self.vfa_IDs = list(vfa_IDs)
        self.product_ID = product_ID
        self.conversion = float(conversion)
        self.product_yield_kg_per_kg_vfa_consumed = float(product_yield_kg_per_kg_vfa_consumed)
        self.biomass_yield_kg_per_kg_vfa_consumed = float(biomass_yield_kg_per_kg_vfa_consumed)
        self.co2_yield_kg_per_kg_vfa_consumed = float(co2_yield_kg_per_kg_vfa_consumed)
        self.oxygen_kg_per_kg_vfa_consumed = float(oxygen_kg_per_kg_vfa_consumed)
        self.target_pH = float(target_pH)
        self.seed_water_kgph = float(seed_water_kgph)
        self.seed_cellmass_kgph = float(seed_cellmass_kgph)
        self.seed = seed = bst.Stream(None, thermo=self.thermo)
        self.ins.append(seed)

    def _run(self):
        # Fresh seed makeup, added as a genuine influent stream instead of
        # via an external pre-built stream + mixer. Recomputed every run so
        # it's empty when both kgph values are 0.
        seed = self.seed
        seed.empty()
        seed.imass["Water"] = self.seed_water_kgph
        seed.imass["CellMass"] = self.seed_cellmass_kgph
        super()._run()

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

        # Yields don't sum to 1, so only this much of vfa_consumed is actually
        # explained by tracked products; the remainder is left as VFA rather
        # than converted into a fictitious byproduct (previously SolubleOrganics).
        accounted = product_formed + biomass_formed + co2_formed

        # Remove accounted VFA mass in proportion to composition
        for cid in self.vfa_IDs:
            if cid in ids:
                m = float(effluent.imass[cid])
                take = accounted * (m / vfa_available) if vfa_available > 0 else 0.0
                effluent.imass[cid] -= take

        # Form product, biomass, CO2
        effluent.imass[self.product_ID] += product_formed
        effluent.imass["CellMass"] += biomass_formed
        effluent.imass["CarbonDioxide"] += co2_formed

        # Consume oxygen supplied by AeratedBioreactor
        available_o2 = float(effluent.imass["Oxygen"])
        effluent.imass["Oxygen"] = max(0.0, available_o2 - o2_demand)

        self.design_results["Target pH"] = self.target_pH
        self.design_results["VFA available (kg/h)"] = vfa_available
        self.design_results["VFA consumed (kg/h)"] = accounted
        self.design_results["Microbial oil formed (kg/h)"] = product_formed
        self.design_results["Biomass formed (kg/h)"] = biomass_formed
        self.design_results["CO2 formed (kg/h)"] = co2_formed
        self.design_results["O2 demand (kg/h)"] = o2_demand


class FermentationMediumTank(bst.Tank):
    """
    Simple medium adjustment tank for pH control and nutrient dosing.
    Reuses BioSTEAM's own conventional stainless-steel mix-tank
    purchase-cost algorithm. Only `_run()` (dosing) is overridden;
    `_design()` and `_cost()` are bst.Tank's own, unmodified.

    Parameters
    ----------
    ins : tuple[stream, stream, stream, stream, stream]
        Broth, ammonia, phosphate, base, and magnesium sulfate feeds.
    outs : stream
        Conditioned broth.
    ammonia_dose_kg_per_m3 : float
        Ammonia dose per m3 of broth (mass basis); 0 disables dosing.
    phosphate_dose_kg_per_m3 : float
        KH2PO4 dose per m3 of broth; 0 disables dosing.
    base_dose_kg_per_m3 : float
        NaOH dose per m3 of broth; 0 disables dosing.
    magnesium_sulfate_dose_kg_per_m3 : float
        MagnesiumSulfate dose per m3 of broth; 0 disables dosing.
    target_pH : float
        Target pH; recorded in `design_results`, not enforced.
    tau : float
        Tank residence time [hr] (bst.Tank's own attribute name).
    mixing_kW_per_m3 : float
        Mixing power intensity.
    **kwargs
        Forwarded to `bst.Tank.__init__`.

    See Also
    --------
    Refer to data/fermentation.yaml for the default values and references.
    """

    _N_ins = 5
    _N_outs = 1
    purchase_cost_algorithms = mix_tank_purchase_cost_algorithms
    _default_vessel_type = 'Conventional'
    _default_tau = 1
    _default_V_wf = 0.8
    _default_vessel_material = 'Stainless steel'

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        ammonia_dose_kg_per_m3: float = _FERMENTATION_MEDIUM_TANK["ammonia_dose_kg_per_m3"],
        phosphate_dose_kg_per_m3: float = _FERMENTATION_MEDIUM_TANK["phosphate_dose_kg_per_m3"],
        base_dose_kg_per_m3: float = _FERMENTATION_MEDIUM_TANK["base_dose_kg_per_m3"],
        magnesium_sulfate_dose_kg_per_m3: float = _FERMENTATION_MEDIUM_TANK["magnesium_sulfate_dose_kg_per_m3"],
        target_pH: float = _FERMENTATION_MEDIUM_TANK["target_pH"],
        tau: float = _FERMENTATION_MEDIUM_TANK["tau"],
        mixing_kW_per_m3: float = _FERMENTATION_MEDIUM_TANK["mixing_kW_per_m3"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.ammonia_dose_kg_per_m3 = float(ammonia_dose_kg_per_m3)
        self.phosphate_dose_kg_per_m3 = float(phosphate_dose_kg_per_m3)
        self.base_dose_kg_per_m3 = float(base_dose_kg_per_m3)
        self.magnesium_sulfate_dose_kg_per_m3 = float(magnesium_sulfate_dose_kg_per_m3)
        self.target_pH = float(target_pH)
        self.tau = float(tau)
        self.kW_per_m3 = float(mixing_kW_per_m3)

    def _run(self):
        broth, ammonia, phosphate, base, mgso4 = self.ins
        out = self.outs[0]
        out.copy_like(broth)
        out.phase = "l"

        vol_m3ph = broth.F_vol
        ammonia.empty(); phosphate.empty(); base.empty(); mgso4.empty()

        chem_ids = set(out.chemicals.IDs)
        if self.ammonia_dose_kg_per_m3 > 0 and "Ammonia" not in chem_ids:
            raise RuntimeError("Ammonia dose specified but 'Ammonia' is not in thermo.")
        if self.phosphate_dose_kg_per_m3 > 0 and "KH2PO4" not in chem_ids:
            raise RuntimeError("Phosphate dose specified but 'KH2PO4' is not in thermo.")
        if self.base_dose_kg_per_m3 > 0 and "NaOH" not in chem_ids:
            raise RuntimeError("Base dose specified but 'NaOH' is not in thermo.")
        if self.magnesium_sulfate_dose_kg_per_m3 > 0 and "MagnesiumSulfate" not in chem_ids:
            raise RuntimeError("MgSO4 dose specified but 'MagnesiumSulfate' is not in thermo.")

        if "Ammonia" in chem_ids and self.ammonia_dose_kg_per_m3 > 0:
            ammonia.imass["Ammonia"] = self.ammonia_dose_kg_per_m3 * vol_m3ph
            out.imass["Ammonia"] += ammonia.imass["Ammonia"]
        if "KH2PO4" in chem_ids and self.phosphate_dose_kg_per_m3 > 0:
            phosphate.imass["KH2PO4"] = self.phosphate_dose_kg_per_m3 * vol_m3ph
            out.imass["KH2PO4"] += phosphate.imass["KH2PO4"]
        if "NaOH" in chem_ids and self.base_dose_kg_per_m3 > 0:
            base.imass["NaOH"] = self.base_dose_kg_per_m3 * vol_m3ph
            out.imass["NaOH"] += base.imass["NaOH"]
        if "MagnesiumSulfate" in chem_ids and self.magnesium_sulfate_dose_kg_per_m3 > 0:
            mgso4.imass["MagnesiumSulfate"] = self.magnesium_sulfate_dose_kg_per_m3 * vol_m3ph
            out.imass["MagnesiumSulfate"] += mgso4.imass["MagnesiumSulfate"]

        self.design_results["Target pH"] = self.target_pH
        self.design_results["Broth flow (m3/h)"] = vol_m3ph

# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
VFA fermentation train: medium conditioning and the Yarrowia lipolytica
fermenter itself.
"""
from typing import Optional
import biosteam as bst
from biosteam.units.tank import Tank
from biosteam.units.design_tools.tank_design import mix_tank_purchase_cost_algorithms

__all__ = ('YarrowiaLipidFermenter', 'FermentationMediumTank')


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

    Sources
    -------
    target_pH:
        Gao et al. 2020, Biotechnol Biofuels -- initial pH 8 reported as
        optimal under high-VFA alkaline cultivation.
    product_yield_kg_per_kg_vfa_consumed:
        Gao et al. 2020, Biotechnol Biofuels -- Y_L/S at pH 8.
    conversion, biomass_yield_kg_per_kg_vfa_consumed,
    co2_yield_kg_per_kg_vfa_consumed, oxygen_kg_per_kg_vfa_consumed,
    residence_time_h:
        No literature source found in assumptions.yaml or elsewhere in this
        codebase (checked both biorefineries/sabre and sabre_ad). Current
        scoping assumptions for process modeling. Should be sensitivity-tested.
    pathway feasibility:
        Fontanille et al. 2012 (Yarrowia lipolytica assimilates VFAs and
        accumulates intracellular lipids); Pereira et al. 2023 (tested on
        food-waste-derived VFAs); Frontiers 2021 review (VFA-to-microbial-oil
        pathway concept).

    References
    ----------
    [1] Gao, R.; Li, Z.; Zhou, X.; Bao, W.; Cheng, S.; Zheng, L. Enhanced Lipid Production by Yarrowia Lipolytica Cultured with Synthetic and Waste-Derived High-Content Volatile Fatty Acids under Alkaline Conditions. Biotechnol Biofuels 2020, 13 (1), 3. https://doi.org/10.1186/s13068-019-1645-y.
    """

    V_max_default = 150

    def _init(
        self,
        vfa_IDs: Optional[list[str]] = None,
        product_ID: str = "MicrobialOil",
        conversion: float = 0.85,  # no literature source found; scoping assumption
        product_yield_kg_per_kg_vfa_consumed: float = 0.144,  # Gao et al. 2020, Y_L/S at pH 8
        biomass_yield_kg_per_kg_vfa_consumed: float = 0.40,  # scoping assumption
        co2_yield_kg_per_kg_vfa_consumed: float = 0.20,  # scoping assumption
        oxygen_kg_per_kg_vfa_consumed: float = 0.80,  # scoping assumption
        residence_time_h: float = 48.0,  # scoping assumption
        dT_hx_loop: float = 8.0,
        Q_O2_consumption: float = -460240.0,
        batch: bool = True,
        broth_density_kg_per_m3: float = 1000.0,
        target_pH: float = 8.0,  # Gao et al. 2020, optimal under high-VFA alkaline cultivation

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


class FermentationMediumTank(Tank):
    """
    Simple medium adjustment tank for pH control and nutrient dosing.
    Reuses BioSTEAM's own conventional stainless-steel mix-tank
    purchase-cost algorithm by inheriting from `bst.units.tank.Tank`
    directly (rather than reimplementing it) -- `_cost()` is Tank's own,
    unmodified. That means the bare-module factor (2.3x) and the
    stainless-steel material-factor normalization Tank._cost() applies
    (dividing the material factor out of the raw correlation cost, then
    reapplying it via `F_M`) both come from BioSTEAM directly. Previously
    this was a `bst.Unit` subclass with a hand-rolled `_cost()` that
    duplicated the cost correlation but hardcoded F_BM=1.0 -- silently
    skipping both factors and understating installed cost.

    Only `_run()` (dosing) and `_design()` (this unit's own
    residence-time/broth-density volume basis, not Tank's generic
    tau * F_vol_out / V_wf) are overridden; `_design()` writes into
    Tank's own `'Total volume'`/`'Residence time'` design_results keys so
    the inherited `_cost()` can read them directly.
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
        ammonia_dose_kg_per_m3: float = 0.0,
        phosphate_dose_kg_per_m3: float = 0.0,
        base_dose_kg_per_m3: float = 0.0,
        magnesium_sulfate_dose_kg_per_m3: float = 0.0,
        broth_density_kg_per_m3: float = 1000.0,
        target_pH: float = 8.0,
        residence_time_h: float = 0.5,
        mixing_kW_per_m3: float = 0.05,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.ammonia_dose_kg_per_m3 = float(ammonia_dose_kg_per_m3)
        self.phosphate_dose_kg_per_m3 = float(phosphate_dose_kg_per_m3)
        self.base_dose_kg_per_m3 = float(base_dose_kg_per_m3)
        self.magnesium_sulfate_dose_kg_per_m3 = float(magnesium_sulfate_dose_kg_per_m3)
        self.broth_density_kg_per_m3 = float(broth_density_kg_per_m3)
        self.target_pH = float(target_pH)
        self.residence_time_h = float(residence_time_h)
        # Tank._init() (run inside super().__init__() above) already set
        # self.kW_per_m3 = 0.0 (Tank's own default); overwrite with ours.
        self.kW_per_m3 = float(mixing_kW_per_m3)

    def _run(self):
        broth, ammonia, phosphate, base, mgso4 = self.ins
        out = self.outs[0]
        out.copy_like(broth)
        out.phase = "l"

        vol_m3ph = broth.F_mass / self.broth_density_kg_per_m3 if self.broth_density_kg_per_m3 > 0 else 0.0
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

    def _design(self):
        # Deliberately not Tank's own generic tau * F_vol_out / V_wf --
        # this unit sizes off the dosed broth's constant-density flow
        # basis instead. Writes into Tank's own design_results keys
        # ('Total volume', 'Residence time') so the inherited _cost()
        # (Tank._cost(), not overridden here) can read them directly.
        vol_m3h = float(self.design_results.get("Broth flow (m3/h)", 0.0))
        self.design_results["Residence time"] = self.residence_time_h
        self.design_results["Total volume"] = vol_m3h * self.residence_time_h

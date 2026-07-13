# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
VFA fermentation-to-microbial-oil system builder for the SaBRe flowsheets.
"""
from __future__ import annotations

import biosteam as bst
import flexsolve as flx

from biorefineries.sabre.units import YarrowiaLipidFermenter, OilExtractionPlaceholder
from biosteam.units.design_tools.tank_design import (
    mix_tank_purchase_cost_algorithms,
    compute_number_of_tanks_and_purchase_cost,
)

__all__ = ('create_vfa_fermentation_system', 'VFAMicrofilter', 'FermentationMediumTank')


class VFAMicrofilter(bst.Unit):
    """
    Split-based representation of a VFA-rich permeate step.
    Includes first-pass power draw and area-based membrane cost.
    """
    _N_ins = 1
    _N_outs = 2
    _F_BM_default = {"Microfilter": 1.0}

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
        membrane_cost_usd_per_m2: float = 200.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.vfa_IDs = tuple(vfa_IDs or [
            "AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"
        ])
        self.solids_IDs = tuple(solids_IDs or [
            "Ash", "Protein", "Lignin", "Glucan", "Xylan", "Mannan", "Galactan",
            "Arabinan", "Alginate", "Fucoidan", "Mannitol", "OtherSolids", "CellMass"
        ])
        self.vfa_to_permeate_frac = float(vfa_to_permeate_frac)
        self.water_to_permeate_frac = float(water_to_permeate_frac)
        self.solids_to_permeate_frac = float(solids_to_permeate_frac)
        self.dissolved_other_to_permeate_frac = float(dissolved_other_to_permeate_frac)
        self.broth_density_kg_per_m3 = float(broth_density_kg_per_m3)
        self.SEC_kWh_per_m3_feed = float(SEC_kWh_per_m3_feed)
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)
        self.membrane_cost_usd_per_m2 = float(membrane_cost_usd_per_m2)

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

    def _cost(self):
        A = float(self.design_results.get("Membrane area (m2)", 0.0))
        self.baseline_purchase_costs["Microfilter"] = max(A, 0.0) * self.membrane_cost_usd_per_m2

class FermentationMediumTank(bst.Unit):
    """
    Simple medium adjustment tank for pH control and nutrient dosing.
    """

    _N_ins = 5
    _N_outs = 1
    _F_BM_default = {"Medium tank": 1.0}

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
        self.mixing_kW_per_m3 = float(mixing_kW_per_m3)

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
        vol_m3h = float(self.design_results.get("Broth flow (m3/h)", 0.0))
        tank_vol = vol_m3h * self.residence_time_h
        self.design_results["Tank residence time (h)"] = self.residence_time_h
        self.design_results["Tank volume (m3)"] = tank_vol
        self.power_utility.consumption = self.mixing_kW_per_m3 * tank_vol

    def _cost(self):
        V_total = max(float(self.design_results.get("Tank volume (m3)", 0.0)), 0.0)
        if V_total <= 0.0:
            self.baseline_purchase_costs["Medium tank"] = 0.0
            return

        alg = mix_tank_purchase_cost_algorithms["Conventional"]
        N, Cp_each = compute_number_of_tanks_and_purchase_cost(V_total, alg)
        V_each = V_total / N

        self.design_results["Number of tanks"] = N
        self.design_results["Tank volume per tank (m3)"] = V_each
        self.baseline_purchase_costs["Medium tank"] = N * Cp_each


def create_vfa_fermentation_system(
    vfa_broth,
    *,
    product_ID: str = "MicrobialOil",
    vfa_IDs: list[str] | None = None,
    conversion: float = 0.85,
    product_yield_kg_per_kg_vfa_consumed: float = 0.144,
    biomass_yield_kg_per_kg_vfa_consumed: float = 0.40,
    co2_yield_kg_per_kg_vfa_consumed: float = 0.20,
    oxygen_kg_per_kg_vfa_consumed: float = 0.80,
    residence_time_h: float = 48.0,
    broth_density_kg_per_m3: float = 1000.0,
    target_pH: float = 8.0,
    ammonia_dose_kg_per_m3: float = 0.0,
    phosphate_dose_kg_per_m3: float = 0.0,
    base_dose_kg_per_m3: float = 0.0,
    magnesium_sulfate_dose_kg_per_m3: float = 0.0,
    seed_water_kgph: float = 0.0,
    seed_cellmass_kgph: float = 0.0,

    vfa_to_permeate_frac: float = 0.98,
    water_to_permeate_frac: float = 0.97,
    solids_to_permeate_frac: float = 0.05,
    dissolved_other_to_permeate_frac: float = 0.90,
    microfilter_SEC_kWh_per_m3_feed: float = 0.08,
    microfilter_design_flux_L_m2_h: float = 20.0,
    microfilter_membrane_cost_usd_per_m2: float = 200.0,

    medium_tank_residence_time_h: float = 0.5,
    medium_tank_mixing_kW_per_m3: float = 0.05,

    target_oil_and_solids_content: float = 60.0,
    target_wastewater_concentration: float = 60.0,
    backend_oil_recovery: float = 0.99,
    backend_oil_water_split: float = 0.0001,

    recycle_total_fraction: float = 0.10,
    recycle_cellmass_wt_frac: float = 0.10,

    oil_extraction_ref_dry_biomass_tph: float = 10.0,
    oil_extraction_ref_installed_cost_usd: float = 7_848_000.0,
    oil_extraction_homogenization_kWh_per_kg: float = 0.203,
    oil_extraction_scale_exponent: float = 0.6,
):

    chem_ids = set(bst.settings.thermo.chemicals.IDs)
    for req in (product_ID, "CarbonDioxide", "CellMass", "Water"):
        if req not in chem_ids:
            raise RuntimeError(f"Missing required chemical in thermo: {req}")

    if vfa_IDs is None:
        vfa_IDs = ["AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"]

    ammonia = bst.Stream("fermentation_ammonia")
    phosphate = bst.Stream("fermentation_phosphate")
    base = bst.Stream("fermentation_base")
    mgso4 = bst.Stream("fermentation_mgso4")
    fresh_seed = bst.Stream(
        "fresh_seed_makeup",
        Water=seed_water_kgph,
        CellMass=seed_cellmass_kgph,
        units="kg/hr",
    )
    recycle_biomass = bst.Stream("recycle_biomass")
    M602 = bst.Mixer("M602", ins=(fresh_seed, recycle_biomass), outs=("seed_to_fermenter",))

    # -------------------------------------------------
    # Upstream conditioning
    # -------------------------------------------------
    MF = VFAMicrofilter(
        "MF",
        ins=vfa_broth,
        outs=("vfa_permeate", "vfa_retentate"),
        vfa_IDs=vfa_IDs,
        vfa_to_permeate_frac=vfa_to_permeate_frac,
        water_to_permeate_frac=water_to_permeate_frac,
        solids_to_permeate_frac=solids_to_permeate_frac,
        dissolved_other_to_permeate_frac=dissolved_other_to_permeate_frac,
        broth_density_kg_per_m3=broth_density_kg_per_m3,
        SEC_kWh_per_m3_feed=microfilter_SEC_kWh_per_m3_feed,
        design_flux_L_m2_h=microfilter_design_flux_L_m2_h,
        membrane_cost_usd_per_m2=microfilter_membrane_cost_usd_per_m2,
    )

    T601 = FermentationMediumTank(
        "T601",
        ins=(MF - 0, ammonia, phosphate, base, mgso4),
        outs=("conditioned_vfa_broth",),
        ammonia_dose_kg_per_m3=ammonia_dose_kg_per_m3,
        phosphate_dose_kg_per_m3=phosphate_dose_kg_per_m3,
        base_dose_kg_per_m3=base_dose_kg_per_m3,
        magnesium_sulfate_dose_kg_per_m3=magnesium_sulfate_dose_kg_per_m3,
        broth_density_kg_per_m3=broth_density_kg_per_m3,
        target_pH=target_pH,
        residence_time_h=medium_tank_residence_time_h,
        mixing_kW_per_m3=medium_tank_mixing_kW_per_m3,
    )

    # -------------------------------------------------
    # Aerated fermenter (Yarrowia lipolytica)
    # -------------------------------------------------
    R601 = YarrowiaLipidFermenter(
        "R601",
        ins=(T601 - 0, M602 - 0),
        outs=("fermentation_vent", "fermentation_broth"),
        vfa_IDs=vfa_IDs,
        product_ID=product_ID,
        conversion=conversion,
        product_yield_kg_per_kg_vfa_consumed=product_yield_kg_per_kg_vfa_consumed,
        biomass_yield_kg_per_kg_vfa_consumed=biomass_yield_kg_per_kg_vfa_consumed,
        co2_yield_kg_per_kg_vfa_consumed=co2_yield_kg_per_kg_vfa_consumed,
        oxygen_kg_per_kg_vfa_consumed=oxygen_kg_per_kg_vfa_consumed,
        residence_time_h=residence_time_h,
        broth_density_kg_per_m3=broth_density_kg_per_m3,
        target_pH=target_pH,
        V_wf=0.8,
        V_max=150.0,
        kW_per_m3=0.06,
        tau_0=3.0,
    )

    # -------------------------------------------------
    # Post-fermentation oil separation
    # -------------------------------------------------
    V605 = bst.MixTank("V605", ins=R601 - 1, outs=("mixed_fermentation_broth",))
    P606 = bst.Pump("P606", ins=V605 - 0, outs=("pumped_fermentation_broth",))

    Ev607 = bst.MultiEffectEvaporator(
        "Ev607",
        ins=P606 - 0,
        outs=("fermentation_concentrate", "evaporator_vapor"),
        P=(101325, 69682, 47057, 30953),
        V=0.90,
        V_definition="First-effect",
        thermo=(R601.outs[1].thermo.ideal()),
        flash=False,
    )
    Ev607.target_oil_and_solids_content = target_oil_and_solids_content
    Ev607.remove_evaporators = True

    P_original = tuple(Ev607.P)
    Pstart = P_original[0]
    Plast = P_original[-1]
    N = len(P_original)

    def concentration_objective(V):
        Ev607.V = V
        Ev607.run()
        effluent = Ev607.outs[0]
        total = effluent.F_mass
        if total <= 0:
            return 0.0
        water = effluent.imass["Water"]
        nonwater_conc = 1000.0 * (1.0 - water / total)
        return Ev607.target_oil_and_solids_content - nonwater_conc

    @Ev607.add_specification(run=False)
    def adjust_evaporation():
        V_last = Ev607.V
        x0 = 0.0
        x1 = 0.5
        Ev607.P = P_original
        Ev607._reload_components = True

        y0 = concentration_objective(x0)
        if y0 <= 0.0:
            Ev607.V = x0
            return
        else:
            Ev607._load_components()
            for i in range(1, N):
                if concentration_objective(1e-6) < 0.0:
                    Ev607.P = tuple(__import__("numpy").linspace(Pstart, Plast, N - 1))
                    Ev607._reload_components = True
                else:
                    break
            y1 = concentration_objective(x1)
            Ev607.V = flx.IQ_interpolation(
                concentration_objective,
                x0, x1, y0, y1,
                x=V_last,
                ytol=1e-5,
                xtol=1e-6,
            )

    P607 = bst.Pump("P607", ins=Ev607 - 0, outs=("pumped_concentrate",), P=101325.0)

    # -------------------------------------------------
    # Cell disruption + oil extraction placeholder
    # Slots between P607 and C603_2.
    # Pass-through unit: carries capital + electricity costs for
    # high-pressure homogenization and lipid extraction.
    # Capital anchor: NREL/TP-5100-55431 (2012), Davis et al.
    # -------------------------------------------------
    OE = OilExtractionPlaceholder(
        "OE",
        ins=P607 - 0,
        outs=("extracted_broth",),
        product_ID=product_ID,
        cellmass_ID="CellMass",
        homogenization_kWh_per_kg_dry_biomass=oil_extraction_homogenization_kWh_per_kg,
        ref_dry_biomass_tph=oil_extraction_ref_dry_biomass_tph,
        ref_installed_cost_usd=oil_extraction_ref_installed_cost_usd,
        scale_exponent=oil_extraction_scale_exponent,
        F_BM=1.0,
    )

    C603_2 = bst.LiquidsSplitCentrifuge(
        "C603_2",
        ins=OE - 0,   # now takes from OE, not P607
        outs=("backend_oil", "cellmass_plus_aqueous"),
        split={product_ID: backend_oil_recovery, "Water": backend_oil_water_split},
    )

    S602 = bst.MockSplitter(
        "S602",
        ins=C603_2 - 1,
        outs=(recycle_biomass, "residual_purge"),
    )
    S602.recycle_total_fraction = recycle_total_fraction
    S602.recycle_cellmass_wt_frac = recycle_cellmass_wt_frac

    @S602.add_specification(run=True)
    def adjust_biomass_recycle():
        feed = S602.ins[0]
        recycle, purge = S602.outs

        recycle.empty()
        purge.copy_like(feed)

        cellmass_available = float(feed.imass["CellMass"]) if "CellMass" in feed.chemicals else 0.0
        water_available = float(feed.imass["Water"]) if "Water" in feed.chemicals else 0.0

        if feed.F_mass <= 0 or cellmass_available <= 0:
            recycle.empty()
            purge.copy_like(feed)
            return

        target_recycle_mass = S602.recycle_total_fraction * feed.F_mass
        target_cellmass = target_recycle_mass * S602.recycle_cellmass_wt_frac
        cellmass_recycle = min(cellmass_available, target_cellmass)
        water_recycle = cellmass_recycle * (1.0 / S602.recycle_cellmass_wt_frac - 1.0)
        water_recycle = min(water_available, water_recycle)

        recycle.imass["CellMass"] = cellmass_recycle
        recycle.imass["Water"] = water_recycle
        purge.imass["CellMass"] -= cellmass_recycle
        purge.imass["Water"] -= water_recycle
        purge.mol.remove_negatives()
        recycle.T = purge.T = feed.T

    S601 = bst.Splitter(
        "S601",
        ins=Ev607 - 1,
        outs=("condensate_to_wastewater", "evaporator_condensate"),
        split=0.5,
    )

    M601 = bst.Mixer(
        "M601",
        ins=(S601 - 0, S602 - 1),
        outs=("fermentation_wastewater",),
    )

    M601.target_wastewater_concentration = target_wastewater_concentration

    @M601.add_specification(run=True, impacted_units=[S601])
    def adjust_wastewater_concentration():
        concentrated_wastewater = M601.ins[1]
        waste = concentrated_wastewater.F_mass - concentrated_wastewater.imass["Water"]
        if concentrated_wastewater.F_vol <= 0 or waste <= 0:
            S601.split[:] = 0.0
            return

        current_concentration = waste / concentrated_wastewater.F_vol
        required_water = (
            (1.0 / M601.target_wastewater_concentration) - (1.0 / current_concentration)
        ) * waste * 1000.0

        F_mass = S601.ins[0].F_mass
        if F_mass > 0:
            split = required_water / F_mass
            split = min(max(split, 0.0), 1.0)
            S601.split[:] = split

    # OE added to system path between P607 and C603_2
    sys = bst.System(
        "VFA_FER_sys",
        path=(MF, T601, M602, R601, V605, P606, Ev607, P607, OE, C603_2, S602, S601, M601)
    )

    key_streams = {
        "vfa_broth": vfa_broth,
        "vfa_permeate": MF.outs[0],
        "vfa_retentate": MF.outs[1],
        "conditioned_vfa_broth": T601.outs[0],
        "fresh_seed_makeup": fresh_seed,
        "recycle_biomass": recycle_biomass,
        "seed": M602.outs[0],
        "vent": R601.outs[0],
        "fermentation_broth": R601.outs[1],
        "extracted_broth": OE.outs[0],
        "backend_oil": C603_2.outs[0],
        "residual_slurry": C603_2.outs[1],
        "residual_purge": S602.outs[1],
        "evaporator_condensate": S601.outs[1],
        "fermentation_wastewater": M601.outs[0],
        "ammonia": ammonia,
        "phosphate": phosphate,
        "base": base,
        "mgso4": mgso4,
    }

    units = {
        "microfilter": MF,
        "medium_tank": T601,
        "seed_mixer": M602,
        "fermenter": R601,
        "post_mix_tank": V605,
        "post_feed_pump": P606,
        "evaporator": Ev607,
        "backend_oil_pump": P607,
        "oil_extraction": OE,
        "oil_centrifuge": C603_2,
        "biomass_recycle_splitter": S602,
        "condensate_splitter": S601,
        "wastewater_mixer": M601,
        # legacy aliases
        "centrifuge": C603_2,
        "lipid_recovery": C603_2,
    }

    return sys, key_streams, units

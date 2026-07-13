# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Acidogenic (arrested) AD unit operation for the VFA platform.
"""
from __future__ import annotations

import math
from typing import Dict, Iterable, Optional

import biosteam as bst

from biorefineries.sabre.units._ad import GAL_TO_M3, interp_capex

__all__ = ('AcidogenicDigester',)


class AcidogenicDigester(bst.Unit):
    """
    Acidogenic / arrested AD unit for VFA platform modeling.

    Outputs:
    0) offgas
    1) acidogenic_broth

    """

    _N_ins = 1
    _N_outs = 2  # offgas, acidogenic_broth

    F_BM = {"Acidogenic digester": 1.0}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        vs_destruction: float = 0.50,
        vfa_kg_per_kg_vs: float = 0.80,
        vfa_split: Optional[Dict[str, float]] = None,
        digestible_IDs: Optional[Iterable[str]] = None,
        produce_offgas_co2: bool = True,
        hrt_days: float = 15.0,
        slurry_density_kg_per_m3: float = 1000.0,
        headspace_frac: float = 0.2,
        max_single_digester_volume_MG: float = 1.5,
        mixing_W_per_m3: float = 5.0,
        influent_temperature_K: float = 298.15,
        target_temperature_K: float = 308.15,
        cp_kJ_per_kgK: float = 4.18,
        enable_heat_shock: bool = False,
        hs_target_temperature_K: float = 338.15,
        hs_events_per_day: float = 1.0 / 7.0,
        hs_heated_fraction_of_liquid: float = 0.10,
        hs_duration_min: float = 15.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.vs_destruction = float(vs_destruction)
        self.vfa_kg_per_kg_vs= float(vfa_kg_per_kg_vs)
        self.vfa_split = dict(vfa_split) if vfa_split is not None else None
        self.digestible_IDs = tuple(digestible_IDs) if digestible_IDs is not None else (
            "Glucan",
            "Xylan",
            "Mannan",
            "Galactan",
            "Arabinan",
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
            "Cellulose",
        )
        self.produce_offgas_co2 = bool(produce_offgas_co2)

        self.hrt_days = float(hrt_days)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)
        self.headspace_frac = float(headspace_frac)
        self.max_single_digester_volume_m3 = float(max_single_digester_volume_MG) * 1e6 * GAL_TO_M3

        self.mixing_W_per_m3 = float(mixing_W_per_m3)
        self.influent_temperature_K = float(influent_temperature_K)
        self.target_temperature_K = float(target_temperature_K)
        self.cp_kJ_per_kgK = float(cp_kJ_per_kgK)

        self.enable_heat_shock = bool(enable_heat_shock)
        self.hs_target_temperature_K = float(hs_target_temperature_K)
        self.hs_events_per_day = float(hs_events_per_day)
        self.hs_heated_fraction_of_liquid = float(hs_heated_fraction_of_liquid)
        self.hs_duration_min = float(hs_duration_min)

        self.F_BM = dict(self.F_BM)

    def _resolve_vfa_split(self):
        chems = bst.settings.thermo.chemicals
        ids = set(chems.IDs)

        if self.vfa_split is not None:
            split = dict(self.vfa_split)
            missing = [k for k in split.keys() if k not in ids]
            if missing:
                raise RuntimeError(
                    f"VFA split includes chemicals not in thermo: {missing}. "
                    "Use thermo IDs that actually exist in your project."
                )
            s = sum(split.values())
            if s <= 0:
                raise RuntimeError("vfa_split sums to 0; provide positive fractions.")
            return {k: v / s for k, v in split.items()}

        common = ["AceticAcid", "PropionicAcid", "ButyricAcid", "ValericAcid", "HexanoicAcid"]
        if all(c in ids for c in common):
            return {
                "AceticAcid": 0.40,
                "PropionicAcid": 0.10,
                "ButyricAcid": 0.30,
                "ValericAcid": 0.05,
                "HexanoicAcid": 0.15,
            }

        if "VFA" in ids:
            return {"VFA": 1.0}

        raise RuntimeError(
            "No valid VFA representation found. Add the acid IDs to thermo or add a pseudo-component named 'VFA'."
        )

    def _available_digestible_pool(self, broth) -> Dict[str, float]:
        ids = set(bst.settings.thermo.chemicals.IDs)
        avail: Dict[str, float] = {}
        for cid in self.digestible_IDs:
            if cid not in ids:
                continue
            m = float(broth.imass[cid])
            if m > 0:
                avail[cid] = m
        return avail

    def _run(self):
        feed = self.ins[0]
        offgas, broth = self.outs

        offgas.empty()
        broth.copy_like(feed)

        offgas.phase = "g"
        broth.phase = "l"

        water = float(broth.imass["Water"]) if "Water" in broth.chemicals else 0.0
        ash = float(broth.imass["Ash"]) if "Ash" in broth.chemicals else 0.0
        TS = max(broth.F_mass - water, 0.0)
        VS = max(TS - ash, 0.0)
        VS_destroyed_target = max(0.0, self.vs_destruction) * VS
        if VS_destroyed_target <= 1e-12:
            return

        avail = self._available_digestible_pool(broth)
        pool = sum(avail.values())
        if pool <= 1e-12:
            return

        remove = min(VS_destroyed_target, pool)
        for cid, m in avail.items():
            take = remove * (m / pool)
            broth.imass[cid] -= take

        split = self._resolve_vfa_split()
        vfa_total = max(0.0, self.vfa_kg_per_kg_vs) * remove
        vfa_total = min(vfa_total, remove)

        for chem_id, frac in split.items():
            broth.imass[chem_id] += vfa_total * frac

        residual_destroyed = max(0.0, remove - vfa_total)
        if residual_destroyed > 0.0:
            ids = set(bst.settings.thermo.chemicals.IDs)
            if self.produce_offgas_co2 and "CarbonDioxide" in ids:
                offgas.imass["CarbonDioxide"] += residual_destroyed
            elif "Water" in ids:
                broth.imass["Water"] += residual_destroyed

        self.design_results["VS in (kg/h)"] = VS
        self.design_results["VS destroyed target (kg/h)"] = VS_destroyed_target
        self.design_results["Digestible pool (kg/h)"] = pool
        self.design_results["Destroyed digestible mass (kg/h)"] = remove
        self.design_results["Total VFA produced (kg/h)"] = vfa_total
        self.design_results["Residual destroyed mass to closure (kg/h)"] = residual_destroyed
        for chem_id, frac in split.items():
            self.design_results[f"{chem_id} produced (kg/h)"] = vfa_total * frac

    def _design(self):
        feed = self.ins[0]
        slurry_m3_per_hr = feed.F_mass / self.slurry_density_kg_per_m3
        V_liquid = slurry_m3_per_hr * 24.0 * self.hrt_days

        hf = min(max(self.headspace_frac, 0.0), 0.95)
        V_total = V_liquid / (1.0 - hf)

        V_max = self.max_single_digester_volume_m3
        N = max(1, math.ceil(V_total / V_max))
        V_each = V_total / N

        self.design_results["Slurry flow (m3/hr)"] = slurry_m3_per_hr
        self.design_results["HRT (days)"] = self.hrt_days
        self.design_results["Headspace frac"] = hf
        self.design_results["Total digester volume (m3)"] = V_total
        self.design_results["Number of digesters"] = N
        self.design_results["Digester volume each (m3)"] = V_each

        mixing_kW = (self.mixing_W_per_m3 * V_liquid) / 1000.0
        self.design_results["Mixing power (kW)"] = mixing_kW
        self.power_utility.consumption = mixing_kW
        self.power_utility.production = 0.0

        T_in = self.influent_temperature_K
        T_base = self.target_temperature_K
        dT = max(0.0, T_base - T_in)
        m_dot_kgph = feed.F_mass
        Q_base_kJph = m_dot_kgph * self.cp_kJ_per_kgK * dT

        self.design_results["Influent T (K)"] = T_in
        self.design_results["Base target T (K)"] = T_base
        self.design_results["Base heating duty (kJ/h)"] = Q_base_kJph

        Q_hs_kJph = 0.0
        if self.enable_heat_shock:
            frac = min(max(self.hs_heated_fraction_of_liquid, 0.0), 1.0)
            events_per_day = max(0.0, self.hs_events_per_day)
            dT_hs = max(0.0, self.hs_target_temperature_K - T_base)
            m_liq_kg = V_liquid * self.slurry_density_kg_per_m3
            Q_event_kJ = (m_liq_kg * frac) * self.cp_kJ_per_kgK * dT_hs
            Q_day_kJ = Q_event_kJ * events_per_day
            Q_hs_kJph = Q_day_kJ / 24.0

            self.design_results["HS enabled"] = True
            self.design_results["HS target T (K)"] = self.hs_target_temperature_K
            self.design_results["HS events/day"] = events_per_day
            self.design_results["HS heated fraction"] = frac
            self.design_results["HS duration (min)"] = self.hs_duration_min
            self.design_results["HS avg duty (kJ/h)"] = Q_hs_kJph
        else:
            self.design_results["HS enabled"] = False

        Q_total_kJph = Q_base_kJph + Q_hs_kJph
        self.design_results["Total heating duty (kJ/h)"] = Q_total_kJph
        if Q_total_kJph > 0:
            try:
                self.add_heat_utility(
                    Q_total_kJph,
                    T_in,
                    max(T_base, self.hs_target_temperature_K if self.enable_heat_shock else T_base),
                )
            except TypeError:
                self.add_heat_utility(Q_total_kJph, T_in)

    def _cost(self):
        V_each = self.design_results["Digester volume each (m3)"]
        N = int(self.design_results["Number of digesters"])
        self.F_BM["Acidogenic digester"] = 1.0

        C_each = interp_capex(V_each)
        C_total = N * C_each

        self.design_results["ADBC_capex_each_$"] = C_each
        self.baseline_purchase_costs["Acidogenic digester"] = C_total

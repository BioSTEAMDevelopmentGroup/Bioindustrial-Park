# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Anaerobic digestion unit operation (AD)

Purpose:
- Convert a biodegradable fraction of Sargassum organics into raw biogas
  using component-specific biodegradability factors
- Produce a residual digestate stream containing undegraded material
- Size the digester using HRT and headspace
- Enforce a maximum single-digester size and model parallel digesters
- Estimate CAPEX using ADBC-based interpolation

Modeling notes:
- Lumped methanogenic AD model
- Methane production is based on:
    feed VS * ch4_kg_per_kg_vs_fed
- Different components can have different biodegradability factors
- Biodegradability + vs_destruction determine digestate solids reduction
- Raw biogas composition is imposed from target mole fractions
"""

from __future__ import annotations

import math

import biosteam as bst

__all__ = ('AnaerobicDigester', 'GAL_TO_M3', 'NM3_PER_KMOL', 'interp_capex')

# Conversion factors and ADBC data
GAL_TO_M3 = 0.003785411784  # US gallons -> m3
NM3_PER_KMOL = 22.414

# ADBC-based digester installed cost anchor points (m3, USD)
ADBC_VOL_M3 = [878, 1755, 2633, 3510, 5265, 8775]
ADBC_CAPEX = [1720964, 1750201, 1779439, 1808676, 1867151, 1984101]


def interp_capex(volume_m3: float) -> float:
    """Linear interpolation/extrapolation of ADBC digester CAPEX"""
    x = ADBC_VOL_M3
    y = ADBC_CAPEX

    if volume_m3 <= x[0]:
        m = (y[1] - y[0]) / (x[1] - x[0])
        return y[0] + m * (volume_m3 - x[0])

    if volume_m3 >= x[-1]:
        m = (y[-1] - y[-2]) / (x[-1] - x[-2])
        return y[-1] + m * (volume_m3 - x[-1])

    for i in range(len(x) - 1):
        if x[i] <= volume_m3 <= x[i + 1]:
            m = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            return y[i] + m * (volume_m3 - x[i])

    raise RuntimeError("CAPEX interpolation failed")


class AnaerobicDigester(bst.Unit):
    """
    Lumped mesophilic anaerobic digester with parallel-train sizing and
    ADBC-interpolated capital cost.

    Inputs:
        ins[0]: pretreated/diluted feed slurry

    Outputs:
        outs[0]: raw biogas
        outs[1]: digestate

    Sources
    -------
    vs_destruction (0.20):
        assumptions.yaml `ad_performance.vs_destruction`, labeled "Global CSTR
        vs_destruction -- uniform across all pretreatment cases." No named
        literature source; yaml's own note: "Conservative lower bound of
        20-70% range for full-scale mesophilic CSTRs. No continuous CSTR data
        for pretreated pelagic Sargassum exists in literature." Scoping
        assumption, not a citation.
    ch4_kg_per_kg_vs_fed (0.100):
        This is the `press_mill_only` baseline pretreatment case value (yaml
        `ad_pretreatment_cases.press_mill_only.ad_effects.ch4_kg_per_kg_vs_fed`),
        not a universal constant -- the 4 other pretreatment cases in yaml use
        different values (0.123 enzymatic, 0.165 peroxide, 0.277 combined_PE,
        0.261 combined_PTE), each cited to a different row of the same table.
        Source: Chikani-Cabrera et al. 2022, Table 4 -- "224.19 +/- 9.45 L
        CH4/kg VS," scaled to 0.100 kg/kg VS for continuous CSTR (yaml key
        `chikani_cabrera_2022_physical_C`).
    raw_biogas_molfrac (Methane 0.55, CarbonDioxide 0.42, HydrogenSulfide 0.03):
        assumptions.yaml `ad_performance.raw_biogas_molfrac`, labeled "Global
        fallback biogas composition (overridden per pretreatment case below)."
        No named source -- scoping default, overridden per-case in real runs.
    headspace_frac (0.25):
        assumptions.yaml `ad.gas_storage_frac_of_total_volume`. No named
        source.
    maintenance_usd_per_m3_yr (10.0):
        assumptions.yaml `ad_costing.maintenance_usd_per_m3_yr`, sourced to
        the ADBC spreadsheet tool (`ADBCv2-48 (2).xlsm`, sheet
        "AnaerobicDigester", cell AB18).
    digestible_IDs, biodegradability:
        assumptions.yaml `ad_performance.digestible_IDs` /
        `ad_performance.biodegradability`. No named source for the individual
        factors beyond the yaml values themselves.
    hrt_days, slurry_density_kg_per_m3, max_single_digester_volume_MG,
    mixing_W_per_m3, influent_temperature_K, target_temperature_K,
    cp_kJ_per_kgK:
        assumptions.yaml `ad` section, matching values already used at the
        real call site (`systems/_ad_biogas_system.py`). No named literature
        source for these beyond the yaml values themselves.
    ADBC_VOL_M3 / ADBC_CAPEX interpolation table (used by `interp_capex`,
    the actual installed-cost calculation):
        UNVERIFIED PROVENANCE. assumptions.yaml has no entry for this
        six-point table. A *separate*, structurally unused single-point
        anchor (`ad_costing.base_volume_m3`/`base_capex_usd`, stored on this
        unit's `base_volume_m3`/`base_capex_usd` __init__ params but never
        read by `_design`/`_cost`) is cited to the same-sounding "ADBC"
        spreadsheet tool (`ADBCv2-48 (2).xlsm`, sheet "AnaerobicDigester").
        Whether the six-point table below was pulled from the same
        spreadsheet run has not been verified by anyone in this conversion --
        do not repeat "ADBC spreadsheet" as a confirmed source for
        ADBC_VOL_M3/ADBC_CAPEX without checking the actual file.
    """

    _N_ins = 1
    _N_outs = 2  # biogas, digestate

    F_BM = {"Anaerobic digester": 1.0}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        # Performance
        vs_destruction=0.20,
        ch4_kg_per_kg_vs_fed=0.100,
        raw_biogas_molfrac=None,
        digestible_IDs=None,
        biodegradability=None,
        # Sizing
        hrt_days=25.0,
        slurry_density_kg_per_m3=1000.0,
        headspace_frac=0.25,
        max_single_digester_volume_MG=1.5,
        # Costing anchors
        base_volume_m3=None,
        base_capex_usd=None,
        maintenance_usd_per_m3_yr=10.0,
        # Utilities
        mixing_W_per_m3=5.0,
        influent_temperature_K=298.15,
        target_temperature_K=308.15,
        cp_kJ_per_kgK=4.18,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.vs_destruction = float(vs_destruction)
        self.ch4_kg_per_kg_vs_fed = float(ch4_kg_per_kg_vs_fed)

        self.raw_biogas_molfrac = raw_biogas_molfrac or {
            "Methane": 0.55,
            "CarbonDioxide": 0.42,
            "HydrogenSulfide": 0.03,
        }

        self.digestible_IDs = tuple(digestible_IDs) if digestible_IDs is not None else (
            "Glucan",
            "Arabinan",
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
        )

        self.biodegradability = biodegradability or {
            "Mannitol": 0.95,
            "Glucan": 0.75,
            "Arabinan": 0.70,
            "Protein": 0.65,
            "Alginate": 0.45,
            "Fucoidan": 0.35,
            "OtherSolids": 0.30,
        }

        self.mixing_W_per_m3 = float(mixing_W_per_m3)
        self.influent_temperature_K = float(influent_temperature_K)
        self.target_temperature_K = float(target_temperature_K)
        self.cp_kJ_per_kgK = float(cp_kJ_per_kgK)

        self.hrt_days = float(hrt_days)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)
        self.headspace_frac = float(headspace_frac)
        self.max_single_digester_volume_m3 = (
            float(max_single_digester_volume_MG) * 1e6 * GAL_TO_M3
        )

        self.base_volume_m3 = base_volume_m3
        self.base_capex_usd = base_capex_usd
        self.maintenance_usd_per_m3_yr = maintenance_usd_per_m3_yr

        self.F_BM = dict(self.F_BM)

    def _available_digestible_pool(self, stream):
        thermo_ids = set(bst.settings.thermo.chemicals.IDs)
        avail = {}

        for cid in self.digestible_IDs:
            if cid not in thermo_ids:
                continue

            mass = float(stream.imass[cid])
            if mass <= 0:
                continue

            factor = float(self.biodegradability.get(cid, 0.0))
            factor = min(max(factor, 0.0), 1.0)

            avail[cid] = {
                "mass": mass,
                "biodegradable_mass": factor * mass,
                "factor": factor,
            }

        return avail

    def _run(self):
        feed = self.ins[0]
        biogas, digestate = self.outs

        biogas.empty()
        digestate.copy_like(feed)

        biogas.phase = "g"
        digestate.phase = "l"

        chems = bst.settings.thermo.chemicals
        ids = set(chems.IDs)

        water_in = float(feed.imass["Water"]) if "Water" in ids else 0.0
        ash_in = float(feed.imass["Ash"]) if "Ash" in ids else 0.0
        TS_in = max(float(feed.F_mass) - water_in, 0.0)
        VS_in = max(TS_in - ash_in, 0.0)

        avail = self._available_digestible_pool(digestate)
        if not avail:
            return

        biodegradable_pool = sum(v["biodegradable_mass"] for v in avail.values())
        if biodegradable_pool <= 1e-12:
            return

        biodegradable_destroyed_target = self.vs_destruction * biodegradable_pool
        biodegradable_destroyed_actual = min(
            biodegradable_destroyed_target, biodegradable_pool
        )

        if biodegradable_destroyed_actual <= 1e-12:
            return

        for cid, info in avail.items():
            biodegradable_mass_i = info["biodegradable_mass"]
            if biodegradable_mass_i <= 0:
                continue

            take = biodegradable_destroyed_actual * (
                biodegradable_mass_i / biodegradable_pool
            )
            take = min(take, float(digestate.imass[cid]))
            digestate.imass[cid] -= take

        required = {"Methane", "CarbonDioxide"}
        missing = required - ids
        if missing:
            raise RuntimeError(
                f"AD unit missing required gas components in thermo: {missing}"
            )

        y = self.raw_biogas_molfrac
        y_ch4 = float(y.get("Methane", 0.0))
        y_co2 = float(y.get("CarbonDioxide", 0.0))
        y_h2s = float(y.get("HydrogenSulfide", 0.0))

        y_sum = y_ch4 + y_co2 + y_h2s
        if abs(y_sum - 1.0) > 1e-6:
            raise ValueError(
                f"Raw biogas mole fractions must sum to 1. Got {y_sum:.6f}"
            )
        if y_ch4 <= 0:
            raise ValueError("Methane mole fraction must be > 0.")

        ch4 = chems["Methane"]
        ch4_mass = self.ch4_kg_per_kg_vs_fed * VS_in
        n_ch4 = ch4_mass / ch4.MW

        if n_ch4 <= 1e-12:
            return

        n_total = n_ch4 / y_ch4

        biogas.imol["Methane"] = n_ch4
        biogas.imol["CarbonDioxide"] = y_co2 * n_total

        if "HydrogenSulfide" in ids and y_h2s > 0:
            biogas.imol["HydrogenSulfide"] = y_h2s * n_total

        # --- simple overall mass closure ---
        # The current model removes destroyed biodegradable solids from the digestate
        # and separately generates biogas from imposed methane yield/composition.
        # To prevent "more gas, same digestate" behavior in reporting, adjust the
        # digestate water inventory so total output mass closes.
        gas_mass = float(biogas.F_mass)
        water_adjust = biodegradable_destroyed_actual - gas_mass

        if "Water" in ids and abs(water_adjust) > 1e-9:
            current_water = float(digestate.imass["Water"])
            new_water = current_water + water_adjust
            if new_water < 0:
                raise RuntimeError(
                    f"AD water balance became negative for unit {self.ID}. "
                    f"Current water={current_water:.3f}, adjustment={water_adjust:.3f}."
                )
            digestate.imass["Water"] = new_water

        self._TS_in_kgph = TS_in
        self._VS_in_kgph = VS_in
        self._biodegradable_pool_kgph = biodegradable_pool
        self._biodegradable_destroyed_kgph = biodegradable_destroyed_actual
        self._methane_kgph = ch4_mass
        self._co2_kgph = float(biogas.imass["CarbonDioxide"]) if "CarbonDioxide" in ids else 0.0
        self._h2s_kgph = float(biogas.imass["HydrogenSulfide"]) if "HydrogenSulfide" in ids else 0.0
        self._biogas_mass_kgph = gas_mass
        self._biogas_Nm3ph = float(biogas.F_mol) * NM3_PER_KMOL
        self._digestate_mass_kgph = float(digestate.F_mass)

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
        self.design_results["Headspace frac"] = self.headspace_frac
        self.design_results["Total digester volume (m3)"] = V_total
        self.design_results["Max single digester (m3)"] = V_max
        self.design_results["Number of digesters"] = N
        self.design_results["Digester volume each (m3)"] = V_each

        mixing_kW = (self.mixing_W_per_m3 * V_liquid) / 1000.0

        T_in = self.influent_temperature_K
        T_target = self.target_temperature_K
        dT = max(0.0, T_target - T_in)

        m_dot_kgph = feed.F_mass
        Q_kJph = m_dot_kgph * self.cp_kJ_per_kgK * dT

        self.design_results["Mixing power (kW)"] = mixing_kW
        self.design_results["Influent T (K)"] = T_in
        self.design_results["Target T (K)"] = T_target
        self.design_results["Heating duty (kJ/h)"] = Q_kJph

        TS_in = getattr(self, "_TS_in_kgph", 0.0)
        VS_in = getattr(self, "_VS_in_kgph", 0.0)
        biodegradable_pool = getattr(self, "_biodegradable_pool_kgph", 0.0)
        biodegradable_destroyed = getattr(self, "_biodegradable_destroyed_kgph", 0.0)
        methane_kgph = getattr(self, "_methane_kgph", 0.0)

        self.design_results["AD inlet TS (kg/hr)"] = TS_in
        self.design_results["AD inlet VS (kg/hr)"] = VS_in
        self.design_results["Biodegradable pool (kg/hr)"] = biodegradable_pool
        self.design_results["Biodegradable destroyed (kg/hr)"] = biodegradable_destroyed
        self.design_results["Methane production (kg/hr)"] = methane_kgph
        self.design_results["CO2 production (kg/hr)"] = getattr(self, "_co2_kgph", 0.0)
        self.design_results["H2S production (kg/hr)"] = getattr(self, "_h2s_kgph", 0.0)
        self.design_results["Raw biogas total (kg/hr)"] = getattr(self, "_biogas_mass_kgph", 0.0)
        self.design_results["Raw biogas total (Nm3/hr)"] = getattr(self, "_biogas_Nm3ph", 0.0)
        self.design_results["Digestate outflow (kg/hr)"] = getattr(self, "_digestate_mass_kgph", 0.0)

        if V_liquid > 0:
            self.design_results["Organic loading proxy (kg/m3-h)"] = (
                biodegradable_pool / V_liquid
            )
            self.design_results["Methane productivity (kg/m3-h)"] = (
                methane_kgph / V_liquid
            )

        if VS_in > 0:
            self.design_results["CH4 yield (kg/kg VS fed)"] = methane_kgph / VS_in

        self.power_utility.consumption = mixing_kW
        self.power_utility.production = 0.0

        if Q_kJph > 0:
            try:
                self.add_heat_utility(Q_kJph, T_in, T_target)
            except TypeError:
                self.add_heat_utility(Q_kJph, T_in)

    def _cost(self):
        V_each = self.design_results["Digester volume each (m3)"]
        N = int(self.design_results["Number of digesters"])

        self.F_BM["Anaerobic digester"] = 1.0

        C_each = interp_capex(V_each)
        C_total = N * C_each

        self.design_results["ADBC CAPEX each ($)"] = C_each
        self.design_results["ADBC CAPEX total ($)"] = C_total
        self.baseline_purchase_costs["Anaerobic digester"] = C_total

        if self.maintenance_usd_per_m3_yr is not None:
            total_volume = self.design_results["Total digester volume (m3)"]
            annual_maintenance = self.maintenance_usd_per_m3_yr * total_volume
            self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

            if annual_maintenance > 0:
                operating_hours_per_yr = 330.0 * 24.0
                self.add_OPEX = {
                    "Digester maintenance": annual_maintenance / operating_hours_per_yr
                }

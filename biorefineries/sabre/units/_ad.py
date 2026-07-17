# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Anaerobic digestion (AD) and the surrounding biogas train: methanogenic and
acidogenic digestion, biogas upgrading, H2S removal, and digestate
dewatering (screw press / decanter centrifuge).

AnaerobicDigester
-----------------
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

import math
from typing import Dict, Iterable, Optional

import biosteam as bst
from biosteam.units.decorators import cost

__all__ = (
    'AnaerobicDigester', 'GAL_TO_M3', 'NM3_PER_KMOL', 'interp_capex',
    'AcidogenicDigester', 'BiogasUpgrading', 'H2SRemoval',
    'DigestateDecanterCentrifuge', 'DigestateScrewPress',
)

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
        data/ad.yaml `ad_performance.methanogenic.vs_destruction`, labeled "Global CSTR
        vs_destruction -- uniform across all pretreatment cases." No named
        literature source; yaml's own note: "Conservative lower bound of
        20-70% range for full-scale mesophilic CSTRs. No continuous CSTR data
        for pretreated pelagic Sargassum exists in literature." Scoping
        assumption, not a citation.
    ch4_kg_per_kg_vs_fed (0.100):
        This is the `press_mill_only` baseline pretreatment case value
        (data/pretreatment.yaml `pretreatment_ad.press_mill_only.ad_effects.ch4_kg_per_kg_vs_fed`),
        not a universal constant -- the 4 other pretreatment cases in yaml use
        different values (0.123 enzymatic, 0.165 peroxide, 0.277 combined_PE,
        0.261 combined_PTE), each cited to a different row of the same table.
        Source: Chikani-Cabrera et al. 2022, Table 4 -- "224.19 +/- 9.45 L
        CH4/kg VS," scaled to 0.100 kg/kg VS for continuous CSTR (yaml key
        `chikani_cabrera_2022_physical_C`).
    raw_biogas_molfrac (Methane 0.55, CarbonDioxide 0.42, HydrogenSulfide 0.03):
        data/ad.yaml `ad_performance.methanogenic.raw_biogas_molfrac`, labeled "Global
        fallback biogas composition (overridden per pretreatment case below)."
        No named source -- scoping default, overridden per-case in real runs.
    headspace_frac (0.25):
        data/ad.yaml `ad.gas_storage_frac_of_total_volume`. No named
        source.
    maintenance_usd_per_m3_yr (10.0):
        data/ad.yaml `ad.cost.maintenance_usd_per_m3_yr`, sourced to
        the ADBC spreadsheet tool (`ADBCv2-48 (2).xlsm`, sheet
        "AnaerobicDigester", cell AB18).
    digestible_IDs, biodegradability:
        data/ad.yaml `ad_performance.digestible_IDs` (shared with the
        acidogenic mode) / `ad_performance.methanogenic.biodegradability`.
        No named source for the individual factors beyond the yaml values
        themselves.
    hrt_days, slurry_density_kg_per_m3, max_single_digester_volume_MG,
    mixing_W_per_m3, influent_temperature_K, target_temperature_K,
    cp_kJ_per_kgK:
        data/ad.yaml `ad` section, matching values already used at the
        real call site (`systems/_ad_biomethane_system.py`). No named literature
        source for these beyond the yaml values themselves.
    ADBC_VOL_M3 / ADBC_CAPEX interpolation table (used by `interp_capex`,
    the actual installed-cost calculation):
        UNVERIFIED PROVENANCE. data/ad.yaml has no entry for this
        six-point table. A *separate*, structurally unused single-point
        anchor (`ad.cost.base_volume_m3`/`base_capex_usd`, stored on this
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


class AcidogenicDigester(bst.Unit):
    """
    Acidogenic / arrested AD unit for VFA platform modeling. Structurally
    parallel to `AnaerobicDigester` (multi-train, ADBC-interpolated capital
    cost via the same `interp_capex`/`ADBC_VOL_M3`/`ADBC_CAPEX` table) but
    produces a VFA-rich broth instead of biogas, and has no biodegradability
    weighting per component -- all mass in `digestible_IDs` is treated as
    equally available for VS destruction.

    Outputs:
    0) offgas
    1) acidogenic_broth

    Sources
    -------
    vs_destruction (0.50):
        data/ad.yaml `ad_performance.acidogenic.cases.seaweed_arrested_fitted.
        vs_destruction`. No named literature source in yaml for this value.
    vfa_kg_per_kg_vs (0.55):
        data/ad.yaml `ad_performance.acidogenic.cases.seaweed_arrested_fitted.
        vfa_kg_per_kg_vs`. No named literature source in yaml for this
        value (yaml's `vfa_fermentation.sources.seaweed_arrested_anchor`
        note describes a related but distinct downstream observation --
        "Peak ~14.5 g/L total VFA... biogas suppression 96%" -- from a
        preprint on arrested AD; it is not cited as the source of this
        specific 0.55 yield number, so not repeated as one here).
    vfa_split (AceticAcid 0.648, PropionicAcid 0.186, ButyricAcid 0.084,
    ValericAcid 0.082, HexanoicAcid 0.000):
        data/ad.yaml `ad_performance.acidogenic.cases.seaweed_arrested_fitted.
        vfa_split`. No named literature source in yaml. This is the
        __init__ default *object* (None, deferring to
        `_resolve_vfa_split()`'s own separate internal fallback split
        below) -- not literally hardcoded as this __init__'s default value.
        `_resolve_vfa_split()`'s internal fallback
        (AceticAcid 0.40/PropionicAcid 0.10/ButyricAcid 0.30/ValericAcid
        0.05/HexanoicAcid 0.15) is left as-is: it is a generic placeholder
        for *any* thermo with the five common VFA IDs, not a stale copy of
        this case's fitted split, so it is not "corrected" to the
        Sargassum-specific numbers above -- doing so would conflate a
        general fallback with a case-specific fitted result.
    digestible_IDs:
        data/ad.yaml `ad_performance.digestible_IDs` (shared with the
        methanogenic mode). The
        pre-conversion default additionally listed `Xylan`, `Mannan`,
        `Galactan` (real thermo chemicals, just absent from this yaml list)
        and `Cellulose` (not a chemical defined in `_chemicals.py` at all --
        a stale leftover reference). All four dropped to match yaml exactly.
    hrt_days (15.0):
        data/ad.yaml `ad.vfa.hrt_days`, sourced to `aad_review_window`:
        "15 d retained as a midpoint design basis within the 10-20 d AAD
        operating window; not a direct fitted value for this exact
        Sargassum system."
    influent_temperature_K, target_temperature_K:
        assumptions.yaml `vfa_ad` section, sourced to `mesophilic_operation`:
        "Mesophilic operation near 35 C."
    headspace_frac (0.25), slurry_density_kg_per_m3,
    max_single_digester_volume_MG, mixing_W_per_m3, cp_kJ_per_kgK:
        assumptions.yaml `vfa_ad` section, matching values already used at
        the real call site (`systems/_ad_vfa_system.py`). No named source
        beyond the yaml values themselves.
    ADBC_VOL_M3 / ADBC_CAPEX interpolation table (defined above in this
    module, used by `interp_capex`, the actual installed-cost
    calculation):
        Same unverified provenance caveat as `AnaerobicDigester` -- see
        that unit's docstring and CONVERSION_NOTES.md. Not re-verified here.
    enable_heat_shock and the hs_* parameters:
        No assumptions.yaml section covers heat-shock operation at all --
        these are pure code-level scoping parameters, off by default, not
        grounded in any cited source.
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
        vfa_kg_per_kg_vs: float = 0.55,
        vfa_split: Optional[Dict[str, float]] = None,
        digestible_IDs: Optional[Iterable[str]] = None,
        produce_offgas_co2: bool = True,
        hrt_days: float = 15.0,
        slurry_density_kg_per_m3: float = 1000.0,
        headspace_frac: float = 0.25,
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
            "Arabinan",
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
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


class BiogasUpgrading(bst.Unit):
    """
    Membrane biogas upgrading to pipeline-quality biomethane.

    Inputs:
        ins[0]: raw biogas (post H2S removal)

    Outputs:
        outs[0]: biomethane
        outs[1]: offgas

    Kept as custom _design/_cost rather than an @cost decorator: annual
    maintenance OPEX here is a fraction of this unit's own installed cost,
    which @cost only computes inside its decorator-owned _cost() -- after
    _design() has already run. There's no clean way to read that
    not-yet-computed value early without duplicating the cost formula (a
    latent-bug risk if the two copies ever drift), so this unit keeps an
    explicit _cost() instead, unlike H2SRemoval (whose OPEX depends only on
    _design()-computed flow, not its own CAPEX).

    Sources
    -------
    co2_removal (0.95), electricity_kwh_per_Nm3_raw (0.25),
    maintenance_frac_of_capex_per_yr (0.035):
        assumptions.yaml `biogas_upgrading` section, sourced to IEA
        Bioenergy Task 37, "Biomethane status and factors affecting market
        development" (2014). Electricity: "typically 0.2-0.3 kWh/Nm3 raw
        biogas for mature upgrading technologies (incl. membranes)."
        Maintenance: "annual maintenance cost for membranes at ~3-4% of
        investment cost."
    ch4_recovery (0.99):
        assumptions.yaml `biogas_upgrading.ch4_recovery`. yaml's
        `sources.methane_slip` note (also IEA Bioenergy Task 37, 2014)
        states "manufacturers guarantee methane losses below ~0.5-2% for
        new plants" -- consistent with a 0.99 (1%) recovery/slip choice, but
        yaml does not cite a number specifically for ch4_recovery itself.
    capex_usd_per_Nm3ph_raw (2200):
        assumptions.yaml `biogas_upgrading.sources.capex`: IEA-ETSAP
        E-TechDS: Biogas and Bio-syngas Production (Dec 2013). "Investment
        cost for biogas-to-biomethane upgrading units reported ~USD
        1950-2600 per Nm3/h for larger raw gas capacities" -- 2200 falls
        within that range.
    """

    _N_ins = 1
    _N_outs = 2  # biomethane, offgas

    COST_ITEM = "Membrane upgrading skid"
    F_BM = {COST_ITEM: 1.0}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        ch4_recovery=0.99,
        co2_removal=0.95,
        electricity_kwh_per_Nm3_raw=0.25,
        capex_usd_per_Nm3ph_raw=2200.0,
        maintenance_frac_of_capex_per_yr=0.035,
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.ch4_recovery = float(ch4_recovery)
        self.co2_removal = float(co2_removal)
        self.electricity_kwh_per_Nm3_raw = float(electricity_kwh_per_Nm3_raw)
        self.capex_usd_per_Nm3ph_raw = float(capex_usd_per_Nm3ph_raw)
        self.maintenance_frac_of_capex_per_yr = float(maintenance_frac_of_capex_per_yr)

        # make sure instance has the same F_BM mapping
        self.F_BM = dict(type(self).F_BM)

        if not (0.0 <= self.ch4_recovery <= 1.0):
            raise ValueError("ch4_recovery must be between 0 and 1.")
        if not (0.0 <= self.co2_removal <= 1.0):
            raise ValueError("co2_removal must be between 0 and 1.")
        if self.electricity_kwh_per_Nm3_raw < 0:
            raise ValueError("electricity_kwh_per_Nm3_raw must be >= 0.")
        if self.capex_usd_per_Nm3ph_raw < 0:
            raise ValueError("capex_usd_per_Nm3ph_raw must be >= 0.")
        if not (0.0 <= self.maintenance_frac_of_capex_per_yr <= 1.0):
            raise ValueError("maintenance_frac_of_capex_per_yr must be between 0 and 1.")

    def _run(self):
        raw = self.ins[0]
        biomethane, offgas = self.outs

        biomethane.empty()
        offgas.empty()

        biomethane.phase = "g"
        offgas.phase = "g"

        ch4_in = float(raw.imol["Methane"])
        co2_in = float(raw.imol["CarbonDioxide"])

        ch4_to_bm = self.ch4_recovery * ch4_in
        biomethane.imol["Methane"] = ch4_to_bm
        offgas.imol["Methane"] = ch4_in - ch4_to_bm

        co2_to_off = self.co2_removal * co2_in
        offgas.imol["CarbonDioxide"] = co2_to_off
        biomethane.imol["CarbonDioxide"] = co2_in - co2_to_off

        for cid in raw.chemicals.IDs:
            if cid in ("Methane", "CarbonDioxide"):
                continue
            n = float(raw.imol[cid])
            if n > 0:
                offgas.imol[cid] = n

    def _design(self):
        raw = self.ins[0]
        biomethane, offgas = self.outs

        # Use dry raw biogas basis for sizing and power
        dry_gas_IDs = ["Methane", "CarbonDioxide", "HydrogenSulfide", "Nitrogen", "Oxygen"]
        n_kmolph_dry = sum(float(raw.imol[i]) for i in dry_gas_IDs if i in raw.chemicals.IDs)

        Q_Nm3ph_dry = 22.414 * n_kmolph_dry
        self.design_results["Raw biogas flow (Nm3/h, dry)"] = Q_Nm3ph_dry

        kW = self.electricity_kwh_per_Nm3_raw * Q_Nm3ph_dry
        self.design_results["Electricity demand (kW)"] = kW

        # Explicit overwrite
        self.power_utility.consumption = kW
        self.power_utility.production = 0.0

        if biomethane.F_mol > 0:
            self.design_results["Biomethane CH4 mol%"] = (
                100.0 * float(biomethane.imol["Methane"]) / float(biomethane.F_mol)
            )
        else:
            self.design_results["Biomethane CH4 mol%"] = 0.0

        self.design_results["Methane slip (kmol/h)"] = float(offgas.imol["Methane"])

        ch4_in = float(raw.imol["Methane"])
        if ch4_in > 0:
            self.design_results["Methane recovery actual"] = (
                float(biomethane.imol["Methane"]) / ch4_in
            )
        else:
            self.design_results["Methane recovery actual"] = 0.0

    def _cost(self):
        Q_Nm3ph = self.design_results.get("Raw biogas flow (Nm3/h, dry)")
        if Q_Nm3ph is None:
            raw = self.ins[0]
            dry_gas_IDs = ["Methane", "CarbonDioxide", "HydrogenSulfide", "Nitrogen", "Oxygen"]
            n_kmolph_dry = sum(float(raw.imol[i]) for i in dry_gas_IDs if i in raw.chemicals.IDs)
            Q_Nm3ph = 22.414 * n_kmolph_dry
            self.design_results["Raw biogas flow (Nm3/h, dry)"] = Q_Nm3ph

        capex = self.capex_usd_per_Nm3ph_raw * float(Q_Nm3ph)
        self.baseline_purchase_costs[self.COST_ITEM] = capex

        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Membrane upgrading maintenance": annual_maintenance / operating_hours_per_yr
            }


@cost('Raw biogas flow (Nm3/h)', 'H2S removal (iron sponge)', units='Nm3/h',
      CE=567.5, cost=450_000.0, S=1700.0, n=0.7, BM=1.0)
class H2SRemoval(bst.Unit):
    """
    Iron sponge H2S removal unit for raw biogas desulfurization.

    Removes hydrogen sulfide from raw biogas before membrane upgrading.
    H2S must be removed because it degrades polymer membranes rapidly and
    is toxic above ~200 ppm in pipeline-quality biomethane.

    Technology: Iron sponge (dry oxidation) -- the standard, lowest-cost
    option for biogas desulfurization at this scale.

    Modeling approach:
    - Pass-through for CH4 and CO2 (not affected by iron sponge)
    - H2S captured to near-zero in treated gas
    - Capital cost scaled to raw biogas flow (Nm3/h) via a power-law @cost item
    - Reagent cost (iron sponge replacement) as add_OPEX, computed in _design
      since it depends on raw biogas flow, not on the unit's own installed cost

    Inputs:
        ins[0]: raw biogas (CH4 + CO2 + H2S)

    Outputs:
        outs[0]: treated biogas (H2S removed to near-zero)
        outs[1]: spent media / captured H2S (solid waste, negligible mass)

    Sources
    -------
    ref_installed_cost_usd ($450,000 at 1,700 Nm3/hr, n=0.7):
        assumptions.yaml `h2s_removal.sources.capex`: Diaz, I.; Ramos, I.;
        Fdz-Polanco, M. Bioresour. Technol. 2015, 192, 280-286. Original
        anchor from Abatzoglou & Boivin (2009); scale exponent 0.7 per
        Green & Perry (2008) convention.
    h2s_removal_efficiency (0.99):
        assumptions.yaml `h2s_removal.h2s_removal_efficiency` -- this is the
        modeled value, chosen conservatively. The yaml citation (Choudhury
        et al. Energies 2019, 12, 4605) reports the cited iron-sponge media
        as achieving >99.9% removal in practice -- 0.99 is not a literal
        transcription of that number, it's a deliberately conservative
        modeling choice.
    reagent_cost_usd_per_Nm3_raw (0.002):
        assumptions.yaml `h2s_removal.sources.reagent`: IEA Bioenergy Task 37
        (2014), yaml note gives a $0.001-0.003/Nm3 range; 0.002 is the
        midpoint.
    """

    _N_ins = 1
    _N_outs = 2   # treated_biogas, spent_media
    _units = {'Raw biogas flow (Nm3/h)': 'Nm3/h'}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        h2s_removal_efficiency: float = 0.99,
        reagent_cost_usd_per_Nm3_raw: float = 0.002,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.h2s_removal_efficiency = float(h2s_removal_efficiency)
        self.reagent_cost_usd_per_Nm3_raw = float(reagent_cost_usd_per_Nm3_raw)

    def _run(self):
        raw = self.ins[0]
        treated, spent = self.outs

        treated.empty()
        spent.empty()
        treated.phase = "g"
        spent.phase = "l"  # negligible solid waste

        ids = set(raw.chemicals.IDs)

        # Pass CH4 and CO2 entirely to treated gas
        for cid in ("Methane", "CarbonDioxide"):
            if cid in ids:
                treated.imol[cid] = float(raw.imol[cid])

        # Remove H2S by efficiency
        if "HydrogenSulfide" in ids:
            h2s_in = float(raw.imol["HydrogenSulfide"])
            h2s_removed = self.h2s_removal_efficiency * h2s_in
            h2s_out = h2s_in - h2s_removed
            treated.imol["HydrogenSulfide"] = h2s_out
            # Spent media is mostly iron sulfide solids —> negligible mass in model
            # Not explicitly tracked as a chemical stream

        # Pass any other gas components to treated
        for cid in ids:
            if cid in ("Methane", "CarbonDioxide", "HydrogenSulfide"):
                continue
            n = float(raw.imol[cid])
            if n > 0:
                treated.imol[cid] = n

    def _design(self):
        raw = self.ins[0]

        # Raw biogas flow in Nm3/h (ideal gas at STP)
        n_kmolph = float(raw.F_mol)
        Q_Nm3ph = 22.414 * n_kmolph

        self.design_results["Raw biogas flow (Nm3/h)"] = Q_Nm3ph
        self.design_results["H2S removal efficiency"] = self.h2s_removal_efficiency

        if "HydrogenSulfide" in raw.chemicals.IDs:
            h2s_in_ppm = (
                float(raw.imol["HydrogenSulfide"]) / float(raw.F_mol) * 1e6
                if raw.F_mol > 0 else 0.0
            )
            treated = self.outs[0]
            h2s_out_ppm = (
                float(treated.imol["HydrogenSulfide"]) / float(treated.F_mol) * 1e6
                if treated.F_mol > 0 else 0.0
            )
            self.design_results["H2S inlet (ppm mol)"] = h2s_in_ppm
            self.design_results["H2S outlet (ppm mol)"] = h2s_out_ppm

        # Iron sponge is passive —> no electricity
        self.power_utility.consumption = 0.0

        # Reagent (iron sponge media replacement) as add_OPEX.
        # Computed here (not in _cost) because it scales with raw biogas
        # flow, not with this unit's own installed cost -- @cost owns _cost.
        reagent_usd_per_hr = self.reagent_cost_usd_per_Nm3_raw * Q_Nm3ph
        if reagent_usd_per_hr > 0:
            self.add_OPEX = {
                "Iron sponge media replacement": reagent_usd_per_hr
            }
            self.design_results["Reagent cost ($/hr)"] = reagent_usd_per_hr
            self.design_results["Reagent cost ($/yr)"] = reagent_usd_per_hr * 330.0 * 24.0


class DigestateDecanterCentrifuge(bst.Unit):
    """
    Post-AD digestate decanter centrifuge (solid-liquid separation).

    An alternative to DigestateScrewPress for digestate dewatering: a decanter
    centrifuge achieves higher TS capture than a screw press (SYSTEMIC D3.2
    report: centrifuge ~59+-17% vs screw press ~33+-14%) at a different cost
    and electricity profile. Not wired into any of the default system builders
    today — kept available as a selectable alternative, same as the
    pretreatment options in systems._ad_biomethane_system.

    Splits digestate into:
        - cake (soil_amendment): captured solids + entrained water to hit cake_moisture_frac
        - centrate (liquid_digestate): remaining liquid + uncaptured solids

    Assumptions:
        - "Solids" are explicitly defined chemical IDs (default: Cellulose, Ash)
        - Everything else is treated as liquid and starts in centrate

    Sizing:
        - Throughput_tph = inlet_mass_kgph / 1000
        - N_parallel = ceil(Throughput_tph / capacity_tph_each)

    Costing (USD):
        - Total purchased cost = N_parallel * centrifuge_purchase_cost_usd_each

    Not currently constructed by any system builder in this codebase
    (confirmed by grep -- only this file and units/__init__.py reference
    the class name) -- so every parameter below is only ever exercised via
    its own __init__ default, never overridden by a real call site. Unlike
    `DigestateScrewPress`, `solids_IDs` here is an active whitelist (see
    _run): only chemicals in this list are captured to cake at
    `ts_capture_frac`; everything else goes entirely to centrate.

    Sources
    -------
    solids_IDs (Cellulose, Ash), ts_capture_frac (0.78),
    cake_moisture_frac (0.75), capacity_tph_each (50.0),
    centrifuge_purchase_cost_usd_each ($297,500):
        assumptions.yaml `digestate_decanter_centrifuge` section (no
        `sources:` sub-block in yaml for this section, so no named
        literature citation for any of these values beyond the yaml values
        themselves). `solids_IDs` matches yaml exactly (`["Cellulose",
        "Ash"]`) -- unlike `DigestateScrewPress`, this default was not
        stale, but note "Cellulose" is not a chemical defined anywhere in
        `_chemicals.py`, so in practice only Ash is ever captured to cake
        by this whitelist; not changed since it matches assumptions.yaml.
        `centrifuge_purchase_cost_usd_each` was `None` (i.e. `_cost()`
        returns early, giving $0 installed cost) -- corrected to yaml's
        297500.0 so a standalone instantiation of this currently-unwired
        unit doesn't silently cost nothing.
    ts_capture_frac comparison to DigestateScrewPress (module docstring,
    pre-existing text, not yaml):
        SYSTEMIC D3.2 report (2021): "centrifuge ~59+-17% vs screw press
        ~33+-14%" -- yaml's modeled 0.78 for this unit is higher than that
        reported ~59% average, not verified further.
    Electricity intensity (1.0 kWh/m3, hardcoded in _design, not an
    __init__ param):
        Tchobanoglous et al., Wastewater Engineering, 5th ed. (pre-existing
        module docstring citation, not independently verified in this
        conversion).
    """

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID="", ins=None, outs=(),
                 solids_IDs=("Cellulose", "Ash"),
                 ts_capture_frac=0.78,
                 cake_moisture_frac=0.75,
                 capacity_tph_each=50.0,
                 centrifuge_purchase_cost_usd_each=297_500.0,
                 F_BM=1.0,
                 **kwargs):
        super().__init__(ID, ins, outs, **kwargs)

        self.solids_IDs = tuple(solids_IDs)
        self.ts_capture_frac = float(ts_capture_frac)
        self.cake_moisture_frac = float(cake_moisture_frac)

        self.capacity_tph_each = float(capacity_tph_each)
        self.centrifuge_purchase_cost_usd_each = centrifuge_purchase_cost_usd_each

        self.F_BM_default = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        cake, centrate = self.outs

        cake.empty()
        centrate.empty()
        cake.phase = "l"
        centrate.phase = "l"

        cap = min(max(self.ts_capture_frac, 0.0), 1.0)

        # split defined solids between cake and centrate
        for sid in self.solids_IDs:
            m = float(feed.imass[sid])
            m_cake = cap * m
            cake.imass[sid] = m_cake
            centrate.imass[sid] = m - m_cake

        # send everything else to centrate
        for chem_id in feed.chemicals.IDs:
            if chem_id in self.solids_IDs:
                continue
            centrate.imass[chem_id] = float(feed.imass[chem_id])

        # entrained water to cake to meet target cake moisture
        TS_cake = sum(float(cake.imass[sid]) for sid in self.solids_IDs)
        if TS_cake > 0 and "Water" in feed.chemicals.IDs:
            mfrac = self.cake_moisture_frac
            if not (0.0 < mfrac < 1.0):
                raise ValueError("cake_moisture_frac must be between 0 and 1")

            # moisture = water/(water+TS) => water = mfrac/(1-mfrac)*TS
            water_needed = (mfrac / (1.0 - mfrac)) * TS_cake
            water_available = float(centrate.imass["Water"])
            water_to_cake = min(water_needed, water_available)

            cake.imass["Water"] += water_to_cake
            centrate.imass["Water"] -= water_to_cake

    def _design(self):
        F_kgph = float(self.ins[0].F_mass)
        F_tph = F_kgph / 1000.0  # metric ton/h

        cap = float(self.capacity_tph_each)
        if cap <= 0.0:
            raise ValueError("capacity_tph_each must be > 0")

        N = int(math.ceil(F_tph / cap))

        self.design_results["Throughput_kgph"] = F_kgph
        self.design_results["Throughput_tph"] = F_tph
        self.design_results["Capacity_tph_each"] = cap
        self.design_results["N_centrifuges_parallel"] = N

        # Electricity: 1.0 kWh/m3 digestate for decanter centrifuge
        # Source: Tchobanoglous et al., Wastewater Engineering 5th ed.
        F_m3ph = float(self.ins[0].F_vol)
        kWh_per_m3 = 1.0
        power_kW = kWh_per_m3 * F_m3ph
        self.power_utility.consumption = power_kW
        self.design_results["Feed flow (m3/h)"] = F_m3ph
        self.design_results["Power (kW)"] = power_kW

    def _cost(self):
        if self.centrifuge_purchase_cost_usd_each is None:
            return

        N = int(self.design_results.get("N_centrifuges_parallel", 1))
        u_cost = float(self.centrifuge_purchase_cost_usd_each)

        key = "Decanter centrifuge (digestate)"
        self.baseline_purchase_costs[key] = N * u_cost
        self.F_BM[key] = self.F_BM_default


class DigestateScrewPress(bst.Unit):
    """
    Post-AD digestate screw press (solid-liquid separation).

    Splits digestate into:
        - cake (soil_amendment): captured solids + entrained water to hit cake_moisture_frac
        - pressate (liquid_digestate): remaining liquid + uncaptured solids

    Assumptions:
        - "Solids" are explicitly defined chemical IDs (default: Cellulose, Ash)
        - Everything else is treated as liquid and starts in pressate
        - Default DM (TS) capture to solids is lower than centrifuge:
            SE_DM screw press ~33+-14% vs centrifuge ~59+-17% (SYSTEMIC D3.2 report, 2021)
        - Solid fraction DM after screw press ~23+-8% => moisture ~77% (SYSTEMIC database)

    Sizing:
        - Throughput_tph = inlet_mass_kgph / 1000
        - N_parallel = ceil(Throughput_tph / capacity_tph_each)

    Energy:
        - Electricity intensity: 0.67 kWh/m3 digestate treated (SYSTEMIC database; n=13)

    Costing:
        - Purchased cost per press from SYSTEMIC Table 2-11 (CAPEX vs ton/h), interpolated
        - Total purchased cost = N_parallel * cost_per_press
        - Optional polymer dosing unit per press (12k-50k EUR) can be included

    Constructed twice in this codebase (systems/_ad_biomethane_system.py's "SP"
    and systems/_ad_vfa_system.py's "SP_VFA"), both routing their own
    digestate streams through the same class with different
    assumptions.yaml-sourced parameters.

    Sources
    -------
    ts_capture_frac (0.40), cake_moisture_frac (0.45), kWh_per_m3 (0.67):
        assumptions.yaml `digestate_screw_press` section. No `sources:`
        sub-block in yaml for this section -- citations below are from this
        module's own pre-existing docstring text, not yaml: SYSTEMIC D3.2
        report (2021), "SE_DM screw press ~33+-14% vs centrifuge ~59+-17%"
        (yaml's modeled ts_capture_frac of 0.40 is close to but not
        identical to that reported ~33% average); SYSTEMIC database
        (n=13) for both cake solids fraction (~23% DM -> yaml's 0.45
        moisture is a rounder number than the ~77% moisture that ~23% DM
        implies -- not verified further) and the 0.67 kWh/m3 electricity
        intensity.
    capex_eur_table (SYSTEMIC Table 2-11, 10 (capacity_tph, capex_eur)
    points):
        assumptions.yaml's `digestate_screw_press` section has no
        `capex_eur_table` key, so this hardcoded table -- not the yaml
        section -- is what's actually used at every real call site. Cited
        to this module's pre-existing docstring text: "SYSTEMIC Table
        2-11 (CAPEX vs ton/h)."
    eur_to_usd (1.19):
        assumptions.yaml `digestate_screw_press.eur_to_usd`. No named
        source in yaml.
    solids_IDs -- accepted and stored but NOT used by _run/_design/_cost:
        TS for capture purposes is instead computed dynamically in _run()
        as "everything except Water and dissolved_IDs," not as this
        explicit whitelist. This parameter is dead code regardless of its
        value -- left at its pre-existing default (which itself references
        "Cellulose," not a chemical defined anywhere in `_chemicals.py`)
        rather than "corrected" to assumptions.yaml's `digestate_screw_press
        .solids_IDs` list, since doing so would misleadingly imply the
        parameter has an effect it doesn't have.
    """

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID="", ins=None, outs=(),
                 solids_IDs=("Cellulose", "Ash"),
                 dissolved_IDs=None,            # chemicals treated as dissolved — always route to pressate
                 ts_capture_frac=0.40,          # assumptions.yaml digestate_screw_press.ts_capture_frac
                 cake_moisture_frac=0.45,       # assumptions.yaml digestate_screw_press.cake_moisture_frac
                 capacity_tph_each=6.0,         # aligns with reported 6 ton/h energy datapoint
                 kWh_per_m3=0.67,               # SYSTEMIC avg electricity intensity
                 capex_eur_table=None,          # override if you want your own table
                 eur_to_usd=1.19,               # assumptions.yaml digestate_screw_press.eur_to_usd
                 include_polymer_dosing=False,
                 polymer_dosing_cost_eur_each=0.0,  # set within 12k–50k if included
                 F_BM=1.0,
                 **kwargs):
        super().__init__(ID, ins, outs, **kwargs)

        self.solids_IDs = tuple(solids_IDs)

        # Dissolved chemicals pass entirely to pressate regardless of ts_capture_frac
        # VFAs are soluble acids — they must not be captured in the cake
        # Defaults cover the standard VFA set; override via dissolved_IDs if needed
        if dissolved_IDs is None:
            dissolved_IDs = (
                "AceticAcid", "PropionicAcid", "ButyricAcid",
                "ValericAcid", "HexanoicAcid",
                "CarbonDioxide", "Ammonia",
            )
        self.dissolved_IDs = tuple(dissolved_IDs)

        self.ts_capture_frac = float(ts_capture_frac)
        self.cake_moisture_frac = float(cake_moisture_frac)

        self.capacity_tph_each = float(capacity_tph_each)

        self.kWh_per_m3 = float(kWh_per_m3)

        # SYSTEMIC Table 2-11 (capacity_tph, capex_eur)
        # (2, 25625), (3, 36187), (4, 30416), (5, 39571),
        # (6.25, 28750), (8, 34583), (9, 44400), (10, 24062),
        # (12, 52500), (15.5, 17000)
        self.capex_eur_table = capex_eur_table or [
            (2.0, 25625.0),
            (3.0, 36187.0),
            (4.0, 30416.0),
            (5.0, 39571.0),
            (6.25, 28750.0),
            (8.0, 34583.0),
            (9.0, 44400.0),
            (10.0, 24062.0),
            (12.0, 52500.0),
            (15.5, 17000.0),
        ]

        self.eur_to_usd = float(eur_to_usd)

        self.include_polymer_dosing = bool(include_polymer_dosing)
        self.polymer_dosing_cost_eur_each = float(polymer_dosing_cost_eur_each)

        self.F_BM_default = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        cake, pressate = self.outs

        cake.empty()
        pressate.empty()
        cake.phase = "l"
        pressate.phase = "l"

        cap = min(max(self.ts_capture_frac, 0.0), 1.0)

        water_id = "Water"
        ids = feed.chemicals.IDs

        # Dissolved chemicals (VFAs, CO2, Ammonia etc.) pass entirely to pressate
        # They are soluble and cannot be mechanically captured in a screw press cake
        dissolved_set = set(self.dissolved_IDs)

        # Define TS as everything except Water AND dissolved chemicals
        ts_ids = [i for i in ids if i != water_id and i not in dissolved_set]

        # Route dissolved chemicals entirely to pressate
        for i in dissolved_set:
            if i in ids:
                pressate.imass[i] = float(feed.imass[i])

        # --- split TS by capture fraction (based on true solids only) ---
        TS_total = sum(float(feed.imass[i]) for i in ts_ids)
        TS_to_cake = cap * TS_total

        if TS_total > 0:
            for i in ts_ids:
                m = float(feed.imass[i])
                m_cake = (m / TS_total) * TS_to_cake
                cake.imass[i] = m_cake
                pressate.imass[i] = m - m_cake

        # all water starts in pressate
        if water_id in ids:
            pressate.imass[water_id] = float(feed.imass[water_id])

        # --- add entrained water to cake to meet target cake moisture ---
        TS_cake = sum(float(cake.imass[i]) for i in ts_ids)
        if TS_cake > 0 and water_id in ids:
            mfrac = self.cake_moisture_frac
            if not (0.0 < mfrac < 1.0):
                raise ValueError("cake_moisture_frac must be between 0 and 1")

            # moisture = water/(water+TS) => water = mfrac/(1-mfrac)*TS
            water_needed = (mfrac / (1.0 - mfrac)) * TS_cake
            water_available = float(pressate.imass[water_id])
            water_to_cake = min(water_needed, water_available)

            cake.imass[water_id] += water_to_cake
            pressate.imass[water_id] -= water_to_cake

    def _design(self):
        F_kgph = float(self.ins[0].F_mass)
        F_tph = F_kgph / 1000.0  # metric ton/h

        cap = float(self.capacity_tph_each)
        if cap <= 0.0:
            raise ValueError("capacity_tph_each must be > 0")

        N = int(math.ceil(F_tph / cap))

        self.design_results["Throughput_kgph"] = F_kgph
        self.design_results["Throughput_tph"] = F_tph
        self.design_results["Capacity_tph_each"] = cap
        self.design_results["N_screw_presses_parallel"] = N

        # Power: kW = (kWh/m3) * (m3/h)
        F_m3ph = float(self.ins[0].F_vol)
        kW = self.kWh_per_m3 * F_m3ph
        self.power_utility.consumption = kW

        self.design_results["Throughput_m3ph"] = F_m3ph
        self.design_results["kWh_per_m3"] = self.kWh_per_m3
        self.design_results["Power_kW"] = kW

    def _interp_capex_eur(self, cap_tph):
        pts = sorted(self.capex_eur_table, key=lambda x: x[0])
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]

        if cap_tph <= xs[0]:
            return ys[0]
        if cap_tph >= xs[-1]:
            return ys[-1]

        for i in range(len(xs) - 1):
            x0, x1 = xs[i], xs[i + 1]
            if x0 <= cap_tph <= x1:
                y0, y1 = ys[i], ys[i + 1]
                # linear interpolation
                t = (cap_tph - x0) / (x1 - x0)
                return y0 + t * (y1 - y0)

        return ys[-1]

    def _cost(self):
        N = int(self.design_results.get("N_screw_presses_parallel", 1))
        cap_each = float(self.capacity_tph_each)

        capex_eur_each = float(self._interp_capex_eur(cap_each))
        capex_usd_each = capex_eur_each * self.eur_to_usd

        total_usd = N * capex_usd_each

        if self.include_polymer_dosing:
            dosing_usd_each = self.polymer_dosing_cost_eur_each * self.eur_to_usd
            total_usd += N * dosing_usd_each

        key = "Screw press (digestate)"
        self.baseline_purchase_costs[key] = total_usd
        self.F_BM[key] = self.F_BM_default

        self.design_results["Capex_eur_each"] = capex_eur_each
        self.design_results["Capex_usd_each"] = capex_usd_each
        if self.include_polymer_dosing:
            self.design_results["Polymer_dosing_cost_eur_each"] = self.polymer_dosing_cost_eur_each
            self.design_results["Polymer_dosing_cost_usd_each"] = self.polymer_dosing_cost_eur_each * self.eur_to_usd
        self.design_results["N_parallel"] = N
        self.design_results["Total_purchase_cost_usd"] = total_usd

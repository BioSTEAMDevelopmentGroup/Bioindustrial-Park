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
"""

import math
from typing import Dict, Iterable, Optional

import biosteam as bst
from biosteam.units.decorators import cost

from biorefineries.sabre.utils import load_assumptions, OPERATING_HOURS_PER_YEAR

__all__ = (
    'AnaerobicDigester', 'AcidogenicDigester',
    'H2SRemoval', 'BiogasUpgrading',
    'DigestateDecanterCentrifuge', 'DigestateScrewPress',
)


# Conversion factors and ADBC data
GAL_TO_M3 = 0.003785411784  # US gallons -> m3
NM3_PER_KMOL = 22.414

# Loaded once at import time -- data/ad.yaml is the single source of truth
# for AnaerobicDigester.__init__'s default values below, so a yaml edit
# doesn't need a matching literal edit here too.
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]
_AD_METHANOGENIC = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
_ADP_METHANOGENIC = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
_AD_COST = _AD_SHARED["cost"]

# ADBC-based digester installed cost anchor points (m3, USD)
ADBC_VOL_M3 = [878, 1755, 2633, 3510, 5265, 8775]
ADBC_CAPEX = [1720964, 1750201, 1779439, 1808676, 1867151, 1984101]

def interp_ad_capex(volume_m3: float) -> float:
    """
    Linear interpolation/extrapolation of digester CAPEX based on the
    Anaerobic Digester Budget Calculator (ADBC).
    See https://csanr.wsu.edu/resources/enterprise-budget-calculator/
    """
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
    Convert a biodegradable fraction of Sargassum organics into raw biogas
    using component-specific biodegradability factors. The digestate contains
    undegraded material.

    Modeling notes:
    - Lumped methanogenic AD model
    - Methane production is based on: feed VS * ch4_kg_per_kg_vs_fed
    - Different components can have different biodegradability factors
    - Biodegradability + vs_destruction determine digestate solids reduction
    - Raw biogas composition is imposed from target mole fractions
    - Paralle-train, reactor sized using HRT and headspace
    - Enforce a maximum single-digester size and model parallel digesters
    - Estimate CAPEX using ADBC-based interpolation

    Parameters
    ----------
    ins : stream
        Pretreated/diluted feed slurry.
    outs : tuple[stream, stream]
        Raw biogas and digestate.
    vs_destruction : float
        Fraction of the biodegradability-weighted VS pool destroyed
        during digestion.
    ch4_kg_per_kg_vs_fed : float
        Methane yield per kg of VS fed (not VS destroyed).
    raw_biogas_molfrac : dict[str, float]
        Target raw biogas composition (mol fractions of Methane,
        CarbonDioxide, HydrogenSulfide); must sum to 1.
    digestible_IDs : Iterable[str]
        Chemical IDs eligible for VS destruction.
    biodegradability : dict[str, float]
        Per-component biodegradability factor (0-1) for IDs in
        `digestible_IDs`; components not listed default to 0
        (unavailable).
    hrt_days : float
        Hydraulic retention time, used to size total liquid digester
        volume.
    slurry_density_kg_per_m3 : float
        Feed slurry density, used to convert mass flow to volumetric
        flow for sizing.
    headspace_frac : float
        Fraction of total digester volume reserved for gas headspace.
    max_single_digester_volume_MG : float
        Maximum volume (million US gallons) of a single digester train;
        larger required volumes are split across parallel trains.
    maintenance_usd_per_m3_yr : float, optional
        Annual maintenance cost per m3 of total digester volume, added
        to `add_OPEX`. If None, no maintenance OPEX is added.
    mixing_W_per_m3 : float
        Mixing power intensity per m3 of liquid volume.
    influent_temperature_K : float
        Feed slurry inlet temperature, used to size heating duty.
    target_temperature_K : float
        Digester operating temperature, used to size heating duty.
    cp_kJ_per_kgK : float
        Feed slurry specific heat capacity, used to size heating duty.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/ad.yaml and data/pretreatment.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2  # biogas, digestate

    F_BM = {"Anaerobic digester": 1.0}

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        *,
        # Performance
        vs_destruction=_ADP_METHANOGENIC["vs_destruction"],
        ch4_kg_per_kg_vs_fed,
        raw_biogas_molfrac=_ADP_METHANOGENIC["raw_biogas_molfrac"],
        digestible_IDs=_ADP_METHANOGENIC["digestible_IDs"],
        biodegradability=_ADP_METHANOGENIC["biodegradability"],
        # Sizing
        hrt_days=_AD_METHANOGENIC["hrt_days"],
        slurry_density_kg_per_m3=_AD_METHANOGENIC["slurry_density_kg_per_m3"],
        headspace_frac=_AD_METHANOGENIC["gas_storage_frac_of_total_volume"],
        max_single_digester_volume_MG=_AD_METHANOGENIC["max_single_digester_volume_MG"],
        # Costing anchors
        maintenance_usd_per_m3_yr=_AD_COST.get("maintenance_usd_per_m3_yr"),
        # Utilities
        mixing_W_per_m3=_AD_METHANOGENIC["mixing_W_per_m3"],
        influent_temperature_K=_AD_METHANOGENIC["influent_temperature_K"],
        target_temperature_K=_AD_SHARED["temperature_regimes"]["mesophilic"]["temperature_K"],
        cp_kJ_per_kgK=_AD_METHANOGENIC["cp_kJ_per_kgK"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.vs_destruction = float(vs_destruction)
        self.ch4_kg_per_kg_vs_fed = float(ch4_kg_per_kg_vs_fed)
        self.raw_biogas_molfrac = dict(raw_biogas_molfrac)
        self.digestible_IDs = tuple(digestible_IDs)
        self.biodegradability = dict(biodegradability)
        self.mixing_W_per_m3 = float(mixing_W_per_m3)
        self.influent_temperature_K = float(influent_temperature_K)
        self.target_temperature_K = float(target_temperature_K)
        self.cp_kJ_per_kgK = float(cp_kJ_per_kgK)
        self.hrt_days = float(hrt_days)
        self.slurry_density_kg_per_m3 = float(slurry_density_kg_per_m3)
        self.headspace_frac = float(headspace_frac)
        self.max_single_digester_volume_m3 = float(max_single_digester_volume_MG) * 1e6 * GAL_TO_M3
        self.maintenance_usd_per_m3_yr = maintenance_usd_per_m3_yr
        self.F_BM = dict(self.F_BM)

    def _available_digestible_pool(self, stream):
        thermo_ids = set(self.chemicals.IDs)
        avail = {}

        for cid in self.digestible_IDs:
            if cid not in thermo_ids:
                continue

            mass = float(stream.imass[cid])

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

        chems = self.chemicals
        ids = set(chems.IDs)

        water_in = float(feed.imass["Water"])
        ash_in = float(feed.imass["Ash"])
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
        y_ch4 = float(y['Methane'])
        y_co2 = float(y['CarbonDioxide'])
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

        if abs(water_adjust) > 1e-9:
            current_water = float(digestate.imass["Water"])
            new_water = current_water + water_adjust
            if new_water < 0:
                raise RuntimeError(
                    f"AD water balance became negative for unit {self.ID}. "
                    f"Current water={current_water:.3f}, adjustment={water_adjust:.3f}."
                )
            digestate.imass["Water"] = new_water

    def _design(self):
        feed = self.ins[0]
        biogas, digestate = self.outs

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

        ids = set(self.chemicals.IDs)

        water_in = float(feed.imass["Water"])
        ash_in = float(feed.imass["Ash"])
        TS_in = max(float(feed.F_mass) - water_in, 0.0)
        VS_in = max(TS_in - ash_in, 0.0)

        avail = self._available_digestible_pool(feed)
        biodegradable_pool = sum(v["biodegradable_mass"] for v in avail.values())
        biodegradable_destroyed = min(self.vs_destruction * biodegradable_pool, biodegradable_pool)

        methane_kgph = float(biogas.imass["Methane"]) if "Methane" in ids else 0.0
        co2_kgph = float(biogas.imass["CarbonDioxide"]) if "CarbonDioxide" in ids else 0.0
        h2s_kgph = float(biogas.imass["HydrogenSulfide"]) if "HydrogenSulfide" in ids else 0.0
        biogas_mass_kgph = float(biogas.F_mass)
        biogas_Nm3ph = float(biogas.F_mol) * NM3_PER_KMOL
        digestate_mass_kgph = float(digestate.F_mass)

        self.design_results["AD inlet TS (kg/hr)"] = TS_in
        self.design_results["AD inlet VS (kg/hr)"] = VS_in
        self.design_results["Biodegradable pool (kg/hr)"] = biodegradable_pool
        self.design_results["Biodegradable destroyed (kg/hr)"] = biodegradable_destroyed
        self.design_results["Methane production (kg/hr)"] = methane_kgph
        self.design_results["CO2 production (kg/hr)"] = co2_kgph
        self.design_results["H2S production (kg/hr)"] = h2s_kgph
        self.design_results["Raw biogas total (kg/hr)"] = biogas_mass_kgph
        self.design_results["Raw biogas total (Nm3/hr)"] = biogas_Nm3ph
        self.design_results["Digestate outflow (kg/hr)"] = digestate_mass_kgph
        self.design_results["Organic loading proxy (kg/m3-h)"] = biodegradable_pool / V_liquid
        self.design_results["Methane productivity (kg/m3-h)"] = methane_kgph / V_liquid
        self.design_results["CH4 yield (kg/kg VS fed)"] = methane_kgph / VS_in

        self.power_utility.consumption = mixing_kW

        if Q_kJph > 0:
            try:
                self.add_heat_utility(Q_kJph, T_in, T_target)
            except TypeError:
                self.add_heat_utility(Q_kJph, T_in)

    def _cost(self):
        V_each = self.design_results["Digester volume each (m3)"]
        N = int(self.design_results["Number of digesters"])

        self.F_BM["Anaerobic digester"] = 1.0

        C_each = interp_ad_capex(V_each)
        C_total = N * C_each

        self.design_results["ADBC CAPEX each ($)"] = C_each
        self.design_results["ADBC CAPEX total ($)"] = C_total
        self.baseline_purchase_costs["Anaerobic digester"] = C_total

        if self.maintenance_usd_per_m3_yr is not None:
            total_volume = self.design_results["Total digester volume (m3)"]
            annual_maintenance = self.maintenance_usd_per_m3_yr * total_volume
            self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

            if annual_maintenance > 0:
                self.add_OPEX = {
                    "Digester maintenance": annual_maintenance / OPERATING_HOURS_PER_YEAR
                }


class AcidogenicDigester(bst.Unit):
    """
    Acidogenic / arrested AD unit for VFA platform modeling. Structurally
    parallel to `AnaerobicDigester` but produces a VFA-rich broth instead of biogas.
    There is no biodegradability weighting per component, all mass in `digestible_IDs`
    is treated as equally available for VS destruction.

    Parameters
    ----------
    ins : stream
        Feed slurry (or milled biomass stream, in integrated mode).
    outs : tuple[stream, stream]
        Offgas and acidogenic broth.
    vs_destruction : float
        Fraction of VS destroyed during digestion (no per-component
        weighting -- all mass in `digestible_IDs` is equally available).
    vfa_kg_per_kg_vs : float
        VFA yield per kg of VS destroyed.
    vfa_split : dict[str, float], optional
        Mass split of produced VFA among chemical IDs; renormalized to
        sum to 1. If not given, resolved at runtime by
        `_resolve_vfa_split()` (a thermo-dependent fallback).
    digestible_IDs : Iterable[str], optional
        Chemical IDs eligible for VS destruction. Falls back to a
        built-in default list if not given.
    produce_offgas_co2 : bool
        If True, VS destroyed but not converted to VFA is sent to the
        offgas as CarbonDioxide; if False (or CarbonDioxide absent from
        thermo), it's returned to the broth as Water instead.
    hrt_days : float
        Hydraulic retention time, used to size total liquid digester
        volume.
    slurry_density_kg_per_m3 : float
        Feed slurry density, used to convert mass flow to volumetric
        flow for sizing.
    headspace_frac : float
        Fraction of total digester volume reserved for gas headspace.
    max_single_digester_volume_MG : float
        Maximum volume (million US gallons) of a single digester train;
        larger required volumes are split across parallel trains.
    mixing_W_per_m3 : float
        Mixing power intensity per m3 of liquid volume.
    influent_temperature_K : float
        Feed slurry inlet temperature, used to size base heating duty.
    target_temperature_K : float
        Digester operating temperature, used to size base heating duty.
    cp_kJ_per_kgK : float
        Feed slurry specific heat capacity, used to size heating duty.
    enable_heat_shock : bool
        If True, adds a periodic heat-shock duty on top of the base
        heating duty.
    hs_target_temperature_K : float
        Peak temperature reached during each heat-shock event.
    hs_events_per_day : float
        Frequency of heat-shock events.
    hs_heated_fraction_of_liquid : float
        Fraction of the liquid volume heated during each event.
    hs_duration_min : float
        Duration of each heat-shock event; recorded in design_results,
        not used in the duty calculation itself.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/ad.yaml for the default values and references.
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
        chems = self.chemicals
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
        ids = set(self.chemicals.IDs)
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
            ids = set(self.chemicals.IDs)
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

        C_each = interp_ad_capex(V_each)
        C_total = N * C_each

        self.design_results["ADBC_capex_each_$"] = C_each
        self.baseline_purchase_costs["Acidogenic digester"] = C_total


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


    Parameters
    ----------
    ins : stream
        Raw biogas.
    outs : tuple[stream, stream]
        Treated biogas (H2S removed to near-zero) and spent media /
        captured H2S.
    h2s_removal_efficiency : float
        Fraction of inlet H2S removed (0-1).
    reagent_cost_usd_per_Nm3_raw : float
        Iron sponge media replacement cost per Nm3 of raw biogas
        treated, added to `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`. Installed cost itself is set
        by the class-level `@cost` decorator (power-law on raw biogas
        flow), not by an __init__ parameter.

    See Also
    --------
    Refer to data/ad.yaml for the default values and references.
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
            self.design_results["Reagent cost ($/yr)"] = reagent_usd_per_hr * OPERATING_HOURS_PER_YEAR


class BiogasUpgrading(bst.Unit):
    """
    Membrane biogas upgrading to pipeline-quality biomethane.

    Parameters
    ----------
    ins : stream
        Raw biogas (post H2S removal).
    outs : tuple[stream, stream]
        Biomethane and offgas.
    ch4_recovery : float
        Fraction of inlet methane recovered to the biomethane product
        (0-1).
    co2_removal : float
        Fraction of inlet CO2 removed to offgas (0-1).
    electricity_kwh_per_Nm3_raw : float
        Electricity intensity per Nm3 of dry raw biogas processed.
    capex_usd_per_Nm3ph_raw : float
        Installed cost per Nm3/h of dry raw biogas capacity.
    maintenance_frac_of_capex_per_yr : float
        Annual maintenance cost as a fraction of installed cost, added
        to `add_OPEX`.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/ad.yaml for the default values and references.
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
            self.add_OPEX = {
                "Membrane upgrading maintenance": annual_maintenance / OPERATING_HOURS_PER_YEAR
            }



class DigestateDecanterCentrifuge(bst.Unit):
    """
    Post-AD digestate decanter centrifuge (solid-liquid separation).

    An alternative to DigestateScrewPress for digestate dewatering: a decanter
    centrifuge achieves higher TS capture than a screw press (SYSTEMIC D3.2
    report: centrifuge ~59+-17% vs screw press ~33+-14%) at a different cost
    and electricity profile.

    Splits digestate into:
        - cake (soil_amendment): captured solids + entrained water to hit cake_moisture_frac
        - centrate (liquid_digestate): remaining liquid + uncaptured solids

    Assumptions:
        - "Solids" are explicitly defined chemical IDs
        - Everything else is treated as liquid and starts in centrate

    Sizing:
        - Throughput_tph = inlet_mass_kgph / 1000
        - N_parallel = ceil(Throughput_tph / capacity_tph_each)

    Costing (USD):
        - Total purchased cost = N_parallel * centrifuge_purchase_cost_usd_each

    Parameters
    ----------
    ins : stream
        AD digestate.
    outs : tuple[stream, stream]
        Cake (soil_amendment) and centrate (liquid_digestate).
    solids_IDs : Iterable[str]
        Chemical IDs treated as solids; only these are captured to cake
        (at `ts_capture_frac`), everything else goes entirely to
        centrate.
    ts_capture_frac : float
        Fraction of solids-ID mass captured to cake (0-1).
    cake_moisture_frac : float
        Target moisture fraction of the cake; entrained water is added
        to meet it.
    capacity_tph_each : float
        Throughput capacity (metric ton/h) of a single centrifuge;
        larger throughputs are split across parallel units.
    centrifuge_purchase_cost_usd_each : float, optional
        Purchase cost per centrifuge. If None, `_cost()` returns early
        and the unit is costed at $0.
    F_BM : float
        Bare-module cost factor applied to the centrifuge purchase
        cost.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/ad.yaml for the default values and references.
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


# DigestateScrewPress per-unit installed cost anchor points (SYSTEMIC Table
# 2-11, CAPEX vs ton/h), pre-converted from EUR to USD at a fixed 1.19
# EUR/USD rate.
SCREW_PRESS_CAPACITY_TPH = [2.0, 3.0, 4.0, 5.0, 6.25, 8.0, 9.0, 10.0, 12.0, 15.5]
SCREW_PRESS_CAPEX_USD = [
    30493.75,  # 25625.00 EUR at 2.00 tph
    43062.53,  # 36187.00 EUR at 3.00 tph
    36195.04,  # 30416.00 EUR at 4.00 tph
    47089.49,  # 39571.00 EUR at 5.00 tph
    34212.50,  # 28750.00 EUR at 6.25 tph
    41153.77,  # 34583.00 EUR at 8.00 tph
    52836.00,  # 44400.00 EUR at 9.00 tph
    28633.78,  # 24062.00 EUR at 10.00 tph
    62475.00,  # 52500.00 EUR at 12.00 tph
    20230.00,  # 17000.00 EUR at 15.50 tph
]


def interp_screw_press_capex_usd(capacity_tph: float) -> float:
    """
    Linear interpolation/extrapolation of DigestateScrewPress per-unit
    installed cost (USD) from SCREW_PRESS_CAPACITY_TPH/SCREW_PRESS_CAPEX_USD.
    """
    xs = SCREW_PRESS_CAPACITY_TPH
    ys = SCREW_PRESS_CAPEX_USD

    if capacity_tph <= xs[0]:
        return ys[0]
    if capacity_tph >= xs[-1]:
        return ys[-1]

    for i in range(len(xs) - 1):
        if xs[i] <= capacity_tph <= xs[i + 1]:
            t = (capacity_tph - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + t * (ys[i + 1] - ys[i])

    return ys[-1]


class DigestateScrewPress(bst.Unit):
    """
    Post-AD digestate screw press (solid-liquid separation).

    Splits digestate into:
        - cake (soil_amendment): captured solids + entrained water to hit cake_moisture_frac
        - pressate (liquid_digestate): remaining liquid + uncaptured solids

    Assumptions:
        - "Solids" are everything except Water and dissolved_IDs (dynamically computed)
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
        - Purchased cost per press interpolated from a fixed SYSTEMIC Table
          2-11 (CAPEX vs ton/h) anchor table, pre-converted to USD -- not a
          user-configurable input, see interp_screw_press_capex_usd() above
        - Total purchased cost = N_parallel * cost_per_press
        - Optional polymer dosing unit per press (~12k-50k EUR, pre-converted
          to USD) can be included


    Parameters
    ----------
    ins : stream
        AD digestate.
    outs : tuple[stream, stream]
        Cake (soil_amendment) and pressate (liquid_digestate).
    dissolved_IDs : Iterable[str], optional
        Chemical IDs always routed entirely to pressate regardless of
        `ts_capture_frac` (e.g. VFAs, which must not be captured in the
        cake). Falls back to the standard VFA set plus CarbonDioxide and
        Ammonia if not given.
    ts_capture_frac : float
        Fraction of dynamically-computed TS (everything except Water
        and `dissolved_IDs`) captured to cake (0-1).
    cake_moisture_frac : float
        Target moisture fraction of the cake; entrained water is added
        to meet it.
    capacity_tph_each : float
        Throughput capacity (metric ton/h) of a single press; larger
        throughputs are split across parallel units.
    kWh_per_m3 : float
        Electricity intensity per m3 of digestate treated.
    include_polymer_dosing : bool
        If True, adds an optional polymer dosing cost per press.
    polymer_dosing_cost_usd_each : float
        Polymer dosing cost per press, in USD; only applied if
        `include_polymer_dosing` is True.
    F_BM : float
        Bare-module cost factor applied to the press purchase cost.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/ad.yaml for the default values and references.
    """

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID="", ins=None, outs=(),
                 dissolved_IDs=None,            # chemicals treated as dissolved — always route to pressate
                 ts_capture_frac=0.40,          # assumptions.yaml digestate_screw_press.ts_capture_frac
                 cake_moisture_frac=0.45,       # assumptions.yaml digestate_screw_press.cake_moisture_frac
                 capacity_tph_each=6.0,         # aligns with reported 6 ton/h energy datapoint
                 kWh_per_m3=0.67,               # SYSTEMIC avg electricity intensity
                 include_polymer_dosing=False,
                 polymer_dosing_cost_usd_each=35000, 
                 F_BM=1.0,
                 **kwargs):
        super().__init__(ID, ins, outs, **kwargs)

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

        self.include_polymer_dosing = bool(include_polymer_dosing)
        self.polymer_dosing_cost_usd_each = float(polymer_dosing_cost_usd_each)

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

    def _cost(self):
        N = int(self.design_results.get("N_screw_presses_parallel", 1))
        cap_each = float(self.capacity_tph_each)

        capex_usd_each = interp_screw_press_capex_usd(cap_each)
        total_usd = N * capex_usd_each

        if self.include_polymer_dosing:
            dosing_usd_each = self.polymer_dosing_cost_usd_each
            total_usd += N * dosing_usd_each

        key = "Screw press (digestate)"
        self.baseline_purchase_costs[key] = total_usd
        self.F_BM[key] = self.F_BM_default

        self.design_results["Capex_usd_each"] = capex_usd_each
        if self.include_polymer_dosing:
            self.design_results["Polymer_dosing_cost_usd_each"] = dosing_usd_each
        self.design_results["N_parallel"] = N
        self.design_results["Total_purchase_cost_usd"] = total_usd

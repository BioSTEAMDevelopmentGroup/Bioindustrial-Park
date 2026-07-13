# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Post-AD digestate screw press (solid-liquid separation)

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
"""

import math
import biosteam as bst

__all__ = ('DigestateScrewPress',)


class DigestateScrewPress(bst.Unit):
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID="", ins=None, outs=(),
                 solids_IDs=("Cellulose", "Ash"),
                 dissolved_IDs=None,            # chemicals treated as dissolved — always route to pressate
                 ts_capture_frac=0.33,          # SYSTEMIC avg SE_DM for screw press
                 cake_moisture_frac=0.77,       # ~23% DM solid fraction
                 capacity_tph_each=6.0,         # aligns with reported 6 ton/h energy datapoint
                 kWh_per_m3=0.67,               # SYSTEMIC avg electricity intensity
                 capex_eur_table=None,          # override if you want your own table
                 eur_to_usd=1.0,                # set in your config if you want USD
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

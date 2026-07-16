# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Post-AD digestate decanter centrifuge (solid-liquid separation)

An alternative to DigestateScrewPress for digestate dewatering: a decanter
centrifuge achieves higher TS capture than a screw press (SYSTEMIC D3.2
report: centrifuge ~59+-17% vs screw press ~33+-14%) at a different cost
and electricity profile. Not wired into any of the default system builders
today — kept available as a selectable alternative, same as the
pretreatment options in systems._ad_biogas_system.

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
"""

import math
import biosteam as bst

__all__ = ('DigestateDecanterCentrifuge',)


class DigestateDecanterCentrifuge(bst.Unit):
    """
    Alternative digestate dewatering unit, not currently constructed by any
    system builder in this codebase (confirmed by grep -- only this file and
    units/__init__.py reference the class name) -- so every parameter below
    is only ever exercised via its own __init__ default, never overridden by
    a real call site. Unlike `DigestateScrewPress`, `solids_IDs` here is an
    active whitelist (see _run): only chemicals in this list are captured to
    cake at `ts_capture_frac`; everything else goes entirely to centrate.

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

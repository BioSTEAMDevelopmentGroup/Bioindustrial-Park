# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Pressate membrane concentrator for the biostimulant side-stream.
"""
from __future__ import annotations

import biosteam as bst

__all__ = ('PressateConcentrator',)


class PressateConcentrator(bst.Unit):
    """
    Membrane concentrator for the pressate-derived biostimulant side-stream.
    Constructed in three places in this codebase
    (`systems/_ad_biogas_system.py`, `systems/_vfa_ad_system.py`,
    `systems/_integrated_system.py`), all with the same
    `assumptions.yaml`-sourced parameters.

    Kept as custom _design/_cost rather than an @cost decorator, same
    reasoning as `BiogasUpgrading`: annual maintenance OPEX here is a
    fraction of this unit's own installed cost (`capex_usd_per_m2 *
    membrane_area_m2`), computed inside `_cost()` itself -- @cost's
    decorator-owned `_cost()` only computes that value after `_design()`
    has already run, so there's no clean way to read it early without
    duplicating the cost formula.

    Sources
    -------
    retained_solute_IDs, retained_solute_recovery_to_concentrate (0.95),
    design_flux_L_m2_h (35.0), operating_pressure_bar (5.0):
        assumptions.yaml `pressate_biostimulant.concentrator` section.
        `design_flux_L_m2_h` sourced to `sources.flux`: "Seaweed/algal
        membrane concentration literature; Sievers et al. (2017);
        Diaz-Reinoso and Dominguez (2020) -- UF fluxes for bioprocess and
        seaweed-derived streams commonly fall in the ~15-60 LMH range; 35
        LMH selected as a central design value." This is the citation that
        was previously mis-attributed to `VFAMicrofilter` during the VFA
        fermentation conversion (see the "cross-contaminated citation" item
        earlier in this file) -- it belongs here, on this unit's matching
        35 LMH value.
    water_recovery_to_permeate (0.93), electricity_kWh_per_m3_feed (2.5):
        assumptions.yaml `pressate_biostimulant.concentrator` section.
        `electricity_kWh_per_m3_feed` sourced to `sources.electricity`:
        "Membrane concentration literature for algal and bioprocess
        streams -- UF electricity commonly falls in the ~0.5-2.5 kWh/m3
        range."
    capex_usd_per_m2 (500.0):
        assumptions.yaml `pressate_biostimulant.concentrator.
        capex_usd_per_m2`. UNVERIFIED CITATION -- assumptions.yaml's
        `sources.capex` block has duplicate `citation`/`note` keys (a yaml
        authoring bug): the raw file text shows "Sethi and Wiesner (2000);
        Nguyen and Yoshikawa (2020); Cheryan (1998)" with the note "UF
        membrane/module cost commonly reported around $200-500/m2; $800/m2
        for conservative assumption" -- immediately followed by a second
        citation/note pair about maintenance fraction, using the exact
        same key names. Standard YAML parsing keeps only the *last* of
        duplicate keys, so `yaml.safe_load()` on this file actually
        returns the maintenance-fraction text for `sources.capex`, not the
        Sethi/Wiesner/Nguyen/Cheryan citation -- confirmed directly in
        this conversion by parsing the file, not just reading it. The
        modeled value (500.0) does fall inside the range the orphaned
        note text states ($200-500/m2), not the "$800/m2 conservative"
        figure also in that same orphaned text. Not fixable here (yaml is
        out of scope for this conversion); flagged in CONVERSION_NOTES.md.
    maintenance_frac_of_capex_per_yr (0.035):
        assumptions.yaml `pressate_biostimulant.concentrator.
        maintenance_frac_of_capex_per_yr`. Its own citation is the same
        orphaned/overwritten yaml text described above.
    Non-target-solute permeate/concentrate split (hardcoded 0.70 inline in
    _run(), not an __init__ param):
        No assumptions.yaml counterpart at all -- an internal scoping
        assumption for chemicals outside `retained_solute_IDs`, not
        corrected here since there is no yaml value to align it to.
    """

    _N_ins = 1
    _N_outs = 2  # concentrated_biostimulant, permeate

    def __init__(
        self,
        ID="",
        ins=None,
        outs=(),
        retained_solute_IDs=(
            "Alginate",
            "Fucoidan",
            "Mannitol",
            "Protein",
            "OtherSolids",
        ),
        water_recovery_to_permeate=0.93,
        retained_solute_recovery_to_concentrate=0.95,
        design_flux_L_m2_h=35.0,
        operating_pressure_bar=5.0,
        electricity_kWh_per_m3_feed=2.5,
        capex_usd_per_m2=500.0,
        maintenance_frac_of_capex_per_yr=0.035,
        F_BM=1.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)

        self.retained_solute_IDs = tuple(retained_solute_IDs)
        self.water_recovery_to_permeate = float(water_recovery_to_permeate)
        self.retained_solute_recovery_to_concentrate = float(
            retained_solute_recovery_to_concentrate
        )
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)
        self.operating_pressure_bar = float(operating_pressure_bar)
        self.electricity_kWh_per_m3_feed = float(electricity_kWh_per_m3_feed)
        self.capex_usd_per_m2 = float(capex_usd_per_m2)
        self.maintenance_frac_of_capex_per_yr = float(
            maintenance_frac_of_capex_per_yr
        )

        self.F_BM = dict(getattr(self, "F_BM", {}))
        self.F_BM["Pressate concentrator"] = float(F_BM)

    def _run(self):
        feed = self.ins[0]
        conc, perm = self.outs

        conc.empty()
        perm.empty()

        conc.phase = "l"
        perm.phase = "l"

        ids = set(feed.chemicals.IDs)

        water = float(feed.imass["Water"]) if "Water" in ids else 0.0
        water_to_perm = self.water_recovery_to_permeate * water
        water_to_conc = water - water_to_perm

        if "Water" in ids:
            conc.imass["Water"] = water_to_conc
            perm.imass["Water"] = water_to_perm

        retained = set(self.retained_solute_IDs)

        for cid in feed.chemicals.IDs:
            if cid == "Water":
                continue

            m = float(feed.imass[cid])
            if m <= 0:
                continue

            if cid in retained:
                m_to_conc = self.retained_solute_recovery_to_concentrate * m
                m_to_perm = m - m_to_conc
            else:
                # base case: non-target solubles mostly follow permeate/water
                m_to_perm = 0.70 * m
                m_to_conc = m - m_to_perm

            conc.imass[cid] = m_to_conc
            perm.imass[cid] = m_to_perm

        self._feed_m3ph = float(feed.F_vol)
        self._conc_m3ph = float(conc.F_vol)
        self._perm_m3ph = float(perm.F_vol)

    def _design(self):
        feed = self.ins[0]

        feed_m3ph = float(feed.F_vol)
        flux_m3_m2_h = self.design_flux_L_m2_h / 1000.0  # LMH -> m3/m2/h
        membrane_area_m2 = feed_m3ph / flux_m3_m2_h if flux_m3_m2_h > 0 else 0.0

        power_kW = self.electricity_kWh_per_m3_feed * feed_m3ph

        self.design_results["Feed flow (m3/h)"] = feed_m3ph
        self.design_results["Design flux (L/m2-h)"] = self.design_flux_L_m2_h
        self.design_results["Operating pressure (bar)"] = self.operating_pressure_bar
        self.design_results["Membrane area (m2)"] = membrane_area_m2
        self.design_results["Electricity demand (kW)"] = power_kW

        self.power_utility(power_kW)

    def _cost(self):
        membrane_area_m2 = self.design_results.get("Membrane area (m2)", 0.0)
        capex = membrane_area_m2 * self.capex_usd_per_m2

        self.baseline_purchase_costs["Pressate concentrator"] = capex

        # Annual maintenance -> add_OPEX so it flows through
        # SABREBaselineTEA.VOC -> solve_product_msp correctly
        annual_maintenance = self.maintenance_frac_of_capex_per_yr * capex
        self.design_results["Annual maintenance ($/yr)"] = annual_maintenance

        if annual_maintenance > 0:
            operating_hours_per_yr = 330.0 * 24.0
            self.add_OPEX = {
                "Pressate concentrator maintenance": annual_maintenance / operating_hours_per_yr
            }

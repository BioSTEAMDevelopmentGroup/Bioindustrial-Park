# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
H2S Removal Unit
==========================================
Removes hydrogen sulfide from raw biogas before membrane upgrading.
H2S must be removed because it degrades polymer membranes rapidly and
is toxic above ~200 ppm in pipeline-quality biomethane

Technology: Iron sponge (dry oxidation) — the standard, lowest-cost
option for biogas desulfurization at this scale

Modeling approach:
- Pass-through for CH4 and CO2 (not affected by iron sponge)
- H2S captured to near-zero in treated gas
- Capital cost scaled to raw biogas flow (Nm3/h)
- Reagent cost (iron sponge replacement) as add_OPEX
"""

from __future__ import annotations
import biosteam as bst

__all__ = ('H2SRemoval',)


class H2SRemoval(bst.Unit):
    """
    Iron sponge H2S removal unit for raw biogas desulfurization.

    Inputs:
        ins[0]: raw biogas (CH4 + CO2 + H2S)

    Outputs:
        outs[0]: treated biogas (H2S removed to near-zero)
        outs[1]: spent media / captured H2S (solid waste, negligible mass)
    """

    _N_ins = 1
    _N_outs = 2   # treated_biogas, spent_media
    _F_BM_default = {"H2S removal (iron sponge)": 1.0}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        h2s_removal_efficiency: float = 0.99,
        ref_flow_Nm3ph: float = 1700.0,
        ref_installed_cost_usd: float = 450_000.0,
        scale_exponent: float = 0.7,
        reagent_cost_usd_per_Nm3_raw: float = 0.005,
        F_BM: float = 1.0,
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.h2s_removal_efficiency = float(h2s_removal_efficiency)
        self.ref_flow_Nm3ph = float(ref_flow_Nm3ph)
        self.ref_installed_cost_usd = float(ref_installed_cost_usd)
        self.scale_exponent = float(scale_exponent)
        self.reagent_cost_usd_per_Nm3_raw = float(reagent_cost_usd_per_Nm3_raw)
        self.F_BM = {"H2S removal (iron sponge)": float(F_BM)}

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

    def _cost(self):
        Q_Nm3ph = self.design_results.get("Raw biogas flow (Nm3/h)", 0.0)
        ref = self.ref_flow_Nm3ph
        ref_cost = self.ref_installed_cost_usd

        if Q_Nm3ph <= 0 or ref <= 0:
            installed_cost = ref_cost
        else:
            installed_cost = ref_cost * (Q_Nm3ph / ref) ** self.scale_exponent

        self.baseline_purchase_costs["H2S removal (iron sponge)"] = installed_cost
        self.design_results["Installed cost ($)"] = installed_cost

        # Reagent (iron sponge media replacement) as add_OPEX
        # Basis: $0.002/Nm3 raw biogas
        reagent_usd_per_hr = self.reagent_cost_usd_per_Nm3_raw * Q_Nm3ph
        if reagent_usd_per_hr > 0:
            self.add_OPEX = {
                "Iron sponge media replacement": reagent_usd_per_hr
            }
            self.design_results["Reagent cost ($/hr)"] = reagent_usd_per_hr
            self.design_results["Reagent cost ($/yr)"] = reagent_usd_per_hr * 330.0 * 24.0

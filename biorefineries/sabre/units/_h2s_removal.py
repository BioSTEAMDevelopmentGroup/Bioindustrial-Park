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
- Capital cost scaled to raw biogas flow (Nm3/h) via a power-law @cost item
- Reagent cost (iron sponge replacement) as add_OPEX, computed in _design
  since it depends on raw biogas flow, not on the unit's own installed cost
"""

from __future__ import annotations
import biosteam as bst
from biosteam.units.decorators import cost

__all__ = ('H2SRemoval',)


@cost('Raw biogas flow (Nm3/h)', 'H2S removal (iron sponge)', units='Nm3/h',
      CE=567.5, cost=450_000.0, S=1700.0, n=0.7, BM=1.0)
class H2SRemoval(bst.Unit):
    """
    Iron sponge H2S removal unit for raw biogas desulfurization.

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

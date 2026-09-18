# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam.units.decorators import cost

from biorefineries.sabre.utils import load_assumptions, get_solids_group_IDs

__all__ = ('VFAMicrofilter',)

# Loaded yaml assumptions (vfa_microfilter now lives in vfa.yaml;
# downstream_processing.yaml was split into vfa.yaml/microbial_oil.yaml).
_VFA_YAML = load_assumptions("vfa.yaml")
_VFA_MICROFILTER = _VFA_YAML["vfa_microfilter"]


@cost('Membrane area (m2)', 'Microfilter', units='m2',
      CE=567.5, cost=_VFA_MICROFILTER["membrane_cost_usd_per_m2"], S=1., n=1., BM=_VFA_MICROFILTER["F_BM"])
class VFAMicrofilter(bst.Unit):
    """
    Split-based representation of a VFA-rich permeate step.
    Includes first-pass power draw and area-based membrane cost.

    Parameters
    ----------
    ins : stream
        VFA-rich broth feed.
    outs : tuple[stream, stream]
        Permeate and retentate.
    vfa_IDs : Iterable[str]
        Chemical IDs treated as VFA.
    solids_IDs : Iterable[str]
        Chemical IDs treated as solids for the split. Defaults to the
        "solids" chemical group when not given (see
        `utils.get_solids_group_IDs`).
    vfa_to_permeate_frac : float
        Fraction of VFA mass routed to permeate (0-1).
    water_to_permeate_frac : float
        Fraction of water mass routed to permeate (0-1).
    solids_to_permeate_frac : float
        Fraction of solids-ID mass routed to permeate (0-1).
    dissolved_other_to_permeate_frac : float
        Fraction of remaining dissolved mass routed to permeate (0-1).
    SEC_kWh_per_m3_feed : float
        Specific electricity consumption per m3 of feed.
    design_flux_L_m2_h : float
        Design membrane flux, used to size membrane area.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/downstream_processing.yaml for the default values and references.
    """
    _N_ins = 1
    _N_outs = 2
    _units = {'Membrane area (m2)': 'm2'}

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        vfa_IDs: list[str] = _VFA_MICROFILTER["vfa_IDs"],
        solids_IDs=None,  # if not given, defaults to the "solids" chemical group (see utils.get_solids_group_IDs)
        vfa_to_permeate_frac: float = _VFA_MICROFILTER["vfa_to_permeate_frac"],
        water_to_permeate_frac: float = _VFA_MICROFILTER["water_to_permeate_frac"],
        solids_to_permeate_frac: float = _VFA_MICROFILTER["solids_to_permeate_frac"],
        dissolved_other_to_permeate_frac: float = _VFA_MICROFILTER["dissolved_other_to_permeate_frac"],
        SEC_kWh_per_m3_feed: float = _VFA_MICROFILTER["SEC_kWh_per_m3_feed"],
        design_flux_L_m2_h: float = _VFA_MICROFILTER["design_flux_L_m2_h"],
        **kwargs,
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.vfa_IDs = tuple(vfa_IDs)
        if solids_IDs is None:
            solids_IDs = get_solids_group_IDs(self.chemicals)
        self.solids_IDs = tuple(solids_IDs)
        self.vfa_to_permeate_frac = float(vfa_to_permeate_frac)
        self.water_to_permeate_frac = float(water_to_permeate_frac)
        self.solids_to_permeate_frac = float(solids_to_permeate_frac)
        self.dissolved_other_to_permeate_frac = float(dissolved_other_to_permeate_frac)
        self.SEC_kWh_per_m3_feed = float(SEC_kWh_per_m3_feed)
        self.design_flux_L_m2_h = float(design_flux_L_m2_h)

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
        feed_m3h = feed.F_vol
        membrane_area_m2 = 0.0
        if self.design_flux_L_m2_h > 0:
            membrane_area_m2 = feed_m3h * 1000.0 / self.design_flux_L_m2_h

        self.design_results["Feed flow (kg/h)"] = feed.F_mass
        self.design_results["Feed flow (m3/h)"] = feed_m3h
        self.design_results["Permeate flow (kg/h)"] = self.outs[0].F_mass
        self.design_results["Retentate flow (kg/h)"] = self.outs[1].F_mass
        self.design_results["Membrane area (m2)"] = membrane_area_m2
        self.power_utility(self.SEC_kWh_per_m3_feed * feed_m3h)

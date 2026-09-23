# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
import thermosteam as tmo
from biosteam.units.design_tools import PressureVessel, cylinder_diameter_from_volume

from biorefineries.sabre.utils import load_assumptions

__all__ = ('DiluteAcidReactor',)

_M_TO_FT = 3.28084
_PA_TO_PSI = 0.000145038

# Loaded assumptions
_PRETREATMENT = load_assumptions("3hp.yaml")["cake_dilute_acid_pretreatment"]


class DiluteAcidReactor(PressureVessel, bst.Unit):
    """
    Dilute-acid pretreatment pressure reactor for the 3-HP cake
    configurations: a sabre-local adaptation of
    biorefineries.cellulosic's `PretreatmentReactorSystem` (corn stover),
    re-implemented here rather than imported for the same reason as
    `HPFermentation` (sabre's dependency boundary, and cellulosic's reaction
    set is written for lignocellulose chemistry that sabre's chemical set
    does not carry).

    The feed is the slurry leaving the steam mixer, already at reaction
    temperature and pressure, so the reactor is isothermal: it applies the
    hydrolysis reactions and sizes the pressure vessel from the residence
    time. Both reactions are exactly atom-balanced:

    - `Glucan + Water -> Glucose`
    - `Alginate + Water -> AlginateMonomer`

    Degradation products (HMF, furfural) are not modeled: they are not in
    sabre's chemical set.

    Parameters
    ----------
    ins : stream
        Steam-mixed cake slurry, already at reaction T and P.
    outs : stream
        Pretreated slurry.
    glucan_conversion : float
        Fraction of Glucan hydrolyzed to Glucose.
    alginate_conversion : float
        Fraction of Alginate hydrolyzed to AlginateMonomer.
    tau : float
        Residence time [hr].
    V_wf : float
        Working volume fraction of the vessel.
    length_to_diameter : float
        Vessel length-to-diameter ratio.
    vessel_material : str
        Pressure vessel material.
    vessel_type : str
        'Horizontal' or 'Vertical'.
    **kwargs
        Forwarded to `bst.Unit.__init__`.

    See Also
    --------
    Refer to data/3hp.yaml (`cake_dilute_acid_pretreatment`) for the default
    values and references.
    """

    _N_ins = 1
    _N_outs = 1
    _units = {**PressureVessel._units,
              'Residence time': 'hr',
              'Reactor volume': 'm3'}

    def __init__(
        self, ID="", ins=None, outs=(),
        glucan_conversion=_PRETREATMENT["glucan_to_glucose_conversion"],
        alginate_conversion=_PRETREATMENT["alginate_to_alginate_monomer_conversion"],
        tau=_PRETREATMENT["reactor_tau_hr"],
        V_wf=_PRETREATMENT["reactor_V_wf"],
        length_to_diameter=_PRETREATMENT["reactor_length_to_diameter"],
        vessel_material=_PRETREATMENT["reactor_vessel_material"],
        vessel_type=_PRETREATMENT["reactor_vessel_type"],
        **kwargs
    ):
        super().__init__(ID, ins, outs, **kwargs)
        self.glucan_conversion = float(glucan_conversion)
        self.alginate_conversion = float(alginate_conversion)
        self.tau = float(tau)
        self.V_wf = float(V_wf)
        self.length_to_diameter = float(length_to_diameter)
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

        self.reactions = tmo.ParallelReaction([
            tmo.Reaction('Glucan + Water -> Glucose', reactant='Glucan',
                         X=self.glucan_conversion),
            tmo.Reaction('Alginate + Water -> AlginateMonomer', reactant='Alginate',
                         X=self.alginate_conversion),
        ])

    def _run(self):
        feed, = self.ins
        effluent, = self.outs
        effluent.copy_like(feed)
        self.reactions(effluent)

    def _design(self):
        V_reactor = self.F_vol_in * self.tau / self.V_wf  # m3
        diameter = cylinder_diameter_from_volume(V_reactor, self.length_to_diameter) * _M_TO_FT
        length = diameter * self.length_to_diameter
        pressure = self.outs[0].P * _PA_TO_PSI

        design = self.design_results
        design['Residence time'] = self.tau
        design['Reactor volume'] = V_reactor
        design.update(self._vessel_design(float(pressure), float(diameter), float(length)))

    def _cost(self):
        design = self.design_results
        self.baseline_purchase_costs.update(
            self._vessel_purchase_cost(design['Weight'], design['Diameter'], design['Length'])
        )

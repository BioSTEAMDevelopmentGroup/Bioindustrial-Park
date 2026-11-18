# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
import thermosteam as tmo
from biosteam.units import BatchBioreactor

from biorefineries.sabre.utils import load_assumptions

__all__ = ('HPFermentation', 'CaHPCrystallizer', 'CrystalCentrifuge', 'HydrolysisTank')

# Loaded assumptions
_3HP_YAML = load_assumptions("3hp.yaml")
_FERMENTATION = _3HP_YAML["fermentation"]
_SUBSTRATES = _FERMENTATION["substrates"]
_CRYSTALLIZATION = _3HP_YAML["crystallization"]
_CRYSTAL_SEPARATOR = _3HP_YAML["crystal_separator"]


class HPFermentation(BatchBioreactor):
    """
    Sargassum-derived-sugar fermentation to 3-hydroxypropionic acid (HP),
    neutralized in situ with Ca(OH)2 to calcium 3-hydroxypropionate
    (Ca3HP2). A sabre-local adaptation of biorefineries.HP's
    BatchCoFermentation (`biorefineries/HP/units.py`), re-implemented here
    rather than imported, to keep sabre's dependency boundary clean.

    For each substrate in `substrates`, two independent reactions draw from
    the same original available pool (mirroring how HP's own
    ParallelReaction handles "Glucose -> 2 HP" and "Glucose -> 6 FermMicrobe
    + 2.4 H2O" as two separate reactions sharing the Glucose reactant):

    - HP formation: `hp_conversion` fraction of available substrate is
      consumed, forming HP at `hp_yield_kg_per_kg_consumed` (mass HP per
      mass substrate consumed by this reaction).
    - Biomass formation: `biomass_conversion` fraction of available
      substrate is independently consumed, forming CellMass at
      `biomass_yield_kg_per_kg_consumed`.

    Glucose -> 2 HP is exactly mass-balanced (2x HP's MW == Glucose's MW),
    but Mannitol and Uronics are not exact multiples of HP's MW.
    thermosteam.Reaction does not enforce or warn on an unbalanced
    reaction -- confirmed by direct testing that it silently destroys real
    mass -- so all three substrates are handled with the same imperative,
    mass-conserving-by-construction arithmetic (consumed * yield) instead
    of thermosteam Reaction objects, mirroring sabre's own
    YarrowiaLipidFermenter._run_reactions pattern. This is the "mass-yield
    scaling against MW" simplification flagged in the design spec (secs.
    7, 12).

    Neutralization (`2 HP + CalciumDihydroxide -> Ca3HP2 + 2 Water`) IS an
    exactly atom-balanced reaction (LHS/RHS MW both 254.24856), so it is
    implemented as a real `tmo.Reaction`, dosed to the exact stoichiometric
    demand of the HP formed plus a small excess.

    Parameters
    ----------
    ins : stream
        Substrate feed, and lime (Ca(OH)2, dosed automatically -- pass an
        empty placeholder stream).
    outs : tuple[stream, stream]
        Vent and fermentation broth (Ca3HP2 already formed in situ).
    substrates : dict[str, dict]
        Per-substrate `hp_conversion`, `hp_yield_kg_per_kg_consumed`,
        `biomass_conversion`, `biomass_yield_kg_per_kg_consumed`.
    product_ID : str
        Chemical ID of the free-acid fermentation product (HP) -- exists
        only transiently within `_run`, since neutralization converts
        (almost) all of it to `product_salt_ID` before `_run` returns.
    cellmass_ID : str
        Chemical ID of cell mass.
    lime_ID : str
        Chemical ID of the neutralization base (Ca(OH)2).
    product_salt_ID : str
        Chemical ID of the neutralized product (Ca3HP2).
    neutralization_conversion : float
        Fraction of HP neutralized to the calcium salt (bounded < 1, since
        `tmo.Reaction` requires X < 1).
    lime_excess_frac : float
        Excess lime dosed beyond the exact stoichiometric demand.
    productivity_g_per_L_per_h : float
        Target volumetric productivity. `tau` (batch time) is computed in
        `_run` as the simulated titer (HP formed per effluent volume, i.e.
        HP-equivalent since neutralization to Ca3HP2 is ~complete) divided
        by this productivity, so it always tracks the actual broth
        concentration rather than an independently-specified titer.
    T : float
        Operating temperature [K] (`BatchBioreactor`'s own attribute name).
    P : float
        Operating pressure [Pa] (`BatchBioreactor`'s own attribute name).
    V : float
        Reactor volume [m3], passed to `BatchBioreactor` as `V_max`.
    **kwargs
        Forwarded to `biosteam.units.BatchBioreactor.__init__`.

    See Also
    --------
    Refer to data/3hp.yaml for the default values and references.
    """

    _N_ins = 2
    _N_outs = 2

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        substrates: dict = _SUBSTRATES,
        product_ID: str = _FERMENTATION["product_ID"],
        cellmass_ID: str = _FERMENTATION["cellmass_ID"],
        lime_ID: str = _FERMENTATION["lime_ID"],
        product_salt_ID: str = _FERMENTATION["product_salt_ID"],
        neutralization_conversion: float = _FERMENTATION["neutralization_conversion"],
        lime_excess_frac: float = _FERMENTATION["lime_excess_frac"],
        productivity_g_per_L_per_h: float = _FERMENTATION["productivity_g_per_L_per_h"],
        T: float = _FERMENTATION["T_K"],
        P: float = _FERMENTATION["P_Pa"],
        V: float = _FERMENTATION["V_m3"],
        **kwargs,
    ):
        # biosteam's batch bioreactor takes its reactor volume as `V_max`
        # (it has no `V` argument); `V` is kept as this unit's own name.
        kwargs.setdefault("V_max", V)
        super().__init__(ID, ins, outs, T=T, P=P, tau=None, **kwargs)

        self.substrates = {k: dict(v) for k, v in substrates.items()}
        self.product_ID = product_ID
        self.cellmass_ID = cellmass_ID
        self.lime_ID = lime_ID
        self.product_salt_ID = product_salt_ID
        self.neutralization_conversion = float(neutralization_conversion)
        self.lime_excess_frac = float(lime_excess_frac)
        self.productivity_g_per_L_per_h = float(productivity_g_per_L_per_h)

        self.neutralization_rxn = tmo.Reaction(
            f"2 {product_ID} + {lime_ID} -> {product_salt_ID} + 2 Water",
            reactant=product_ID,
            X=self.neutralization_conversion,
        )

    def _run(self):
        feed, lime = self.ins
        vent, effluent = self.outs

        lime.empty()
        vent.empty()

        effluent.copy_like(feed)
        effluent.phase = "l"

        self._run_reactions(effluent)
        self._run_neutralization(effluent, lime)

        effluent.T = vent.T = self.T
        effluent.P = vent.P = self.P

        # Batch time from the simulated titer (kg HP/h over m3 broth/h =
        # g/L), not an independently-specified one -- see the class
        # docstring's `productivity_g_per_L_per_h` entry.
        titer_g_per_L = self.design_results["HP formed (kg/h)"] / effluent.F_vol
        self.design_results["Titer (g HP-equivalent/L)"] = titer_g_per_L
        self.tau = titer_g_per_L / self.productivity_g_per_L_per_h

    def _run_reactions(self, effluent):
        ids = set(self.chemicals.IDs)
        required = (self.product_ID, self.cellmass_ID)
        missing = [i for i in required if i not in ids]
        if missing:
            raise RuntimeError(f"Missing required chemicals in thermo: {missing}")

        hp_formed_total = 0.0
        biomass_formed_total = 0.0
        substrate_available_total = 0.0

        for substrate_ID, params in self.substrates.items():
            if substrate_ID not in ids:
                continue
            available = float(effluent.imass[substrate_ID])
            if available <= 1e-12:
                continue
            substrate_available_total += available

            hp_consumed = params["hp_conversion"] * available
            hp_formed = params["hp_yield_kg_per_kg_consumed"] * hp_consumed

            biomass_consumed = params["biomass_conversion"] * available
            biomass_formed = params["biomass_yield_kg_per_kg_consumed"] * biomass_consumed

            # Mass-conserving by construction (mirrors sabre's own
            # YarrowiaLipidFermenter._run_reactions "yields don't sum to 1"
            # pattern): only the mass actually formed (hp_formed +
            # biomass_formed) is removed from the substrate pool, not the
            # nominal *_consumed amounts -- the gap between "consumed" (via
            # conversion X) and "accounted" (via the <1 mass yield) simply
            # stays as unreacted substrate rather than vanishing.
            accounted = hp_formed + biomass_formed
            effluent.imass[substrate_ID] -= accounted

            hp_formed_total += hp_formed
            biomass_formed_total += biomass_formed

        effluent.imass[self.product_ID] += hp_formed_total
        effluent.imass[self.cellmass_ID] += biomass_formed_total

        self.design_results["Substrate available (kg/h)"] = substrate_available_total
        self.design_results["HP formed (kg/h)"] = hp_formed_total
        self.design_results["Biomass formed (kg/h)"] = biomass_formed_total

    def _run_neutralization(self, effluent, lime):
        rxn = self.neutralization_rxn
        hp_mol = float(effluent.imol[self.product_ID])
        lime_mol_needed = 0.5 * hp_mol * (1.0 + self.lime_excess_frac)
        lime.imol[self.lime_ID] = lime_mol_needed
        effluent.mol += lime.mol

        rxn(effluent)

        self.design_results["Lime dosed (kg/h)"] = lime.F_mass


class CaHPCrystallizer(bst.BatchCrystallizer):
    """
    Batch crystallizer for calcium 3-hydroxypropionate (Ca3HP2), following
    the structural pattern of biorefineries.succinic's
    SuccinicAcidCrystallizer / biorefineries.TAL's TALCrystallizer (a
    single downstream `bst.units.SolidsCentrifuge` does the actual
    solid/liquid mechanical separation), but using a constant
    `target_recovery` rather than a temperature-solubility correlation:
    the patent reports fixed recovery/purity at specified conditions (room
    temperature, ~300 rpm stirring), not a solubility-vs-temperature
    curve, so a fixed-recovery split is the fidelity level the design spec
    (sec. 8) calls for. Inherits `bst.BatchCrystallizer`'s batch-vessel
    sizing/costing unchanged.

    Unlike those two references, this unit's outlet is NOT tagged as a
    genuine 2-phase ('l', 's') MultiStream -- confirmed by direct testing
    that doing so propagates badly through the downstream
    `bst.units.SolidsCentrifuge` -> `bst.units.DrumDryer` chain (the
    dryer's own internal 'g'-phase auxiliary streams raise
    `UndefinedPhase` when fed a MultiStream upstream of them). The actual
    separation is enforced entirely by the downstream `CrystalCentrifuge`'s
    own `product_recovery` (set to this unit's `target_recovery` by whoever
    wires the system script -- see data/3hp.yaml `crystal_separator`'s
    comment); `target_recovery` here only determines what this unit
    reports in `design_results`.

    Parameters
    ----------
    ins : stream
        Concentrated broth (product dissolved, from the upstream
        evaporator).
    outs : stream
        Effluent (same composition as the feed; the split downstream is
        not yet applied here -- see Notes above).
    product_ID : str
        Chemical ID of the crystallized product (Ca3HP2).
    target_recovery : float
        Fixed fraction of dissolved product mass recovered to the solid
        phase, reported in `design_results` only (see Notes above).
    **kwargs
        Forwarded to `bst.BatchCrystallizer.__init__`.

    See Also
    --------
    Refer to data/3hp.yaml for the default values and references.
    """

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        product_ID: str = _CRYSTALLIZATION["product_ID"],
        target_recovery: float = _CRYSTALLIZATION["target_recovery"],
        T: float = _CRYSTALLIZATION["T_K"],
        tau: float = _CRYSTALLIZATION["tau"],
        V: float = _CRYSTALLIZATION["V_m3"],
        **kwargs,
    ):
        kwargs.setdefault("V", V)
        super().__init__(ID, ins, outs, T=T, tau=tau, **kwargs)
        self.product_ID = product_ID
        self.target_recovery = float(target_recovery)

    def _run(self):
        feed, = self.ins
        effluent, = self.outs

        # Kept as a single-phase stream rather than a genuine 2-phase
        # ('l', 's') MultiStream (unlike succinic/TAL's crystallizers,
        # which this unit otherwise mirrors structurally): confirmed by
        # direct testing that a 2-phase MultiStream leaving this unit
        # propagates badly through the downstream
        # bst.units.SolidsCentrifuge -> bst.units.DrumDryer chain (the
        # dryer's own internal 'g'-phase auxiliary streams raise
        # UndefinedPhase when fed a MultiStream upstream of them). The
        # actual solid/liquid separation is enforced downstream by
        # crystal_separator's own `split` fraction (see data/3hp.yaml
        # crystal_separator's comment), not by phase tags on this stream,
        # so nothing physical is lost by not tagging phases here --
        # target_recovery is still reported below for documentation.
        effluent.copy_like(feed)
        effluent.T = self.T

        product_mol = float(feed.imol[self.product_ID])
        solid_mol = self.target_recovery * product_mol

        self.design_results["Product recovered to solids (kg/h)"] = (
            solid_mol * feed.chemicals[self.product_ID].MW
        )


class CrystalCentrifuge(bst.units.SolidsCentrifuge):
    """
    Solid/liquid separation of the crystal slurry, in which the wet cake
    carries mother liquor along with the crystals (so the dried product is
    not 100% pure).

    `bst.units.SolidsCentrifuge`'s own `moisture_content` handling only
    moves pure water into the cake, leaving every dissolved species
    (sugars, salts, enzyme, ...) in the liquid-rich outlet, which makes any
    product downstream come out essentially 100% pure by construction. This
    subclass instead sends `product_recovery` of the product to the cake as
    crystals, then entrains just enough *whole* mother liquor (all
    remaining species, water included, at the liquor's own composition) for
    the cake's water content to equal `moisture_content` (a wet-basis mass
    fraction, as in `SolidsCentrifuge`). The impurities dissolved in that
    entrained liquor stay with the crystals through the downstream dryer,
    which removes only water. Costing and sizing are inherited from
    `bst.units.SolidsCentrifuge` unchanged.

    The cake's liquor fraction `phi` follows from the water balance
    `phi * W = moisture_content * (crystals + phi * L)`, where W is the
    feed's water and L its non-crystal mass, so it adjusts automatically to
    the feed composition (e.g. through a recycle loop) instead of being a
    fixed input. The entrained liquor also carries its small share of
    dissolved product, so the cake's product is slightly above
    `product_recovery` * feed.

    Parameters
    ----------
    ins : stream
        Crystal slurry from `CaHPCrystallizer`.
    outs : tuple[stream, stream]
        Crystal cake (crystals + entrained mother liquor) and mother
        liquor.
    product_ID : str
        Chemical ID of the crystallized product (Ca3HP2).
    product_recovery : float
        Fraction of the feed's product recovered to the cake as crystals
        (before the small extra amount carried in the entrained liquor).
    moisture_content : float
        Water mass fraction of the wet cake.
    **kwargs
        Forwarded to `bst.units.SolidsCentrifuge.__init__`.

    See Also
    --------
    Refer to data/3hp.yaml (`crystal_separator`) for the default values and
    references.
    """

    def __init__(
        self,
        ID: str = "",
        ins=None,
        outs=(),
        *,
        product_ID: str = _CRYSTALLIZATION["product_ID"],
        product_recovery: float = _CRYSTALLIZATION["target_recovery"],
        moisture_content: float = _CRYSTAL_SEPARATOR["moisture_content"],
        **kwargs,
    ):
        super().__init__(
            ID, ins, outs,
            split={product_ID: product_recovery},
            moisture_content=moisture_content,
            **kwargs,
        )
        self.product_ID = product_ID
        self.product_recovery = float(product_recovery)

    def _run(self):
        feed, = self.ins
        cake, liquor = self.outs

        cake.empty()
        liquor.empty()
        cake.phase = liquor.phase = "l"
        if feed.isempty():
            return

        rho = self.product_recovery
        m = self.moisture_content
        product_ID = self.product_ID

        crystals_kg = rho * float(feed.imass[product_ID])
        water_kg = float(feed.imass["Water"])
        liquor_feed_kg = float(feed.F_mass) - crystals_kg

        denominator = water_kg - m * liquor_feed_kg
        if denominator <= 0:
            raise RuntimeError(
                f"{self.ID}: feed has too little water ({water_kg:.1f} kg/h of "
                f"{liquor_feed_kg:.1f} kg/h non-crystal mass) to form a cake at "
                f"moisture_content={m}."
            )
        phi = min(m * crystals_kg / denominator, 1.0)

        cake.mol[:] = phi * feed.mol
        cake.imol[product_ID] = float(feed.imol[product_ID]) * (rho + phi * (1.0 - rho))
        liquor.mol[:] = feed.mol - cake.mol
        cake.T = liquor.T = feed.T
        cake.P = liquor.P = feed.P

        self.design_results["Entrained liquor fraction"] = phi
        self.design_results["Entrained liquor (kg/h)"] = phi * liquor_feed_kg


class HydrolysisTank(bst.Unit):
    """
    Enzymatic hydrolysis holding tanks, costed like the continuous
    saccharification tanks of `biorefineries.cellulosic`
    (`ContinuousPresaccharification`, from Humbird et al. 2011, NREL/TP-5100-47764)
    and `biorefineries.HP` (`Saccharification`): large tanks whose cost scales
    with volume to the 0.7 power from the reference tank, rather than many small
    stirred tanks. The reaction itself is applied by a specification on the unit
    (see systems/_3hp_system.py); `_run` only passes the slurry through.

    Parameters
    ----------
    tau : float
        Residence time [h].
    V_max : float
        Largest single tank volume [m3]; the number of tanks is
        ceil(total volume / V_max). Default 3,785 m3 (1,000,000 gal), the
        maximum vessel size of `biorefineries.cellulosic.Saccharification`.
    V_wf : float
        Working volume fraction.
    T : float
        Temperature [K].

    Notes
    -----
    Reference tank: 250,000 gal (946 m3) at $480,000 (3,840,000 / 8 tanks, CE
    522), bare-module factor 2.0, exponent 0.7. Transfer pumps: $47,200 and
    58 kW each at 352 gpm (80 m3/h), exponent 0.8, bare-module factor 2.3
    (cellulosic basis). Sizes beyond the reference are extrapolated. The
    tanks carry no agitation power in that basis.
    """
    _N_ins = 1
    _N_outs = 1
    _F_BM_default = {"Tanks": 2.0, "Transfer pumps": 2.3}
    _units = {"Tank volume": "m3", "Flow rate": "m3/hr"}

    _V_ref = 250e3 * 0.00378541  # m3
    _tank_cost_ref = 3840e3 / 8
    _pump_cost_ref = 47200.0
    _pump_kW_ref = 58.0
    _pump_flow_ref = 352 * 0.2271247  # gpm -> m3/hr
    _CE_ref = 522.0

    def _init(self, tau=96.0, V_max=3785.4118, V_wf=0.9, T=48 + 273.15):
        self.tau = tau
        self.V_max = V_max
        self.V_wf = V_wf
        self.T = T

    def _run(self):
        feed, = self.ins
        effluent, = self.outs
        effluent.copy_like(feed)
        effluent.T = self.T

    def _design(self):
        import math
        v_0 = self.ins[0].F_vol
        total = v_0 * self.tau / self.V_wf
        N = max(1, math.ceil(total / self.V_max))
        N_pumps = max(1, math.ceil(v_0 / self._pump_flow_ref))
        D = self.design_results
        D["Number of tanks"] = N
        D["Tank volume"] = total / N
        D["Total volume"] = total
        D["Number of pumps"] = N_pumps
        D["Flow rate"] = v_0 / N_pumps

    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        scale = bst.CE / self._CE_ref
        C["Tanks"] = D["Number of tanks"] * self._tank_cost_ref * scale * (D["Tank volume"] / self._V_ref) ** 0.7
        q = D["Flow rate"] / self._pump_flow_ref
        C["Transfer pumps"] = D["Number of pumps"] * self._pump_cost_ref * scale * q ** 0.8
        self.power_utility(D["Number of pumps"] * self._pump_kW_ref * q)

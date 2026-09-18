# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Acidogenic-AD-to-microbial-oil system builder for the SaBRe flowsheets.
"""
from pathlib import Path

import biosteam as bst
import flexsolve as flx

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.units import (
    YarrowiaLipidFermenter,
    OilExtraction,
    FermentationMediumTank,
)
from biorefineries.sabre.systems._ad_vfa_system import create_ad_vfa_system
from biorefineries.sabre.systems._biostimulant_system import BIOSTIMULANT_UNIT_IDS
from biorefineries.sabre._tea import create_tea

__all__ = ('create_ad_fermentation_system', 'price_ad_fermentation_system')

# Load assumptions
_FERMENTATION_YAML = load_assumptions("fermentation.yaml")
_VFA_FERM = _FERMENTATION_YAML["vfa"]
_VFA_CASE = _VFA_FERM["cases"][_VFA_FERM["case"]]
_DOWNSTREAM_PROCESSING_YAML = load_assumptions("downstream_processing.yaml")
_VFA_DOWNSTREAM = _DOWNSTREAM_PROCESSING_YAML["oil_extraction"]
_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def _create_vfa_fermentation_system(vfa_broth):
    product_ID = _VFA_CASE["product_ID"]

    ammonia = bst.Stream("fermentation_ammonia")
    ammonia.price = _TEA_PRICE["ammonia"]["baseline"]
    phosphate = bst.Stream("fermentation_phosphate")
    phosphate.price = _TEA_PRICE["phosphate"]["baseline"]
    base = bst.Stream("fermentation_base")
    base.price = _TEA_PRICE["naoh"]["baseline"]
    mgso4 = bst.Stream("fermentation_mgso4")
    mgso4.price = _TEA_PRICE["mgso4"]["baseline"]
    recycle_biomass = bst.Stream("recycle_biomass")

    # -------------------------------------------------
    # Upstream conditioning
    # -------------------------------------------------
    T601 = FermentationMediumTank(
        "T601",
        ins=(vfa_broth, ammonia, phosphate, base, mgso4),
        outs=("conditioned_vfa_broth",),
    )

    # -------------------------------------------------
    # Aerated fermenter (Yarrowia lipolytica)
    # -------------------------------------------------
    # outs[0] (fermentation_vent) is a gas waste stream -- no price needed.
    R601 = YarrowiaLipidFermenter(
        "R601",
        ins=(T601 - 0, recycle_biomass),
        outs=("fermentation_vent", "fermentation_broth"),
    )

    # -------------------------------------------------
    # Post-fermentation oil separation
    # -------------------------------------------------
    V605 = bst.MixTank("V605", ins=R601 - 1, outs=("mixed_fermentation_broth",))
    P606 = bst.Pump("P606", ins=V605 - 0, outs=("pumped_fermentation_broth",))

    Ev607 = bst.MultiEffectEvaporator(
        "Ev607",
        ins=P606 - 0,
        outs=("fermentation_concentrate", "evaporator_vapor"),
        P=(101325, 69682, 47057, 30953),
        V=0.90,
        V_definition="First-effect",
        thermo=(R601.outs[1].thermo.ideal()),
        flash=False,
    )
    Ev607.target_oil_and_solids_content = _VFA_DOWNSTREAM["target_oil_and_solids_content_g_per_L"]
    Ev607.remove_evaporators = True

    P_original = tuple(Ev607.P)
    Pstart = P_original[0]
    Plast = P_original[-1]
    N = len(P_original)

    def concentration_objective(V):
        Ev607.V = V
        Ev607.run()
        effluent = Ev607.outs[0]
        total = effluent.F_mass
        if total <= 0:
            return 0.0
        water = effluent.imass["Water"]
        nonwater_conc = 1000.0 * (1.0 - water / total)
        return Ev607.target_oil_and_solids_content - nonwater_conc

    @Ev607.add_specification(run=False)
    def adjust_evaporation():
        V_last = Ev607.V
        x0 = 0.0
        x1 = 0.5
        Ev607.P = P_original
        Ev607._reload_components = True

        y0 = concentration_objective(x0)
        if y0 <= 0.0:
            Ev607.V = x0
            return
        else:
            Ev607._load_components()
            for i in range(1, N):
                if concentration_objective(1e-6) < 0.0:
                    Ev607.P = tuple(__import__("numpy").linspace(Pstart, Plast, N - 1))
                    Ev607._reload_components = True
                else:
                    break
            y1 = concentration_objective(x1)
            Ev607.V = flx.IQ_interpolation(
                concentration_objective,
                x0, x1, y0, y1,
                x=V_last,
                ytol=1e-5,
                xtol=1e-6,
            )

    P607 = bst.Pump("P607", ins=Ev607 - 0, outs=("pumped_concentrate",), P=101325.0)

    # -------------------------------------------------
    # Cell disruption + oil extraction
    # Slots between P607 and C603_2.
    # Pass-through unit: carries capital + electricity costs for
    # high-pressure homogenization and lipid extraction.
    # Capital anchor: NREL/TP-5100-55431 (2012), Davis et al.
    # -------------------------------------------------
    OE = OilExtraction(
        "OE",
        ins=P607 - 0,
        outs=("extracted_broth",),
        product_ID=product_ID,
        cellmass_ID="CellMass",
    )

    C603_2 = bst.LiquidsSplitCentrifuge(
        "C603_2",
        ins=OE - 0,   # now takes from OE, not P607
        outs=("microbial_oil", "cellmass_plus_aqueous"),
        split={
            product_ID: _VFA_DOWNSTREAM["oil_recovery"],
            "Water": _VFA_DOWNSTREAM["oil_water_split"],
        },
    )
    C603_2.outs[0].price = _TEA_PRICE["microbial_oil"]["baseline"]

    S602 = bst.MockSplitter(
        "S602",
        ins=C603_2 - 1,
        outs=(recycle_biomass, "residual_purge"),
    )
    S602.recycle_total_fraction = _VFA_DOWNSTREAM["recycle_total_fraction"]
    S602.recycle_cellmass_wt_frac = _VFA_DOWNSTREAM["recycle_cellmass_wt_frac"]

    @S602.add_specification(run=True)
    def adjust_biomass_recycle():
        feed = S602.ins[0]
        recycle, purge = S602.outs

        recycle.empty()
        purge.copy_like(feed)

        cellmass_available = float(feed.imass["CellMass"]) if "CellMass" in feed.chemicals else 0.0
        water_available = float(feed.imass["Water"]) if "Water" in feed.chemicals else 0.0

        if feed.F_mass <= 0 or cellmass_available <= 0:
            recycle.empty()
            purge.copy_like(feed)
            return

        target_recycle_mass = S602.recycle_total_fraction * feed.F_mass
        target_cellmass = target_recycle_mass * S602.recycle_cellmass_wt_frac
        cellmass_recycle = min(cellmass_available, target_cellmass)
        water_recycle = cellmass_recycle * (1.0 / S602.recycle_cellmass_wt_frac - 1.0)
        water_recycle = min(water_available, water_recycle)

        recycle.imass["CellMass"] = cellmass_recycle
        recycle.imass["Water"] = water_recycle
        purge.imass["CellMass"] -= cellmass_recycle
        purge.imass["Water"] -= water_recycle
        purge.mol.remove_negatives()
        recycle.T = purge.T = feed.T

    M601 = bst.Mixer(
        "M601",
        ins=(Ev607 - 1, S602 - 1),
        outs=("wastewater",),
    )
    M601.outs[0].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

    # OE added to system path between P607 and C603_2
    sys = bst.System(
        "_vfa_fermentation_sys",
        path=(T601, R601, V605, P606, Ev607, P607, OE, C603_2, S602, M601)
    )

    key_streams = {
        "vfa_broth": vfa_broth,
        "conditioned_vfa_broth": T601.outs[0],
        "recycle_biomass": recycle_biomass,
        "vent": R601.outs[0],
        "fermentation_broth": R601.outs[1],
        "extracted_broth": OE.outs[0],
        "microbial_oil": C603_2.outs[0],
        "residual_slurry": C603_2.outs[1],
        "residual_purge": S602.outs[1],
        "wastewater": M601.outs[0],
        "ammonia": ammonia,
        "phosphate": phosphate,
        "base": base,
        "mgso4": mgso4,
    }

    units = {
        "medium_tank": T601,
        "fermenter": R601,
        "post_mix_tank": V605,
        "post_feed_pump": P606,
        "evaporator": Ev607,
        "fermentation_broth_pump": P607,
        "oil_extraction": OE,
        "oil_centrifuge": C603_2,
        "biomass_recycle_splitter": S602,
        "wastewater_mixer": M601,
    }

    return sys, key_streams, units


def create_ad_fermentation_system(
    feedstock: str | bst.Stream = "pelagic",
    biostimulant_price: float | None = None,
):
    """
    Build the full feedstock-to-product system: raw Sargassum (or an
    already-milled biomass stream) -> acidogenic AD -> VFA broth -> VFA
    fermentation -> microbial oil.

    Parameters
    ----------
    feedstock : str or stream
        Forwarded to create_ad_vfa_system() -- see its docstring for the
        str-vs-stream distinction.
    biostimulant_price : float, optional
        Forwarded to create_ad_vfa_system() -> create_biostimulant_system()
        -- only used when feedstock is a str.

    Returns
    -------
    sys : bst.System
        The full feedstock -> fermentation-product system. Key streams
        are accessible via `sys.flowsheet.stream`, including the AD-VFA
        subsystem's 'vfa_broth' (post-microfiltration permeate),
        'vfa_cake', 'acidogenic_solid_digestate', and the final product
        'microbial_oil'.
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    # add_product_splitter=False: this system needs the raw, unsplit
    # vfa_broth as its own fermentation feed, not ad_vfa_sys's standalone
    # pure_vfa/vfa_disposal split.
    ad_vfa_sys = create_ad_vfa_system(
        feedstock=feedstock, add_product_splitter=False,
        biostimulant_price=biostimulant_price,
    )
    feed = feedstock if not isinstance(feedstock, str) else ad_vfa_sys.feeds[0]
    vfa_broth = ad_vfa_sys.flowsheet.stream.vfa_broth
    vfa_cake = ad_vfa_sys.flowsheet.stream.vfa_cake
    acidogenic_solid_digestate = ad_vfa_sys.flowsheet.stream.acidogenic_solid_digestate

    fer_sys, fer_streams, fer_units = _create_vfa_fermentation_system(vfa_broth=vfa_broth)

    # Fold in both subsystems' units, but not their own HXN facilities -- this
    # system gets its own HXN below, scoped to all units visible here, so nesting
    # either subsystem's narrower one would double-count already-optimized utilities.
    combined_units = [
        u for u in (*ad_vfa_sys.units, *fer_sys.units)
        if not isinstance(u, bst.HeatExchangerNetwork)
    ]
    HXN = bst.HeatExchangerNetwork("HXN", units=tuple(combined_units))
    combined_units.append(HXN)

    sys = bst.System.from_units(
        "ad_fermentation_sys",
        units=combined_units,
    )
    create_tea(sys)

    return sys

# Not used now, kept as a reference
# def _patch_ev607():
#     """
#     Replace Ev607's cost + utility with a low-duty placeholder when V < 0.02.
#     BioSTEAM's MultiEffectEvaporator cost correlation produces nonsensical
#     vessel geometry at near-zero evaporation duty (see
#     legacy_analyses/vfa_fermentation_tea.py::_patch_ev607).
#     """
#     try:
#         ev607 = bst.main_flowsheet.unit["Ev607"]
#     except Exception:
#         return
#     v = getattr(ev607, "V", None)
#     if v is None or v >= 0.02:
#         return
#     feed = ev607.ins[0]
#     feed_m3h = max(feed.F_mass / 1000.0, 1.0)
#     placeholder_usd = 50000.0 * (feed_m3h ** 0.6)
#     for k in list(ev607.baseline_purchase_costs.keys()):
#         ev607.baseline_purchase_costs[k] = 0.0
#     ev607.baseline_purchase_costs["Evaporator (low-duty placeholder)"] = placeholder_usd
#     ev607.heat_utilities.clear()


def price_ad_fermentation_system(credit_tipping_fee: bool = False) -> dict:
    """
    Parameters
    ----------
    credit_tipping_fee : bool
        If False (default): the fixed data/tea.yaml assumption basis --
        biostimulant is priced at its flat `price.biostimulant.baseline`
        assumption, and microbial oil's price is solved on a single
        combined TEA covering every unit in this system (biostimulant's
        Press/PressateConcentrator/Evaporator included), so microbial oil's
        price also recovers a share of biostimulant's own capital.

        If True: biostimulant is instead priced at its own standalone MSP
        (still target IRR), and microbial oil's price is solved against an
        AD-only TEA scope that excludes biostimulant's capital entirely --
        which biostimulant's own standalone price already recovers on its
        own -- plus a tipping-fee credit equal to data/tea.yaml
        `price.disposal_solid.baseline` on the pressed_cake mass taken in,
        since biostimulant would otherwise have paid that same fee to
        dispose of it (whether to a landfill or here).
    """
    from biorefineries.sabre._tea import (
        solve_product_msp,
        disposal_avoided_credit_usd_per_yr, apply_revenue_credit,
    )

    tipping_fee_usd_per_yr = 0.0

    if not credit_tipping_fee:
        # Fixed assumption basis: biostimulant priced at its flat tea.yaml
        # baseline, combined TEA over every unit (biostimulant's own
        # capital included).
        biostimulant_price = _TEA_PRICE["biostimulant"]["baseline"]

        bst.main_flowsheet.clear()
        sys = create_ad_fermentation_system(feedstock="pelagic", biostimulant_price=biostimulant_price)
        sys.simulate()
        # _patch_ev607() not used now, kept as a reference

        product = sys.flowsheet.stream.microbial_oil
        msp = solve_product_msp(tea=sys.TEA, product_stream=product)
    else:
        from biorefineries.sabre.systems._biostimulant_system import price_biostimulant_system

        # Biostimulant's own standalone MSP (target IRR) -- set on the
        # biostimulant_product stream purely for reporting/consistency. It
        # plays no role in this pathway's own product price below.
        biostimulant_price = price_biostimulant_system()["msp_usd_per_kg"]

        bst.main_flowsheet.clear()
        sys = create_ad_fermentation_system(feedstock="pelagic", biostimulant_price=biostimulant_price)
        sys.simulate()
        # _patch_ev607() not used now, kept as a reference

        # AD-only TEA scope: same simulated units, minus the embedded
        # biostimulant subsystem's units, so microbial oil's solved price
        # recovers only this pathway's own capital (Mill, acidogenic AD,
        # screw press, microfilter, fermentation train, HXN) -- not
        # biostimulant's Press/PressateConcentrator/Evaporator, which
        # biostimulant's own price already covers in its own standalone
        # system.
        ad_specific_units = [u for u in sys.units if u.ID not in BIOSTIMULANT_UNIT_IDS]
        ad_specific_sys = bst.System.from_units("ad_fermentation_specific_sys", units=ad_specific_units)
        ad_specific_tea = create_tea(ad_specific_sys)

        product = sys.flowsheet.stream.microbial_oil
        msp = solve_product_msp(tea=ad_specific_tea, product_stream=product)

        tipping_fee_usd_per_yr = disposal_avoided_credit_usd_per_yr(
            sys.flowsheet.stream.pressed_cake, _TEA_PRICE["disposal_solid"]["baseline"], ad_specific_tea,
        )
        msp = apply_revenue_credit(msp, tipping_fee_usd_per_yr)

    return {
        "label": "AD-fermentation",
        "product_desc": "crude microbial oil",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "biostimulant_price_usd_per_kg": biostimulant_price,
        "credit_tipping_fee": credit_tipping_fee,
        "tipping_fee_usd_per_yr": tipping_fee_usd_per_yr,
        'sys': sys,
    }


if __name__ == '__main__':
    results = price_ad_fermentation_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")
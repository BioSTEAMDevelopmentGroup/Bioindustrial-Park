# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
VFA-producing acidogenic AD system builder for the SaBRe flowsheets.
"""
from pathlib import Path

import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import load_assumptions
from biorefineries.sabre.units import (
    AcidogenicDigester,
    DigestateScrewPress,
    Mill,
    VFAMicrofilter,
)
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system, BIOSTIMULANT_UNIT_IDS
from biorefineries.sabre._tea import create_tea

__all__ = ('create_ad_vfa_system', 'price_ad_vfa_system')

_VFA_PRODUCT_SPLITTER = load_assumptions("downstream_processing.yaml")["vfa_product_splitter"]
_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def create_ad_vfa_system(
    feedstock: str | bst.Stream = "pelagic",
    add_product_splitter: bool = True,
    biostimulant_price: float | None = None,
):
    """
    Build the VFA acidogenic AD system: VFA_AD -> SP_VFA -> MF (-> SP_PRODUCT).

    Parameters
    ----------
    feedstock : str or stream
        If a str, it's a feedstock type (data/feedstock.yaml
        `feedstock_type`) and create_biostimulant_system() is called to
        build Press -> PressateConcentrator -> BiostimulantEvaporator,
        followed by a Mill on the pressed cake, to get a milled feed.
        If a stream, it's used directly as the already-milled feed (e.g.
        a splitter-derived stream from systems._ad_integrated_system,
        which builds its own shared preprocessing once).
    add_product_splitter : bool
        If True (default), adds SP_PRODUCT: a splitter on vfa_broth that
        recovers `vfa_product_splitter.vfa_recovery_frac` (data/
        downstream_processing.yaml) of the VFA mass to a `pure_vfa`
        product stream, sending the remainder (unrecovered VFA + everything
        else) to a `vfa_disposal` stream. Set to False when this system is
        built as a component of a larger system (e.g.
        systems._ad_fermentation_system, which needs the raw, unsplit
        `vfa_broth` as its own feed) so that system is unaffected.
    biostimulant_price : float, optional
        Forwarded to create_biostimulant_system() -- only used when
        feedstock is a str (i.e. this system builds its own biostimulant
        subsystem).

    Key outputs (accessible via flowsheet):
        - acidogenic_offgas
        - vfa_broth (post-microfiltration permeate; ready for fermentation)
        - vfa_cake
        - acidogenic_solid_digestate
        - pure_vfa, vfa_disposal (only if add_product_splitter is True)
        - biostimulant_membrane_concentrate (standalone mode only)
        - pressate_permeate (standalone mode only)
        - biostimulant_product (standalone mode only)
    """
    try: bst.settings.get_chemicals()
    except Exception: create_chemicals()
    path = []

    if isinstance(feedstock, str):
        bio_sys = create_biostimulant_system(feedstock=feedstock, biostimulant_price=biostimulant_price)
        # Fold in the biostimulant subsystem's units, but not its own HXN facility --
        # this system gets its own HXN below, scoped to all units visible here, so
        # nesting the subsystem's narrower one would double-count already-optimized
        # utilities.
        path.extend(u for u in bio_sys.units if not isinstance(u, bst.HeatExchangerNetwork))

        # milling_losses: no price -- pure mass loss, not a disposed waste stream.
        ML = Mill("ML", ins=bio_sys.flowsheet.stream.pressed_cake, outs=("milled_biomass", "milling_losses"))
        path.append(ML)
        ad_inlet = ML - 0
    else:
        ad_inlet = feedstock

    # AD.outs[0] (acidogenic_offgas) is a gas waste stream -- no price needed.
    AD = AcidogenicDigester("VFA_AD", ins=ad_inlet, outs=("acidogenic_offgas", "acidogenic_broth"))
    SP = DigestateScrewPress(ID="SP_VFA", ins=AD - 1, outs=("acidogenic_solid_digestate", "raw_vfa_broth"))
    SP.outs[0].price = _TEA_PRICE["disposal_solid"]["baseline"]
    MF = VFAMicrofilter("MF", ins=SP - 1, outs=("vfa_broth", "vfa_cake"))
    # vfa_cake: solid waste (screw-pressed VFA-microfilter retentate).
    MF.outs[1].price = _TEA_PRICE["disposal_solid"]["baseline"]
    path.extend([AD, SP, MF])

    if add_product_splitter:
        vfa_recovery_frac = float(_VFA_PRODUCT_SPLITTER["vfa_recovery_frac"])
        split = {cid: vfa_recovery_frac for cid in MF.vfa_IDs}
        SP_PRODUCT = bst.Splitter(
            "SP_PRODUCT", ins=MF - 0, outs=("pure_vfa", "vfa_disposal"), split=split,
        )
        SP_PRODUCT.outs[0].price = _TEA_PRICE["vfa"]["baseline"]
        SP_PRODUCT.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]
        path.append(SP_PRODUCT)

    HXN = bst.HeatExchangerNetwork("HXN", units=tuple(path))
    path.append(HXN)

    sys = bst.System("ad_vfa_sys", path=tuple(path))
    create_tea(sys)

    return sys


def price_ad_vfa_system(credit_tipping_fee: bool = False) -> dict:
    """
    Parameters
    ----------
    credit_tipping_fee : bool
        If False (default): the fixed data/tea.yaml assumption basis --
        biostimulant is priced at its flat `price.biostimulant.baseline`
        assumption, and VFA's price is solved on a single combined TEA
        covering every unit in this system (biostimulant's Press/
        PressateConcentrator/Evaporator included), so VFA's price also
        recovers a share of biostimulant's own capital.

        If True: biostimulant is instead priced at its own standalone MSP
        (still target IRR), and VFA's price is solved against an AD-only
        TEA scope that excludes biostimulant's capital entirely -- which
        biostimulant's own standalone price already recovers on its own --
        plus a tipping-fee credit equal to data/tea.yaml
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
        sys = create_ad_vfa_system(feedstock="pelagic", biostimulant_price=biostimulant_price)
        sys.simulate()

        product = sys.flowsheet.stream.pure_vfa
        msp = solve_product_msp(tea=sys.TEA, product_stream=product)
    else:
        from biorefineries.sabre.systems._biostimulant_system import price_biostimulant_system

        # Biostimulant's own standalone MSP (target IRR) -- set on the
        # biostimulant_product stream purely for reporting/consistency. It
        # plays no role in this pathway's own product price below.
        biostimulant_price = price_biostimulant_system()["msp_usd_per_kg"]

        bst.main_flowsheet.clear()
        sys = create_ad_vfa_system(feedstock="pelagic", biostimulant_price=biostimulant_price)
        sys.simulate()

        # AD-only TEA scope: same simulated units, minus the embedded
        # biostimulant subsystem's units, so VFA's solved price recovers
        # only this pathway's own capital (Mill, AD, screw press,
        # microfilter, product splitter, HXN) -- not biostimulant's Press/
        # PressateConcentrator/Evaporator, which biostimulant's own price
        # already covers in its own standalone system.
        ad_specific_units = [u for u in sys.units if u.ID not in BIOSTIMULANT_UNIT_IDS]
        ad_specific_sys = bst.System.from_units("ad_vfa_specific_sys", units=ad_specific_units)
        ad_specific_tea = create_tea(ad_specific_sys)

        product = sys.flowsheet.stream.pure_vfa
        msp = solve_product_msp(tea=ad_specific_tea, product_stream=product)

        tipping_fee_usd_per_yr = disposal_avoided_credit_usd_per_yr(
            sys.flowsheet.stream.pressed_cake, _TEA_PRICE["disposal_solid"]["baseline"], ad_specific_tea,
        )
        msp = apply_revenue_credit(msp, tipping_fee_usd_per_yr)

    return {
        "label": "AD-VFA",
        "product_desc": "VFA broth (total-VFA basis)",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "biostimulant_price_usd_per_kg": biostimulant_price,
        "credit_tipping_fee": credit_tipping_fee,
        "tipping_fee_usd_per_yr": tipping_fee_usd_per_yr,
        'sys': sys,
    }


if __name__ == '__main__':
    results = price_ad_vfa_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")
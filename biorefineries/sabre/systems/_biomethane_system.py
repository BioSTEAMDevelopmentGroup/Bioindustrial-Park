# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BIOSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
AD/biomethane system builder for the SaBRe (Sargassum Biorefinery) flowsheets.
"""

from pathlib import Path

import biosteam as bst

from biorefineries.sabre._chemicals import create_chemicals
from biorefineries.sabre.utils import load_assumptions, get_ad_temperature_K
from biorefineries.sabre.units import (
    MethanogenicAD, BiogasUpgrading, H2SRemoval, DigestateScrewPress,
    Mill, HeatingPretreatment, EnzymaticPretreatment, PeroxidePretreatment,
)
from biorefineries.sabre.systems._biostimulant_system import create_biostimulant_system, BIOSTIMULANT_UNIT_IDS
from biorefineries.sabre._tea import create_tea, usd_per_mmbtu_to_usd_per_kg

__all__ = ('create_ad_biomethane_system', 'price_ad_biomethane_system')


# Load assumptions
_PRETREATMENT_AD = load_assumptions("pretreatment.yaml")["pretreatment_ad"]
_AD_YAML = load_assumptions("ad.yaml")
_AD_SHARED = _AD_YAML["ad"]
_AD_PERFORMANCE = _AD_YAML["ad_performance"]
_TEA_PRICE = load_assumptions("tea.yaml")["price"]


def create_ad_biomethane_system(
    feedstock: str | bst.Stream = "pelagic",
    pretreatment_case: str = 'press_mill_only',
    biostimulant_price: float | None = None,
):
    """
    Build the AD/biomethane system: [optional pretreatment] -> AD -> H2S
    removal -> biogas upgrading -> digestate screw press.

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
    pretreatment_case : str
        data/pretreatment.yaml `pretreatment_ad` case name.
    biostimulant_price : float, optional
        Forwarded to create_biostimulant_system() -- only used when
        feedstock is a str (i.e. this system builds its own biostimulant
        subsystem).
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
        ad_feed = ML - 0
    else:
        ad_feed = feedstock

    adS = {**_AD_SHARED, **_AD_SHARED.get("methanogenic", {})}
    adS["temperature_K"] = get_ad_temperature_K(_AD_SHARED, "mesophilic")
    adp = {**_AD_PERFORMANCE, **_AD_PERFORMANCE.get("methanogenic", {})}
    pretreatments = _PRETREATMENT_AD

    pt_units = []
    HT = EZ = PX = None

    if pretreatment_case == "press_mill_only":
        pass

    elif pretreatment_case == "enzymatic":
        EZ = EnzymaticPretreatment("EZ", ins=ad_feed, outs=("enzyme_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.append(EZ)

    elif pretreatment_case == "peroxide":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        ad_feed = PX - 0
        pt_units.append(PX)

    elif pretreatment_case == "combined_PE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=PX - 0, outs=("combined_PE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, EZ])

    elif pretreatment_case == "combined_PTE":
        PX = PeroxidePretreatment("PX", ins=ad_feed, outs=("peroxide_treated_biomass",))
        HT = HeatingPretreatment("HT", ins=PX - 0, outs=("heated_biomass",))
        EZ = EnzymaticPretreatment("EZ", ins=HT - 0, outs=("combined_PTE_treated_biomass",))
        ad_feed = EZ - 0
        pt_units.extend([PX, HT, EZ])

    else:
        raise ValueError(f"Unknown pretreatment_case '{pretreatment_case}'")

    pt_case = pretreatments[pretreatment_case]
    ad_effects = pt_case["ad_effects"]

    AD = MethanogenicAD(
        "AD", ins=ad_feed, outs=("biogas", "digestate"),
        vs_destruction=float(ad_effects["vs_destruction"]),
        ch4_kg_per_kg_vs_fed=float(ad_effects["ch4_kg_per_kg_vs_fed"]),
        raw_biogas_molfrac=dict(ad_effects["raw_biogas_molfrac"]),
        biodegradability=dict(adp["biodegradability"]),
        target_temperature_K=adS["temperature_K"],
    )

    H2SR = H2SRemoval("H2SR", ins=AD - 0, outs=("treated_biogas", "spent_h2s_media"))
    H2SR.outs[1].price = _TEA_PRICE["disposal_solid"]["baseline"]
    # UP.outs[1] (methanogenic_offgas) is a gas waste stream -- no price needed.
    UP = BiogasUpgrading("UP", ins=H2SR - 0, outs=("biomethane", "methanogenic_offgas"))
    # data/tea.yaml `biomethane_mmbtu` is priced on an energy basis; convert
    # to a flat $/kg via mmbtu_per_kg (see _tea.py's usd_per_mmbtu_to_usd_per_kg).
    UP.outs[0].price = usd_per_mmbtu_to_usd_per_kg(
        _TEA_PRICE["biomethane_mmbtu"]["baseline"],
        _TEA_PRICE["biomethane_mmbtu"]["mmbtu_per_kg"],
    )
    SP = DigestateScrewPress(ID="SP", ins=AD - 1, outs=("methanogenic_solid_digestate", "liquid_digestate"))
    SP.outs[0].price = _TEA_PRICE["disposal_solid"]["baseline"]
    SP.outs[1].price = _TEA_PRICE["disposal_wastewater"]["baseline"]

    path.extend([*pt_units, AD, H2SR, UP, SP])

    HXN = bst.HeatExchangerNetwork("HXN", units=tuple(path))
    path.append(HXN)

    sys = bst.System("ad_biomethane_sys", path=path)
    create_tea(sys)

    return sys


def price_ad_biomethane_system(
    pretreatment_case: str = 'press_mill_only',
    credit_tipping_fee: bool = False,
) -> dict:
    """
    Parameters
    ----------
    credit_tipping_fee : bool
        If False (default): the fixed data/tea.yaml assumption basis --
        biostimulant is priced at its flat `price.biostimulant.baseline`
        assumption, and biomethane's price is solved on a single combined
        TEA covering every unit in this system (biostimulant's Press/
        PressateConcentrator/Evaporator included), so biomethane's price
        also recovers a share of biostimulant's own capital.

        If True: biostimulant is instead priced at its own standalone MSP
        (still target IRR), and biomethane's price is solved against an
        AD-only TEA scope that excludes biostimulant's capital entirely --
        which biostimulant's own standalone price already recovers on its
        own -- plus a tipping-fee credit equal to data/tea.yaml
        `price.disposal_solid.baseline` on the pressed_cake mass taken in,
        since biostimulant would otherwise have paid that same fee to
        dispose of it (whether to a landfill or here).
    """
    from biorefineries.sabre._tea import (
        solve_product_msp, CH4_MMBTU_PER_KG,
        disposal_avoided_credit_usd_per_yr, apply_revenue_credit,
    )

    tipping_fee_usd_per_yr = 0.0

    if not credit_tipping_fee:
        # Fixed assumption basis: biostimulant priced at its flat tea.yaml
        # baseline, combined TEA over every unit (biostimulant's own
        # capital included).
        biostimulant_price = _TEA_PRICE["biostimulant"]["baseline"]

        bst.main_flowsheet.clear()
        sys = create_ad_biomethane_system(
            feedstock="pelagic", pretreatment_case=pretreatment_case,
            biostimulant_price=biostimulant_price,
        )
        sys.simulate()

        product = sys.flowsheet.stream.biomethane
        msp = solve_product_msp(
            tea=sys.TEA, product_stream=product,
            energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
        )
    else:
        from biorefineries.sabre.systems._biostimulant_system import price_biostimulant_system

        # Biostimulant's own standalone MSP (target IRR) -- set on the
        # biostimulant_product stream purely for reporting/consistency. It
        # plays no role in this pathway's own product price below.
        biostimulant_price = price_biostimulant_system()["msp_usd_per_kg"]

        bst.main_flowsheet.clear()
        sys = create_ad_biomethane_system(
            feedstock="pelagic", pretreatment_case=pretreatment_case,
            biostimulant_price=biostimulant_price,
        )
        sys.simulate()

        # AD-only TEA scope: same simulated units, minus the embedded
        # biostimulant subsystem's units, so biomethane's solved price
        # recovers only this pathway's own capital (Mill, pretreatment, AD,
        # H2S removal, biogas upgrading, screw press, HXN) -- not
        # biostimulant's Press/PressateConcentrator/Evaporator, which
        # biostimulant's own price already covers in its own standalone
        # system.
        ad_specific_units = [u for u in sys.units if u.ID not in BIOSTIMULANT_UNIT_IDS]
        ad_specific_sys = bst.System.from_units("ad_biomethane_specific_sys", units=ad_specific_units)
        ad_specific_tea = create_tea(ad_specific_sys)

        product = sys.flowsheet.stream.biomethane
        msp = solve_product_msp(
            tea=ad_specific_tea, product_stream=product,
            energy_content_mmbtu_per_kg=CH4_MMBTU_PER_KG,
        )

        tipping_fee_usd_per_yr = disposal_avoided_credit_usd_per_yr(
            sys.flowsheet.stream.pressed_cake, _TEA_PRICE["disposal_solid"]["baseline"], ad_specific_tea,
        )
        msp = apply_revenue_credit(msp, tipping_fee_usd_per_yr)

    return {
        "label": "AD-biomethane",
        "pretreatment_case": pretreatment_case,
        "product_desc": "biomethane (whole-stream basis)",
        "msp_usd_per_kg": msp["usd_per_kg"],
        "annual_product_kg": msp["annual_product_kg"],
        "msp_usd_per_mmbtu": msp["usd_per_mmbtu"],
        "annual_product_mmbtu": msp["annual_product_mmbtu"],
        "biostimulant_price_usd_per_kg": biostimulant_price,
        "credit_tipping_fee": credit_tipping_fee,
        "tipping_fee_usd_per_yr": tipping_fee_usd_per_yr,
        'sys': sys,
    }


if __name__ == '__main__':
    results = price_ad_biomethane_system()
    sys = results['sys']

    figures_dir = Path(__file__).resolve().parent.parent / "results" / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    diagram_path = figures_dir / f"{sys.ID}.svg"
    sys.diagram(file=str(figures_dir / sys.ID), format="svg")
    print(f"System diagram saved to: {diagram_path}")
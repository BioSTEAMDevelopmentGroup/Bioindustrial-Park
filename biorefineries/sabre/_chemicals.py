# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""""
Creating the chemical list used for the SaBRe systems
"""""
import biosteam as bst
import thermosteam as tmo

from biorefineries.sabre.utils import load_assumptions

__all__ = ('create_chemicals',)


def _pseudo_solid(ID: str, MW: float = 1.0, formula: str = None):
    if formula is not None:
        return bst.Chemical(ID, search_db=False, default=True, phase="s", formula=formula)
    return bst.Chemical(ID, search_db=False, default=True, phase="s", MW=MW)

# Structural-carbohydrate constants reused from biorefineries.cellulosic.chemicals,
# which sources Hf from the Humbird et al. 2011 NREL report
# (https://www.nrel.gov/docs/fy11osti/47764.pdf) and Cp from
# https://doi.org/10.1007/s10853-013-7815-6 (lignin/cellulose/hemicellulose
# heat capacities are approximately equal near 350 K).
_cal2joule = 4.184
_Cp_structural = 1.364  # J/g/K
_rho_solids = 1540  # kg/m3, same value used in cellulosic/cane/microalgae


def _structural_solid(ID: str, formula: str, Hf_cal: float, Cp: float = _Cp_structural, rho: float = _rho_solids):
    chemical = bst.Chemical(ID, search_db=False, default=True, phase="s",
                             formula=formula, Hf=Hf_cal * _cal2joule)
    chemical.Cn.add_model(Cp * chemical.MW, top_priority=True)
    chemical.V.add_model(tmo.functional.rho_to_V(rho, chemical.MW), top_priority=True)
    return chemical


def create_chemicals(set_thermo: bool = True):
    Water = bst.Chemical("Water")

    # Sargassum components
    Ash = _pseudo_solid("Ash")
    # Cp and density from biorefineries.cane._chemicals (create_sugarcane_chemicals);
    # Ash is kept as a MW=1 mixture placeholder since it has no single formula.
    Ash.Cn.add_model(0.09 * _cal2joule * Ash.MW, top_priority=True)
    Ash.V.add_model(tmo.functional.rho_to_V(_rho_solids, Ash.MW), top_priority=True)

    # Hf/Cp from biorefineries.cellulosic.chemicals (Humbird et al. 2011 NREL report)
    Glucan = _structural_solid("Glucan", "C6H10O5", -233200)
    Xylan = _structural_solid("Xylan", "C5H8O4", -182100)
    Mannan = Glucan.copy("Mannan")
    Galactan = Glucan.copy("Galactan")
    Arabinan = Xylan.copy("Arabinan")
    Lignin = _structural_solid("Lignin", "C8H8O3", -108248)  # vanillin used as a structural proxy, per biorefineries.cane
    # Same elemental formula/Hf as biorefineries.cellulosic.chemicals; Cp from
    # Leow et al., Green Chem. 2015, 17, 3584-3599 (as in biorefineries.microalgae)
    Protein = _structural_solid("Protein", "CH1.57O0.31N0.29S0.007", -17618, Cp=1.25)

    # Alginic acid uronic-acid repeat unit, dehydrated (C6H10O7 - H2O); no database
    # entry or literature Hf found, so Hf falls back to the pseudo-chemical default (0).
    Alginate = _pseudo_solid("Alginate", formula="C6H8O6")
    # Fucose repeat unit, dehydrated (C6H12O5 - H2O); ignores fucoidan's sulfate
    # substitution. Hf likewise falls back to 0 (no formula-only estimation available).
    Fucoidan = _pseudo_solid("Fucoidan", formula="C6H10O4")
    Mannitol = bst.Chemical("Mannitol")
    OtherSolids = _pseudo_solid("OtherSolids")

    # Gases
    CH4 = bst.Chemical("Methane", phase="g")
    CO2 = bst.Chemical("CarbonDioxide", phase="g")
    H2S = bst.Chemical("HydrogenSulfide", phase="g")
    Oxygen = bst.Chemical("Oxygen", phase="g")
    Nitrogen = bst.Chemical("Nitrogen", phase="g")

    # VFAs
    AceticAcid = bst.Chemical("AceticAcid", phase="l")
    PropionicAcid = bst.Chemical("PropionicAcid", phase="l")
    ButyricAcid = bst.Chemical("ButyricAcid", phase="l")
    ValericAcid = bst.Chemical("ValericAcid", phase="l")
    HexanoicAcid = bst.Chemical("HexanoicAcid", phase="l")

    # VFA fermentation / microbial oil pathway.
    # MicrobialOil: real properties of triolein (TAG), the standard
    # microbial/single-cell-oil proxy used across biorefineries.cane,
    # biorefineries.HP, and biorefineries.OHFA (oleaginous-yeast fermentation
    # yields are tracked directly as TAG in those biorefineries). The
    # database entry for triolein itself has no Hf (confirmed by testing
    # bst.Chemical('Triolein').Hf is None), so biorefineries.cane doesn't
    # rely on a lookup either -- it hardcodes the same literature value
    # (Hf_triolein) used here, in create_acyl_olein().
    MicrobialOil = bst.Chemical("MicrobialOil", search_ID="Triolein")
    MicrobialOil.Hf = -1776e3

    # CellMass: matches biorefineries.actag's 'Cells' chemical (chemicals.yaml),
    # which is explicitly a Yarrowia lipolytica stand-in (synonym
    # 'YarrowiaLipolytica') -- same organism as sabre's fermentation step.
    # Its formula is the generic-yeast composition also used for
    # biorefineries.cane's 'Yeast' chemical, while its Hf is the same
    # Humbird et al. 2011 value used for Z_mobilis (-31169.39 cal/mol);
    # actag pairs that Hf with the yeast formula rather than Z_mobilis's
    # own bacterial formula, which is a better match for Yarrowia than the
    # bacterium-based placeholder used here previously.
    CellMass = bst.Chemical("CellMass", search_db=False, default=True, phase="s",
                             formula="CH1.61O0.56N0.16", Hf=-31169.39 * _cal2joule)
    CellMass.V.add_model(tmo.functional.rho_to_V(_rho_solids, CellMass.MW), top_priority=True)

    # Conditioner / nutrient additions -- real database chemicals. Ammonia's
    # natural phase_ref is gas; forced to liquid since it's dosed as an
    # aqueous nutrient stream (matches the NH3 handling in
    # biorefineries.cornstover/cellulosic).
    Ammonia = bst.Chemical("Ammonia")
    Ammonia.at_state("l")
    KH2PO4 = bst.Chemical("KH2PO4")
    NaOH = bst.Chemical("NaOH")
    MagnesiumSulfate = bst.Chemical("MagnesiumSulfate")

    # These are dosed as dilute solutes in aqueous nutrient streams, so their
    # own volume is treated as negligible (rho=1e5 kg/m3 proxy), matching
    # biorefineries.cane's treatment of HCl/NaOH. This also works around
    # KH2PO4 and MagnesiumSulfate's database liquid-volume correlations
    # ('NEGLECT_P') raising at ambient T/P -- confirmed by testing that
    # bst.Chemical('KH2PO4').V('l', 298.15, 101325) and the MagnesiumSulfate
    # equivalent both fail without this override.
    for _nutrient in (NaOH, KH2PO4, MagnesiumSulfate):
        _nutrient.V.l.add_model(tmo.functional.rho_to_V(1e5, _nutrient.MW), top_priority=True)

    # KH2PO4 and MagnesiumSulfate's database liquid-viscosity correlations
    # ('NEGLECT_P') are likewise invalid at process T/P -- confirmed by
    # testing bst.Chemical('KH2PO4').mu('l', 298.15, 101325) and the
    # MagnesiumSulfate equivalent. As dilute solutes, approximate their
    # solution viscosity contribution with Water's (NaOH's own liquid mu
    # model is valid, so it's left alone).
    for _nutrient in (KH2PO4, MagnesiumSulfate):
        _nutrient.copy_models_from(Water, ["mu"])

    chems = bst.Chemicals([
        Water,
        Ash, Protein, Lignin, Glucan, Xylan, Mannan, Galactan, Arabinan,
        Alginate, Fucoidan, Mannitol, OtherSolids,
        CH4, CO2, H2S, Oxygen, Nitrogen,
        AceticAcid, PropionicAcid, ButyricAcid, ValericAcid, HexanoicAcid,
        MicrobialOil, CellMass,
        Ammonia, KH2PO4, NaOH, MagnesiumSulfate,
    ])
    chems.compile()

    solids_IDs = load_assumptions("feedstock.yaml").get("solids_IDs", [])
    solids_IDs = [i for i in solids_IDs if i in chems]
    if solids_IDs:
        chems.define_group("solids", solids_IDs)

    if set_thermo:
        thermo = tmo.Thermo(chems)
        bst.settings.set_thermo(thermo)
    return chems
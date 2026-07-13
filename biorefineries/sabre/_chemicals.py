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
from __future__ import annotations

import biosteam as bst
import thermosteam as tmo

__all__ = ('create_chemicals', 'set_thermo')


def _pseudo_chemical(ID: str, phase: str, MW: float = 1.0):
    return bst.Chemical(ID, search_db=False, default=True, phase=phase, MW=MW)


def _pseudo_solid(ID: str, MW: float = 1.0):
    return _pseudo_chemical(ID, phase="s", MW=MW)


def _pseudo_liquid(ID: str, MW: float = 1.0):
    return _pseudo_chemical(ID, phase="l", MW=MW)


def _pseudo_gas(ID: str, MW: float = 1.0):
    return _pseudo_chemical(ID, phase="g", MW=MW)


def create_chemicals():
    Water = bst.Chemical("Water")

    # Sargassum components
    Ash = _pseudo_solid("Ash")
    Protein = _pseudo_solid("Protein")
    Lignin = _pseudo_solid("Lignin")
    Glucan = _pseudo_solid("Glucan")
    Xylan = _pseudo_solid("Xylan")
    Mannan = _pseudo_solid("Mannan")
    Galactan = _pseudo_solid("Galactan")
    Arabinan = _pseudo_solid("Arabinan")
    Alginate = _pseudo_solid("Alginate")
    Fucoidan = _pseudo_solid("Fucoidan")
    Mannitol = _pseudo_solid("Mannitol")
    OtherSolids = _pseudo_solid("OtherSolids")

    # HTL lumped pseudo-components
    HTLBiocrude = _pseudo_liquid("HTLBiocrude")
    HTLAqueousOrg = _pseudo_liquid("HTLAqueousOrg")
    HTLGas = _pseudo_gas("HTLGas")
    HTLChar = _pseudo_solid("HTLChar")

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

    # VFA fermentation / microbial oil pathway
    MicrobialOil = _pseudo_liquid("MicrobialOil")
    CellMass = _pseudo_solid("CellMass", MW=24.0)
    SolubleOrganics = _pseudo_liquid("SolubleOrganics")

    # Conditioner / nutrient additions
    Ammonia = _pseudo_liquid("Ammonia", MW=17.031)
    KH2PO4 = _pseudo_solid("KH2PO4", MW=136.086)
    NaOH = _pseudo_solid("NaOH", MW=40.0)
    MagnesiumSulfate = _pseudo_solid("MagnesiumSulfate", MW=120.366)

    chems = bst.Chemicals([
        Water,
        Ash, Protein, Lignin, Glucan, Xylan, Mannan, Galactan, Arabinan,
        Alginate, Fucoidan, Mannitol, OtherSolids,
        HTLBiocrude, HTLAqueousOrg, HTLGas, HTLChar,
        CH4, CO2, H2S, Oxygen, Nitrogen,
        AceticAcid, PropionicAcid, ButyricAcid, ValericAcid, HexanoicAcid,
        MicrobialOil, CellMass, SolubleOrganics,
        Ammonia, KH2PO4, NaOH, MagnesiumSulfate,
    ])
    chems.compile()

    thermo = tmo.Thermo(chems)
    bst.settings.set_thermo(thermo)
    tmo.settings.set_thermo(thermo)
    return chems


def set_thermo():
    return create_chemicals()

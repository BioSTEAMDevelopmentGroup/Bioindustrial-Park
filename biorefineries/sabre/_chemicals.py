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

__all__ = ('create_chemicals',)


# Structural-carbohydrate constants reused from biorefineries.cellulosic.chemicals,
# which sources Hf from the Humbird et al. 2011 NREL report
# (https://www.nrel.gov/docs/fy11osti/47764.pdf) and Cp from
# https://doi.org/10.1007/s10853-013-7815-6 (lignin/cellulose/hemicellulose
# heat capacities are approximately equal near 350 K).
_cal2joule = 4.184
_Cp_structural = 1.364  # J/g/K
_rho_solids = 1540  # kg/m3, same value used in cellulosic/cane/microalgae


def _structural_solid(ID: str, formula: str, Hf_cal: float,
                       Cp: float = _Cp_structural, rho: float = _rho_solids):
    chemical = bst.Chemical(ID, search_db=False, default=True, phase="s",
                             formula=formula, Hf=Hf_cal * _cal2joule)
    chemical.Cn.add_model(Cp * chemical.MW, top_priority=True)
    chemical.V.add_model(tmo.functional.rho_to_V(rho, chemical.MW), top_priority=True)
    return chemical


def create_chemicals(set_thermo: bool = True, include_hp3: bool = False):
    """
    Parameters
    ----------
    include_hp3 : bool
        If True, also construct and include the chemicals used only by
        the EnzymaticPress/3-HP-from-Sargassum pathway (Glucose,
        AlginateMonomer, Enzyme, ...). False by default so every existing
        system (biostimulant, Biomethane, VFA, Microbial oil)
        keeps exactly the chemical roster -- and therefore the exact
        floating-point behavior -- it always has. Adding chemicals to a
        shared Chemicals object changes NumPy's pairwise-summation array
        length for every mixture-property calculation package-wide, which
        can nudge iterative solvers (recycle convergence, MSP root-finding)
        by a small but nonzero amount even for chemicals with zero flow in
        a given system -- confirmed by testing that this is insensitive to
        *where* in the chemical list the new species are inserted. Keeping
        the 3-HP-only chemicals out of the default roster avoids that
        drift for every pathway that doesn't need them.
    """
    Water = bst.Chemical("Water")

    # Sargassum components
    # Ash modeled as CaO (following biorefineries.lactic's precedent for
    # generic biomass ash). CaO's database liquid-viscosity correlation
    # ('NEGLECT_P') is invalid at process T/P, the same issue this file
    # already works around for KH2PO4/MagnesiumSulfate below,
    # approximated with Water's viscosity for the same reason.
    #
    # HHV=0/LHV=0 passed to the constructor do NOT actually take effect
    # for a real database chemical -- CaO.Hf/HHV/LHV all still come back
    # as its real (nonzero) database value of -634900 J/mol (confirmed by
    # direct testing). This was dormant/harmless everywhere Ash's moles
    # are only ever split/mixed, never reacted (every non-3hp sabre
    # pathway), since Hf cancels out across any non-reactive step
    # regardless of its value. It stopped being harmless once
    # sabre/systems/_3hp_system.py added a real BoilerTurbogenerator:
    # thermosteam auto-generates a combustion "reaction" for every
    # chemical (even non-combustible ones), and for CaO that comes out as
    # the nonsensical `Ash -> 0.5 Oxygen` (thermosteam's default
    # combustion-reaction heuristic assumes any formula oxygen in excess
    # of what's needed to oxidize C/H/S/N is released as O2, which isn't
    # physically real for a metal oxide) -- with CaO's real, large
    # negative Hf, BT's simulated energy balance saw this as a
    # substantial *endothermic* heat sink (~318 GJ/h at this pathway's
    # scale) that isn't physically real. Explicitly zeroing Hf (not just
    # the ineffective HHV/LHV kwargs) fixes both: Ash becomes a true
    # energy-inert pass-through regardless of what nonsense reaction
    # thermosteam invents for it.
    CaO = bst.Chemical("CaO", phase="s", HHV=0, LHV=0)
    CaO.Hf = 0.0
    CaO.copy_models_from(Water, ["mu"])
    Ash = CaO.copy(ID="Ash")

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

    # Alginic acid uronic-acid repeat unit, dehydrated (C6H10O7 - H2O); MW is
    # computed from the formula. No database entry or literature Hf found for
    # alginic acid or its mannuronic/guluronic acid monomers (checked NIST
    # WebBook and general literature), so Hf borrows Glucan's structural-
    # carbohydrate value as the least-bad available proxy; Cn/V likewise
    # reuse Glucan's Cp/density basis rather than bst.Chemical's generic
    # (and here badly wrong -- implied density ~1e6 kg/m3) default estimates.
    Alginate = _structural_solid("Alginate", "C6H8O6", -233200)
    # Fucose repeat unit, dehydrated (C6H12O5 - H2O); ignores fucoidan's sulfate
    # substitution. No literature Hf found for fucoidan or its L-fucose
    # monomer either, so Hf and Cn/V all borrow Glucan's basis for the same
    # reason as Alginate.
    Fucoidan = _structural_solid("Fucoidan", "C6H10O4", -233200)
    Mannitol = bst.Chemical("Mannitol")
    # No known formula/composition for this lumped catch-all; treated as
    # generically Glucan-like (same MW/Cn/V/Hf), consistent with how
    # Mannan/Galactan already borrow Glucan's properties above.
    OtherSolids = Glucan.copy("OtherSolids")

    # EnzymaticPress/3-HP-only chemicals (sabre/units/_preprocessing.py).
    # Gated behind include_hp3 -- see create_chemicals()'s docstring for
    # why these aren't unconditionally in the default roster.
    hp3_chemicals = []
    if include_hp3:
        # Real database chemical -- soluble hydrolysis product of Glucan
        # (Glucan + Water -> Glucose is atom-balanced against Glucan's
        # C6H10O5 formula above).
        Glucose = bst.Chemical("Glucose")

        # Free monomer form of Alginate (mannuronic/guluronic acid, hydrated:
        # Alginate's C6H8O6 repeat unit + H2O -> C6H10O7), formed by
        # enzymatic hydrolysis in EnzymaticPress. No literature Hf found
        # for the free monomer either, so Hf/Cp borrow the same generic
        # structural-carbohydrate basis already used for Alginate above.
        # Unlike Alginate, this is modeled as phase="l" (dissolved, not
        # solid-locked) since it must follow the "solubles" split path
        # once hydrolyzed, not the "solids" path -- see
        # utils.get_solids_group_IDs(). Liquid density is a rough
        # concentrated-sugar-acid-solution proxy (1200 kg/m3, vs. Glucan's
        # solid-phase 1540 kg/m3), and viscosity borrows Water's model for
        # the same NEGLECT_P-workaround reason as the nutrient salts below.
        AlginateMonomer = bst.Chemical("AlginateMonomer", search_db=False, default=True,
                                        phase="l", formula="C6H10O7", Hf=-233200 * _cal2joule)
        AlginateMonomer.Cn.add_model(_Cp_structural * AlginateMonomer.MW, top_priority=True)
        AlginateMonomer.V.add_model(tmo.functional.rho_to_V(1200, AlginateMonomer.MW), top_priority=True)
        AlginateMonomer.copy_models_from(Water, ["mu"])

        # Generic-protein proxy for the hydrolysis enzyme cocktail (cellulase +
        # alginate lyase, dosed in EnzymaticPress) -- same formula/Hf
        # convention biorefineries.cellulosic uses for its own 'Enzyme'
        # chemical (CH1.59O0.42N0.24S0.01, Hf shared with Protein-like
        # chemicals), also standing in for alginate lyase since no separate
        # property data exists for it. Modeled as phase="l" (dosed as a dilute
        # aqueous enzyme prep, not solid-locked), with Water-like density/
        # viscosity as a dilute-solute proxy (same technique used for
        # KH2PO4/MagnesiumSulfate below).
        Enzyme = bst.Chemical("Enzyme", search_db=False, default=True, phase="l",
                               formula="CH1.59O0.42N0.24S0.01", Hf=-17618 * _cal2joule)
        Enzyme.Cn.add_model(1.25 * Enzyme.MW, top_priority=True)
        Enzyme.V.add_model(tmo.functional.rho_to_V(1000, Enzyme.MW), top_priority=True)
        Enzyme.copy_models_from(Water, ["mu"])

        # Real database chemical: 3-hydroxypropionic acid, the fermentation
        # product (HPFermentation, sabre/units/_3hp.py).
        HP = bst.Chemical("HP", search_ID="3-Hydroxypropionic acid")

        # Real database chemical: fermentation neutralization base, forming
        # Ca3HP2 in situ (HPFermentation). Dosed into a liquid-phase
        # broth, but its database liquid-volume/viscosity correlations
        # ('NEGLECT_P') are invalid at process T/P -- same issue this file
        # already works around for KH2PO4/MagnesiumSulfate/NaOH below, fixed
        # the same way (dilute-solute density/viscosity proxy).
        CalciumDihydroxide = bst.Chemical("CalciumDihydroxide", search_ID="Calcium hydroxide")
        CalciumDihydroxide.V.l.add_model(
            tmo.functional.rho_to_V(1e5, CalciumDihydroxide.MW), top_priority=True
        )
        CalciumDihydroxide.copy_models_from(Water, ["mu"])

        # Calcium 3-hydroxypropionate, the crystallized final product
        # (CaHPCrystallizer, sabre/units/_3hp.py). 3-HP and lactic acid are
        # structural isomers (both C3H6O3), so their calcium salts share the
        # exact same molecular formula (CaC6H10O6) -- properties borrowed
        # from the database's real "Calcium lactate" entry as the best
        # available proxy (same borrowing technique used for Alginate/
        # Fucoidan above, but via a formula-identical real analog here
        # rather than an approximate one). Locked to phase="s" (the
        # crystallized/dried state this chemical is always tracked in);
        # the database entry has no Hf, so one is estimated assuming the
        # neutralization reaction (2 HP + Ca(OH)2 -> Ca3HP2 + 2 H2O) is
        # thermoneutral (no literature heat of reaction found):
        # Hf(Ca3HP2) = 2*Hf(HP) + Hf(Ca(OH)2) - 2*Hf(Water). The database
        # entry also has no solid-phase molar volume model, so density
        # borrows the same generic structural-solid proxy (_rho_solids)
        # used for Glucan/Alginate/etc. above. Its database liquid-
        # viscosity correlation ('NEGLECT_P') is invalid outside
        # atmospheric pressure (confirmed by testing: a downstream Pump
        # operating this stream under vacuum raised a RuntimeError even
        # though this chemical is locked to phase="s") -- same NEGLECT_P
        # issue this file already works around for KH2PO4/MagnesiumSulfate/
        # NaOH/CalciumDihydroxide above, fixed the same way.
        Ca3HP2 = bst.Chemical("Ca3HP2", search_ID="Calcium lactate", phase="s")
        Ca3HP2.Hf = 2 * HP.Hf + CalciumDihydroxide.Hf - 2 * Water.Hf
        Ca3HP2.V.add_model(tmo.functional.rho_to_V(_rho_solids, Ca3HP2.MW), top_priority=True)
        Ca3HP2.copy_models_from(Water, ["mu"])

        # Combustion product of the sulfur atoms in Protein/Enzyme's
        # formulas (both already S-bearing), needed only because
        # bst.BoilerTurbogenerator (sabre/systems/_3hp_system.py's BT,
        # burning pressed_cake/milling_losses/cell_mass) calls
        # chemicals.get_combustion_reactions() at simulation time, which
        # raises UndefinedChemicalAlias('SO2') without this -- confirmed
        # by direct testing. Gated behind include_hp3 like the rest of
        # this block (not the default roster) even though Protein/Enzyme's
        # sulfur isn't itself 3-HP-specific, since no other sabre system
        # has a BT yet and adding any chemical to the default roster
        # drifts the other pathways' regression tests (see this function's
        # own docstring).
        SO2 = bst.Chemical("SO2", phase="g")

        # Same issue, for the phosphorus in KH2PO4's formula (KH2PO4 is
        # already in the default roster below, for nutrient dosing, but
        # this combustion byproduct is only ever needed by BT, so it's
        # gated here rather than added alongside KH2PO4 itself). Never
        # actually flows in this pathway (HPFermentation doesn't dose
        # KH2PO4 at all -- only biorefineries.sabre's microbial_oil
        # pathway does), but get_combustion_reactions() builds a reaction
        # for every chemical in the compiled roster regardless of whether
        # it's ever present in a real stream, and the database's real
        # entry is missing Psat/Tb/Hvap (confirmed by testing) -- since
        # the exact values are irrelevant for a chemical that never flows,
        # Water's own models are borrowed as a placeholder rather than
        # sourcing real P4O10 property data.
        P4O10 = bst.Chemical("P4O10", phase="s")
        P4O10.copy_models_from(Water, ["Psat", "Hvap", "V", "mu"])
        P4O10.Tb = Water.Tb

        # BT's flue-gas desulfurization reaction
        # (SO2 + Ca(OH)2 + 0.5 O2 -> CaSO4 + H2O) needs gypsum (CaSO4) as
        # a product chemical too; same missing-property gaps as P4O10,
        # fixed the same way (never flows in bulk, so the borrowed values
        # are inconsequential).
        CaSO4 = bst.Chemical("CaSO4", phase="s")
        CaSO4.copy_models_from(Water, ["Psat", "Hvap", "mu"])
        CaSO4.Tb = Water.Tb

        hp3_chemicals = [
            Glucose, AlginateMonomer, Enzyme, HP, CalciumDihydroxide, Ca3HP2,
            SO2, P4O10, CaSO4,
        ]

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

    # Conditioner / nutrient additions. Ammonia's
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
    # ('NEGLECT_P') raising at ambient T/P.
    for _nutrient in (NaOH, KH2PO4, MagnesiumSulfate):
        _nutrient.V.l.add_model(tmo.functional.rho_to_V(1e5, _nutrient.MW), top_priority=True)

    # KH2PO4 and MagnesiumSulfate's database liquid-viscosity correlations
    # ('NEGLECT_P') are likewise invalid at process T/P.
    # As dilute solutes, approximate their
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
        *hp3_chemicals,
    ])
    chems.compile()

    # "solids" group derived directly from which chemicals are actually
    # locked to phase='s' above -- the single source of truth for which
    # chemicals count as solids in unit simulations (Press,
    # VFAMicrofilter, DigestateDecanterCentrifuge via
    # utils.get_solids_group_IDs()) is the chemical models themselves, not
    # a separately-maintained ID list that could drift out of sync with them.
    solid_IDs = [c.ID for c in chems if c.locked_state == "s"]
    if solid_IDs:
        chems.define_group("solids", solid_IDs)

    if set_thermo:
        thermo = tmo.Thermo(chems)
        bst.settings.set_thermo(thermo)
    return chems